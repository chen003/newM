#include "sp_shell.h"
struct JschemeConfig { //存单粒子轨道
  int j;
  int l;
};
struct MschemeConfig {
  int m;
  int JschemeID; //从 1 开始
};
struct JschemeConfig Jscheme[n_orbtJ], JschemeH[n_orbtJH];
struct MschemeConfig Mscheme[n_orbtM], MschemeH[n_orbtMH];

void setMJID(int n_J, int n_M, FILE *fp, struct JschemeConfig *Jscheme,
             struct MschemeConfig *Mscheme) {
  fseek(fp, 0, SEEK_SET);
  for (int i = 0; i < n_J; i++) {
    fscanf(fp, "%d ", &(Jscheme[i].j));
  }
  for (int i = 0; i < n_J; i++) {
    fscanf(fp, "%d ", &(Jscheme[i].l));
  }
  int NinJ = 0;
  int Jid = 0;
  for (int i = 0; i < n_M; i++) {
    Mscheme[i].JschemeID = Jid + 1;
    Mscheme[i].m = Jscheme[Jid].j - 2 * NinJ;
    NinJ++;
    if (NinJ == Jscheme[Jid].j + 1) {
      Jid++;
      NinJ = 0;
    }
  }
  fseek(fp, 0, SEEK_SET);
}

void setMschemeID2JschemeID(
    FILE *fp,
    FILE *fpH) { // 从 fp 和 fpH 读取轨道角动量等数据, 让 jfromA 等能运行
  setMJID(n_orbtJ, n_orbtM, fp, Jscheme, Mscheme);
  setMJID(n_orbtJH, n_orbtMH, fpH, JschemeH, MschemeH);
}

int JschemeJ(int type, int a) { // type 是粒子类型, A 是 J-scheme 轨道标号
  if (type == 1 || type == 2) {
    return Jscheme[a - 1].j;
  } else if (type == 3) {
    return JschemeH[a - 1].j;
  } else {
    fprintf(stderr, "JschemeJ type error.\n");
    exit(1);
  }
}

int jfromA(int type, int A) { // type 是粒子类型, A 是 M-scheme 轨道标号
  if (type == typepp || type == typenn) {
    return Jscheme[Mscheme[A - 1].JschemeID - 1].j;
  } else if (type == typell) {
    return JschemeH[MschemeH[A - 1].JschemeID - 1].j;
  } else {
    fprintf(stderr, "jfromA type error.\n");
    exit(1);
  }
}

int mfromA(int type, int a) {
  if (type == typepp || type == typenn) {
    return Mscheme[a - 1].m;
  } else if (type == typell) {
    return MschemeH[a - 1].m;
  } else {
    fprintf(stderr, "mfromA type error.\n");
    exit(1);
  }
}

int lfromA(int type, int a) {
  if (type == typepp || type == typenn) {
    return Jscheme[Mscheme[a - 1].JschemeID - 1].l;
  } else if (type == typell) {
    return JschemeH[MschemeH[a - 1].JschemeID - 1].l;
  } else {
    fprintf(stderr, "lfromA type error.\n");
    exit(1);
  }
}

int typeAfromAB(int type, int lr) {
  switch (type) {
  case 4:          // pn 相互作用
    if (lr == 1) { // lr 1 代表左, 2 代表右
      return 1;
    } else if (lr == 2) {
      return 2;
    }
    break;
  case 1:
    return 1;
  case 2:
    return 2;
  case 5: // nΛ 相互作用, 等价于 pΛ
    if (lr == 1) {
      return 1;
    } else if (lr == 2) {
      return 3;
    }
    break;
  case 6: // Λp 相互作用
    if (lr == 1) {
      return 3;
    } else if (lr == 2) {
      return 1;
    }
    break;
  default:
    fprintf(stderr, "%d type doesn't exist in typeAfromAB.\n", type);
    exit(1);
    break;
  }
}

int numfromA(
    int A,
    int type) { // 从 Mschmeme 单粒子编号 得到 Jscheme 编号. type 是粒子类型
  if (type == 1 || type == 3) {
    return Mscheme[A - 1].JschemeID;
  } else if (type == 2) {
    return Mscheme[A - 1].JschemeID + n_orbtJ;
  } else {
    fprintf(stderr, "numfromA error");
    exit(1);
  }
}

int moveToNextLine(FILE *fp) {
  for (int x = 0; x < 1; x++) {
    char temp[255 + 1];          //一行最长1024个字符, 255随手写的
    if (!fgets(temp, 255, fp)) { // fgets 会在读到新的一行时停下来.
      break;
    }
  }
  return 0;
}

// base 进制转换
void itoa(int n, int base, char *buf) {
  div_t dv = {.quot = n};
  char *p = buf;
  do {
    dv = div(dv.quot, base);
    *p++ = "0123456789abcdef"[abs(dv.rem)];
  } while (dv.quot);
  if (n < 0)
    *p++ = '-';
  *p-- = '\0';
  while (buf < p) {
    char c = *p;
    *p-- = *buf;
    *buf++ = c;
  } // reverse
}
