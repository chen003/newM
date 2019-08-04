#include "sp_shell.h"
#include <math.h>
#include <stdio.h>

#define nlines_T 100

#define PATHint                                                                \
  "../newM/gogny" // JT表象 TBME 存放的文件夹, 生成的文件也在此路径, 与 main.c
                  // 中的 PATHgogny 相同

// replace A3Z1

struct I4 {
  int i1;
  int i2;
  int i3;
  int i4;
  int J;
  int T;
  double me;
};

void getJTME(FILE *fp, struct I4 *p) { //从文件读取一行 JT 表象的矩阵元
  fscanf(fp, "%d  %d  %d  %d  %d  %d  %lf\n", &p->i1, &p->i2, &p->i3, &p->i4,
         &p->J, &p->T, &p->me);
}

bool ConfigEqualQ(struct I4 p, struct I4 q) {
  if (p.i1 == q.i1 && p.i2 == q.i2 && p.i3 == q.i3 && p.i4 == q.i4 &&
      p.J == q.J)
    return true;
  else
    return false;
}

// 转换为 pp 和 nn TBME 的标记
void transVal(struct I4 *p, struct I4 *q) {
  p->i1 = q->i1;
  p->i2 = q->i2;
  p->i3 = q->i3;
  p->i4 = q->i4;
  p->J = q->J;
  p->T = q->T;
  p->me = q->me;
}

void transValn(struct I4 *p, struct I4 *q) {
  p->i1 = q->i1 + n_orbtJ;
  p->i2 = q->i2 + n_orbtJ;
  p->i3 = q->i3 + n_orbtJ;
  p->i4 = q->i4 + n_orbtJ;
  p->J = q->J;
  p->T = q->T;
  p->me = q->me;
}

// JT 表象的 TBME 有四种交换方式(有对称关系), 相应可以生成四种可能的 pn 表象.
// 记为 pn, pnex, pnex1, pnex2
void transValpn(struct I4 *p, struct I4 *q) {
  p->i1 = q->i1;
  p->i2 = q->i2 + n_orbtJ;
  p->i3 = q->i3;
  p->i4 = q->i4 + n_orbtJ;
  p->J = q->J;
  p->T = 0;
  p->me = q->me;
}

void transValpnex1(struct I4 *p, struct I4 *q) {
  p->i1 = q->i2;
  p->i2 = q->i1 + n_orbtJ;
  p->i3 = q->i4;
  p->i4 = q->i3 + n_orbtJ;
  p->J = q->J;
  p->T = 0;
  p->me = q->me;
}

void transValpnex2(struct I4 *p, struct I4 *q) {
  p->i1 = q->i2;
  p->i2 = q->i1 + n_orbtJ;
  p->i3 = q->i3;
  p->i4 = q->i4 + n_orbtJ;
  p->J = q->J;
  p->T = 0;
  p->me = q->me;
}

void transValpnex(struct I4 *p,
                  struct I4 *q) { // pp nn 是全同粒子, 在算M-scheme TBME
  // 时可以用对称关系, 但 pn 不方便直接用,
  // 所以在这里转换表象的时候用
  p->i1 = q->i1;
  p->i2 = q->i2 + n_orbtJ;
  p->i3 = q->i4;
  p->i4 = q->i3 + n_orbtJ;
  p->J = q->J;
  p->T = 0;
  p->me = q->me;
}

// 输出 pn 表象的 TBME. T = 0 表示来自两项求和.
void print2file(FILE *fpw, struct I4 p) {
  fprintf(fpw, "%d  %d  %d  %d  %d  %d  %f\n", p.i1, p.i2, p.i3, p.i4, p.J, p.T,
          p.me);
}

void factorize(struct I4 *p) {
  if (p->i1 == p->i2 - n_orbtJ && p->i3 == p->i4 - n_orbtJ) {
    p->me = p->me;
  } else if ((p->i1 != p->i2 - n_orbtJ && p->i3 == p->i4 - n_orbtJ) ||
             (p->i1 == p->i2 - n_orbtJ && p->i3 != p->i4 - n_orbtJ)) {
    p->me = sqrt(1. / 2) * p->me;
  } else if (p->i1 != p->i2 - n_orbtJ && p->i3 != p->i4 - n_orbtJ) {
    p->me = p->me / 2;
  }
}

void convert1line(FILE *fpw, struct I4 *p, int *j) {
  struct I4 r;
  if (p->T == 1) { // pp 和 nn 情形
    print2file(fpw, *p);
    transValn(&r, p);
    print2file(fpw, r);
  }
  double tot = p->me;
  transValpn(&r, p);
  r.me = tot;
  factorize(&r);
  print2file(fpw, r);
  if (p->i3 != p->i4) {
    transValpnex(&r, p); // 交换 34
    r.me = -((((j[p->i3 - 1] + j[p->i4 - 1]) / 2 - p->J + 1 - p->T) % 2 == 0)
                 ? 1
                 : -1) *
           tot;
    factorize(&r);
    print2file(fpw, r);
  }
  if (p->i1 != p->i2) { // 1234 换 2143. 当 p 和 n 在不同 k 轨道时交换后输出的
    if (p->i1 != p->i3 || p->i2 != p->i4) { // km 态不同, 是不能交换出来的.
      transValpnex2(&r, p);
      r.me = -((((j[p->i1 - 1] + j[p->i2 - 1]) / 2 - p->J + 1 - p->T) % 2 == 0)
                   ? 1
                   : -1) *
             tot;
      factorize(&r);
      print2file(fpw, r);
    }
    if (p->i3 != p->i4) {
      if (p->i1 != p->i4 || p->i2 != p->i3) {
        transValpnex1(&r, p);
        r.me = ((((j[p->i1 - 1] + j[p->i2 - 1]) / 2 - p->J + 1 - p->T) % 2 == 0)
                    ? 1
                    : -1) *
               ((((j[p->i3 - 1] + j[p->i4 - 1]) / 2 - p->J + 1 - p->T) % 2 == 0)
                    ? 1
                    : -1) *
               tot;
        factorize(&r);
        print2file(fpw, r);
      }
    }
  }
}

void convert2lines(FILE *fpw, struct I4 *p, struct I4 *q, int *j) {
  struct I4 r;
  if (p->T == 1) { // pp 和 nn 情形
    print2file(fpw, *p);
    transValn(&r, p);
    print2file(fpw, r);
  } else if (q->T == 1) {
    print2file(fpw, *q);
    transValn(&r, q);
    print2file(fpw, r);
  }
  double tot = p->me + q->me; // pn 情形
  double dif =
      -((((j[p->i3 - 1] + j[p->i4 - 1]) / 2 - p->J + 1 - p->T) % 2 == 0) ? 1
                                                                         : -1) *
          p->me -
      ((((j[p->i3 - 1] + j[p->i4 - 1]) / 2 - p->J + 1 - q->T) % 2 == 0) ? 1
                                                                        : -1) *
          q->me;
  // p q 只有 T 不一样.
  transValpn(&r, p);
  r.me = tot;
  factorize(&r);
  print2file(fpw, r);
  if (p->i3 != p->i4) {
    transValpnex(&r, p);
    r.me = dif;
    factorize(&r);
    print2file(fpw, r);
  }
  if (p->i1 != p->i2) { // 1234 换 2143. 当 p 和 n 在不同 k 轨道时交换后输出的
    if (p->i1 != p->i3 || p->i2 != p->i4) { // km 态不同, 是不能交换出来的.
      transValpnex2(&r, p);
      r.me = -((((j[p->i1 - 1] + j[p->i2 - 1]) / 2 - p->J + 1 - p->T) % 2 == 0)
                   ? 1
                   : -1) *
                 p->me -
             ((((j[p->i1 - 1] + j[p->i2 - 1]) / 2 - p->J + 1 - q->T) % 2 == 0)
                  ? 1
                  : -1) *
                 q->me;
      factorize(&r);
      print2file(fpw, r);
    }
    if (p->i3 != p->i4) {
      if (p->i1 != p->i4 || p->i2 != p->i3) {
        transValpnex1(&r, p);
        r.me =
            ((((j[p->i1 - 1] + j[p->i2 - 1]) / 2 - p->J + 1 - p->T) % 2 == 0)
                 ? 1
                 : -1) *
                ((((j[p->i3 - 1] + j[p->i4 - 1]) / 2 - p->J + 1 - p->T) % 2 ==
                  0)
                     ? 1
                     : -1) *
                p->me +
            ((((j[q->i1 - 1] + j[q->i2 - 1]) / 2 - q->J + 1 - q->T) % 2 == 0)
                 ? 1
                 : -1) *
                ((((j[q->i3 - 1] + j[q->i4 - 1]) / 2 - q->J + 1 - q->T) % 2 ==
                  0)
                     ? 1
                     : -1) *
                q->me;
        factorize(&r);
        print2file(fpw, r);
      }
    }
  }
}

void swapVal(struct I4 *p, struct I4 *q) {
  struct I4 temp;
  transVal(&temp, p);
  transVal(p, q);
  transVal(q, &temp);
}

void copyline(FILE *fpw, FILE *fp, int n) {
  char sentence[256];
  for (int i = 0; i < n; i++) {
    fgets(sentence, 256, fp);
    fputs(sentence, fpw);
  }
}

int main() {
  char path[50], pathex[50], pathpn[50], name[50];
  int A = 11, Z = 5;
  int run = 1;
  while (run != 0) {
    printf("please input &A &Z:\n");
    scanf("%d %d", &A, &Z); // 不能写"%d %d\n" ?
    sprintf(name, "A%dZ%d", A, Z);
    sprintf(path, "%s/%s.int", PATHint, name);
    sprintf(pathex, "%s/%sex.int", PATHint, name);
    sprintf(pathpn, "%s/%spn.int", PATHint, name);
    //    sprintf(name, "%s", "ckii");
    //    sprintf(path, "%s/%s.int", PATHint, name);
    //    sprintf(pathex, "%s/%sex.int", PATHint, name);
    //    sprintf(pathpn, "%s/%spn.int", PATHint, name);

    FILE *fpr;
    if (!(fpr = fopen(path, "rt"))) {
      printf("Can't find file 1!\n");
      exit(1);
    }
    FILE *fp;
    if (!(fp = fopen(pathex, "w+"))) {
      printf("Can't find file 2!\n");
      exit(1);
    }
    copyline(fp, fpr, 2);
    struct I4 read[nlines_T];
    int nlines;
    for (nlines = 0; nlines < nlines_T && !feof(fpr); nlines++) {
      getJTME(fpr, read + nlines);
    }
    if (nlines == nlines_T) {
      printf("not enough lines.\n");
      exit(1);
    }
    for (int i = 0; i < nlines; i++) {
      for (int j = i + 1; j < nlines; j++) {
        if (ConfigEqualQ(read[i], read[j])) {
          swapVal(read + i + 1, read + j); // 把只有 T 不同的矩阵元换到相邻行
          i++;
          break;
        }
      }
    }

    for (int i = 0; i < nlines; i++) {
      print2file(fp, read[i]);
    }
    fprintf(fp, "0  0  0  0  0  0   0\n"); //最后两行 0 用来停止转换程序
    fprintf(fp, "0  0  0  0  0  0   0\n");
    fseek(fp, 0, SEEK_SET);
    FILE *fpw = fopen(pathpn, "w+");
    FILE *fpsp;
    FILE *fpspH;
    if (!(fpsp = fopen(PATHsp, "r"))) {
      printf("Can't find sp file!\n");
      exit(1);
    }
    if (!(fpspH = fopen(PATHspH, "r"))) {
      printf("Can't find sp file!\n");
      exit(1);
    }
    setMschemeID2JschemeID(fpsp, fpspH);
    int *j = (int *)malloc(n_orbtJ * sizeof(int));
    for (int i = 1; i < n_orbtJ + 1; i++) {
      j[i - 1] = JschemeJ(1, i); // 1 代表质子
      fprintf(fpw, "%d  ", j[i - 1]);
    }
    fprintf(fpw, "\n");
    moveToNextLine(fp);
    copyline(fpw, fp, 1);
    struct I4 p, q;
    getJTME(fp, &p);
    getJTME(fp, &q);
    while (!feof(fp) ||
           q.me != 0) { //一次操作两行 config 相同的, 下次读两行. 如果不相同,
      //这次只操作一行, 下次只读一行.
      //并且为了方便, 在 testt.int 的结尾填上两行 0  0  0  0  0  0   0
      if (ConfigEqualQ(p, q)) {
        convert2lines(fpw, &p, &q, j);
        getJTME(fp, &p);
        getJTME(fp, &q);
      } else { //只操作一行
        convert1line(fpw, &p, j);
        transVal(&p, &q);
        getJTME(fp, &q);
      }
    }
    if (fabs(p.me) > 1e-10) { //处理多余的行.
      convert1line(fpw, &p, j);
    }
    fclose(fp);
    fclose(fpw);
    free(j);
    printf("type 1 for continue, 0 for exit:\n");
    scanf("%d", &run);
  }
  return 0;
}
