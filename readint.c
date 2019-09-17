#include "sp_shell.h"

int **getdatac(FILE *fp) {
  fseek(fp, 0, SEEK_SET);
  int **p, j;
  p = (int **)malloc(num_lines * sizeof(int *)); //存TBME 的 type
                                                 //和四个单粒子态. 1~3是质子的
                                                 // s1/2, p3/2, p1/2, 4~6是中子
  moveToNextLine(fp);
  moveToNextLine(fp); //跳过开头两行没用的.
  for (int i = 0; i < num_lines; i++) {
    p[i] = (int *)malloc(5 * sizeof(int));
    for (j = 0; j < 5; j++) {
      //     fscanf(fp, "  %c", &p[i][j]);
      fscanf(fp, "%d ", &p[i][j]); //字符转成数字
    }
    moveToNextLine(fp);
  }
  return p;
}

double *getdatad(FILE *fp) {
  fseek(fp, 0, SEEK_SET);
  double *p;
  p = (double *)malloc(num_lines * sizeof(double)); //存 TBME 的值
  moveToNextLine(fp);
  moveToNextLine(fp); //跳过开头两行没用的.
  for (int i = 0; i < num_lines; i++) {
    for (int j = 0; j < 6;
         j++) { //若 pn 表象的 .int 文件有一列 T' 就改为 6, 多读一个 int
      int f;
      fscanf(fp, "%c ", &f);
    }
    fscanf(fp, "%lf", &p[i]);
    moveToNextLine(fp);
  }
  return p;
}

//以下计算超子-核子相互作用*
int **Hgetdatac(FILE *fp) {
  fseek(fp, 0, SEEK_SET);
  int **p, i, j;
  p = (int **)malloc(numlinesH * sizeof(int *));
  for (i = 0; i < numlinesH; i++) {
    p[i] = (int *)malloc(6 * sizeof(int));
    for (j = 0; j < 6; j += 2) {
      fscanf(fp, "%d  ", &p[i][j]);
      p[i][j + 1] = 1; //在偶数位补上 1, 代表 Λ 只处在 s1/2.
    }
    moveToNextLine(fp);
  }
  return p;
}

double *Hgetdatad(FILE *fp) {
  fseek(fp, 0, SEEK_SET);
  int i;
  double *p;
  p = (double *)malloc(numlinesH * sizeof(double));
  for (i = 0; i < numlinesH; i++) {
    for (int j = 0; j < 3; j++) {
      int f;
      fscanf(fp, "%c  ", &f);
    }
    fscanf(fp, "%lf", &p[i]);
    moveToNextLine(fp);
  }
  return p;
}

double hw(int A, int Z) {
  //  switch (Z) {
  //  case 2: {
  //    switch (A) {
  //    case 6:
  //      return 14.2;
  //    case 7:
  //      return 14.7;
  //    case 8:
  //      return 14.3;
  //    default:
  //      break;
  //    }
  //  }
  //  case 3: {
  //    switch (A) {
  //    case 6:
  //      return 14.9;
  //    case 7:
  //      return 15.1;
  //    case 8:
  //      return 15.2;
  //    case 9:
  //      return 13.0;
  //    default:
  //      break;
  //    }
  //  }
  //  case 4: {
  //    switch (A) {
  //    case 8:
  //      return 17.0;
  //    case 9:
  //      return 13.8;
  //    case 10:
  //      return 16.6;
  //    default:
  //      break;
  //    }
  //  }
  //  case 5: {
  //    switch (A) {
  //    case 9:
  //      return 14.8;
  //    case 10:
  //      return 16.6;
  //    case 11:
  //      return 14.9;
  //    case 12:
  //      return 14.4;
  //    case 13:
  //      return 12.6;
  //    default:
  //      break;
  //    }
  //  }
  //  case 6: {
  //    switch (A) {
  //    case 11:
  //      return 15.1;
  //    case 12:
  //      return 14.4;
  //    case 13:
  //      return 14.1;
  //    case 14:
  //      return 14.0;
  //    default:
  //      break;
  //    }
  //  }
  //  case 7: {
  //    switch (A) {
  //    case 14:
  //      return 14.1;
  //    case 15:
  //      return 13.1;
  //    default:
  //      break;
  //    }
  //  }
  //  case 8: {
  //    switch (A) {
  //    case 15:
  //      return 13.7;
  //    default:
  //      break;
  //    }
  //  }
  //  }
  //  printf("hw for A Z doesn't exist! return empirical\n");
  return 45 * pow(A, -1 / 3.) - 25 * pow(A, -2 / 3.);
}
