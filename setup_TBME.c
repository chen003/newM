//
//  main.c
//  test
//
//  Created by 1400011329 on 29/03/2017.
//  Copyright © 2017 Peking University. All rights reserved.
//
#include "sp_shell.h"
#include <gsl/gsl_sf_coupling.h>

// static double OTBMEsav[5184],HTBMEsav[1296]; 24576
static double TBMEsav[18144];

//用 CG 系数和对称关系 从 J-scheme 得到 m-scheme TBME. pn 情形对称关系不成立,
//但由于两边同时交换, 对结果不影响.
double OTBMEm(int type, int A, int B, int a, int b, int **od1, double *od2,
              int lines) //读取函数并返回矩阵元
{
  int j1ld, j2ld, j1rd, j2rd, l1l, l2l, l1r, l2r; // l 左矢, r 右矢
  int m1ld, m2ld, m1rd, m2rd;
  int type1 = typeAfromAB(type, 1);
  int type2 = typeAfromAB(type, 2);
  j1ld = jfromA(type1, A);
  j2ld = jfromA(type2, B);
  j1rd = jfromA(type1, a);
  j2rd = jfromA(type2, b);
  l1l = lfromA(type1, A);
  l2l = lfromA(type2, B);
  l1r = lfromA(type1, a);
  l2r = lfromA(type2, b);
  m1ld = mfromA(type1, A);
  m2ld = mfromA(type2, B);
  m1rd = mfromA(type1, a);
  m2rd = mfromA(type2, b);
  int J, Jmax;
  double
      tbme = 0,
      ele = 0,
      inc;            // tbme是输出的值，ele是J-scheme出的矩阵元(加上对称关系得到的符号)，inc是连加号后的每一项（ele*两个CG系数）
  int n1, n2, n3, n4; // J-scheme的前四个编号，1~3质子s1/2,p3/2,p1/2，4~6中子
  int i;
  Jmax = (j1rd + j2rd) / 2;
  // J-scheme矩阵元转到m-scheme。 对J和M求和
  for (J = abs((j1rd - j2rd) / 2); J <= Jmax; J++) { //对J求和
    n1 = numfromA(A, type1);                         //读取
    n2 = numfromA(B, type2);
    n3 = numfromA(a, type1);
    n4 = numfromA(b, type2);
    for (i = 0; i < lines; i++) { //按各种交换方式查找, 利用对称条件得到系数.
      if (od1[i][0] == n1 && od1[i][1] == n2 && od1[i][2] == n3 &&
          od1[i][3] == n4 && od1[i][4] == J) {
        ele = od2[i];
        break;
      } else if (od1[i][0] == n1 && od1[i][1] == n2 && od1[i][2] == n4 &&
                 od1[i][3] == n3 && od1[i][4] == J && type < 3) {
        ele = -od2[i] * ((((j1rd + j2rd) / 2 - J) % 2) == 0 ? 1 : -1);
        break;
      } else if (od1[i][0] == n2 && od1[i][1] == n1 && od1[i][2] == n4 &&
                 od1[i][3] == n3 && od1[i][4] == J) {
        ele = od2[i] * ((((j1ld + j2ld + j1rd + j2rd) / 2 - 2 * J) % 2 == 0)
                            ? 1
                            : -1); //左右都交换
        break;
      } else if (od1[i][0] == n2 && od1[i][1] == n1 && od1[i][2] == n3 &&
                 od1[i][3] == n4 && od1[i][4] == J && type < 3) {
        ele = -od2[i] * ((((j1ld + j2ld) / 2 - J) % 2 == 0 ? 1 : -1));
        break;
      } else if (od1[i][0] == n3 && od1[i][1] == n4 && od1[i][2] == n1 &&
                 od1[i][3] == n2 && od1[i][4] == J) {
        ele = od2[i];
        break;
      } else if (od1[i][0] == n4 && od1[i][1] == n3 && od1[i][2] == n1 &&
                 od1[i][3] == n2 && od1[i][4] == J && type < 3) {
        ele = -od2[i] * ((((j1rd + j2rd) / 2 - J) % 2 == 0 ? 1 : -1));
        break;
      } else if (od1[i][0] == n4 && od1[i][1] == n3 && od1[i][2] == n2 &&
                 od1[i][3] == n1 && od1[i][4] == J) {
        ele = od2[i] *
              ((((j1ld + j2ld + j1rd + j2rd) / 2 - 2 * J) % 2 == 0) ? 1 : -1);
        break;
      } else if (od1[i][0] == n3 && od1[i][1] == n4 && od1[i][2] == n2 &&
                 od1[i][3] == n1 && od1[i][4] == J && type < 3) {
        ele = -od2[i] * ((((j1ld + j2ld) / 2 - J) % 2 == 0 ? 1 : -1));
        break;
      }
    }
    if (i == lines) {
      ele = 0;
    }
    int M2 = (m1ld + m2ld);
    inc = (((j1ld - j2ld + j1rd - j2rd) / 2 + M2) % 2 == 0 ? 1 : -1) *
          (2 * (double)J + 1) *
          gsl_sf_coupling_3j(j1ld, j2ld, (2 * J), m1ld, m2ld, -(M2)) *
          gsl_sf_coupling_3j(j1rd, j2rd, (2 * J), m1rd, m2rd, -(M2)) *
          ele; // 这里 M2 放在第二个 cg 系数里充当了δ函数的作用
    tbme += inc;
  }
  if (type < 3) { //全同粒子时 J-scheme TBME 会带上归一化系数
    double tbme1 = sqrt((1 + (n1 == n2)) * (1 + (n3 == n4))) * tbme;
    return tbme1;
  } else {
    return tbme;
  }
}

double HTBMEm(int A, int B, int a, int b, int **od1, double *od2,
              int linesH) //读取函数并返回矩阵元
{
  return OTBMEm(5, A, B, a, b, od1, od2, linesH);
}

void makeTBME(FILE *fp, FILE *fpH,
              int l_or_h) { // fp是核子核子相互作用, 跟 A 有关.

  FILE *fpsp;
  FILE *fpspH; // j,m是两倍
  if (!(fpsp = fopen(PATHsp, "r"))) {
    printf("Can't find sp file!\n");
    exit(1);
  }
  if (!(fpspH = fopen(PATHspH, "r"))) {
    printf("Can't find sp file!\n");
    exit(1);
  }
  setMschemeID2JschemeID(fpsp, fpspH);
  int **od1 = getdatac(fp);   //读前面六列
  double *od2 = getdatad(fp); //读最后一列
                              //  int **od3 = Hgetdatac(fpH);  //读前面六列
  //  double *od4 = Hgetdatad(fpH); //读最后一列
  for (int type = 1; type < n_type + 1 + 8; type++) {
    for (int a = 1; a < n_orbtM + 1; a++) {
      for (int b = 1; b < n_orbtM + 1; b++) {
        for (int c = 1; c < n_orbtM + 1; c++) {
          for (int d = 1; d < n_orbtM + 1; d++) {
            //存到静态数组 TBMEsav 中, 访问时用两题相互作用的 type 和左右矢的
            // m-scheme 单粒子态. -1555因为 type 和 abcd 都是从1开始的.
            // 1555 = 1296 + 216 + 36 + 6 + 1.
            if (type < 3 || type == 4) {
              TBMEsav[power(n_orbtM, 4) * (type - 1) +
                      power(n_orbtM, 3) * (a - 1) +
                      power(n_orbtM, 2) * (b - 1) + n_orbtM * (c - 1) + d - 1] =
                  OTBMEm(type, a, b, c, d, od1, od2, num_lines);
            } else if (type == 3) {
              TBMEsav[power(n_orbtM, 4) * (type - 1) +
                      power(n_orbtM, 3) * (a - 1) +
                      power(n_orbtM, 2) * (b - 1) + n_orbtM * (c - 1) + d - 1] =
                  0;
            } else {
              TBMEsav[power(n_orbtM, 4) * (type - 1) +
                      power(n_orbtM, 3) * (a - 1) +
                      power(n_orbtM, 2) * (b - 1) + n_orbtM * (c - 1) + d - 1] =
                  Mtbme_T(type, a, b, c, d, l_or_h);
            }
          }
        }
      }
    }
  }
  fclose(fpsp);
}

double TBME(int type, int A, int B, int a,
            int b) //从static一维数组中提取矩阵元
{
  double tbme =
      TBMEsav[power(n_orbtM, 4) * (type - 1) + power(n_orbtM, 3) * (A - 1) +
              power(n_orbtM, 2) * (B - 1) + n_orbtM * (a - 1) + b - 1];
  if (type == 5 && B > 2 && b > 2) {
    tbme = 0;
  } else if (type == 6 && B > 2 && b > 2) {
    tbme = 0;
  }
  return tbme;
}

//输出结果, 根据对称性和直觉手动检查
void writeTBME(void) {
  FILE *fp;
  if (!(fp = fopen("result.txt", "wt"))) {
    printf("Can't write file!\n");
    exit(1);
  }
  for (int type = 1; type < n_type + 1; type++) {
    for (int a = 1; a < n_orbtM + 1; a++) {
      for (int b = 1; b < n_orbtM + 1; b++) {
        for (int c = 1; c < n_orbtM + 1; c++) {
          for (int d = 1; d < n_orbtM + 1; d++) {
            if (TBME(type, a, b, c, d) != 0) {
              fprintf(fp, "%d  %d  %d  %d  %d  %lf\n", type, a, b, c, d,
                      TBME(type, a, b, c, d));
            }
          }
        }
      }
    }
  }
  fclose(fp);
}
