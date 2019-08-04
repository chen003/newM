//
//  main.c
//  nuclear_A2B1
//
//  Created by xingyan chen on 09/03/2017.
//  Copyright © 2017 xingyan chen. All rights reserved.
//
#include "eigen_calc.h"
#include "sp_shell.h"
#include <complex.h>
#include <math.h>
#include <stdio.h>
#define dim_tot power(n_orbtM, 3)
#define dim_r 21
//#define dim 57
// #define PATH "../M-scheme/int/test2pn.int"
#define PATH "../newM/int/p-shellpn.int"
#define PATHh "../newM/int/p-lam.int"

// file name A3Z1pn.int
int main() {
  //  freopen("out.txt", "w", stdout);

  FILE *fp;
  if (!(fp = fopen(PATH, "rt"))) {
    printf("Can't find file 1!\n");
    exit(1);
  }
  FILE *fpH;
  if (!(fpH = fopen(PATHh, "rt"))) {
    printf("Can't find file 2!\n");
    exit(1);
  }

  //**************单独使用 jfromA 等. ***********
  //  FILE *fpsp;
  //  FILE *fpspH; // j,m是两倍
  //  if (!(fpsp = fopen(PATHsp, "r"))) {
  //    printf("Can't find sp file!\n");
  //    exit(1);
  //  }
  //  if (!(fpspH = fopen(PATHspH, "r"))) {
  //    printf("Can't find sp file!\n");
  //    exit(1);
  //  }
  //  setMschemeID2JschemeID(fpsp, fpspH);
  //*******************************************

  makeTBME(fp, fpH, 0); // 这一步后才能使用 mfromA 函数
  //**************查看TBME*********************
  //    printf("input:\n");
  //    int t,a,b,c,d;
  //    scanf("%d %d %d %d %d",&t, &a,&b,&c,&d);
  //    printf("%f\n",TBME(t,a,b,c,d));
  //  double me;
  //  for (int type = 1; type <= 14; type++) {
  //    for (int a = 1; a <= 6; a++) {
  //      for (int b = 1; b <= 6; b++) {
  //        for (int c = 1; c <= 6; c++) {
  //          for (int d = 1; d <= 6; d++) {
  //            me = TBME(type, a, b, c, d);
  //            if (fabs(me) > 1e-6) {
  //              printf("%d\t%d\t%d\t%d\t%d\t%f\n", type, a, b, c, d, me);
  //            }
  //          }
  //        }
  //      }
  //    }
  //  }
  //*******************************************

  //  fclose(stdout);
  return 0;
}
