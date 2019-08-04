//
//  main.c
//  Mtbme_T
//
//  Created by 1400011329 on 25/08/2017.
//  Copyright © 2017 Peking University. All rights reserved.
//

#include "sp_shell.h"
#include <gsl/gsl_sf_coupling.h>
#include <math.h>

int powminus1(int i) { //-1的指数
  if (i % 2 == 0) {
    return 1;
  } else
    return -1;
}

double CGcoefficient(int J1times2, int J2times2, int Jtimes2, int M1times2,
                     int M2times2, int Mtimes2) { // CG系数
  double cg = 1;
  cg *= powminus1((J1times2 - J2times2 + Mtimes2) / 2);
  cg *= sqrt(Jtimes2 + 1);
  cg *= gsl_sf_coupling_3j(J1times2, J2times2, Jtimes2, M1times2, M2times2,
                           -Mtimes2);
  return cg;
}

int *ttimes2tztimes2fromtype(int type) { //从type提取四个t和tz
  if (type <= 4 || type > 14) {
    fprintf(stderr, "ERROR in Mtbme_T's type input.\n");
    exit(1);
  }
  int *p = (int *)malloc(9 * sizeof(int));
  if (type == 5) { // lambda,p-lambda,p
    p[1] = 0;
    p[2] = 0;
    p[3] = 1;
    p[4] = -1;
    p[5] = 0;
    p[6] = 0;
    p[7] = 1;
    p[8] = -1;
  } else if (type == 6) { // lambda,n-lambda,n
    p[1] = 0;
    p[2] = 0;
    p[3] = 1;
    p[4] = 1;
    p[5] = 0;
    p[6] = 0;
    p[7] = 1;
    p[8] = 1;
  } else if (type == 7) {
    p[1] = 2;
    p[2] = -2;
    p[3] = 1;
    p[4] = -1;
    p[5] = 2;
    p[6] = -2;
    p[7] = 1;
    p[8] = -1;
  } else if (type == 8) {
    p[1] = 2;
    p[2] = -2;
    p[3] = 1;
    p[4] = 1;
    p[5] = 2;
    p[6] = -2;
    p[7] = 1;
    p[8] = 1;
  } else if (type == 9) {
    p[1] = 2;
    p[2] = 0;
    p[3] = 1;
    p[4] = -1;
    p[5] = 2;
    p[6] = 0;
    p[7] = 1;
    p[8] = -1;
  } else if (type == 10) {
    p[1] = 2;
    p[2] = 0;
    p[3] = 1;
    p[4] = -1;
    p[5] = 2;
    p[6] = -2;
    p[7] = 1;
    p[8] = 1;
  } else if (type == 11) {
    p[1] = 0;
    p[2] = 0;
    p[3] = 1;
    p[4] = -1;
    p[5] = 2;
    p[6] = 0;
    p[7] = 1;
    p[8] = -1;
  } else if (type == 12) {
    p[1] = 0;
    p[2] = 0;
    p[3] = 1;
    p[4] = -1;
    p[5] = 2;
    p[6] = -2;
    p[7] = 1;
    p[8] = 1;
  } else if (type == 13) {
    p[1] = 0;
    p[2] = 0;
    p[3] = 1;
    p[4] = 1;
    p[5] = 2;
    p[6] = 0;
    p[7] = 1;
    p[8] = 1;
  } else if (type == 14) {
    p[1] = 0;
    p[2] = 0;
    p[3] = 1;
    p[4] = 1;
    p[5] = 2;
    p[6] = 2;
    p[7] = 1;
    p[8] = -1;
  }
  return p;
}

double Vpara(double Vpar, double Delta, double Splus, double Sminus, double T,
             int J, int j2times2,
             int j4times2) { //根据J,j2times2,j4times2决定公式
  if (J == 0 && j2times2 == 1 && j4times2 == 1) {
    return Vpar + Delta / 4 - 2 * Splus - 6 * T;
  } else if (J == 1 && j2times2 == 1 && j4times2 == 1) {
    return Vpar - Delta / 12 - 2 * Splus / 3 + 2 * T + 4 * Sminus / 3;
  } else if ((J == 1 && j2times2 == 3 && j4times2 == 1) ||
             (J == 1 && j2times2 == 1 && j4times2 == 3)) {
    return sqrt(2) / 3 * (Delta - Splus + 3 * T - Sminus);
  } else if (J == 1 && j2times2 == 3 && j4times2 == 3) {
    return Vpar - 5 * Delta / 12 - Splus / 3 + T - 4 * Sminus / 3;
  } else if (J == 2 && j2times2 == 3 && j4times2 == 3) {
    return Vpar + Delta / 4 + Splus - 3 * T / 5;
  } else {
    return 0;
  }
}

// 两组参数对应 p 壳轻核和重核
double Vtjtypej2j4(int Ttimes2, int J, int type, int j2times2, int j4times2,
                   int l_or_h) { //根据type和T决定代值
  double vjt = 0;
  double Vpar, Delta, Splus, Sminus, T;
  double VparCSB = 0.024864, DeltaCSB = 0.052128;
  VparCSB = 0;
  DeltaCSB = 0;
  if (type == 5) {
    if (l_or_h == 0) {
      Vpar = -1.06 - VparCSB;
      Delta = 0.430 - DeltaCSB;
      Splus = -0.2025;
      Sminus = 0.1875;
      T = 0.030;
      vjt = Vpara(Vpar, Delta, Splus, Sminus, T, J, j2times2, j4times2);
    } else if (l_or_h == 1) {
      Vpar = -1.06 - VparCSB;
      Delta = 0.330 - DeltaCSB;
      Splus = -0.1825;
      Sminus = 0.1675;
      T = 0.0239;
      vjt = Vpara(Vpar, Delta, Splus, Sminus, T, J, j2times2, j4times2);
    }
  } else if (type == 6) {
    if (l_or_h == 0) {
      Vpar = -1.06 + VparCSB;
      Delta = 0.430 + DeltaCSB;
      Splus = -0.2025;
      Sminus = 0.1875;
      T = 0.030;
      vjt = Vpara(Vpar, Delta, Splus, Sminus, T, J, j2times2, j4times2);
    } else if (l_or_h == 1) {
      Vpar = -1.06 + VparCSB;
      Delta = 0.330 + DeltaCSB;
      Splus = -0.1825;
      Sminus = 0.1675;
      T = 0.0239;
      vjt = Vpara(Vpar, Delta, Splus, Sminus, T, J, j2times2, j4times2);
    }
  } else if (type >= 7 && type <= 10) {
    if (Ttimes2 == 1) {
      Vpar = 1.0100;
      Delta = -7.2150;
      Splus = -0.0010;
      Sminus = 0.0000;
      T = -0.3640;
      vjt = Vpara(Vpar, Delta, Splus, Sminus, T, J, j2times2, j4times2);
    } else if (Ttimes2 == 3) {
      Vpar = -1.1070;
      Delta = 2.2750;
      Splus = -0.2680;
      Sminus = 0.0000;
      T = 0.1870;
      vjt = Vpara(Vpar, Delta, Splus, Sminus, T, J, j2times2, j4times2);
    }
  } else if (type >= 11 & type <= 14) {
    Vpar = 1.4500;
    Delta = 3.0400;
    Splus = -0.0850;
    Sminus = 0.0000;
    T = 0.1570;
    vjt = Vpara(Vpar, Delta, Splus, Sminus, T, J, j2times2, j4times2);
  }
  return vjt;
}

double Mtbme_T(int type, int k1m1, int k2m2, int k3m3, int k4m4,
               int l_or_h) { //算m-scheme矩阵元, l_or_h(0:light;1:heavy)
  double mtbme_T = 0;
  int j1times2, j2times2, j3times2, j4times2, l1, l2, l3, l4, m1times2,
      m2times2, m3times2, m4times2;
  j1times2 = jfromA(typepp, k1m1);
  j2times2 = jfromA(typell, k2m2);
  j3times2 = jfromA(typepp, k3m3);
  j4times2 = jfromA(typell, k4m4);
  l1 = lfromA(typepp, k1m1);
  l2 = lfromA(typell, k2m2);
  l3 = lfromA(typepp, k3m3);
  l4 = lfromA(typell, k4m4);
  m1times2 = mfromA(typepp, k1m1);
  m2times2 = mfromA(typell, k2m2);
  m3times2 = mfromA(typepp, k3m3);
  m4times2 = mfromA(typell, k4m4);
  if (m1times2 + m2times2 != m3times2 + m4times2) {
    return 0;
  }
  if (l_or_h != 0 && l_or_h != 1) {
    fprintf(stderr, "invalid l_or_h!\n");
    exit(1);
  }
  for (int J = abs((j1times2 - j2times2) / 2); J <= (j1times2 + j2times2) / 2;
       J++) { //对J求和
    double sum_T = 0;
    int *p = ttimes2tztimes2fromtype(type);
    int t1times2 = p[3], tz1times2 = p[4], t2times2 = p[1], tz2times2 = p[2],
        t3times2 = p[7], tz3times2 = p[8], t4times2 = p[5], tz4times2 = p[6];
    int Ttimes2max =
        0; //这里赋值为0是为了在type输入不合法的时候下面对T的求和就取消了(0<1)
    if (type >= 7 && type <= 10) {
      Ttimes2max = 3;
    } else if (type >= 11 & type <= 14 || type == 5 || type == 6) {
      Ttimes2max = 1;
    }
    for (int Ttimes2 = 1; Ttimes2 <= Ttimes2max; Ttimes2 += 2) { //对T求和
      double inc = 1;
      inc *= CGcoefficient(t1times2, t2times2, Ttimes2, tz1times2, tz2times2,
                           tz1times2 + tz2times2);
      inc *= CGcoefficient(t3times2, t4times2, Ttimes2, tz3times2, tz4times2,
                           tz3times2 + tz4times2);
      inc *= Vtjtypej2j4(Ttimes2, J, type, j1times2, j3times2, l_or_h);
      sum_T += inc;
      /*if (J==1&&j2times2==3&&j4times2==3) {
          printf("inc=%lf\n",inc);
      }*/
    }
    /*if (J==1&&j2times2==3&&j4times2==3) {
        printf("aaaa%lf\n",sum_T);
    }*/
    sum_T *= CGcoefficient(j1times2, j2times2, J * 2, m1times2, m2times2,
                           m1times2 + m2times2) *
             CGcoefficient(j3times2, j4times2, J * 2, m3times2, m4times2,
                           m3times2 + m4times2);
    mtbme_T += sum_T;
  }
  //  if (fabs(mtbme_T) > 1e-3) {
  //    printf("%d %d %d %d %d %d %d %d %d %d %d %d\n", j1times2, j2times2,
  //           j3times2, j4times2, l1, l2, l3, l4, m1times2, m2times2, m3times2,
  //           m4times2);
  //  }
  return mtbme_T;
}
// int main(int argc, const char * argv[]) {
//     int type, k1m1, k2m2, k3m3, k4m4;
//     //int Ttimes2=3, J=1, type=8, j2times2=3, j4times2=3;
//     //printf("%f\n",Vtjtypej2j4(Ttimes2, J, type, j2times2, j4times2));
//     printf("Input:(type,k1m1,k2m2,k3m3,k4m4)\n");
//     scanf("%d,%d,%d,%d,%d", &type, &k1m1, &k2m2, &k3m3, &k4m4);
//     printf("%lf\n", Mtbme_T(type, k1m1, k2m2, k3m3, k4m4));
//     return 0;
// }
