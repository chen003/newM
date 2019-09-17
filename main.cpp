#include "sp_shell.h"
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sf_coupling.h>
// #include <gsl/gsl_matrix.h>
#include "eigen_calc.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#define dim_occup 500

const char *isotope[7] = {"He", "Li", "Be", "B", "C", "N", "O"};

#define PATH                                                                   \
  "../newM/int/ckipn.int" // 手动改成其他 pn.int 文件. 比如从 nushell
                          // 拷贝过来的
// cki(6-16) 和 ckii(8-16)
#define PATHh "../newM/int/p-lam.int" // 现在不需要了.
#define PATHgogny                                                              \
  "../newM/gogny" // 有两组, 一组是自己算的 spe, 另一组是实验 spe

#define PATHPLOT "./plot"
#define PATHPLOTEXP "../newM/plotfiles"

/*
 * 程序考虑了: (后两者默认关闭, 因为精度尚未达到要求)
 * 1. p 壳轻重核的两套参数.
 * 2. Gal 2015 PLB 中提到的 CSB
 * 3. Σ+-0 的质量差
 * 可调参数:
 * 1. Σ 单粒子能 para
 * 2. 若使用 gogny 力, 可在 int 文件中调 hw (不再本程序内)
 * 可计算:
 * 1. 能谱
 * 2. Λ 和 ΛΛ 结合能
 */

// 画图用
#define NUM_LEVELS 6 //  画几条
#define NUM_POINTS 5
#define NUM_COMMANDS 2
const double width =
    0.04; // 根据图片打印出来的实际尺寸和12号字体算出来的相对间隔
// ********************************************************************************1
int fromdoubletoint(double d) {
  if (d <= -0.1) {
    return -1;
  }
  for (int i = 0; i < 32767; i++) {
    if (d - i < 0.1 && i - d < 0.1) {
      return i;
    }
  }
  return -1;
}
// ********************************************************************************1
void CreateConfig(struct config *cfg, struct config *cfgJplus, int n_proton,
                  int n_neutron, int n_lam, int *dim, int *dimJplus, int M2) {
  int n1 = 0, n2 = 0;
  int cfgmaxp = power(2, n_orbtM) - power(2, n_orbtM - n_proton);
  int cfgmaxn = power(2, n_orbtM) - power(2, n_orbtM - n_neutron);
  int cfgmaxl =
      power(2, n_orbtMH) - power(2, n_orbtMH - n_lam); // 三种粒子的最大位数
  for (int p_cfg = power(2, n_proton) - 1; p_cfg <= cfgmaxp;
       p_cfg = nextbitpermut(p_cfg)) {
    for (int n_cfg = power(2, n_neutron) - 1; n_cfg <= cfgmaxn;
         n_cfg = nextbitpermut(n_cfg)) {
      for (int l_cfg = power(2, n_lam) - 1; l_cfg <= cfgmaxl;
           l_cfg = nextbitpermut(l_cfg)) { //  用算法按字典序得到下一组态
        if (mfromPermut(typepp, p_cfg) + mfromPermut(typenn, n_cfg) +
                mfromPermut(typell, l_cfg) ==
            M2) {
          cfg[n1].p = p_cfg;
          cfg[n1].n = n_cfg;
          cfg[n1].lam = l_cfg;
          n1++;
        } else if (mfromPermut(typepp, p_cfg) + mfromPermut(typenn, n_cfg) +
                       mfromPermut(typell, l_cfg) ==
                   M2 + 2) { //  升算符作用后的态由这些基组成
          cfgJplus[n2].p = p_cfg;
          cfgJplus[n2].n = n_cfg;
          cfgJplus[n2].lam = l_cfg;
          n2++;
        }
      }
    }
  }
  *dim = n1;
  *dimJplus = n2;
  if (*dim > dim_occup) {
    fprintf(stderr, "out of dimension.\n"); //  预先分配的内存不足
    exit(1);
  }
}

//  非对角元
//  Λ-Σ0 项
void CreateHamiltonianOffDiagDif0(double **h, int beginright, int beginleft,
                                  int dimleft, int dimright,
                                  struct config *cfgleft,
                                  struct config *cfgright, int typeNL,
                                  int typePL) {
  double me = 0;
  for (int i = 0; i < dimleft; i++) {
    for (int j = 0; j < dimright; j++) {
      me = 0;
      int difp =
          difcount(cfgleft[i].p, cfgright[j].p); // 根据缩并情况讨论相互作用类型
      int difn = difcount(cfgleft[i].n, cfgright[j].n);
      //  五种相互作用. (除ΛΛ以外)
      me += (difp == 0 ? TBOME2(typeNL, cfgleft[i].n, cfgleft[i].lam,
                                cfgright[j].n, cfgright[j].lam)
                       : 0);
      me += (difn == 0 ? TBOME2(typePL, cfgleft[i].p, cfgleft[i].lam,
                                cfgright[j].p, cfgright[j].lam)
                       : 0);

      h[i + beginleft][j + beginright] = h[j + beginright][i + beginleft] = me;
    }
  }
}

int signDifn(int p,
             int m) { //  返回反对易关系得到的+-1, 方法是求 p 在第 m
                      //  位之后有多少占据粒子, 作用在质子部分
  int sign = 1;
  p = p & ~(~0 << (m - 1)); //  屏蔽 m 位之后的
  if (__builtin_popcount(p) % 2 != 0) {
    sign = -1;
  }
  return sign;
}

int signDifp(int p, int m) { //  返回反对易关系得到的+-1, 类比 beta
  //  衰变的算符(同位旋升降算符). 求 p 在第 m
  //  位之前有多少占据粒子, 作用在中子部分, 两个函数联合使用
  int sign = 1;
  p = p & (~0 << m); //  屏蔽 m 位之前的
  if (__builtin_popcount(p) % 2 != 0) {
    sign = -1;
  }
  return sign;
}

void CreateHamiltonianOffDiagDif1(double **h, int beginright, int beginleft,
                                  int dimleft, int dimright,
                                  struct config *cfgleft,
                                  struct config *cfgright, int type) {
  double me = 0;
  for (int i = 0; i < dimleft; i++) {
    for (int j = 0; j < dimright; j++) {
      me = 0;
      int difp =
          difcount(cfgleft[i].p, cfgright[j].p); // 根据缩并情况讨论相互作用类型
      int difn = difcount(cfgleft[i].n, cfgright[j].n);
      int sign, m1, m3;
      if (difp == 1 && difn == 1) {
        if (type == typeLampSig1n) {
          m1 = leaddifbit(cfgleft[i].p, cfgright[j].p);
          m3 = leaddifbit(cfgright[j].n, cfgleft[i].n);
          sign = signDifp(cfgright[j].p, m1) * signDifn(cfgleft[i].n, m3);
          me = sign * TBME(type, m1, cfgleft[i].lam, m3, cfgright[j].lam);
        } //  使用 leaddifbit 目的是当粒子数不对时会出负号.
        else if (type == typeLamnSig3p) {
          m1 = leaddifbit(cfgleft[i].n, cfgright[j].n);
          m3 = leaddifbit(cfgright[j].p, cfgleft[i].p);
          sign = signDifn(cfgright[j].n, m1) * signDifp(cfgleft[i].p, m3);
          me = sign * TBME(type, m1, cfgleft[i].lam, m3, cfgright[j].lam);
        } else if (type == typeSig2pSig1n) { //  这里要左右换个位置, 因为
                                             //  sig1,2,3 有顺序
          m1 = leaddifbit(cfgright[j].p, cfgleft[i].p);
          m3 = leaddifbit(cfgleft[i].n, cfgright[j].n);
          sign = signDifp(cfgleft[i].p, m1) * signDifn(cfgright[j].n, m3);
          me = sign * TBME(type, m1, cfgright[j].lam, m3, cfgleft[i].lam);
        }
      }
      h[i + beginleft][j + beginright] = h[j + beginright][i + beginleft] = me;
    }
  }
}

double massdif(int type, double para) {
  double mass = 75;
  //   switch (type) {
  //   case -1: //  用 typeNL - typePL 来区分 Σ+-0
  //     mass = 81.766;
  //     break;
  //   case 0:
  //     mass = 76.959;
  //     break;
  //   case 1:
  //     mass = 73.687;
  //     break;
  //   default:
  //     fprintf(stderr, "mass type error.\n");
  //     exit(1);
  //     break;
  //   }
  return mass + para; //  para 用来整体平移, 相当于调 Σ 的单粒子能
}

void CreateHamiltonianDiag(double **h, int begin, int dim, struct config *cfg,
                           int typeNL, int typePL) {
  for (int i = 0; i < dim; i++) {
    for (int j = i; j < dim; j++) {
      double me =
          (i == j ? singlePermutE(cfg[i].p, typepp) +
                        singlePermutE(cfg[i].n, typenn) +
                        singlePermutE(cfg[i].lam, typell) +
                        (typeNL == typenl ? 0 : massdif(typeNL - typePL, -10))
                  : 0);                        //  单粒子能.
      int difp = difcount(cfg[i].p, cfg[j].p); // 根据缩并情况讨论相互作用类型
      int difn = difcount(cfg[i].n, cfg[j].n);
      int difl = difcount(cfg[i].lam, cfg[j].lam);
      //  五种相互作用. (除ΛΛ以外)
      me +=
          (difp == 0 ? TBOME(typeNL, cfg[i].n, cfg[i].lam, cfg[j].n, cfg[j].lam)
                     : 0);
      me +=
          (difn == 0 ? TBOME(typePL, cfg[i].p, cfg[i].lam, cfg[j].p, cfg[j].lam)
                     : 0);
      me += (difl == 0 ? TBOME(typepn, cfg[i].p, cfg[i].n, cfg[j].p, cfg[j].n)
                       : 0);
      me += ((difp + difl) == 0 ? TBOMEnn(typenn, cfg[i].n, cfg[j].n) : 0);
      me += ((difn + difl) == 0 ? TBOMEnn(typepp, cfg[i].p, cfg[j].p) : 0);
      h[i + begin][j + begin] = h[j + begin][i + begin] = me;
    }
  }
}

void CreateJsqr(double *eigen_v, double *eigen_v2, struct config *cfgJplus,
                struct config *cfg, int dim, int dimJplus) {
  double *eigen_vJplus =
      (double *)malloc(dimJplus * sizeof(double)); // 升算符作用后
  for (int i = 0; i < dimJplus; i++) {
    eigen_vJplus[i] = 0;
  }
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dimJplus; j++) {
      int difp = difcount(cfgJplus[j].p, cfg[i].p);
      int difn = difcount(cfgJplus[j].n, cfg[i].n);
      int difl = difcount(cfgJplus[j].lam, cfg[i].lam);
      if ((difp + difn + difl) != 2) { // 粒子数在生成时就保证正确了
        continue;
      } else if (difp == 2) {
        eigen_vJplus[j] += vJplus(typepp, cfgJplus[j].p, cfg[i].p, eigen_v[i]);
      } else if (difn == 2) {
        eigen_vJplus[j] += vJplus(typenn, cfgJplus[j].n, cfg[i].n, eigen_v[i]);
      } else if (difl == 2) {
        eigen_vJplus[j] +=
            vJplus(typell, cfgJplus[j].lam, cfg[i].lam, eigen_v[i]);
      }
    }
  }
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dimJplus; j++) {
      int difp = difcount(cfgJplus[j].p, cfg[i].p);
      int difn = difcount(cfgJplus[j].n, cfg[i].n);
      int difl = difcount(cfgJplus[j].lam, cfg[i].lam);
      if ((difp + difn + difl) != 2) { // 粒子数在生成时就保证正确了
        continue;
      } else if (difp == 2) {
        eigen_v2[i] +=
            vJminus(typepp, cfgJplus[j].p, cfg[i].p, eigen_vJplus[j]);
      } else if (difn == 2) {
        eigen_v2[i] +=
            vJminus(typenn, cfgJplus[j].n, cfg[i].n, eigen_vJplus[j]);
      } else if (difl == 2) {
        eigen_v2[i] +=
            vJminus(typell, cfgJplus[j].lam, cfg[i].lam, eigen_vJplus[j]);
      }
    }
  }
  free(eigen_vJplus);
}

void makeBasis(int A, int Z, int n_lam, int *dim, struct config *cfg) {
  int n_proton;
  int n_neutron;
  n_proton = Z - 2;
  n_neutron = A - Z - 2 - n_lam;
  // makeTBME(fp, fpH, l_or_h); //  这一步后才能使用 mfromA 函数

  int M2 =
      ((n_proton + n_neutron + n_lam) % 2 == 0 ? 0 : -1); //  最低角动量 z 分量

  int n1 = 0;
  int cfgmaxp = power(2, n_orbtM) - power(2, n_orbtM - n_proton);
  int cfgmaxn = power(2, n_orbtM) - power(2, n_orbtM - n_neutron);
  int cfgmaxl =
      power(2, n_orbtMH) - power(2, n_orbtMH - n_lam); // 三种粒子的最大位数
  for (int p_cfg = power(2, n_proton) - 1; p_cfg <= cfgmaxp;
       p_cfg = nextbitpermut(p_cfg)) {
    for (int n_cfg = power(2, n_neutron) - 1; n_cfg <= cfgmaxn;
         n_cfg = nextbitpermut(n_cfg)) {
      for (int l_cfg = power(2, n_lam) - 1; l_cfg <= cfgmaxl;
           l_cfg = nextbitpermut(l_cfg)) { //  用算法按字典序得到下一组态
        if (mfromPermut(typepp, p_cfg) + mfromPermut(typenn, n_cfg) +
                mfromPermut(typell, l_cfg) ==
            M2) {
          cfg[n1].p = p_cfg;
          cfg[n1].n = n_cfg;
          cfg[n1].lam = l_cfg;
          n1++;
        }
      }
    }
  }
  *dim = n1;
}

//  操作的能级类型
struct level {
  double e;
  double j;
  bool prty;
  double *v;
};
//  计算单个核
int shell(int A, int Z, int n_lam, struct level *level, int interaction) {
  int n_proton;
  int n_neutron;
  n_proton = Z - 2;
  n_neutron = A - Z - 2 - n_lam;
  char path[50], pathcreate[50], pathh[50] = PATHh;
  if (interaction == 0) {
    sprintf(path, "%s/A%dZ%dpn.int", PATHgogny, A - n_lam, Z);
  } else if (interaction == 1) {
    sprintf(path, PATH);
  } else {
    printf("interaction error!\n");
    return 0;
  }
  //             sprintf(pathh, "../newMBack/int/p-lam%d_45_MBE.int", A);
  if (n_lam == 0) {
    sprintf(pathcreate, "result/%d%s.txt", A, isotope[Z - 2]);
  } else {
    sprintf(pathcreate, "result/%d%sΛ.txt", A, isotope[Z - 2]);
  }

  FILE *fp;
  if (!(fp = fopen(path, "rt"))) {
    fprintf(stderr, "Can't find file 1 for A%dZ%d!\n", A, Z);
    return 0;
  }
  FILE *fpcreate;
  if (!(fpcreate = fopen(pathcreate, "w"))) {
    fprintf(stderr, "Can't find file create!\n");
    exit(1);
  }
  FILE *fpH;
  if (!(fpH = fopen(pathh, "rt"))) {
    fprintf(stderr, "Can't find file 2!\n");
    exit(1);
  }
  int l_or_h = 1;
  if (A <= 10) {
    l_or_h = 0;
  }
  makeTBME(fp, fpH, l_or_h); //  这一步后才能使用 mfromA 函数

  int M2 =
      ((n_proton + n_neutron + n_lam) % 2 == 0 ? 0 : -1); //  最低角动量 z 分量

  struct config *cfgl =
      (struct config *)malloc(dim_occup * sizeof(struct config));
  struct config *cfglJplus = (struct config *)malloc(
      dim_occup * sizeof(struct config)); // 存升算符作用后的结果.
  int diml, dimlJplus, dimsig1, dimsig1Jplus, dimsig2, dimsig2Jplus, dimsig3,
      dimsig3Jplus; //  sig1,2,3 对应 Σ+0-
  CreateConfig(cfgl, cfglJplus, n_proton, n_neutron, n_lam, &diml, &dimlJplus,
               M2);
  struct config cfgsig2[diml], cfgsig2Jplus[dimlJplus];
  CreateConfig(cfgsig2, cfgsig2Jplus, n_proton, n_neutron, n_lam, &dimsig2,
               &dimsig2Jplus, M2);

  struct config *cfgsig1 =
      (struct config *)malloc(dim_occup * sizeof(struct config));
  struct config *cfgsig1Jplus = (struct config *)malloc(
      dim_occup * sizeof(struct config)); // 存升算符作用后的结果.
  struct config *cfgsig3 =
      (struct config *)malloc(dim_occup * sizeof(struct config));
  struct config *cfgsig3Jplus = (struct config *)malloc(
      dim_occup * sizeof(struct config)); // 存升算符作用后的结果.

  CreateConfig(cfgsig1, cfgsig1Jplus, n_proton - 1, n_neutron + 1, n_lam,
               &dimsig1, &dimsig1Jplus, M2);
  CreateConfig(cfgsig3, cfgsig3Jplus, n_proton + 1, n_neutron - 1, n_lam,
               &dimsig3, &dimsig3Jplus, M2);

  const int dim = diml + dimsig1 + dimsig2 + dimsig3;
  //   const int dimJplus = dimlJplus + dimsig1Jplus + dimsig2Jplus +
  //   dimsig3Jplus;

  double **h = (double **)malloc(sizeof(double *) * dim);
  {
    double *p = (double *)malloc(sizeof(double) * dim * dim);
    for (int i = 0; i < dim; i++) {
      h[i] = p + i * dim;
    }
  }

  fp = fopen(path, "rt"); //  singleE 需要读取文件 pn.int
  set_singleE(fp);

  for (int i = 0; i < dim; i++) {
    for (int j = i; j < dim; j++) {
      h[i][j] = h[j][i] = 0;
    }
  }

  // 计算哈密顿矩阵, 对于两个核子直接调用 TBME 即可
  //  利用对称性只需要处理上三角.
  CreateHamiltonianDiag(h, 0, diml, cfgl, typenl, typepl);
  CreateHamiltonianDiag(h, diml, dimsig1, cfgsig1, typeSig1n, typeSig1p);
  CreateHamiltonianDiag(h, diml + dimsig1, dimsig2, cfgsig2, typeSig2p,
                        typeSig2p);
  CreateHamiltonianDiag(h, diml + dimsig1 + dimsig2, dimsig3, cfgsig3,
                        typeSig1p, typeSig1n);

  CreateHamiltonianOffDiagDif1(h, diml, 0, diml, dimsig1, cfgl, cfgsig1,
                               typeLampSig1n);
  CreateHamiltonianOffDiagDif1(h, diml + dimsig1 + dimsig2, 0, diml, dimsig3,
                               cfgl, cfgsig3, typeLamnSig3p);
  CreateHamiltonianOffDiagDif1(h, diml + dimsig1, diml, dimsig1, dimsig2,
                               cfgsig1, cfgsig2,
                               typeSig2pSig1n); //  已经在内部调整了位置
  CreateHamiltonianOffDiagDif1(h, diml + dimsig1 + dimsig2, diml + dimsig1,
                               dimsig2, dimsig3, cfgsig2, cfgsig3,
                               typeSig2pSig1n);

  CreateHamiltonianOffDiagDif0(h, diml + dimsig1, 0, diml, dimsig2, cfgl,
                               cfgsig2, typeLamnSig2n, typeLampSig2p);
  /*
    char hmat[50];
    sprintf(hmat, "../A%dZ%dΛ%dhmat.txt", A, Z, n_lam);
    FILE *fpmat = fopen(hmat, "w+");
    //  打印哈密顿矩阵
    for (int i = 0; i < dim; i++) {
      for (int j = 0; j < dim; j++) {
        fprintf(fpmat, "%f\t", h[i][j]);
      }
      fprintf(fpmat, "\n");
    }

    // ******************* 打印哈密顿矩阵对应的态, 检查
    char s1[n_orbtM + 1], s2[n_orbtM + 1], s3[n_orbtMH + 1], t1[n_orbtM + 1],
        t2[n_orbtM + 1], t3[n_orbtMH + 1];
    for (int i = 0; i < dimsig1; i++) {
      itoa(cfgsig1[i].p, 2, s1);
      itoa(cfgsig1[i].n, 2, s2);
      itoa(cfgsig1[i].lam, 2, s3);
      for (int j = 0; j < dimsig2; j++) {
        itoa(cfgsig2[j].p, 2, t1);
        itoa(cfgsig2[j].n, 2, t2);
        itoa(cfgsig2[j].lam, 2, t3);
        fprintf(fpmat, "(%s,%s,%s)(%s,%s,%s)\t", s1, s2, s3, t1, t2, t3);
      }
      fprintf(fpmat, "\n");
    }
    for (int i = 0; i < dimsig2; i++) {
      itoa(cfgsig2[i].p, 2, s1);
      itoa(cfgsig2[i].n, 2, s2);
      itoa(cfgsig2[i].lam, 2, s3);
      for (int j = 0; j < dimsig3; j++) {
        itoa(cfgsig3[j].p, 2, t1);
        itoa(cfgsig3[j].n, 2, t2);
        itoa(cfgsig3[j].lam, 2, t3);
        fprintf(fpmat, "(%s,%s,%s)(%s,%s,%s)\t", s1, s2, s3, t1, t2, t3);
      }
      fprintf(fpmat, "\n");
    }
    fclose(fpmat);
    */
  // ************************************

  // 再求本征值, h 是 double** 哈密顿矩阵, 维数 dim,
  //  eigen_values 存算出来的本征值(排序), eigen_vec
  // 本征向量
  double _Complex eigen_values[dim];
  double _Complex eigen_vec[dim][dim];

  double h1[dim][dim];
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      h1[i][j] = h[i][j];
    }
  }
  GetEigenValuesAndVectors(dim, (double *)h1, eigen_values,
                           (double _Complex *)eigen_vec);

  //     double eigen_values[dim];
  //     double eigen_vec[dim][dim];
  //  gsl 库算不对特征向量, 不知道为什么. 只能用
  //  eigen3 库
  //   gsl_matrix *hmat = gsl_matrix_alloc(dim,
  //   dim); for (int i = 0; i < dim; i++) {
  //     for (int j = 0; j < dim; j++) {
  //       gsl_matrix_set(hmat, i, j, h[i][j]);
  //     }
  //   }

  //   gsl_vector *eval = gsl_vector_alloc(dim);
  //   gsl_matrix *evec = gsl_matrix_alloc(dim,
  //   dim);

  //   gsl_eigen_symmv_workspace *workspace =
  //   gsl_eigen_symmv_alloc(dim);
  //   gsl_eigen_symmv(hmat, eval, evec, workspace);
  //   gsl_eigen_symmv_sort(eval, evec,
  //   GSL_EIGEN_SORT_VAL_ASC);

  //   for (int i = 0; i < dim; i++) {
  //     eigen_values[i] = gsl_vector_get(eval, i);
  //     for (int j = 0; j < dim; j++) {
  //       eigen_vec[i][j] = gsl_matrix_get(evec, i,
  //       j);
  //     }
  //   }

  //   gsl_vector_free(eval);
  //   gsl_matrix_free(hmat);
  //   gsl_matrix_free(evec);
  //   gsl_eigen_symmv_free(workspace);

  fprintf(fpcreate, "E(MeV)\tEx(Mev)\t2J\n");

  // ********************************************************************************2
  double y[NUM_LEVELS], yn[NUM_LEVELS];
  int pm[NUM_LEVELS], kplotmax = 0, Jtimes2[NUM_LEVELS];
  // ********************************************************************************2
  // 角动量
  for (int k = 0; k < dim; k++) {
    if (k == 4) {
      printf("\n");
    }
    double eigen_v[dim];
    double eigen_v2[dim]; //  存升降算符作用后的.
    for (int i = 0; i < dim; i++) {
      eigen_v[i] = creal(eigen_vec[k][i]);
      eigen_v2[i] = 0;
    }

    CreateJsqr(eigen_v, eigen_v2, cfglJplus, cfgl, diml, dimlJplus);
    CreateJsqr(eigen_v + diml, eigen_v2 + diml, cfgsig1Jplus, cfgsig1, dimsig1,
               dimsig1Jplus);
    CreateJsqr(eigen_v + diml + dimsig1, eigen_v2 + diml + dimsig1,
               cfgsig2Jplus, cfgsig2, dimsig2, dimsig2Jplus);
    CreateJsqr(eigen_v + diml + dimsig1 + dimsig2,
               eigen_v2 + diml + dimsig1 + dimsig2, cfgsig3Jplus, cfgsig3,
               dimsig3, dimsig3Jplus);

    double J_sqr = (M2 * M2 / 4.) + M2 / 2.;
    for (int i = 0; i < dim; i++) { //  eigen_v2 与 eigen_v 内积
      J_sqr += eigen_v[i] * eigen_v2[i];
    }
    bool prty;
    for (int i = 0; i < diml; i++) {
      if (fabs(eigen_v[i]) < 1 / sqrt(dim))
        continue;
      prty = (lfromPermut(typepp, cfgl[i].p) + lfromPermut(typenn, cfgl[i].n) +
              lfromPermut(typell, 0)) %
                         2 ==
                     0
                 ? true
                 : false; // 暂时不需要计 Λ
      break;
    }

    level[k].e = creal(eigen_values[k]);
    level[k].j = -1 + sqrt(1 + 4 * J_sqr);
    level[k].prty = prty;

    for (int i = 0; i < dim; i++) {
      level[k].v[i] = eigen_v[i];
    }

    /*
    // ****************fpmat
    if (k >= 0) { //  打印特征向量的分量
      fprintf(fpmat, "\n%d\t%f\n", k, level[k].j);
      for (int i = 0; i < dim; i++) {
        fprintf(fpmat, "%f\t%f\n", eigen_v[i], eigen_v2[i]);
      }
      fprintf(fpmat, "\n");
    }
      fclose(fpmat);
*/

    if (prty) {
      fprintf(fpcreate, "%.3f\t%.3f\t%.3f+\n", creal(eigen_values[k]),
              creal(eigen_values[k]) - creal(eigen_values[0]),
              -1 + sqrt(1 + 4 * J_sqr));

      // ********************************************************************************3
      if (k < NUM_LEVELS) {
        y[k] = creal(eigen_values[k]) - creal(eigen_values[0]);
        yn[k] = y[k];
        Jtimes2[k] = fromdoubletoint(-1 + sqrt(1 + 4 * J_sqr));
        pm[k] = 1; // "+"
        kplotmax = k;
      }
      // ********************************************************************************3

    } else {
      fprintf(fpcreate, "%.3f\t%.3f\t%.3f-\n", creal(eigen_values[k]),
              creal(eigen_values[k]) - creal(eigen_values[0]),
              -1 + sqrt(1 + 4 * J_sqr));
      // ********************************************************************************3
      if (k < NUM_LEVELS) {
        y[k] = creal(eigen_values[k]) - creal(eigen_values[0]);
        yn[k] = y[k];
        Jtimes2[k] = fromdoubletoint(-1 + sqrt(1 + 4 * J_sqr));
        pm[k] = 0; // "-"
        kplotmax = k;
      }
      // ********************************************************************************3
    }
  }
  // ********************************************************************************5
  char name[50];
  if (n_lam == 0) {
    sprintf(name, "%d%s", A, isotope[Z - 2]);
  } else {
    sprintf(name, "%d%sΛ", A, isotope[Z - 2]);
  }
  FILE *gnuplotPipe = popen("gnuplot", "w");

  fprintf(gnuplotPipe, "set encoding utf8 \n");
  fprintf(gnuplotPipe, "set terminal svg enhanced font \"Times,15\"\n");
  fprintf(gnuplotPipe, "set output \"result/%s.svg\"\n", name);
  fprintf(gnuplotPipe, "emin=0 \n");
  fprintf(gnuplotPipe, "emax=%.3f \n", y[kplotmax]);
  fprintf(gnuplotPipe, "xmin=-0.2 \n");
  fprintf(gnuplotPipe, "xmax=3.6 \n");
  fprintf(gnuplotPipe, "yspan=emax-emin \n");
  fprintf(gnuplotPipe, "ymin=emin-(yspan*0.1) \n");
  fprintf(gnuplotPipe, "ymax=emax+(yspan*0.1) \n");
  fprintf(gnuplotPipe, "unset xtics \n");
  fprintf(gnuplotPipe, "unset key \n");
  fprintf(gnuplotPipe, "set title \"%s\"\n",
          name); //  pathcreate
                 //  前三个字符是上层路径,
                 //  不在图中输出
  for (int k = 0; k < kplotmax; k++) {
    if (yn[k + 1] - yn[k] < y[kplotmax] * width) {
      int j;
      for (j = k + 1; j < kplotmax; j++) {
        if (yn[j + 1] - yn[j] >= y[kplotmax] * width) {
          break;
        }
      }
      for (int l = k; l <= j; l++) {
        yn[l] = (y[j] + y[k]) / 2 -
                y[kplotmax] * width * (((double)(j - k)) / 2 - (double)(l - k));
      }
      k = j;
    }
  }
  for (int k = 0; k <= kplotmax; k++) {
    fprintf(gnuplotPipe,
            "set arrow from 1.2,%.3f to 2.2,%.3f "
            "nohead \n",
            y[k], y[k]);
    fprintf(gnuplotPipe, "set label \"%.3f\" at 0.8,%.3f \n", y[k], yn[k]);
    if (Jtimes2[k] % 2 != 0) {
      if (pm[k] == 0) {
        fprintf(gnuplotPipe, "set label \"%d/2-\" at 2.3,%.3f \n", Jtimes2[k],
                yn[k]);
      } else {
        fprintf(gnuplotPipe, "set label \"%d/2+\" at 2.3,%.3f \n", Jtimes2[k],
                yn[k]);
      }
    } else {
      if (pm[k] == 0) {
        fprintf(gnuplotPipe, "set label \"%d-\" at 2.3,%.3f \n", Jtimes2[k] / 2,
                yn[k]);
      } else {
        fprintf(gnuplotPipe, "set label \"%d+\" at 2.3,%.3f \n", Jtimes2[k] / 2,
                yn[k]);
      }
    }
  }
  fprintf(gnuplotPipe, "plot [xmin:xmax][ymin:ymax] NaN\n");
  //   fprintf(gnuplotPipe, "set terminal postscript\n");
  //   fprintf(gnuplotPipe, "set output \"%s.eps\"\n", pathcreate);
  //   fprintf(gnuplotPipe, "replot\n");
  //   fprintf(gnuplotPipe, "set terminal x11\n");
  //   fprintf(gnuplotPipe, "set output\n");

  pclose(gnuplotPipe);

  char cmd[200];
  sprintf(cmd, "svg2pdf result/%s.svg result/%s.pdf", name, name);
  system(cmd);
  // ********************************************************************************5

  //  在命令行也输出一次
  printf("\t\t\t%s\n", name);
  printf("E(MeV)\tEx(Mev)\t2J\tprty\n");
  for (int i = 0; i < 10; i++) {
    printf("%f\t%f\t%f\t%d\n", level[i].e, level[i].e - level[0].e, level[i].j,
           level[i].prty);
  }

  free(h[0]);
  free(h);
  free(cfgl);
  free(cfgsig1);
  free(cfgsig3);
  free(cfglJplus);
  free(cfgsig1Jplus);
  free(cfgsig3Jplus);

  fclose(fpcreate);
  fclose(fp);
  fclose(fpH);
  return 0;
}

bool shellBlam1(int A, int Z, int n_lam, double *val, struct level *energylevel,
                int inter) {
  double e1;
  bool lam = true;
  if (n_lam != 1) {
    printf("n_lam must be 1!\n");
    lam = false;
    return lam;
  }
  shell(A - 1, Z, 0, energylevel, inter);
  e1 = energylevel[0].e;
  shell(A, Z, 1, energylevel, inter);
  double e2 = energylevel[0].e;
  printf("Λ binding energy of %d%sΛ is %f.\n", A, isotope[Z - 2], e1 - e2);
  *val = e1 - e2;
  return lam;
}

bool shellplot(int A, int Z, int n_lam, int numberoflevels,
               struct level *level1, struct level *level2) {
  bool lam = true;
  char name[50];
  if (n_lam == 0) {
    sprintf(name, "%d%s", A, isotope[Z - 2]);
  } else {
    sprintf(name, "%d%sΛ", A, isotope[Z - 2]);
  }
  FILE *fpplot;
  if (!(fpplot = fopen((std::string("plot/") + name + ".plt").c_str(), "w+"))) {
    printf(".plt cannot open!\n");
    lam = false;
    return lam;
  }
  fprintf(fpplot, "name=%s\n", name);
  fprintf(fpplot, "num_columns=2\n");
  shell(A, Z, n_lam, level1, 0);
  shell(A, Z, n_lam, level2, 1);
  double emax = std::max(level1[numberoflevels - 1].e - level1[0].e,
                         level2[numberoflevels - 1].e - level2[0].e);
  fprintf(fpplot, "emax=%f\n", emax);
  fprintf(fpplot, "Gogny %d\n", numberoflevels);
  for (int i = 0; i < numberoflevels; i++) {
    if (level1[i].prty) {
      if (fromdoubletoint(level1[i].j) % 2 != 0) {
        fprintf(fpplot, "%f %d/2+\n", level1[i].e - level1[0].e,
                fromdoubletoint(level1[i].j));
      } else {
        fprintf(fpplot, "%f %d+\n", level1[i].e - level1[0].e,
                fromdoubletoint(level1[i].j) / 2);
      }
    } else {
      if (fromdoubletoint(level1[i].j) % 2 != 0) {
        fprintf(fpplot, "%f %d/2-\n", level1[i].e - level1[0].e,
                fromdoubletoint(level1[i].j));
      } else {
        fprintf(fpplot, "%f %d-\n", level1[i].e - level1[0].e,
                fromdoubletoint(level1[i].j) / 2);
      }
    }
  }
  fprintf(fpplot, "(8-16)2BME %d\n", numberoflevels);
  for (int i = 0; i < numberoflevels; i++) {
    if (level2[i].prty) {
      if (fromdoubletoint(level1[i].j) % 2 != 0) {
        fprintf(fpplot, "%f %d/2+\n", level2[i].e - level2[0].e,
                fromdoubletoint(level2[i].j));
      } else {
        fprintf(fpplot, "%f %d+\n", level2[i].e - level2[0].e,
                fromdoubletoint(level2[i].j) / 2);
      }
    } else {
      if (fromdoubletoint(level1[i].j) % 2 != 0) {
        fprintf(fpplot, "%f %d/2-\n", level2[i].e - level2[0].e,
                fromdoubletoint(level2[i].j));
      } else {
        fprintf(fpplot, "%f %d-\n", level2[i].e - level2[0].e,
                fromdoubletoint(level2[i].j) / 2);
      }
    }
  }
  fclose(fpplot);
  return lam;
}

void fullplot(FILE *fp, FILE *fpexp) {
  char name[50];
  int num_columns;
  double emax, t;
  fscanf(fp, "name=%s", name);
  moveToNextLine(fp);
  fscanf(fp, "num_columns=%d", &num_columns);
  num_columns += 2;
  moveToNextLine(fp);
  fscanf(fp, "emax=%lf", &emax);
  moveToNextLine(fp);
  fscanf(fpexp, "emax=%lf", &t);
  moveToNextLine(fpexp);
  if (t > emax) {
    emax = t;
  }
  FILE *gnuplotPipe;
  if (!(gnuplotPipe = popen("gnuplot", "w"))) {
    throw std::ios_base::failure("gnuplot");
  }
  fprintf(gnuplotPipe, "set encoding utf8 \n");
  fprintf(gnuplotPipe,
          "set terminal svg size 640,480 enhanced font \"Times,15\"\n");
  fprintf(gnuplotPipe, "set output \"./fullplot/%s.svg\"\n", name);
  fprintf(gnuplotPipe, "emin=0 \n");
  fprintf(gnuplotPipe, "emax=%.3f \n", emax);
  fprintf(gnuplotPipe, "xmin=0 \n");
  fprintf(gnuplotPipe, "xmax=%d \n", 3 * num_columns);
  fprintf(gnuplotPipe, "yspan=emax-emin \n");
  fprintf(gnuplotPipe, "ymin=emin-(yspan*0.2) \n");
  fprintf(gnuplotPipe, "ymax=emax+(yspan*0.1) \n");
  fprintf(gnuplotPipe, "unset xtics \n");
  fprintf(gnuplotPipe, "unset key \n");
  fprintf(gnuplotPipe, "set title \"%s\"\n", name);
  for (int i = 0; i < num_columns; i++) {
    char columnname[50];
    int num_levels;
    double y[20], yn[20];
    char J[20][10];
    if (i == 0 || i == 1) {
      fscanf(fpexp, "%s", columnname);
      fscanf(fpexp, "%d", &num_levels);
      moveToNextLine(fpexp);
      for (int k = 0; k < num_levels; k++) {
        fscanf(fpexp, "%lf", &y[k]);
        yn[k] = y[k];
        fscanf(fpexp, "%s", J[k]);
        moveToNextLine(fpexp);
      }
    } else {
      fscanf(fp, "%s", columnname);
      fscanf(fp, "%d", &num_levels);
      moveToNextLine(fp);
      for (int k = 0; k < num_levels; k++) {
        fscanf(fp, "%lf", &y[k]);
        yn[k] = y[k];
        fscanf(fp, "%s", J[k]);
        moveToNextLine(fp);
      }
    }
    // 方案: 记录字体大小, 检测碰撞. 如果相邻有碰撞则分开至恰好不碰撞为止.
    // 如果整体碰撞, 那就改变图片高度, 纵向整体放大.
    for (int k = 0; k < num_levels - 1; k++) {
      if (yn[k + 1] - yn[k] < emax * width) {
        yn[k + 1] = (y[k + 1] + y[k]) / 2 + emax * width / 2;
        yn[k] = (y[k + 1] + y[k]) / 2 - emax * width / 2;
      }
    }
    fprintf(gnuplotPipe, "set label \"%s\" at %.1f,%.3f \n", columnname,
            1.3 + 3 * i, 0.0 - emax * 0.1);
    for (int k = 0; k < num_levels; k++) {
      fprintf(gnuplotPipe, "set arrow from %.1f,%.3f to %.1f,%.3f nohead \n",
              1.1 + 3 * i, y[k], 2.1 + 3 * i, y[k]);
      fprintf(gnuplotPipe, "set label \"%.3f\" at %.1f,%.3f \n", y[k],
              0.2 + 3 * i, yn[k]);
      fprintf(gnuplotPipe, "set label \"%s\" at %.1f,%.3f \n", J[k],
              2.3 + 3 * i, yn[k]);
    }
  }
  fprintf(gnuplotPipe, "plot [xmin:xmax][ymin:ymax] NaN\n");
  pclose(gnuplotPipe);

  char cmd[200];
  sprintf(cmd, "svg2pdf ./fullplot/%s.svg ./fullplot/%s.pdf", name, name);
  system(cmd);
}

bool shellBlam2(int A, int Z, int n_lam, double *val, struct level *energylevel,
                int inter) {
  double e1;
  bool lam = true;
  if (n_lam != 2) {
    printf("n_lam must be 2!\n");
    lam = false;
    return lam;
  }
  shell(A - 2, Z, 0, energylevel, inter);
  e1 = energylevel[0].e;
  double j1 = energylevel[0].j;
  shell(A - 1, Z, 1, energylevel, inter);
  double b;
  if (fabs(j1) < 1e-3) {
    b = -(energylevel[0].e - e1);
  } else {
    //  角动量2J+1加权平均, level.j 已经是两倍了
    b = -(energylevel[0].e * (energylevel[0].j + 1) +
          energylevel[1].e * (energylevel[1].j + 1)) /
            (energylevel[0].j + energylevel[1].j + 2) +
        e1;
  }
  printf("double-Λ binding energy of %d%sΛΛ is %f.\n", A, isotope[Z - 2],
         2 * b + 0.67);
  *val = 2 * b + 0.67;
  return lam;
}

int main() {
  int n_lam = 1;
  int A = 7;
  int Z = 3;

  int ext = 1;
  struct level *energylevel;
  struct level *energylevel2;
  energylevel = (struct level *)malloc(
      sizeof(struct level) *
      dim_occup); //  这样如果出错就会跟 shell 里面同时报错

  // 存本征值, 设向量维数不超过 dim_occup
  for (int i = 0; i < dim_occup; i++) {
    energylevel[i].v = (double *)malloc(sizeof(double) * dim_occup);
  }
  struct config *cfgl;
  cfgl = (struct config *)malloc(sizeof(struct config) * dim_occup);
  int dim_cfgl = 0;

  int interaction = 0;
  system("mkdir result"); //  新建一个路径存结果
  double val;
  while (ext != 0) {
    printf("input 1 for shell model, 2 for Λ binding energy, "
           "3 for double-Λ binding energy, "
           "4 for all p-shell Λ binding energies, "
           "5 for all p-shell ΛΛ binding energies, "
           "6 for creating .plt file for energy level plot, "
           "7 for BE2 between the first 5 states for a given shell, "
           "0 exit\n");

    int numberoflevels;
    int choose = 7;
    scanf("%d", &choose);
    switch (choose) {
    case 7: {
      printf("input &A &Z &n_lam:\n");
      scanf("%d %d %d", &A, &Z, &n_lam);
      printf("input interaction type:(0 for gogny, 1 for fit)\n");
      scanf("%d", &interaction);

      shell(A, Z, n_lam, energylevel, interaction);
      makeBasis(A, Z, n_lam, &dim_cfgl, cfgl);

      //      for (int i = 0; i < dim_cfgl; i++) {
      //        printf("%d %d\n", cfgl[i].p, cfgl[i].n);
      //      }

      double be2mat[dim_cfgl][dim_cfgl];
      for (int i = 0; i < dim_cfgl; i++) {
        for (int j = i; j < dim_cfgl; j++) {
          if (j == i) {
            be2mat[i][i] = be2diagfromPermut(typenn, cfgl[i].n) +
                           be2diagfromPermut(typepp, cfgl[i].p);
            //            printf("be2nondiag i = %d, be2 for %d %d:%f\n", i,
            //            cfgl[i].n,
            //                   cfgl[i].p, be2mat[i][j]);
          } else if (difcount(cfgl[i].lam, cfgl[j].lam) == 0) {

            if (difcount(cfgl[i].p, cfgl[j].p) == 0 &&
                difcount(cfgl[i].n, cfgl[j].n) == 2) {
              be2mat[i][j] = be2mat[j][i] =
                  be2nondiagfromPermut(typenn, cfgl[i].n, cfgl[j].n);
              //              printf("be2nondiag i = %d, j = %d, be2 for %d %d:
              //              %f\n", i, j,
              //                     cfgl[i].n, cfgl[j].n, be2mat[i][j]);

            } else if (difcount(cfgl[i].p, cfgl[j].p) == 2 &&
                       difcount(cfgl[i].n, cfgl[j].n) == 0) {
              be2mat[i][j] = be2mat[j][i] =
                  be2nondiagfromPermut(typepp, cfgl[i].p, cfgl[j].p);
            } // should only have in total (p & n) one difference.
            else {
              be2mat[i][j] = be2mat[j][i] = 0;
            }
          } else {
            be2mat[i][j] = be2mat[j][i] = 0;
          }
        }
      }

      for (int i = 0; i < dim_cfgl; i++) {
        for (int j = 0; j < dim_cfgl; j++) {
          printf("%.2f\t", be2mat[i][j]);
        }
        printf("\n");
      }

      for (int i = 0; i < dim_cfgl; i++) {
        printf("p ");
        bin(cfgl[i].p);
        printf("\tn ");
        bin(cfgl[i].n);
        printf("\tlam ");
        bin(cfgl[i].lam);
        printf("\n");
      }

      int n_ini, n_fin;
      int next = 0;
      while (next == 0) {

        do {
          printf("input initial and final state, &n_ini &n_fin:\n");
          scanf("%d %d", &n_ini, &n_fin);
        } while (!(n_fin < n_ini && n_fin >= 0));

        double be2 = 0;
        for (int i = 0; i < dim_cfgl; i++) {
          for (int j = 0; j < dim_cfgl; j++) {
            be2 += energylevel[n_ini].v[i] * be2mat[i][j] *
                   energylevel[n_fin].v[j];
          }
        }

        for (int i = 0; i < dim_cfgl; i++) {
          printf("%f\t%f\n", energylevel[n_ini].v[i], energylevel[n_fin].v[i]);
        }

        int M2 = 0;
        if (A % 2 == 1) {
          M2 = -1;
        }

        be2 = be2 * 41.471 / hw(A, Z) /
              gsl_sf_coupling_3j((int)round(energylevel[n_ini].j), 4,
                                 (int)round(energylevel[n_fin].j), M2, 0, -M2);
        be2 = be2 * be2 / (energylevel[n_ini].j + 1); // 平方
        printf("BE2 from %d to %d equals:\n%f\n", n_ini, n_fin, be2);

        printf("Continue 0 or not 1:\n");
        scanf("%d", &next);
      }

      break;
    }
    case 1: {
      printf("input &A &Z &n_lam:\n");
      scanf("%d %d %d", &A, &Z, &n_lam);
      printf("input interaction type:(0 for gogny, 1 for fit)\n");
      scanf("%d", &interaction);
      shell(A, Z, n_lam, energylevel, interaction);
      break;
    }
    case 2: {
      do {
        printf("input &A &Z &n_lam:\n");
        scanf("%d %d %d", &A, &Z, &n_lam);
        printf("input interaction type:(0 for gogny, 1 for fit)\n");
        scanf("%d", &interaction);
      } while (!shellBlam1(A, Z, n_lam, &val, energylevel, interaction));
      break;
    }
    case 3: {
      do {
        printf("input &A &Z &n_lam:\n");
        scanf("%d %d %d", &A, &Z, &n_lam);
        printf("input interaction type:(0 for gogny, 1 for fit)\n");
        scanf("%d", &interaction);
      } while (!shellBlam2(A, Z, n_lam, &val, energylevel, interaction));
      break;
    }
    case 4: {
      FILE *fpB1;
      if (!(fpB1 = fopen("result/BΛ.txt", "w+"))) {
        fprintf(stderr, "Can't find file BΛ!\n");
        exit(1);
      }
      printf("input interaction type:(0 for gogny, 1 for fit)\n");
      scanf("%d", &interaction);
      if (interaction != 0 && interaction != 1) {
        printf("interaction error in case 4.\n");
        break;
      }
      fprintf(fpB1, "\tB(Λ)\n");
      for (int Z = 2; Z <= 8; Z++) {
        for (int N = 2; N <= 8; N++) {
          int A = Z + N + 1;
          shellBlam1(A, Z, 1, &val, energylevel, interaction);
          fprintf(fpB1, "%d%sΛ\t%f\n", A, isotope[Z - 2], val);
        }
      }
      fclose(fpB1);
    } break;
    case 5: {
      FILE *fpB2;
      if (!(fpB2 = fopen("result/BΛΛ.txt", "w+"))) {
        fprintf(stderr, "Can't find file BΛΛ!\n");
        exit(1);
      }
      printf("input interaction type:(0 for gogny, 1 for fit)\n");
      scanf("%d", &interaction);
      if (interaction != 0 && interaction != 1) {
        printf("interaction error in case 5.\n");
        break;
      }
      fprintf(fpB2, "\tB(ΛΛ)\n");
      for (int Z = 2; Z <= 8; Z++) {
        for (int N = 2; N <= 8; N++) {
          int A = Z + N + 2;
          shellBlam2(A, Z, 2, &val, energylevel, interaction);
          fprintf(fpB2, "%d%sΛΛ\t%f\n", A, isotope[Z - 2], val);
        }
      }
      fclose(fpB2);
    } break;
    case 6: {
      system("mkdir plot");
      system("mkdir fullplot");
      energylevel2 = (struct level *)malloc(sizeof(struct level) * dim_occup);
      do {
        do {
          printf("input &A &Z &n_lam, n_lam must be 1:\n");
          if (scanf("%d %d %d", &A, &Z, &n_lam) != 3) {
            printf("input A Z n_lam\n");
            continue;
          }
        } while (n_lam != 1);
        printf("input number of energy levels to plot:\n");
        scanf("%d", &numberoflevels);
      } while (
          !shellplot(A, Z, n_lam, numberoflevels, energylevel, energylevel2));
      free(energylevel2);

      char plotfilename[FILENAME_MAX], plotfileexpname[FILENAME_MAX];
      sprintf(plotfilename, "%s/%d%sΛ.plt", PATHPLOT, A, isotope[Z - 2]);
      sprintf(plotfileexpname, "%s/%d%sΛexp.plt", PATHPLOTEXP, A,
              isotope[Z - 2]);
      FILE *fp;
      if (!(fp = fopen(plotfilename, "r"))) {
        throw std::ios_base::failure("fullplot filename error.");
      }

      FILE *fpexp;
      if (!(fpexp = fopen(plotfileexpname, "r"))) {
        throw std::ios_base::failure("fullplot expfilename error.");
      }
      fullplot(fp, fpexp);
      fclose(fp);
      fclose(fpexp);
    } break;
    case 0:
      ext = 0;
      break;
    default:
      break;
    }
  }
  for (int i = 0; i < dim_occup; i++) {
    free(energylevel[i].v);
  }
  free(energylevel);

  free(cfgl);
  return 0;
}

// 上面的
// “//
// ********************************************************************************1,2,3,4,5”
// 表示之间的是与画图相关的代码

/*
int main() {
  //   // *****************在程序内 Debug shell 使用*************

  //   int n_lam = 1;
  //   int A = 7, Z = 3;

  struct level *energylevel;
  energylevel = malloc(sizeof(struct level) *
                       dim_occup); //  这样如果出错就会跟 shell 里面同时报错
  //   shell(A, Z, n_lam, energylevel);
  //   // *******************************************

  const char *isotope[7] = {"He", "Li", "Be", "B", "C", "N", "O"};
  double val;
  FILE *fpB1;
  if (!(fpB1 = fopen("../BΛ.txt", "wb+"))) {
    fprintf(stderr, "Can't find file BΛ!\n");
    exit(1);
  }
  fprintf(fpB1, "\tB(Λ)\n");
  for (int Z = 2; Z <= 8; Z++) {
    for (int N = 2; N <= 8; N++) {
      int A = Z + N + 1;
      shellBlam1(A, Z, 1, &val, energylevel);
      fprintf(fpB1, "%d%sΛ\t%f\n", A, isotope[Z - 2], val);
    }
  }
  fclose(fpB1);
  return 0;
}
*/
