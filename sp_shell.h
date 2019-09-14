//
//  sp'_shell.h
//  nuclear_A1B1
//
//  Created by xingyan chen on 09/03/2017.
//  Copyright © 2017 xingyan chen. All rights reserved.
//

#ifndef sp__shell_h
#define sp__shell_h

#include <complex.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#define typepn 4
#define typepl 5
#define typenl 6
#define typepp 1
#define typenn 2
#define typell 3
#define typeSig1p                                                              \
  7 // 实际调用时方便起见, 让核子在超子前面. 比如Σ+p 或 Σ-n
    //, 实际上是 p 和 n 在 Σ前面
#define typeSig1n 8       // Σ+n 或 Σ-p
#define typeSig2p 9       // Σ0p 或 Σ0n
#define typeSig2pSig1n 10 // Σ0p-Σ+n 或 Σ-p-Σ0n
#define typeLampSig2p 11  // Λp-Σ0p
#define typeLampSig1n 12  // Λp-Σ+n
#define typeLamnSig2n 13  // Λn-Σ0n
#define typeLamnSig3p 14  // Λn-Σ-p

// 除了核子数和 int 相互作用文件外全部输入集中在这里
#define n_orbtM 6  //  p 或 n 的 M-scheme 轨道个数
#define n_orbtJ 2  // J-scheme 轨道个数
#define n_orbtMH 2 // Λ 的 M-scheme 轨道个数
#define n_orbtJH 1
#define n_type 6                   //相互作用类型. pp nn ΛΛ pn nΛ Λp
#define num_lines 34               // pn.int 文件
#define numlinesH 5                // 超核的相互作用文件行数
#define PATHsp "../newM/sp/p.sp"   // sp 壳文件
#define PATHspH "../newM/sp/sH.sp" // Λ 的 s 壳文件

#define efchargeP 1.3 // BE2 有效电荷
#define efchargeN 0.5

#endif /* sp__shell_h */

#if defined(__cplusplus)
extern "C" {
#endif

// 二进制转换
void itoa(int n, int base, char *buf);

// print
void bin(unsigned n);

// 生成 TBME, 存入数组. fp 是  pn.int 文件
void makeTBME(FILE *fp, FILE *fpH, int l_or_h);

// setup_TBME 里调用的, 用于读 int 文件
int **getdatac(FILE *fp);
double *getdatad(FILE *fp);
int **Hgetdatac(FILE *fp);
double *Hgetdatad(FILE *fp);

// setup_TBME 里求和用的来自 readsp.c 的函数
int moveToNextLine(FILE *fp);
void setMschemeID2JschemeID(
    FILE *fp,
    FILE *fpH); //使用 mfromA, lfromA, jfromA, JschemeJ, numfromA 先运行
int mfromA(int type, int A);       //从 M-scheme 轨道得到 m 值(两倍)
int lfromA(int type, int A);       //一倍
int jfromA(int type, int A);       //两倍
int JschemeJ(int type, int a);     // 得到 Jscheme 角动量
int typeAfromAB(int type, int lr); // lr 1 是左 2 是右
int numfromA(int A, int type);

// BE2
double efcharge(int type);                 // effective charge from type
double be2diagfromPermut(int type, int w); // p 壳的 BE2
double be2nondiagfromPermut(int type, int w, int v);

double hw(int A, int Z);

int power(int a, int b); //整数次幂

struct config {
  int p;
  int n;
  int lam; //用__builtin 函数加位运算, 更高效.
};

//用于生成组态, bit 位和数组相互转换. bitConfig.c
unsigned int nextbitpermut(unsigned int v);
int mfromPermut(int type, int w);
int lfromPermut(int type, int w);
double singlePermutE(int w, int type);
int difcount(int p, int q);
int leaddifbit(int p, int q);
int leastdifbit(int p, int q);

// Two-Body Operator, act on different kinds
double TBOME(int type, int p1, int q1, int p2, int q2);
double TBOME2(int type, int p1, int q1, int p2, int q2);
// act on the same kind
double TBOMEnn(int type, int p, int q);

double vJplus(int type, int cfgJplus, int cfg, double eigen_v);
double vJminus(int type, int cfgJplus, int cfg, double eigen_vJplus);

//总角动量升降算符的作用系数, a 标记升降前的态
static inline double c_Jplus(int type, int a) {
  return sqrt(jfromA(type, a) * (jfromA(type, a) + 2) -
              mfromA(type, a) * (mfromA(type, a) + 2)) /
         2;
}

static inline double c_Jminus(int type, int a) {
  return sqrt(jfromA(type, a) * (jfromA(type, a) + 2) -
              mfromA(type, a) * (mfromA(type, a) - 2)) /
         2;
}

//单粒子能级 a 的能量
void set_singleE(FILE *fp);
double singleE(int type, int a);

//二体矩阵元<a',b'|V_type|a,b>_nas, type 标记相互作用类型(包括对角元,
//但不包括单粒子能级), AA\BB\CC\AB\BC\AC 分别对应type 1~6, 方便用 type%3. 当是
// AA 时, 应考虑费米子统计导致的二体算符有一个1/2以及算符的対易关系. 避免重复,
//规定c1<=c2.
double TBME(int type, int c1, int c2, int c3, int c4);

double Mtbme_T(int type, int k1m1, int k2m2, int k3m3, int k4m4, int lorh);

//测试用
void writeTBME();

#if defined(__cplusplus)
}
#endif
