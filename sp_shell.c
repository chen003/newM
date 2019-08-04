//
//  sp'_shell.c
//  nuclear_A1B1
//
//  Created by xingyan chen on 09/03/2017.
//  Copyright © 2017 xingyan chen. All rights reserved.
//

#include "sp_shell.h"

int power(int a, int b) {
  int n = 1;
  for (int i = 0; i < b; i++) {
    n = n * a;
  }
  return n;
}
static double e[n_orbtJ];

void set_singleE(FILE *fp) { //把 pn 的单粒子能读到数组 e 中, 让 singleE 能访问
  fseek(fp, 0, SEEK_SET);
  moveToNextLine(fp);
  fscanf(fp, "%*d ");
  for (int i = 0; i < n_orbtJ; i++) {
      fscanf(fp, "%lf ",e + i);
  }
  fseek(fp, 0, SEEK_SET);
}

double singleE(int type, int a) {
  double E = 0.;
  switch (type) {
  case 1:
  case 2:
    E = e[numfromA(a, 1) - 1];
    break;
  case 3:
    E = -3.12;  // 暂且这样处理, 以后文件输入
    break;
  }
  return E;
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//测试 main 函数求哈密顿矩阵的算法, 就让这些矩阵元函数 return 参数之和.
//原来的函数见 setup_TBME.c  double TBME(int type, int c1, int c2, int
// c3,
// int c4 ){
//    if (type==4)
//        return c1+c2+(c3+c4);
//    else return 0;
//}

double MEV3(int typeAB, int A1, int A2, int B, int a1, int a2, int b) {
  double ME = 0.;
  //只计入 type AB的TBME, 另外两类:一类是 type AA 的TBME, 还有单粒子能级.
  //分开算

  if (A1 == a1)
    ME = ME + TBME(typeAB, A2, B, a2, b);
  if (A2 == a2)
    ME = ME + TBME(typeAB, A1, B, a1, b);
  if (A1 == a2)
    ME = ME - TBME(typeAB, A2, B, a1, b);
  if (A2 == a1)
    ME = ME - TBME(typeAB, A1, B, a2, b);

  return ME;
}

//计算 MEV4用到
//用来存考察的那一对下标是否相等 AB 两部分单粒子算符返回的下标
struct para {
  bool equal;
  int m_;
  int m;
};

//数组的序号存储相等那一对参数的位置, 可以返回矩阵元的符号.
//序号0代表a1'=a1.  1代表a1'=a2. 2代表a2'=a2. 3代表a2'=1.

//计算用到
void MEV4_getpara(int A1, int A2, int a1, int a2, struct para paraA[4]);

//计算形如<a1',a2',b1',b2'|V_typeAB|a1,a2,b1,b2>的矩阵元
double MEV4(int typeAB, int A1, int A2, int B1, int B2, int a1, int a2, int b1,
            int b2) {
  double ME = 0.;
  //数组的序号存储相等那一对参数的位置, 可以返回矩阵元的符号.
  //序号0代表a1'=a1.  1代表a1'=a2. 2代表a2'=a2. 3代表a2'=a1.
  struct para paraA[4], paraB[4];

  MEV4_getpara(A1, A2, a1, a2, paraA);

  MEV4_getpara(B1, B2, b1, b2, paraB);

  for (int i = 0; i < 4; i++) {
    if (paraA[i].equal == false) {
      continue;
    }
    for (int j = 0; j < 4; j++) {
      if (paraB[j].equal == false) {
        continue;
      }
      ME = ME +
           ((i + j) % 2 == 0 ? 1 : -1) *
               TBME(typeAB, paraA[i].m_, paraB[j].m_, paraA[i].m, paraB[j].m);
    }
  }
  return ME;
}

//计算 MEV4单粒子算符返回的参数
void MEV4_getpara(int A1, int A2, int a1, int a2, struct para paraA[4]) {
  //初始化
  for (int i = 0; i < 4; i++) {
    paraA[i].equal = false;
  }
  //四对
  if (A1 == a1) {
    paraA[0].equal = true;
    paraA[0].m_ = A2;
    paraA[0].m = a2;
  }
  if (A1 == a2) {
    paraA[1].equal = true;
    paraA[1].m_ = A2;
    paraA[1].m = a1;
  }
  if (A2 == a2) {
    paraA[2].equal = true;
    paraA[2].m_ = A1;
    paraA[2].m = a1;
  }
  if (A2 == a1) {
    paraA[3].equal = true;
    paraA[3].m_ = A1;
    paraA[3].m = a2;
  }
}

int max(int i, int j) { return i < j ? j : i; }

//用于按顺序生成 pn 表象 J-scheme TBME. 如 text.int 中不重复也不少的顺序
void printConfig(int type, int num_orbits) {
  switch (type) {
  case 1:
    for (int i1 = 1; i1 <= num_orbits; i1++) {
      for (int i2 = i1; i2 <= num_orbits; i2++) {
        for (int i3 = i1; i3 <= num_orbits; i3++) {
          for (int i4 = max(i3, i2); i4 <= num_orbits; i4++) {
            printf("%d %d %d %d\n", i1, i2, i3, i4);
          }
        }
      }
    }

    break;
  case 2:
    for (int i1 = 1; i1 <= num_orbits; i1++) {
      for (int i2 = i1 + num_orbits; i2 <= 2 * num_orbits; i2++) {
        for (int i3 = i1; i3 <= num_orbits; i3++) {
          for (int i4 = i3 + num_orbits; i4 <= 2 * num_orbits; i4++) {
            printf("%d %d %d %d\n", i1, i2, i3, i4);
          }
        }
      }
    }
    break;
  case 3:
    for (int i1 = 1 + num_orbits; i1 <= 2 * num_orbits; i1++) {
      for (int i2 = i1; i2 <= 2 * num_orbits; i2++) {
        for (int i3 = i1; i3 <= 2 * num_orbits; i3++) {
          for (int i4 = max(i3, i2); i4 <= 2 * num_orbits; i4++) {
            printf("%d %d %d %d\n", i1, i2, i3, i4);
          }
        }
      }
    }
  default:
    break; //不能嵌套定义吗?
  }
}
