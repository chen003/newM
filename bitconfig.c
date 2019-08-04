#include "sp_shell.h"

#define max(a, b) (a > b ? a : b)
#define min(a, b) (a < b ? a : b)

unsigned int nextbitpermut(
    unsigned int
        v) { // https://graphics.stanford.edu/~seander/bithacks.html#NextBitPermutation
  unsigned int t = v | (v - 1); // t gets v's least significant 0 bits set to 1
  // Next set to 1 the most significant bit to change,
  // set to 0 the least significant ones, and add the necessary 1 bits.
  unsigned int w = (t + 1) | (((~t & -~t) - 1) >> (__builtin_ctz(v) + 1));
  return w;
}

int mfromPermut(int type, int w) {  //type 种类粒子的总 m 角动量
  int n_nucleon = __builtin_popcount(w);
  int m = 0;
  int n;
  for (int i = 0; i < n_nucleon; i++) {
    n = __builtin_ctz(w) + 1;
    m += mfromA( type,n);
    w = w & (~0 << n);
  }
  return m;
}

int lfromPermut(int type, int w) { //type 种类粒子的总轨道角动量
  int n_nucleon = __builtin_popcount(w);
  int m = 0;
  int n;
  for (int i = 0; i < n_nucleon; i++) {
    n = __builtin_ctz(w) + 1;
    m += lfromA(type, n);
    w = w & (~0 << n);
  }
  return m;
}

double singlePermutE(int w, int type) { // 从 int 得到多粒子波函数的单粒子能,
                                        // type 标记粒子种类, 暂无区别.
  int n_nucleon = __builtin_popcount(w);
  double e = 0;
  int n;
  for (int i = 0; i < n_nucleon; i++) {
    n = __builtin_ctz(w) + 1;
    e += singleE(type, n);
    w = w & (~0 << n);
  }
  return e;
}

int difcount(int p, int q) { return __builtin_popcount((p ^ q)); } //不同位数

int msbit(int p) { return __builtin_clrsb(0) - __builtin_clrsb(p); } //最大位数

int lsbit(int p) { return __builtin_ctz(p); } // 最小位数

int leaddifbit(int p, int q) { //最高位不同的, 位数. 换成数组要-1
  int frst;
  int s = (p ^ q);
  frst = msbit(s);
  if ((p >> (frst - 1)) % 2 == 0) {
    return -frst;
  } else {
    return frst;
  }
}

int leastdifbit(int p, int q) { // 与 lead 一致, 换成数组要 -1
  int lst;
  int s = (p ^ q);
  lst = __builtin_ctz(s) + 1;
  if ((p >> (lst - 1)) % 2 == 0) {
    return -lst; // 负号表示该位置被 q 占据, 正号表示被 p 占据.
  } else {
    return lst;
  }
}

void setTBO2(int *m1, int *m3, int p1,
             int p2) { //只用于非同种粒子 difcount == 2 的情形, 得到编号分配(属于左矢还是右矢)
  if ((*m1 = leaddifbit(p1, p2)) < 0) {
    *m3 = -*m1;
    *m1 = abs(leastdifbit(p1, p2));
  } else {
    *m3 = abs(leastdifbit(p1, p2));
  }
}

int signTBO2(int m1, int m3, int p) { // 返回反对易关系得到的+-1, 方法是求 p 在 m1 和 m3 位之间有多少占据粒子
  int sign = 1;
  p = p & (~0 << min(m1, m3));
  p = p & (~(~0 << (max(m1, m3) - 1))); // 屏蔽两端
  if (__builtin_popcount(p) % 2 != 0) {
    sign = -1;
  }
  return sign;
}

double addTBO1(int type, int m1, int m3, int w, int lr) { //当每一位都相同时, 对一种粒子的缩并情况求和.
  double add = 0;
  int n_nucleon = __builtin_popcount(w);
  int n;
  if (lr == 1) { // 两种情况 13 left 和 24 right
    for (int i = 0; i < n_nucleon; i++) {
      n = __builtin_ctz(w) + 1;
      add +=
          TBME(type, m1, n, m3, n); // 两个 n 不产生符号. m1 m3 的符号在外面考虑
      w = w & (~0 << n);
    }
  } else if (lr == 2) {
    for (int i = 0; i < n_nucleon; i++) {
      n = __builtin_ctz(w) + 1;
      add += TBME(type, n, m1, n, m3);
      w = w & (~0 << n);
    }
  }
  return add;
}

// Two-Body Operator, act on different kinds. e.g. ΛΝ-ΛΝ
double TBOME(int type, int p1, int q1, int p2, int q2) {
  double me = 0;
  int m1, m2, m3, m4;
  if (difcount(p1, p2) == 2) { //  difcount == 2 意味着 p 只有一种缩并可能
    setTBO2(&m1, &m3, p1, p2);
    int sign13 = signTBO2(m1, m3, p1);
    if (difcount(q1, q2) == 2) {
      setTBO2(&m2, &m4, q1, q2);
      int sign24 = signTBO2(m2, m4, q1);
      me += sign13 * sign24 * TBME(type, m1, m2, m3, m4);
    } else if (difcount(q1, q2) == 0) {
      me += sign13 * addTBO1(type, m1, m3, q1, 1);
    }
  } else if (difcount(p1, p2) == 0) { // 对每一个缩并求和
    if (difcount(q1, q2) == 2) {
      setTBO2(&m2, &m4, q1, q2);
      int sign24 = signTBO2(m2, m4, q1);
      me += sign24 * addTBO1(type, m2, m4, p1, 2);
    } else if (difcount(q1, q2) == 0) {
      int n_nucleon = __builtin_popcount(p1);
      int n, w = p1;
      for (int i = 0; i < n_nucleon; i++) {
        n = __builtin_ctz(w) + 1;
        me += addTBO1(type, n, n, q1, 1);
        w = w & (~0 << n);
      }
    }
  }
  return me;
}

// Two-Body Operator, act on three kinds of particles. e.g. ΛΝ-ΣΝ. no one-boy operator dexchange sign
double TBOME2(int type, int p1, int q1, int p2, int q2) {
  double me = 0;
  int m1, m2, m3, m4;
  if (difcount(p1, p2) == 2) { //  difcount == 2 意味着 p 只有一种缩并可能
    setTBO2(&m1, &m3, p1, p2);
    int sign13 = signTBO2(m1, m3, p1);
    if (difcount(q1, q2) == 2) {
      setTBO2(&m2, &m4, q1, q2);
//      int sign24 = signTBO2(m2, m4, q1);
      me += sign13 * TBME(type, m1, m2, m3, m4);
    } else if (difcount(q1, q2) == 0) {
      me += sign13 * addTBO1(type, m1, m3, q1, 1);
    }
  } else if (difcount(p1, p2) == 0) { // 对每一个缩并求和
    if (difcount(q1, q2) == 2) {
      setTBO2(&m2, &m4, q1, q2);
//      int sign24 = signTBO2(m2, m4, q1);
      me += addTBO1(type, m2, m4, p1, 2);
    } else if (difcount(q1, q2) == 0) {
      int n_nucleon = __builtin_popcount(p1);
      int n, w = p1;
      for (int i = 0; i < n_nucleon; i++) {
        n = __builtin_ctz(w) + 1;
        me += addTBO1(type, n, n, q1, 1);
        w = w & (~0 << n);
      }
    }
  }
  return me;
}

//  two-body operator acting on the same kind
double TBOMEnn(int type, int p, int q) {
  double me = 0;
  int m1 = 0, m2 = 0, m3 = 0, m4 = 0, s;
  int sign12 = 1;
  int sign34 = 1;
  int sign = 1;
  if (difcount(p, q) == 4) {
    s = (p ^ q);
    int n_nucleon = __builtin_popcount(s);
    int n;
    for (int i = 0; i < n_nucleon; i++) { // 从低位开始, 逐个讨论缩并情况
      n = __builtin_ctz(s) + 1;
      if ((p >> (n - 1)) % 2 == 0) {
        if (m3 == 0) {
          m3 = n;
        } else {
          m4 = n;
        }
      } else {
        if (m1 == 0) {
          m1 = n;
        } else {
          m2 = n;
        }
      }
      s = s & (~0 << n);
    }
    sign12 = signTBO2(m1, m2, p);  // 反对易关系产生的符号
    sign34 = signTBO2(m3, m4, q);
    me += sign12 * sign34 * TBME(type, m1, m2, m3, m4);
  } else if (difcount(p, q) == 2) {
    s = (p ^ q);                           //不同位
    int w = p & q;                         //相同位
    int n_nucleon = __builtin_popcount(w); // 对缩并掉的粒子求和, 求和次数 == 可能的缩并数
    int n;
    int lead = leaddifbit(p, q);
    int least = leastdifbit(p, q);        // difcount == 2 时这两者异号.
    for (int i = 0; i < n_nucleon; i++) { // 从低位开始
      n = __builtin_ctz(w) + 1;
      if (n < abs(least)) {  // 缩并的位置小于 least, 这里画个图就清楚了
        m1 = m3 = n;
        if (lead < 0) {
          m4 = -lead;
          m2 = abs(least);
        } else {
          m2 = lead;
          m4 = abs(least);
        }
      } else if (n > abs(lead)) { //缩并位置大于 lead
        m2 = m4 = n;
        if (lead < 0) {
          m3 = -lead;
          m1 = abs(least);
        } else {
          m1 = lead;
          m3 = abs(least);
        }
      } else { // 夹在中间
        if (lead < 0) {
          m2 = m3 = n;
          m4 = -lead;
          m1 = least;
        } else {
          m1 = m4 = n;
          m2 = lead;
          m3 = -least;
        }
      }
      sign = signTBO2(m1,m2,p) * signTBO2(m3,m4,q);
      me += sign * TBME(type, m1, m2, m3, m4);
      w = w & (~0 << n);
    }
  } else if (difcount(p, q) == 0) { // 无符号
    int n_nucleon = __builtin_popcount(p);
    int n, w = p;
    int id[n_nucleon]; // 存下标(不用+1)
    for (int i = 0; i < n_nucleon; i++) {
      n = __builtin_ctz(w) + 1;
      id[i] = n;
      w = w & (~0 << n);
    }
    for (int i = 0; i < n_nucleon; i++) {
      for (int j = i; j < n_nucleon; j++) {
        me += TBME(type, id[i], id[j], id[i], id[j]);
      }
    }
  }
  return me;
}

double vJplus(int type, int cfgJplus, int cfg, double eigen_v) {
    int mpl = leastdifbit(cfgJplus,cfg); // 升之后, 下标变小
    int m = -leaddifbit(cfgJplus,cfg);
    if (mpl > 0 && (m - mpl) == 1 && (mfromA(type,m) + 1) <= jfromA(type,m)) {
            return eigen_v * c_Jplus(type,m);
          }
    else
        return 0;
}

double vJminus(int type, int cfgJplus, int cfg, double eigen_vJplus) {
    int mpl = leastdifbit(cfgJplus,cfg); // 升之后, 下标变小
    int m = -leaddifbit(cfgJplus,cfg);
            if (mpl > 0 && (m - mpl) == 1 && (mfromA(type,m) + 1) <= jfromA(type,m)) {
                return c_Jminus(type,mpl) * eigen_vJplus;
            }
    else
        return 0;
}

void setJminus(int type, int norbt, int n_proton, double *eigen_v2i, double *eigen_vp,
              int cfgi_p) {
    int cfgmax = power(2, norbt) - power(2, norbt - n_proton);
  for (int p = power(2, n_proton) - 1;
       p <= cfgmax;
       p = nextbitpermut(p)) {
    int mpl, m; // 降之前是 mpl(指标小的, 属于 p), 降之后是 m
    if (difcount(p, cfgi_p) == 2) {
        mpl = leastdifbit(p, cfgi_p);
        m = -leaddifbit(p, cfgi_p);
        if (mpl > 0 && (m - mpl) == 1 && (mfromA(type,m) + 1) <= jfromA(type,m)) {
            *eigen_v2i += c_Jminus(type,mpl) * eigen_vp[p];
        }
    }
  }
}
