#include "sp_shell.h"
#define nlines 1288
// 比较 writeTBME 的输出结果.
int main() {
  FILE *fp, *fpr;
  if (!(fp = fopen("../M-scheme/int/A4Z2_45pn.int", "rt"))) {
    printf("Can't find file 1!\n");
    exit(1);
  }
  makeTBME(fp);
  int p, q, r, s, t;
  double me[nlines];
  if (!(fpr = fopen("./resultzz.txt", "rt"))) {
    printf("Can't find file 2!\n");
    exit(1);
  }
  int n = 0;
  while (!feof(fpr) || n < nlines) {
    fscanf(fpr, "%d  %d  %d  %d  %d  %lf\n", &p, &q, &r, &s, &t, &me[n]);
    n++;
  }
  n = 0;
  for (int type = 1; type < n_type + 1; type++) {
    for (int a = 1; a < n_orbtM + 1; a++) {
      for (int b = 1; b < n_orbtM + 1; b++) {
        for (int c = 1; c < n_orbtM + 1; c++) {
          for (int d = 1; d < n_orbtM + 1; d++) {
            if (TBME(type, a, b, c, d) != 0) {
              if (fabs(me[n] - TBME(type, a, b, c, d)) > 1e-6) {
                printf("error at %d %d %d %d %d %f %f\n", type, a, b, c, d,
                       TBME(type, a, b, c, d), me[n]);
              }
              n++;
            }
          }
        }
      }
    }
  }
  fclose(fp);
  fclose(fpr);
  return 0;
}
