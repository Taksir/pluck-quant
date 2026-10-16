// Minimal sparse-probe EM (SQUAREM) over equivalence classes. Input binary files written by dump.py.
// Usage: ./em DIR  -> prints iterations, seconds; writes DIR/X.bin
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
static double *rd(const char *d, const char *f, long *n, int isint) {
  char p[512]; snprintf(p, sizeof p, "%s/%s", d, f); FILE *fp = fopen(p, "rb"); if (!fp) { perror(p); exit(1); }
  fseek(fp, 0, SEEK_END); long sz = ftell(fp); fseek(fp, 0, SEEK_SET);
  void *buf = malloc(sz); if (fread(buf, 1, sz, fp) != (size_t)sz) exit(2); fclose(fp);
  *n = sz / 8; return (double *)buf;
}
static long nec, ne, ntx, ng; static long *eoff, *et, *gof; static double *ea, *cec, *A;
static void step(const double *X, double *out) {
  memset(out, 0, ntx * sizeof(double));
  for (long k = 0; k < nec; k++) {
    if (cec[k] <= 0) continue;
    double mu = 0; for (long q = eoff[k]; q < eoff[k + 1]; q++) mu += ea[q] * X[et[q]];
    if (mu <= 0) continue;
    double f = cec[k] / mu; for (long q = eoff[k]; q < eoff[k + 1]; q++) out[et[q]] += ea[q] * X[et[q]] * f;
  }
  for (long t = 0; t < ntx; t++) out[t] = A[t] > 0 ? out[t] / A[t] : 0;
}
int main(int argc, char **argv) {
  long n; const char *d = argv[1]; double tol = argc > 2 ? atof(argv[2]) : 1e-3;
  eoff = (long *)rd(d, "eoff.bin", &n, 1); nec = n - 1;
  et = (long *)rd(d, "et.bin", &ne, 1); ea = rd(d, "ea.bin", &n, 0); cec = rd(d, "cec.bin", &n, 0);
  A = rd(d, "A.bin", &ntx, 0); gof = (long *)rd(d, "gof.bin", &n, 1);
  ng = 0; for (long t = 0; t < ntx; t++) if (gof[t] + 1 > ng) ng = gof[t] + 1;
  double *X = malloc(ntx * 8), *X1 = malloc(ntx * 8), *X2 = malloc(ntx * 8), *Xn = malloc(ntx * 8), *g = calloc(ng, 8), *gp = calloc(ng, 8);
  struct timespec t0, t1; clock_gettime(CLOCK_MONOTONIC, &t0);
  for (long t = 0; t < ntx; t++) { X[t] = A[t] > 0; gp[gof[t]] += X[t]; }
  int it;
  for (it = 0; it < 3000; it++) {
    step(X, X1); step(X1, X2);
    double rr = 0, vv = 0;
    for (long t = 0; t < ntx; t++) { double r = X1[t] - X[t], v = X2[t] - X1[t] - r; rr += r * r; vv += v * v; }
    double a = vv > 0 ? -sqrt(rr) / sqrt(vv) : -1; if (a > -1) a = -1;
    for (long t = 0; t < ntx; t++) { double r = X1[t] - X[t], v = X2[t] - X1[t] - r; double x = X[t] - 2 * a * r + a * a * v; X2[t] = x > 0 ? x : 0; }
    step(X2, Xn);
    memset(g, 0, ng * 8); for (long t = 0; t < ntx; t++) g[gof[t]] += Xn[t];
    double dmax = 0; for (long j = 0; j < ng; j++) if (g[j] > 10) { double dd = fabs(g[j] - gp[j]) / g[j]; if (dd > dmax) dmax = dd; }
    memcpy(X, Xn, ntx * 8); memcpy(gp, g, ng * 8);
    if (dmax < tol) break;
  }
  clock_gettime(CLOCK_MONOTONIC, &t1);
  printf("iters %d em_steps %d seconds %.3f nec %ld entries %ld ntx %ld\n", it + 1, 3 * (it + 1), (t1.tv_sec - t0.tv_sec) + 1e-9 * (t1.tv_nsec - t0.tv_nsec), nec, ne, ntx);
  char p[512]; snprintf(p, sizeof p, "%s/X.bin", d); FILE *fo = fopen(p, "wb"); fwrite(X, 8, ntx, fo); fclose(fo);
  return 0;
}
