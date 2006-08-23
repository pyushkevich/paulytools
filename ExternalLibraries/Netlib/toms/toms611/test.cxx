#include <cstdio>

extern "C" {
  void sumsl_(
    int &n,
    double *d,
    double *x,
    void (*calcf)(int &, double *, int &, double &, int *, double *, void *),
    void (*calcg)(int &, double *, int &, double *, int *, double *, void *),
    int *iv,
    int &liv,
    int &lv,
    double *v,
    int *uiparm,
    double *urparm,
    void (*ufparm)());
}

// Define the function (something simple)
void myf(int &n, double *x, int &nf, double &f, int *, double *, void *)
{
  f = x[0] * x[0] + x[1] * x[1];
}

void mydf(int &n, double *x, int &nf, double *g, int *, double *, void *)
{
  g[0] = 2 * x[0];
  g[1] = 2 * x[1];
}

int main()
{
  int n = 2;
  double dScale[] = {1.0, 1.0};
  double x[] = {0.565, 0.245};
  int liv = 60;
  int lv = 71+n*(n+15)/2;
  int *iv = new int[liv];
  double *v = new double[lv];

  sumsl_(n, dScale, x, &myf, &mydf, iv, liv, lv, v, NULL, NULL, NULL);
  return 0;
}
     

