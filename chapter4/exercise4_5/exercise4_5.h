#ifndef EXERCISE4_5_H
#define EXERCISE4_5_h

#include <math.h>
#include <iostream>

#define F77NAME(x) x##_

extern "C"
{
    void F77NAME(daxpy)(const int &n, const double &alpha, const double *x, const int &incx, double *y, const int &incy);
    double F77NAME(ddot)(const int &n, const double *x, const int &incx, const double *y, const int &incy);
    void F77NAME(dgbmv)(const char &trans, const int &m, const int &n, const int &kl, const int &ku, const double &alpha,
                        const double *a, const int &LDA, const double *x, const int &incx, const double &beta, double *y, const int &incy);
    double F77NAME(dnrm2)(const int &n, const double *x, const int &incx);
}

void matvecmulband(double *a, double *b, double *c, int n);
double vecmult(int n, double *x, double *y);
void bandRowMaj(double a, double b, int n, double *returnable);
void forcfunc(double *a, int n);
void conjgradsolve(int n, double *a, double *b, double *x0);

#endif