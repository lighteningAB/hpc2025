#ifndef EXERCISE4_4_H //include guard
#define EXERCISE4_4_H

#include <math.h>
#include <iostream>

#define F77NAME(x) x##_
// helmholtz
extern "C"
{
    void F77NAME(daxpy)(const int &n, const double &alpha, const double *x, const int &incx, double *y, const int &incy);
    double F77NAME(ddot)(const int &n, const double *x, const int &incx, const double *y, const int &incy);
    double F77NAME(dsymv)(const char &uplo, const int &n, const double &alpha, const double *a, const int &lda,
                          const double *x, const int &incx, const double &beta, double *y, const int &incy);
    double F77NAME(dnrm2)(const int &n, const double *x, const int &incx);
}

void matvecmul(double *a, double *b, double *c, int n);
double vecmult(int n, double *x, double *y);\
void symmetricRowMaj(double a, double b, int n, double *returnable);
void symmetricColMaj(double a, double b, int n, double *returnable);
void forcfunc(double *a, int n);
void forcfunc2(double *a, int n);
void conjgradsolve(int n, double *a, double *b, double *x0);

#endif