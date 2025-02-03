#ifndef EXERCISE4_3_H // include guard
#define EXERCISE4_3_H

#include <iostream>
#include <random>

#define F77NAME(x) x##_
extern "C"
{
    void F77NAME(dgemm)(const char &transa, const char &transb, const int &m, const int &n, const int &k,
                        const double &alpha, const double *a, const int &lda, const double *b, const int &ldb,
                        const double &beta, double *c, const int &ldc);
    void F77NAME(dgemv)(const char &transa, const int &m, const int &n, const double &alpha, const double *a, const int &lda,
                        const double *x, const int &incx, const double &beta, double *y, const int &incy);
    void F77NAME(daxpy)(const int &n, const double &alpha, const double *x, const int &incx, double *y, const int &incy);
    void F77NAME(dger)(const int &n, const int &m, const double &alpha, const double *x, const int &incx, const double *y,
                       const int &incy, double *a, const int &LDA);
    double F77NAME(ddot)(const int &n, const double *x, const int &incx, const double *y, const int &incy);
    double F77NAME(dnrm2)(const int &n, const double *x, const int &incx);
}

void matvecmul(double *a, double *b, double *c, int n);
double vecmult(int n, double *x, double *y);

#endif