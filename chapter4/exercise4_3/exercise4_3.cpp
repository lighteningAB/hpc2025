// conjugate gradient method
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
// 1
void squarematgen(int n, double *a)
{
    for (int i = 0; i < n * n; i++)
    {
        *(a + i) = (double)rand() * 2 / RAND_MAX - 1;
    }
}

void vecgen(int n, double *a)
{
    for (int i = 0; i < n; i++)
    {
        *(a + i) = (double)rand() * 2 / RAND_MAX - 1;
    }
}

// 2
void inplacesym(double *a, int n)
{
    double *b = new double[n * n];
    F77NAME(dgemm)('T', 'N', n, n, n, 1.0, a, n, a, n, 0.0, b, n);
    a = b;
    delete[] b;
}

// 3
void matvecmul(double *a, double *b, double *c, int n)
{
    F77NAME(dgemv)('T', n, n, 1.0, a, n, b, 1, 0.0, c, 1);
}

double vecmult(int n, double *x, double *y)
{
    return F77NAME(ddot)(n, x, 1, y, 1);
}
// 4
int main()
{
    int n = 5;
    double *a = new double[n * n];
    squarematgen(n, a);
    // generated m matrix
    double *x = new double[n];
    vecgen(n, x);
    // generated x matrix
    inplacesym(a, n);
    // in place conversion of m to A
    double *b = new double[n];
    matvecmul(a, x, b, n); // b = Ax
    // fill in vector b with vector product of A and x
    double *x0 = new double[n];
    for (int i = 0; i <= n; i++)
    {
        x0[i] = 0;
    }
    // generate x0 of zeroes

    double *r = new double[n];
    F77NAME(daxpy)(n, 1, b, 1, r, 1); // r = b
    // use vector sum to make vector r equal to b
    double *c = new double[n];
    // c is helper for us to get A x0
    matvecmul(a, x0, c, n);            // Ax0
    F77NAME(daxpy)(n, -1, c, 1, r, 1); // r = r0
    // create r0 where r0 = b-Ax0
    double *p = new double[n];
    F77NAME(daxpy)(n, 1, r, 1, p, 1); // p0 = r0
    // cloning r0
    int k = 0;
    while (0 == 0)
    { // replace with abs r_k+1 < eta break condition
        double *bottomhelpervec = new double[n];
        matvecmul(a, p, bottomhelpervec, n);
        double a_k = vecmult(n, r, r) / vecmult(n, p, bottomhelpervec);
        // a_k = r_k^t * r_k / p_k^t*A*p_k
        F77NAME(daxpy)(n, a_k, p, 1, x0, 1);
        // x_k+1 = x_k + a_k*p_k
        double *helpervec_2 = new double[n];
        matvecmul(a, p, helpervec_2, n);
        double *r_prev = new double[n];
        for (int i = 0; i < n; i++)
        {
            r_prev[i] = r[i];
        }
        F77NAME(daxpy)(n, a_k * -1, helpervec_2, 1, r, 1);
        // r_k+1 = r_k - a_k*a*p_k
        double stop = F77NAME(dnrm2)(n, r, 1);
        if (stop<0.01){
            break;
        }
        double B_k = vecmult(n, r, r) / vecmult(n, r_prev, r_prev);
        // B_k = r_k+1^t * r_k+1/ r_k^t * r_k
        double *p_prev = new double[n];
        for (int i = 0; i < n; i++)
        {
            p_prev[i] = p[i];
        }
        F77NAME(daxpy)(n, a_k, helpervec_2, 1, p, 1); // p_k+1 = r_k+1
        F77NAME(daxpy)(n, B_k, p_prev, 1, p, 1);
        // p_k+1 = r_k+1 + B_k*p_k
        k += 1;
        delete[] bottomhelpervec;
        delete[] helpervec_2;
        delete[] r_prev;
        delete[] p_prev;
    }

    for (int i = 1; i < n; i++)
    {
        std::cout << x0[i] << " ";
    }
    std::cout << std::endl;

    delete[] a;
    delete[] x;
    delete[] b;
    delete[] x0;
    delete[] r;
    delete[] c;
    delete[] p;
}