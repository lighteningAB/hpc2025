#include <complex>
#include <cmath>
#include <iostream>

#define F77NAME(x) x##_

extern "C"
{
    void F77NAME(zsymv)(const char &uplo, const int &n, const std::complex<double> &alpha, const std::complex<double> *a, const int &lda,
                        const std::complex<double> *x, const int &incx, const std::complex<double> &beta, std::complex<double> *y, const int &incy);
    double F77NAME(dznrm2)(const int &n, const std::complex<double> *a, const int &incx);
}

void genW(std::complex<double> *a, int n)
{
    const std::complex<double> i(0.0, 1.0);
    for (int k = 0; k < n; k++)
    {
        for (int j = 0; j < n; j++)
        {
            a[k * n + j] = pow(exp((-2.0 / n * M_PI * i)), j * k) * 1.0 / sqrt(double(n));
        }
    }
}

void genX(int n, std::complex<double> *x)
{
    double h = 1.0 / (double(n) - 1);
    for (int i = 0; i < n; i++)
    {
        x[i] = std::complex<double>(cos(2 * M_PI * i * h) + cos(6 * M_PI * i * h) + 1.0, 0);
    }
}

void multCompSym(int n, std::complex<double> *a, std::complex<double> *b, std::complex<double> *c)
{
    std::complex<double> alpha = 1.0;
    std::complex<double> beta = 0.0;
    F77NAME(zsymv)('L', n, alpha, a, n, b, 1, beta, c, 1);
}

int main()
{
    int n = 16;
    std::complex<double> *w = new std::complex<double>[n * n];
    genW(w, n);
    /*
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            std::cout << w[i * n + j] << " ";
        }
        std::cout << std::endl;
    }
    */
    std::complex<double> *x = new std::complex<double>[n];
    std::complex<double> *bigX = new std::complex<double>[n];
    genX(n, x);
    multCompSym(n, w, x, bigX);
    for (int i = 0; i < n / 2; i++)
    {
        std::cout << bigX[i] << " ";
    }
    std::cout << std::endl;
    for (int i = 0; i < 3; i++)
    {
        std::copy(bigX, bigX + n, x);
        multCompSym(n, w, x, bigX);
    }
    std::complex<double> *y = new std::complex<double>[n];
    double h = 1.0 / (double(n) - 1);
    for (int i = 0; i < n; i++)
    {
        y[i] = bigX[i] - ((2 * M_PI * i * h) + cos(6 * M_PI * i * h) + 1.0);
    }
    double diff = F77NAME(dznrm2)(n, y, 1);
    std::cout << diff << std::endl;

    delete[] w;
    delete[] x;
    delete[] bigX;
    delete[] y;
}