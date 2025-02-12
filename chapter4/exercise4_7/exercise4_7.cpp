#include <iostream>
#include <random>
#include <complex>

#define F77NAME(x) x##_

extern "C"
{
    void F77NAME(zgemm)(const char &transa, const char &transb, const int &m, const int &n, const int &k,
                        const std::complex<double> &alpha, const std::complex<double> *a, const int &lda, const std::complex<double> *b,
                        const int &ldb, const std::complex<double> &beta, const std::complex<double> *c, const int &ldc);
}

// a and b need to be size n*n
void genAB(int n, std::complex<double> *a, std::complex<double> *b)
{
    double real;
    double imag;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            real = ((double)rand() / RAND_MAX) * (20.0) - 10.0;
            imag = ((double)rand() / RAND_MAX) * (20.0) - 10.0;
            a[i * n + j] = std::complex<double>(real, -imag);
            b[i + n * j] = std::complex<double>(real, imag);
        }
    }
}

void multCompSquare(int n, std::complex<double> *a, std::complex<double> *b, std::complex<double> *c)
{
    F77NAME(zgemm)('N', 'N', n, n, n, 1, a, n, b, n, 1, c, n);
}

bool hermitianVerify(int n, std::complex<double> *a)
{
    bool returnable = true;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            a[i * n + j] == std::conj(a[i + j * n]) ? returnable = returnable : returnable = false;
        }
    }
    return returnable;
}

int main()
{
    int n = 2;
    std::complex<double> *a = new std::complex<double>[n * n];
    std::complex<double> *b = new std::complex<double>[n * n];
    genAB(n, a, b);
    /*
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            std::cout << a[i * n + j] << " ";
        }
        std::cout << std::endl;
    }
    */

    /*
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            std::cout << b[i * n + j] << " ";
        }
        std::cout << std::endl;
    }
    */
    std::complex<double> *c = new std::complex<double>[n * n];
    multCompSquare(n, a, b, c);
    ///*
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            std::cout << c[i * n + j] << " ";
        }
        std::cout << std::endl;
    }
    //*/
    std::cout << ((hermitianVerify(n, c) == 1)? "true":"false") << std::endl;
    delete[] a;
    delete[] b;
    delete[] c;    
}