// multiply two nxn square matrices where n has to be a power of 2 greater than 2x2
// do with blas and not blas
// use blas to compare results

#include <iostream>
#include <random>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#define F77NAME(x) x##_

using namespace boost;
using namespace boost::program_options;

extern "C"
{
    void F77NAME(dgemm)(const char &transa, const char &transb, const int &m, const int &n, const int &k,
                        const double &alpha, const double *a, const int &lda, const double *b, const int &ldb,
                        const double &beta, double *c, const int &ldc);
    void F77NAME(daxpy)(const int &n, const double &alpha, const double *x, const int &incx, double *y,
                        const int &incy);
}

// blas array multiplication
void multarrays(double *a, double *b, double *c, int n)
{
    F77NAME(dgemm)
    ('T', 'T', n, n, n, 1.0, a, n, b, n, 0.0, c, n);
}

// manual array multiplication
void manualmultarrays(double *a, double *b, double *c, int n)
{
    // accessing each row of a
    for (int i = 0; i < n; i++)
    {
        // accessing each column of b
        for (int j = 0; j < n; j++)
        {
            double sum = 0;
            // going down b's column
            for (int k = 0; k < n; k++)
            {
                sum += a[i * n + k] * b[k * n + j];
            }
            c[i * n + j] = sum;
        }
    }
}

void subtractArrays(double *a, double *b, double *c, int n)
{
    F77NAME(daxpy)
    (n, -1.0, b, 1, a, 1);
}

int main(int argc, char *argv[])
{
    options_description desc{"Options"};
    desc.add_options()("size,n", value<int>()->default_value(3), "Size of the array 2^N");

    variables_map vm;
    store(parse_command_line(argc, argv, desc), vm);
    notify(vm);

    int n = vm["size"].as<int>();
    int size = pow(2, n);
    double *a = new double[size * size];
    double *b = new double[size * size];
    double *c = new double[size * size];
    double *c2 = new double[size * size];
    double *c3 = new double[size * size];

    double max = 1.0;
    for (int i = 0; i < size * size; i++)
    {
        a[i] = (double)rand() / RAND_MAX * max;
        b[i] = (double)rand() / RAND_MAX * max;
    }

    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            std::cout << a[i * size + j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "_________" << std::endl;

    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            std::cout << b[i * size + j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "_________" << std::endl;

    multarrays(a, b, c, size);

    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            std::cout << c[i * size + j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "_________" << std::endl;

    manualmultarrays(a, b, c2, size);

    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            std::cout << c2[i * size + j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "_________" << std::endl;

    subtractArrays(c, c2, c3, n);

    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            std::cout << c3[i * size + j] << " ";
        }
        std::cout << std::endl;
    }

    delete[] a;
    delete[] b;
    delete[] c;

    return 0;
}