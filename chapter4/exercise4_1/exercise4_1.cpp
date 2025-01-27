// generate an array of N random doubles between 0 and 1
#include <random>
#include <iostream>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#define F77NAME(x) x##_

using namespace boost;
using namespace boost::program_options;

extern "C"
{
    double F77NAME(dasum)(const int &n,
                          const double *x, const int &incx);
    double F77NAME(dnrm2)(const int &n,
                          const double *x, const int &incx);
    double F77NAME(idamax)(const int &n,
                            const double *x, const int &incx);
}

double calcs(double *a, int n)
{
    double returnable = 0.0;
    return returnable = F77NAME(dasum)(n, a, 1);
}

double manualCalcS(double *a, int n)
{
    double returnable = 0.0;
    for (int i = 0; i < n; i++)
    {
        returnable += fabs(a[i]);
    }
    return returnable;
}

double calcd(double *a, int n)
{
    double returnable = 0.0;
    return returnable = F77NAME(dnrm2)(n, a, 1);
}

double manualCalcD(double *a, int n)
{
    double returnable = 0.0;
    for (int i = 0; i < n; i++)
    {
        returnable += a[i] * a[i];
    }
    return returnable = sqrt(returnable);
}

double calcm(double *a, int n)
{
    double returnable = 0.0;
    return returnable = F77NAME(idamax)(n, a, 1);
}

double manualCalcM(double *a, int n)
{
    double returnable = 0.0;
    for (int i = 0; i < n; i++)
    {
        if (fabs(a[i]) > returnable)
        {
            returnable = fabs(a[i]);
        }
    }
    return returnable;
}

int main(int argc, char *argv[])
{

    options_description desc{"Options"};
    desc.add_options()("size", value<int>()->default_value(5));
    variables_map vm;
    try
    {
        store(parse_command_line(argc, argv, desc), vm);
        notify(vm);
    }
    catch (const error &ex)
    {
        std::cerr << ex.what() << std::endl;
        return 1;
    }
    int n = vm["size"].as<int>();
    double *array = new double[n];
    for (int i = 0; i < n; i++)
    {
        array[i] = (double)rand() / RAND_MAX;
    }
    for (int j = 0; j < n; j++)
    {
        std::cout << array[j] << std::endl;
    }
    double s = calcs(array, n);
    double d = calcd(array, n);
    double m = calcm(array, n);
    double s2 = manualCalcS(array, n);
    double d2 = manualCalcD(array, n);
    double m2 = manualCalcM(array, n);

    std::cout << "s: " << s << " s2: " << s2 << std::endl;
    std::cout << "d: " << d << " d2: " << d2 << std::endl;
    std::cout << "m: " << m << " m2: " << m2 << std::endl;

    delete[] array;
    return 0;
}