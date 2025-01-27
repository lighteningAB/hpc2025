// insertion sort takes command link args that specify random array size, range of values, and ascending or descending
// i.e. ./insertion-sort --size=100 --min=0 --max=100 --descending

#include <iostream>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <random>

using namespace boost;
using namespace boost::program_options;
// p is index of last sorted element
// i is where the insertion element is
double *insert(double a[], double x, int p)
{
    int i = p + 1;
    while (i >= 0 && a[i - 1] > x)
    {
        if (i == 0)
        {
            a[i] = x;
            return a;
        }
        a[i] = a[i - 1];
        i--;
    }
    a[i] = x;
    return a;
}

double *insertion_sort(double a[], int n)
{
    for (int i = 1; i <= n - 1; i++)
    {
        a = insert(a, a[i], i - 1);
    }
    return a;
}

int main(int argc, char *argv[])
{
    // default values

    // options
    options_description desc{"Options"};
    desc.add_options()("size", value<int>()->default_value(10), "Size of the array")("min", value<double>()->default_value(0.0), "Minimum value")("max", value<double>()->default_value(5.0), "Maximum value")("ascending", "Sort in ascending order")("descending", "Sort in descending order");

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

    int size = vm["size"].as<int>();
    double min = vm["min"].as<double>();
    double max = vm["max"].as<double>();
    double *a = new double[size];
    for (int i = 0; i < size; i++)
    {
        a[i] = (max - min) * (rand() / (double)RAND_MAX) + min;
    }
    double *b = insertion_sort(a, size);
    if (vm.count("descending"))
    {
        for (int i = 0; i < size / 2; i++)
        {
            double temp = b[i];
            b[i] = b[size - i - 1];
            b[size - i - 1] = temp;
        }
    }
    for (int i = 0; i < size; i++)
    {
        std::cout << b[i] << " ";
    }
    delete[] a;
    return 0;
}