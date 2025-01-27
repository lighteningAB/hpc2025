//generate an array of N random doubles between 0 and 1
#include <random>
#include <iostream>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

using namespace boost;
using namespace boost::program_options;

int main(int argc, char *argv[]){

    options_description desc{"Options"};
    desc.add_options()
        ("size", value<int>()->default_value(5));
    variables_map vm;
    try{
        store(parse_command_line(argc, argv, desc), vm);
        notify(vm);
    }
    int n = 10;
    double* array = new double[n];
    for (int i = 0; i < n; i++)
    {
        array[i] = (double)rand() / RAND_MAX;
    }
}