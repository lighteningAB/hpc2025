#include <iostream>
#include <cmath>
#include <cstdlib>
#include <vector>

std::vector<double> vecadd(std::vector<double> a, std::vector<double> b)
{
    // allocating dynamic memory
    std::vector<double> returnable;
    for (int i = 0; i < a.size(); i++)
    {
        returnable.push_back(a.at(i) + b.at(i));
    }
    return returnable;
}

std::vector<double> dot(std::vector<double> a, std::vector<double> b)
{
    std::vector<double> returnable;
    for (int i = 0; i < a.size(); i++)
    {
        returnable.push_back(a.at(i) * b.at(i));
    }
    return returnable;
}

int main()
{
    // test for vecadd
    return 0;
}