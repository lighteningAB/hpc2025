#include <cstdlib>
#include <iostream>
#include <cmath>
#include <boost/timer/timer.hpp>

double pi_4_approx(int n)
{
    double returnable = 0.0;
    for (int i = 1; i < n + 1; i++)
    {
        int k = i - 1;
        returnable += pow(-1,k) / (2.0 * k + 1.0);
    }
    return returnable;
}

double pi_2_approx(int n)
{
    double returnable = 1.0;
    for (int i = 1; i < n + 1; i++)
    {
        returnable = returnable * (4.0 * pow(i,2)) / (4.0 * pow(i,2) - 1.0);
    }
    return returnable;
}

int main(){
    boost::timer::auto_cpu_timer t;
    double prev = 100.0;
    double curr = 1000.0;
    int n = 1;
    while (fabs(prev-curr)>0.00001){
        prev = curr;
        curr = pi_4_approx(n);
        n+=1;
    }
    std::cout << "The value of pi is " << curr*4 << " using the first approximation, found in "<< n << " iterations" <<std::endl;
    prev = 100;
    curr = 1000;
    n = 1;
    while (fabs(prev-curr)>0.00001){
        prev = curr;
        curr = pi_2_approx(n);
        n+=1;
    }
    std::cout << "The value of pi is " << curr*2 << " using the second approximation, found in "<< n << " iterations" <<std::endl;
    return 0;
}