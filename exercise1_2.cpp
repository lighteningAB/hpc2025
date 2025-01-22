#include <cstdlib>
#include <iostream>

double pi_4_approx(int n)
{
    double returnable = 0;
    for (int i = 1; i < n + 1; i++)
    {
        int k = i - 1;
        returnable += (-1) ^ k / (2 * k + 1);
    }
    return returnable;
}

double pi_2_approx(int n)
{
    double returnable = 1;
    for (int i = 1; i < n + 1; i++)
    {
        returnable = returnable * (4 * i ^ 2) / (4 * i ^ 2 - 1);
    }
    return returnable;
}

void main(){
    int prev = 100;
    int curr = 1000;
    int n = 1;
    while (!(abs(prev-curr)<0.00001)){
        prev = curr;
        curr = pi_4_approx(n);
        n++;
    }
    std::cout << "The value of pi is " << curr*4 << " using the first approximation, found in "<< n << " iterations" <<std::endl;
    int prev = 100;
    int curr = 1000;
    int n = 1;
    while (!(abs(prev-curr)<0.00001)){
        prev = curr;
        curr = pi_2_approx(n);
        n++;
    }
    std::cout << "The value of pi is " << curr*2 << " using the second approximation, found in "<< n << " iterations" <<std::endl;
}