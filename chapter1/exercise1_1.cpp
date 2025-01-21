#include <iostream>
//recursive function to calculate fibonnaci numbers

int fib(int n, int one = 0, int two = 1){
    if (n==0){
        return two;
    }
    else{
        int helper = two;
        two = two + one;
        one = helper;
        return fib(n-1, one, two);
    }
}
