#include <iostream>
//recursive function to calculate fibonnaci numbers

int fib(int n, int prev = 0, int curr = 1){
    if (n<=0){
        return 0;
    }
    else if (n==1){
        return curr;
    }
    else{
        int helper = curr;
        curr = prev + curr;
        prev = helper;
        return fib(n-1, prev, curr);
    }
}

int main() {
    int n;
    std::cout << "Enter a number: ";
    std::cin >> n;
    std::cout << "Fibonacci(" << n << ") = " << fib(n) << std::endl;
    return 0;
}
