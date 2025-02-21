#include <iostream>
#include "fibonacci.h"

int main() {
    int n;
    std::cout << "Enter an integer: ";
    std::cin >> n;
    std::cout << "Fibonacci(" << n << ") = " << fibonacci(n) << std::endl;
    return 0;
}