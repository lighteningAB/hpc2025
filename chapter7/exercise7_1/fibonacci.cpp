#include "fibonacci.h"

int fibonacci(int n, int prev, int curr) {
    if (n <= 0) {
        return 0;
    }
    else if (n == 1) {
        return curr;
    }
    else {
        int helper = curr;
        curr = prev + curr;
        prev = helper;
        return fibonacci(n-1, prev, curr);
    }
}