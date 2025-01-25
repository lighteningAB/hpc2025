//insertion sort

#include <iostream>
//p is index of last sorted element
//i is where the insertion element is
double* insert(double a[], double x, int p){
    int i = p+1;
    while (i >= 0 && a[i-1] > x){
        if (i == 0){
            a[i] = x;
            return a;
        }
        a[i] = a[i-1];
        i--;
    }
    a[i] = x;
    return a;
}

double* insertion_sort(double a[], int n){
    for (int i = 1; i <= n-1; i++){
        a = insert(a, a[i], i-1);
    }
    return a;
}

int main(){
    double a[] = {3, 2, 1, 4, 5};
    int n = 5;
    double* b = insertion_sort(a, n);
    for (int i = 0; i < n; i++){
        std::cout << b[i] << " ";
    }
    return 0;
}