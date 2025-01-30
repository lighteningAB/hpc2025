//conjugate gradient method
#include <iostream>
#include <random>

#define F77NAME(x) x##_

extern "C"
{
    void F77NAME(dgemm)(const char &transa, const char &transb, const int &m, const int &n, const int &k,
                        const double &alpha, const double *a, const int &lda, const double *b, const int &ldb,
                        const double &beta, double *c, const int &ldc);
    void F77NAME(dgemv)(const char &transa, const int &m, const int &n, const double &alpha, const double *a, const int &lda,
                        const double *x, const int &incx, const double &beta, double *y, const int &incy);
}
//1
void squarematgen(int n, double *a){
    for (int i = 0; i<n*n; i++){
        *(a+i) = (double)rand()*2/RAND_MAX - 1;
    }
}

void vecgen(int n, double* a){
    for (int i = 0; i<n; i++){
        *(a+i) = (double)rand()*2/RAND_MAX - 1;
    }
}

//2
void inplacesym(double* a, int n){
    double * b = new double[n*n];
    F77NAME(dgemm)('T', 'N', n, n, n, 1.0, a, n, a, n, 0.0, b, n);
    a = b;
    delete[] b;
}

//3 
void matvecmul(double* a, double* b, double* c, int n){
    F77NAME(dgemv)('T', n, n, 1.0, a, n, b, 1, 0.0, c, 1);
}

//4
int main(){
    int n = 5;
    double* m = new double[n*n];
    squarematgen(n, m);
    double* x = new double[n];
    vecgen(n, x);
    inplacesym(m, n);
    double* b = new double[n];
    matvecmul(m, x, b, n);

    //do later
}