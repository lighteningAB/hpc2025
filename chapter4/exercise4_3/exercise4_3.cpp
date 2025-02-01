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
    void F77NAME(daxpy)(const int &n , const double &alpha, const double *x, const int &incx, double *y, const int &incy);
    void F77NAME(dger)(const int &n, const int&m, const double&alpha, const double* x, const int &incx, const double*y, 
                       const int &incy, double* a, const int &LDA)
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
    double* a = new double[n*n];
    squarematgen(n, a);
    //generated m matrix
    double* x = new double[n];
    vecgen(n, x);
    //generated x matrix
    inplacesym(a, n);
    //in place conversion of m to A
    double* b = new double[n];
    matvecmul(a, x, b, n); //b = Ax
    //fill in vector b with vector product of A and x
    double* x0 = new double[n];
    for (int i = 0; i<=n; i++){
        x0[i] = 0;
    }
    //generate x0 of zeroes

    double *r = new double[n];
    F77NAME(daxpy)(n, 1, b, 1, r, 1); //r = b
    //use vector sum to make vector r equal to b
    double * c = new double[n];
    //c is helper for us to get A x0
    matvecmul(a, x0, c, n); //Ax0
    F77NAME(daxpy)(n, -1, c, 1, r, 1); //r = r0
    //create r0 where r0 = b-Ax0
    double * p = new double[n];
    F77NAME(daxpy)(n, 1, r, 1, p, 1); //p0 = r0
    //cloning r0
    int k = 0;
    while (0==0 ){ //replace with abs r_k+1 < eta break condition
        double * a_k = new double [n*n];
        F77NAME(dger) (n, n, 1.0, r, 1, r, 1, a_k, 1); //top part of fraction
        //a_k = r_k^t * r_k / p_k^t*A*p_k
        //x_k = x_k + 
    }
}