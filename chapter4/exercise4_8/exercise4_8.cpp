//use blas to operate on submatrices to recursively multiply two square matrices of dimensions nxn

#define F77NAME(x) x##_
#include <iostream>
#include <random>

extern "C"{
    void F77NAME(dgemm)(const char &transa, const char &transb, const int &m, 
                        const int &n, const int &k, const double &alpha, 
                        const double* a, const int &lda, const double* b, 
                        const int &ldb, const double &beta, double* c, const int &ldc);
}

//function has input of two nxn arrays where n is a multiple of 2 and returns the product
double* horizontalcat(double* a, double* b, int n){
    double * returnable = new double [2*n*n];
    int index = 0; 
    //changes which row
    for(int i = 0; i<n; i++){
        //switches from a to b
        for (int s = 0; s<2; s++){
        //iterates through a row
            for (int j = 0; j<n; j++){
                if(s==0){
                    *(returnable+index)=(a[i*n+j]);
                    index++;
                }
                else{
                    *(returnable+index)=(b[i*n+j]);
                    index++;
                }
            }
        }
    }
    return returnable;
}

//concatenate two rectangular arrays of size (2n, n)
double* vertcat(double* a, double* b, int n){
    double * returnable = new double [4*n*n];
    int index = 0;
    for(int i = 0; i<2*n*n; i++){
        *(returnable+index)=(a[i]);
        index++;
    }
    for(int i = 0; i<2*n*n; i++){
        *(returnable+index)=(b[i]);
        index++;
    }
    return returnable;
}

//add two arrays of size (n,n)
double* vecadd(double * a, double * b, int n){
    double * returnable = new double [n*n];
    for(int i = 0; i<n*n; i++){
        *(returnable+i)=*(a+i)+*(b+i);
    }
    return returnable;
}

/**
 * @brief Recursively compute the product of two square arrays A and B
 *
 * @param a The pointer to the first item of the first array A
 * @param b The pointer to the first item of the second array B
 * @param matsize the size of the square arrays A and B
 * @param ld the ld that we use to find the next values in A and B
 * @return The product of a and b
 */
double* blassen(double* a, double* b, int matsize, int ld){
    //base case, multiply 2x2 matrix in blas
    if (matsize == 2){
        double* result = new double[4];
        F77NAME(dgemm)('T', 'T', matsize, matsize, matsize, 1.0, a, ld, b, ld, 1.0, result, matsize);
        return result;
    }
    else{
        int halfsize = matsize/2;
        double * a11b11 = new double[halfsize*halfsize]();
        double * a12b21 = new double[halfsize*halfsize]();
        double * a11b12 = new double[halfsize*halfsize]();
        double * a12b22 = new double[halfsize*halfsize]();
        double * a21b11 = new double[halfsize*halfsize]();
        double * a22b21 = new double[halfsize*halfsize]();
        double * a21b12 = new double[halfsize*halfsize]();
        double * a22b22 = new double[halfsize*halfsize]();

        a11b11 = blassen(a, b, halfsize, ld);
        a12b21 = blassen(a+halfsize, b+halfsize*ld, halfsize, ld);
        a11b12 = blassen(a, b+halfsize, halfsize, ld);
        a12b22 = blassen(a+halfsize, b+halfsize*(ld+1), halfsize, ld);
        a21b11 = blassen(a+halfsize * ld, b, halfsize, ld);
        a22b21 = blassen(a+halfsize*(ld+1), b+halfsize*ld, halfsize, ld);
        a21b12 = blassen(a+halfsize*ld, b+halfsize, halfsize, ld);
        a22b22 = blassen(a+halfsize*(ld+1), b+halfsize*(ld+1), halfsize, ld);

        return vertcat(
            horizontalcat(
                vecadd(a11b11, a12b21, halfsize),
                vecadd(a21b11, a22b21, halfsize),
                halfsize
            ),
            horizontalcat(
                vecadd(a11b12, a12b22, halfsize),
                vecadd(a21b12, a22b22, halfsize),
                halfsize
            ),
            halfsize
        );
    }
}

int main(){
    int n = 16;
    int max = 3;
    double* a = new double [n*n];
    double* b = new double [n*n];
    for (int i = 0; i<n; i++){
        for(int j = 0; j<n; j++){
            a[i*n+j] = rand()/(RAND_MAX/(max+1));
            b[i*n+j] = rand()/(RAND_MAX/(max+1));
        }
    }
    std::cout << "Matrix A:\n";
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            std::cout << a[i * n + j] << " ";
        }
        std::cout << "\n";
    }
    std::cout << "Matrix B:\n";
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            std::cout << b[i * n + j] << " ";
        }
        std::cout << "\n";
    }

    double* c = new double [n*n];
    c = blassen(a, b, n, n);
    std::cout << "Resulting Matrix:\n";
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            std::cout << c[i * n + j] << " ";
        }
        std::cout << "\n";
    }

    double* d = new double [n*n];
    F77NAME(dgemm)('T', 'T', n, n, n, 1.0, a, n, b, n, 0.0, d, n);
    std::cout << "Resulting Matrix True:\n";
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            std::cout << d[i * n + j] << " ";
        }
        std::cout << "\n";
    }
    delete[] a;
    delete[] b;
    delete[] c;
    delete[] d;
}