//solve 4.4 with lapack, i.e. solve matrix equations of form Ax=b

#include "exercise5_1.h"

void convertToFull(double* packed, double* full, int n) {
    // Zero out full matrix first
    for (int i = 0; i < n * n; i++) {
        full[i] = packed[i];
    }

    for (int i = 0; i<n; i++){
        for (int j = i+1; j<n; j++){
            full[j*n+i] = packed[i*n+j];
        }
    }
}

int main(){
    int n = 19;
    double *symmat = new double[n * n];
    double alpha = -2.0/pow((.1),2)-1.0;
    double beta = 1.0/(pow(.1,2));
    symmetricColMaj(alpha, beta, 19, symmat);

    for(int i = 0; i<n; i++){
        for(int j = 0; j<n; j++){
            std::cout<<symmat[i*n+j]<<" ";
        }
        std::cout<<std::endl;
    }

    double *b = new double[n];
    for (int i = 0; i<n;i++){
        b[i] = (i+1)*2.0/20.0;
    }

    forcfunc(b, n);

    std::cout << std::endl;
    double *x0 = new double[n];
    conjgradsolve(n, symmat, b, x0);
    //function solution
    for (int i = 0; i < n; i++)
    {
        std::cout << x0[i] << " ";
    }
    std::cout << std::endl;
    std::cout << std::endl;

    //solution
    for (int i = 0; i<n; i++){
        std::cout << sin(M_PI*((i+1)*0.1)) << " ";
    }
    std::cout<<std::endl;

    //compare
    for (int i = 0; i<n; i++){
        std::cout <<sin(M_PI*((i+1)*0.1)) - (x0[i]) << " ";
    }
    std::cout<<std::endl;   

    std::cout<<"using lapack"<<std::endl;

    double* fullMatrix = new double[n * n];
    convertToFull(symmat, fullMatrix, n);

    int* ipiv = new int[n]; // Vector for pivots
    int info = 0;
    F77NAME(dgesv)(n, 1, fullMatrix, n, ipiv, b, n, info);

    //lapack solution
    for (int i = 0; i < n; i++)
    {
        std::cout << b[i] << " ";
    }
    std::cout << std::endl;
    std::cout << std::endl;

    //compare
    for (int i = 0; i<n; i++){
        std::cout <<sin(M_PI*((i+1)*0.1)) - (b[i]) << " ";
    }
    std::cout<<std::endl;   

    delete[] ipiv;
    delete[] symmat;
    delete[] b;
    delete[] x0;
}
