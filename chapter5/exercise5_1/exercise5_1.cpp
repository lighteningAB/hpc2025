//solve 4.4 with lapack, i.e. solve matrix equations of form Ax=b

#include "exercise5_1.h"

void conjgradsolvelapack(int n, double *a, double* b){
    int info = 0;
    F77NAME(dposv)('L', n, 1, a, n, b, 1, info);
}

int main(){
    int n = 19;
    double *symmat = new double[n * n];
    double alpha = -2.0/pow((.1),2)-1.0;
    double beta = 1.0/(pow(.1,2));
    symmetricColMaj(alpha, beta, 19, symmat);

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
    conjgradsolvelapack(n, symmat, b);

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


    delete[] symmat;
    delete[] b;
    delete[] x0;
}
