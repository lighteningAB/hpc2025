#include <iostream>
#include <cmath>
#include <fstream>

#define F77NAME(x) x##_

extern "C"
{
    void F77NAME(dsteqr)(const char &compz, const int &n, double *d,
                         double *e, double *z, const int &ldz,
                         double *work, int &info);
    void F77NAME(dgtsv)(const int &n, const int &nrhs, double *dl,
                        double *d, double *du, double* b, int &ldb, int &info);
}

/**
 * @brief generates explicit operator in tridiagonal format
 * for forward euler time integration
 *
 * @param v double, dt/(dx^2), needs to be <0.5 for stability
 * @param returnable double*, 3n-2 vector where explicit operator is stored
 * @param n int, width and height of explicit operator square matrix
 */
void genA(double v, double *returnable, int n)
{
    int position = 0;
    int info = 0;
    double *fill = new double[1];
    for(int i = 0; i<n-1; i++){
        returnable[position] = v;
        position+=1;
    } 
    for(int i = 0; i<n; i++){
        returnable[position] = 1-2*v;
        position+=1;
    } 
    for(int i = 0; i<n-1; i++){
        returnable[position] = v;
        position+=1;
    } 
    double *copy = new double[2 * n - 1];
    std::copy(returnable, returnable + 2 * n - 2, copy);
    // check if stable
    F77NAME(dsteqr)('N', n, copy+n-1, copy, fill, 1, fill, info);
    int unstable = 0;
    for (int i = 0; i < n; i++)
    {
        //std::cout<<returnable[i]<<" ";
        if (abs(copy[i]) >= 1)
        {
            unstable = 1;
        }
    }
    //unstable == 0 ? std::cout << "the matrix A is stable" : std::cout << "the matrix A is unstable";
    //std::cout << std::endl;
    delete[] fill;
    delete[] copy;
}
/**
 * @brief repeatedly solve the equation Ax=U for tridiagonal matrix A
 * 
 * @param A double*, tridiagonal matrix in order N-1 upper, N diagonal, N-1 lower
 * @param U double*, vector U, on succesful completion will contain solution X
 * @param n int, dimension n of square matrix A
 */
void repeatSolve(double * A, double * U, int n){
    int info = 0;
    F77NAME(dgtsv)(n, 1, A, A+n-1, A+2*n-1, U, n, info);
}

//generate vector
void initialVec(int gridp, double * x){
    for(int i = 0; i<gridp; i++){
        x[i] = sin(M_PI*double(i)/double(gridp-1));
    }
}

int main()
{
    //1
    int n = 5;
    double *returnable = new double[3 * n - 2];
    genA(0.4, returnable, n);
    genA(20, returnable, n);
    //2
    int gridp;
    double ss;
    double ts;
    std::cout<<"how many grid points to use?: ";
    std::cin>>gridp;  
    std::cout<<std::endl;
    std::cout<<"what stepsize?: ";
    std::cin>>ss;
    std::cout<<std::endl;
    std::cout<<"amount of timesteps?: ";
    std::cin>>ts;
    std::cout<<std::endl;

    double*A = new double[gridp*3-2];
    double*U = new double[gridp];
    double dx = 1.0/double(gridp);
    double v = -1*(ss/dx)/dx;
    std::cout<<v<<" "<<std::endl;
    initialVec(gridp, U);
    U[0] = 0;
    U[gridp-1] = 0;
    std::ofstream file("vector_evolution.csv");
    for(int i = 0; i<ts; i++){
        for (int j = 0; j < gridp; j++) {
            file << U[j] << (j < gridp-1 ? "," : "\n");
        }
        genA(v, A, gridp);
        A[gridp-2] = 0;
        A[2*gridp-1] = 0;
        repeatSolve(A, U, gridp);
        U[0] = 0;
        U[gridp-1] = 0;
    }
    file.close();
    std::cout << "Data saved to vector_evolution.csv\n";
    delete[] returnable;
    delete[] A;
    delete[] U;
}
