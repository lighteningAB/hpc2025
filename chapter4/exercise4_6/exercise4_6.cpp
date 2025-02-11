#include <iostream>
#include <cmath>
#include <fstream>

#define F77NAME(x) x##_

extern "C"{
    void F77NAME(dgbmv)(const char &trans, const int &m, const int &n, const int &kl, const int &ku, const double &alpha,
        const double *a, const int &LDA, const double *x, const int &incx, const double &beta, double *y, const int &incy);
}

void matvecmulband(double *a, double *b, double *c, int n, int m, int kl, int ku)
{
    F77NAME(dgbmv)('N', n, n, kl, ku, 1.0, a, 1+kl+ku, b, 1, 0.0, c, 1);
}

//generate banded matrix
//banded matrix 3*n
void genA(double ts, double dx, double * returnable, int n){
    double v = ts/dx/dx;
    int position = 0;
    for (int i = 0; i < n; i++)
    {
        returnable[position] = v;
        returnable[position + 1] = 1.0-2.0*v;
        returnable[position + 2] = v;
        position += 3;
    }
}

//generate vector
void initialVec(double ss, int gridp, double * x){
    for(int i = 0; i<gridp; i++){
        x[i] = sin(M_PI*i*ss);
    }
}

int main(){
    //getting inputs
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

    double * A = new double[3*gridp];
    double * x = new double[gridp];
    double * x1 = new double[gridp];
    genA(ts, ss, A, gridp);
    initialVec(ss, gridp, x);
    int tries = 100;
    std::ofstream file("vector_evolution.csv");
    for(int i = 0; i<tries; i++){
        for (int i = 0; i < gridp; i++) {
            file << x[i] << (i < gridp - 1 ? "," : "\n");
        }
        matvecmulband(A, x, x1, gridp, 3, 1, 1);
        std::copy(x, x+gridp, x1);
    }
    file.close();
    std::cout << "Data saved to vector_evolution.csv\n";
}

