// use banded matrices for 4_4
#include "exercise4_5.h"

void matvecmulband(double *a, double *b, double *c, int n, int m, int kl, int ku)
{
    F77NAME(dgbmv)('N', n, n, kl, ku, 1.0, a, 3, b, 1, 0.0, c, 1);
}

double vecmult(int n, double *x, double *y)
{
    return F77NAME(ddot)(n, x, 1, y, 1);
}

// returnable should be of size 3 x n for our purposes
void bandcolMaj(double a, double b, int n, double *returnable)
{
    int position = 0;
    for(int i = 0; i<n; i++){
        returnable[position] = b;
        returnable[position+1] = a;
        returnable[position+2] = b;
        position+=3;
    }
}

// forcing func calculation
void forcfunc(double *a, int n)
{
    // function is -(theta+pi^2)sin(pi*x)
    for (int i = 0; i < n; i++)
    {
        a[i] = -1 * (1 + pow(M_PI, 2)) * sin(M_PI * a[i]);
    }
}

void conjgradsolve(int n, double *a, double *b, double *x0)
{
    std::fill(x0, x0 + n, 0.0);
    double *r = new double[n];
    std::copy(b, b+n, r);
    double *c = new double[n];
    matvecmulband(a, x0, c, n, 3, 1, 1); //Ax0
    F77NAME(daxpy)(n, -1, c, 1, r, 1); // r=r0 ___ r0 = b-Ax0
    double *p = new double[n];
    std::copy(r, r+n, p); //p0 = r0, clone r0
    int k = 0;
    int max_iter = 100;
    double stop = 40;
    double *bottomhelpervec = new double[n];
    double *helpervec2 = new double[n];
    double *r_prev = new double[n];
    double *p_prev = new double[n];
    while (k < max_iter)
    {
        matvecmulband(a, p , bottomhelpervec, n, 3, 1, 1); 
        double a_k = vecmult(n, r, r)/ vecmult(n, p, bottomhelpervec);
        F77NAME(daxpy)(n, a_k, p, 1, x0, 1);
        std::copy(r, r+n, r_prev);
        matvecmulband(a, p, helpervec2, n, 3 , 1, 1);
        F77NAME(daxpy)(n, a_k * -1, helpervec2, 1, r, 1); // r_k+1 = r_k - a_k*a*p_k
        stop = F77NAME(dnrm2)(n, r, 1);
        //std::cout<<stop<<std::endl;
        if (stop < 0.0000000000000000001){
            break;
        }
        double B_k = vecmult(n, r, r)/ vecmult(n, r_prev, r_prev);
        std::copy(p, p+n, p_prev);
        std::copy(r, r + n, p);                  // p = r_k+1
        F77NAME(daxpy)(n, B_k, p_prev, 1, p, 1); // p_k+1 = r_k+1 + B_k*p_k
        k += 1;
    }

    delete[] r;
    delete[] c;
    delete[] p;
    delete[] bottomhelpervec;
    delete[] helpervec2;
    delete[] r_prev;
    delete[] p_prev;
    std::cout << "found in " << k << " iterations" << std::endl;
}

int main()
{
    int n = 19;
    double *bandmat = new double[n * 3];
    double alpha = -2.0/pow((.1),2)-1.0;
    double beta = 1.0/(pow(.1,2));
    bandcolMaj(alpha, beta, n, bandmat);
    for (int i = 0; i<3; i++){
        for(int j = 0; j<n; j++){
            std::cout<<bandmat[i*n+j]<<" ";
        }
        std::cout<<std::endl;
    }
    double *b = new double[n];
    for(int i = 0; i<n; i++){
        b[i] = (i+1.0)*2.0/20.0;
    }
    forcfunc(b, n);
    double *x0 = new double[n];
    conjgradsolve(n, bandmat, b, x0);
    for (int i = 0; i < n; i++)
    {
        std::cout << x0[i] << " ";
    }
    std::cout << std::endl;
    for (int i = 0; i<n; i++){
        std::cout << sin(M_PI*((i+1)*0.1)) << " ";
    }
    std::cout<<std::endl;

    for (int i = 0; i<n; i++){
        std::cout <<sin(M_PI*((i+1)*0.1)) - (x0[i]) << " ";
    }
    std::cout<<std::endl;
}