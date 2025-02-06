#include "exercise4_4.h"

void matvecmul(double *a, double *b, double *c, int n)
{
    F77NAME(dsymv)('L', n, 1.0, a, n, b, 1, 0.0, c, 1);
}

double vecmult(int n, double *x, double *y)
{
    return F77NAME(ddot)(n, x, 1, y, 1);
}
// generate reduced matrix symmetric n=21
// column major format
// structure of B a B repeating, row 0 in positions -1 (not in) 0 , 1
// row 1: 0,1,2
// row 18: 17,18, 19 (not in)
//  B = 2 a = 0

// n is the dimensions of the square array
// returnable should be of appropriate size n*n
void symmetricRowMaj(double a, double b, int n, double *returnable)
{
    int counter = 0;
    // for each row
    for (int i = 0; i < n; i++)
    {
        // for each row we want to generate the requisite number of 0s then populate in b,a where a is on the diagonal
        // then add in the remaining 0s

        // we will need j zeroes, if j is negative it should be equivalent to 0
        for (int j = i - 1; j > 0; j--)
        {
            returnable[counter] = 0;
            counter += 1;
        }

        // we will need to populate the next elements of the array with beta alpha, the edge case exists on
        // the first row, maybe there is a smart way to do it but for now i will hardcode them
        if (i == 0)
        { // first row
            returnable[counter] = a;
            counter += 1;
        }
        else
        {
            returnable[counter] = b;
            returnable[counter + 1] = a;
            counter += 2;
        }

        for (int j = 0; j < n - 2 - i; j++)
        {
            returnable[counter] = 0;
            counter += 1;
        }
    }
}

// same as above but column major
void symmetricColMaj(double a, double b, int n, double *returnable)
{
    int counter = 0;
    for (int i = 0; i < n; i++)
    { // each column
        for (int j = 0; j < i; j++)
        { // initial 0s
            returnable[counter] = 0;
            counter += 1;
        }
        returnable[counter] = a;
        if (i != n - 1)
        {
            returnable[counter + 1] = b;
        }
        counter += 2;
        for (int j = 0; j < n - 2 - i; j++)
        {
            returnable[counter] = 0;
            counter += 1;
        }
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

void forcfunc2(double *a, int n)
{
    // function is -(theta+pi^2)cos(pi*x)
    for (int i = 0; i < n; i++)
    {
        a[i] = -1 * (1 + pow(M_PI, 2)) * cos(M_PI * a[i]);
    }
}

void conjgradsolve(int n, double *a, double *b, double *x0)
{
    std::fill(x0, x0 + n, 0.0);
    double *r = new double[n];
    std::copy(b, b + n, r);
    double *c = new double[n];
    matvecmul(a, x0, c, n);            // Ax0
    F77NAME(daxpy)(n, -1, c, 1, r, 1); // r=r0 ___ r0 = b-Ax0
    double *p = new double[n];
    std::copy(r, r + n, p); // p0 = r0, cloning r0
    int k = 0;
    int max_iter = 10000; // can make these inputs to the function
    double stop = 40;
    double *bottomhelpervec = new double[n]; // allocate memory
    double *helpervec2 = new double[n];
    double *r_prev = new double[n];
    double *p_prev = new double[n];
    while (k < max_iter)
    {
        matvecmul(a, p, bottomhelpervec, n);
        double a_k = vecmult(n, r, r) / vecmult(n, p, bottomhelpervec); // a_k = r_k^t * r_k / (p_k^t*A*p_k)
        F77NAME(daxpy)(n, a_k, p, 1, x0, 1);                            // x_k+1 = x_k + a_k*p_k
        std::copy(r, r + n, r_prev);
        matvecmul(a, p, helpervec2, n);
        F77NAME(daxpy)(n, a_k * -1, helpervec2, 1, r, 1); // r_k+1 = r_k - a_k*a*p_k
        stop = F77NAME(dnrm2)(n, r, 1);
        if (stop < 0.00000000000000001)
        {
            break;
        }
        double B_k = vecmult(n, r, r) / vecmult(n, r_prev, r_prev); // B_k = r_k+1^t * r_k+1 / r_k^t * r_k
        std::copy(p, p + n, p_prev);
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
    double *symmat = new double[n * n];
    double alpha = -2.0/pow((.1),2)-1.0;
    double beta = 1.0/(pow(.1,2));
    symmetricColMaj(alpha, beta, 19, symmat);
    /*
    for (int i = 0; i<n; i++){
        for(int j = 0; j<n; j++){
            std::cout<<symmat[i*n+j]<<" ";
        }
        std::cout<<std::endl;
    }
    */
    double *b = new double[n];
    for (int i = 0; i<n;i++){
        b[i] = (i+1)*2.0/20.0;
    }
    /*
    for (int i = 0; i < n; i++)
    {
        std::cout << b[i] << " ";
    }
    std::cout << std::endl;
    */
    forcfunc(b, n);
    /*
    for (int i = 0; i < n; i++)
    {
        std::cout << b[i] << " ";
    }
    */
    std::cout << std::endl;
    double *x0 = new double[n];
    conjgradsolve(n, symmat, b, x0);
    for (int i = 0; i < n; i++)
    {
        std::cout << x0[i] << " ";
    }
    std::cout << std::endl;
    std::cout << std::endl;
    for (int i = 0; i<n; i++){
        std::cout << sin(M_PI*((i+1)*0.1)) << " ";
    }
    std::cout<<std::endl;

    for (int i = 0; i<n; i++){
        std::cout <<sin(M_PI*((i+1)*0.1)) - (x0[i]) << " ";
    }
    std::cout<<std::endl;

    int m = 29;
    double *symmat2 = new double[m * m];
    alpha = -2.0/pow((.1),2)-1.0;
    beta = 1.0/(pow(.1,2));
    symmetricColMaj(alpha, beta, 29, symmat2);
    /*
    for (int i = 0; i<n; i++){
        for(int j = 0; j<n; j++){
            std::cout<<symmat[i*n+j]<<" ";
        }
        std::cout<<std::endl;
    }
    */
    double *c = new double[m];
    for (int i = 0; i<m;i++){
        c[i] = (i+1)*.1;
    }
    /*
    for (int i = 0; i < m; i++)
    {
        std::cout << c[i] << " ";
    }
    */
    std::cout << std::endl;
    forcfunc2(c, m);
    c[0] += -1/pow(.1,2);
    c[m-1] += 1/pow(.1,2);
    /*
    for (int i = 0; i < n; i++)
    {
        std::cout << c[i] << " ";
    }
    */
    std::cout << std::endl;
    double *x02 = new double[m];
    conjgradsolve(m, symmat2, c, x02);
    for (int i = 0; i < m; i++)
    {
        std::cout << x02[i] << " ";
    }
    std::cout << std::endl;
    

    for (int i = 0; i<m; i++){
        std::cout << cos(M_PI*((i+1)*0.1)) << " ";
    }
    std::cout<<std::endl;

    for (int i = 0; i<m; i++){
        std::cout << cos(M_PI*((i+1)*0.1)) - (x02[i]) << " ";
    }
    std::cout<<std::endl;

    delete[] symmat;
    delete[] symmat2;
    delete[] b;
    delete[] c;
    delete[] x0;
    delete[] x02;
}