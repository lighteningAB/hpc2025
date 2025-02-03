#include <math.h>
// helmholtz

// generate reduced matrix symmetric n=21
// structure of B a B repeating, row 0 in positions -1 (not in) 0 , 1
// row 1: 0,1,2
// row 18: 17,18, 19 (not in)
//  B = 2 a = 0

// n is the dimensions of the square array
// returnable should be of appropriate size
//  [(1+n)*n]/2
void symmetric(double a, double b, int n, double *returnable)
{
    int counter = 0;
    // for each row
    for (int i = 0; i < n; i++)
    {
        // for each row we want to generate the requisite number of 0s then populate in b,a where a is on the diagonal

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

void conjgradsolve(int n, double* a, double *b){
    
}
