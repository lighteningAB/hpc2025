//parallel program, two 1024 floating point numbers, calculate dot product of both vectors and the norm of
//each vector, print result on process 0

#include <iostream>
#include <mpi.h>

#define F77NAME(x) x##_
extern "C"
{
    double F77NAME(ddot)(const int&n, const double* dx, const int&incx, const double* dy, const int& incy);
}