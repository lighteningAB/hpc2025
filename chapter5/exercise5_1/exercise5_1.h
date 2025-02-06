#ifndef EXERCISE5_1_H //include guard
#define EXERCISE5_1_H

#include "exercise4_4.h"

//include lapack function to solve linear system
extern "C"{
    void F77NAME(dgesv)(const int& n, const int& nrhs, const double * A,
                        const int& lda, int * ipiv, double * B,
                        const int& ldb, int& info);
}

void convertToFull(double* packed, double* full, int n);

#endif