#ifndef EXERCISE5_1_H //include guard
#define EXERCISE5_1_H

#include "../../chapter4/exercise4_4/exercise4_4.h"

//include lapack function to solve linear system
#define F77NAME(x) x##_
extern "C"{
    void F77NAME(dposv)(const char &uplo, const int &n, const int &nrhs, double *a, const int &LDA, double *b,
                        const int &LDB, int info);
}

void conjgradsolvelapack(double *a, double* b, int n);
#endif