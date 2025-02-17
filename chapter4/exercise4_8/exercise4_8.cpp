//use blas to operate on submatrices to recursively multiply two square matrices of dimensions nxn

#define F77NAME(x) x##_

extern "C"{
    void F77NAME(dgemm)(char &transa, char &transb)
}

//function has input of two nxn arrays where n is a multiple of 2 and returns the product
