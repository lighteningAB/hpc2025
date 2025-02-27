// parallel program, two 1024 floating point numbers, calculate dot product of both vectors and the norm of
// each vector, print result on process 0

#include <iostream>
#include <random>
#include <mpi.h>

#define F77NAME(x) x##_
extern "C"
{
    double F77NAME(ddot)(const int &n, const double *dx, const int &incx, const double *dy, const int &incy);
    double F77NAME(dnrm2)(const int &n, const double *x, const int &incx);
}

int main(int argc, char *argv[])
{
    double norm1;
    double norm2;
    double dotprod;
    int rank, size, retval_rank, retval_size;
    MPI_Init(&argc, &argv);
    retval_rank = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    retval_size = MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (retval_rank == MPI_ERR_COMM || retval_size == MPI_ERR_COMM)
    {
        std::cout << "Invalid communicator" << std::endl;
        return 1;
    }

    if (size > 32 || !(size && (!(size & (size - 1)))))
    {
        std::cout << "process number must be positive power of 2 and less than 32" << std::endl;
        return 1;
    }

    double *vec1 = new double[1024 / size];
    double *vec2 = new double[1024 / size];

    for (int i = 0; i < 1024 / size; i++)
    {
        vec1[i] = (double)rand() / RAND_MAX * 10.0;
        vec2[i] = (double)rand() / RAND_MAX * 10.0;
    }

    dotprod = F77NAME(ddot)(1024 / size, vec1, 1, vec2, 1);
    norm1 = F77NAME(dnrm2)(1024 / size, vec1, 1);
    norm2 = F77NAME(dnrm2)(1024 / size, vec2, 1);

    // send from all slave nodes
    if (rank != 0)
    {
        MPI_Send(&dotprod, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        MPI_Send(&norm1, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
        MPI_Send(&norm2, 1, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD);
    }
    // receive at root, combine, and show

    else
    {
        double dotprodsum = dotprod;
        double norm1sum = norm1;
        double norm2sum = norm2;
        for (int i = 1; i < size; i++)
        {
            MPI_Recv(&dotprod, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            dotprodsum += dotprod;
            MPI_Recv(&norm1, 1, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            norm1sum += norm1;
            MPI_Recv(&norm2, 1, MPI_DOUBLE, i, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            norm2sum += norm2;
        }
        std::cout << "the norm of vector 1 is: " << norm1sum << std::endl;
        std::cout << "the norm of vector 2 is: " << norm2sum << std::endl;
        std::cout << "the dot product is: " << dotprodsum << std::endl;
    }

    delete[] vec1;
    delete[] vec2;
    MPI_Finalize();

    return 0;
}