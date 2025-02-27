#include <iostream>
#include <mpi.h>

int main(int argc, char* argv[]){
    int input;
    int rank, size, retval_rank, retval_size;
    MPI_Init(&argc, &argv);
    retval_rank = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    retval_size = MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (retval_rank == MPI_ERR_COMM || retval_size == MPI_ERR_COMM) {
        std::cout << "Invalid communicator" << std::endl;
        return 1;
    }

    if (rank == 0){
        std::cout<<"enter an integer: ";
        std::cin>>input;
        std::cout<<std::endl;
        for (int i = 1; i<size; i++){
            MPI_Send(&input, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
        }
    }
    else{
        MPI_Recv(&input, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    std::cout<<"process of rank "<< rank <<", user integer is: "<<input<<std::endl;
    MPI_Finalize();
    return 0;
}