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
    }
    else{
        MPI_Recv(&input, 1, MPI_INT, rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        std::cout<<"process of rank "<< rank <<", user integer is: "<<input<<" received from "<< rank-1 <<std::endl;
    }

    if (rank!=size-1){
        MPI_Send(&input, 1, MPI_INT, rank+1, 0, MPI_COMM_WORLD);
    }
    else{
        MPI_Send(&input, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }
    
    if(rank == 0){
        std::cout<<"process of rank "<< rank <<", user integer is: "<<input<<" received from "<< size <<std::endl;
    }


    MPI_Finalize();
    return 0;
}