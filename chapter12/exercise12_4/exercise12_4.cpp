//trapezium pi expansion
//take in an integer amount of intervals to evaluate integral of 1/(1+x^2) from 0 to 1
//split these among the procceses
//profit

#include <iostream>
#include <mpi.h>

int main(int argc, char*argv[]){
    int t_intervals; //total intervals
    int p_intervals; //amount of intervals assigned to the process
    double localsum = 0.0;
    double dx;
    double startpoint;
    double total = 0.0;
    int rank, size, retval_rank, retval_size;
    MPI_Init(&argc, &argv);
    retval_rank = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    retval_size = MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (retval_rank == MPI_ERR_COMM || retval_size == MPI_ERR_COMM)
    {
        std::cout << "Invalid communicator" << std::endl;
        return 1;
    }
    
    //root process input
    if(rank == 0)
    {
        std::cout<<"Enter the number of intervals: ";
        std::cin>>t_intervals;
        std::cout<<std::endl;
    }
    MPI_Bcast(&t_intervals, 1, MPI_INT, 0, MPI_COMM_WORLD);

    //allocating intervals
    if (rank == 0){
        p_intervals = std::max(0, t_intervals-(size-1)*(t_intervals/size+1));
    }
    else{
        p_intervals = t_intervals/size+1;
    }
    //each process summation
    
    dx = 1.0/double(t_intervals);
    startpoint = std::max(0, t_intervals-(size-1)*(t_intervals/size+1))*std::min(1,rank)+(std::max(0,rank-1))*(p_intervals);
    startpoint = dx*double(startpoint);
    for(int i = 0; i<p_intervals; i++){
        localsum+=(1.0/(1.0+startpoint*startpoint) + 1.0/(1.0+(startpoint+dx)*(startpoint+dx)))/2.0*dx;   
        startpoint+=dx;
    }

    //collect vals
    rank == 0? (total+=localsum) : total = total;
    if (rank == 0){
        for(int i = 1; i<size; i++){
            MPI_Recv(&localsum, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            total+=localsum;
        }
    }
    else{
        MPI_Send(&localsum, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }

    if(rank == 0)
    { 
        std::cout<<"pi approxmated as "<<total*4<<std::endl;
    }
    //

    MPI_Finalize();
}
