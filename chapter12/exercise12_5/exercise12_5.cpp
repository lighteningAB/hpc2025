//Write a parallel MPI program which generates a random matrix and vector, distributed evenly
//amongst the processes, and performs a matrix-vector multiplication in parallel, printing the
//result on the root process. Use a matrix/vector size of 2^n and run on 2^m processes, where
//n >> m.