#include <mpi.h>
#include <bits/stdc++.h>
#include "QuickSort.h"

using namespace std;

int main(int argc, char** argv) {
    srand (time (NULL));
  
    // Initialize the MPI environment
    MPI_Init(NULL, NULL);

    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    
    int N, *data, *pdata; 
    if(world_rank == 0) {
      cin >> N;
      data = (int *) malloc(N*sizeof(int));
      for(int i = 0; i < N; i++) {
	cin >> data[i];
      }
    }
    
    MPI_Bcast(&N,1,MPI_INT,0,MPI_COMM_WORLD);
    int split_size = ceil(static_cast<double>(N)/world_size);
    pdata = (int *) malloc(split_size*sizeof(int));
    MPI_Scatter(data,split_size,MPI_INT,pdata,split_size,MPI_INT,0,MPI_COMM_WORLD);

    QuickSort<int> *qs = new QuickSort<int>(world_rank,world_size,split_size,pdata);
    qs->quickSort(0,world_size-1);
    
    // Finalize the MPI environment.
    MPI_Finalize();
}