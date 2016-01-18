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
    
    int N, split_size, *data, *pdata; 
    if(world_rank == 0) {
      cin >> N;
      split_size = ceil(static_cast<double>(N)/world_size);
      data = (int *) malloc(split_size*world_size*sizeof(int));
      for(int i = 0; i < N; i++) {
	cin >> data[i];
      }
      for(int i = N; i < split_size*world_size; i++)
	data[i] = INT_MAX;
    }
    
    MPI_Bcast(&N,1,MPI_INT,0,MPI_COMM_WORLD);
    split_size = ceil(static_cast<double>(N)/world_size);
    pdata = (int *) malloc(split_size*sizeof(int));
    MPI_Scatter(data,split_size,MPI_INT,pdata,split_size,MPI_INT,0,MPI_COMM_WORLD);

    QuickSort<int> *qs = new QuickSort<int>(split_size,pdata);
    qs->quickSort(0,0,N, MPI_COMM_WORLD);
    int* d = qs->getData();
    if(world_rank == 0) {
      for(int i = 0; i < N; i++)
	cout << data[i] << " ";
      cout << endl;
    }
    MPI_Gather(d,split_size,MPI_INT,data,split_size,MPI_INT,0,MPI_COMM_WORLD);
    if(world_rank == 0) {
      for(int i = 0; i < N; i++)
	cout << data[i] << " ";
      cout << endl;
    }
    
    // Finalize the MPI environment.
    MPI_Finalize();
}