#include <mpi.h>
#include <bits/stdc++.h>

using namespace std;

int main(int argc, char** argv) {
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
      cout << N << endl;
      data = (int *) malloc(N*sizeof(int));
      for(int i = 0; i < N; i++) {
	cin >> data[i];
      }
    }
    
    MPI_Bcast(&N,1,MPI_INT,0,MPI_COMM_WORLD);
    int split_size = ceil(static_cast<double>(N)/world_size);
    pdata = (int *) malloc(split_size*sizeof(int));
    MPI_Scatter(data,split_size,MPI_INT,pdata,split_size,MPI_INT,0,MPI_COMM_WORLD);

    // Get the name of the processor
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    // Print off a hello world message
    cout << "Elements of processor " << world_rank << ": " << endl;
    for(int i = 0; i < split_size; i++) {
      cout << pdata[i] << " ";
    }
    cout << endl;

    // Finalize the MPI environment.
    MPI_Finalize();
}