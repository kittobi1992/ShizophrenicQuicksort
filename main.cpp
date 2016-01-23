#include <mpi.h>
#include <bits/stdc++.h>
#include <chrono>
#include "QuickSort.h"

using namespace std;
using HighResClockTimepoint = std::chrono::time_point<std::chrono::high_resolution_clock>;

typedef int DATATYPE;

MPI_Datatype mpi_type = GenericDatatype<DATATYPE>::getMPIDatatype();

bool isSorted(DATATYPE* data, int N) {
  for(int i = 1; i < N; i++)
    if(data[i-1] > data[i])
      return false;
  return true;
}

int countWrongIndicies(DATATYPE* data, int N) {
  int idx = 0;
  for(int i = 1; i < N; i++)
    if(data[i-1] > data[i])
      idx++;
  return idx;
} 

int main(int argc, char** argv) {
    ios::sync_with_stdio(false);
    srand (time (NULL));
  
    // Initialize the MPI environment
    MPI_Init(NULL, NULL);

    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    
    int N, split_size;
    DATATYPE *data, *pdata; 
    if(world_rank == 0) {
      cin >> N;
      split_size = ceil(static_cast<double>(N)/world_size);
      data = (DATATYPE *) malloc(split_size*world_size*sizeof(DATATYPE));
      for(int i = 0; i < N; i++) {
	cin >> data[i];
      }
      for(int i = N; i < split_size*world_size; i++)
	data[i] = numeric_limits<DATATYPE>::max();
    }
    
    MPI_Bcast(&N,1,MPI_INT,0,MPI_COMM_WORLD);
    split_size = ceil(static_cast<double>(N)/world_size);
    pdata = (DATATYPE *) malloc(split_size*sizeof(DATATYPE));
    MPI_Scatter(data,split_size,mpi_type,pdata,split_size,mpi_type,0,MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    
    QuickSort<DATATYPE> *qs = new QuickSort<DATATYPE>(split_size,pdata);
    QSInterval<DATATYPE> ival(0,0,N,MPI_COMM_WORLD);
    HighResClockTimepoint start = std::chrono::high_resolution_clock::now();
    qs->quickSort(ival);
    HighResClockTimepoint end = std::chrono::high_resolution_clock::now();
    DATATYPE* d = qs->getData();

    std::chrono::duration<double> start_time = start.time_since_epoch();
    std::chrono::duration<double> end_time = end.time_since_epoch();
    double d_start = start_time.count();
    double d_end = end_time.count();
    double min_start, max_end;
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Reduce(&d_start,&min_start,1,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);
    MPI_Reduce(&d_end,&max_end,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
    
    double *runningTimes = (double *)malloc(world_size*sizeof(double));
    double runningTime = (d_end-d_start);
    MPI_Gather(&runningTime,1,MPI_DOUBLE,runningTimes,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    
    MPI_Gather(d,split_size,mpi_type,data,split_size,mpi_type,0,MPI_COMM_WORLD);
    if(world_rank == 0) {
      bool isSort = isSorted;
      cout << "RESULT "
	   << "algo=ShizophrenicQuicksort "
	   << "time=" << (max_end-min_start) << " "
	   << "isSorted=" << isSort << " ";
      for(int i = 0; i < world_size; i++)
	cout << "timePE" << i << "=" << runningTimes[i] << " ";
      cout << endl;
    }
    
    
    // Finalize the MPI environment.
    MPI_Finalize();
}