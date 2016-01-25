#include <mpi.h>
#include <bits/stdc++.h>
#include "io/QuickSortExecuter.cpp"

using namespace std;

void executeQuickSort(int N, string datatype, MPI_Comm comm, int seed, string &benchmark, ifstream *stream = nullptr) {
  int rank;
  MPI_Comm_rank(comm,&rank);
  SORT_TYPE type = SortingDatatype<void>::getSortType(datatype,rank);
  switch(type) {
    case SORT_TYPE::Int: {
      QuickSortExecuter<int> *qs_exec = new QuickSortExecuter<int>(N,comm,stream);
      qs_exec->execute(seed,benchmark); }
      break;
    case SORT_TYPE::Long: {
      QuickSortExecuter<long long> *qs_exec = new QuickSortExecuter<long long>(N,comm,stream);
      qs_exec->execute(seed,benchmark); }
      break;
    case SORT_TYPE::Double: {
      QuickSortExecuter<double> *qs_exec = new QuickSortExecuter<double>(N,comm,stream);
      qs_exec->execute(seed,benchmark); }
      break;
    case SORT_TYPE::Float: {
      QuickSortExecuter<float> *qs_exec = new QuickSortExecuter<float>(N,comm,stream);
      qs_exec->execute(seed,benchmark); }
      break;
  }
}

int main(int argc, char** argv) {
    ios::sync_with_stdio(false);
    int seed = atoi(argv[2]);
    srand (seed);
    
    // Initialize the MPI environment
    MPI_Init(NULL, NULL);

    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    
    int N, split_size;
    string benchmark = argv[1];
    if(world_rank == 0) {
      string datatype;
      ifstream stream(benchmark);
      if(stream.is_open()) {
	stream >> N >> datatype;
      }
      int datatype_string_size = datatype.size();
      MPI_Bcast(&N,1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(&datatype_string_size,1,MPI_INT,0,MPI_COMM_WORLD);
      char *c_datatype = (char *) datatype.c_str();
      MPI_Bcast(c_datatype,datatype_string_size,MPI_CHAR,0,MPI_COMM_WORLD);
      executeQuickSort(N,datatype,MPI_COMM_WORLD,seed,benchmark,&stream);
    } else {
      int datatype_string_size;
      MPI_Bcast(&N,1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(&datatype_string_size,1,MPI_INT,0,MPI_COMM_WORLD);
      char *c_datatype = (char *) malloc(datatype_string_size*sizeof(char));
      MPI_Bcast(c_datatype,datatype_string_size,MPI_CHAR,0,MPI_COMM_WORLD);
      string datatype(c_datatype,datatype_string_size);
      executeQuickSort(N,datatype,MPI_COMM_WORLD,seed,benchmark);
    }
    
    // Finalize the MPI environment.
    MPI_Finalize();
}