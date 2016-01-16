#include <mpi.h>
#include <bits/stdc++.h>

using namespace std;

/* random-number generator used for values
 * which are identical on every PE
 * The constants can be found Press Et Al.
 * Numerical Recipes in C, 2nd edition, page 284 
 */
unsigned long globalRandState;
#define sGlobalRand(s) (globalRandState = (s))
#define globalRand() (\
   globalRandState = ((1664525L*globalRandState + 1013904223L) & 0xffffffffL))
#define globalRandInt(n) ((globalRand() >> 10) % (n))

template<typename T>
class QuickSort {
 
public:
  QuickSort(int rank, int world_size, int N, T *data) : rank(rank), world_size(world_size), N(N), data(data) { }
  
  void quickSort(int i, int j) {
    int pivot = getPivot();
    int largeIdx = 0;
    for(int i = 0; i < N; i++) {
      if(data[i] <= pivot) {
	swap(data[i],data[largeIdx]);
	largeIdx++;
      }
    }
    int prefix_sum;
    MPI_Scan(&largeIdx,&prefix_sum,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    prefix_sum -= largeIdx;
    int *new_data = (int *) malloc(N*sizeof(int));
    for(int i = 0; i < N; i++)
      new_data[i] = 0;
    
    for(int i = 0; i < N; i++) {
      cout << data[i] << "("<<rank<<") ";
    }
    cout << endl;
    MPI_Barrier(MPI_COMM_WORLD);
    
    sendData(0,prefix_sum,data,new_data,largeIdx);
    prefix_sum += largeIdx;
    int small_elements = prefix_sum;
    if(rank == world_size-1)
      cout << "Small Elements: " << small_elements << endl;
    MPI_Bcast(&small_elements,1,MPI_INT,world_size-1,MPI_COMM_WORLD);
    int large_elements = N-largeIdx;
    MPI_Scan(&large_elements,&prefix_sum,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    prefix_sum -= large_elements;
    sendData(small_elements,prefix_sum,data+largeIdx,new_data,large_elements);
    
    swap(data,new_data);
    for(int i = 0; i < N; i++) {
      cout << data[i] << "("<<rank<<") ";
    }
    cout << endl;
    
  }
  
private:
  
  int getPivot() {
    int pivot = data[rand() % N];
    int pivotPE = globalRandInt(world_size);
    MPI_Bcast(&pivot,1,MPI_INT,pivotPE,MPI_COMM_WORLD);
    if(rank == pivotPE)
      cout << "Pivot: "<< pivot << endl;
    return pivot;
  }
  
  void sendData(int start_idx, int prefix_sum, int *data, int *new_data, int length) {
    pair<int,int> process_idx = getProcessIdx(start_idx+prefix_sum);
    cout << process_idx.first << " " << process_idx.second << " " << length << endl;
    int *data_send = (int *) malloc(world_size*sizeof(int));
    int *data_length_send = (int *) malloc(world_size*sizeof(int));
    int *data_recv = (int *) malloc(world_size*sizeof(int));
    int *data_length_recv = (int *) malloc(world_size*sizeof(int));
    for(int i = 0; i < world_size; i++) {
      data_length_send[i] = 0;
      if(i == process_idx.first) {
	data_send[i] = process_idx.second;
	data_length_send[i] = min(N-data_send[i],length);
	prefix_sum += data_length_send[i];
	length -= data_length_send[i];
	if(length == 0)
	  process_idx.first = -1;
	else {
	  process_idx = getProcessIdx(start_idx+prefix_sum);
	}
      } else 
	data_send[i] = -1;
    } 
    MPI_Alltoall(data_send,1,MPI_INT,data_recv,1,MPI_INT,MPI_COMM_WORLD);
    MPI_Alltoall(data_length_send,1,MPI_INT,data_length_recv,1,MPI_INT,MPI_COMM_WORLD);
    
    int idx = 0;
    for(int i = 0; i < world_size; i++) {
	if(data_send[i] != -1) {
	  MPI_Send(data+idx,data_length_send[i],MPI_INT,i,42,MPI_COMM_WORLD);
	  idx += data_length_send[i];
	}
    }
    for(int i = 0; i < world_size; i++) {
	MPI_Status status;
	if(data_recv[i] != -1)
	  MPI_Recv(new_data+data_recv[i],data_length_recv[i],MPI_INT,i,42,MPI_COMM_WORLD,&status);
    }
    
  }
  
  void printArray(int N, int *array) {
    for(int i = 0; i < N; i++)
      cout << array[i] << "("<<rank<<") ";
    cout << endl;
  }
  
  pair<int,int> getProcessIdx(int idx) {
    return make_pair(idx/N,idx-(idx/N)*N);
  }
  
  int rank, world_size, N;
  T *data;
  
};