#include <mpi.h>
#include <bits/stdc++.h>
#include <future>

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

//Macros
#define V(X) cout << #X << "=" << X << endl
#define W(X) #X << "=" << X

void c(int i) {
  cout << i << endl;
}

template<typename T>
class QuickSort {
 
public:
  QuickSort(int N, T *data) : N(N), data(data) { }
  
  void quickSort(int start_pid, int start, int end, MPI_Comm comm, bool less_equal = true) {
    int world_size, rank;
    MPI_Comm_size(comm,&world_size);
    MPI_Comm_rank(comm, &rank);
    pair<int,int> startIdx = getProcessIdx(start_pid,start), endIdx = getProcessIdx(start_pid,end-1);
    if(world_size == 0)
      return;
    else if(world_size == 1) {
      sort(data+startIdx.second,data+endIdx.second+1);
      return;
    }
    
    
    int largeIdx = 0;
    if(rank != startIdx.first)
      startIdx.second = 0;
    if(rank != endIdx.first)
      endIdx.second = N-1;
        
    int pivot = getPivot(rank,world_size,startIdx.second,endIdx.second+1,comm);
    for(int i = startIdx.second; i < endIdx.second+1; i++) {
      if(less_equal) {
	if(data[i] <= pivot) {
	  swap(data[i],data[startIdx.second+largeIdx++]);
	}
      } else {
	if(data[i] < pivot) {
	  swap(data[i],data[startIdx.second+largeIdx++]);
	}	
      }
    }

    MPI_Barrier(comm);
    
    int prefix_sum;
    MPI_Scan(&largeIdx,&prefix_sum,1,MPI_INT,MPI_SUM,comm);
    prefix_sum -= largeIdx;
    int *new_data = (int *) malloc(N*sizeof(int));
    for(int i = 0; i < N; i++)
      new_data[i] = 0;
    
    //printArray(N,rank,data);
    
    sendData(start,start_pid,prefix_sum,data+startIdx.second,new_data,largeIdx,comm);
    prefix_sum += largeIdx;
    int small_elements = prefix_sum;
    MPI_Bcast(&small_elements,1,MPI_INT,world_size-1,comm);
    
    if(small_elements == end-start) {
      quickSort(start_pid,start,end,comm,false);
      return;
    } else if(small_elements == 0 && !less_equal)
      return;
    
    int large_elements = (endIdx.second-startIdx.second+1)-largeIdx;
    MPI_Scan(&large_elements,&prefix_sum,1,MPI_INT,MPI_SUM,comm);
    prefix_sum -= large_elements;
    sendData(start+small_elements,start_pid,prefix_sum,data+startIdx.second+largeIdx,new_data,large_elements,comm);
    
    
    for(int i = startIdx.second; i < endIdx.second+1; i++)
      data[i] = new_data[i];
    
    /*if(rank == 0) {
      V(pivot);
      V(small_elements);
      V(start);
      V(end);
    }*/
    
    //printArray(N,start_pid+rank,data);
    
    MPI_Comm left, right;
    int shizophrenicPID = -1;
    pair<int,int> new_start_pid = createNewCommunicators(start_pid,start,start+small_elements,end,comm,shizophrenicPID,&left,&right);
    future<void> fut1,fut2;
    if(MPI_COMM_NULL != left)
      fut1 = async(launch::async,&QuickSort::quickSort,this,new_start_pid.first+start_pid,start,start+small_elements,left,true);
    else
      fut1 = async([]() {return;});

    if(shizophrenicPID == rank)
      fut1.get();
    
    if(MPI_COMM_NULL != right)
      fut2 = async(launch::async,&QuickSort::quickSort,this,new_start_pid.second+start_pid,start+small_elements,end,right,true);
    else
      fut2 = async([]() {return;});
    
    if(shizophrenicPID == -1)
      fut1.get();
    fut2.get();
  }
  
   T* getData() {
    return data;
  }
  
private:
  
  int getPivot(int rank, int world_size, int start, int end, MPI_Comm &comm) {
    int c = 5;
    int pivot = 0;
    for(int i = 0; i < c; i++) {
      pivot += data[start + rand() % (end-start)];
    }
    pivot /= c;
    int pivotPE;
    if(rank == 0)
      pivotPE = rand() % world_size;
    MPI_Request *req; MPI_Status *status;
    MPI_Bcast(&pivotPE,1,MPI_INT,0,comm);
    MPI_Bcast(&pivot,1,MPI_INT,pivotPE,comm);
    return pivot;
  }
  
  void sendData(int start_idx, int start_pid, int prefix_sum, int *data, int *new_data, int length, MPI_Comm &comm) {
    
    int world_size, rank;
    MPI_Comm_size(comm,&world_size);
    MPI_Comm_rank(comm, &rank);
    
    pair<int,int> process_idx = getProcessIdx(start_pid,start_idx+prefix_sum);
    int *data_send = (int *) malloc(world_size*sizeof(int));
    int *data_length_send = (int *) malloc(world_size*sizeof(int));
    int *data_recv = (int *) malloc(world_size*sizeof(int));
    int *data_length_recv = (int *) malloc(world_size*sizeof(int));
    int send_count = 0, recv_count = 0;
    for(int i = 0; i < world_size; i++) {
      data_length_send[i] = 0;
      if(i == process_idx.first) {
	if(i != rank)
	  send_count++;
	data_send[i] = process_idx.second;
	data_length_send[i] = min(N-data_send[i],length);
	prefix_sum += data_length_send[i];
	length -= data_length_send[i];
	if(length == 0)
	  process_idx.first = -1;
	else {
	  process_idx = getProcessIdx(start_pid,start_idx+prefix_sum);
	}
      } else 
	data_send[i] = -1;
    } 
    
    MPI_Barrier(comm);
    
    MPI_Alltoall(data_send,1,MPI_INT,data_recv,1,MPI_INT,comm);
    MPI_Alltoall(data_length_send,1,MPI_INT,data_length_recv,1,MPI_INT,comm);
    
    for(int i = 0; i < world_size; i++) {
	if(data_recv[i] != -1 && i != rank)
	  recv_count++;
    }
   
    MPI_Request *send_req = (MPI_Request *) malloc(send_count*sizeof(MPI_Request));
    MPI_Request *recv_req = (MPI_Request *) malloc(recv_count*sizeof(MPI_Request));
    int send_idx = 0, recv_idx = 0;
    for(int i = 0; i < world_size; i++) {
	if(rank == i) {
	  int idx = 0;
	  for(int j = 0; j < world_size; j++) {
	    if(data_send[j] != -1 && j != rank) {
	      MPI_Isend(data+idx,data_length_send[j],MPI_INT,j,42,comm,&send_req[send_idx++]);
	      idx += data_length_send[j];
	    } else if(rank == j && data_send[j] != -1) {
	      for(int k = 0; k < data_length_send[rank]; k++)
		*(new_data+(data_recv[rank]+k)) = data[k+idx];
	      idx += data_length_send[rank];
	    }
	  }
	}
	else {
	  MPI_Status status;
	  if(data_recv[i] != -1)
	    MPI_Irecv(new_data+data_recv[i],data_length_recv[i],MPI_INT,i,42,comm,&recv_req[recv_idx++]);
	}
    }
    
    MPI_Status *send_status = (MPI_Status *) malloc(send_count*sizeof(MPI_Status));
    MPI_Status *recv_status = (MPI_Status *) malloc(recv_count*sizeof(MPI_Status));
    MPI_Waitall(send_count,send_req,send_status);
    MPI_Waitall(recv_count,recv_req,recv_status);
  }
  
  pair<int,int> createNewCommunicators(int start_pid, int start, int middle, int end, MPI_Comm comm, int &shizophrenicPID, MPI_Comm *left, MPI_Comm *right) {
      MPI_Group group;
      MPI_Comm_group(comm,&group);
      pair<int,int> startIdx = getProcessIdx(start_pid,start), middleIdx1 = getProcessIdx(start_pid,middle-1), 
		    middleIdx2 = getProcessIdx(start_pid,middle),endIdx = getProcessIdx(start_pid,end-1);
      int size_left = (middleIdx1.first-startIdx.first+1), size_right = (endIdx.first-middleIdx2.first+1);
      int *ranks_left = (int *) malloc(size_left*sizeof(int));
      int *ranks_right = (int *) malloc(size_right*sizeof(int));
      for(int i = 0; i < size_left; i++)
	ranks_left[i] = startIdx.first+i;
      for(int i = 0; i < size_right; i++)
	ranks_right[i] = middleIdx2.first+i;
      if(middleIdx1.first == middleIdx2.first)
	shizophrenicPID = middleIdx1.first;
      MPI_Group group_left, group_right;
      MPI_Group_incl(group,size_left,ranks_left,&group_left);
      MPI_Group_incl(group,size_right,ranks_right,&group_right);
      MPI_Comm_create(comm,group_left,left);
      MPI_Comm_create(comm,group_right,right);
      return make_pair(startIdx.first,middleIdx2.first);
  }
  
  void printArray(int N, int rank, int *array) {
    for(int i = 0; i < N; i++)
      cout << array[i] << "("<<rank<<") ";
    cout << endl;
  }
  
  
  pair<int,int> getProcessIdx(int start_pid,int idx) {
    return make_pair(idx/N-start_pid,idx%N);
  }
  
  int N;
  T *data;
  
};