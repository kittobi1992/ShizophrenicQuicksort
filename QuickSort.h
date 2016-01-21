#include <mpi.h>
#include <bits/stdc++.h>
#include <future>

using namespace std;

//Macros
#define V(X) cout << #X << "=" << X << endl
#define W(X) #X << "=" << X

struct QSInterval {
  
  QSInterval(int start_pid, int global_start, int global_end, MPI_Comm comm, bool less_equal = true) :
  start_pid(start_pid), N(0), start(0), end(0),
  global_start(global_start),global_end(global_end),p(0),rank(0),
  comm(comm),
  less_equal(less_equal) {
    
    MPI_Comm_size(comm,&p);
    MPI_Comm_rank(comm,&rank);
    
  }
  
  void toString() {
    cout << "GlobalRank=" << (start_pid+rank) << ", " << W(global_start) << ", " << W(global_end) << ", "
	 << W(start) << ", " << W(end) << ", " << W(p) << ", " << W(rank) << endl;
  }
  
  
  void setDataSize(int size) {
    N = size;
    pair<int,int> startIdx = make_pair(global_start/N-start_pid,global_start%N);
    pair<int,int> endIdx = make_pair((global_end-1)/N-start_pid,(global_end-1)%N);
    if(rank != startIdx.first)
      startIdx.second = 0;
    if(rank != endIdx.first)
      endIdx.second = N-1;
    start = startIdx.second;
    end = endIdx.second;
  }
  
  void calculateExchangeData(int bound) {
    small_elements = bound; large_elements = (end-start+1) - bound;
    MPI_Scan(&small_elements,&presum_small,1,MPI_INT,MPI_SUM,comm);
    MPI_Scan(&large_elements,&presum_large,1,MPI_INT,MPI_SUM,comm);
    if(rank == p-1) {
      small_size = presum_small;
      large_size = presum_large;
    }
    MPI_Bcast(&small_size,1,MPI_INT,p-1,comm); 
    MPI_Bcast(&small_size,1,MPI_INT,p-1,comm); 
    presum_small -= small_elements; presum_large -= large_elements;
  }
  
  int start_pid, global_start,global_end, start, end;
  int N, p, rank;
  int presum_small, presum_large, small_elements, large_elements, small_size, large_size;
  MPI_Comm comm;
  bool less_equal;
};

template<typename T>
class QuickSort {
 
public:
  QuickSort(int N, T *data) : N(N), data(data) {
    buffer = (int *) malloc(N*sizeof(int));
    for(int i = 0; i < N; i++)
      buffer[i] = 0;
  }
  
  void quickSort(QSInterval ival) {
    ival.setDataSize(N);
    //ival.toString();

    if(ival.p == 1) {
      sort(data+ival.start,data+ival.end+1);
      return;
    }
        
    int bound = partition_data(ival);
    
    ival.calculateExchangeData(bound);
    
    if(ival.small_size == ival.global_end-ival.global_start) {
      ival.less_equal = false;
      quickSort(ival);
      return;
    } else if(ival.small_size == 0 && !ival.less_equal)
      return;
    
    //Exchange small elements (all elements <= pivot)
    exchangeData(ival,ival.global_start,bound,ival.presum_small,data+ival.start,buffer);
    //Exchange large elements (all elements > pivot)
    exchangeData(ival,ival.global_start+ival.small_size,ival.large_elements,ival.presum_large,data+ival.start+ival.small_elements,buffer);
    
    
    for(int i = ival.start; i < ival.end+1; i++)
      data[i] = buffer[i];
    
    MPI_Comm left, right;
    int shizophrenicPID = -1;
    pair<int,int> new_start_pid = createNewCommunicators(ival,ival.global_start+ival.small_size,shizophrenicPID,&left,&right);
    future<void> fut1,fut2;
    if(MPI_COMM_NULL != left) {
      QSInterval left_ival(new_start_pid.first+ival.start_pid,ival.global_start,ival.global_start+ival.small_size,left);
      fut1 = async(launch::async,&QuickSort::quickSort,this,left_ival); 
    }
    else
      fut1 = async([]() {return;});

    if(shizophrenicPID == ival.rank)
      fut1.get();
    
    if(MPI_COMM_NULL != right) {
      QSInterval right_ival(new_start_pid.second+ival.start_pid,ival.global_start+ival.small_size,ival.global_end,right);
      fut2 = async(launch::async,&QuickSort::quickSort,this,right_ival);
    }
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
  
  int getPivot(QSInterval &ival) {
    int c = 5;
    int pivot = 0;
    for(int i = 0; i < c; i++) {
      pivot += data[ival.start + rand() % (ival.end+1-ival.start)];
    }
    pivot /= c;
    int pivotPE;
    if(ival.rank == 0)
      pivotPE = rand() % ival.p;
    MPI_Request *req; MPI_Status *status;
    MPI_Bcast(&pivotPE,1,MPI_INT,0,ival.comm);
    MPI_Bcast(&pivot,1,MPI_INT,pivotPE,ival.comm);
    return pivot;
  }
  
  int partition_data(QSInterval &ival) {
    int bound = 0;
    int pivot = getPivot(ival);
    for(int i = ival.start; i < ival.end+1; i++) {
      if(ival.less_equal) {
	if(data[i] <= pivot) {
	  swap(data[i],data[ival.start+bound++]);
	}
      } else {
	if(data[i] < pivot) {
	  swap(data[i],data[ival.start+bound++]);
	}	
      }
    }
    return bound;
  }
  
  void exchangeData(QSInterval &ival, int start_idx, int length, int prefix_sum, int *data, int *buffer) {
    
    pair<int,int> process_idx = getProcessIdx(ival.start_pid,start_idx+prefix_sum);
    int *data_send = (int *) malloc(ival.p*sizeof(int));
    int *data_length_send = (int *) malloc(ival.p*sizeof(int));
    int *data_recv = (int *) malloc(ival.p*sizeof(int));
    int *data_length_recv = (int *) malloc(ival.p*sizeof(int));
    int send_count = 0, recv_count = 0;
    for(int i = 0; i < ival.p; i++) {
      data_length_send[i] = 0;
      if(i == process_idx.first) {
	if(i != ival.rank)
	  send_count++;
	data_send[i] = process_idx.second;
	data_length_send[i] = min(N-data_send[i],length);
	prefix_sum += data_length_send[i];
	length -= data_length_send[i];
	if(length == 0)
	  process_idx.first = -1;
	else {
	  process_idx = getProcessIdx(ival.start_pid,start_idx+prefix_sum);
	}
      } else 
	data_send[i] = -1;
    } 
    
    MPI_Barrier(ival.comm);
    
    MPI_Alltoall(data_send,1,MPI_INT,data_recv,1,MPI_INT,ival.comm);
    MPI_Alltoall(data_length_send,1,MPI_INT,data_length_recv,1,MPI_INT,ival.comm);
    
    for(int i = 0; i < ival.p; i++) {
	if(data_recv[i] != -1 && i != ival.rank)
	  recv_count++;
    }
   
    MPI_Request *send_req = (MPI_Request *) malloc(send_count*sizeof(MPI_Request));
    MPI_Request *recv_req = (MPI_Request *) malloc(recv_count*sizeof(MPI_Request));
    int send_idx = 0, recv_idx = 0;
    for(int i = 0; i < ival.p; i++) {
	if(ival.rank == i) {
	  int idx = 0;
	  for(int j = 0; j < ival.p; j++) {
	    if(data_send[j] != -1 && j != ival.rank) {
	      MPI_Isend(data+idx,data_length_send[j],MPI_INT,j,42,ival.comm,&send_req[send_idx++]);
	      idx += data_length_send[j];
	    } else if(ival.rank == j && data_send[j] != -1) {
	      for(int k = 0; k < data_length_send[ival.rank]; k++)
		*(buffer+(data_recv[ival.rank]+k)) = data[k+idx];
	      idx += data_length_send[ival.rank];
	    }
	  }
	}
	else {
	  MPI_Status status;
	  if(data_recv[i] != -1)
	    MPI_Irecv(buffer+data_recv[i],data_length_recv[i],MPI_INT,i,42,ival.comm,&recv_req[recv_idx++]);
	}
    }
    
    MPI_Status *send_status = (MPI_Status *) malloc(send_count*sizeof(MPI_Status));
    MPI_Status *recv_status = (MPI_Status *) malloc(recv_count*sizeof(MPI_Status));
    MPI_Waitall(send_count,send_req,send_status);
    MPI_Waitall(recv_count,recv_req,recv_status);
  }
  
  pair<int,int> createNewCommunicators(QSInterval &ival, int middle, int &shizophrenicPID, MPI_Comm *left, MPI_Comm *right) {
      MPI_Group group;
      MPI_Comm_group(ival.comm,&group);
      pair<int,int> startIdx = getProcessIdx(ival.start_pid,ival.global_start), middleIdx1 = getProcessIdx(ival.start_pid,middle-1), 
		    middleIdx2 = getProcessIdx(ival.start_pid,middle),endIdx = getProcessIdx(ival.start_pid,ival.global_end-1);
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
      MPI_Comm_create(ival.comm,group_left,left);
      MPI_Comm_create(ival.comm,group_right,right);
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
  T *data, *buffer;
  
};