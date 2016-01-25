#include <mpi.h>
#include <bits/stdc++.h>
#include "../io/SortingDatatype.h"

using namespace std;

//Macros
#define V(X) cout << #X << "=" << X << endl
#define W(X) #X << "=" << X

template<typename T>
struct QSInterval {
  
  QSInterval(int start_pid, int global_start, int global_end, MPI_Comm comm, bool less_equal = true) :
  start_pid(start_pid), N(0), start(0), end(0),
  global_start(global_start),global_end(global_end),p(0),rank(0),
  comm(comm),
  less_equal(less_equal) {
    
    MPI_Comm_size(comm,&p);
    MPI_Comm_rank(comm,&rank);
    
  }
  
  void toString(bool shizophren) {
    cout << "GlobalRank=" << (start_pid+rank) << ", " << W(global_start) << ", " << W(global_end) << ", "
	 << W(start) << ", " << W(end) << ", " << W(p) << ", " << W(rank) << ", "
	 << W(pivot) << ", " << W(shizophren) << endl;
  }
  
  void sortData(T *data) {
    sort(data+start,data+end+1);    
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
    MPI_Bcast(&large_size,1,MPI_INT,p-1,comm); 
    presum_small -= small_elements; presum_large -= large_elements;
  }
  
  bool allDataAreEqualOnInterval(T *data) {
      if(small_size == (global_end-global_start)) {
	int i = end+1;
	while(--i > start && data[i] == data[start]) {}
	bool ivalEqual = (i == start); bool allEqual;
	MPI_Reduce(&ivalEqual,&allEqual,1,MPI_C_BOOL,MPI_LAND,0,comm);
	MPI_Bcast(&allEqual,1,MPI_C_BOOL,0,comm);
	if(allEqual) {
	  T d = data[start];
	  T min_d, max_d;
	  MPI_Reduce(&d,&min_d,1,SortingDatatype<T>::getMPIDatatype(),MPI_MIN,0,comm);
	  MPI_Reduce(&d,&max_d,1,SortingDatatype<T>::getMPIDatatype(),MPI_MAX,0,comm);
	  if(rank == 0) {
	    allEqual = (min_d == max_d);
	    MPI_Bcast(&allEqual,1,MPI_C_BOOL,0,comm);
	  }
	  else
	    MPI_Bcast(&allEqual,1,MPI_C_BOOL,0,comm);
	  if(allEqual) {
	    return true;
	  }
	}
	else
	  return false;
      }
      else
	return false;
  }
  
  int start_pid, global_start,global_end, start, end;
  int N, p, rank;
  T pivot;
  int presum_small, presum_large, small_elements, large_elements, small_size, large_size;
  MPI_Comm comm;
  bool less_equal;
};

template<typename T>
class QuickSort {
 
public:
  QuickSort(int N, T *data) : N(N), data(data) {
    buffer = (T *) malloc(N*sizeof(T));
    for(int i = 0; i < N; i++)
      buffer[i] = 0;
    mpi_type = SortingDatatype<T>::getMPIDatatype();
  }
  
  void quickSort(QSInterval<T> &ival) {
    ival.setDataSize(N);
    
    if(ival.global_end-ival.global_start < 10) {
      global_sort(ival);
      return;
    }

    if(ival.p == 1) {
      ival.sortData(data);
      return;
    }

    int bound = partition_data(ival);
    
    ival.calculateExchangeData(bound);
    
    if(ival.allDataAreEqualOnInterval(data)) {
      return;
    }
    
        
    ival.toString(false);
        
    
    //Exchange small elements (all elements <= pivot)
    exchangeData(ival,ival.global_start,
		 ival.small_elements,
		 ival.presum_small,
		 data+ival.start);
    //Exchange large elements (all elements > pivot)
    exchangeData(ival,ival.global_start+ival.small_size,
		 ival.large_elements,
		 ival.presum_large,
		 data+ival.start+ival.small_elements);
    
    for(int i = ival.start; i < ival.end+1; i++)
      data[i] = buffer[i];
    
    
    MPI_Comm left, right;
    pair<int,int> new_start_pid = createNewCommunicators(ival,ival.global_start+ival.small_size,&left,&right);
    
    if(MPI_COMM_NULL != left && MPI_COMM_NULL != right) {
      QSInterval<T> left_ival(new_start_pid.first+ival.start_pid,ival.global_start,ival.global_start+ival.small_size,left);
      QSInterval<T> right_ival(new_start_pid.second+ival.start_pid,ival.global_start+ival.small_size,ival.global_end,right);
      shizophrenicQuickSort(left_ival,right_ival);
    } else if(MPI_COMM_NULL != left) {
      QSInterval<T> left_ival(new_start_pid.first+ival.start_pid,ival.global_start,ival.global_start+ival.small_size,left);
      quickSort(left_ival);
    } else if(MPI_COMM_NULL != right) {
      QSInterval<T> right_ival(new_start_pid.second+ival.start_pid,ival.global_start+ival.small_size,ival.global_end,right);
      quickSort(right_ival);
    }


  }
  
   T* getData() {
    return data;
  }
  
private:
  

  void shizophrenicQuickSort(QSInterval<T> &left_ival, QSInterval<T> &right_ival) {
     left_ival.setDataSize(N);
     right_ival.setDataSize(N);

    if(left_ival.p == 1) {
      left_ival.sortData(data);
      if(right_ival.p > 1) {
	quickSort(right_ival);
	return;
      } else {
	right_ival.sortData(data);
	return;
      }
    }
    if(right_ival.p == 1) {
      right_ival.sortData(data);
      quickSort(left_ival);
      return;
    }
    
    if(left_ival.global_end-left_ival.global_start < 10) {
      global_sort(left_ival); 
      quickSort(right_ival);
      return;
    }
    
    if(right_ival.global_end-right_ival.global_start < 10) {
      global_sort(right_ival); 
      quickSort(left_ival);
      return;
    }
    
    
    
        
    int bound1 = partition_data(left_ival);
    int bound2 = partition_data(right_ival);
    
    left_ival.calculateExchangeData(bound1);
    right_ival.calculateExchangeData(bound2);
    
    bool dataEqualLeft = left_ival.allDataAreEqualOnInterval(data);
    bool dataEqualRight = right_ival.allDataAreEqualOnInterval(data);
    
    left_ival.toString(true);
    right_ival.toString(true);
    
    
    //Exchange small elements for left Interval (all elements <= pivot)
    if(!dataEqualLeft)
      exchangeData(left_ival, left_ival.global_start, 
		  left_ival.small_elements,
		  left_ival.presum_small,
		  data+left_ival.start);
    //Exchange small elements for right Interval (all elements <= pivot)
    if(!dataEqualRight)
      exchangeData(right_ival, right_ival.global_start, 
		  right_ival.small_elements,
		  right_ival.presum_small,
		  data+right_ival.start);
    //Exchange large elements for left Interval (all elements > pivot)
    if(!dataEqualLeft)
      exchangeData(left_ival,left_ival.global_start+left_ival.small_size, 
		  left_ival.large_elements,
		  left_ival.presum_large,
		  data+left_ival.start+left_ival.small_elements);
    //Exchange large elements (all elements > pivot)
    if(!dataEqualRight)
      exchangeData(right_ival,right_ival.global_start+right_ival.small_size, 
		  right_ival.large_elements,
		  right_ival.presum_large,
		  data+right_ival.start+right_ival.small_elements);
    
    if(!dataEqualLeft)
      for(int i = left_ival.start; i < left_ival.end+1; i++)
	data[i] = buffer[i];
    if(!dataEqualRight)
      for(int i = right_ival.start; i < right_ival.end+1; i++)
	data[i] = buffer[i];
      
    MPI_Comm left1, right1, left2, right2;
    pair<int,int> new_start_pid1, new_start_pid2;
    if(!dataEqualLeft)
      new_start_pid1 = createNewCommunicators(left_ival,
					      left_ival.global_start + left_ival.small_size,
					      &left1,&right1);
    if(!dataEqualRight)
      new_start_pid2 = createNewCommunicators(right_ival,
					      right_ival.global_start + right_ival.small_size,
					      &left2,&right2);
      
      
    if(dataEqualLeft && dataEqualRight) {
      return;
    } else if(dataEqualLeft) {
      if(left2 != MPI_COMM_NULL && right2 != MPI_COMM_NULL) {
	QSInterval<T> left_ival2(new_start_pid2.first+right_ival.start_pid,right_ival.global_start,right_ival.global_start+right_ival.small_size,left2);
	QSInterval<T> right_ival2(new_start_pid2.second+right_ival.start_pid,right_ival.global_start+right_ival.small_size,right_ival.global_end,right2);
	shizophrenicQuickSort(left_ival2,right_ival2);
      } else if(left2 == MPI_COMM_NULL && right2 != MPI_COMM_NULL) {
	QSInterval<T> right_ival2(new_start_pid2.second+right_ival.start_pid,right_ival.global_start+right_ival.small_size,right_ival.global_end,right2);
	quickSort(right_ival2);
      } else if(left2 != MPI_COMM_NULL && right2 == MPI_COMM_NULL) {
	QSInterval<T> left_ival2(new_start_pid2.first+right_ival.start_pid,right_ival.global_start,right_ival.global_start+right_ival.small_size,left2);
	quickSort(left_ival2);
      }
      return;
    } else if(dataEqualRight) {
      if(left1 != MPI_COMM_NULL && right1 != MPI_COMM_NULL) {
	QSInterval<T> left_ival1(new_start_pid1.first+left_ival.start_pid,left_ival.global_start,left_ival.global_start+left_ival.small_size,left1);
	QSInterval<T> right_ival1(new_start_pid1.second+left_ival.start_pid,left_ival.global_start+left_ival.small_size,left_ival.global_end,right1);
	shizophrenicQuickSort(left_ival1,right_ival1);
      } else if(left1 == MPI_COMM_NULL && right1 != MPI_COMM_NULL) {
	QSInterval<T> right_ival1(new_start_pid1.second+left_ival.start_pid,left_ival.global_start+left_ival.small_size,left_ival.global_end,right1);
	quickSort(right_ival1);
      } else if(left1 != MPI_COMM_NULL && right1 == MPI_COMM_NULL) {
	QSInterval<T> left_ival1(new_start_pid1.first+left_ival.start_pid,left_ival.global_start,left_ival.global_start+left_ival.small_size,left1);
	quickSort(left_ival1);
      }  
      return;
    }

     //cout << (left1 != MPI_COMM_NULL) << " " << (right1 != MPI_COMM_NULL) << " " <<  (left2 != MPI_COMM_NULL) << " " << (right2 != MPI_COMM_NULL) << endl;
    
    if(left1 != MPI_COMM_NULL && right1 != MPI_COMM_NULL && left2 != MPI_COMM_NULL && right2 != MPI_COMM_NULL) {
      QSInterval<T> left_ival1(new_start_pid1.first+left_ival.start_pid,left_ival.global_start,left_ival.global_start+left_ival.small_size,left1);
      QSInterval<T> right_ival1(new_start_pid1.second+left_ival.start_pid,left_ival.global_start+left_ival.small_size,left_ival.global_end,right1);
      QSInterval<T> left_ival2(new_start_pid2.first+right_ival.start_pid,right_ival.global_start,right_ival.global_start+right_ival.small_size,left2);
      QSInterval<T> right_ival2(new_start_pid2.second+right_ival.start_pid,right_ival.global_start+right_ival.small_size,right_ival.global_end,right2);
      shizophrenicQuickSort(right_ival1,left_ival2);
      shizophrenicQuickSort(left_ival1,right_ival2);
      return;
    } else if(left1 == MPI_COMM_NULL && right1 != MPI_COMM_NULL && left2 != MPI_COMM_NULL && right2 != MPI_COMM_NULL) {
      QSInterval<T> right_ival1(new_start_pid1.second+left_ival.start_pid,left_ival.global_start+left_ival.small_size,left_ival.global_end,right1);
      QSInterval<T> left_ival2(new_start_pid2.first+right_ival.start_pid,right_ival.global_start,right_ival.global_start+right_ival.small_size,left2);
      QSInterval<T> right_ival2(new_start_pid2.second+right_ival.start_pid,right_ival.global_start+right_ival.small_size,right_ival.global_end,right2);
      quickSort(left_ival2);
      shizophrenicQuickSort(right_ival1,right_ival2);
      return;
    } else if(left1 != MPI_COMM_NULL && right1 != MPI_COMM_NULL && left2 != MPI_COMM_NULL && right2 == MPI_COMM_NULL) {
      QSInterval<T> left_ival1(new_start_pid1.first+left_ival.start_pid,left_ival.global_start,left_ival.global_start+left_ival.small_size,left1);
      QSInterval<T> right_ival1(new_start_pid1.second+left_ival.start_pid,left_ival.global_start+left_ival.small_size,left_ival.global_end,right1);
      QSInterval<T> left_ival2(new_start_pid2.first+right_ival.start_pid,right_ival.global_start,right_ival.global_start+right_ival.small_size,left2);
      quickSort(right_ival1);
      shizophrenicQuickSort(left_ival1,left_ival2);
      return;
    } else if(left1 == MPI_COMM_NULL && right1 != MPI_COMM_NULL && left2 != MPI_COMM_NULL && right2 == MPI_COMM_NULL) {
      QSInterval<T> right_ival1(new_start_pid1.second+left_ival.start_pid,left_ival.global_start+left_ival.small_size,left_ival.global_end,right1);
      QSInterval<T> left_ival2(new_start_pid2.first+right_ival.start_pid,right_ival.global_start,right_ival.global_start+right_ival.small_size,left2);
      shizophrenicQuickSort(right_ival1,left_ival2);
      return;
    }

  }
  
  T getPivot(QSInterval<T> &ival) {
    double c = 5.0;
    double pivot = 0.0;
    for(int i = 0; i < c; i++) {
      pivot += static_cast<double>((data[ival.start + rand() % (ival.end+1-ival.start)]))/c;
    }
    pivot /= ival.p;
    double global_pivot;
    MPI_Reduce(&pivot,&global_pivot,1,MPI_DOUBLE,MPI_SUM,0,ival.comm);
    MPI_Bcast(&global_pivot,1,MPI_DOUBLE,0,ival.comm);
    ival.pivot = global_pivot;
    return static_cast<T>(global_pivot);
  }
  
  int partition_data(QSInterval<T> &ival) {
    int bound = 0;
    T pivot = getPivot(ival);
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
  
  void exchangeData(QSInterval<T> &ival, int start_idx, int length, int prefix_sum, T *data) {
    
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
	      MPI_Isend(data+idx,data_length_send[j],mpi_type,j,42,ival.comm,&send_req[send_idx++]);
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
	    MPI_Irecv(buffer+data_recv[i],data_length_recv[i],mpi_type,i,42,ival.comm,&recv_req[recv_idx++]);
	}
    }
    
    MPI_Status *send_status = (MPI_Status *) malloc(send_count*sizeof(MPI_Status));
    MPI_Status *recv_status = (MPI_Status *) malloc(recv_count*sizeof(MPI_Status));
    MPI_Waitall(send_count,send_req,send_status);
    MPI_Waitall(recv_count,recv_req,recv_status);
  }
  
  pair<int,int> createNewCommunicators(QSInterval<T> &ival, int middle, MPI_Comm *left, MPI_Comm *right) {
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
      MPI_Group group_left, group_right;
      MPI_Group_incl(group,size_left,ranks_left,&group_left);
      MPI_Group_incl(group,size_right,ranks_right,&group_right);
      MPI_Comm_create(ival.comm,group_left,left);
      MPI_Comm_create(ival.comm,group_right,right);
      return make_pair(startIdx.first,middleIdx2.first);
  }
  
  void global_sort(QSInterval<T> &ival) {
    ival.sortData(data);
    int idx = ival.global_start;
    int local_idx = ival.start;
    T *data_recv = (T *) malloc(ival.p*sizeof(T));
    while(idx < ival.global_end) {
      T d = (local_idx <= ival.end ? data[local_idx] : numeric_limits<T>::max());
      pair<int,int> cur_proc = getProcessIdx(ival.start_pid,idx);
      MPI_Gather(&d,1,mpi_type,data_recv,1,mpi_type,cur_proc.first,ival.comm);
      int min_proc = -1;
      if(ival.rank == cur_proc.first) {
	T min_d = numeric_limits<T>::max();
	for(int i = 0; i < ival.p; i++) {
	  if(min_d > data_recv[i]) {
	   min_d = data_recv[i];
	   min_proc = i;
	  }
	}
	buffer[cur_proc.second] = min_d;
      }
      MPI_Bcast(&min_proc,1,mpi_type,cur_proc.first,ival.comm);
      idx++;
      if(ival.rank == min_proc)
	local_idx++;
    }
    for(int i = ival.start; i < ival.end+1; i++)
      data[i] = buffer[i];
  }
  
  void printArray(int rank, T *array) {
    for(int i = 0; i < N; i++)
      cout << array[i] << "("<<rank<<") ";
    cout << endl;
  }
  
  
  pair<int,int> getProcessIdx(int start_pid,int idx) {
    return make_pair(idx/N-start_pid,idx%N);
  }
  
  int N;
  T *data, *buffer;
  MPI_Datatype mpi_type;
  
};