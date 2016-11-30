/*****************************************************************************
 * This file is part of Project ShizophrenicQuicksort 
 * https://github.com/kittobi1992/ShizophrenicQuicksort.git
 * 
 * Copyright (c) 2016-2017, Tobias Heuer <tobias.heuer@gmx.net>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#include <mpi.h>
#include <bits/stdc++.h>
#include <chrono>
#include <type_traits>
#include "InputReader.cpp"
#include "../sort/QuickSort.h"

using namespace std;
using HighResClockTimepoint = std::chrono::time_point<std::chrono::high_resolution_clock>;

//Macros
#define V(X) cout << #X << "=" << X << endl
#define W(X) #X << "=" << X

template<typename T>
class QuickSortExecuter {
 
public:
  
  QuickSortExecuter(int N, MPI_Comm comm, ifstream *stream = nullptr) : N(N), comm(comm) {
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    mpi_type = SortingDatatype<T>::getMPIDatatype();
    if(rank == 0) {
      data = InputReader<T>::readInput(N,size,*stream);
    }
    int split_size = ceil(static_cast<double>(N)/size);
    pdata = (T *) malloc(split_size*sizeof(T));
    MPI_Scatter(data,split_size,mpi_type,pdata,split_size,mpi_type,0,MPI_COMM_WORLD);
  }
  
  void execute(int seed, string benchmark) {
    int split_size = ceil(static_cast<double>(N)/size);
    QuickSort<T> *qs = new QuickSort<T>(split_size,pdata);
    QSInterval<T> ival(0,0,N,comm);
    HighResClockTimepoint start = std::chrono::high_resolution_clock::now();
    qs->quickSort(ival);
    HighResClockTimepoint end = std::chrono::high_resolution_clock::now();
    T* d = qs->getData();  
    std::chrono::duration<double> start_time = start.time_since_epoch();
    std::chrono::duration<double> end_time = end.time_since_epoch();
    double d_start = start_time.count();
    double d_end = end_time.count();
    double min_start, max_end;
    MPI_Barrier(comm);
    MPI_Reduce(&d_start,&min_start,1,MPI_DOUBLE,MPI_MIN,0,comm);
    MPI_Reduce(&d_end,&max_end,1,MPI_DOUBLE,MPI_MAX,0,comm);
    
    double *runningTimes = (double *)malloc(size*sizeof(double));
    double runningTime = (d_end-d_start);
    MPI_Gather(&runningTime,1,MPI_DOUBLE,runningTimes,1,MPI_DOUBLE,0,comm);
    
    MPI_Gather(d,split_size,mpi_type,data,split_size,mpi_type,0,comm);
    if(rank == 0) {
      bool isSort = isSorted(data,N);
      cout << "RESULT "
	   << "benchmark=" << benchmark << " "
	   << "P=" << size << " "
	   << "N=" << N << " "
	   << "algo=ShizophrenicQuicksort "
	   << "seed=" << seed << " "
	   << "datatype=" << typeid(T).name() << " "
	   << "isSorted=" << isSort << " "
	   << "time=" << (max_end-min_start) << " ";
      for(int i = 0; i < size; i++)
	cout << "timePE" << i << "=" << runningTimes[i] << " ";
      cout << endl;
    }
  }
  
private:
  bool isSorted(T *data, int N) {
    for(int i = 1; i < N; i++)
      if(data[i-1] > data[i])
	return false;
    return true;
  }
  
  
  int N, size, rank;
  T *data, *pdata;
  MPI_Comm comm;
  MPI_Datatype mpi_type;
  
};