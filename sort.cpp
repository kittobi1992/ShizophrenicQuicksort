/*****************************************************************************
 * This file is part of Project ShizophrenicQuicksort 
 * https://github.com/kittobi1992/ShizophrenicQuicksort.git
 * 
 * Copyright (c) 2016-2017, Tobias Heuer <tobias.heuer@gmx.net>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#include <bits/stdc++.h>
#include <chrono>

using namespace std;
using HighResClockTimepoint = std::chrono::time_point<std::chrono::high_resolution_clock>;

bool isSorted(int *data, int N) {
    for(int i = 1; i < N; i++)
      if(data[i-1] > data[i])
	return false;
    return true;
}

int main(int argc, char** argv) {
    ios::sync_with_stdio(false);
    string benchmark(argv[1]);
    ifstream stream(benchmark);
    int N; stream >> N;
    string datatype; stream >> datatype;
    int *data = (int *) malloc(N*sizeof(int));
    for(int i = 0; i < N; i++)
      stream >> data[i];
    stream.close();
    HighResClockTimepoint start = std::chrono::high_resolution_clock::now();
    sort(data,data+N);
    HighResClockTimepoint end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    cout << "RESULT benchmark=" << benchmark
	 << " algo=std_sort" 
	 << " N=" << N
	 << " datatype=" << typeid(int).name()
         << " time=" << elapsed_seconds.count() 
	 << " isSorted=" << isSorted(data,N) << endl;
}