#include <mpi.h>
#include <bits/stdc++.h>
#include <type_traits>

using namespace std;

//Macros
#define V(X) cout << #X << "=" << X << endl
#define W(X) #X << "=" << X

template<typename T>
class InputReader {
 
public:
  
  static T * readInput(int &N, int world_size, ifstream &stream) {
    T *data; 
    int split_size;
    if(stream.is_open()) {
      split_size = ceil(static_cast<double>(N)/world_size);
      data = (T *) malloc(split_size*world_size*sizeof(T));
      for(int i = 0; i < N; i++) {
	stream >> data[i];
      }
    }
    stream.close();
    for(int i = N; i < split_size*world_size; i++)
	data[i] = numeric_limits<T>::max();
    return data;
  }
};