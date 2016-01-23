#include <mpi.h>
#include <bits/stdc++.h>
#include <type_traits>

using namespace std;

//Macros
#define V(X) cout << #X << "=" << X << endl
#define W(X) #X << "=" << X


template<typename T>
class GenericDatatype {
 
public:
  
  static MPI_Datatype getMPIDatatype() {
    if(is_same<T,int>::value)
      return MPI_INT;
    else if(is_same<T,short>::value)
      return MPI_SHORT;
    else if(is_same<T,char>::value)
      return MPI_CHAR;
    else if(is_same<T,long long>::value)
      return MPI_LONG_LONG;  
    else if(is_same<T,float>::value)
      return MPI_FLOAT;
    else if(is_same<T,double>::value)
      return MPI_DOUBLE;
  }
  
};