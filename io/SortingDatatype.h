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
#include <type_traits>

using namespace std;

enum class SORT_TYPE : uint8_t {Int, Double, Long, Float};

template<typename T>
class SortingDatatype {
 
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
  
  static SORT_TYPE getSortType(string datatype) {
    if(datatype.compare("int") == 0)
      return SORT_TYPE::Int;
    else if(datatype.compare("double") == 0)
      return SORT_TYPE::Double;
    else if(datatype.compare("long") == 0)
      return SORT_TYPE::Long;
    else if(datatype.compare("float") == 0)
      return SORT_TYPE::Float;
    else
      return SORT_TYPE::Int;
  }
  
};