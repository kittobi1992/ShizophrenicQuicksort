#include <bits/stdc++.h>

using namespace std;

int main() {
  
  int N = 1000; int MAX_NUM = 500;
  cout << N << endl;
  for(int i = 0; i < N; i++) {
    cout << rand() % MAX_NUM << (i == N ? "" : " ");
  }
  cout << endl;
  
  return 0;
}