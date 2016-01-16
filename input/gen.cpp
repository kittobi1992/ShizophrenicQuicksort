#include <bits/stdc++.h>

using namespace std;

int main() {
  
  int N = 40; int MAX_NUM = 20;
  cout << N << endl;
  for(int i = 0; i < N; i++) {
    cout << rand() % MAX_NUM << (i == N ? "" : " ");
  }
  cout << endl;
  
  return 0;
}