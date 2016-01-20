#include <bits/stdc++.h>

using namespace std;

int main(int argc, char** argv) {
    int N; cin >> N;
    vector<int> data(N);
    for(int i = 0; i < N; i++)
      cin >> data[i];
    sort(data.begin(),data.end());
}