#include <bits/stdc++.h>
#include <chrono>

using namespace std;
using HighResClockTimepoint = std::chrono::time_point<std::chrono::high_resolution_clock>;

int main(int argc, char** argv) {
    ios::sync_with_stdio(false);
    int N; cin >> N;
    vector<int> data(N);
    for(int i = 0; i < N; i++)
      cin >> data[i];
    HighResClockTimepoint start = std::chrono::high_resolution_clock::now();
    sort(data.begin(),data.end());
    HighResClockTimepoint end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    cout << elapsed_seconds.count() << "s" << endl;
}