#include <cassert>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <stdlib.h>

using namespace std;

int main(int argc, char *argv[]) {
  if (argc != 5) {
    cout << "error in parameter: input_graph_path output_graph_path n m" << endl;
    return 0;
  }

  ifstream input(argv[1]);

  if (!input.is_open()) {
    cout << "error in opening file " << argv[1] << endl;
    return 0;
  }

  int n;
  long long m;
  n = atoi(argv[3]);
  m=atoll(argv[4]);
  //input >> n >> m;
  cout << "n = " << n << ", m = " << m << endl;

  vector<vector<int>> graph;
  graph.resize(n);
  int s, t;

  long long verify_m = 0;
  vector<int> indeg(n, 0);

  while (input >> s >> t) {
    assert(s < n);
    assert(t < n);
    graph[s].push_back(t);
    indeg[t]++;
    verify_m++;
  }
  assert(verify_m == m);

  vector<double> prob(n);
  for (int node = 0; node < n; node++) {
    prob[node] = 1.0 / indeg[node];
  }

  FILE *ofile = fopen(argv[2], "w");
  if (!ofile)
  {
    cout << "error in opening file " << argv[2] << endl;
    return 0;
  }

  ofstream text("test.txt");
  if (!text.is_open())
  {
    cout << "error in opening file " << "test.txt" << endl;
    return 0;
  }

  for (long long i = 0; i < (long long)graph.size(); i++) {
    auto nbrs = graph[i];
    if (i % 10000 == 0) {
      cout << i << endl;
    }
    for (int j = 0; j < (int)nbrs.size(); j++)
    {
      int target = nbrs[j];
      fwrite(&i, sizeof(int), 1, ofile);
      fwrite(&target, sizeof(int), 1, ofile);
      double p = prob[target];
      fwrite(&p, sizeof(double), 1, ofile);
      //text << i << " " << target << " " << p << endl;
    }
  }

  text.close();
  fclose(ofile);
  return 0;
}