#pragma once
#include "CommonStruc.h"
#include <fstream>
#include <iostream>
#include <tuple>
#include <algorithm>
using std::tuple;
using namespace std;
class GraphBase
{
public:

	static void load_graph_directly_nbr_sorted(const string filename, Graph &O_graph, Graph &R_graph)
	{
		size_t numV, numE;
		uint32_t srcId, dstId;
		ifstream infile(filename.c_str());
		if (!infile.is_open())
		{
			cout << "The file \"" + filename + "\" can NOT be opened\n";
			return;
		}
		infile >> numV >> numE;
		O_graph.resize(numV);
		R_graph.resize(numV);
		for (auto i = numE; i--;)
		{
			infile >> srcId >> dstId;
			auto pos_src=upper_bound(O_graph[srcId].begin(), O_graph[srcId].end(), dstId);
			auto pos_dst=upper_bound(R_graph[dstId].begin(), R_graph[dstId].end(), srcId);
			if(pos_src==O_graph[srcId].end())
			{
				O_graph[srcId].push_back(dstId);
			}
			else
			{
				O_graph[srcId].insert(pos_src, dstId);
			}
			if(pos_dst==R_graph[dstId].end())
			{
				R_graph[dstId].push_back(srcId);
			}
			else
			{
				R_graph[dstId].insert(pos_dst, srcId);
			}
		}
	}
}; // cls