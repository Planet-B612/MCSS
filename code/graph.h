#pragma once
#include "FileCtrl.h"
#include "serialize.h"
#include "CommonStruc.h"
#include <fstream>
#include <iostream>
#include <tuple>
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

	/// Format the input for future computing, which is much faster for loading. Vector serialization is used.
	static void format_graph(const string filename, const bool isReverse)
	{
		size_t numV, numE;
		uint32_t srcId, dstId;
		// float weight = 0.0;
		ifstream infile(filename.c_str());
		if (!infile.is_open())
		{
			cout << "The file \"" + filename + "\" can NOT be opened\n";
			return;
		}
		infile >> numV >> numE;
		if (isReverse)
		{
			Graph vecGRev(numV);
			vector<size_t> vecInDeg(numV);
			for (auto i = numE; i--;)
			{
				infile >> srcId >> dstId;
				// vecGRev[dstId].push_back(Edge(srcId, weight));
				vecGRev[dstId].push_back(srcId);
			}
			infile.close();
			for (size_t idx = 0; idx < numV; idx++)
			{
				vecInDeg[idx] = vecGRev[idx].size();
			}
			// auto idx = 0;
			// for (auto &inNbrs : vecGRev)
			// {
			// 	if (inNbrs.empty())
			// 		continue; // Skip if there is no in-neighbor.
			// 	weight = (float)1.0 / vecInDeg[idx++];
			// 	for (auto &inNbr : inNbrs)
			// 	{
			// 		get<1>(inNbr) = weight;
			// 	}
			// }
			TIO::save_graph_struct(filename, vecGRev, true);
		} // if Reverse
		else // not reverse
		{
			Graph vecG(numV);
			vector<size_t> vecInDeg(numV);
			for (auto i = numE; i--;)
			{
				infile >> srcId >> dstId;
				// pair<uint32_t, float> e = make_pair(dstId, weight);
				vecG[srcId].push_back(dstId);
				vecInDeg[dstId] += 1;
			}
			infile.close();
			// auto idx = 0;
			// for (auto &outNbrs : vecG)
			// {
			// 	if (outNbrs.empty())
			// 		continue; // Skip if there is no in-neighbor.
			// 	weight = (float)1.0 / vecInDeg[idx++];
			// 	for (auto &outNbr : outNbrs)
			// 	{
			// 		get<1>(outNbr) = weight;
			// 		//get<1>(outNbr) = (float)1.0 / vecInDeg[get<0>(outNbr)];
			// 	}
			// }
			TIO::save_graph_struct(filename, vecG, false);
		}

		cout << "The graph is formatted!" << endl;
	}

/*
	/// Load graph via vector deserialization.
	static Graph load_graph(const string graphName, const bool isReverse)
	{
		Graph graph;
		TIO::load_graph_struct(graphName, graph, isReverse);

		if (isReverse)
		{
			// Reverse graph
			for (auto &nbrs : graph)
			{
				for (auto &nbr : nbrs)
				{
					get<1>(nbr) = float(1.0 / nbrs.size());
				}
			}
		}
		else
		{
			// Forward graph
			vector<uint32_t> vecInDeg(graph.size());
			for (auto &nbrs : graph)
			{
				for (auto &nbr : nbrs)
				{
					vecInDeg[get<0>(nbr)]++;
				}
			}
			for (auto &nbrs : graph)
			{
				for (auto &nbr : nbrs)
				{
					get<1>(nbr) = (float)1.0 / vecInDeg[get<0>(nbr)];
				}
			}
		}
		return graph;
	}
*/

}; // cls