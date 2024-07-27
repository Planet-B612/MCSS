#pragma once
#include "./dSFMT/dSFMT.h"
#include "CommonStruc.h"
#include <memory>
#include <queue>
#include "CommonFunc.h"
using namespace std;
class RRcollection
{
private:
	/// __numV: number of nodes in the graph.
	uint32_t __numV;
	/// __numE: number of edges in the graph.
	size_t __numE = 0;
	/// __numRRsets: number of RR sets.
	size_t __numRRsets = 0;
	std::vector<bool> __vecVisitBool;
	Nodelist __vecVisitNode;
	Nodelist __roots;
	uint8_t __mode=0;

	/// Initialization
	void init_RRcollection()
	{
		__numV = (uint32_t)_graph.size();  
		for (auto& nbrs : _graph) __numE += nbrs.size();
		_FRsets = FRsets(__numV);
		__vecVisitBool = std::vector<bool>(__numV);
		__vecVisitNode = Nodelist(__numV);
	}

public:
	Graph & _graph; 
	FRsets _FRsets;
	RRsets _RRsets;
	FRcollection _MFFcol;
	FRcollection _GFFcol;

	CascadeModel _cascadeModel = IC;

	explicit RRcollection(Graph & graph) : _graph(graph) 
	{
		init_RRcollection();
	}

	/// Set estimation mode
	void set_model_mode(const CascadeModel model, const uint8_t mode)
	{
		_cascadeModel=model;
		__mode=mode;
	}
	/// Set cascade model
	void set_cascade_model(const CascadeModel model)
	{
		_cascadeModel = model;
	}

	/// Returns the number of nodes in the graph.
	uint32_t get_nodes() const
	{
		return __numV;
	}

	/// Returns the number of edges in the graph.
	size_t get_edges() const
	{
		return __numE;
	}

	/// Returns the number of RR sets in the graph.
	size_t get_RR_sets_size() const
	{
		return __numRRsets;
	}

	/// Get out degree in the original graph
	std::vector<size_t> get_out_degree() const
	{
		std::vector<size_t> outDeg(__numV);
		for (auto& nbrs : _graph)
		{
			for (auto& nbr : nbrs)
			{
				outDeg[get<0>(nbr)]++;
			}
		}
		return outDeg;
	}

	/// Generate a set of n RR sets
	void build_n_RRsets(const size_t numSamples)
	{
		if (numSamples > SIZE_MAX)
		{
			std::cout << "Error: the number of RRsets you need is too large" << std::endl;
			exit(1);
		}
		const auto prevSize = __numRRsets;
		__numRRsets = __numRRsets > numSamples ? __numRRsets : numSamples;
		for (auto i = prevSize; i < numSamples; i++)
		{
			build_one_RRset(dsfmt_gv_genrand_uint32_range(__numV), i);
		}
	}

	/// Generate one mRR set
	void build_one_RRset(const uint32_t uStart, const size_t RRsetIdx)
	{
		size_t numVisitNode = 0, currIdx = 0;  
		_FRsets[uStart].push_back(RRsetIdx);
		__vecVisitNode[numVisitNode++] = uStart;
		__vecVisitBool[uStart] = true;
		while (currIdx < numVisitNode)
		{
			const auto expand = __vecVisitNode[currIdx++];
			if (_cascadeModel == IC)
			{
				for (auto& nbr : _graph[expand])
				{
					const auto nbrId = get<0>(nbr);
					if (__vecVisitBool[nbrId])
						continue;
					const auto randDouble = dsfmt_gv_genrand_open_close();
					if (randDouble > get<1>(nbr))
						continue;
					__vecVisitNode[numVisitNode++] = nbrId;
					__vecVisitBool[nbrId] = true;
					_FRsets[nbrId].push_back(RRsetIdx);
				}
			}
			else if (_cascadeModel == LT)
			{
				if (_graph[expand].size() == 0)
					continue;
				int index = dsfmt_gv_genrand_uint32_range(_graph[expand].size());
				int v=get<0>(_graph[expand][index]);
				if(__vecVisitBool[v])
					continue;
				__vecVisitNode[numVisitNode++] = v;
				__vecVisitBool[v] = true;
				_FRsets[v].push_back(RRsetIdx);
			}
		}
		for (size_t i = 0; i < numVisitNode; i++) 
		{
			__vecVisitBool[__vecVisitNode[i]] = false;
		}
		_RRsets.push_back(RRset(__vecVisitNode.begin(), __vecVisitNode.begin() + numVisitNode));
	}

	void generate_MFFcol_20231029(size_t realizations, size_t numRR_real)
	{
		uint32_t nodeId;
		ASSERTT(realizations>0 && numRR_real>0, "The number of realizations or roots is wrong!");
		for(uint32_t phi=1; phi<realizations+1; phi++)
		{
			queue<uint32_t> Que;
			vector<uint32_t> roots;
			FRsets _FRsets_real(__numV);
			uint32_t RRIdx=0; 
			vector<uint32_t> selected_edge(__numV,__numV+1);
			vector<uint32_t> phi_round(__numV, 0);
			for (uint32_t i = 1; i < numRR_real+1; i++)
			{
				size_t source = dsfmt_gv_genrand_uint32_range(__numV);
				roots.push_back(source);
			}
			for(auto root:roots)
			{
				_FRsets_real[root].push_back(RRIdx);  
				vector<bool> added(__numV, false);  
				added[root]=true;
				Que.push(root);
				// BFS traversal
				if (_cascadeModel == IC)
				{
					while (!Que.empty())
					{
						nodeId = Que.front();
						Que.pop();
						for (auto& nbr : _graph[nodeId])
						{
							if(get<2>(nbr)==phi) 
							{
								if(get<3>(nbr)==phi)
								{
									if(added[get<0>(nbr)])   
										{   continue;   }
									else
									{
										_FRsets_real[get<0>(nbr)].push_back(RRIdx);  
										Que.push(get<0>(nbr));
										added[get<0>(nbr)]=true;
									}
								}
								else  // The edge status is dead.
								{
									continue;
								}
							}
							else  // The edge has not been tested. Thus, we test it here.
							{
								if (dsfmt_gv_genrand_open_close() <= get<1>(nbr))  // This edge is tested to be live, and the neighbor is activated.
								{
									get<3>(nbr)=phi;
									if(!added[get<0>(nbr)]) 
									{
										_FRsets_real[get<0>(nbr)].push_back(RRIdx);
										Que.push(get<0>(nbr));
										added[get<0>(nbr)]=true;
									}
								}
								get<2>(nbr)=phi; 
							}
						}
					}
				}
				else if (_cascadeModel == LT) 
				{
					while (!Que.empty())
					{
						nodeId = Que.front();
						Que.pop();
						
						if (phi_round[nodeId]==phi) // the node has selected its neighbor in this realization
						{
							auto nbr=selected_edge[nodeId];
							if(added[nbr]) continue;							
							Que.push(nbr);
							_FRsets_real[nbr].push_back(RRIdx);
							added[nbr]=true;
						}
						else
						{
							if (_graph[nodeId].size() == 0)  continue;
							uint32_t index = dsfmt_gv_genrand_uint32_range(_graph[nodeId].size());
							auto nbr=get<0>(_graph[nodeId][index]);
							Que.push(nbr);
							_FRsets_real[nbr].push_back(RRIdx);
							added[nbr]=true;
							phi_round[nodeId]=phi;
							selected_edge[nodeId]=nbr;
						}
					}
				}
				vector<bool>().swap(added);
				++RRIdx;
			}
			_MFFcol.push_back(_FRsets_real);
			for (auto i = 0; i < _FRsets_real.size();i--)
			{
				FRset().swap(_FRsets_real[i]);
			}			
			FRsets().swap(_FRsets_real);
			vector<uint32_t>().swap(roots);
		}
	}

	void generate_GRR(size_t numGRR, size_t num_gRR)
	{
		ASSERTT(numGRR>0 && num_gRR>0, "The number of realizations or roots is wrong!");
		for (size_t i = 0; i < numGRR; i++)
		{
			FRsets GFFsets(__numV); // Equivalent to a realization
			uint32_t RRIdx=0;
			for (size_t j = 0; j < num_gRR; j++)
			{
				size_t numVisitNode = 0, currIdx = 0; 
				uint32_t root=dsfmt_gv_genrand_uint32_range(__numV);
				GFFsets[root].push_back(RRIdx);
				__vecVisitNode[numVisitNode++] = root;
				__vecVisitBool[root] = true;
				while (currIdx < numVisitNode)
				{
					const auto expand = __vecVisitNode[currIdx++];
					if (_cascadeModel == IC)
					{
						for (auto& nbr : _graph[expand])
						{
							const auto nbrId = get<0>(nbr);
							if (__vecVisitBool[nbrId])
								continue;
							const auto randDouble = dsfmt_gv_genrand_open_close();
							if (randDouble > get<1>(nbr))
								continue;
							__vecVisitNode[numVisitNode++] = nbrId;
							__vecVisitBool[nbrId] = true;
							GFFsets[nbrId].push_back(RRIdx);
						}
					}
					else if (_cascadeModel == LT)
					{
						if (_graph[expand].size() == 0)
							continue;
						int index = dsfmt_gv_genrand_uint32_range(_graph[expand].size());
						int v=get<0>(_graph[expand][index]);
						if(__vecVisitBool[v])
							continue;
						__vecVisitNode[numVisitNode++] = v;
						__vecVisitBool[v] = true;
						GFFsets[v].push_back(RRIdx);
					}
				}
				//cout<<endl;
				for (size_t i = 0; i < numVisitNode; i++) 
				{
					__vecVisitBool[__vecVisitNode[i]] = false;
				}
				RRIdx++;
			}
			_GFFcol.push_back(GFFsets);
			for(auto k=0;k<GFFsets.size();k++)
			{
				FRset().swap(GFFsets[k]);
			}
			FRsets().swap(GFFsets);
		}
	}

	/// Refresh the RRsets
	void refresh_RRsets()
	{
		for (auto i = __numRRsets; i--;)
		{
			RRset().swap(_RRsets[i]);
		}
		RRsets().swap(_RRsets);
		for (auto i = __numV; i--;)
		{
			FRset().swap(_FRsets[i]);
		}
		__numRRsets = 0;
	}

	/// Refresh MRR-sets
	void refresh_MFF()
	{
		for(auto i=_MFFcol.size();i--;)
		{
			for(auto j=_MFFcol[i].size();j--;)
			{
				FRset().swap(_MFFcol[i][j]);
			}
			FRsets().swap(_MFFcol[i]);
		}
		FRcollection().swap(_MFFcol);
	}

	void refresh_GFF()
	{
		for(auto i=_GFFcol.size();i--;)
		{
			for(auto j=_GFFcol[i].size();j--;)
			{
				FRset().swap(_GFFcol[i][j]);
			}
			FRsets().swap(_GFFcol[i]);
		}
		FRcollection().swap(_GFFcol);
	}

	/// Release memory
	void release_memory()
	{
		refresh_RRsets();
		refresh_MFF();
		refresh_GFF();
		std::vector<bool>().swap(__vecVisitBool);
		Nodelist().swap(__vecVisitNode);
		FRsets().swap(_FRsets);
		Nodelist().swap(__roots);
	}

	/// Generate one node with probabilities according to their weights for the LT cascade model
	static inline size_t gen_random_node_by_weight_LT(const Edgelist& edges)
	{
		const double weight = dsfmt_gv_genrand_open_close();
		size_t minIdx = 0, maxIdx = edges.size() - 1;
		if (weight < get<1>(edges.front())) return 0; // First element
		if (weight > get<1>(edges.back())) return edges.size() + 1; // No element
		while (maxIdx > minIdx)
		{
			const size_t meanIdx = (minIdx + maxIdx) / 2;
			const auto meanWeight = get<1>(edges[meanIdx]);
			if (weight <= meanWeight) maxIdx = meanIdx;
			else minIdx = meanIdx + 1;
		}
		return maxIdx;
	}

	/// Efficiently estimate the influence spread with sampling error epsilon within probability 1-delta
	double effic_inf_valid_algo(const Nodelist& vecSeed, const double delta = 1e-3, const double eps = 0.01)
	{
		const double c = 2.0 * (exp(1.0) - 2.0);
		const double LambdaL = 1.0 + 2.0 * c * (1.0 + eps) * log(2.0 / delta) / (eps * eps);
		size_t numHyperEdge = 0;
		size_t numCoverd = 0;
		std::vector<bool> vecBoolSeed(__numV);
		for (auto seed : vecSeed) vecBoolSeed[seed] = true;

		while (numCoverd < LambdaL)
		{
			numHyperEdge++;
			size_t numVisitNode = 0, currIdx = 0;
			const auto uStart = dsfmt_gv_genrand_uint32_range(__numV);
			if (vecBoolSeed[uStart])
			{
				// Stop, this sample is covered
				numCoverd++;
				continue;
			}
			__vecVisitNode[numVisitNode++] = uStart;
			__vecVisitBool[uStart] = true;
			while (currIdx < numVisitNode)
			{
				const auto expand = __vecVisitNode[currIdx++];
				if (_cascadeModel == IC)
				{
					for (auto& nbr : _graph[expand])
					{
						const auto nbrId = get<0>(nbr);
						if (__vecVisitBool[nbrId])
							continue;
						const auto randDouble = dsfmt_gv_genrand_open_close();
						if (randDouble > get<1>(nbr))
							continue;
						if (vecBoolSeed[nbrId])
						{
							// Stop, this sample is covered
							numCoverd++;
							goto postProcess;
						}
						__vecVisitNode[numVisitNode++] = nbrId;
						__vecVisitBool[nbrId] = true;
					}
				}
				else if (_cascadeModel == LT)
				{
					if (_graph[expand].size() == 0)
						continue;
					const auto nextNbrIdx = gen_random_node_by_weight_LT(_graph[expand]);
					if (nextNbrIdx >= _graph[expand].size()) break; // No element activated
					const auto nbrId = get<0>(_graph[expand][nextNbrIdx]);
					if (__vecVisitBool[nbrId]) break; // Stop, no further node activated
					if (vecBoolSeed[nbrId])
					{
						// Stop, this sample is covered
						numCoverd++;
						goto postProcess;
					}
					__vecVisitNode[numVisitNode++] = nbrId;
					__vecVisitBool[nbrId] = true;
				}
			}
		postProcess:
			for (size_t i = 0; i < numVisitNode; i++)
				__vecVisitBool[__vecVisitNode[i]] = false;
		}
		return 1.0 * numCoverd * __numV / numHyperEdge;
	}
	
};

using TRRcollection = RRcollection;
using PRRcollection = std::shared_ptr<TRRcollection>;
