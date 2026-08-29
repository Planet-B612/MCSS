#pragma once
#include "../dSFMT/dSFMT.h"
#include "CommonStruc.h"
#include <memory>
#include <queue>
#include "CommonFunc.h"
#include <cmath>
#include <set>
#include "Argument.h"
#include <algorithm>
#include <iomanip>
#include <climits>
#include <random>
#include <immintrin.h>
using namespace std;

class mRRcollection
{
	private:
	/// __numV: number of nodes in the graph.
	int __numV;
	/// __numE: number of edges in the graph.
	size_t __numE = 0;
	/// _num_mRRsets: number of RR sets.
	vector<int> __vecVisitBool;
	vint __vecTree;
	// vint __vecSeq;
	vint __vecNewTree;
	vint __vecVisitNode;
	float rand_div = 1.0;

	public:
	Graph PO;
	FRsets _FRsets;
	FRsets _FRsets_veri;
	mRRsets _mRRsets;
	mRRsets vec_mRR_layer;
	// #ifdef DEBUG
	vsint vec_hash_FR;
	vsint vec_hash_mRR;
	// #endif // !NDEBUG
	vint vecRoot_num;
	vint __vecSeq;
	ulint _num_mRRsets = 0, _num_mRRsets_veri = 0;
	int pre_root_num = 0;
	Argument *__arg;
	string model;
	string result;
	int floor_root_RR = 0;
	int ceil_root_RR = 0;
	int num_update_this_round=0;
	int num_add_root_this_round = 0;
	int num_delete_root_this_round = 0;
	double build_mRRset_time = 0.0;
	double revise_mRRset_time = 0.0;
	float __regen_threshold = 0.15; 
	int __root_num_bound = 250; // try to delete unnecessary mRR-sets when the number of roots in an mRR exceeds this value. sample:0, facebook: 25, dblp: 250
	vvint vv_polluted_nodes;
	vector<vint_aligned> vv_virtual_roots;
	vvint vv_next_mRRnode;
	std::random_device rd; // initialize random number generator

	double mRR_traversal_time = 0.0;

	explicit mRRcollection(Argument &arg)
	{
		__arg = &arg;
		__numV = arg.numV;
		_FRsets = FRsets(__numV);
		_FRsets_veri = FRsets(__numV);
		__regen_threshold = arg.regen_threshold;
		__root_num_bound = arg.root_num_bound;
		// #ifdef DEBUG
		// vec_hash_FR.resize(__numV);
		// #endif
		model = arg.model;
		__vecVisitBool = std::vector<int>(__numV, 0);
		__vecTree = std::vector<int>(__numV, -1);
		if(model != "IC")
		{
			__vecSeq = std::vector<int>(__numV, -1);
		}
		__vecNewTree = std::vector<int>(__numV, -1);
		__vecVisitNode = vint(__numV);
		result = arg.result_dir;
		PO.resize((__numV), vint_aligned());
		if (arg.real_time_pw == true)
		{
			generate_possible_world();
		}
		else // load previously generated PO
		{
			string pw_path;
			if(model=="IC")
			{
				pw_path = arg.pw_path + arg.dataset[arg.dataset_No] + "_pw_ic" + to_string(arg.times) + ".txt";
			}
			else
			{
				pw_path = arg.pw_path + arg.dataset[arg.dataset_No] + "_pw_lt" + to_string(arg.times) + ".txt";
				if(arg.times==2)
				{
					pw_path=arg.pw_path + arg.dataset[arg.dataset_No] + "_pw_lt6.txt";
					// cout<<"i = "<<arg.times<<" using pw file: "<< pw_path<<endl;
				}
			}
			cout << "used PO path: " + pw_path << endl;
			// pw_path+="_pw_ic.txt";
			ifstream load_pw;
			load_pw.open(pw_path);
			assert(!load_pw.fail());
			int i, nbr;
			while (!load_pw.eof())
			{
				load_pw >> i >> nbr;
				PO[i].push_back(nbr);
			}
			PO[i].pop_back(); // the last row is empty, due to the mechanism of eof, a duplicated nbr will be added. Thus, we need to pop_back here
		}
	}

	/// Genrerate a possible world, PO.
	void generate_possible_world()
	{
		for (int i = 0; i < (__numV); i++)
		{
			auto nbrs = (O_graph)[i];
			for (auto nbr : nbrs)
			{
				if ((dsfmt_gv_genrand_open_close() / rand_div) < Inv_inDeg[nbr])
				{
					PO[i].push_back(nbr);
				}
			}
		}
	}

	int realization(vint seeds)
	{
		int curr_Node = 0, numVisitNode = 0;
		int counter_real = 0; // local counter not used
		for (auto seed : seeds)
		{
			++counter_real;
			(__Activated)[seed] = 1;
			__vecVisitNode[numVisitNode++] = seed;
			for (const auto &rrid : _FRsets[seed])
			{
				if (rrid >= static_cast<int>(_num_mRRsets))
					continue;
				vv_polluted_nodes[rrid].push_back(seed);
			}
		}
		while (curr_Node < numVisitNode)
		{
			int expand = __vecVisitNode[curr_Node++];
			for (auto v : PO[expand])
			{
				if ((__Activated)[v])
					continue;
				__vecVisitNode[numVisitNode++] = v;
				++counter_real;
				(__Activated)[v] = 1;
				for (const auto &rrid : _FRsets[v])
				{
					if (rrid >= static_cast<int>(_num_mRRsets))
						continue;
					vv_polluted_nodes[rrid].push_back(v);
				}
			}
		}
		return counter_real;
	}

	int realization_fresh_vec(vint seeds)
	{
		int curr_Node = 0, numVisitNode = 0;
		int counter_real = 0; // local counter not used
		for (auto seed : seeds)
		{
			++counter_real;
			(__Activated)[seed] = 1;
			__vecVisitNode[numVisitNode++] = seed;
		}
		while (curr_Node < numVisitNode)
		{
			int expand = __vecVisitNode[curr_Node++];
			for (auto v : PO[expand])
			{
				if ((__Activated)[v])
					continue;
				__vecVisitNode[numVisitNode++] = v;
				++counter_real;
				(__Activated)[v] = 1;
			}
		}
		return counter_real;
	}

	#include "test_ic.h"

	/// Generate a set of n mRR sets
	void build_n_mRRsets_tree(const ulint numSamples, const int pre_theta)
	{
		floor_root_RR = 0;
		ceil_root_RR = 0;   // the number mRR-sets with root number root_num+1 in the previous revisable mRR-sets
		ulint prevSize = _num_mRRsets;								   // previous total number of mRR-sets
		vector<int> mRR_mark(prevSize, -1);								   // false indicates this mRR is not directly reused.
		ulint num_revise_RR = (prevSize > numSamples ? numSamples : prevSize); // the number of mRR-sets that will be revised (revisable mRR-sets)
		vint vec_rootnum_RRid, vec_rootnum_1_RRid;
		vec_rootnum_RRid.reserve(num_revise_RR);
		vec_rootnum_1_RRid.reserve(num_revise_RR);
        std::mt19937 gen(rd());
        std::binomial_distribution<int> dist(num_revise_RR, residual);
        ceil_root_RR=dist(gen);
		int num_polluted=0, root_num_1=root_num+1;

		auto single_start = std::chrono::high_resolution_clock::now();
		for(ulint i=pre_theta;i<num_revise_RR;i++) 
		{
			mRR_mark[i]=vv_polluted_nodes[i].size();
			if(mRR_mark[i]>0)
			{
				num_polluted++;
			}
		}
		size_t count = std::count_if(mRR_mark.begin(), mRR_mark.end(), [](int x) { return x > 0; });
		if(count<5 && _num_mRRsets>0)
		{
			cout<<"Too few polluted mRR-sets, please check the correctness of the input seeds."<<endl;
			exit(1);
		}		
		if(pre_theta>100 && (1.0*(num_polluted)/(numSamples-pre_theta)>__regen_threshold))  // brute regen is needed
		{
			refresh_FRmRRsets(pre_theta);
			num_revise_RR=0;
			prevSize=pre_theta;
		}
		else // sensible regen can be uses
		{
			for(ulint i=pre_theta;i<num_revise_RR;i++)  // build the basic information of previous mRR-sets, and update these mRR-sets
			{
				if(i%2000==0)
				{
					cout<<"updating mRR-sets: "<<i<<"/"<<num_revise_RR<<endl;
				}
				if (mRR_mark[i] > 0)
				{
					if(model=="IC")
					{
						if(ablation_update)
						{
							naive_mRR_update(i, vv_polluted_nodes[i]);
						}
						else
						{
							num_update++;
							mRR_update(i, vv_polluted_nodes[i]);
						}
						// naive_mRR_update(i, vv_polluted_nodes[i]);
					}
					else
					{
						mRR_update_lt(i, vv_polluted_nodes[i]);
					}
					vv_polluted_nodes[i].clear();
					num_update_this_round++;
				}
				// The root info should be recorded after the mRR-sets are updated.
				if (vecRoot_num[i] == root_num)
				{
					vec_rootnum_RRid.push_back(i);
				}
				else if (vecRoot_num[i] == root_num_1)
				{
					vec_rootnum_1_RRid.push_back(i);
				}
			}
		}
		floor_root_RR = num_revise_RR - ceil_root_RR;
		for (auto mRRid : vec_rootnum_RRid) // directly reuse updated previous mRR-sets with root number root_num, if there is any such mRR-sets
		{
			if (floor_root_RR > 0) // if still need floor_root_RR
			{
				mRR_mark[mRRid] = INT_MAX;  // re-usable
				floor_root_RR--;
			}
		}
		for (auto mRRid : vec_rootnum_1_RRid) // directly reuse updated previous mRR-sets with root number root_num+1
		{
			if (ceil_root_RR > 0)
			{
				mRR_mark[mRRid] = INT_MAX;  // re-usable
				ceil_root_RR--;
			}
		}
		int root_diff=0;
		for(ulint i=pre_theta;i<num_revise_RR;i++)
		{
			if(i%2000==0)
			{
				cout<<"adding roots to mRR-sets: "<<i<<"/"<<num_revise_RR<<endl;
			}
			if (mRR_mark[i] <__numV) // for mRR-sets that have not been directly reused
			{
				if (floor_root_RR > 0) // derive floor-rooted mRR first
				{
					root_diff = vecRoot_num[i] - root_num;
					if (root_diff > 0)
					{
						delete_root(i, root_diff);
						num_delete_root_this_round++;
					}
					else
					{
						if(model=="IC")
						{
							// auto single_start = std::chrono::high_resolution_clock::now();
							// add_root(i, -root_diff);
							if(ablation_add_root)
							{
								// naive_add_root(i, -root_diff);
								more_naive_add_root(i, -root_diff);
							}
							else
							{
								num_add_root++;
								add_root(i, -root_diff);
							}
							// naive_add_root(i, -root_diff);
							// more_naive_add_root(i, -root_diff);
						}
						else
						{
							add_root_lt(i, -root_diff);
						}
						num_add_root_this_round++;
					}
					floor_root_RR--;
				}
				else if (ceil_root_RR > 0)
				{
					root_diff = vecRoot_num[i] - root_num - 1;
					if (root_diff > 0)
					{
						delete_root(i, root_diff);
						num_delete_root_this_round++;
					}
					else
					{
						if(model=="IC")
						{
							if(ablation_add_root)
							{
								// naive_add_root(i, -root_diff);
								more_naive_add_root(i, -root_diff);
							}
							else
							{
								add_root(i, -root_diff);
							}
							// add_root(i, -root_diff);
							// naive_add_root(i, -root_diff);
							// more_naive_add_root(i, -root_diff);
						}
						else
						{
							add_root_lt(i, -root_diff);
						}
						num_add_root_this_round++;
					}
					ceil_root_RR--;
				}
			}
		}
		auto single_end = std::chrono::high_resolution_clock::now();
		revise_mRRset_time += std::chrono::duration<double>(single_end - single_start).count();
		
		if (prevSize < numSamples)
		{
			_num_mRRsets = numSamples;
			vv_virtual_roots.resize(numSamples);
			_mRRsets.resize(numSamples);
		// #ifdef DEBUG
		// 	vec_hash_mRR.resize(numSamples);
		// #endif
			vv_polluted_nodes.resize(numSamples);
			if(model=="IC")
			{
				vec_mRR_layer.resize(numSamples);
			}
			// else
			// {
			// 	vv_next_mRRnode.resize(numSamples);
			// }
			vecRoot_num.resize(numSamples);
		}

		// single_start = std::chrono::high_resolution_clock::now();
		for (auto i = prevSize; i < numSamples; i++)  // if the number of previous mRR-sets is not enough, new mRR-sets will be generated
		{
			build_one_mRRset_tree(i, root_num, residual);
		}
		// single_end = std::chrono::high_resolution_clock::now();
		// build_mRRset_time += std::chrono::duration<double>(single_end - single_start).count();
	}

	int build_one_mRRset_tree(int mRRid, int root_num, double residual)
	// Each adj_list in in the form of adjacency list, so that the first node of each entry automatically constitutes the original __vecVisitNode
	{
		int root;
		root_num += (dsfmt_gv_genrand_open_close() <= residual);
		vecRoot_num[mRRid] = root_num;
		mRRset &mRR = _mRRsets[mRRid];
		// vint &vec_next_mRRnode=vv_next_mRRnode[mRRid];
		// #ifdef DEBUG
		// sint &mRR_hash = vec_hash_mRR[mRRid];
		// #endif
		mRR.resize(root_num);
		auto &vec_RR_layer = vec_mRR_layer[mRRid];		
		vec_RR_layer.resize(root_num);
		for (int i = 0; i < root_num; i++) // roots should be independent, and thus are selected in advance, while the diffusion from them is dependent
		{
			root = dsfmt_gv_genrand_uint32_range(__numV);
			// if(numVisitNode>=__numV_left)  break;  // In the last round, the number of users visited by previous roots may be the whole network
			while ((__Activated)[root] || __vecVisitBool[root])
			{
				root = dsfmt_gv_genrand_uint32_range(__numV);
			}
			__vecVisitBool[root] = 1; // only record the state of roots, but do not push the root into the queue, since we are not going to diffuse here.
			_FRsets[root].push_back(mRRid);
			mRR[i].push_back(root);
// #ifdef DEBUG
// 			mRR_hash.insert(root);
// 			vec_hash_FR[root].insert(mRRid);
// #endif
		}
		if(model=="IC")
		{
			for (int i = 0; i < root_num;i++)
			{
				auto &RR = mRR[i];
				int layer_start = 0, layer_end = 1, expand;
				while (layer_start < layer_end)
				{
					vec_RR_layer[i].push_back(layer_start);
					for (int j = layer_start; j < layer_end; j++)
					{
						expand = RR[j];
						auto prob=Inv_inDeg[expand];
						for (const auto &nbrId : (R_graph)[expand])
						{
							if(__vecVisitBool[nbrId] || (__Activated)[nbrId] || dsfmt_gv_genrand_open_close() > prob)
								continue;
							RR.push_back(nbrId);
							__vecVisitBool[nbrId] = 1; // mark the node as in the new mRR
							#ifndef PREFETCH // be careful to here, this line should be applied only when _FRsets is prefetched below.
							_FRsets[nbrId].push_back(mRRid);
							#endif
						}
						#ifdef PREFETCH  // usefull for 5% acceleration
						if(RR.size()-j>8)
						{
							_mm_prefetch(R_graph[j+8].data(), _MM_HINT_T0);
						}
						#endif
					}
					layer_start = layer_end;
					layer_end = RR.size();	 // update the end index of the next layer
				}
				#ifdef PREFETCH
				for(int j=1;j<layer_end;j++)
				{
					if(j+2<layer_end)
					{
						auto &frset=_FRsets[RR[j+2]];
						_mm_prefetch(frset.data()+frset.size(), _MM_HINT_T0);
					}
					_FRsets[RR[j]].push_back(mRRid);
				}
				#endif
			}
		}
		else // LT model
		{
			// vec_next_mRRnode.reserve(root_num);
			for (int i = 0; i < root_num; i++)
			{
				auto &RR = mRR[i];
				int node= RR[0];
				while(true)
				{
					auto &nbrs= (R_graph)[node];
					ulint nbrs_size = nbrs.size();
					if (nbrs_size == 0)
					{
						// vec_next_mRRnode.push_back(-1);
						break;
					}
					int nbrId = nbrs[dsfmt_gv_genrand_uint32_range(nbrs_size)];
					if ((__Activated)[nbrId])
					{
						// vec_next_mRRnode.push_back(-1);						
						break;
					}
					if (__vecVisitBool[nbrId])
					{
						// vec_next_mRRnode.push_back(nbrId);
						break;
					}
					__vecVisitBool[nbrId] = 1;
					#ifndef PREFETCH
					_FRsets[nbrId].push_back(mRRid);
					#endif
					RR.push_back(nbrId);
					node = nbrId;
					// #ifdef DEBUG
					// 	mRR_hash.insert(nbrId);
					// 	vec_hash_FR[nbrId].insert(mRRid);
					// #endif
				}
				#ifdef PREFETCH
				ulint RR_size=RR.size();
				for(ulint j=1;j<RR_size;j++)
				{
					if(j+2<RR_size)
					{
						auto &frset=_FRsets[RR[j+2]];
						_mm_prefetch(frset.data()+frset.size(), _MM_HINT_T0);
					}
					_FRsets[RR[j]].push_back(mRRid);
				}
				#endif
			}
		}
		for (const auto &RR : mRR)
		{
			ulint RR_size = RR.size(), i=0;
			for(;i+16<=RR_size;i+=16)
			{
				__m512i indices = _mm512_load_epi32((const void*)(RR.data()+i));
				_mm512_i32scatter_epi32((void*)__vecVisitBool.data(), indices, zeros512, 4);
			}
			for (;i<RR_size;i++)
			{
				__vecVisitBool[RR[i]] = 0;
			}
		}
		// #ifdef DEBUG
		// if(vec_value_check(__vecVisitBool, false, 1, string(__func__) + "in the " + to_string(mRRid) +"-th mRR-set, __vecVisitBool includes TRUE values."))
		// {
		// 	print_single_mRRset(mRRid);
		// 	exit(1);
		// }
		// if(FR_sorted_check(__func__))
		// {
		// 	out_FR();
		// 	out_mRRset();
		// 	exit(1);
		// }
		// #endif
		return 0;
	}

	void build_n_mRRsets_fresh_vec(int theta)
	{
		for (int i = _num_mRRsets; i < theta; i++) // if the number of previous mRR-sets is not enough, new mRR-sets will be generated
		{
			build_one_mRRset_fresh_vec(i, root_num, residual);
		}
		_num_mRRsets = theta;
	}

	int build_one_mRRset_fresh_vec(int mRRid, int root_num, double residual)
	// Each adj_list in in the form of adjacency list, so that the first node of each entry automatically constitutes the original __vecVisitNode
	{
		int numVisitNode = 0, currNode = 0, root;
		root_num += (dsfmt_gv_genrand_open_close() <= residual);
		vint roots;
		for (int i = 0; i < root_num; i++) // roots should be independent, and thus are selected in advance, while the diffusion from them is dependent
		{
			root = dsfmt_gv_genrand_uint32_range(__numV);
			while ((__Activated)[root] || __vecVisitBool[root])
			{
				root = dsfmt_gv_genrand_uint32_range(__numV);
			}
			__vecVisitBool[root] = 1;
			_FRsets[root].push_back(mRRid);
			roots.push_back(root);
		}
		if(model=="IC")
		{
			for (int i = 0; i < root_num; i++)
			{
				root = roots[i]; // Take out a root
				__vecVisitNode[numVisitNode++] = root;
				while (currNode < numVisitNode)
				{
					const auto expand = __vecVisitNode[currNode];
					currNode++;
					for (auto &nbrId : (R_graph)[expand])
					{
						if (__vecVisitBool[nbrId] || (__Activated)[nbrId])
							continue;
						if (dsfmt_gv_genrand_open_close() > Inv_inDeg[expand])
							continue;
						__vecVisitNode[numVisitNode++] = nbrId;
						__vecVisitBool[nbrId] = 1;
						_FRsets[nbrId].push_back(mRRid);
					}
				}
			}
		}
		else
		{
			for (int i = 0; i < root_num; i++)
			{
				int node = roots[i];
				__vecVisitNode[numVisitNode++] = node;
				while(true)
				{
					auto &nbrs= (R_graph)[node];
					ulint nbrs_size = nbrs.size();
					if (nbrs_size == 0)
					{
						continue;
					}
					int nbrId = nbrs[dsfmt_gv_genrand_uint32_range(nbrs_size)];
					if (__vecVisitBool[nbrId] || (__Activated)[nbrId])
						continue;
					__vecVisitBool[nbrId] = 1;
					__vecVisitNode[numVisitNode++] = nbrId;
					_FRsets[nbrId].push_back(mRRid);
					node = nbrId;
				}
			}
		}
		for (int i = 0; i < numVisitNode; i++)
		{
			__vecVisitBool[__vecVisitNode[i]] = 0;
		}
		return 0;
	}

	int build_veri_mRRset_fresh(ulint theta, bool veri=false)
	{
		FRsets *FR;
		ulint *num;
		if(veri)
		{
			FR=&_FRsets_veri;
			num= &_num_mRRsets_veri;
		}
		else
		{
			FR=&_FRsets;
			num=&_num_mRRsets;
		}
		for (ulint mRRid = *num; mRRid < theta; mRRid++) 
		{
			int numVisitNode = 0, currNode = 0, root;
			root = dsfmt_gv_genrand_uint32_range(__numV);
			while ((__Activated)[root] || __vecVisitBool[root])
			{
				root = dsfmt_gv_genrand_uint32_range(__numV);
			}
			__vecVisitBool[root] = 1;
			(*FR)[root].push_back(mRRid);
			if(model=="IC")
			{
				__vecVisitNode[numVisitNode++] = root;
				while (currNode < numVisitNode)
				{
					const auto expand = __vecVisitNode[currNode];
					currNode++;
					for (auto &nbrId : (R_graph)[expand])
					{
						if (__vecVisitBool[nbrId] || (__Activated)[nbrId])
							continue;
						if (dsfmt_gv_genrand_open_close() > Inv_inDeg[expand])
							continue;
						__vecVisitNode[numVisitNode++] = nbrId;
						__vecVisitBool[nbrId] = 1;
						(*FR)[nbrId].push_back(mRRid);
					}
				}
			}
			else
			{
				__vecVisitNode[numVisitNode++] = root;
				while(true)
				{
					auto &nbrs= (R_graph)[root];
					ulint nbrs_size = nbrs.size();
					if (nbrs_size == 0)
					{
						continue;
					}
					int nbrId = nbrs[dsfmt_gv_genrand_uint32_range(nbrs_size)];
					if (__vecVisitBool[nbrId] || (__Activated)[nbrId])
						continue;
					__vecVisitBool[nbrId] = 1;
					__vecVisitNode[numVisitNode++] = nbrId;
					(*FR)[nbrId].push_back(mRRid);
					root = nbrId;
				}
			}
			for (int i = 0; i < numVisitNode; i++)
			{
				__vecVisitBool[__vecVisitNode[i]] = 0;
			}
		}
		*num = theta;		
		return 0;
	}

	void mRR_update(int mRRid, vint &del_nodes)
	{
		mRRset &mRR = _mRRsets[mRRid];
		mRRset &mRR_layer = vec_mRR_layer[mRRid];
		#ifdef DEBUG
		mRRset mRR_original = mRR, mRR_layer_original = mRR_layer;
		vint_aligned v_roots_original= vv_virtual_roots[mRRid];
		#endif
		int mRR_size = static_cast<int>(mRR.size());
		vint_aligned &v_roots = vv_virtual_roots[mRRid];
		ulint v_roots_size = v_roots.size();
		vint roots, del_roots;
		roots.reserve(v_roots_size + mRR_size);
		del_roots.reserve(mRR_size);

		ulint i = 0;
		vint del_v_roots;  del_v_roots.reserve(v_roots_size);
		for(; i+16<=v_roots_size; i+=16)
		{
			__m512i idx = _mm512_load_epi32((const void*)(v_roots.data() + i));
			__m512i active = _mm512_i32gather_epi32(idx, __Activated.data(), 4);
			__mmask16 mask = _mm512_cmpeq_epi32_mask(active, ones512);
			if (mask == 0)
				continue;
			while (mask) 
			{
				int bit = __builtin_clz((unsigned)mask);
				del_v_roots.push_back(i + bit);
				mask &= mask - 1;
        	}
		}
		for(; i<v_roots_size; i++)
		{
			if (__Activated[v_roots[i]]) // is a del_node
			{
				del_v_roots.push_back(i);
			}			
		}
		for(int i = static_cast<int>(del_v_roots.size()) - 1; i > -1; i--)
		{
			v_roots.erase(v_roots.begin() + i);
		}
		i=0;
		v_roots_size = v_roots.size();
		for (; i + 16 <= v_roots_size; i += 16) 
		{
			__m512i idx = _mm512_load_epi32((const void*)(v_roots.data()+i));
			_mm512_i32scatter_epi32(__vecNewTree.data(), idx, _mm512_set1_epi32(mRR_size), 4);
    	}
		for(; i < v_roots_size; i++)
		{
			__vecNewTree[v_roots[i]] = mRR_size;
		}
		#ifdef DEBUG
		if (v_roots_check(mRRid, string(__func__) + " end"))
		{
			exit(1);
		}
		#endif

		int min_tree = __numV, affected_layer_idx = __numV, first_del_idx = __numV;
		ulint min_tree_RR_size=__numV;
		bool find_del = false;
		for (int i = 0; i < mRR_size; i++) // mark previous roots, and traverse the mRRset. Traversing from the end is not necessary, since we need to know whether a node us already in the mRR if regenerating.
		{
			auto &RR = mRR[i];
			min_tree_RR_size = (RR.size());
			const __m512i treeId512 = _mm512_set1_epi32(i);
			ulint j = 0;
			for (; j+16 <= min_tree_RR_size; j+=16)
			{
				__m512i idx = _mm512_load_epi32((const void*)(RR.data()+j));
				_mm512_i32scatter_epi32(__vecTree.data(), idx, treeId512, 4);
				__m512i act = _mm512_i32gather_epi32(idx, __Activated.data(), 4);
				__mmask16 mask = _mm512_cmpeq_epi32_mask(act, ones512);
        		if (mask != 0) 
				{
					first_del_idx = __builtin_clz((unsigned)mask)+j; // count from high bit.
					find_del=true;					
					min_tree = i;
					break;
				}
			}
			for (; j < min_tree_RR_size; j++)
			{
				__vecTree[RR[j]] = i;
				if (__Activated[RR[j]])
				{
					first_del_idx = j;
					min_tree = i;
					find_del = true;
					break;
				}
			}
			if (find_del)
			{
				break;
			}
		}
		auto &min_tree_layer = mRR_layer[min_tree];
		affected_layer_idx = upper_bound(min_tree_layer.begin(), min_tree_layer.end(), first_del_idx) - min_tree_layer.begin() - 1;
		vint_aligned &min_tree_RR = mRR[min_tree];
		int affected_layer_beg = min_tree_layer[affected_layer_idx];
		int affected_next_layer_beg=0, min_tree_layer_size = static_cast<int>(min_tree_layer.size());
		if (affected_layer_idx == min_tree_layer_size - 1)
		{
			affected_next_layer_beg = min_tree_RR_size;
		}
		else
		{
			affected_next_layer_beg = min_tree_layer[affected_layer_idx + 1];
		}
		i=affected_next_layer_beg;
		auto i_aligned = min((((ulint)i + 15) / 16) * 16, min_tree_RR_size);
		for (; i < i_aligned; i++)
		{
			__vecTree[min_tree_RR[i]] = __numV;
		}
		for (; i+16 <= min_tree_RR_size; i+=16)
		{
			__m512i idx = _mm512_load_epi32((const void*)(min_tree_RR.data()+i));
			_mm512_i32scatter_epi32(__vecTree.data(), idx, numV512, 4);
		}
		for (; i < min_tree_RR_size; i++)
		{
			__vecTree[min_tree_RR[i]] = __numV;
		}
		mRRset mRR_copy(mRR_size - min_tree);
#ifdef DEBUG
		if (affected_next_layer_beg > min_tree_RR_size)
		{
			cout << "Error: affected_next_layer_beg > min_tree_RR_size in mRR_update, mRRid=" << mRRid << ", affected_next_layer_beg=" << affected_next_layer_beg << ", min_tree_RR_size=" << min_tree_RR_size << endl;
			exit(1);
		}
#endif
		mRR_copy[0].reserve(min_tree_RR_size - affected_next_layer_beg);
		mRR_copy[0] = vint_aligned(std::make_move_iterator(min_tree_RR.begin() + affected_next_layer_beg), std::make_move_iterator(min_tree_RR.end()));
		
		min_tree_RR.resize(affected_next_layer_beg);
		min_tree_layer.resize(affected_layer_idx); // do not contain the beggin of the affected layer.

		for (int i = min_tree + 1; i < mRR_size; i++) // traverse trees after the min_tree
		{
			auto &RR = mRR[i];
			if (!__Activated[RR[0]])
			{
				roots.push_back(RR[0]);
				__vecNewTree[RR[0]] = mRR_size; // should not be added into some tree during the new exploration
			}
			ulint j=0;
			const __m512i treeId512 = _mm512_set1_epi32(i);
			for(; j+16 <= RR.size(); j+=16)
			{
				__m512i idx = _mm512_load_epi32((const void*)(RR.data()+j));
				_mm512_i32scatter_epi32(__vecTree.data(), idx, treeId512, 4);
			}
			for (;j<RR.size();j++)
			{
				__vecTree[RR[j]] = i;
			}
			RR.swap(mRR_copy[i - min_tree]);
			mRR_layer[i].clear();
		}

		for (int i = static_cast<int>(v_roots_size - 1); i > -1; i--)
		{
			int root = v_roots[i];
			if (__vecTree[root] > min_tree) // delete realized v_roots
			{
				roots.push_back(root);
				v_roots.erase(v_roots.begin() + i);
			}
		}
		for (int i = affected_next_layer_beg - 1; i >= affected_layer_beg; i--)
		{
			int node = min_tree_RR[i];
			if (__Activated[node]) // only consider the activated nodes in this layer now
			{
				__vecTree[node] = -1;						// necessary to make __vecTree all -1
				min_tree_RR.erase(min_tree_RR.begin() + i); // remove del_nodes
				affected_next_layer_beg--;
			}
			else
			{
				__vecTree[node] = min_tree;
				__vecNewTree[node] = min_tree;
			}
		}
		int layer_start = affected_layer_beg, layer_end = affected_next_layer_beg, expand;
		// min_tree_layer.pop_back();
		while (layer_start < layer_end)
		{
			min_tree_layer.push_back(layer_start);
			for (int j = layer_start; j < layer_end; j++)
			{
				expand = min_tree_RR[j];
				auto prob=Inv_inDeg[expand];
				for (const auto &nbrId : (R_graph)[expand])
				{
					if (__Activated[nbrId] || (__vecTree[nbrId] > -1 && __vecTree[nbrId] <= min_tree) || __vecNewTree[nbrId] > -1 || dsfmt_gv_genrand_open_close() > prob)
						continue;
					min_tree_RR.push_back(nbrId);
					__vecNewTree[nbrId] = min_tree; // mark the node as in the new mRR
					#ifdef PREFETCH
					if (__vecTree[nbrId] < 0)		// nbrId was not in this mRR previously
					{
						auto &frset = _FRsets[nbrId];
						auto it = lower_bound(frset.begin(), frset.end(), mRRid);
						frset.insert(it, mRRid);
					}
					#endif
				}
				#ifdef PREFETCH  // usefull for 5% acceleration
				if(min_tree_RR.size()-j>8)
				{
					_mm_prefetch(R_graph[j+8].data(), _MM_HINT_T0);
				}
				#endif
			}
			layer_start = layer_end; // update the start index of the next layer
			layer_end = min_tree_RR.size();
		}
		#ifdef PREFETCH
		for(int j=layer_start;j<layer_end;j++)
		{
			if(j+2<layer_end)
			{
				auto &frset=_FRsets[min_tree_RR[j+2]];
				_mm_prefetch(frset.data()+frset.size(), _MM_HINT_T0);
			}
			if (__vecTree[min_tree_RR[j]] < 0)
			{
				auto &frset = _FRsets[min_tree_RR[j]];
				auto it = lower_bound(frset.begin(), frset.end(), mRRid);
				frset.insert(it, mRRid);
			}
		}
		#endif
		if (min_tree_RR.size() < 1) // make sure it is not empty, // if the last RR root is a del_node
		{
			mRR.resize(min_tree); // remove the last RR
			mRR_layer.resize(min_tree);
			mRR_size = min_tree;
		}
		else
		{
			mRR_size = min_tree + 1;
		}
		ulint num_new_roots = roots.size();
		mRR.resize(mRR_size + num_new_roots);
		mRR_layer.resize(mRR_size + num_new_roots);
		for (ulint i = 0; i < num_new_roots; i++)
		{
			int mRR_size_plus_i= mRR_size + i;
			auto &RR = mRR[mRR_size_plus_i];
			auto &RR_layer = mRR_layer[mRR_size_plus_i];
			RR.push_back(roots[i]);
			int layer_start = 0, layer_end = 1, node;
			while (layer_start < layer_end)
			{
				RR_layer.push_back(layer_start);
				for (int j = layer_start; j < layer_end; j++)
				{
					node = RR[j];
					double prob=Inv_inDeg[node];
					for (const auto &nbrId : (R_graph)[node])
					{
						if (__Activated[nbrId] || (__vecTree[nbrId] > -1 && __vecTree[nbrId] <= min_tree) || __vecNewTree[nbrId] > -1 || dsfmt_gv_genrand_open_close() > prob)
							continue;
						RR.push_back(nbrId);
						__vecNewTree[nbrId] = mRR_size_plus_i; // cannot use SIMD here, because this state may be immediately used.
						#ifdef PREFETCH
						if (__vecTree[nbrId] < 0) // nbrId was not in this mRR previously
						{
							auto &frset = _FRsets[nbrId];
							auto it = lower_bound(frset.begin(), frset.end(), mRRid);
							frset.insert(it, mRRid);
						}
						#endif
					}
					#ifdef PREFETCH  // usefull for 5% acceleration
					if(RR.size()-j>8)
					{
						_mm_prefetch(R_graph[j+8].data(), _MM_HINT_T0);
					}
					#endif
				}
				layer_start = layer_end; // update the start index of the next layer
				layer_end = RR.size();
			}
			#ifdef PREFETCH  // prefetch here may be unuseful, since the position to be inserted is unknown.
			for(int j=1;j<layer_end;j++)
			{
				if(j+2<layer_end)
				{
					auto &frset=_FRsets[RR[j+2]];
					_mm_prefetch(frset.data()+frset.size(), _MM_HINT_T0);
					_mm_prefetch(frset.data(), _MM_HINT_T0);
				}
				if (__vecTree[RR[j]] < 0)
				{
					auto &frset = _FRsets[RR[j]];
					auto it = lower_bound(frset.begin(), frset.end(), mRRid);
					frset.insert(it, mRRid);
				}
			}
			#endif
		}
		for(ulint i=0;i<mRR_copy.size();i++)
		{
			auto &RR = mRR_copy[i];
			ulint j=0, RR_size=RR.size();
			for (; j + 16 <= RR_size; j += 16)
        	{
				__m512i idx = _mm512_load_si512(reinterpret_cast<const void*>(RR.data() + j));
				__m512i vals = _mm512_i32gather_epi32(idx, __vecNewTree.data(), 4);
				__mmask16 mask = _mm512_cmplt_epi32_mask(vals, zeros512);
				_mm512_i32scatter_epi32(__vecTree.data(), idx, minus_one512, 4);  // set vecTree

				alignas(64) int nodes[16];
            	_mm512_store_si512(reinterpret_cast<void*>(nodes), idx);

				while (mask)
				{
					int node = nodes[__builtin_ctz(mask)];
					auto &frset = _FRsets[node];
					auto it = lower_bound(frset.begin(), frset.end(), mRRid);
						frset.erase(it);
					mask &= mask - 1;
				}
			}
			for (; j < RR_size; ++j)
			{
				int node = RR[j];
				if (__vecNewTree[node] < 0)
				{
				    auto &frset = _FRsets[node];
					auto it = lower_bound(frset.begin(), frset.end(), mRRid);
					frset.erase(it);
				}
				__vecTree[node] = -1;
			}
		}
		mRR_size = mRR.size();
		int safe_size = min_tree;
		if (min_tree >= mRR_size)
		{
			safe_size = mRR_size - 1;
		}
		for (int i = 0; i <= safe_size; i++)
		{
			auto &RR = mRR[i];
			ulint j=0, RR_size=RR.size();
			for (; j + 16 <= RR_size; j += 16)
        	{
				__m512i idx = _mm512_load_si512(reinterpret_cast<const void*>(RR.data() + j));
				_mm512_i32scatter_epi32(__vecTree.data(), idx, minus_one512, 4);
			}
			for (;j<RR_size;j++)
			{
				__vecTree[RR[j]] = -1; // reset the tree id
			}
		}

		for (int i = min_tree; i < mRR_size; i++)
		{
			auto &RR = mRR[i];
			ulint j=0, RR_size=RR.size();
			for(;j+16<=RR_size;j+=16)
			{
				__m512i idx = _mm512_load_si512(reinterpret_cast<const void*>(RR.data() + j));
				_mm512_i32scatter_epi32(__vecNewTree.data(), idx, minus_one512, 4);
			}
			for (;j<RR_size;j++)
			{
				__vecNewTree[RR[j]] = -1; // reset the tree id
			}
		}
		for (const auto &root : v_roots) // reset the v_roots
		{
			__vecNewTree[root] = -1; // reset the tree id
		}
		vecRoot_num[mRRid] = mRR_size + v_roots.size();
		#ifdef DEBUG
		if (synthetic_check(mRRid, string(__func__) + " end", 0, 0, 0, 0, 0, 0, 0) || identical_element_check(mRRid, string(__func__) + " end")) //previous: 1, 1, 0, 0, 0, 0, 0
		{
			out_mRRset(mRR_original);
			out_layer(mRR_layer_original);
			out_vec(v_roots_original);

			out_mRRset(mRR);
			out_layer(mRR_layer);
			out_vec(v_roots);
			exit(1);
		}
		#endif
	}

	void naive_mRR_update(int mRRid, vint &del_nodes)
	{
		mRRset &mRR = _mRRsets[mRRid];
		auto &vec_RR_layer = vec_mRR_layer[mRRid];
		auto &v_roots = vv_virtual_roots[mRRid];
		vint roots;
		for(const auto del_node: del_nodes)
		{
			__vecVisitBool[del_node] = 1; // is deleted nodes
		}
		for(const auto &RR:mRR)
		{
			int root=RR[0];
			if(__vecVisitBool[root]==0)
			{
				roots.push_back(root);
				__vecVisitBool[root]=1;
			}
			for(const auto &node:RR)
			{
				auto &frset = _FRsets[node];
				auto it = lower_bound(frset.begin(), frset.end(), mRRid);
				if (it != frset.end() && *it == mRRid)
				{
					frset.erase(it);
					// #ifdef DEBUG
					// vec_hash_FR[node].erase(mRRid);
					// mRR_hash.erase(node);
					// #endif // !NDEBUG
				}
				else
				{
					cout << __LINE__ << ", Error: mRRid is not in the FRset of " << node << endl;
				}
			}
			auto &frset = _FRsets[root];
			auto it = lower_bound(frset.begin(), frset.end(), mRRid);
			// #ifdef DEBUG
			// if (it != frset.end() && *it == mRRid)
			// {
			// 	cout << __LINE__ << ": Error: mRRid=" << mRRid << ", nbrId=" << root << " already in _FRsets." << endl;
			// 	exit(1);
			// }
			// #endif
			frset.insert(it, mRRid);
		}
		for(const auto &v_root:v_roots)
		{
			if(__vecVisitBool[v_root]==0)
			{
				roots.push_back(v_root);
				__vecVisitBool[v_root]=1;
				auto &frset = _FRsets[v_root];
				auto it = lower_bound(frset.begin(), frset.end(), mRRid);
				// #ifdef DEBUG
				// if (it != frset.end() && *it == mRRid)
				// {
				// 	cout << __LINE__ << ": Error: mRRid=" << mRRid << ", nbrId=" << v_root << " already in _FRsets." << endl;
				// 	exit(1);
				// }
				// #endif
				frset.insert(it, mRRid);
			}
		}
		mRR.clear();
		vec_RR_layer.clear();
		mRR.resize(roots.size());
		vec_RR_layer.resize(roots.size());
		if(model=="IC")
		{
			for (ulint i = 0; i < roots.size(); i++)
			{
				auto &RR = mRR[i];
				RR.push_back(roots[i]);
				ulint layer_start = 0, layer_end = 1;
				while (layer_start < layer_end)
				{
					vec_RR_layer[i].push_back(layer_start);
					for (ulint j = layer_start; j < layer_end; j++)
					{
						int node = RR[j];
						for (const auto &nbrId : (R_graph)[node])
						{
							if (__vecVisitBool[nbrId] || (__Activated)[nbrId])
								continue;
							if (dsfmt_gv_genrand_open_close() > Inv_inDeg[node])
								continue;
							RR.push_back(nbrId);
						// #ifdef DEBUG
						// 	mRR_hash.insert(nbrId);
						// 	vec_hash_FR[nbrId].insert(mRRid);
						// #endif
							__vecVisitBool[nbrId] = 1;
							auto &frset = _FRsets[nbrId];
							auto it = lower_bound(frset.begin(), frset.end(), mRRid);
							// #ifdef DEBUG
							// if (it != frset.end() && *it == mRRid)
							// {
							// 	cout << __LINE__ << ": Error: mRRid=" << mRRid << ", nbrId=" << root << " already in _FRsets." << endl;
							// 	exit(1);
							// }
							// #endif
							frset.insert(it, mRRid);
						}
					}
					layer_start = layer_end; // update the start index of the next layer
					layer_end = RR.size();	 // update the end index of the next layer
				}
			}
		}
		for (const auto &RR : mRR)
		{
			for(const auto &node:RR)
			{
				__vecVisitBool[node] = 0;
			}
		}
		for(const auto node:roots)
		{
			__vecVisitBool[node] = 0;
		}
		for(const auto del_node: del_nodes)
		{
			__vecVisitBool[del_node] = 0;
		}
		return;
	}

	void mRR_update_lt(int mRRid, vint &del_nodes)
	{
		mRRset &mRR = _mRRsets[mRRid];
		// vint &vec_next_mRRnode = vv_next_mRRnode[mRRid];
		// #ifdef DEBUG
		// mRRset pre_mRR = mRR;
		// sint &mRR_hash = vec_hash_mRR[mRRid];
		// if (synthetic_check(mRRid, string(__func__) + "beg", 1, 1, 0, 0, 1, 1))
		// {
		// 	exit(1);
		// }
		// #endif
		int mRR_size = static_cast<int>(mRR.size());
		auto &v_roots = vv_virtual_roots[mRRid];
		ulint v_roots_size = v_roots.size();
		vint roots, del_roots;
		// vvint vv_del_idx(mRR_size);
		// vint extendable_chains;
		// extendable_chains.reserve(mRR_size);
		roots.reserve(v_roots_size + mRR_size);
		del_roots.reserve(mRR_size);

		ulint i = 0;
		vint del_v_roots;  del_v_roots.reserve(v_roots_size);
		for(; i+16<=v_roots_size; i+=16)
		{
			__m512i idx = _mm512_load_epi32((const void*)(v_roots.data()+i));
			__m512i active = _mm512_i32gather_epi32(idx, __Activated.data(), 4);
			__mmask16 mask = _mm512_cmpeq_epi32_mask(active, ones512);
			if (mask == 0)
				continue;
			while (mask) 
			{
				int bit = __builtin_clz((unsigned)mask);
				del_v_roots.push_back(i + bit);
				mask &= mask - 1;
        	}
		}
		for(; i<v_roots_size; i++)
		{
			if (__Activated[v_roots[i]]) // is a del_node
			{
				del_v_roots.push_back(i);
			}			
		}
		for(int i = static_cast<int>(del_v_roots.size()) - 1; i > -1; i--)
		{
			v_roots.erase(v_roots.begin() + i);
		}
		i=0;
		v_roots_size = v_roots.size();
		const __m512i mRR_size512 = _mm512_set1_epi32(mRR_size);
		for (; i + 16 <= v_roots_size; i += 16) 
		{
			__m512i idx = _mm512_load_epi32((const void*)(v_roots.data()+i));
			_mm512_i32scatter_epi32(__vecNewTree.data(), idx, mRR_size512, 4);
    	}
		for(; i < v_roots_size; i++)
		{
			__vecNewTree[v_roots[i]] = mRR_size;
		}

		int min_tree = __numV, first_del_idx = __numV; 
		ulint min_tree_RR_size=__numV;
		bool find_del = false;
		for (int i = 0; i < mRR_size; i++)
		{
			auto &RR = mRR[i];
			min_tree_RR_size = (RR.size());
			ulint j = 0;
			for (; j+16 <= min_tree_RR_size; j+=16)
			{
				__m512i idx = _mm512_load_epi32((const void*)(RR.data()+j));
				_mm512_i32scatter_epi32(__vecTree.data(), idx, _mm512_set1_epi32(i), 4);
				__m512i act = _mm512_i32gather_epi32(idx, __Activated.data(), 4);
				__mmask16 mask = _mm512_cmpeq_epi32_mask(act, ones512);
        		if (mask != 0) 
				{
					first_del_idx = __builtin_clz((unsigned)mask)+j;
					find_del=true;					
					min_tree = i;
					break;
				}
			}
			for (; j < min_tree_RR_size; j++)
			{
				__vecTree[RR[j]] = i;
				if (__Activated[RR[j]])
				{
					first_del_idx = j;
					min_tree = i;
					find_del = true;
					break;
				}
			}
			if (find_del)
			{
				break;
			}
		}

		auto &min_tree_RR=mRR[min_tree];
		ulint first_del_idx_1=first_del_idx+1;
		i=first_del_idx_1;
		auto i_aligned = min((((ulint)i + 15) / 16) * 16, min_tree_RR_size);
		for (; i < i_aligned; i++)
		{
			__vecTree[min_tree_RR[i]] = __numV;
		}
		for (; i+16 <= min_tree_RR_size; i+=16)
		{
			__m512i idx = _mm512_load_epi32((const void*)(min_tree_RR.data()+i));
			_mm512_i32scatter_epi32(__vecTree.data(), idx, numV512, 4);
		}
		for (; i < min_tree_RR_size; i++)
		{
			__vecTree[min_tree_RR[i]] = __numV;
		// #ifdef DEBUG
		// 			mRR_hash.erase(min_tree_RR[i]);
		// #endif
		}
		mRRset mRR_copy(mRR_size - min_tree);  // mRR_copy should contain the truncated part of the min_tree_RR
		if(first_del_idx_1< min_tree_RR_size)
		{
			mRR_copy[0] = vint_aligned(std::make_move_iterator(min_tree_RR.begin() + first_del_idx_1), std::make_move_iterator(min_tree_RR.end()));
		}
		min_tree_RR.resize(first_del_idx);
		for (int i = min_tree + 1; i < mRR_size; i++)
		{
			auto &RR = mRR[i];
			if (__Activated[RR[0]]==0)
			{
				roots.push_back(RR[0]);
				__vecNewTree[RR[0]] = mRR_size; // the root should not be added into some tree during the new exploration
			}
			ulint j=0;
			const __m512i treeId512 = _mm512_set1_epi32(i);
			for(; j+16 <= RR.size(); j+=16)
			{
				__m512i idx = _mm512_load_epi32((const void*)(RR.data()+j));
				_mm512_i32scatter_epi32(__vecTree.data(), idx, treeId512, 4);
			}
			for (;j<RR.size();j++)
			{
				__vecTree[RR[j]] = i;
			}
			RR.swap(mRR_copy[i - min_tree]);
		}
		for (int i = static_cast<int>(v_roots_size - 1); i > -1; i--)
		{
			int root = v_roots[i];
			if (__vecTree[root] > min_tree) // delete realized v_roots
			{
				roots.push_back(root);
				v_roots.erase(v_roots.begin() + i);
			}
		}
		if(mRR[min_tree].size() < 1) 
		{
			mRR.resize(min_tree);
			mRR_size = min_tree; 
		}
		else
		{
			mRR.resize(min_tree + 1);
			mRR_size = min_tree + 1;
		}
		ulint num_new_roots = roots.size();
		mRR.resize(mRR_size + num_new_roots);
		for(ulint i = 0; i < num_new_roots; i++)
		{
			auto &RR = mRR[mRR.size() - num_new_roots + i];
			RR.push_back(roots[i]);
		// #ifdef DEBUG
		// 	if(mRR_hash.find(roots[i]) == mRR_hash.end())
		// 		mRR_hash.insert(roots[i]);
		// #endif // !NDEBUG
			int node= RR[0], mRR_size_plus_i = mRR_size + i;
			while(true)
			{
				auto &nbrs= (R_graph)[node];
				ulint nbrs_size = nbrs.size();
				if (nbrs_size == 0)
				{
					break;
				}
				int nbrId = nbrs[dsfmt_gv_genrand_uint32_range(nbrs_size)];
				if (__Activated[nbrId] || (__vecTree[nbrId] > -1 && __vecTree[nbrId] <= min_tree) || __vecNewTree[nbrId] > -1)
					break;
				__vecNewTree[nbrId] = mRR_size_plus_i;
				// #ifdef DEBUG
				// 	mRR_hash.insert(nbrId);
				// #endif
				if (__vecTree[nbrId] < 0) // nbrId was not in this mRR previously
				{
					auto &frset = _FRsets[nbrId];
					auto it = lower_bound(frset.begin(), frset.end(), mRRid);
					if (it != frset.end() && *it == mRRid)
					{
						cout <<__LINE__<< "Error: mRRid=" << mRRid << ", nbrId=" << nbrId << " already in _FRsets." << endl;
						exit(1);
					}
					frset.insert(it, mRRid);
				}
				RR.push_back(nbrId);
				node = nbrId;
			}
		}
		for(ulint i=0;i<mRR_copy.size();i++)
		{
			auto &RR = mRR_copy[i];
			ulint j=0, RR_size=RR.size();
			for (; j + 16 <= RR_size; j += 16)
        	{
				__m512i idx = _mm512_load_si512(reinterpret_cast<const void*>(RR.data() + j));
				__m512i vals = _mm512_i32gather_epi32(idx, __vecNewTree.data(), 4);
				__mmask16 mask = _mm512_cmplt_epi32_mask(vals, zeros512);
				_mm512_i32scatter_epi32(__vecTree.data(), idx, minus_one512, 4);  // set vecTree

				alignas(64) int nodes[16];
            	_mm512_store_si512(reinterpret_cast<void*>(nodes), idx);

				while (mask)
				{
					int node = nodes[__builtin_ctz(mask)];
					auto &frset = _FRsets[node];
					auto it = lower_bound(frset.begin(), frset.end(), mRRid);
					// if (it != frset.end() && *it == mRRid)
					// {
						frset.erase(it);
		// #ifdef DEBUG
		// 						vec_hash_FR[node].erase(mRRid);
		// #endif // !NDEBUG
					// }
					mask &= mask - 1;
				}
			}
			for (; j < RR_size; ++j)
			{
				int node = RR[j];
				if (__vecNewTree[node] < 0)
				{
				    auto &frset = _FRsets[node];
					auto it = lower_bound(frset.begin(), frset.end(), mRRid);
					frset.erase(it);
				}
				__vecTree[node] = -1;
			}
		}
		// #ifdef DEBUG
		// 	if (v_roots_check(mRRid, string(__func__) + " end"))
		// 	{
		// 		exit(1);
		// 	}
		// #endif
		mRR_size = mRR.size();
		int safe_size = min_tree;
		if (min_tree >= mRR_size)
		{
			safe_size = mRR_size - 1;
		}
		for (int i = 0; i <= safe_size; i++)
		{
			auto &RR = mRR[i];
			ulint j=0, RR_size=RR.size();
			for (; j + 16 <= RR_size; j += 16)
        	{
				__m512i idx = _mm512_load_si512(reinterpret_cast<const void*>(RR.data() + j));
				_mm512_i32scatter_epi32(__vecTree.data(), idx, minus_one512, 4);
			}
			for (;j<RR_size;j++)
			{
				__vecTree[RR[j]] = -1; // reset the tree id
			}
		}

		for (int i = min_tree; i < mRR_size; i++)
		{
			auto &RR = mRR[i];
			ulint j=0, RR_size=RR.size();
			for (; j + 16 <= RR_size; j += 16)
        	{
				__m512i idx = _mm512_load_si512(reinterpret_cast<const void*>(RR.data() + j));
				_mm512_i32scatter_epi32(__vecNewTree.data(), idx, minus_one512, 4);
			}
			for (;j<RR_size;j++)
			{
				__vecNewTree[RR[j]] = -1; // reset the tree id
			}
		}
		for (const auto &root : v_roots) // reset the v_roots
		{
			__vecNewTree[root] = -1; // reset the tree id
		}
		vecRoot_num[mRRid] = mRR_size + v_roots.size();
		// #ifdef DEBUG
		// if (synthetic_check(mRRid, string(__func__) + " end", 1, 1, 0, 0, 1, 1))
		// {
		// 	exit(1);
		// }
		// #endif
	}

	void add_root(int mRRid, int num)
	{
		vecRoot_num[mRRid] += num;
		mRRset &mRR = _mRRsets[mRRid];
		auto &vec_RR_layer = vec_mRR_layer[mRRid];
		#ifdef DEBUG
		mRRset mRR_original = mRR, mRR_layer_original = vec_RR_layer;
		vint_aligned v_roots_original= vv_virtual_roots[mRRid];
		#endif
		vint_aligned roots, new_roots;
		roots.reserve(root_num);
		new_roots.reserve(num);
		mRR.reserve(root_num + num);
		auto &v_roots = vv_virtual_roots[mRRid];
		for (const auto &RR : mRR) // mark previous roots, and traverse the mRRset
		{
			int root = RR[0];
			roots.push_back(root);
			__vecTree[root] = 1024;
			ulint RR_size=RR.size(), j=0;
			for (; j+16 <= RR_size; j+=16)
			{
				__m512i idx = _mm512_load_epi32((const void*)(RR.data()+j));
				_mm512_i32scatter_epi32(__vecVisitBool.data(), idx, ones512, 4);
			}
			for (;j<RR_size;j++)
			{
				__vecVisitBool[RR[j]] = 1;
			}
		}
		for (auto root : v_roots)
		{
			__vecTree[root] = 1024; // a root
		}
		#ifdef DEBUG
		if (v_roots_check(mRRid, string(__func__) + " end"))
		{
			exit(1);
		}
		#endif
		for (int j = 0; j < num; j++) // generate new roots
		{
			int root = dsfmt_gv_genrand_uint32_range(__numV);
			while (__Activated[root] || __vecTree[root] > 0)
			{
				root = dsfmt_gv_genrand_uint32_range(__numV);
			}
			__vecTree[root] = 1024;
			new_roots.push_back(root);
		}
		for (const auto &root : new_roots) // traverse the mRRset
		{
			if (__vecVisitBool[root])
			{
				v_roots.push_back(root);
				continue;
			}
			else
			{
				auto mRR_size_1 = mRR.size() + 1;
				mRR.resize(mRR_size_1);
				vec_RR_layer.resize(mRR_size_1);
				auto &RR = mRR[mRR_size_1 - 1];
				RR.push_back(root);
			// #ifdef DEBUG
			// 	mRR_hash.insert(root);
			// #endif
				__vecVisitBool[root] = 1;
				auto &frset = _FRsets[root];
				auto it = lower_bound(frset.begin(), frset.end(), mRRid);
			// #ifdef DEBUG
			// 	if (it != frset.end() && *it == mRRid)
			// 	{
			// 		cout << __LINE__ << ": Error: mRRid=" << mRRid << ", nbrId=" << root << " already in _FRsets." << endl;
			// 		exit(1);
			// 	}
			// #endif
				frset.insert(it, mRRid);
			// #ifdef DEBUG
			// 	mRR_hash.insert(root);
			// 	vec_hash_FR[root].insert(mRRid);
			// #endif // !NDEBUG

				int layer_start = 0, layer_end = 1;
				while (layer_start < layer_end)
				{
					vec_RR_layer[mRR_size_1 - 1].push_back(layer_start);
					for (int j = layer_start; j < layer_end; j++)
					{
						int node = RR[j];
						auto prob=Inv_inDeg[node];
						for (const auto &nbrId : (R_graph)[node])
						{
							if (__vecVisitBool[nbrId] || (__Activated)[nbrId] || dsfmt_gv_genrand_open_close() > prob)
								continue;
							RR.push_back(nbrId);
			// #ifdef DEBUG
			// 				mRR_hash.insert(root);
			// #endif
							__vecVisitBool[nbrId] = 1;
							auto &frset = _FRsets[nbrId];
							auto it = lower_bound(frset.begin(), frset.end(), mRRid);
			// #ifdef DEBUG
			// 				if (it != frset.end() && *it == mRRid)
			// 				{
			// 					cout << __LINE__ << ", Error: mRRid " << mRRid << " already in the FRset of node " << node << endl;
			// 					exit(1);
			// 				}
			// #endif
							frset.insert(it, mRRid);
			// #ifdef DEBUG
			// 				vec_hash_FR[nbrId].insert(mRRid);
			// #endif
						}
					}
					layer_start = layer_end; // update the start index of the next layer
					layer_end = RR.size();	 // update the end index of the next layer
				}
			}
		}
		// #ifdef DEBUG
		// if (v_roots_check(mRRid, string(__func__) + " end"))
		// {
		// 	exit(1);
		// }
		// mRRset a = {{}};
		// if (FR_reverse_check_hash(mRRid, a, "add_root"))
		// {
		// 	cout << __LINE__ << ", Error: FR_reverse_check_hash" << endl;
		// 	exit(1);
		// }
		// #endif
		for (const auto &RR : mRR)
		{
			ulint j=0, RR_size=RR.size();
			for (; j + 16 <= RR_size; j += 16)
			{
				__m512i idx = _mm512_load_si512(reinterpret_cast<const void*>(RR.data() + j));
				_mm512_i32scatter_epi32(__vecVisitBool.data(), idx, zeros512, 4);
			}
			for (; j < RR_size; ++j)
			{
				__vecVisitBool[RR[j]] = 0;
			}
		}
		ulint all_roots_size = roots.size()+v_roots.size()+new_roots.size(), k=0;
		roots.insert(roots.end(),
         std::make_move_iterator(v_roots.begin()),
         std::make_move_iterator(v_roots.end()));
		roots.insert(roots.end(),
         std::make_move_iterator(new_roots.begin()),
         std::make_move_iterator(new_roots.end()));
		for(;k+16<=all_roots_size;k+=16)
		{
			__m512i idx = _mm512_load_si512(reinterpret_cast<const void*>(roots.data() + k));
			_mm512_i32scatter_epi32(__vecTree.data(), idx, minus_one512, 4);
		}
		for (; k < all_roots_size; ++k)
		{
			__vecTree[roots[k]] = -1;
		}
		#ifdef DEBUG
		if (synthetic_check(mRRid, string(__func__) + " end", 0, 0, 0, 0, 0, 0, 0) || identical_element_check(mRRid, string(__func__) + " end")) //previous: 1, 1, 0, 0, 0, 0, 0
		{
			out_mRRset(mRR_original);
			out_layer(mRR_layer_original);
			out_vec(v_roots_original);

			out_mRRset(mRR);
			out_layer(vec_RR_layer);
			out_vec(v_roots);
			exit(1);
		}
		#endif
		return;
	}

	void more_naive_add_root(int mRRid, int num)
	{
		vecRoot_num[mRRid] += num;
		mRRset &mRR = _mRRsets[mRRid];
		auto &vec_RR_layer = vec_mRR_layer[mRRid];
		vint roots;
		for (const auto &RR : mRR) // mark previous roots, and traverse the mRRset
		{
			int root = RR[0];
			roots.push_back(root);
			__vecVisitBool[root] = 1;
			for(const auto &node:RR)
			{
				auto &frset = _FRsets[node];
				auto it = lower_bound(frset.begin(), frset.end(), mRRid);
				if (it != frset.end() && *it == mRRid)
				{
					frset.erase(it);
					// #ifdef DEBUG
					// vec_hash_FR[node].erase(mRRid);
					// mRR_hash.erase(node);
					// #endif // !NDEBUG
				}
				else
				{
					cout << __LINE__ << ", Error: mRRid is not in the FRset of " << node << endl;
				}
			}
			auto &frset = _FRsets[root];
			auto it = lower_bound(frset.begin(), frset.end(), mRRid);
			// #ifdef DEBUG
			// if (it != frset.end() && *it == mRRid)
			// {
			// 	cout << __LINE__ << ": Error: mRRid=" << mRRid << ", nbrId=" << root << " already in _FRsets." << endl;
			// 	exit(1);
			// }
			// #endif
			frset.insert(it, mRRid);
		}
		for (int j = 0; j < num; j++) // generate new roots
		{
			int root = dsfmt_gv_genrand_uint32_range(__numV);
			while (__Activated[root] || __vecVisitBool[root] == 1)
			{
				root = dsfmt_gv_genrand_uint32_range(__numV);
			}
			roots.push_back(root);
			__vecVisitBool[root] = 1;
			auto &frset = _FRsets[root];
			auto it = lower_bound(frset.begin(), frset.end(), mRRid);
			// #ifdef DEBUG
			// if (it != frset.end() && *it == mRRid)
			// {
			// 	cout << __LINE__ << ": Error: mRRid=" << mRRid << ", nbrId=" << root << " already in _FRsets." << endl;
			// 	exit(1);
			// }
			// #endif
			frset.insert(it, mRRid);
		}
		mRR.clear();
		vec_RR_layer.clear();
		mRR.resize(roots.size());
		vec_RR_layer.resize(roots.size());
		if(model=="IC")
		{
			vec_RR_layer.resize(roots.size());
			for (ulint i = 0; i < roots.size(); i++)
			{
				auto &RR = mRR[i];
				RR.push_back(roots[i]);
				int layer_start = 0, layer_end = 1;
				while (layer_start < layer_end)
				{
					vec_RR_layer[i].push_back(layer_start);
					for (int j = layer_start; j < layer_end; j++)
					{
						int node = RR[j];
						for (const auto &nbrId : (R_graph)[node])
						{
							if (__vecVisitBool[nbrId] || (__Activated)[nbrId])
								continue;
							if (dsfmt_gv_genrand_open_close() > Inv_inDeg[node])
								continue;
							RR.push_back(nbrId);
						// #ifdef DEBUG
						// 	mRR_hash.insert(nbrId);
						// 	vec_hash_FR[nbrId].insert(mRRid);
						// #endif
							__vecVisitBool[nbrId] = 1;
							auto &frset = _FRsets[nbrId];
							auto it = lower_bound(frset.begin(), frset.end(), mRRid);
							// #ifdef DEBUG
							// if (it != frset.end() && *it == mRRid)
							// {
							// 	cout << __LINE__ << ": Error: mRRid=" << mRRid << ", nbrId=" << root << " already in _FRsets." << endl;
							// 	exit(1);
							// }
							// #endif
							frset.insert(it, mRRid);
						}
					}
					layer_start = layer_end; // update the start index of the next layer
					layer_end = RR.size();	 // update the end index of the next layer
				}
			}
		}
		for (const auto &RR : mRR)
		{
			for(const auto &node:RR)
			{
			__vecVisitBool[node] = 0;
			}
		}
		for(const auto node:roots)
		{
			__vecVisitBool[node] = 0;
		}
		return;

	}

	void naive_add_root(int mRRid, int num)
	{
		vecRoot_num[mRRid] += num;
		mRRset &mRR = _mRRsets[mRRid];
		auto &vec_RR_layer = vec_mRR_layer[mRRid];
// #ifdef DEBUG
// 		sint &mRR_hash = vec_hash_mRR[mRRid];
// #endif // !NDEBUG
		vint roots, new_roots;
		roots.reserve(root_num);
		new_roots.reserve(num);
		mRR.reserve(root_num + num);
		for (const auto &RR : mRR) // mark previous roots, and traverse the mRRset
		{
			int root = RR[0];
			roots.push_back(root);
			__vecTree[root] = 1024;
			for (const auto &node : RR)
			{
				__vecVisitBool[node] = 1;
			}
		}
		bool virtual_root = false;
		for (int j = 0; j < num; j++) // generate new roots
		{
			int root = dsfmt_gv_genrand_uint32_range(__numV);
			while (__Activated[root] || __vecTree[root] > 0)
			{
				root = dsfmt_gv_genrand_uint32_range(__numV);
			}
			if (__vecVisitBool[root])
			{
				virtual_root = true;
			}
			__vecTree[root] = 1024;
			new_roots.push_back(root);
			auto &frset = _FRsets[root];
			auto it = lower_bound(frset.begin(), frset.end(), mRRid);
			// #ifdef DEBUG
			// if (it != frset.end() && *it == mRRid)
			// {
			// 	cout << __LINE__ << ": Error: mRRid=" << mRRid << ", nbrId=" << root << " already in _FRsets." << endl;
			// 	exit(1);
			// }
			// #endif
			frset.insert(it, mRRid);
		}
		roots.insert(roots.end(), new_roots.begin(), new_roots.end());
		ulint roots_size = roots.size();
		if(virtual_root) // generate from fresh
		{
			mRRset mRR_copy; mRR_copy.swap(mRR);
			mRR.clear();
			vec_RR_layer.clear();
			mRR.resize(roots_size);
			vec_RR_layer.resize(roots_size);
			for(ulint i = 0; i < roots_size; i++)
			{
				auto &RR = mRR[i];
				auto &layer_info=vec_RR_layer[i];
				RR.push_back(roots[i]);
// #ifdef DEBUG
// 				mRR_hash.insert(roots[i]);
// 				vec_hash_FR[roots[i]].insert(mRRid);
// #endif

				int layer_start = 0, layer_end = 1;
				while (layer_start < layer_end)
				{
					layer_info.push_back(layer_start);
					for (int j = layer_start; j < layer_end; j++)
					{
						int node = RR[j];
						for (const auto &nbrId : (R_graph)[node])
						{
							if (__vecTree[nbrId] > -1 || (__Activated)[nbrId])
								continue;
							if (dsfmt_gv_genrand_open_close() > Inv_inDeg[node])
								continue;
							RR.push_back(nbrId);
						// #ifdef DEBUG
						// 	mRR_hash.insert(nbrId);
						// 	vec_hash_FR[nbrId].insert(mRRid);
						// #endif
							if(__vecVisitBool[nbrId]==0)
							{
								auto &frset = _FRsets[nbrId];
								auto it = lower_bound(frset.begin(), frset.end(), mRRid);
								// if (it != frset.end() && *it == mRRid)
								// {
								// 	cout << "Error: mRRid=" << mRRid << ", nbrId=" << nbrId << " already in _FRsets." << endl;
								// 	exit(1);
								// }
								frset.insert(it, mRRid);
							}
							__vecTree[nbrId]=1024;
						}
					}
					layer_start = layer_end; // update the start index of the next layer
					layer_end = RR.size();	 // update the end index of the next layer
				}
			}
			for(const auto &RR:mRR_copy)
			{
				for(const auto node:RR)
				{
					__vecVisitBool[node]=0;
					if(__vecTree[node]<0)
					{
						auto &frset = _FRsets[node];
						auto it = lower_bound(frset.begin(), frset.end(), mRRid);
						if (it != frset.end() && *it == mRRid)
						{
							frset.erase(it);
							// #ifdef DEBUG
							// vec_hash_FR[node].erase(mRRid);
							// mRR_hash.erase(node);
							// #endif // !NDEBUG
						}
						else
						{
							cout << __LINE__ << ", Error: mRRid is not in the FRset of " << node << endl;
						}
					}
				}
			}
			for(const auto &RR:mRR)
			{
				for(const auto node:RR)
				{
					__vecTree[node]=-1;
				}
			}
		}
		else // simply add
		{
			for(const auto &root:new_roots)
			{
				auto mRR_size_1 = mRR.size() + 1;
				mRR.resize(mRR_size_1);
				vec_RR_layer.resize(mRR_size_1);
				auto &RR = mRR[mRR_size_1 - 1];
				RR.push_back(root);
				__vecVisitBool[root] = 1;
// #ifdef DEBUG
// 				mRR_hash.insert(root);
// 				vec_hash_FR[root].insert(mRRid);
// #endif // !NDEBUG
				int layer_start = 0, layer_end = 1;
				while (layer_start < layer_end)
				{
					vec_RR_layer[mRR_size_1 - 1].push_back(layer_start);
					for (int j = layer_start; j < layer_end; j++)
					{
						int node = RR[j];
						for (const auto &nbrId : (R_graph)[node])
						{
							if (__vecVisitBool[nbrId] || (__Activated)[nbrId] || __vecTree[nbrId]>-1)
								continue;
							if (dsfmt_gv_genrand_open_close() > Inv_inDeg[node])
								continue;
							RR.push_back(nbrId);
							// #ifdef DEBUG
							// mRR_hash.insert(nbrId);
							// vec_hash_FR[nbrId].insert(mRRid);
							// #endif
							__vecVisitBool[nbrId] = 1;
							auto &frset = _FRsets[nbrId];
							auto it = lower_bound(frset.begin(), frset.end(), mRRid);
							// #ifdef DEBUG
							// if (it != frset.end() && *it == mRRid)
							// {
							// 	cout << __LINE__ << ", Error: mRRid " << mRRid << " already in the FRset of node " << node << endl;
							// 	exit(1);
							// }
							// #endif
							frset.insert(it, mRRid);
						}
					}
					layer_start = layer_end; // update the start index of the next layer
					layer_end = RR.size();	 // update the end index of the next layer
				}
			}
			for (const auto &RR : mRR)
			{
				for (auto &node : RR)
				{
					__vecVisitBool[node] = 0;
				}
			}
			for(const auto root:roots)
			{
				__vecTree[root]=-1;
			}
		}
// #ifdef DEBUG
// 		mRRset a = {{}};
// 		if (FR_reverse_check_hash(mRRid, a, "add_root"))
// 		{
// 			cout << __LINE__ << ", Error: FR_reverse_check_hash" << endl;
// 			exit(1);
// 		}
// #endif
// #ifdef DEBUG
// 		if (synthetic_check(mRRid, string(__func__) + " end", 0, 0, 0, 0, 0, 1))
// 		{
// 			cout<<"num of new roots is "<<num<<endl;
// 			cout<<"the original mRR is "<<endl;
// 			for(const auto &RR:mRR_debug)
// 			{
// 				for(const auto node:RR)				
// 				{
// 					cout<<node<<" ";
// 				}
// 				cout<<endl;
// 			}	
// 			cout<<"the new mRR is "<<endl;
// 			for(const auto &RR:mRR)
// 			{
// 				for(const auto node:RR)				
// 				{
// 					cout<<node<<" ";
// 					auto &frset = _FRsets[node];
// 					auto it = lower_bound(frset.begin(), frset.end(), mRRid);
// 					if (it == frset.end())
// 					{
// 						cout << __LINE__ << ", Error: mRRid " << mRRid << " is not in the FRset of node " << node << endl;
// 						cout<<"the FRset of node "<<node<<" is: ";
// 						for(const auto i:frset)
// 						{
// 							cout<<i<<" ";
// 						}
// 						cout<<endl;
// 						cout<<"the hash FRset of node "<<node<<" is: ";
// 						for(const auto i:vec_hash_FR[node])
// 						{
// 							cout<<i<<" ";
// 						}
// 						cout<<endl;
// 						cout<<"the mRRhash is ";
// 						for(const auto i:mRR_hash)
// 						{
// 							cout<<i<<" ";
// 						}
// 						cout<<endl;
// 					}
// 				}
// 				cout<<endl;			
// 			}
// 			cout<<"virtual_root = "<<virtual_root<<endl;
// 			exit(1);
// 		}
// #endif
		return;
	}

	void add_root_lt(int mRRid, int num)
	{
		vecRoot_num[mRRid] += num;
		mRRset &mRR = _mRRsets[mRRid];
		// #ifdef DEBUG
		// 	if (synthetic_check(mRRid, string(__func__) + "beg", 1, 1, 0, 0, 1, 1))
		// 	{
		// 		exit(1);
		// 	}
		// 	sint &mRR_hash = vec_hash_mRR[mRRid];
		// #endif // !NDEBUG
		vint roots, new_roots;
		roots.reserve(root_num);
		new_roots.reserve(num);
		mRR.reserve(root_num + num);
		auto &v_roots = vv_virtual_roots[mRRid];
		for (const auto &RR : mRR) // mark previous roots, and traverse the mRRset
		{
			int root = RR[0];
			roots.push_back(root);
			__vecTree[root] = 1024;
			ulint RR_size=RR.size(), j=0;
			for (; j+16 <= RR_size; j+=16)
			{
				__m512i idx = _mm512_load_epi32((const void*)(RR.data()+j));
				_mm512_i32scatter_epi32(__vecVisitBool.data(), idx, ones512, 4);
			}
			for (;j<RR_size;j++)
			{
				__vecVisitBool[RR[j]] = 1;
			}
		}
		// #ifdef DEBUG
		// 	if (v_roots_check(mRRid, string(__func__) + " end"))
		// 	{
		// 		exit(1);
		// 	}
		// #endif
		for (auto root : v_roots)
		{
			__vecTree[root] = 1024; // a root
		}
		for (int j = 0; j < num; j++) // generate new roots
		{
			int root = dsfmt_gv_genrand_uint32_range(__numV);
			while (__Activated[root] || __vecTree[root] > 0)
			{
				root = dsfmt_gv_genrand_uint32_range(__numV);
			}
			__vecTree[root] = 1024;
			new_roots.push_back(root);
		}
		for (const auto &root : new_roots)
		{
			if (__vecVisitBool[root])
			{
				v_roots.push_back(root);
				continue;
			}
			auto mRR_size = mRR.size();
			mRR.resize(mRR_size+1);
			auto &RR = mRR[mRR_size];
			RR.push_back(root);
			// #ifdef DEBUG
			// 	mRR_hash.insert(root);
			// #endif
			__vecVisitBool[root] = 1;
			auto &frset = _FRsets[root];
			auto it = lower_bound(frset.begin(), frset.end(), mRRid);
			if (it != frset.end() && *it == mRRid)
			{
				cout << __LINE__ << ": Error: mRRid=" << mRRid << ", nbrId=" << root << " already in _FRsets." << endl;
				exit(1);
			}
			frset.insert(it, mRRid);
			// #ifdef DEBUG
			// 	mRR_hash.insert(root);
			// 	vec_hash_FR[root].insert(mRRid);
			// #endif // !NDEBUG
			int node=root;
			while(true)
			{
				auto &nbrs= (R_graph)[node];
				ulint nbrs_size = nbrs.size();
				if (nbrs_size == 0)
				{
					break;
				}
				int nbrId = nbrs[dsfmt_gv_genrand_uint32_range(nbrs_size)];
				if (__vecVisitBool[nbrId] || (__Activated)[nbrId])
					break;
				auto &frset = _FRsets[nbrId];
				auto it = lower_bound(frset.begin(), frset.end(), mRRid);
				if (it != frset.end() && *it == mRRid)
				{
					cout << "Error: mRRid=" << mRRid << ", nbrId=" << nbrId << " already in _FRsets." << endl;
					exit(1);
				}
				frset.insert(it, mRRid);
				// #ifdef DEBUG
				// 	mRR_hash.insert(root);
				// 	vec_hash_FR[nbrId].insert(mRRid);
				// #endif // !NDEBUG
				__vecVisitBool[nbrId] = 1;
				RR.push_back(nbrId);
				node = nbrId;
			}
		}
		// #ifdef DEBUG
		// 	if (v_roots_check(mRRid, string(__func__) + " end"))
		// 	{
		// 		exit(1);
		// 	}
		// 	mRRset a = {{}};
		// 	if (FR_reverse_check_hash(mRRid, a, "add_root"))
		// 	{
		// 		cout << __LINE__ << ", Error: FR_reverse_check_hash" << endl;
		// 		exit(1);
		// 	}
		// #endif
		for (const auto &RR : mRR)
		{
			ulint j=0, RR_size=RR.size();
			for (; j + 16 <= RR_size; j += 16)
			{
				__m512i idx = _mm512_load_si512(reinterpret_cast<const void*>(RR.data() + j));
				_mm512_i32scatter_epi32(__vecVisitBool.data(), idx, zeros512, 4);
			}
			for (; j < RR_size; ++j)
			{
				__vecVisitBool[RR[j]] = 0;
			}
		}
		ulint all_roots_size = roots.size()+v_roots.size()+new_roots.size(), k=0;
		roots.insert(roots.end(),
         std::make_move_iterator(v_roots.begin()),
         std::make_move_iterator(v_roots.end()));
		roots.insert(roots.end(),
         std::make_move_iterator(new_roots.begin()),
         std::make_move_iterator(new_roots.end()));
		for(;k+16<=all_roots_size;k+=16)
		{
			__m512i idx = _mm512_load_si512(reinterpret_cast<const void*>(roots.data() + k));
			_mm512_i32scatter_epi32(__vecTree.data(), idx, minus_one512, 4);
		}
		for (; k < all_roots_size; ++k)
		{
			__vecTree[roots[k]] = -1;
		}
		// #ifdef DEBUG
		// 	if (synthetic_check(mRRid, string(__func__) + " end", 0, 1, 1, 0, 1, 1))
		// 	{
		// 		exit(1);
		// 	}
		// #endif
	}

	void delete_root(int mRRid, int num_del_roots)
	{
		vecRoot_num[mRRid] -= num_del_roots;
		mRRset &mRR = _mRRsets[mRRid];
		ulint mRR_size = mRR.size();
		vint_aligned &v_roots = vv_virtual_roots[mRRid], roots;
		ulint v_roots_size = v_roots.size();
		if (v_roots_size >= static_cast<ulint>(num_del_roots))
		{
			v_roots.resize(v_roots_size - static_cast<ulint>(num_del_roots));
			num_del_roots = 0;
		}
		else
		{
			v_roots.clear();
			num_del_roots -= v_roots_size;
		}
		for (int i = 0; i < num_del_roots; i++)
		{
			auto &RR = mRR[mRR_size - 1 - i];
			for (const auto &node : RR)
			{
				auto &frset = _FRsets[node];
				auto it = lower_bound(frset.begin(), frset.end(), mRRid);
				if (it != frset.end() && *it == mRRid)
				{
					frset.erase(it);
// #ifdef DEBUG
// 					vec_hash_mRR[mRRid].erase(node);
// 					vec_hash_FR[node].erase(mRRid);
// #endif // !NDEBUG
				}
				else
				{
					cout << __LINE__ << ", Error: " << node << " is not in " << mRRid << ", when deleting it." << endl;
				}
			}
		}
		mRR.resize(mRR_size - num_del_roots);
		if(model=="IC")
		{
			vec_mRR_layer[mRRid].resize(mRR_size - num_del_roots);
		}
	}

	/// Refresh the RRsets
	void refresh_RRsets()
	{
		for (auto &fr:_FRsets)
		{
			fr.clear();
		}
		_mRRsets.clear();
		if(model=="IC")
		{
			vec_mRR_layer.clear();
		}
		vv_virtual_roots.clear();
		vecRoot_num.clear();
		vv_polluted_nodes.clear();
		_num_mRRsets = 0; // important	
	}

	void refresh_FRmRRsets(int max_size)
	{
		for (auto i = 0; i < __numV; i++)
		{
			auto &frset = _FRsets[i];
			auto it = lower_bound(frset.begin(), frset.end(), max_size);
			auto k = it - frset.begin();
			frset.resize(k);
		}
		for (ulint i = max_size; i < _num_mRRsets; i++)
		{
			mRRset().swap(_mRRsets[i]);
		}
		if(model=="IC")
		{
			for (ulint i = max_size; i < _num_mRRsets; i++)
			{
				mRRset().swap(vec_mRR_layer[i]);
			}
 			vec_mRR_layer.resize(max_size);
		}
		_mRRsets.resize(max_size);
// #ifdef DEBUG
// 		vec_hash_mRR.resize(max_size);
// #endif // !NDEBUG
		for (ulint i = max_size; i < _num_mRRsets; i++)
		{
			vint_aligned().swap(vv_virtual_roots[i]);
		}
		vv_virtual_roots.resize(max_size);

		vecRoot_num.resize(max_size);
		vv_polluted_nodes.resize(max_size);
		_num_mRRsets = max_size; // important
	}

	void refresh_mRRFRsets()
	{
		for (int i = 0; i < __numV; i++)
		{
			FRset().swap(_FRsets[i]);
			FRset().swap(_FRsets_veri[i]);
		}
		for (auto &vec : vv_virtual_roots)
		{
			vint_aligned().swap(vec);
		}
		// no need to refresh mRRsets, since it is never recorded in ending rounds
		_num_mRRsets = 0;
		_num_mRRsets_veri=0;
	}

	/// Release memory
	void release_memory()
	{
		refresh_RRsets();
		std::vector<int>().swap(__vecVisitBool);
		vint().swap(__vecVisitNode);
		FRsets().swap(_FRsets);
		vector<int>().swap(__vecTree);
		vector<vint_aligned>().swap(vv_virtual_roots);
	}

	/// Set cascade model
	// void set_cascade_model(const string model)
	// {
	// 	model = model;
	// }

	void out_PO()
	{
		std::ofstream po("/data/fc/graphInfo/new/sample_po.txt");
		for (long unsigned int i = 0; i < PO.size(); i++)
		{
			if (PO[i].size() > 0)
				po << i << "; ";
			for (auto node : PO[i])
			{
				po << node << " ";
			}
			po << endl;
		}
		po.close();
	}

	void out_graph(Graph g)
	{
		std::ofstream out_g("/data/fc/graphInfo/new/test_graph.txt");
		for (long unsigned int k = 0; k < g.size(); k++)
		{
			out_g << k << ": ";
			for (auto &node : g[k])
			{
				out_g << node << ", ";
			}
			out_g << endl;
		}
	}

	void print_single_mRRset(ulint mRRid)
	{
		for(auto &RR:_mRRsets[mRRid])
		{
			for(auto node:RR)
			{
				std::cout << node << " ";
			}
			std::cout << endl;
		}
	}

	void out_mRRset(mRRset &mRR)
	{
		std::ofstream out_mRR(mRR_output_path, std::ios::app);
		// for (ulint k = 0; k < _mRRsets.size(); k++)
		// {
			out_mRR <<"The mRR-set: "<<endl<<"--";
			for (auto &RR : mRR)
			{
				for(auto node:RR)
				{
					out_mRR << node << " ";
				}
				out_mRR << endl<<"--";
			}
			// out_mRR << endl;
		// }
	}

	void out_layer(mRRset &mRR_layer)
	{
		std::ofstream out_layer(layer_output_path, std::ios::app);
		// for(ulint k=0;k<vec_mRR_layer.size();k++)
		// {
			out_layer <<"The mRR-layer: "<<endl<<"--";
			for (auto &layer : mRR_layer)
			{
				for(auto node:layer)
				{
					out_layer << node << " ";
				}
				out_layer << endl<<"--";
			}
			// out_layer << endl;
		// }
	}

	void out_FR()
	{
		std::ofstream out_FR(FR_output_path, std::ios::app);
		for(ulint k=0;k<_FRsets.size();k++)
		{
			out_FR <<"The " <<k << "-th FR-set: ";
			for (auto &node : _FRsets[k])
			{
				out_FR << node << " ";
			}
			out_FR << endl;
		}
	}

	template<typename T>
	void out_vec(T vec)
	{
		std::ofstream out_vec("../vec.txt", std::ios::app);
		for(auto &node:vec)
		{
			out_vec << node << " ";
		}
		out_vec << endl;
	}
};

using TmRRcollection = mRRcollection;
using PmRRcollection = std::shared_ptr<TmRRcollection>;