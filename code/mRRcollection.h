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
using namespace std;

class mRRcollection
{
	private:
	/// __numV: number of nodes in the graph.
	int __numV;
	/// __numE: number of edges in the graph.
	size_t __numE = 0;
	/// _num_mRRsets: number of RR sets.
	vector<bool> __vecVisitBool;
	vector<bool> __vecAffected;
	vector<tuple<int, int, int, bool>> __possible_pars;
	vint __vecTree;
	vint __vecSeq;
	vint __vecNewTree;
	vint __vecNewSeq;
	vvint vv_out_nbr_q;
	vector<int> __vecPreState;
	Nodelist __vecVisitNode;
	Nodelist __vecParentNode;
	vector<vector<int>> vec_adj_size;
	vector<vector<int>> vv_virtual_roots;
	float rand_div=1.0;

	public:
	vector<vector<int>> PO;
	FRsets _FRsets;
	mRRsets _mRRsets;
	mRRsets _mRRpars;
	Mchains _Mchns;
	int num_renew_mRR=0;
	size_t _num_mRRsets = 0;
	double decimal=1.0;
	int root_num = 1;
	double residual=0.0;
	int pre_root_num=0;
	int_vec_patchmap root_num_map;
	Argument *__arg;
	string _cascadeModel;
	vector<double> __Inv_inDeg;
	string result;
	float __q_ratio;
	int __linear_search_thr;
	int num_update=0;
	int num_add_root=0;
	int num_delete=0;
	vector<vector<int>> vv_roots;
	vint vec_round;


	double mRR_traversal_time = 0.0;

	explicit mRRcollection(Argument & arg)
	{
		__arg=&arg;
		__Inv_inDeg=arg.Inv_inDeg;
		__numV = arg.numV;
		_FRsets = FRsets(__numV);
		__vecVisitBool = std::vector<bool>(__numV, false);
		__vecAffected = std::vector<bool>(__numV, false);
		__possible_pars = std::vector<tuple<int, int, int, bool>>(__numV, make_tuple(INT_MAX, -1, -1, false));
		__vecTree = std::vector<int>(__numV, -1);
		__vecSeq = std::vector<int>(__numV, -1);
		__vecNewTree = std::vector<int>(__numV, -1);
		__vecNewSeq = std::vector<int>(__numV, -1);
		__vecPreState = std::vector<int>(__numV, -1);
		vv_out_nbr_q = vector<vector<int>>(__numV, vector<int>());
		#ifdef debug
		__vecVisitNode = Nodelist(5*__numV);
		__vecParentNode = Nodelist(5*__numV);
		#else
		__vecVisitNode = Nodelist(__numV);
		__vecParentNode = Nodelist(__numV);
		#endif
		_cascadeModel=arg.model;
		result=arg.result_dir;
		__linear_search_thr=arg.linear_search_thr;
		PO.resize((__numV), vector<int>());
		if(arg.real_time_pw==true)
		{
			generate_possible_world();
		}
		else  // load previously generated PO
		{
			string pw_path=arg.pw_path + arg.dataset[arg.dataset_No] + "_pw_ic" + to_string(arg.times) + ".txt";
			cout << "used PO path: " + pw_path <<endl;
			// pw_path+="_pw_ic.txt";
			ifstream load_pw;
			load_pw.open(pw_path);
			assert(!load_pw.fail());
			int i, nbr;
			while(!load_pw.eof())
			{
				load_pw>>i>>nbr;
				PO[i].push_back(nbr);
			}
			PO[i].pop_back();  // the last row is empty, due to the mechanism of eof, a duplicated nbr will be added. Thus, we need to pop_back here
		}
		// __q_ratio = arg.q_ratio;
	}

	/// Genrerate a possible world, PO.
	void generate_possible_world()
	{

		for(int i=0;i<(__numV);i++)
		{
			auto nbrs=(O_graph)[i];
			for(auto nbr:nbrs)
			{
				if((dsfmt_gv_genrand_open_close()/rand_div)<__Inv_inDeg[nbr])
				{
					PO[i].push_back(nbr);
				}
			}
		}		
	}

	int realization(Nodelist seeds)
	{
		int curr_Node=0, numVisitNode = 0; 
		int counter_real = 0;  // local counter not used
		for (auto seed : seeds)
		{		
		 	++counter_real;  // only one seed
		 	(__Activated)[seed] = true;
		 	__vecVisitNode[numVisitNode++]=seed;
		}
		while (curr_Node<numVisitNode)
		{
			int expand = __vecVisitNode[curr_Node++];
			for (auto v : PO[expand])
			{
				if ((__Activated)[v])continue;
				__vecVisitNode[numVisitNode++]=v;
				++counter_real;
				(__Activated)[v] = true;
			}
		}
		Nodelist temp(__vecVisitNode.begin(), __vecVisitNode.begin()+numVisitNode); 
		activated_nodes.insert((activated_nodes).end(),temp);  
		return counter_real;
	}

 #include "test_ic.h"

	/// Generate a set of n mRR sets
	void build_n_mRRsets_tree(const size_t numSamples, const int pre_theta)
	{
		root_num_map.clear();  // the map of previous revisable mRR-sets
		int floor_root_RR=0;
		int ceil_root_RR=0;  // the number mRR-sets with root number root_num+1 in the previous revisable mRR-sets
		const auto prevSize = _num_mRRsets;  // previous total number of mRR-sets
		vector<bool> mRR_mark(prevSize, false); // false indicates this mRR is not directly reused.
		decimal = 1.0 * (__numV_left) / (__eta_left);
		root_num = floor(decimal);
		residual = decimal - root_num;  // in (0,1)
		int num_revise_RR=(prevSize>numSamples?numSamples:prevSize);  // the number of mRR-sets that will be revised (revisable mRR-sets)
		vector<int> vecRoot_num(num_revise_RR,0);
		int curr_round;
		if(prevSize<numSamples)
		{
			vv_virtual_roots.resize(numSamples);
			vv_roots.resize(numSamples);
			vec_round.resize(numSamples);
			_mRRpars.resize(numSamples);
			_mRRsets.resize(numSamples);
		}
		for(int i=0;i<num_revise_RR;i++)  // build the basic information of previous mRR-sets, and update these mRR-sets
		{
			Nodelist polluted_nodes; polluted_nodes.reserve(100);
			curr_round=vec_round[i];
			for(auto j=curr_round;j<round_num;j++)  // append all previously polluted nodes, which have not take effect.  CHECK the for limits
			{
				polluted_nodes.insert(polluted_nodes.end(), (activated_nodes)[j].begin(), (activated_nodes)[j].end());
			}
			if(curr_round<round_num)
			{
				Nodelist p_nodes; p_nodes.reserve(polluted_nodes.size());
				for(auto p_node:polluted_nodes)
				{
					// if(_FRsets[p_node].find(i)!=_FRsets[p_node].end())
					auto &frset= _FRsets[p_node];
					auto it=lower_bound(frset.begin(), frset.end(), i);
					if(it != frset.end() && *it == i)
					{
						p_nodes.push_back(p_node);
					}
				}
				vec_round[i]=round_num;  // the mRR is also updated even there is no p_node for it
				if(p_nodes.size()>0)
				{
					mRR_update(i, p_nodes);
				}
			}
			// The root info should be recorded after the mRR-sets are updated.
			pre_root_num=vv_roots[i].size()+vv_virtual_roots[i].size();  // do not use mRR size, since there exist virtual vv_roots that do not add an adj_list but add a root to rnd_roots
			root_num_map[pre_root_num].emplace_back(i);
			vecRoot_num[i]=pre_root_num;
			if(dsfmt_gv_genrand_open_close() <= residual)  ceil_root_RR++;  // Get the number of ceil_root_num to be generated BTW
		}
		floor_root_RR=num_revise_RR-ceil_root_RR;
		for(auto mRRid:root_num_map[root_num])  // directly reuse updated previous mRR-sets with root number root_num, if there is any such mRR-sets
		{
			if(floor_root_RR>0)  // if still need floor_root_RR
			{
				mRR_mark[mRRid]=true;
				floor_root_RR--;
			}
		}
		for(auto mRRid:root_num_map[root_num+1])  // directly reuse updated previous mRR-sets with root number root_num+1
		{
			if(ceil_root_RR>0)
			{
				mRR_mark[mRRid]=true;
				ceil_root_RR--;
			}
		}
		int root_diff=0;
		for(int i=0;i<num_revise_RR;i++)
		{
			if(mRR_mark[i]==false)  // for mRR-sets that have not been directly reused
			{
				if(floor_root_RR>0)  // derive floor-rooted mRR first
				{
					root_diff=vecRoot_num[i]-root_num;
					if(root_diff>0)
					{
						delete_root(i,root_diff);
					}
					else
					{
						// output_info(i,true);
						add_root(i,-root_diff);
						// output_info(i,false);
					}
					floor_root_RR--;
				}
				else if(ceil_root_RR>0)
				{
					root_diff=vecRoot_num[i]-root_num-1;
					if(root_diff>0)
					{
						// output_info(i,true);
						delete_root(i,root_diff);
						// output_info(i,false);
					}
					else
					{
						// output_info(i,true);
						add_root(i,-root_diff);
						// output_info(i,false);
					}
					ceil_root_RR--;
				}
			}
		}
		for (auto i = prevSize; i < numSamples; i++)  // if the number of previous mRR-sets is not enough, new mRR-sets will be generated
		{
			build_one_mRRset_tree(i, root_num, residual);
		}
		if(prevSize<numSamples)
		{
			_num_mRRsets=numSamples;
		}
	}

	#include "test_ic.h"

	int build_one_mRRset_tree(int mRRid, int root_num, double residual)
	// Each adj_list in in the form of adjacency list, so that the first node of each entry automatically constitutes the original __vecVisitNode
	{
		int root;
		root_num += (dsfmt_gv_genrand_open_close() <= residual);
		mRRset &mRR=_mRRsets[mRRid];
		mRR.resize(root_num);
		// auto mRR_size = mRR.size();
		// vint &roots=vv_roots[mRRid];
		vec_round[mRRid]=round_num;
		for (int i = 0; i < root_num; i++) // roots should be independent, and thus are selected in advance, while the diffusion from them is dependent
		{
			root = dsfmt_gv_genrand_uint32_range(__numV);
			// if(numVisitNode>=__numV_left)  break;  // In the last round, the number of users visited by previous roots may be the whole network
			while ((__Activated)[root] || __vecVisitBool[root])
			{
				root = dsfmt_gv_genrand_uint32_range(__numV);
			}
			__vecVisitBool[root] = true; // only record the state of roots, but do not push the root into the queue, since we are not going to diffuse here.
			_FRsets[root].push_back(mRRid);
			mRR[i].push_back(root);
		}
		int idx=0;
		for (auto &RR:mRR)
		{
			int numVisitNode = 1, currNode = 0;
			while(currNode<numVisitNode)
			{
				const auto expand=RR[currNode++];
				for (const auto &nbrId : (R_graph)[expand])
				{
					if (__vecVisitBool[nbrId] || (__Activated)[nbrId])
						continue;
                    if (dsfmt_gv_genrand_open_close() > __Inv_inDeg[expand])
						continue;
					RR.push_back(nbrId);
					numVisitNode++;
					__vecVisitBool[nbrId] = true;
					_FRsets[nbrId].push_back(mRRid);
				}
			}
		}
		for(const auto &RR:mRR)
		{
			for (const auto &expand : RR)
			{
				__vecVisitBool[expand] = false;
			}
		}
		vec_value_check(__vecVisitBool, false, 1, string(__func__) + "beg="+to_string(0)+" __vecVisitBool includes TRUE values.");
		FR_sorted_check(__func__);
		return 0;
	}

	void build_n_mRRsets_fresh_vec(int theta)
	{
		decimal = 1.0 * (__numV_left) / (__eta_left);
		root_num = floor(decimal);
		residual = decimal - root_num;  // in (0,1)
		for (int i = _num_mRRsets; i < theta; i++)  // if the number of previous mRR-sets is not enough, new mRR-sets will be generated
		{
			build_one_mRRset_fresh_vec(i, root_num, residual);
		}
		_num_mRRsets=theta;
	}

	int build_one_mRRset_fresh_vec(int mRRid, int root_num, double residual)
	// Each adj_list in in the form of adjacency list, so that the first node of each entry automatically constitutes the original __vecVisitNode
	{
		int numVisitNode = 0, currNode = 0;
		int root;
		root_num += (dsfmt_gv_genrand_open_close() <= residual);
		vector<int> rnd_roots;
		for (int i = 0; i < root_num; i++) // roots should be independent, and thus are selected in advance, while the diffusion from them is dependent
		{
			root = dsfmt_gv_genrand_uint32_range(__numV);
			while ((__Activated)[root] || __vecVisitBool[root])
			{
				root = dsfmt_gv_genrand_uint32_range(__numV);
			}
			__vecVisitBool[root] = true;
			_FRsets[root].push_back(mRRid);
			rnd_roots.push_back(root);
		}
		for (int i = 0; i < root_num; i++) // skip the round element
		{
			root = rnd_roots[i];						   // Take out a root
			__vecVisitNode[numVisitNode++] = root; 
			while (currNode < numVisitNode)
			{
				const auto expand = __vecVisitNode[currNode];
				currNode++;
				
				for (auto &nbrId : (R_graph)[expand])
				{
					if (__vecVisitBool[nbrId] || (__Activated)[nbrId])
						continue;
					double randDouble;
					randDouble = dsfmt_gv_genrand_open_close();
					if (randDouble > __Inv_inDeg[expand])
						continue;
					__vecVisitNode[numVisitNode++] = nbrId;
					__vecVisitBool[nbrId] = true;
					_FRsets[nbrId].push_back(mRRid);
				}
			}
		}
		for (int i = 0; i < numVisitNode; i++)
		{
			__vecVisitBool[__vecVisitNode[i]] = false;
		}
		return 0;
	}


	void mRR_update(int mRRid, Nodelist &del_nodes)
	{
        vec_round[mRRid]=round_num;
		mRRset &mRR=_mRRsets[mRRid];
		ulint mRR_size=mRR.size(), del_nodes_size=del_nodes.size();
		vint &v_roots=vv_virtual_roots[mRRid], q;
		ulint v_roots_size=v_roots.size();
        vint roots; roots.reserve(v_roots_size+mRR_size);
		vvint del_nodes_into_RR(mRR_size);  // record the v_roots in each tree
		// child_end may not be necessary
        mRRset mRR_copy;
		int min_tree=0, first_del_node, first_del_idx;
        bool find_del=false;
		for (const auto &RR:mRR)  // mark previous roots, and traverse the mRRset. Traversing from the end is not necessary, since we need to know whether a node us already in the mRR if regenerating.
		{
            if(!find_del)
            {
                for(const auto &node : RR)
                {
                    __vecTree[node] = min_tree;
                    if(__Activated[node])
                    {
                        first_del_idx=&node-RR.begin();
                        find_del=true;
                        first_del_node=node;
                        break;
                    }
                }
                if(find_del)
                {
                    break;
                }
            }
			min_tree++;
		}
        mRR_copy.resize(mRR_size-min_tree);
        q.reserve(mRR[min_tree].size());
        for(ulint j=min_tree;j<mRR_size;j++)
        {
            auto &RR_j=mRR[j];
            if(!__Activated[RR_j[0]] && j>min_tree)  // only records roots after min_tree
            {
                roots.push_back(RR_j[0]);
            }
            q.reserve(q.capacity()+RR_j.size());
            RR_j.swap(mRR_copy[j-min_tree]);
        }
        mRR.resize(min_tree+1);
		for(ulint i=v_roots_size-1;i>-1;i--)  // I did not record the idx of v_roots here like before
		{
			int root=v_roots[i];
			if(__Activated[root])  // is a del_node
			{
				v_roots.erase(v_roots.begin()+i);
				continue;
			}
			if(__vecTree[root]>=min_tree)  // only realize v_roots that are affected by del_nodes
			{
				roots.push_back(root); // to facilitate the regeneration process
				v_roots.erase(v_roots.begin()+i);  // this v_root will be realized and thus should be removed from the v_roots
                __vecTree[root] = mRR_size; // mark realized v_roots
			}
		}
        auto &last_RR= mRR_copy[0];
        int pre_del_idx=first_del_idx;
        for(ulint i=0;i<first_del_idx;i++)
        {
            __vecSeq[last_RR[i]] = i;  // record the sequence of nodes in the RR
            if(__vecTree[last_RR[i]]==mRR_size)  // a v_root
            {
                first_del_idx=i;
                first_del_node=last_RR[i];
                break;
            }
        }
        if(pre_del_idx!=first_del_idx)
        {
            for(ulint i=first_del_idx;i<pre_del_idx;i++)
            {
                __vecTree[last_RR[i]] = -1;
            }
        }
        int par_node_seq=-1;
        bool unique_par=true;
        const auto &nbr_first_del_node=O_graph[first_del_node];
        for(const auto &node : nbr_first_del_node)
        {
            if(__vecTree[node]==min_tree && __vecSeq[node]<first_del_idx)
            {
                if(!unique_par)  // there are multiple parents, we can not tell which influenced the first del_node
                {
                    par_node_seq=-1;
                    first_del_idx=1;
                    break;
                }
                par_node_seq=__vecSeq[node];
                unique_par=false;
            }
        }        
        auto &new_RR=mRR[min_tree];
        new_RR.insert(new_RR.begin(), last_RR.begin(), last_RR.begin()+first_del_idx);  // check the range
        if(par_node_seq>-1)  // the first del_node has a unique parent
        {
            auto it=std::upper_bound(nbr_first_del_node.begin(), nbr_first_del_node.end(), first_del_node); // the first nbr larger than first_del_node
            for(ulint i=it-nbr_first_del_node.begin();i<nbr_first_del_node.size();i++)
            {
                const int node=nbr_first_del_node[i];
                for (const auto &nbrId : (R_graph)[node])
				{
					if (__vecVisitBool[nbrId] || (__Activated)[nbrId])
						continue;
					if (dsfmt_gv_genrand_open_close() > __Inv_inDeg[node])
						continue;
                    new_RR.push_back(nbrId);
					__vecTree[nbrId] = min_tree;
            }
        }
        for(ulint i=first_del_idx+1;i<new_RR.size();i++)
        {
            const int node=new_RR[i];
            for(const auto &nbrId : (R_graph)[node])
            {
                if(__vecTree[nbrId]>-1 || __Activated[nbrId])
                {
                    continue;
                }
                if(dsfmt_gv_genrand_open_close() > __Inv_inDeg[node])
                {
                    continue;
                }
                new_RR.push_back(nbrId);
                __vecTree[nbrId] = min_tree;
            }
        }
        ulint num_new_roots=roots.size();
        mRR.resize(min_tree+1+num_new_roots);
        for(ulint i=0;i<num_new_roots;i++)
        {
            auto &RR=mRR[min_tree+1+i];
            RR.push_back(roots[i]);
            __vecTree[roots[i]] = min_tree+1+i;
            int numVisitNode = 1, currNode = 0;
			while(currNode<numVisitNode)
			{
				const auto expand=RR[currNode++];
				for (const auto &nbrId : (R_graph)[expand])
				{
					if (__vecTree[nbrId]>-1 || (__Activated)[nbrId])
						continue;
                    if (dsfmt_gv_genrand_open_close() > __Inv_inDeg[expand])
						continue;
					RR.push_back(nbrId);
					numVisitNode++;
					__vecTree[nbrId] = min_tree+1+i;
				}
			}
        }

        

                
        for(ulint j=i;j<first_del_idx;j++)
        {
            __vecTree[last_RR[j]] = -1;
        }
		for(ulint i=0;i<mRR_size;i++)
		{
			const auto &del_nodes_i=del_nodes_into_RR[i];
			if(!del_nodes_i.empty())
			{
				const auto &RR=mRR[i];
				const auto &RRpar=mRRpar[i];
				const auto RR_size_i=RR.size();
				auto &child_start_i=vv_child_start[i]; //, &child_end_i=vv_child_end[i];
				auto &del_intvl_beg_i=del_interval_beg[i];
				auto &del_interval_end_i=del_interval_end[i];
				int k=0; // number of del_nodes that have added intervals
				for(const auto &node:del_nodes_i)
				{
					int par=RRpar[node];
					bool flag=false;
					while(par!=-1)
					{
						if(__Activated[par])  // meet another del_node
						{
							flag=true;
							break;
						}
						else
						{
							par=RRpar[par];
						}
					}
					if(flag)
					{
						continue;
					}
					q.push_back(node);
					int seq=__vecSeq[node];
					int head=seq, tail=seq, head_1=head+1, tail_1=seq+1, q_size=1;
					auto this_pos=del_intvl_beg_i.begin();
					auto intvl_end_pos=del_interval_end_i.begin();
					if(k==0)
					{
						del_intvl_beg_i.push_back(seq);	
						del_interval_end_i.push_back(seq);
					}
					else
					{
						int idx=ascend_sorted_insert(del_intvl_beg_i, seq, this_pos, k);
						// this_pos = std::upper_bound(this_pos, this_pos+k, seq);
						// del_interval_i.insert(this_pos, seq);
						del_interval_end_i.insert(intvl_end_pos+idx, seq);
					}
					while(true)
					{
						bool find_head=false;
						ulint j=tail_1;
						for(;j<RR_size_i;j++)
						{
							if(child_start_i[j]<0)
							{
								continue;
							}
							else
							{
								tail=j-1;
								tail_1=j;
								break;
							}
						}
						if(j==RR_size_i)
						{
							tail=j-1;
							tail_1=j;
						}
						for(ulint idx=head;idx<=tail;idx++)
						{
							if(child_start_i[idx]<0)
							{
								continue;
							}
							else
							{
								head=idx;
								find_head=true;
								break;
							}
						}					
						if(find_head)  // implies both true
						{
							head=child_start_i[head]; // the new head	
							q.insert(q.end(), RR.begin()+head, RR.begin()+tail+1);  // check the range
							// here, it should not be push_back directly, since there can be multiple del_nodes.
							// Instead, for the first del_node, push_back is okay. For subsequent nodes, the first interval should use binary seach to find the position. For subsequent intervals, we only need to do binary search from the last position to the position+k, where k is the number of del_nodes that have added intervals.
							int last_size=tail-head+1;
							if(k==0)
							{
								del_intvl_beg_i.push_back(head);
								del_interval_end_i.push_back(tail);
							}
							else
							{
								int idx=ascend_sorted_insert(del_intvl_beg_i, head, this_pos, k);
								del_interval_end_i.insert(del_interval_end_i.begin()+idx, tail);
							}
						}
						else
						{
							break;
						}
					}
					++k;
				}
			}
		}
		total_affected_num=q.size();
		// for(const auto &node : del_nodes)
		// {
		// 	int tree_id=__vecTree[node];
		// 	ulint now_num=vec_RR_size[tree_id]-__vecSeq[node]; // # nodes affected by this del_node
		// 	int &pre_num=node_num_Tree[tree_id].first;
		// 	if(now_num>pre_num)
		// 	{
		// 		total_affected_num+=(now_num-pre_num);
		// 		pre_num=now_num;
		// 		if(pre_num<1)  // have not been counted as affected
		// 		{
		// 			regen_num+=vec_RR_size[tree_id];
		// 		}
		// 	}
		// 	if(tree_id<min_del_tree)
		// 	{
		// 		min_del_tree=tree_id;
		// 	}
		// }
		if(total_affected_num>regen_num*__q_ratio)
		{
			// regenerate the mRRset
			mRRset mRR_copy(mRR_size-min_del_tree);
			for(int i=min_del_tree;i<mRR_size;i++)
			{
				mRR[i].swap(mRR_copy[i-min_del_tree]);
			}
			mRR.resize(min_del_tree); //check the size value!
			for(int i=min_del_tree;i<mRR_size;i++)
			{
				// 1. regen from the original root
				auto &roots=vv_roots[i];
				ulint roots_size=roots.size();
				ulint now_size=mRR.size();
				mRR.resize(now_size+roots_size);
				mRRpar.resize(now_size+roots_size);
				for(ulint j=0;j<roots_size;j++)				
				{
					vint &RR=mRR[now_size+j];
					RR.push_back(roots[j]);  // preserve the root
					vint &RRpar=mRRpar[now_size+j];
					RRpar.push_back(-1);
					int numVisitNode = 1, currNode = 0;
					while(currNode<numVisitNode)
					{
						const auto expand=RR[currNode++];
						for (const auto &nbrId : R_graph[expand])
						{
							if (__vecVisitBool[nbrId] || (__Activated)[nbrId] || (__vecTree[nbrId]>-1)&&(__vecTree[nbrId]<min_del_tree))
								continue;
							double randDouble;
							randDouble = dsfmt_gv_genrand_open_close();
							if (randDouble > __Inv_inDeg[expand])
								continue;
							RR.push_back(nbrId);
							RRpar.push_back(expand);
							__vecVisitBool[nbrId] = true;				
							numVisitNode++;
							if(__vecTree[nbrId]<0)  // nbrId was not in this mRR previously
							{
								auto &frset = _FRsets[nbrId];
								auto it=lower_bound(frset.begin(), frset.end(), mRRid);
								frset.insert(it, mRRid);
							}
						}
					}
				}
			}
			for(const auto &RR:mRR_copy)
			{
				for(const auto & node:RR)
				{
					__vecTree[node] = -1;  // reset the tree id
					if(__vecVisitBool[node]==false)  // have not been added into mRR
					{
						auto &frset = _FRsets[node];
						auto it=lower_bound(frset.begin(), frset.end(), mRRid);
						if (it != frset.end() && *it == mRRid) 
						{
							frset.erase(it); 
						}
					}
				}
			}
			for(const auto &roots:vv_roots)
			{
				for(const auto &root:roots)
				{
					__vecVisitBool[root] = false;
				}
			}
			ulint now_size=mRR.size();
			for(ulint i=min_del_tree;i<mRR_size;i++)
			{
				auto &RR=mRR[i];
				for(const auto &node:RR)
				{
					__vecVisitBool[node] = false;
				}
			}
		}
		else
		{
			// 2. update the mRRset
			// we do not delete the affected nodes immediately. Instead, we preserve them to facilitate the look up. After the nodes finish finding new parents, we carry out deletion and addition together.			
			vvint prefix_deleted(mRR_size);
			vvint insert_intvl_pos(mRR_size);
			vvint insert_intvl_size(mRR_size);
			vector<int> del_nodes_idx; del_nodes_idx.reserve(del_nodes_size);
			for(ulint i=0;i<total_affected_num;i++)
			{
				int node=q[i];
				if(__Activated[node])  //  should not be a del_node (included in __vecAffected)
				{
					del_nodes_idx.push_back(i);
					continue;
				}
				__vecAffected[node]=true;  // only mark addable nodes as affected
			}
			#ifndef NDEBUG
			if(del_nodes_size!=del_nodes_idx.size())
			{
				cout<<__func__<<": del_nodes index record error, del_nodes_size="<<del_nodes_size<<", del_nodes_idx.size()="<<del_nodes_idx.size()<<endl;
			}
			#endif
			for(ulint j=del_nodes_size-1;j>-1;j--)
			{
				q.erase(q.begin()+del_nodes_idx[j]);
			}
			ulint q_size=q.size();
			vvint vv_in_nbr;vv_in_nbr.resize(q_size);
			for(ulint i=0;i<q_size;i++)
			{
				int node=q[i];
				int tree_id=__vecTree[node], seq=__vecSeq[node], pre_par=mRRpar[tree_id][seq], par_seq=__vecSeq[pre_par];
				for(const auto &nbrId : O_graph[node])
				{
					int nbr_tree=__vecTree[nbrId], nbr_seq=__vecSeq[nbrId];
					if(nbr_tree<tree_id || (nbr_tree==tree_id && nbr_seq<par_seq))  // including the case where nbrId is never in the mRR. nbrId should not be earlier than node's previous parent, or an affected node. seq is compared with <, to allow the previous parent
					{
						continue;
					}
					if(get<0>(__possible_pars[node])<nbr_tree)
					{
						continue;
					}
					else if(get<0>(__possible_pars[node])==nbr_tree)
					{
						if(get<1>(__possible_pars[node])<nbr_seq)
						{
							continue;
						}		
					}
					if(nbrId==pre_par)  // add it directly, since it is the previous parent
					{
						vv_in_nbr[i].push_back(nbrId);
						continue;
					}
					// for other nbrs, we should check whether they can influence this node
					if (dsfmt_gv_genrand_open_close() > __Inv_inDeg[nbrId])
					{
						continue;
					}
					if(__vecAffected[nbrId] && __Activated[nbrId]==false)  // nbrId is an affected node, but not a del_node
					{
						vv_in_nbr[i].push_back(nbrId);
						vv_out_nbr_q[nbrId].push_back(node);
						continue;
					}
					__possible_pars[node]=make_tuple(nbr_tree, nbr_seq, nbrId, false);
				}
			}
			for(ulint i=mRR_size-1;i>=min_del_tree;i--)  // delete affected nodes
			{
				auto &RR=mRR[i];
				auto &RRpar=mRRpar[i];
				auto &del_interval_beg_i=del_interval_beg[i];
				auto &del_interval_end_i=del_interval_end[i];
				auto &prefix_deleted_i=prefix_deleted[i];
				ulint del_intvl_num=del_interval_beg_i.size();
				prefix_deleted_i.resize(del_intvl_num+1,0);
				for(ulint j=del_intvl_num-1;j>=0;j--)
				{
					auto del_beg=RR.begin()+del_interval_beg_i[j], del_end=RR.begin()+del_interval_end_i[j]+1;
					RR.erase(del_beg, del_end);
					RRpar.erase(del_beg, del_end);
					prefix_deleted_i[j+1]=del_end-del_beg+1;
				}
				for(ulint j=1;j<=del_intvl_num;j++)
				{
					prefix_deleted_i[j]+=prefix_deleted_i[j-1];  // prefix sum
				}
			}
			// try to add back nodes directly
			vint added; added.reserve(q_size);
			for(ulint i=0;i<q_size;i++)
			{
				int node=q[i], par_tree_id=get<0>(__possible_pars[node]), par_seq=get<1>(__possible_pars[node]);
				bool can_add=true;
				for(int uncertained_possible_par:vv_in_nbr[i])
				{
					int possible_tree_id=get<0>(__possible_pars[uncertained_possible_par]);
					if(__vecNewTree[possible_tree_id]<0) // have not been added back
					{
						can_add=false;
						break;  // break here, since the time of nbrs can not be ascertained
					}
					else
					{
						// can the below be simplified to the comparation between tuples??
						if(par_tree_id<possible_tree_id)  // current node's possible parent is earlier
						{
							continue;
						}
						else if(par_tree_id==possible_tree_id)
						{
							if(par_seq<get<1>(__possible_pars[uncertained_possible_par]))   // current node's possible parent is earlier
							{
								continue;
							}
							// the below else deals with the possible case that uncertained_possible_par and __possible_pars[node] want to be added to the same position
							else if(par_seq==get<1>(__possible_pars[uncertained_possible_par])) // node's current possible parent is earlier
							{
								auto &RRpar=mRRpar[possible_tree_id];
								int pre_par_possible_par=get<2>(__possible_pars[node]), uncertained_possible_par_par=get<2>(__possible_pars[uncertained_possible_par]);
								if(__vecSeq[RRpar[__vecSeq[pre_par_possible_par]]]<__vecSeq[RRpar[uncertained_possible_par_par]])
								{
									continue;
								}
							}
						}
						__possible_pars[node]=make_tuple(possible_tree_id, get<1>(__possible_pars[uncertained_possible_par]), uncertained_possible_par, false);  // update the possible parent
					}
				}
				if(can_add)
				{
					// add the node back
					auto &RR=mRR[par_tree_id];
					auto &RRpar=mRRpar[par_tree_id];
					auto &child_start_i=vv_child_start[par_tree_id];
					if(child_start_i[par_seq]>0)
					{
						int insert_pos_del=upper_bound(prefix_deleted[par_tree_id].begin(), prefix_deleted[par_tree_id].end(), child_start_i[par_seq])-prefix_deleted[par_tree_id].begin()-1;
						RR.insert(RR.begin()+child_start_i[par_seq], node);
						RRpar.insert(RRpar.begin()+child_start_i[par_seq], get<2>(__possible_pars[node]));
					}
					__vecNewTree[node]=par_tree_id;  // update the tree id
					added.push_back(i);
					get<3>(__possible_pars[node])=true;
					
				}
			}
			for(ulint j=added.size()-1;j>-1;j--)
			{
				q.erase(q.begin()+added[j]);
				vv_in_nbr.erase(vv_in_nbr.begin()+added[j]);
			}
			// make indexed min_heap which allows us to know the index of the node in the heap
			make_heap(q.begin(), q.end(), compare);
			int tree_id, RR_size, seq;
			for(const auto &node:q)
			{
				tree_id=get<0>(__possible_pars[node]);
				if(tree_id>__numV)
				{
					break;
				}
				const auto &RR=mRR[tree_id];
				const auto &RRpar=mRRpar[tree_id];
				const auto &child_start_i=vv_child_start[tree_id];
				RR_size=vec_RR_size[tree_id];
				seq=get<1>(__possible_pars[node]);
				if(child_start_i[seq]<0) // -1, no child, then, insert the node before the first node's child_start
				{
					ulint j=seq+1;
					for(;j<RR_size;j++)  // note for the case seq=RR_size-1
					{
						if(child_start_i[j]>-1)
						{
							break;
						}
					}						
					get<1>(__possible_pars[node])=j;  // put it here to consider the case that seq=RR_size-1
				}
				else  // >=0, has child, insert it into its children list
				{
					int node_par=RRpar[node];
					for(ulint j=child_start_i[seq];j<RR_size;j++)  // note for the case seq=RR_size-1
					{
						if(RR[j]>node)
						{
							get<1>(__possible_pars[node])=j;  // here, no exception will occur, even if seq=RR_size-1
							break;
						}
						else if(RRpar[j]!=node_par)
						{
							get<1>(__possible_pars[node])=j;
							break;
						}
					}
				}
				get<3>(__possible_pars[node])=true;  // mark that the node has found someone to add back
			}
		}
	}

	void add_root(int mRRid, int num)
	{
		mRRset &mRR=_mRRsets[mRRid];
		vint roots, new_roots; roots.reserve(root_num); new_roots.reserve(num); mRR.reserve(root_num+num);
        auto &v_roots=vv_virtual_roots[mRRid];
		for (const auto &RR:mRR)  // mark previous roots, and traverse the mRRset
		{
			int root = RR[0];
			roots.push_back(root);
			__vecTree[root] = 1024;
			for(const auto &node : RR)
			{
				__vecVisitBool[node] = true;
			}
		}
        for(auto root:v_roots)
		{
			__vecTree[root]=1024; // a root
		}
		for (int j = 0; j < num; j++)  // generate new roots
		{
			int root = dsfmt_gv_genrand_uint32_range(__numV);
			while (__Activated[root] || __vecTree[root]>0)
			{
				root = dsfmt_gv_genrand_uint32_range(__numV);
			}
			__vecTree[root] = 1024;
			new_roots.push_back(root);
		}
		#ifndef NDEBUG
		cout<<"True num : "<<count(__vecVisitBool.begin(), __vecVisitBool.end(), true)<<endl;
		#endif
		for(const auto &root : new_roots)  // traverse the mRRset
		{
			if(__vecVisitBool[root])
			{
				vv_virtual_roots[mRRid].push_back(root);
				continue;
			}
			else
			{
				auto mRR_size_1 = mRR.size()+1;
				mRR.resize(mRR_size_1);
				auto &RR=mRR[mRR_size_1-1];
				RR.push_back(root);

				__vecVisitBool[root] = true;
				auto &frset = _FRsets[root];
				auto it=lower_bound(frset.begin(), frset.end(), mRRid);
				frset.insert(it, mRRid);

				int numVisitNode = 1, currNode = 0;
				while(currNode<numVisitNode)
				{
					const auto expand=RR[currNode++];
					for (const auto &nbrId : R_graph[expand])
					{
						if (__vecVisitBool[nbrId] || (__Activated)[nbrId])
							continue;
						double randDouble;
						randDouble = dsfmt_gv_genrand_open_close();
						if (randDouble > __Inv_inDeg[expand])
							continue;
						RR.push_back(nbrId);
						__vecVisitBool[nbrId] = true;						
						numVisitNode++;
						auto &frset = _FRsets[nbrId];
						auto it=lower_bound(frset.begin(), frset.end(), mRRid);
						frset.insert(it, mRRid);
					}
				}
			}
		}
		for(const auto &RR:mRR)
		{
			for (const auto &expand : RR)
			{
				__vecVisitBool[expand] = false;
			}
		}
		for (const auto &root : roots)
		{
			__vecTree[root] = -1;
		}
		for (const auto &root : new_roots)
		{
			__vecTree[root] = -1;
		}
		vec_value_check(__vecVisitBool, false, 1, string(__func__) + "beg="+to_string(0)+" __vecVisitBool includes TRUE values.");
		vec_value_check(__vecTree, -1, 1, string(__func__) + "beg="+to_string(0)+" __vecTree includes non -1 values.");
		FR_sorted_check(__func__);
		return;
	}

	void delete_root(int mRRid, int i)
	{
		
	}

	/// Refresh the RRsets
	void refresh_RRsets()
	{
		for (size_t i =0; i< _num_mRRsets; i++)
		{
			for(auto &RR:_mRRsets[i])
			{
				RR.clear();
			}
			mRRset().swap(_mRRsets[i]);
		}
		mRRsets().swap(_mRRsets);
		for (auto i = __numV; i--;)
		{
			FRset().swap(_FRsets[i]);
		}
		_num_mRRsets = 0;  // important
		for(auto &vec:vv_roots)
		{
			Nodelist().swap(vec);
		}
		for(auto &vec:vv_virtual_roots)
		{
			Nodelist().swap(vec);
		}
	}

	void refresh_FRmRRsets(int max_size)
	{
		// refresh FRsets
		// for (auto i =max_size; i< _num_mRRsets; i++)
		// {
		// 	for(auto &RR:_mRRsets[i])
		// 	{
		// 		for(auto &entry:RR)
		// 		{
		// 			_FRsets[entry.first].erase(i);
		// 		}
		// 	}
		// }
		for(auto i=0;i<__numV;i++)
		{
			auto &frset= _FRsets[i];
			auto it=lower_bound(frset.begin(), frset.end(), max_size);
			auto k= it- frset.begin();
			frset.resize(k);
		}
		for (auto i =max_size; i< _num_mRRsets; i++)
		{
			for(auto &RR:_mRRsets[i])
			{
				RR.clear();
			}
			mRRset().swap(_mRRsets[i]);
		}
		_mRRsets.resize(max_size);
		for (auto i =max_size; i< _num_mRRsets; i++)
		{
			Nodelist().swap(vv_roots[i]);
		}
		vv_roots.resize(max_size);
		for (auto i =max_size; i< _num_mRRsets; i++)
		{
			Nodelist().swap(vv_virtual_roots[i]);
		}
		vv_virtual_roots.resize(max_size);
		_num_mRRsets = max_size;  // important
	}

	void refresh_mRRFRsets()
	{
		for (int i =0;i< __numV; i++)
		{
			FRset().swap(_FRsets[i]);
		}
		for(auto &vec:vv_virtual_roots)
		{
			Nodelist().swap(vec);
		}
		// no need to refresh mRRsets, since it is never recorded in ending rounds
		_num_mRRsets = 0;
	}

	/// Release memory
	void release_memory()
	{
		refresh_RRsets();
		std::vector<bool>().swap(__vecVisitBool);
		Nodelist().swap(__vecVisitNode);
		FRsets().swap(_FRsets);
		Nodelist().swap(__vecParentNode);
		vector<int>().swap(__vecTree);
		vector<vector<int>>().swap(vv_virtual_roots);
	}

	/// Set cascade model
	void set_cascade_model(const string model)
	{
		_cascadeModel = model;
	}

	void out_PO()
	{
		std::ofstream po("/data/fc/graphInfo/new/sample_po.txt");
		for(long unsigned int i=0;i<PO.size();i++)
		{
			if(PO[i].size()>0)
				po<<i<<"; ";
			for(auto node:PO[i])
			{
				po<<node<<" ";
			} 
			po<<endl;
		}
		po.close();
	}

	void out_graph(Graph g)
	{
		std::ofstream out_g("/data/fc/graphInfo/new/test_graph.txt");
		for(long unsigned int k=0;k<g.size();k++)
		{
			out_g<<k<<": ";
			for(auto &node:g[k])
			{
				out_g<<node<<", ";
			}
			out_g<<endl;
		}
	}

};

using TmRRcollection = mRRcollection;
using PmRRcollection = std::shared_ptr<TmRRcollection>;