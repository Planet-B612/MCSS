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
	vint __vecTree;
	vint __vecSeq;
	vint __vecNewTree;
	Nodelist __vecVisitNode;
	vector<vector<int>> vv_virtual_roots;
	float rand_div=1.0;

	public:
	vector<vector<int>> PO;
	FRsets _FRsets;
	mRRsets _mRRsets;
	mRRsets vec_mRR_layer;
	// #ifdef DEBUG
	vsint vec_hash_FR;
	vsint vec_hash_mRR;
	// #endif // !NDEBUG
    vector<int> vecRoot_num;
	ulint _num_mRRsets = 0;
	double decimal=1.0;
	double residual=0.0;
	int pre_root_num=0;
	Argument *__arg;
	string _cascadeModel;
	vector<double> __Inv_inDeg;
	string result;
	int num_update=0;
	int num_add_root=0;
    int num_delete_root=0;
	vint vec_round;
    vvint vv_polluted_nodes;
    std::random_device rd; // initialize random number generator


	double mRR_traversal_time = 0.0;

	explicit mRRcollection(Argument & arg)
	{
		__arg=&arg;
		__Inv_inDeg=arg.Inv_inDeg;
		__numV = arg.numV;
		_FRsets = FRsets(__numV);
		#ifdef DEBUG
		vec_hash_FR.resize(__numV);
		#endif 

		__vecVisitBool = std::vector<bool>(__numV, false);
		__vecTree = std::vector<int>(__numV, -1);
		__vecSeq = std::vector<int>(__numV, -1);
		__vecNewTree = std::vector<int>(__numV, -1);
		__vecVisitNode = Nodelist(__numV);
		_cascadeModel=arg.model;
		result=arg.result_dir;
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
		 	++counter_real;
		 	(__Activated)[seed] = true;
		 	__vecVisitNode[numVisitNode++]=seed;
            for(const auto &rrid: _FRsets[seed])
            {
				if(rrid>=_num_mRRsets) continue;
                vv_polluted_nodes[rrid].push_back(seed);
            }
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
                for(const auto &rrid: _FRsets[v])
                {
					if(rrid>=_num_mRRsets) continue;
                    vv_polluted_nodes[rrid].push_back(v);
                }
			}
		}
		Nodelist temp(__vecVisitNode.begin(), __vecVisitNode.begin()+numVisitNode); 
		activated_nodes.insert((activated_nodes).end(),temp);
		return counter_real;
	}

 #include "test_ic.h"


	/// Generate a set of n mRR sets
	void build_n_mRRsets_tree(const ulint numSamples, const int pre_theta)
	{
		int floor_root_RR=0;
		int ceil_root_RR=0;  // the number mRR-sets with root number root_num+1 in the previous revisable mRR-sets
		const ulint prevSize = _num_mRRsets;  // previous total number of mRR-sets
		vector<bool> mRR_mark(prevSize, false); // false indicates this mRR is not directly reused.
		ulint num_revise_RR=(prevSize>numSamples?numSamples:prevSize);  // the number of mRR-sets that will be revised (revisable mRR-sets)
        vint vec_rootnum_RRid, vec_rootnum_1_RRid; 
        vec_rootnum_RRid.reserve(num_revise_RR); vec_rootnum_1_RRid.reserve(num_revise_RR);
		if(prevSize<numSamples)
		{
			vv_virtual_roots.resize(numSamples);
			_mRRsets.resize(numSamples);
			vec_mRR_layer.resize(numSamples);
			#ifdef DEBUG
			vec_hash_mRR.resize(numSamples);
			#endif // !NDEBUG
            vv_polluted_nodes.resize(numSamples);
            vecRoot_num.resize(numSamples);
		}
        std::mt19937 gen(rd());
        std::binomial_distribution<int> dist(num_revise_RR, residual);
        ceil_root_RR=dist(gen);
		#ifdef DEBUG
		int original_ceil_root_RR=ceil_root_RR;
		#endif
		for(ulint i=pre_theta;i<num_revise_RR;i++)  // build the basic information of previous mRR-sets, and update these mRR-sets
		{
            vint &polluted_nodes=vv_polluted_nodes[i];
            if(polluted_nodes.size()>0)
            {
                mRR_update(i, polluted_nodes);
                polluted_nodes.clear();
				num_update++;
            }
			// The root info should be recorded after the mRR-sets are updated.
            if(vecRoot_num[i]==root_num)
            {
                vec_rootnum_RRid.push_back(i);
            }
            else if(vecRoot_num[i]==root_num+1)
            {
                vec_rootnum_1_RRid.push_back(i);
            }
		}
		floor_root_RR=num_revise_RR-ceil_root_RR;
		for(auto mRRid:vec_rootnum_RRid)  // directly reuse updated previous mRR-sets with root number root_num, if there is any such mRR-sets
		{
			if(floor_root_RR>0)  // if still need floor_root_RR
			{
				mRR_mark[mRRid]=true;
				floor_root_RR--;
			}
		}
		for(auto mRRid:vec_rootnum_1_RRid)  // directly reuse updated previous mRR-sets with root number root_num+1
		{
			if(ceil_root_RR>0)
			{
				mRR_mark[mRRid]=true;
				ceil_root_RR--;
			}
		}
		int root_diff=0;
		for(ulint i=pre_theta;i<num_revise_RR;i++)
		{
			if(mRR_mark[i]==false)  // for mRR-sets that have not been directly reused
			{
				if(floor_root_RR>0)  // derive floor-rooted mRR first
				{
					root_diff=vecRoot_num[i]-root_num;
					if(root_diff>0)
					{
						delete_root(i,root_diff);
						num_delete_root++;
					}
					else
					{
						// output_info(i,true);
						add_root(i,-root_diff);
						// output_info(i,false);
						num_add_root++;
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
						num_delete_root++;
					}
					else
					{
						// output_info(i,true);
						add_root(i,-root_diff);
						// output_info(i,false);
						num_add_root++;
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

	int build_one_mRRset_tree(int mRRid, int root_num, double residual)
	// Each adj_list in in the form of adjacency list, so that the first node of each entry automatically constitutes the original __vecVisitNode
	{
		int root;
		root_num += (dsfmt_gv_genrand_open_close() <= residual);
        vecRoot_num[mRRid]=root_num;
		mRRset &mRR=_mRRsets[mRRid];
		#ifdef DEBUG
		sint &mRR_hash= vec_hash_mRR[mRRid];
		#endif // !NDEBUG
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
			__vecVisitBool[root] = true; // only record the state of roots, but do not push the root into the queue, since we are not going to diffuse here.
			_FRsets[root].push_back(mRRid);
			mRR[i].push_back(root);
			#ifdef DEBUG
			mRR_hash.insert(root);
			vec_hash_FR[root].insert(mRRid);
			#endif // !NDEBUG
		}
		for (int i = 0; i < root_num; i++)
		{
			auto &RR = mRR[i];
			int layer_start = 0, layer_end = 1;
			while(layer_start < layer_end) 
			{ 
				vec_RR_layer[i].push_back(layer_start);
				for (int j = layer_start; j < layer_end; j++)
				{
					int node = RR[j];
					for(const auto &nbrId : (R_graph)[node])
					{
						if (__vecVisitBool[nbrId] || (__Activated)[nbrId])
							continue;
						if (dsfmt_gv_genrand_open_close() > __Inv_inDeg[node])
							continue;
						RR.push_back(nbrId);
						#ifdef DEBUG
						mRR_hash.insert(nbrId);
						vec_hash_FR[nbrId].insert(mRRid);
						#endif
						__vecVisitBool[nbrId] = true;
						_FRsets[nbrId].push_back(mRRid);
					}
				}
				layer_start = layer_end;  // update the start index of the next layer
				layer_end = RR.size();  // update the end index of the next layer
			}
		}	

		for(const auto &RR:mRR)
		{
			for (const auto &node : RR)
			{
				__vecVisitBool[node] = false;
			}
		}
		// vec_value_check(__vecVisitBool, false, 1, string(__func__) + "beg="+to_string(0)+" __vecVisitBool includes TRUE values.");
		// FR_sorted_check(__func__);
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
		mRRset &mRR=_mRRsets[mRRid];
		mRRset &mRR_layer=vec_mRR_layer[mRRid];
		#ifdef DEBUG
		mRRset pre_mRR=mRR;
		mRRset pre_mRR_layer=mRR_layer;
		sint &mRR_hash=vec_hash_mRR[mRRid];
		#endif // !NDEBUG
		int mRR_size=static_cast<int>(mRR.size());
		vint &v_roots=vv_virtual_roots[mRRid];
		ulint v_roots_size=v_roots.size();
        vint roots, del_roots; roots.reserve(v_roots_size+mRR_size); del_roots.reserve(mRR_size);

		for(int i=static_cast<int>(v_roots_size-1);i>-1;i--)
		{
			int root=v_roots[i];
			if(__Activated[root])  // is a del_node
			{
				v_roots.erase(v_roots.begin()+i);
				continue;
			}
            __vecNewTree[root] = mRR_size; // mark realized v_roots
		}
		v_roots_size=v_roots.size();

		int min_tree=__numV, affected_layer_idx=__numV, first_del_idx=__numV, min_tree_RR_size;
        bool find_del=false;
		for(int i=0;i<mRR_size;i++)  // mark previous roots, and traverse the mRRset. Traversing from the end is not necessary, since we need to know whether a node us already in the mRR if regenerating.
		{
			auto &RR= mRR[i];
			min_tree_RR_size=static_cast<int>(RR.size());
			for(int j=0;j<min_tree_RR_size;j++)
			{
				__vecTree[RR[j]] = i;
				if(!find_del && __Activated[RR[j]])
				{
					first_del_idx=j;
					min_tree=i;
					find_del=true;
				}
				if(find_del)
				{
					break;
				}
			}			
			if(find_del)
			{
				break;
			}         
		}
		auto &min_tree_layer=mRR_layer[min_tree];
		affected_layer_idx=upper_bound(min_tree_layer.begin(), min_tree_layer.end(), first_del_idx)-min_tree_layer.begin()-1;
		vint &min_tree_RR=mRR[min_tree];	
		int affected_layer_beg=min_tree_layer[affected_layer_idx];
		int affected_next_layer_beg, min_tree_layer_size=static_cast<int>(min_tree_layer.size());
		if(affected_layer_idx==min_tree_layer_size-1)
		{
			affected_next_layer_beg=min_tree_RR_size;
		}
		else
		{
			affected_next_layer_beg=min_tree_layer[affected_layer_idx+1];
		}	
		if(affected_next_layer_beg>min_tree_RR_size)
		{
			cout<<"Abnormal resize in Line: "<<__LINE__<<endl;
		}
		for(int i=affected_next_layer_beg; i<min_tree_RR_size; i++)
		{
			__vecTree[min_tree_RR[i]]=__numV;
			#ifdef DEBUG
			mRR_hash.erase(min_tree_RR[i]);
			#endif
		}
        mRRset mRR_copy(mRR_size-min_tree);
		if(affected_next_layer_beg > min_tree_RR_size)
		{
			cout << "Error: affected_next_layer_beg > min_tree_RR_size in mRR_update, mRRid=" << mRRid << ", affected_next_layer_beg=" << affected_next_layer_beg << ", min_tree_RR_size=" << min_tree_RR_size << endl;
			exit(1);
		}
		mRR_copy[0].reserve(min_tree_RR_size-affected_next_layer_beg);
		mRR_copy[0]=vector<int>(std::make_move_iterator(min_tree_RR.begin()+affected_next_layer_beg), std::make_move_iterator(min_tree_RR.end()));
		
		min_tree_RR.resize(affected_next_layer_beg);
		min_tree_layer.resize(affected_layer_idx);

		for(int i=min_tree+1;i<mRR_size;i++)
		{
			auto &RR= mRR[i];
			if(!__Activated[RR[0]])
			{
				roots.push_back(RR[0]);
				__vecNewTree[RR[0]] = mRR_size;  // should not be added into some tree during the new exploration
			}
			for(const auto &node:RR)
			{
				#ifdef DEBUG
				mRR_hash.erase(node);
				#endif // !NDEBUG
				__vecTree[node] = i;
			}
			RR.swap(mRR_copy[i-min_tree]);
			mRR_layer[i].clear();
		}
		
		for(int i=static_cast<int>(v_roots_size-1);i>-1;i--)
		{
			int root=v_roots[i];
			if(__vecTree[root]>min_tree)  // delete realized v_roots
			{
				roots.push_back(root);
				v_roots.erase(v_roots.begin()+i);
			}
		}
		for(int i=affected_next_layer_beg-1;i>=affected_layer_beg;i--)
		{
			int node=min_tree_RR[i];
			if(__Activated[node]) // only consider the activated nodes in this layer now
			{
				__vecTree[node]=-1;  // necessary to make __vecTree all -1
				min_tree_RR.erase(min_tree_RR.begin()+i);  // remove del_nodes
				#ifdef DEBUG
				mRR_hash.erase(node);
				#endif // !NDEBUG
				affected_next_layer_beg--;
			}
			else
			{
				__vecTree[node] = min_tree;
				__vecNewTree[node] = min_tree;
			}
		}
		#ifdef DEBUG
			// if(del_nodes_check(mRRid, del_nodes))
			// {
			// 	cout<<__LINE__<<": Error: del_nodes_check failed in mRR_update, mRRid="<<mRRid<<endl;
			// 	exit(1);
			// }
			if(v_roots_check(mRRid, string(__func__)+" end"))
			{
				cout<<__LINE__<<": Error: v_roots_check failed in mRR_update, mRRid="<<mRRid<<endl;
				exit(1);
			}
            if(min_tree>=mRR_size)
            {
                cout<<__LINE__<<": Error: min_tree > mRR_size in mRR_update, mRRid="<<mRRid<<", min_tree="<<min_tree<<", mRR_size="<<mRR_size<<endl;
                mRR_out(mRRid);
                exit(1);
            }
        #endif
		int layer_start = affected_layer_beg, layer_end=affected_next_layer_beg, expand;
		// min_tree_layer.pop_back();
		while(layer_start < layer_end)
		{
			min_tree_layer.push_back(layer_start);
			for (int j = layer_start; j < layer_end; j++)
			{
				expand = min_tree_RR[j];
				for (const auto &nbrId : (R_graph)[expand])
				{
					if(__Activated[nbrId] || (__vecTree[nbrId]>-1 && __vecTree[nbrId]<=min_tree) || __vecNewTree[nbrId]>-1)
						continue;
					if (dsfmt_gv_genrand_open_close() > __Inv_inDeg[expand])
						continue;
					min_tree_RR.push_back(nbrId);
					#ifdef DEBUG
					mRR_hash.insert(nbrId);
					#endif // !NDEBUG
					__vecNewTree[nbrId] = min_tree;  // mark the node as in the new mRR
					if(__vecTree[nbrId]<0)  // nbrId was not in this mRR previously
					{
						auto &frset = _FRsets[nbrId];
						auto it=lower_bound(frset.begin(), frset.end(), mRRid);
						if(it!=frset.end() && *it==mRRid)
						{
							cout<<"Error: mRRid="<<mRRid<<", nbrId="<<nbrId<<" already in _FRsets."<<endl;
							exit(1);
						}
						frset.insert(it, mRRid);
						#ifdef DEBUG
						vec_hash_FR[nbrId].insert(mRRid);
						#endif // !NDEBUG
					}
				}
			}
			layer_start = layer_end; // update the start index of the next layer
			layer_end = min_tree_RR.size();
		}
		if(min_tree_RR.size()<1)  // make sure it is not empty, // if the last RR root is a del_node
		{
			mRR_size=min_tree;
		}
		else
		{
			mRR_size=min_tree+1;
		}
		mRR.resize(mRR_size);  // remove the last RR
		mRR_layer.resize(mRR_size);
        ulint num_new_roots=roots.size();
		if(mRR_size+num_new_roots<mRR.size())
		{
			cout<<"Abnormal resize in Line: "<<__LINE__<<endl;
		}
        mRR.resize(mRR_size+num_new_roots);
		mRR_layer.resize(mRR_size+num_new_roots);
        for(ulint i=0;i<num_new_roots;i++)
        {
            auto &RR=mRR[mRR_size+i];
			auto &RR_layer=mRR_layer[mRR_size+i];
            RR.push_back(roots[i]);
			#ifdef DEBUG
			mRR_hash.insert(roots[i]);
			#endif // !NDEBUG
			int layer_start = 0, layer_end = 1, node;
			while(layer_start < layer_end)
			{
				RR_layer.push_back(layer_start);
				for (int j = layer_start; j < layer_end; j++)
				{
					node = RR[j];
					for (const auto &nbrId : (R_graph)[node])
					{
						if(__Activated[nbrId] || (__vecTree[nbrId]>-1 && __vecTree[nbrId]<=min_tree) || __vecNewTree[nbrId]>-1)
							continue;
						if (dsfmt_gv_genrand_open_close() > __Inv_inDeg[node])
							continue;
						RR.push_back(nbrId);
						__vecNewTree[nbrId] = mRR_size+i;
						#ifdef DEBUG
						mRR_hash.insert(nbrId);
						#endif
						if(__vecTree[nbrId]<0)  // nbrId was not in this mRR previously
						{
							auto &frset = _FRsets[nbrId];
							auto it=lower_bound(frset.begin(), frset.end(), mRRid);
							if(it!=frset.end() && *it==mRRid)
							{
								cout<<"Error: mRRid="<<mRRid<<", nbrId="<<nbrId<<" already in _FRsets."<<endl;
								exit(1);
							}
							frset.insert(it, mRRid);
							#ifdef DEBUG
							vec_hash_FR[nbrId].insert(mRRid);
							#endif // !NDEBUG
						}
					}
				}
				layer_start = layer_end; // update the start index of the next layer
				layer_end = RR.size();
			}
		}
		#ifdef DEBUG
		if(FR_reverse_check_hash(mRRid, mRR_copy, "mRR_update end"))
		{
			cout<<__LINE__<<", Error: FR_reverse_check_hash"<<endl;
			exit(1);
		}
		#endif

        for(const auto &RR:mRR_copy)
        {
            for(const auto &node:RR)
            {
				if(__vecNewTree[node]<0)  // previously in mRR but now not in mRR
				{
					auto &frset = _FRsets[node];
					auto it=lower_bound(frset.begin(), frset.end(), mRRid);
					if (it != frset.end() && *it == mRRid) 
					{
						frset.erase(it);
						#ifdef DEBUG
						vec_hash_FR[node].erase(mRRid);
						#endif // !NDEBUG
					}
					else
					{
						cout<<__LINE__<<", Error: mRRid is not in the FRset of "<<node<<endl;
					}
				}
				__vecTree[node] =-1;
            }
        }
        // reset all static variables
		#ifdef DEBUG
        if(v_roots_check(mRRid, string(__func__)+" end"))
        {
            exit(1);
        }
		#endif
        mRR_size=mRR.size();
		int safe_size=min_tree;
		if(min_tree>=mRR_size) 
		{
			safe_size=mRR_size-1;
		}
        for(int i=0;i<=safe_size;i++)
        {
            auto &RR=mRR[i];
            for(const auto &node:RR)
            {
				__vecTree[node] = -1;  // reset the tree id
            }
        }
        
        for(int i=min_tree;i<mRR_size;i++)
        {
            auto &RR=mRR[i];
            for(const auto &node:RR)
            {
				__vecNewTree[node] = -1;  // reset the tree id
            }
        }
		for(const auto &root : v_roots)  // reset the v_roots
		{
			__vecNewTree[root] = -1;  // reset the tree id
		}
        vecRoot_num[mRRid]=mRR_size+v_roots.size();
		#ifdef DEBUG
        if(synthetic_check(mRRid, string(__func__)+" end", 1,1,1,0,1,1,1))
        {
            exit(1);
        }
		#endif
    }

	void add_root(int mRRid, int num)
	{
        vecRoot_num[mRRid] += num;
		mRRset &mRR=_mRRsets[mRRid];
		auto &mRR_layer = vec_mRR_layer[mRRid];
		#ifdef DEBUG
		sint &mRR_hash=vec_hash_mRR[mRRid];
		#endif // !NDEBUG
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
		#ifdef DEBUG
        if(v_roots_check(mRRid, string(__func__)+" end"))
        {
            exit(1);
        }
		#endif
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
		for(const auto &root : new_roots)  // traverse the mRRset
		{
			if(__vecVisitBool[root])
			{
				v_roots.push_back(root);
				continue;
			}
			else
			{
				auto mRR_size_1 = mRR.size()+1;
				mRR.resize(mRR_size_1);
				mRR_layer.resize(mRR_size_1);
				auto &RR=mRR[mRR_size_1-1];
				RR.push_back(root);
				#ifdef DEBUG
				mRR_hash.insert(root);
				#endif
				__vecVisitBool[root] = true;
				auto &frset = _FRsets[root];
				auto it=lower_bound(frset.begin(), frset.end(), mRRid);
                if(it!=frset.end() && *it==mRRid)
                {
                    cout<<__LINE__<<": Error: mRRid="<<mRRid<<", nbrId="<<root<<" already in _FRsets."<<endl;
                    exit(1);
                }
				frset.insert(it, mRRid);
				#ifdef DEBUG
				mRR_hash.insert(root);
				vec_hash_FR[root].insert(mRRid);
				#endif // !NDEBUG

				int layer_start = 0, layer_end = 1;
				while(layer_start < layer_end) 
				{
					mRR_layer[mRR_size_1-1].push_back(layer_start);
					for (int j = layer_start; j < layer_end; j++)
					{
						int node = RR[j];
						for(const auto &nbrId : (R_graph)[node])
						{
							if (__vecVisitBool[nbrId] || (__Activated)[nbrId])
								continue;
							if (dsfmt_gv_genrand_open_close() > __Inv_inDeg[node])
								continue;
							RR.push_back(nbrId);
							#ifdef DEBUG
							mRR_hash.insert(root);
							#endif
							__vecVisitBool[nbrId] = true;
							auto &frset=_FRsets[nbrId];
							auto it=lower_bound(frset.begin(), frset.end(),mRRid);
							if(it!=frset.end() && *it==mRRid)
							{
								cout<<__LINE__<<", Error: mRRid "<< mRRid<<" already in the FRset of node "<<node<<endl;
								exit(1);
							}
							else
							{
								frset.insert(it,mRRid);
								#ifdef DEBUG
								vec_hash_FR[nbrId].insert(mRRid);
								#endif
							}
						}
					}
					layer_start = layer_end;  // update the start index of the next layer
					layer_end = RR.size();  // update the end index of the next layer
				}
			}
		}
		#ifdef DEBUG
        if(v_roots_check(mRRid, string(__func__)+" end"))
        {
            exit(1);
        }
		mRRset a={{}};
		if(FR_reverse_check_hash(mRRid, a, "add_root"))
		{
			cout<<__LINE__<<", Error: FR_reverse_check_hash"<<endl;
			exit(1);
		}
		#endif
		for(const auto &RR:mRR)
		{
			for(auto &node : RR)
			{
				__vecVisitBool[node] = false;
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
		
        for(auto root:v_roots)
		{
			__vecTree[root] = -1;
		}
		#ifdef DEBUG
        if(synthetic_check(mRRid, string(__func__)+" end",0,1,1,0,1,1,1))
        {
            exit(1);
        }
		#endif
		return;
	}

	void delete_root(int mRRid, int num_del_roots)
	{
        vecRoot_num[mRRid] -= num_del_roots;
		mRRset &mRR=_mRRsets[mRRid];
		ulint mRR_size=mRR.size();
		vint &v_roots=vv_virtual_roots[mRRid], roots;
		ulint v_roots_size=v_roots.size();
        if(v_roots_size>= num_del_roots)  
        {
            v_roots.resize(v_roots_size-num_del_roots);
            num_del_roots=0;
        }
        else
        {
            v_roots.clear();
            num_del_roots -= v_roots_size;
        }
        for(int i=0;i<num_del_roots;i++)
        {
            auto &RR=mRR[mRR_size-1-i];
            for(const auto &node:RR)
            {
				auto &frset = _FRsets[node];
				auto it=lower_bound(frset.begin(), frset.end(), mRRid);
				if (it != frset.end() && *it == mRRid) 
				{
					frset.erase(it);
					#ifdef DEBUG
					vec_hash_mRR[mRRid].erase(node);
					vec_hash_FR[node].erase(mRRid);
					#endif // !NDEBUG
				}
				else
				{
					cout<<__LINE__<<", Error: "<<node<<" is not in "<<mRRid<<", when deleting it."<<endl;
				}
            }
        }
        mRR.resize(mRR_size-num_del_roots);
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
			for(auto &layer:vec_mRR_layer[i])
			{
				layer.clear();
			}
			mRRset().swap(vec_mRR_layer[i]);
			mRRset().swap(_mRRsets[i]);
		}
		mRRsets().swap(_mRRsets);
		mRRsets().swap(vec_mRR_layer);
		for (auto i = __numV; i--;)
		{
			FRset().swap(_FRsets[i]);
		}
		_num_mRRsets = 0;  // important
		for(auto &vec:vv_virtual_roots)
		{
			Nodelist().swap(vec);
		}
	}

	void refresh_FRmRRsets(int max_size)
	{
		for(auto i=0;i<__numV;i++)
		{
			auto &frset= _FRsets[i];
			auto it=lower_bound(frset.begin(), frset.end(), max_size);
			auto k= it- frset.begin();
			frset.resize(k);
		}
		for (ulint i =max_size; i< _num_mRRsets; i++)
		{
			mRRset().swap(_mRRsets[i]);
			mRRset().swap(vec_mRR_layer[i]);
		}
		_mRRsets.resize(max_size);
		vec_mRR_layer.resize(max_size);
		#ifdef DEBUG
		vec_hash_mRR.resize(max_size);
		#endif // !NDEBUG
		for (ulint i =max_size; i< _num_mRRsets; i++)
		{
			Nodelist().swap(vv_virtual_roots[i]);
		}
        vecRoot_num.resize(max_size);
        vv_polluted_nodes.resize(max_size);
		vv_virtual_roots.resize(max_size);
        vec_round.resize(max_size);
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