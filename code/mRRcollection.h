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
    vector<int> vecRoot_num;
	float rand_div=1.0;

	public:
	vector<vector<int>> PO;
	FRsets _FRsets;
	mRRsets _mRRsets;
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
		__vecVisitBool = std::vector<bool>(__numV, false);
		__vecTree = std::vector<int>(__numV, -1);
		__vecSeq = std::vector<int>(__numV, -1);
		__vecNewTree = std::vector<int>(__numV, -1);
		#ifdef debug
		__vecVisitNode = Nodelist(5*__numV);
		__vecParentNode = Nodelist(5*__numV);
		#else
		__vecVisitNode = Nodelist(__numV);
		#endif
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
	void build_n_mRRsets_tree(const ulint numSamples)
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
            vv_polluted_nodes.resize(numSamples);
            vecRoot_num.resize(numSamples);
		}
        std::mt19937 gen(rd());
        std::binomial_distribution<int> dist(num_revise_RR, residual);
        ceil_root_RR=dist(gen);
		for(ulint i=0;i<num_revise_RR;i++)  // build the basic information of previous mRR-sets, and update these mRR-sets
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
		for(ulint i=0;i<num_revise_RR;i++)
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
		mRR.resize(root_num);
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
			mRR[i].push_back({root});
		}
		for (auto &RR:mRR)
		{
			int layer=0;
			while(RR.size()>layer)
			{
				vint &layer_nodes= RR[layer], new_layer_nodes;
				for(const auto &node:layer_nodes)
				{
					for(const auto &nbrId : (R_graph)[node])
					{
						if (__vecVisitBool[nbrId] || (__Activated)[nbrId])
							continue;
						if (dsfmt_gv_genrand_open_close() > __Inv_inDeg[node])
							continue;
						new_layer_nodes.push_back(nbrId);
						__vecVisitBool[nbrId] = true;
						_FRsets[nbrId].push_back(mRRid);
					}
				}
				if(new_layer_nodes.size()>0)
				{
					RR.emplace_back(std::move(new_layer_nodes));
					layer++;
				}
				else
				{
					break;
				}
			}
		}
		for(const auto &RR:mRR)
		{
			for (const auto &layer : RR)
			{
				for(const auto &expand : layer)
				{
					__vecVisitBool[expand] = false;
				}
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
		int mRR_size=static_cast<int>(mRR.size());
		vint &v_roots=vv_virtual_roots[mRRid];
		ulint v_roots_size=v_roots.size();
        vint roots, del_roots; roots.reserve(v_roots_size+mRR_size); del_roots.reserve(mRR_size);
        mRRset mRR_copy(mRR_size);
		int min_tree=__numV, first_del_node=0, affected_layer=-1;
        bool find_del=false;
		for(int i=0;i<mRR_size;i++)  // mark previous roots, and traverse the mRRset. Traversing from the end is not necessary, since we need to know whether a node us already in the mRR if regenerating.
		{
			auto &RR= mRR[i];
			int RR_size=static_cast<int>(RR.size());
			for(int j=0;i<RR_size;j++)
			{
				auto &layer=RR[j];
				for(const auto &node:layer)
				{
					__vecTree[node] = i;
					if(!find_del && __Activated[node])
					{
						min_tree=i;
						affected_layer=j;
						find_del=true;
					}
				}
			}
            if(i>min_tree)  // only records roots after min_tree
            {
                if(!__Activated[RR[0][0]])
                {
                    roots.push_back(RR[0][0]);
                    __vecNewTree[RR[0][0]] = mRR_size;  // should not be added into some tree during the new exploration
                }
                RR.swap(mRR_copy[i]);
            }
		}
		mRR[min_tree].swap(mRR_copy[min_tree]);
		#ifndef NDEBUG
			if(mRR_copy[min_tree].empty())
			{
				cout<<"Error: mRR_copy[min_tree] is empty in mRR_update, mRRid="<<mRRid<<", min_tree="<<min_tree<<endl;
				exit(1);
			}
			if(del_nodes_check(mRRid, del_nodes))
			{
				cout<<"Error: del_nodes_check failed in mRR_update, mRRid="<<mRRid<<endl;
				exit(1);
			}
			if(v_roots_check(mRRid, string(__func__)+" end"))
			{
				exit(1);
			}
            if(min_tree>mRR_size)
            {
                cout<<"Error: min_tree > mRR_size in mRR_update, mRRid="<<mRRid<<", min_tree="<<min_tree<<", mRR_size="<<mRR_size<<endl;
                mRR_out(mRRid);
                exit(1);
            }
        #endif
        mRR.resize(min_tree);  // delete empty trees
		for(int i=int(v_roots_size-1);i>-1;i--)  // I did not record the idx of v_roots here like before
		{
			int root=v_roots[i];
            #ifndef NDEBUG
                if(__vecTree[root]<0)
                {
                    cout<<__func__<<": root="<<root<<" not in mRR, __vecTree[root]="<<__vecTree[root]<<", mRRid="<<mRRid<<endl;
                    exit(1);
                }
            #endif
			if(__Activated[root])  // is a del_node
			{
				v_roots.erase(v_roots.begin()+i);
				continue;
			}
			if(__vecTree[root]>=min_tree)  // only realize v_roots that are affected by del_nodes
			{
				roots.push_back(root); // to facilitate the regeneration process
				v_roots.erase(v_roots.begin()+i);  // this v_root will be realized and thus should be removed from the v_roots
                __vecNewTree[root] = mRR_size; // mark realized v_roots
			}
		}
        vvint last_RR; last_RR.assign(mRR_copy[min_tree].begin(), mRR_copy[min_tree].begin()+affected_layer+1);  // last_RR is the last RR that will be regenerated
		int pre_affected_layer=affected_layer;
		bool find_new_layer=false;
        for(int i=0;i<=affected_layer;i++)
        {
			if(find_new_layer)
			{
				break;  
			}
			auto &layer_nodes=last_RR[i];
			for(const auto &node:layer_nodes)
			{
				if(__vecNewTree[node]>-1)  // a v_root
				{
					find_new_layer=true;
					affected_layer=i;
					break;
				}
				__vecNewTree[node] = min_tree;
			}
        }
		auto &layer_nodes=last_RR[affected_layer];
		ulint layer_nodes_size=layer_nodes.size();
		for(ulint i=layer_nodes_size-1;i>-1;i--)
		{
			if(__Activated[layer_nodes[i]] || __vecNewTree[layer_nodes[i]]>-1)
			{
				layer_nodes.erase(layer_nodes.begin()+i);  // remove del_nodes and v_roots
			}
			else
			{
				__vecNewTree[layer_nodes[i]] = min_tree;
			}
		}
		while(last_RR.size()>affected_layer)
		{
			vint &layer_nodes=last_RR[affected_layer], new_layer_nodes;
			for(const auto &node:layer_nodes)
			{
				for(const auto &nbrId : (R_graph)[node])
				{
					if(__Activated[nbrId] || (__vecTree[nbrId]>-1 && __vecTree[nbrId]<min_tree) || __vecNewTree[nbrId]>-1)
						continue;
					if(dsfmt_gv_genrand_open_close() > __Inv_inDeg[node])
						continue;
					new_layer_nodes.push_back(nbrId);
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
					}
				}
			}
		}
		if(last_RR[0].size()>0)  // make sure it is not empty
		{
			mRR.emplace_back(std::move(last_RR));
			mRR_size=min_tree+1;
		}
		else
		{			
			mRR_size=min_tree;
		}
        ulint num_new_roots=roots.size();
        mRR.resize(mRR_size+num_new_roots);
        for(ulint i=0;i<num_new_roots;i++)
        {
            auto &RR=mRR[mRR_size+i];
            RR.push_back({roots[i]});
			int layer=0;
			while(RR.size()>layer)
			{
				vint &layer_nodes= RR[layer], new_layer_nodes;
				for(const auto &node:layer_nodes)
				{
					for(const auto &nbrId : (R_graph)[node])
					{
						if(__Activated[nbrId] || (__vecTree[nbrId]>-1 && __vecTree[nbrId]<min_tree) || __vecNewTree[nbrId]>-1)
							continue;
						if (dsfmt_gv_genrand_open_close() > __Inv_inDeg[node])
							continue;
						new_layer_nodes.push_back(nbrId);
						__vecNewTree[nbrId] = mRR_size+i;
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
						}
					}
				}
				if(new_layer_nodes.size()>0)
				{
					RR.emplace_back(std::move(new_layer_nodes));
					layer++;
				}
				else
				{
					break;
				}
			}
		}
        for(const auto &RR:mRR_copy)
        {
            for(const auto &layer:RR)
            {
				for(const auto &node:layer)
				{
					if(__vecNewTree[node]<0)  // previously in mRR but now not in mRR
					{
						auto &frset = _FRsets[node];
						auto it=lower_bound(frset.begin(), frset.end(), mRRid);
						if (it != frset.end() && *it == mRRid) 
						{
							frset.erase(it); 
						}
					}
					__vecTree[node] =-1;
				}
            }
        }
        // reset all static variables
		#ifndef NDEBUG
        if(v_roots_check(mRRid, string(__func__)+" end"))
        {
            exit(1);
        }
		#endif
        for(int i=0;i<min_tree;i++)
        {
            auto &RR=mRR[i];
            for(const auto &layer:RR)
            {
				for(const auto &node:layer)
				{
                	__vecTree[node] = -1;  // reset the tree id
				}
            }
        }
        for(const auto &root:roots)
        {
            __vecTree[root] = -1;
        }
        
        mRR_size=mRR.size();
        for(int i=min_tree;i<mRR_size;i++)
        {
            auto &RR=mRR[i];
            for(const auto &layer:RR)
            {
				for(const auto &node:layer)
				{
                	__vecNewTree[node] = -1;  // reset the tree id
				}
            }
        }
        vecRoot_num[mRRid]=mRR_size+v_roots.size();
		#ifndef NDEBUG
        if(synthetic_check(mRRid, string(__func__)+" end", 1,1,1,1,1,1))
        {
            exit(1);
        }
		#endif
    }

	void add_root(int mRRid, int num)
	{
        vecRoot_num[mRRid] += num;
		mRRset &mRR=_mRRsets[mRRid];
		vint roots, new_roots; roots.reserve(root_num); new_roots.reserve(num); mRR.reserve(root_num+num);
        auto &v_roots=vv_virtual_roots[mRRid];
		for (const auto &RR:mRR)  // mark previous roots, and traverse the mRRset
		{
			int root = RR[0][0];
			roots.push_back(root);
			__vecTree[root] = 1024;
			for(const auto &layer : RR)
			{
				for(const auto &node : layer)
				{
					__vecVisitBool[node] = true;
				}
			}
		}
        for(auto root:v_roots)
		{
			__vecTree[root]=1024; // a root
		}
		#ifndef NDEBUG
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
				auto &RR=mRR[mRR_size_1-1];
				RR.push_back({root});

				__vecVisitBool[root] = true;
				auto &frset = _FRsets[root];
				auto it=lower_bound(frset.begin(), frset.end(), mRRid);
                if(it!=frset.end() && *it==mRRid)
                {
                    cout<<"Error: mRRid="<<mRRid<<", nbrId="<<root<<" already in _FRsets."<<endl;
                    exit(1);
                }
				frset.insert(it, mRRid);

				int layer=0;
				while(RR.size()>layer)
				{
					vint &layer_nodes= RR[layer], new_layer_nodes;
					for(const auto &node:layer_nodes)
					{
						for(const auto &nbrId : (R_graph)[node])
						{
							if (__vecVisitBool[nbrId] || (__Activated)[nbrId])
								continue;
							if (dsfmt_gv_genrand_open_close() > __Inv_inDeg[node])
								continue;
							new_layer_nodes.push_back(nbrId);
							__vecVisitBool[nbrId] = true;
							auto &frset = _FRsets[nbrId];
							auto it=lower_bound(frset.begin(), frset.end(), mRRid);
							if(it!=frset.end() && *it==mRRid)
							{
								cout<<"Error: mRRid="<<mRRid<<", nbrId="<<nbrId<<" already in _FRsets."<<endl;
								exit(1);
							}
						}
					}
					if(new_layer_nodes.size()>0)
					{
						RR.emplace_back(std::move(new_layer_nodes));
						layer++;
					}
					else
					{
						break;
					}
				}
			}
		}
		#ifndef NDEBUG
        if(v_roots_check(mRRid, string(__func__)+" end"))
        {
            exit(1);
        }
		#endif
		for(const auto &RR:mRR)
		{
			for (const auto &layer : RR)
			{
				for(const auto &expand : layer)
				{
					__vecVisitBool[expand] = false;
				}
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
		#ifndef NDEBUG
        if(synthetic_check(mRRid, string(__func__)+" end",0,1,1,0,1,1))
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
            for(const auto &layer:RR)
            {
				for(const auto &node:layer)
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
			mRRset().swap(_mRRsets[i]);
		}
		mRRsets().swap(_mRRsets);
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
		}
		_mRRsets.resize(max_size);
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