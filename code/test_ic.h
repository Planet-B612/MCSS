#pragma once
#include "Argument.h"

void test_mRR()
{
	gene_syn_mRR();
	vint del_nodes={16};
	for(const auto &node:del_nodes)
	{
		__Activated[node]=true;
	}
	vv_virtual_roots[0].push_back(1);
	// add_root(0, 2);
	// delete_root(0, 1);
	mRR_update(0,del_nodes);
}

void gene_syn_mRR()
{
	__Activated[17]=true;
	mRRset mRR;
	mRR.push_back({5,4,21,0,6,22,18,15,13});
	mRR.push_back({20,12,19,10,14,1,16,3,7});
	for(const auto &RR:mRR)
	{
		for(const auto &node:RR)
		{
			_FRsets[node].push_back(0);
		}
	}
	_mRRsets.push_back(mRR);
	vv_virtual_roots.resize(1);
	vv_polluted_nodes.resize(1);
	vecRoot_num={2};
}

bool synthetic_check(int mRRid, string str, int newtree, int tree, int visitBool, int seq, int fr)
{
	bool a=false,b=false,c=false,d=false,e=false;
	if(newtree) a=vec_value_check(__vecNewTree, -1, 1, str + " __vecNewTree includes non -1 values.");
	if(tree) b=vec_value_check(__vecTree, -1, 1, str + " __vecTree includes non -1 values.");
	if(visitBool) c=vec_value_check(__vecVisitBool, false, 1, str + " __vecVisitBool includes TRUE values.");
	if(seq) d=vec_value_check(__vecSeq, -1, 1, str + " __vecSeq includes non -1 values.");
	if(fr) e=FR_check(mRRid, str+" FR_check");
	if(a || b || c || d || e)
	{
		cout<<"Error in synthetic_check of "<<str<<endl;
		return true;
	}
	else
	{
		return false;
	}
}

bool v_roots_check(int mRRid, string str)
{
	for(const auto &root : vv_virtual_roots[mRRid])
	{
		if(__vecTree[root]<0 && __vecNewTree[root]<0 && __vecVisitBool[root]==false)
		{
			cout<<"Error in v_roots_check of "<<str<<", the root "<<root<<" is not in mRR "<<mRRid<<" or new mRR."<<endl;
			return true;
		}
	}
	return false;
}

bool del_nodes_check(int mRRid, vint &del_nodes)
{
	for(const auto &node:del_nodes)
	{
		if(__vecTree[node]<0)
		{
			cout<<"Error in del_nodes_check of mRR "<<mRRid<<", the del_node "<<node<<" is not in mRR or new mRR."<<endl;
			vec_out(del_nodes, "del_nodes: ");
			return true;
		}
	}
	return false;
}

template <typename T, typename T1>
bool vec_value_check(T &vec, T1 val, int equality, string str)  // equality: 1: should be equal to val, -1: shoud not equal to val, 2: should be greater than, -2: should be smaller than
{
	bool flag=false;
	Nodelist vec_ind;
	if(equality==1)
	{
		for(auto i=0;i<vec.size();i++)
		{
			auto this_val=vec[i];
			if(this_val!=val)
			{
				flag=true;
				vec_ind.push_back(i);
			}
		}
		if(flag)
		{
			std::fstream result_bk(result, ios::app);
			assert(!result_bk.fail());
			result_bk<<str<<" vec Value errors: "<<endl;
			for(auto i:vec_ind)
			{
				result_bk<<vec[i]<<", ";
			}
			result_bk<<endl;
			result_bk.close();
			vec_out(vec_ind, "vec_ind: ");
			return true;
		}
		return false;
	}
	else if(equality==-1)
	{
		for(auto i:vec)
		{
			if(i==val)
			{
				cout<<str<<" vec Value error: "<<i<<" = "<<val<<endl;
				vec_out(vec);
				// exit(0);
				return true;
			}
		}
		return false;
	}
	else if(equality==2)
	{
		for(auto i:vec)
		{
			if(i<=val)
			{
				cout<<str<<" vec Value error: "<<i<<" != "<<val<<endl;
				vec_out(vec);
				// exit(0);
				return true;
			}
		}
		return false;
	}
	else if(equality==-2)
	{
		for(auto i:vec)
		{
			if(i>=val)
			{
				cout<<str<<" vec Value error: "<<i<<" != "<<val<<endl;
				vec_out(vec);
				// exit(0);
				return true;
			}
		}
		return false;
	}
	return false;
}

bool FR_check(int rid, string str, Nodelist p_nodes={})
{
	mRRset &mRR=_mRRsets[rid];
	for(auto &adj_list:mRR)
	{
		for(auto &node:adj_list)
		{
			// auto node=entry.first;
			// if(_FRsets[node].find(rid)==_FRsets[node].end())
			auto &frset= _FRsets[node];
			auto it= lower_bound(frset.begin(), frset.end(), rid);
			if( it == frset.end() )
			{
				set_out({node});
				output_info(rid, false, p_nodes);
				cout<<str+" Error in FR_full_check, mRR "<<rid<<" contains the node "<<node<<"; but the mRRid is not in node's _FRsets"<<endl;
				// output_info(rid);
				return true;
			}
		}
	}
	// for(int i=0;i<__numV;i++)
	// {
	// 	// if(_FRsets[i].find(rid)!=_FRsets[i].end())
	// 	auto &frset= _FRsets[i];
	// 	auto it= lower_bound(frset.begin(), frset.end(), rid);
	// 	if( it != frset.end() )
	// 	{
	// 		bool find_it=false;
	// 		for(auto &adj_list:mRR)
	// 		{
	// 			// if(adj_list.find(i)==adj_list.end())
	// 			// {
	// 			// 	continue;
	// 			// }
	// 			// else
	// 			// {
	// 			// 	find_it=true;
	// 			// 	break;
	// 			// }
	// 		}
	// 		if(!find_it)
	// 		{
	// 			set_out({i});
	// 			// output_info(rid);
	// 			cout<<"Error in FR_full_check, _FRset[node] contains mRRid "<<rid<<", but the node "<<i<<" is not in any adj_list of the mRR "<<endl;
	// 			output_info(rid);
	// 			return i;
	// 		}
	// 	}
	// }
	return false;
}

bool FR_sorted_check(string str)
{
	int i=0;
	for(const auto &fr: _FRsets)
	{
		if(!is_sorted(fr.begin(), fr.end()))
		{
			cout<<"Error in FR_sorted_check of "+str+" in checking the "+to_string(i)+"-th FRset. The _FRsets is not sorted."<<endl;
			exit(1);
		}
		i++;
	}
	return 0;
}

bool FR_insert_check(int rid, int node)
{
	// if(_FRsets[node].find(rid)!=_FRsets[node].end())
	auto &frset= _FRsets[node];
	auto it= lower_bound(frset.begin(), frset.end(), rid);
	if (it != frset.end() && *it == rid) 		
	{
		set_out({node});
		// __log_message("INFO", __FILE__, __LINE__, __func__);
		output_info(rid);
		// exit(1);
		return 1;
	}
	return 0;
}

template <typename T>
void vec_out(T &vec, string str="")
{
	std::fstream result_bk(result, ios::app);
	assert(!result_bk.fail());
	result_bk<<str+"vec: ";
	for(auto i:vec)
	{
		result_bk<<i<<", ";
	}
	result_bk<<endl;
	result_bk.close();
	// cout<<str+"vec: ";
	// for(auto i:vec)
	// {
	// 	cout<<i<<", ";
	// }
	// cout<<endl;
}

void mRR_out(int mRRid, string str="")
{
	std::fstream result_bk(result, ios::app);
	assert(!result_bk.fail());
	result_bk<<str+"mRRset: "<<mRRid<<endl;
	auto &mRR=_mRRsets[mRRid];
	for(auto i=0;i<mRR.size();i++)
	{
		result_bk<<i<<"-th RR: ";
		for(auto j:mRR[i])
		{
			result_bk<<j<<", ";
		}
		result_bk<<endl;
	}
	result_bk.close();
}

void set_out(Nodelist p_nodes)
{
	// std::fstream result_bk(result, ios::out);
	(*__arg).result_bk.open(result);
	assert(!(*__arg).result_bk.fail());
	for(auto i:p_nodes)
	{
		(*__arg).result_bk<<i<<"'s hashset values: "<<endl;
		for(auto j:_FRsets[i])
		{
			(*__arg).result_bk<<j<<", ";
		}		
		(*__arg).result_bk<<endl;
	}
	(*__arg).result_bk.close();
}

bool output_info(int mRRid, bool erase=false, Nodelist p_nodes={})
{		
	// cout<<vec_virtual_roots[mRRid].size()<<endl;
	// if(p_nodes.size()==0 && vec_virtual_roots[mRRid].size()==0)
	// {
	// 	return true;
	// }
	if(erase)
		(*__arg).result_bk.open(result);
	else
		(*__arg).result_bk.open(result, ios::app);
	assert(!(*__arg).result_bk.fail());
	(*__arg).result_bk<<"==========================================================="<<endl;
	(*__arg).result_bk<<"The vec_virtual_roots is: "<<endl;
	for(auto i:vv_virtual_roots[mRRid])
	{
		(*__arg).result_bk<<i<<", ";
	}
	(*__arg).result_bk<<endl;
	(*__arg).result_bk<<"The p_nodes are: "<<endl;
	for(auto node:p_nodes)
	{
		(*__arg).result_bk<<node<<", ";
	}
	(*__arg).result_bk<<endl;
	(*__arg).result_bk<<"The mRRid is: "<<mRRid<<endl;
	for(auto i=0;i<__numV;i++)
	{
		(*__arg).result_bk<<__vecVisitBool[i]<<", ";
	}
	(*__arg).result_bk<<endl;
	(*__arg).result_bk<<"The trees are: "<<endl;
	auto k=0;
	for(auto adj_list:_mRRsets[mRRid])
	{
		//if(k==0) 		continue;
		(*__arg).result_bk<<k<<"-th tree is: "<<endl;
		for(auto node:adj_list)
		{
			(*__arg).result_bk<<node<<", ";
		}
		(*__arg).result_bk<<endl;
		k++;
	}
	(*__arg).result_bk.close();
	return false;
}