#pragma once
#include "mRRcollection.h"
// #include "MChainCollection.h"
#include "CommonStruc.h"
#include "CommonFunc.h"
#include <memory>
#include <algorithm>
#include <queue>	//priority_queue
#include <malloc.h>
#include "Memory.h"
#include "MemoryUsage.h"
#include "test_ic.h"
using namespace std;


class Algorithm
{
private:
	size_t num_RRsets = 0;
	float __eta;
	double delta = 0.01;
	float eps = 0.1;
	vector<vector<bool>> __vecCover; // record whether an FRset_real is covered by some see
	int __numV;
	double total_theta=0;
	string  _cascadeModel;
	vector<tuple<int, int, double, int>> ratio_plain, ratio_UB;
	vector<double> vec_UB; // store the upper bound of coverage for each node
	
	int counter=0;  // record the number of nodes being affected in total
	Nodelist seed_batch;
	int seed;
	int batch_size=1;
	int coverage=0;
	double cov_UB=0.0;
	double cov_LB=0.0;
	int deg=0;
	double approx=1.0;
	double theta=0;
	int ending_rnd=300000;  // not needed, if q_ratio is properly set. sample: 4, facebook: 80, dblp: 3000
	const int root_num_bound=250; // try to delete unnecessary mRR-sets when the number of roots in an mRR exceeds this value. sample:0, facebook: 25, dblp: 250
	const int window_size=5;  // sample:2, facebook: 3, dblp: 5
	bool in_ending_rnd=false;
	bool delete_extra_mRR=false;
	vector<int> vec_mRR_size;
	int window_beg=0;
	int max_size_within_window=0;
	// Nodelist num_deg_incremental;
	Nodelist vec_deg;
	vector<bool> RR_Mark;
	double a=0.0;
	int __dataset_No = 0;
	float __q_ratio = 0.0;
	// double total_build_seedset_time = 0;
	double pre_theta = 0.0;
	float __left_num = 600.0;
	float __over_pnodes = 10.0;

public:
	mRRcollection RR;

	Algorithm(Argument &arg):RR(arg)
	{
		_cascadeModel=arg.model;
		__numV=arg.numV;
		eps=arg.eps;
		delta=arg.delta;
		__eta=arg.eta_0;
		batch_size=arg.batch;
		round_num=0;
		vec_mRR_size.resize(1e4,0);
		__q_ratio = arg.q_ratio;
		__dataset_No = arg.dataset_No;
		__left_num = arg.left_num;
		__over_pnodes = arg.over_pnodes;

	}
	
	~Algorithm()
	{
	}

	static double logcnk(int n, int k) 
	{
		double ans = 0;
		for (int i = n - k + 1; i <= n; i++)
		{
			ans += log(i);
		}
		for (int i = 1; i <= k; i++)
		{
			ans -= log(i);
		}
		return ans;
	}

	void build_seedset(int theta)
	{
		vec_deg.assign(__numV,0);
		seed_batch.clear();
		ratio_plain.clear();
		ratio_UB.clear(); 
		RR_Mark.assign(theta,false);
		if(in_ending_rnd)
		{
			if(use_UB)
			{
				ratio_UB.reserve(__numV);
				for (int i = (__numV); i--;) 
				{
					if((__Activated)[i]) continue;  // activated nodes should not be considered
					deg = RR._FRsets[i].size();  //The number of RR-sets covered by i.
					vec_deg[i] = deg;
					ratio_plain.push_back(make_tuple(i,deg,1.0*deg/(cost[i]),0));  // push back the plain ratio first
					double deg_UB=deg+a+sqrt(2.0*a*deg+1.0*a*a);
					vec_deg[i] = deg;
					ratio_UB.push_back(make_tuple(i,deg,1.0*deg_UB/cost[i],0));  // do not push back, otherwise this vector will be very long
				}
			}
			else 
			{
				for (int i = (__numV); i--;)  // initialize the benefit-to-cost ratio of each user
				{
					if((__Activated)[i]) continue;  // activated nodes should not be considered
					deg = RR._FRsets[i].size();  //The number of RR-sets covered by i.
					vec_deg[i] = deg;
					ratio_plain.push_back(make_tuple(i,deg,1.0*deg/cost[i],0));  // do not push back, otherwise this vector will be very long
				}
			}
		}
		else
		{
			for(auto i=0;i<__numV;i++)  // in the same round, only newly updated mRR-sets will contribute to vec_deg
			{
				auto &frset= RR._FRsets[i];
				auto it= lower_bound(frset.begin(), frset.end(), theta);
				auto k=it-frset.begin();  
				vec_deg[i] = k;  // The number of RR-sets covered by i.
			}

			if(use_UB)
			{
				ratio_UB.reserve(__numV);
				for(int i=0;i<(__numV);i++)
				{
					if((__Activated)[i]) continue;  // activated nodes should not be considered
					ratio_plain.push_back(make_tuple(i,vec_deg[i],1.0*vec_deg[i]/(cost)[i],0));  // push back the plain ratio first
					double deg_UB=1.0*(vec_deg[i]+a+sqrt(2.0*a*vec_deg[i]+1.0*a*a));
					ratio_UB.push_back(make_tuple(i,vec_deg[i],deg_UB/(cost)[i],0));
				}
			}
			else 
			{
				for(int i=0;i<(__numV);i++)
				{
					if((__Activated)[i]) continue;  // activated nodes should not be considered
					ratio_plain.push_back(make_tuple(i,vec_deg[i],1.0*vec_deg[i]/(cost)[i],0));
				}
			}
		}

		make_max_heap(ratio_plain);
		coverage=0;  // reset coverage in each trial
		for(int i =0;i<batch_size;i++)  // select k seeds
		{
			while(get<3>(ratio_plain[0])!=i)
			{
				seed=get<0>(ratio_plain[0]);
				// int nodeDeg=get<1>(ratio_plain[0]);
				int nodeDeg=vec_deg[seed];  // RR_Mark is shared throughout the selection of this batch. Thus, true states will be taken account repeatedly. Each time the number of true states is actually the total number.
				auto pre_nodeDeg=nodeDeg;
				for(int RRId: RR._FRsets[seed])
				{
					// if((RR_Mark[RRId]==true)||(RRId>=theta)) 	nodeDeg--;
					if((RR_Mark[RRId]==true) && RRId<theta) 	nodeDeg--;  // may need to check RRId<theta, since acessing a value outside RR_Mark is permitted in C++  // no need to check RRId>=theta, since nodeDeg only counts deg in current mRR-sets
				}
				if(nodeDeg>RR._num_mRRsets) 
				{
					cout<<"error, nodeDeg>num_mRRsets: "<<nodeDeg<<",  "<<RR._num_mRRsets<<". The original node_deg is "<<pre_nodeDeg<<endl; exit(1);
				}
				tuple<int,int, double, int> updated_node=make_tuple(seed,nodeDeg, 1.0*nodeDeg/(cost)[seed],i);
				max_heap_replace_max_value(ratio_plain, updated_node);
			}
			seed=get<0>(ratio_plain[0]);
			seed_batch.push_back(seed);
			coverage+=get<1>(ratio_plain[0]);
			for(auto &rr:RR._FRsets[seed])
			{
				if(RR_Mark[rr] || rr>=theta) continue;  // only mark currently used mRR-sets
				RR_Mark[rr]=true;
			}
			tuple<int, int, double, int> disable_node=make_tuple(seed, 0, -1.0, i);
			max_heap_replace_max_value(ratio_plain, disable_node);
		}
	}

	void build_max_single_seed(int theta){
		// RR.output_info(1);
		vec_deg.assign(__numV,0);
		seed_batch.clear();
		ratio_plain.clear();ratio_UB.reserve(__numV);		
		RR_Mark.assign(theta,false);
		vec_UB.resize(__numV);
		fill(vec_UB.begin(), vec_UB.end(), 0.0);
		if(in_ending_rnd)
		{
			if(use_UB)
			{
				ratio_UB.reserve(__numV);
				for (int i = (__numV); i--;) 
				{
					if((__Activated)[i]) continue;  // activated nodes should not be considered
					deg = RR._FRsets[i].size();  //The number of RR-sets covered by i.
					ratio_plain.push_back(make_tuple(i,deg,1.0*deg/(cost[i]),0));  // push back the plain ratio first
					double deg_UB=deg+a+sqrt(2.0*a*deg+1.0*a*a);
					ratio_UB.push_back(make_tuple(i,deg,1.0*deg_UB/cost[i],0));  // do not push back, otherwise this vector will be very long
				}
			}
			else 
			{
				for (int i = (__numV); i--;)  // initialize the benefit-to-cost ratio of each user. From 0 to numV-1
				{
					if((__Activated)[i]) continue;  // activated nodes should not be considered
					deg = RR._FRsets[i].size();  //The number of RR-sets covered by i.
					double deg_UB=deg+a+sqrt(2.0*a*deg+1.0*a*a);
					vec_UB[i]=deg_UB;
					ratio_plain.push_back(make_tuple(i,deg,1.0*deg/cost[i],0));  // do not push back, otherwise this vector will be very long
				}
			}
		}
		else
		{
			for(auto i=0;i<__numV;i++)  // in the same round, only newly updated mRR-sets will contribute to vec_deg
			{
				auto &frset= RR._FRsets[i];
				auto it= lower_bound(frset.begin(), frset.end(), theta);
				auto k=it-frset.begin();  
				vec_deg[i] = k;  // The number of RR-sets covered by i.
			}
			if(use_UB)
			{
				ratio_UB.reserve(__numV);
				for(int i=0;i<(__numV);i++)
				{
					if((__Activated)[i]) continue;  // activated nodes should not be considered
					ratio_plain.push_back(make_tuple(i,vec_deg[i],1.0*vec_deg[i]/(cost)[i],0));  // push back the plain ratio first
					double deg_UB=1.0*(vec_deg[i]+a+sqrt(2.0*a*vec_deg[i]+1.0*a*a));
					ratio_UB.push_back(make_tuple(i,vec_deg[i],deg_UB/(cost)[i],0));
				}
			}
			else 
			{
				for(int i=0;i<(__numV);i++)
				{
					if((__Activated)[i]) continue;  // activated nodes should not be considered
					deg = vec_deg[i];
					double deg_UB=deg+a+sqrt(2.0*a*deg+1.0*a*a);
					vec_UB[i]=deg_UB;
					ratio_plain.push_back(make_tuple(i,deg,1.0*deg/(cost)[i],0));
				}
			}
		}
		double max_ratio = 0.0;
		int seed = 0;
		coverage=0;
		for (int i = 0; i< ratio_plain.size();i++){
			double cur_ratio = get<2>(ratio_plain[i]);
			if (cur_ratio > max_ratio){
				seed = get<0>(ratio_plain[i]);
				max_ratio = cur_ratio;
				coverage = get<1>(ratio_plain[i]);
			}
		}
		seed_batch.push_back(seed);
	}

	void OneRoundSelect()
	{
		// RR.out_a(__LINE__);
		// double __eps_hat, __eps_prime, __delta=1.0/(__numV_left);
		// __eps_prime=(1-__eps_hat)/(1+__eps_hat);
		// __eps_hat=(eps-__delta)/(1-__delta);
		// double theta_max=ceil(2.0*__numV_left*(1+__eps_prime/3.0)*log(6.0/__delta)/(__eps_prime*__eps_prime*(1-1.0/2.71828)));
		// double T=ceil(log(1.0*__numV_left/(__eps_prime*__eps_prime))/log(2));
		// theta=1.0*theta_max/T;
		// a=log(3.0*T/__delta);
		//===================================
		delta=eps/(100.0*(1-1/2.71828)*(1-eps)*__eta_left);
		double eps_hat=99.0*eps/(100.0-eps);
		const double alpha = sqrt(log(6.0 / delta));
		const double beta = sqrt((logcnk(__numV_left, batch_size) + log(6.0 / delta)) / approx);

		theta = 2 * (alpha + beta)* (alpha + beta);
		const double theta_max = 2 * __numV_left*(alpha + beta)*(alpha + beta) / eps_hat / eps_hat / 1.0 * (batch_size);

		const double i_max = ceil(log(__numV_left / batch_size / eps_hat / eps_hat) / log(2)) + 1;

		const double a1 = log(3 * i_max / delta) + logcnk(__numV_left, batch_size);	
		const double a2 = log(3 * i_max / delta);
		a = a2;
		//===================================
		
		bool set_q = false;
		RR.__q_ratio = __q_ratio;
		while(theta<theta_max)
		{	
		    /* change q_ratio when eta_left < 500 or to much activated nodes in pre round*/
			if (!set_q)
			{
				if(__eta_left < __left_num){
					RR.__q_ratio = 1e-8;
					set_q = true;
				}
			}
			// 	else if ( round_num < 1){
			// 		RR.__q_ratio = __q_ratio; //default q_ratio
			// 		set_q = true;
			// 	}
			// 	else
			// 	{				
			// 		//count interval activated nodes number
			// 		int num_activated_nodes = 0;
			// 		int cur_mRR_size = theta;
			// 		int last_overmRR_index = -1;
			// 		for (int i = round_num -1 ; i >= 0; i--) {
			// 			if (vec_mRR_size[i] > cur_mRR_size) {
			// 				last_overmRR_index = i;
			// 				break;
			// 			}
			// 		}

			// 		if(last_overmRR_index!=-1){
			// 			for(int i = last_overmRR_index + 1; i < activated_nodes.size(); i++){
			// 				num_activated_nodes += activated_nodes[i].size(); //activated_nodes initialize with a empty entry
			// 			}
			// 		}
			// 		if (num_activated_nodes > __over_pnodes){
			// 			RR.__q_ratio = 1e-8;
			// 			set_q = true;
			// 		}
			// 	}		
			// }
			// result << "round " << (round_num+1) <<" q_ratio value " << RR.__q_ratio << " theta " <<theta << endl;
			if(in_ending_rnd)
			{
				RR.build_n_mRRsets_fresh_vec(theta);
			}
			else
			{
				RR.build_n_mRRsets_tree(theta, 0);
			}
			// build_seedset(theta);		
			if (batch_size > 1) build_seedset(theta);
			else  build_max_single_seed(theta);
			bool select=false;
			//===================================
			// cov_UB=1.0*coverage/approx;
			// cov_UB=cov_UB+a+sqrt(2.0*a*cov_UB+a*a);
			// cov_LB=1.0*( coverage+2.0*a/3.0-sqrt(2.0*a*coverage+4.0*a*a/9.0) );
			// if(cov_LB>=(1-__eps_hat)*approx*cov_UB) {select=true;}
			//===================================
			cov_UB=coverage*1.0/approx+a2+sqrt(2.0*a2*coverage/approx+a2*a2);
			cov_LB=1.0*(coverage+2.0*a1/3.0-sqrt(2.0*a1*coverage+4.0*a1*a1/9.0));
			if(cov_LB>=(1-eps_hat)*approx*cov_UB) {select=true;}
			//===================================
			if(!select && batch_size < 2)
			{
				int i=0;
				double ratio = 1.0 * coverage /cost[seed_batch[0]];
				for(;i<__numV_left;i++)
				{
					if(__Activated[i]) continue; // skip already activated nodes
					if(ratio<vec_UB[i]*1.581976707/cost[i])
					{
						if(i!=seed_batch[0])	break;
					}
				}
				if(i==__numV_left)
				{
					select=true;
				}
				// else  // for the batch version, we let the budget in each round be the costs used by the candidate set.
				// {
				// 	if(coverage/approx)
				// }
			}
			if(select)
			{
				(seed_set).insert((seed_set).end(),seed_batch.begin(),seed_batch.end());
				total_theta+=theta;
				pre_theta = theta;
				return;
			}
			if(theta>=RR_thr)
			{
				theta+=RR_step;
			}
			else
			{
				theta*=2;
			}
		}
		// build_seedset(theta);
		if (batch_size > 1) build_seedset(theta);
		else  build_max_single_seed(theta);		
		cov_UB=1.0*coverage/approx;
		cov_UB=cov_UB+a+sqrt(2.0*a*cov_UB+a*a);
		cov_LB=1.0*( coverage+2.0*a/3.0-sqrt(2.0*a*coverage+4.0*a*a/9.0) );
		(seed_set).insert((seed_set).end(),seed_batch.begin(),seed_batch.end());
		total_theta+=theta;
		return;
	}

	tuple<int,int,int,int,double,double> AdaptiveSelect()
	{
		approx=1.0-power((1-1.0/batch_size),batch_size);
		std::ofstream result;
		string file_name = "../results/ourround/round_" + std::to_string(__dataset_No) + "_" + std::to_string(int(__eta*__numV))+ "_" + std::to_string(__q_ratio);
		result.open(file_name, ios::app);
		assert(!result.fail());
		auto single_start = std::chrono::high_resolution_clock::now();
		while((__eta_left)>0)
		{
			auto start = std::chrono::high_resolution_clock::now();
			// if((round_num)/5==0) malloc_trim(0);
			root_num=ceil(1.0*__numV_left/__eta_left);
			
			if(in_ending_rnd==false)
			{
				if((__eta_left)<((__numV)/ending_rnd))
				{
					in_ending_rnd=true;
					cout<<"Entering the ending round"<<endl;
					RR.refresh_RRsets();
					// cout<<"Refreshed the RR-sets complete"<<endl;
				}
			}
			else  // in an ending round, refresh mRR and FR-sets
			{
				RR.refresh_mRRFRsets();
				// cout<<"Refreshed mRRFR-sets in ending rounds success at round: "<<(round_num)<<endl;
			}

			if((__eta_left)<=batch_size)
			{
				seed_batch.clear();
				for(int i=0;i<(__eta_left);i++)
				{
					float min_cost=INT_MAX;
					int node=(__numV);
					for(int j=0;j<(__numV);j++)
					{
						if((__Activated)[j]) continue;
						if(min_cost>(cost)[j])
						{
							min_cost=(cost)[j];
							node=j;
						}
					}
					(__Activated)[node]=true;
					seed_batch.push_back(node);
				}
			}
			else
			{
				OneRoundSelect();
			}
			counter=RR.realization(seed_batch);
			auto now = std::chrono::high_resolution_clock::now();
			std::chrono::duration<double> elapsed = now - start;

			(__numV_left)-= counter;  // mRR-sets need to be updated;
			(__eta_left)-=counter;

			vec_mRR_size[round_num] = theta;

			// if(1.0*(__numV_left)/(__eta_left)>root_num_bound)
			// {
			// 	vec_mRR_size[window_beg+window_size]=theta;
			// 	if(max_size_within_window<=theta)  // if previous max size > theta, and the earliest size in the window is not the previous max_size, no thing needs to be done
			// 	{
			// 		max_size_within_window=theta;
			// 	}
			// 	else if(vec_mRR_size[window_beg]==max_size_within_window)  // the element removed from the window is not the previous maximum one
			// 	{
			// 		max_size_within_window=*max_element(vec_mRR_size.begin()+window_beg+1, vec_mRR_size.begin()+window_beg+1+window_size);
			// 	}
			// 	window_beg++;
			// 	if(RR._num_mRRsets>max_size_within_window)
			// 	{
			// 		// cout << "Truncating the mRR-sets at round "<<(round_num)<<endl;
			// 		RR.refresh_FRmRRsets(max_size_within_window);
			// 	}
			// }
			if(1.0*(__numV_left)/(__eta_left)>root_num_bound)
			{
				if (round_num < window_size){
					window_beg = 0;
				}
				else{
					window_beg = round_num - window_size;
				}
				auto iter = max_element(vec_mRR_size.begin()+window_beg, vec_mRR_size.begin()+round_num);
				max_size_within_window = *iter;
				if(RR._num_mRRsets> max_size_within_window)
				{
					// cout << "Truncating the mRR-sets at round "<<(round_num)<<endl;
					RR.refresh_FRmRRsets(max_size_within_window);
				}
			}
	
			round_num++;
			result<<(round_num)<<", \t"<<counter<<", \t"<<(__eta_left)<<", \t"<<1.0*(__numV_left)/(__eta_left)<<", \t theta = "<<theta<<"; \t"<< RR.num_update<<" = "<<num_gen<<"+"<<num_addback<<", \t"<< RR.num_add_root<<", \t"<<RR.num_delete<<", \t"<< disp_mem_usage()<<" MB, \t"<<elapsed.count() << " 秒"<<endl;  // the round that is currently running

		}
		auto single_end = std::chrono::high_resolution_clock::now();
		std::chrono::duration <double> single_elapsed = single_end - single_start;
		result.close();
		// cout << "build seed set traversal time " << total_build_seedset_time <<endl;
		// cout << "mRR traversal time " <<RR.mRR_traversal_time <<endl;
		cout << "Single time " << single_elapsed.count() << " s" << endl;
		cout << "Single spread " << (-(__eta_left))+ __eta*__numV <<endl;
		return make_tuple(total_theta, RR.num_update, RR.num_add_root, RR.num_delete,single_elapsed.count(),(-(__eta_left))+ __eta*__numV);
	}

	void output_mRR(int mRRid)
	{
		string res="../results/backup.txt";
		std::ofstream result_bk;
		result_bk.open(res, ios::app);
        assert(!result_bk.fail());
		// for(auto &adj_list:RR._mRRsets[mRRid])
		// {
		// 	for(auto &entry:adj_list)
		// 	{
		// 		result_bk<<entry.first<<", ";
		// 		for(auto &nbr:entry.second)
		// 		{
		// 			result_bk<<nbr<<", ";
		// 		}
		// 		result_bk<<endl;
		// 	}
		// }
	}

	void release_memory()
	{
		RR.release_memory();
	}

}; // cls

using TAlg = Algorithm;
using PAlg = std::shared_ptr<TAlg>;