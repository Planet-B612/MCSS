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
#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>
#include <iomanip>
#include <string>
// #include "test_ic.h"
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
	int total_theta = 0;
	string  _cascadeModel;
	vector<tuple<int, int, double, int>> ratio_plain, ratio_UB;
	// vector<double> vec_UB; // store the upper bound of coverage for each node
	// vector<double> vec_LB;
	vint vec_deg;
	string model="IC";
	
	int counter=0;  // record the number of nodes being affected in total
	vint seed_batch;
	int seed;
	int batch_size=1;
	int deg=0;
	double approx=1.0;
	int theta=0;
	int ending_rnd=3000000;  // not needed, if q_ratio is properly set. sample: 4, facebook: 80, dblp: 3000. Let it be a large value, so that it will never enter the ending round to generate fresh mRRsets.
	int root_num_bound=250; // try to delete unnecessary mRR-sets when the number of roots in an mRR exceeds this value. sample:0, facebook: 25, dblp: 250
	const int window_size=5;  // sample:2, facebook: 3, dblp: 5
	bool in_ending_rnd=false;
	bool delete_extra_mRR=false;
	vector<int> vec_mRR_num;
	int window_beg=0;
	ulint max_size_within_window=0;
	// vint num_deg_incremental;
	vector<bool> RR_Mark;
	double a_1=0.0, a_2=0.0;
	int __dataset_No = 0;
	int __eta_left_threshold = 5;  // the threshold of eta_left, below which the algorithm will stop
	int __regen_threshold = 0.2;  // the threshold of the ratio of polluted mRR-sets, below which the algorithm will regenerate mRR-sets
	
	double total_round_time = 0.0;
	double total_realizaiton_time = 0.0;
	double total_seed_selection_time = 0.0;
	float __left_num = 600.0;
	float __over_pnodes = 10.0;
	int _delta_amp=1;

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
		vec_mRR_num.resize(1e4,0);
		__dataset_No = arg.dataset_No;
		__eta_left_threshold = arg.eta_left_threshold;
		__regen_threshold = arg.regen_threshold;
		root_num_bound = arg.root_num_bound;
		// __left_num = arg.left_num;
		// __over_pnodes = arg.over_pnodes;
		// vec_UB= vector<double>(__numV, 0.0);
		// vec_LB= vector<double>(__numV, 0.0);
		vec_deg= vector<int>(__numV, 0);
		model=arg.model;
		_delta_amp=arg.delta_amp;
		_delta_amp=arg.delta_amp;
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

	/// the build_seedset function for the batch mode needs to be further examined. How to get the optimal upper bound set (although the objective function is submodular), and how to get the seed set(at least we can simply select the seed set based on the vanilla coverage).
	bool build_seedset(int theta, double ratio=1.0)
	{
		vec_deg.assign(__numV,0);
		seed_batch.clear();
		ratio_plain.clear();
		// ratio_UB.clear(); 
		RR_Mark.assign(theta,false);
		if(in_ending_rnd)
		{
			// ratio_UB.reserve(__numV);
			for (int i = (__numV); i--;) 
			{
				if((__Activated)[i]) continue;  // activated nodes should not be considered
				vec_deg[i] = RR._FRsets[i].size();  //The number of RR-sets covered by i.
				// double deg_UB=deg+a_2+sqrt(2.0*a_2*deg+1.0*a_2*a_2);
				ratio_plain.push_back(make_tuple(i,deg,1.0*vec_deg[i]/cost[i],0));  // do not push back, otherwise this vector will be very long
			}
		}
		else
		{
			for(auto i=0;i<__numV;i++)  // in the same round, only newly updated mRR-sets will contribute to vec_deg
			{
				if((__Activated)[i]) continue;
				auto &frset= RR._FRsets[i];
				vec_deg[i]= lower_bound(frset.begin(), frset.end(), theta)-frset.begin();
				ratio_plain.push_back(make_tuple(i,vec_deg[i],1.0*vec_deg[i]/(cost)[i],0));
			}
		}

		make_max_heap(ratio_plain);
		int coverage=0;  // reset coverage in each trial
		for(int i =0;i<batch_size;i++)  // select k seeds
		{
			while(get<3>(ratio_plain[0])!=i)
			{
				seed=get<0>(ratio_plain[0]);
				// int nodeDeg=get<1>(ratio_plain[0]);
				int nodeDeg=vec_deg[seed];  // RR_Mark is shared throughout the selection of this batch. Thus, true states will be taken account repeatedly. Each time the number of true states is actually the total number.
				for(int RRId: RR._FRsets[seed])
				{
					// if((RR_Mark[RRId]==true)||(RRId>=theta)) 	nodeDeg--;
					if(RRId>=theta) break;
					if((RR_Mark[RRId]==true)) 	nodeDeg--;  // may need to check RRId<theta, since acessing a value outside RR_Mark is permitted in C++  // no need to check RRId>=theta, since nodeDeg only counts deg in current mRR-sets
				}
				assert(nodeDeg>=0);
				// double nodeDeg_UB=1.0*(nodeDeg+a_2+sqrt(2.0*a_2*nodeDeg+a_2*a_2));
				tuple<int,int, double, int> updated_node=make_tuple(seed,nodeDeg, nodeDeg/(cost)[seed],i);
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
		return (coverage + 2.0 * a_1 / 3.0 - sqrt(2 * a_1 * coverage + 4.0* a_1 *a_1 / 9.0)) >= ratio * (coverage + a_2 + sqrt(2 * a_2 * coverage + a_2 * a_2));  // fix the number of seeds selected in each round, and assume the budget fits this round automatically
	}

	bool build_max_single_seed(int theta, double ratio=1.0)
	{
		double a_1_23 = 2.0 * a_1 / 3.0, a_1_49=4.0* a_1 *a_1 / 9.0, a_2_2=a_2 * a_2;		
		double seed_ratio_LB = 0.0, max_ratio_UB = 0.0, this_LB_ratio=0.0, this_UB_ratio=0.0, ratio_i=0.0;
		int seed = 0;
		seed_batch.clear();
		for(auto i=0;i<__numV;i++)  // in the same round, only newly updated mRR-sets will contribute to vec_deg
		{
			if((__Activated)[i]) continue;
			auto &frset= RR._FRsets[i];
			ratio_i=frset.size()/(cost[i]);
			if(ratio_i<seed_ratio_LB) continue;		
			int k= lower_bound(frset.begin(), frset.end(), theta)-frset.begin();
			this_LB_ratio = (k + a_1_23 - sqrt(2 * a_1 * k + a_1_49))/(cost[i]);
			this_UB_ratio = (k+a_2+sqrt(2*a_2*k+a_2_2))/(cost[i]);
			if(this_LB_ratio>seed_ratio_LB)  // find the max seed ratio
			{
				seed_ratio_LB=this_LB_ratio;
				seed=i;
			}
			if(this_UB_ratio>max_ratio_UB)  // find the max ratio
			{
				max_ratio_UB=this_UB_ratio;
			}
		}
		seed_batch.push_back(seed);
		return seed_ratio_LB>=ratio*max_ratio_UB;  // seed ratio and max ratio
	}

	void OneRoundSelect()
	{
		//===================================
		double theta_max, i_max, ratio;
		int pre_theta=0;
		if(batch_size>1)
		{			
			delta=_delta_amp*eps/(100.0*(1-1/2.71828)*(1-eps)*__eta_left);
			// delta=eps/(100.0*(1-1/2.71828)*(1-eps)*__eta_left);
			double eps_hat=99.0*eps/(100.0-eps);
			const double alpha = sqrt(log(6.0 / delta));
			const double beta = sqrt((logcnk(__numV_left, batch_size) + log(6.0 / delta)) / approx);

			theta = 2 * (alpha + beta)* (alpha + beta);
			theta_max = 2 * __numV_left*(alpha + beta)*(alpha + beta) / eps_hat / eps_hat / 1.0 * (batch_size);

			i_max = ceil(log(__numV_left / batch_size / eps_hat / eps_hat) / log(2)) + 1;

			a_1 = log(3 * i_max / delta) + logcnk(__numV_left, batch_size);
			a_2 = log(3 * i_max / delta);

			ratio=(1-eps_hat)*approx;
		}
		else
		{
			delta=_delta_amp*1.0/__numV_left;
			delta=_delta_amp*1.0/__numV_left;
			double eps_hat=(eps-delta)/(1-delta);
			double eps_prime=(1-eps_hat)/(1+eps_hat);
			theta_max=2*(1+eps_hat/3.0)*__numV_left*log(6.0/delta)/(eps_prime*eps_prime*(1-1/2.71828));
			theta=ceil(theta_max* eps_hat*eps_hat/(2*__numV_left));
			i_max=ceil(log(__numV_left / (eps_hat*eps_hat)) / log(2));

			a_1 = log(3 * i_max / delta) + log(__numV_left);
			a_2 = log(3 * i_max / delta);
			ratio=(1-eps_hat)*approx;
		}


		//===================================
		bool select=false;
		while(theta<theta_max)
		{
			if(in_ending_rnd)
			{
				RR.build_n_mRRsets_fresh_vec(theta);
			}
			else
			{
				RR.build_n_mRRsets_tree(theta, pre_theta);
			}
			auto seed_slection_start = std::chrono::high_resolution_clock::now();
			if (batch_size > 1) select=build_seedset(theta, ratio);
			else  select=build_max_single_seed(theta, ratio);
			auto seed_slection_end = std::chrono::high_resolution_clock::now();
			total_seed_selection_time += std::chrono::duration<double>(seed_slection_end - seed_slection_start).count();
			if(select) 
			{
				(seed_set).insert((seed_set).end(),seed_batch.begin(),seed_batch.end());
				total_theta+=theta;
				return;
			}
			pre_theta=theta;
			if(theta>=RR_thr)  // modify the increase of mRR-sets from double to linear
			{
				theta+=RR_step;
			}
			else
			{
				theta*=2;
			}
		}
		// build_seedset(theta);
		auto seed_slection_start = std::chrono::high_resolution_clock::now();
		if (batch_size > 1) build_seedset(theta);
		else  build_max_single_seed(theta);
		auto seed_slection_end = std::chrono::high_resolution_clock::now();
		total_seed_selection_time += std::chrono::duration<double>(seed_slection_end - seed_slection_start).count();
		(seed_set).insert((seed_set).end(),seed_batch.begin(),seed_batch.end());
		total_theta+=theta;
		return;
	}

	tuple<int,int,int,int,double,double> AdaptiveSelect()
	{
		approx=1.0-power((1-1.0/batch_size),batch_size);
// #ifdef DEBUG
// 		auto now = std::chrono::system_clock::now();
// 		time_t tt = std::chrono::system_clock::to_time_t(now);
// 		tm ltm;
// 		localtime_r(&tt, &ltm);
// 		std::stringstream ss;
// 		ss << std::put_time(&ltm, "%Y%m%d_%H%M%S"); 
// 		string time_str = ss.str();
// 		std::ofstream result;
// 		string file_name = "/home/cfeng/mRR_Regen/code/log_mine/round_" + std::to_string(__dataset_No) + "_" + std::to_string(static_cast<int>(__eta*__numV))+ "_" + std::to_string(batch_size) + "_" + std::to_string(eps)+time_str + ".txt";
// 		result.open(file_name, ios::app);
// 		if (result.fail()) 
// 		{
// 			std::cerr << "Failed to open file: " << strerror(errno) << std::endl;
// 			assert(false);
// 		}
// 		assert(!result.fail());
// 		// result << "start recording at " << std::chrono::system_clock::now() << std::endl;
// #endif
		auto single_start = std::chrono::high_resolution_clock::now();
		while((__eta_left)>0)
		{
			// cout<<"round_num = "<<round_num<<endl;
			auto start = std::chrono::high_resolution_clock::now();
			decimal = 1.0 * (__numV_left) / (__eta_left);
			root_num=floor(decimal);
			residual = decimal - root_num;  // in (0,1)
			
			// if(root_num > 30000)
			// {
			// 	in_ending_rnd=true;
			// 	// cout<<"Entering the ending round"<<endl;
			// 	RR.refresh_RRsets();
			// 	// cout<<"Refreshed the RR-sets complete"<<endl;
			// }
			
			auto round_start = std::chrono::high_resolution_clock::now();
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
					(__Activated)[node]=1;
					seed_batch.push_back(node);
				}
			}
			else
			{
				OneRoundSelect();
			}
			auto round_end = std::chrono::high_resolution_clock::now();
			std::chrono::duration<double> round_elapsed = round_end - round_start;
			total_round_time+=round_elapsed.count();

			auto realization_start = std::chrono::high_resolution_clock::now();
			if(in_ending_rnd)
			{
				counter=RR.realization_fresh_vec(seed_batch);
			}
			else
			{
				counter=RR.realization(seed_batch);
			}

			auto realization_end = std::chrono::high_resolution_clock::now();
			std::chrono::duration<double> realization_elapsed = realization_end - realization_start;
			total_realizaiton_time+=realization_elapsed.count();

			(__numV_left)-= counter;  // mRR-sets need to be updated;
			(__eta_left)-=counter;

			vec_mRR_num[round_num] = theta;

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
			if(decimal>root_num_bound && round_num > window_size && in_ending_rnd==false)
			{
				window_beg = round_num - window_size;
				max_size_within_window = *(max_element(vec_mRR_num.begin()+window_beg, vec_mRR_num.begin()+round_num));
				if(RR._num_mRRsets> max_size_within_window)
				{
					// cout << "Truncating the mRR-sets at round "<<(round_num)<<endl;
					RR.refresh_FRmRRsets(max_size_within_window);
				}
			}
	
			round_num++;
// #ifdef DEBUG
// 			result<<counter<<", \t"<<_eta_left<<", \t"<<1.0*(__numV_left+counter)/(__eta_left+counter)<<", \t θ = "<<theta<<"; \t update: "<< RR.num_update_this_round <<"\t add: "<< RR.num_add_root_this_round<<", \t"<< disp_mem_usage()<<" MB, \t"<<round_elapsed.count() << " s"<<endl;  // the round that is currently running
// #endif
			RR.num_update_this_round=0;
			RR.num_add_root_this_round = 0;
			RR.num_delete_root_this_round = 0;

		}
		auto single_end = std::chrono::high_resolution_clock::now();
		std::chrono::duration <double> single_elapsed = single_end - single_start;
// #ifdef DEBUG
// 		result.close();
// #endif
		// cout << "build seed set traversal time " << total_build_seedset_time <<endl;
		// cout << "mRR traversal time " <<RR.mRR_traversal_time <<endl;
		cout << "Single time " << single_elapsed.count() << " s" << endl;
		cout << "Single spread " << (-(__eta_left))+ __eta*__numV <<endl;
		cout << "Single build mRR set time " << RR.build_mRRset_time << endl;
		cout << "Single revise mRR set time " << RR.revise_mRRset_time << endl;
		cout << "Single round time " << total_round_time << " s" << endl;
		cout << "Single realization time " << total_realizaiton_time << " s" << endl;
		cout << "Single seed selection time " << total_seed_selection_time << " s" << endl;
		return make_tuple(total_theta, RR.num_update_this_round, RR.num_add_root_this_round, RR.num_delete_root_this_round, single_elapsed.count(),(-(__eta_left))+ __eta*__numV);
	}

	bool build_seedset_veri(int theta, double ratio=1.0)
	{
		vec_deg.assign(__numV,0);
		seed_batch.clear();
		ratio_plain.clear();
		RR_Mark.assign(theta,false);
		for (int i = (__numV); i--;) 
		{
			if((__Activated)[i]) continue;  // activated nodes should not be considered
			vec_deg[i] = RR._FRsets[i].size();  //The number of RR-sets covered by i.
			// double deg_UB=deg+a_2+sqrt(2.0*a_2*deg+1.0*a_2*a_2);
			ratio_plain.push_back(make_tuple(i,deg,1.0*vec_deg[i]/cost[i],0));  // do not push back, otherwise this vector will be very long
		}

		make_max_heap(ratio_plain);
		int coverage=0; 
		for(int i =0;i<batch_size;i++)
		{
			while(get<3>(ratio_plain[0])!=i)
			{
				seed=get<0>(ratio_plain[0]);
				int nodeDeg=vec_deg[seed];
				for(int RRId: RR._FRsets[seed])
				{
					if(RRId>=theta) break;
					if((RR_Mark[RRId]==true)) 	nodeDeg--;
				}
				assert(nodeDeg>=0);
				tuple<int,int, double, int> updated_node=make_tuple(seed,nodeDeg, nodeDeg/(cost)[seed],i);
				max_heap_replace_max_value(ratio_plain, updated_node);
			}
			seed=get<0>(ratio_plain[0]);
			seed_batch.push_back(seed);
			coverage+=get<1>(ratio_plain[0]);
			for(auto &rr:RR._FRsets[seed])
			{
				if(RR_Mark[rr] || rr>=theta) continue; 
				RR_Mark[rr]=true;
			}
			tuple<int, int, double, int> disable_node=make_tuple(seed, 0, -1.0, i);
			max_heap_replace_max_value(ratio_plain, disable_node);
		}
		vbool RR_mark_verify(theta,false);
		for(int i=0;i<batch_size;i++)
		{
			for(auto &rr:RR._FRsets_veri[seed_batch[i]])
			{
				RR_mark_verify[rr]=true;
			}
		}
		int coverage_verify=0;
		for(int i=0;i<theta;i++)
		{
			if(RR_mark_verify[i]) coverage_verify++;
		}
		return (coverage_verify + 2.0 * a_1 / 3.0 - sqrt(2 * a_1 * coverage_verify + 4.0* a_1 *a_1 / 9.0)) >= ratio * (coverage_verify + a_2 + sqrt(2 * a_2 * coverage_verify + a_2 * a_2));
	}

	tuple<int,int,int,int,double,double> AdaptiveIM()
	{
		auto single_start = std::chrono::high_resolution_clock::now();
		approx=1.0-power((1-1.0/batch_size),batch_size);
		int round=0;
		while((__eta_left)>0)
		{
			decimal = 1.0 * (__numV_left) / (__eta_left);
			root_num=floor(decimal);
			residual = decimal - root_num;  // in (0,1)
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
					(__Activated)[node]=1;
					seed_batch.push_back(node);
				}
			}
			else
			{
				double theta_max, i_max, ratio;
				delta=_delta_amp*eps/(100.0*(1-1/2.71828)*(1-eps)*__eta_left);
				// delta=eps/(100.0*(1-1/2.71828)*(1-eps)*__eta_left);
				double eps_hat=99.0*eps/(100.0-eps);
				const double alpha = sqrt(log(6.0 / delta));
				const double beta = sqrt((logcnk(__numV_left, batch_size) + log(6.0 / delta)) / approx);

				theta = 2 * (alpha + beta)* (alpha + beta);
				theta_max = 2 * __numV_left*(alpha + beta)*(alpha + beta) / eps_hat / eps_hat / 1.0 * (batch_size);

				i_max = ceil(log(__numV_left / batch_size / eps_hat / eps_hat) / log(2)) + 1;

				a_1 = log(3 * i_max / delta) + logcnk(__numV_left, batch_size);
				a_2 = log(3 * i_max / delta);

				ratio=(1-eps_hat)*approx;
				bool select=false;
				while(theta<theta_max)
				{
					RR.build_veri_mRRset_fresh(theta, true);
					RR.build_veri_mRRset_fresh(theta, false);
					select=build_seedset_veri(theta, ratio);
					if(select) 
					{
						(seed_set).insert((seed_set).end(),seed_batch.begin(),seed_batch.end());
						total_theta+=theta;
						break;
					}
					theta*=2;
				}
				if(!select && theta>=theta_max) 
				{
					build_seedset_veri(theta);
					(seed_set).insert((seed_set).end(),seed_batch.begin(),seed_batch.end());
				}
			}
			counter=RR.realization_fresh_vec(seed_batch);
			// cout<<"round "<<round<<", the number of influenced nodes is "<<counter<<endl;
			(__numV_left)-= counter;  // mRR-sets need to be updated;
			(__eta_left)-=counter;
			RR.refresh_mRRFRsets();
		}
		auto single_end = std::chrono::high_resolution_clock::now();
		std::chrono::duration <double> single_elapsed = single_end - single_start;
		cout << "Single time " << single_elapsed.count() << " s" << endl;
		cout << "Single spread " << (-(__eta_left))+ __eta*__numV <<endl;
		return make_tuple(total_theta, RR.num_update_this_round, RR.num_add_root_this_round, RR.num_delete_root_this_round, single_elapsed.count(),(-(__eta_left))+ __eta*__numV);
	}

	tuple<double,double,double,double,double> accuracy_verification()
	{
		double update_estimation, RR_estimation=0.0, MC_estimation=0.0, inf_LB=0.0, inf_UB=0.0;

		// std::sort(pol_node_num_for_accuracy_verification.begin(), pol_node_num_for_accuracy_verification.end());
		// int max_pol_node_num = pol_node_num_for_accuracy_verification.back();
		// if(pol_node_num_for_accuracy_verification+eta_for_verification>__numV || pol_node_num_for_accuracy_verification<0)
		// {
		// 	cout<<"The number of polluted nodes for accuracy verification is not valid, please check the input."<<endl;
		// 	return make_tuple(0.0, 0.0, 0.0, 0.0, 0.0);
		// }
		if(eta_for_verification<seed_num_for_accuracy_verification)
		{
			cout<<"The number of seeds for accuracy verification is larger than the number of nodes left, please check the input."<<endl;
			return make_tuple(0.0, 0.0, 0.0, 0.0, 0.0);
		}
		vector<bool> selected(__numV, false);
		vint seeds;  
		seeds.reserve(seed_num_for_accuracy_verification);
		for(int i=0;i<seed_num_for_accuracy_verification;i++)
		{
			int node = dsfmt_gv_genrand_uint32_range(__numV);  // generate a random node
			while(selected[node])
			{
				node = dsfmt_gv_genrand_uint32_range(__numV);  // generate a random node
			}
			selected[node] = true;
			seeds.push_back(node);
		}

		int pre_pol_node_num=0;		
		vint pol_nodes;
		for(int pol_node_num_for_accuracy_verification:pol_node_num_for_accuracy_verification_list)
		{
			cout << "pol_nodes: " << pol_node_num_for_accuracy_verification <<endl;
			pol_nodes.reserve(pol_node_num_for_accuracy_verification);
			for(int i=pre_pol_node_num;i<pol_node_num_for_accuracy_verification;i++)
			{
				int node = dsfmt_gv_genrand_uint32_range(__numV);  // generate a random node
				while(selected[node])
				{
					node = dsfmt_gv_genrand_uint32_range(__numV);  // generate a random node
				}
				selected[node] = true;
				pol_nodes.push_back(node);
			}
			__Activated.assign(__numV, 0);
			for(auto pol_node: pol_nodes){
				__Activated[pol_node] = 1;
			}
			__numV_left=__numV-pol_node_num_for_accuracy_verification;
			__eta_left=eta_for_verification;			
			decimal = 1.0 * (__numV_left) / (__eta_left);
			root_num=floor(decimal);
			residual = decimal - root_num;
			cout<<__LINE__<<endl;
			MC_estimation= spread_simulation(seeds, MC_round, pol_node_num_for_accuracy_verification);
			double MC_LB_error=(__eta_left-seed_num_for_accuracy_verification)*sqrt( log(__numV_left)/(2*MC_round) );
			inf_LB=MC_estimation-MC_LB_error;
			double MC_UB_error=MC_LB_error/sqrt(log(2));
			inf_UB=MC_estimation+MC_UB_error;
			double theta=max(2*__numV_left*log(__numV_left)/(inf_LB*eps_for_verification*eps_for_verification), (2+2*eps_for_verification/3)*__numV_left*log(__numV_left)/(inf_LB*eps_for_verification*eps_for_verification));

			// estimate with fresh mRR-sets
			__Activated.assign(__numV, 0);
			for(auto pol_node: pol_nodes){
				__Activated[pol_node] = 1;
			}
			cout<<__LINE__<<": building fresh mRR-sets on residual graph for theta = "<<theta<<endl;
			RR.build_n_mRRsets_fresh_vec(theta);
			double RR_coverage=0;
			vector<bool> RR_mark(theta, false);
			for(int seed:seeds)
			{
				for(int rid:RR._FRsets[seed])
				{
					if(RR_mark[rid]==false)
					{
						RR_mark[rid]=true;
						RR_coverage++;
					}
				}
			}
			RR_estimation = (RR_coverage/theta)*__eta_left;
			RR.refresh_RRsets();

			// estimate with updated mRR-sets
			theta=ceil(theta);
			RR._num_mRRsets = theta;
			RR.vv_virtual_roots.resize(theta);
			RR._mRRsets.resize(theta);
			RR.vec_mRR_layer.resize(theta);
			RR.vv_polluted_nodes.resize(theta);
			RR.vecRoot_num.resize(theta);
			__Activated.assign(__numV, 0);  // reset the __Activated states
			vvint vec_del_nodes(theta);
			__numV_left=__numV;
			__eta_left=eta_for_verification+pol_node_num_for_accuracy_verification;
			decimal = 1.0 * (__numV_left) / (__eta_left);
			root_num=floor(decimal);
			residual = decimal - root_num;
			// generate fresh mRR-sets
			cout<<__LINE__<<": building tree mRR-sets on original graph"<<endl;
			for (auto i = 0; i < theta; i++)  // if the number of previous mRR-sets is not enough, new mRR-sets will be generated
			{
				RR.build_one_mRRset_tree(i, root_num, residual);
			}
			cout<<"Reach here: "<<__LINE__<<endl;
			// activate the nodes
			for(int pol_node:pol_nodes)
			{
				__Activated[pol_node] = 1;  // activate the polluted nodes
				for(int rid:RR._FRsets[pol_node])
				{
					vec_del_nodes[rid].push_back(pol_node);
				}
			}
			__numV_left=__numV-pol_node_num_for_accuracy_verification;
			__eta_left=eta_for_verification;
			decimal = 1.0 * (__numV_left) / (__eta_left);
			root_num=floor(decimal);
			residual = decimal - root_num;
			// update the mRR-sets
			cout<<__LINE__<<": updating the tree mRR-sets"<<endl;
			for (auto i = 0; i < theta; i++)  // if the number of previous mRR-sets is not enough, new mRR-sets will be generated
			{
				if(!vec_del_nodes[i].empty())
				{
					RR.mRR_update(i, vec_del_nodes[i]);
				}
				int root_diff=root_num-RR._mRRsets[i].size()+(dsfmt_gv_genrand_open_close() <= residual);
				if(root_diff>0)
				{
					RR.add_root(i, root_diff);
				}
				else if(root_diff<0)
				{
					RR.delete_root(i, -root_diff);
				}
			}
			RR_coverage=0;
			RR_mark.assign(theta, false);
			for(int seed:seeds)
			{
				for(auto rid:RR._FRsets[seed])
				{
					if(RR_mark[rid]==false)
					{
						RR_mark[rid]=true;
						RR_coverage++;
					}
				}
			}
			update_estimation = RR_coverage/theta*__eta_left;
			RR.refresh_RRsets();
			double approx=1.0-std::exp(-1.0);
			cout << "Accuracy verification results: inf_UB = " << inf_UB << " (1+eps)*inf_UB = "<<(1.0+eps_for_verification)*inf_UB<<", MC_estimation = "<<MC_estimation<<", fresh mRR_estimation = "<<RR_estimation<<", update_estimation = "<<update_estimation<<", (1-1/e)(1-eps)*inf_LB = "<<approx*(1.0-eps_for_verification)*inf_LB<<", inf_LB = "<<inf_LB << endl;
			pre_pol_node_num=pol_node_num_for_accuracy_verification;
		}
		return make_tuple(inf_UB, MC_estimation,RR_estimation, update_estimation, inf_LB);
	}

	double spread_simulation(vint &vec_seed, int MC_round=1000, int pol_node_num=0)
	{
		auto start = std::chrono::high_resolution_clock::now();
		uint32_t nodeId, num_truncate=0;
		double spread=0.0;
		vint vec_visitNode; vec_visitNode.reserve(__numV);  
		uint32_t* visited = (uint32_t *)calloc(__numV, sizeof(uint32_t)); 
		vector<double> vecThr(__numV);  // Threshold of LT
		vector<double> vecActivateWeight(__numV, 0.0);  // total weight of active neighbors in LT
		for (auto seed : vec_seed)
		{
			__Activated[seed] = 1;
		}
		vec_visitNode= vec_seed;
		for (uint32_t i = 0; i < MC_round; i++)
		{
			// if(std::fmod(i+1, 1000) == 0)
			// {
			// 	cout << "Round " << i+1 << " is running..." << endl;
			// }
			vec_visitNode.resize(seed_num_for_accuracy_verification);
			int curIdx=0, numVisit=seed_num_for_accuracy_verification;
			if (model == "IC")
			{
				while (curIdx<numVisit)
				{
					nodeId = vec_visitNode[curIdx++];
					for (const auto& nbr : O_graph[nodeId])
					{
						if (__Activated[nbr]) continue;
						if (dsfmt_gv_genrand_open_close() <= Inv_inDeg[nbr])
						{
							__Activated[nbr] = 1;
							vec_visitNode.push_back(nbr);
							numVisit++;
						}
					}
				}
			}
			else 
			{
				while (curIdx<numVisit)
				{
					nodeId = vec_visitNode[curIdx++];
					for (const auto& nbr : O_graph[nodeId])
					{
						if (__Activated[nbr]) continue;
						if (visited[nbr] < i + 1)
						{
							visited[nbr] = i + 1;
							vecThr[nbr] = dsfmt_gv_genrand_open_close();
							vecActivateWeight[nbr] = 0.0;
						}
						vecActivateWeight[nbr] += Inv_inDeg[nbr];
						if (vecActivateWeight[nbr] >= vecThr[nbr])
						{
							__Activated[nbr] = 1;
							vec_visitNode.push_back(nbr);
							numVisit++;
						}
					}
				}
			}
			int active = count(__Activated.begin(), __Activated.end(), 1)-pol_node_num;
			if(active >__eta_left)
			{
				active=__eta_left;
				num_truncate++;
			}
			spread += active;
			for(int i=seed_num_for_accuracy_verification; i<vec_visitNode.size(); i++)
			{
				__Activated[vec_visitNode[i]] = 0;  // reset the activated nodes
			}
		}
		free(visited);
		auto end = std::chrono::high_resolution_clock::now();
		std::chrono::duration <double> elapsed = end - start;
		cout<<"The time for simulation is "<<elapsed.count() << " s"<<endl;
		cout<<"The percentage of truncated rounds is "<<1.0*num_truncate/MC_round<<endl;
		return 1.0*spread/MC_round;
	}

	void output_mRR(int mRRid)
	{
		string res="../results/backup.txt";
		std::ofstream result_bk;
		result_bk.open(res, ios::app);
        assert(!result_bk.fail());
		result_bk.close();
	}

	void release_memory()
	{
		RR.release_memory();
	}

}; // cls

using TAlg = Algorithm;
using PAlg = std::shared_ptr<TAlg>;