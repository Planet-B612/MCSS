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
	double total_theta=0;
	string  _cascadeModel;
	vector<tuple<int, int, double, int>> ratio_plain, ratio_UB;
	// vector<double> vec_UB; // store the upper bound of coverage for each node
	// vector<double> vec_LB;
	vint vec_deg;
	string model="IC";
	
	int counter=0;  // record the number of nodes being affected in total
	Nodelist seed_batch;
	int seed;
	int batch_size=1;
	int deg=0;
	double approx=1.0;
	int theta=0;
	int ending_rnd=3000000;  // not needed, if q_ratio is properly set. sample: 4, facebook: 80, dblp: 3000. Let it be a large value, so that it will never enter the ending round to generate fresh mRRsets.
	const int root_num_bound=250; // try to delete unnecessary mRR-sets when the number of roots in an mRR exceeds this value. sample:0, facebook: 25, dblp: 250
	const int window_size=5;  // sample:2, facebook: 3, dblp: 5
	bool in_ending_rnd=false;
	bool delete_extra_mRR=false;
	vector<int> vec_mRR_num;
	int window_beg=0;
	ulint max_size_within_window=0;
	// Nodelist num_deg_incremental;
	vector<bool> RR_Mark;
	double a_1=0.0, a_2=0.0;
	int __dataset_No = 0;
	// double total_build_seedset_time = 0;
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
		vec_mRR_num.resize(1e4,0);
		__dataset_No = arg.dataset_No;
		__left_num = arg.left_num;
		__over_pnodes = arg.over_pnodes;
		// vec_UB= vector<double>(__numV, 0.0);
		// vec_LB= vector<double>(__numV, 0.0);
		vec_deg= vector<int>(__numV, 0);
		model=arg.model;
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
			delta=eps/(100.0*(1-1/2.71828)*(1-eps)*__eta_left);
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
			delta=1.0/__numV_left;
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
			if (batch_size > 1) select=build_seedset(theta, ratio);
			else  select=build_max_single_seed(theta, ratio);
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
		if (batch_size > 1) build_seedset(theta);
		else  build_max_single_seed(theta);
		(seed_set).insert((seed_set).end(),seed_batch.begin(),seed_batch.end());
		total_theta+=theta;
		return;
	}

	tuple<int,int,int,int,double,double> AdaptiveSelect()
	{
		approx=1.0-power((1-1.0/batch_size),batch_size);
		std::ofstream result;
		string file_name = "../results/round/round_" + std::to_string(__dataset_No) + "_" + std::to_string(static_cast<int>(__eta*__numV))+ "_" + std::to_string(batch_size) + "_" + std::to_string(eps) + ".txt";
		result.open(file_name, ios::app);
		assert(!result.fail());
		auto single_start = std::chrono::high_resolution_clock::now();
		while((__eta_left)>0)
		{
			auto start = std::chrono::high_resolution_clock::now();
			decimal = 1.0 * (__numV_left) / (__eta_left);
			root_num=floor(decimal);
			residual = decimal - root_num;  // in (0,1)
			
			if(__eta_left<eta_left_threshold)
			{
				in_ending_rnd=true;
				// cout<<"Entering the ending round"<<endl;
				RR.refresh_RRsets();
					// cout<<"Refreshed the RR-sets complete"<<endl;
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
			if(in_ending_rnd)
			{
				counter=RR.realization_fresh_vec(seed_batch);
			}
			else
			{
				counter=RR.realization(seed_batch);
			}
			auto now = std::chrono::high_resolution_clock::now();
			std::chrono::duration<double> elapsed = now - start;

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
			result<<(round_num)<<", \t"<<counter<<", \t"<<(__eta_left)<<", \t"<<1.0*(__numV_left)/(__eta_left)<<", \t theta = "<<theta<<"; \t update: "<< RR.num_update <<"\t add root: "<< RR.num_add_root<<" \t delete root: "<<RR.num_delete_root<<", \t"<< disp_mem_usage()<<" MB, \t"<<elapsed.count() << " 秒"<<endl;  // the round that is currently running
			RR.num_update=0;
			RR.num_add_root=0;
			RR.num_delete_root=0;
			RR.num_update_this_round=0;
		}
		auto single_end = std::chrono::high_resolution_clock::now();
		std::chrono::duration <double> single_elapsed = single_end - single_start;
		result.close();
		// cout << "build seed set traversal time " << total_build_seedset_time <<endl;
		// cout << "mRR traversal time " <<RR.mRR_traversal_time <<endl;
		cout << "Single time " << single_elapsed.count() << " s" << endl;
		cout << "Single spread " << (-(__eta_left))+ __eta*__numV <<endl;
		return make_tuple(total_theta, RR.num_update, RR.num_add_root, RR.num_delete_root, single_elapsed.count(),(-(__eta_left))+ __eta*__numV);
	}

	void accuracy_verification()
	{
		std::sort(pol_node_num_for_accuracy_verification.begin(), pol_node_num_for_accuracy_verification.end());
		int max_pol_node_num = pol_node_num_for_accuracy_verification.back();
		if(max_pol_node_num+__eta_left>__numV || pol_node_num_for_accuracy_verification[0]<0)
		{
			cout<<"The number of polluted nodes for accuracy verification is not valid, please check the input."<<endl;
			return;
		}
		if(__eta_left<seed_num_for_accuracy_verification)
		{
			cout<<"The number of seeds for accuracy verification is larger than the number of nodes left, please check the input."<<endl;
			return;
		}
		vector<bool> selected(__numV, false);
		vint pol_nodes, seeds; pol_nodes.reserve(max_pol_node_num); seeds.reserve(seed_num_for_accuracy_verification);
		for(int i=0;i<max_pol_node_num;i++)
		{
			int node = dsfmt_gv_genrand_uint32_range(__numV);  // generate a random node
			while(selected[node])
			{
				node = dsfmt_gv_genrand_uint32_range(__numV);  // generate a random node
			}
			selected[node] = true;
			pol_nodes.push_back(node);
		}
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

		for(int pol_node_num:pol_node_num_for_accuracy_verification)
		{
			for(int i=0;i<pol_node_num;i++)
			{
				__Activated[pol_nodes[i]] = true;  // activate the polluted nodes
			}
		}

	}

	double spread_simulation(vint &vec_seed, int MC_round=1000)
	{
		auto start = std::chrono::high_resolution_clock::now();
		uint32_t nodeId, num_truncate=0;
		queue<uint32_t> Que;
		double spread=0.0;
		uint32_t* visited = (uint32_t *)calloc(__numV, sizeof(uint32_t)); 
		vector<double> vecThr(__numV);  // Threshold of LT
		vector<double> vecActivateWeight(__numV, 0.0);  // total weight of active neighbors in LT
		for (uint32_t i = 0; i < MC_round; i++)
		{
			__Activated.assign(__numV, false);  // reset the activated nodes
			for (auto seed : vec_seed)
			{
				Que.push(seed);
				__Activated[seed] = true;
			}

			if (model == "IC")
			{
				while (!Que.empty())
				{
					nodeId = Que.front();
					Que.pop();
					for (auto& nbr : O_graph[nodeId])
					{
						if (__Activated[nbr]) continue;
						if (dsfmt_gv_genrand_open_close() <= Inv_inDeg[nbr])
						{
							__Activated[nbr] = true;
							Que.push(nbr);
						}
					}
				}
			}
			else 
			{
				while (!Que.empty())
				{
					nodeId = Que.front();
					Que.pop();
					for (auto& nbr : O_graph[nodeId])
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
							__Activated[nbr] = true;
							Que.push(nbr);
						}
					}
				}
			}
			int active = count(__Activated.begin(), __Activated.end(), true);
			if(active >__eta_left)
			{
				active=__eta_left;
				num_truncate++;
			}
			spread += active;
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