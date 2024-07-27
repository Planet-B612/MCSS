#pragma once
#include "RRcollection.h"
#include "CommonStruc.h"
#include "CommonFunc.h"
//#include "Argument.h"
#include <memory>
#include <algorithm>
#include "Memory.h"
using namespace std;

class Algorithm
{
	private:
		uint32_t __numV=0;
		size_t __numE=0;
		size_t num_RRsets = 0;
		size_t __veri_num_RRsets=0;
		TRRcollection RR, veriRR;
		Nodelist Seeds;
		vector<double> __cost;
		double __eta=1;
		double __eps_MC=0.1;
		double __delta_MC=0.1;
		double __tau=0.1;
		double __delta_1=0.1;
		double __lambda=__eps_MC;
		double __delta_2=0.1;
		uint8_t __mode=0; // 0: no prob estimation, 1: MC, 2: RR
		double __prob=0.1;
		double __epsilon_Inf=0.1;
		double __delta_Inf=0.01;
		double __veri_epsilon_Inf=0.05;
		Graph __F_graph;
		Graph __R_graph;
		CascadeModel __model=IC;
		string __result_dir;
		uint __E_PCG;
		double __RR_ratio=1.0;
		uint32_t __dist=1;
		double __xi=0.002;

		size_t __kappa=0;
		size_t __theta=0;
		vector<vector<bool>> __vecCover;  // record whether an FRset_real is covered by some seed
		vector<uint32_t> __vecCount;  // record how many FRset_real in a realization has been covered by the seeds
		vector<bool> __vecEta;	// record whether this realization has reached eta
		vector<double> __vec_prob_GRR;
		double __freq=0.0;

	public:

		Algorithm(Graph& graph) : RR(graph), veriRR(graph)
		{
			__numV = RR.get_nodes();
			__numE = RR.get_edges();
			__R_graph=graph;
		}
		~Algorithm()
		{ }

bool comp(tuple<uint32_t, uint32_t, double, uint32_t> a, tuple<uint32_t, uint32_t, double, uint32_t> b)
{
	return get<2>(a) > get<2>(b);
}; // sort in descending order.

vector<uint32_t> max_ratio_lazy(const double Q)
{
	vector<tuple<uint32_t, double, double, uint32_t>> ratio(__numV);
	double total_cost=0.0;
	//double maxDeg = 0.0;
	for (uint32_t i = __numV; i--;)
	{
		const uint32_t deg = RR._FRsets[i].size();  // The number of RR-sets covered by i.
		get<0>(ratio[i]) = i;
		get<1>(ratio[i])=deg;
		get<2>(ratio[i])=1.0*deg/__cost[i];
		get<3>(ratio[i])=0;
		//if (deg > maxDeg) maxDeg = deg;
	}
	make_max_heap(ratio);

	size_t total_deg = 0;  // Record the number of _RRsets covered by seed set.
	vector<bool> RR_Mark(num_RRsets, false);  // check if an RR-set is removed
	vector<bool> veriRR_state(num_RRsets, false);

	Seeds.clear();
	double expInf=0.0;
	double est_prob=0.0;
	for(uint32_t i=0;i<__numV;i++) 
	{
		uint32_t nodeId;
		while(get<3>(ratio[0])!=i)
		{
			nodeId=get<0>(ratio[0]);
			uint32_t nodeDeg=RR._FRsets[nodeId].size();
			for(uint32_t RRId: RR._FRsets[nodeId])
			{
				if(RR_Mark[RRId]==true) 	--nodeDeg;
			}
			tuple<uint32_t, double, double, uint32_t> updated_node=make_tuple(nodeId, nodeDeg, 1.0*nodeDeg/__cost[nodeId],i);
			max_heap_replace_max_value(ratio, updated_node);
		}
		nodeId=get<0>(ratio[0]);
		Seeds.push_back( nodeId );
		total_cost=total_cost+__cost[nodeId];

		if(__mode==1)	// 1: MC
		{
			total_deg=total_deg+get<1>(ratio[0]);  // No prob estimation, only evaluate based on expected influence.
			expInf=1.0*__numV*total_deg/num_RRsets;
			if(std::fmod(i, __dist)==0) est_prob=prob_est_MC(Seeds);  
			if(est_prob>=__prob+__eps_MC) 
			{
				//cout<<"The total cost is "<<total_cost<<endl;
				cout<<"The expected inf is: "<<expInf<<endl;
				return Seeds;
			}
		}
		else if (__mode==2)  // 2: RR
		{
			total_deg=total_deg+get<1>(ratio[0]); 
			expInf=1.0*__numV*total_deg/num_RRsets;
			est_prob=prob_est_RR(nodeId); 
			if(est_prob>=__prob+__lambda) 
			{
				cout<<"The expected inf is: "<<expInf<<", the estimated prob is: "<<est_prob<<endl;
				return Seeds;
			}
		}
		else if (__mode==3) // 3: GRR
		{
			total_deg=total_deg+get<1>(ratio[0]);  
			expInf=1.0*__numV*total_deg/num_RRsets;
			est_prob=prob_est_GRR(nodeId); 
			if(est_prob>=__prob+__lambda)  
			{
				cout<<"The expected inf is: "<<expInf<<", the estimated prob is: "<<est_prob<<endl;
				return Seeds;
			}
		}
		else
		{
			total_deg=total_deg+get<1>(ratio[0]);  
			expInf=1.0*__numV*total_deg/num_RRsets;
			if(expInf>=Q)
			{
				cout<<"No estimation is needed. The expInf is "<<expInf<<endl;
				return Seeds;
			}
		}

		for(auto rr:RR._FRsets[nodeId])
		{
			if(RR_Mark[rr]) continue;
			RR_Mark[rr]=true;
		}
		tuple<uint32_t, double, double, uint32_t> disable_node=make_tuple(nodeId, -1.0, -1.0, i);
		max_heap_replace_max_value(ratio, disable_node);
	}
	cout<<"Probably an error exists in seed selection, this line should not appear."<<endl;
	cout<<"The expInf is "<<expInf<<endl;
	return {};

}

vector<uint32_t> max_ratio_lazy_ECG_mine(const double Q)
{
	vector<tuple<uint32_t, double, double, uint32_t>> ratio(__numV);
	double total_cost=0.0;
	double RR_Deg_Thresh=1.0*(1+__epsilon_Inf)*Q*num_RRsets/__numV;
	double veriRR_Deg_Thresh=1.0*(1+__veri_epsilon_Inf)*Q*__veri_num_RRsets/__numV;
	for (uint32_t i = __numV; i--;)
	{
		const uint32_t deg = RR._FRsets[i].size();  // The number of RR-sets covered by i.
		get<0>(ratio[i]) = i;
		get<1>(ratio[i])=deg;
		get<2>(ratio[i])=1.0*deg/__cost[i];
		get<3>(ratio[i])=0;
	}
	make_max_heap(ratio);

	double total_deg = 0.0;  // Record the number of _RRsets covered by seed set.
	double veri_total_deg=0.0;
	double pre_total_deg=0.0;
	double eps_current=0.0;
	double w=log(__numV);
	double c_UB=0.0001;
	double expInf=0.0;
	double c_i=0.0001;
	fstream result_bk(__result_dir, ios::app);
	vector<bool> RR_Mark(num_RRsets, false);  // check if an edge is removed
	vector<bool> veriRR_Mark(__veri_num_RRsets, false);
	Seeds.clear();
	for(uint32_t i=0;i<__eta;i++) 
	{
		uint32_t nodeId;
		while(get<3>(ratio[0])!=i)
		{
			nodeId=get<0>(ratio[0]);
			double nodeDeg=RR._FRsets[nodeId].size();
			for(uint32_t RRId: RR._FRsets[nodeId])
			{
				if(RR_Mark[RRId]==true) 	nodeDeg-=1.0;
			}
			nodeDeg=min(1.0*nodeDeg, RR_Deg_Thresh-1.0*total_deg);
			tuple<uint32_t, double, double, uint32_t> updated_node=make_tuple(nodeId, nodeDeg, 1.0*nodeDeg/__cost[nodeId],i);
			max_heap_replace_max_value(ratio, updated_node);
		}
		nodeId=get<0>(ratio[0]);
		Seeds.push_back( nodeId );
		c_i=__cost[nodeId];
		total_cost=total_cost+c_i;
		total_deg=total_deg+get<1>(ratio[0]);
		expInf=1.0*__numV*total_deg/num_RRsets;
		pre_total_deg=veri_total_deg;
		for(auto rr: veriRR._FRsets[nodeId])
		{
			if(veriRR_Mark[rr]) continue;
			else
			{
				veriRR_Mark[rr]=true;
				veri_total_deg+=1.0;
			}
		}
		veri_total_deg=min(1.0*veri_total_deg, veriRR_Deg_Thresh);
		if((pre_total_deg>2.0*w/3.0)&&(__E_PCG==1)) // To make sure eps_current is positive.
		{
			eps_current = 1.0*pre_total_deg/(power((sqrt(pre_total_deg+2.0*w/9.0))-sqrt(w/2.0),2) - w/18.0)-1.0;
			c_UB=max( c_UB, 1.0*c_i*(__eta*__veri_num_RRsets/__numV-(1.0+eps_current)*pre_total_deg)/((1.0+eps_current)*(total_deg-pre_total_deg)+2.0*eps_current*pre_total_deg/c_i) );
		}

		if(veri_total_deg>=veriRR_Deg_Thresh-0.01)
		{
			cout<<"expInf = "<<expInf<<", total cost = "<<total_cost<<", c_UB = "<<c_UB<<". approx_ratio= "<<1.0*total_cost/c_UB<<endl;
			if(__E_PCG==1) result_bk<<"numV = "<<__numV<<", eta = "<<1.0*__eta/__numV<<", expInf = "<<expInf<<", veri_total_deg = "<<veri_total_deg<<", veriRR_Deg_Thresh = "<<veriRR_Deg_Thresh<<". approx_ratio = "<<1.0*total_cost/c_UB<<endl;
			result_bk.close();
			return Seeds;
		}
		for(auto rr:RR._FRsets[nodeId])
		{
			if(RR_Mark[rr]) continue;
			RR_Mark[rr]=true;
		}
		tuple<uint32_t, double, double, uint32_t> disable_node=make_tuple(nodeId, 0.0, -1.0, i);
		max_heap_replace_max_value(ratio, disable_node);
	}
	cout<<"Probably an error exists in seed selection, this line should not appear."<<endl;
	cout<<"veri_total_deg = "<<veri_total_deg<<", veriRR_Deg_Thresh = "<<veriRR_Deg_Thresh<<endl;
	return {};
}

double max_cover_lazy(const double Q)
{
	vector<double> coverage(__numV, 0);
	size_t maxDeg = 0;
	for (uint32_t i = __numV; i--;)
	{
		const double deg = RR._FRsets[i].size();  // The number of RR-sets covered by i.
		coverage[i] = deg;
		if (deg > maxDeg) maxDeg = deg;
	}
	RRsets degMap(maxDeg + 1);
	for (auto i = __numV; i--;)
	{
		if (coverage[i] == 0) continue;
		degMap[coverage[i]].push_back(i); 
	}

	size_t total_deg = 0;  // Record the number of _RRsets covered by seed set.
	vector<bool> edgeMark(num_RRsets, false);  

	Seeds.clear();
	double expInf=0.0;
	for (auto deg = maxDeg; deg > 0; deg--) 
	{
		auto& vecNode = degMap[deg];  
		for (auto idx = vecNode.size(); idx--;)
		{
			auto argmaxIdx = vecNode[idx];
			const auto currDeg = coverage[argmaxIdx];  
			if (deg > currDeg) 
			{
				degMap[currDeg].push_back(argmaxIdx);
				continue;
			}
			Seeds.push_back(argmaxIdx);
			total_deg += currDeg;
			expInf = 1.0 * total_deg * __numV / num_RRsets;
			double prob=prob_est_MC(Seeds);  // MC by default
			if(prob>=__prob)
			{
				return Seeds.size();
			}
			coverage[argmaxIdx] = 0;
			for (auto edgeIdx : RR._FRsets[argmaxIdx])
			{
				if (edgeMark[edgeIdx]) continue;
				edgeMark[edgeIdx] = true;
				for (auto nodeIdx : RR._RRsets[edgeIdx])
				{
					if (coverage[nodeIdx] == 0) continue; 
					coverage[nodeIdx]--; 
				}
			}
		}
		degMap.pop_back();
	}
	cout<<"???"<<endl;
	return 1.0 * __numV; // All RR sets are covered.
}

double max_cover_topk(const double Q)
{
	FRset coverage(__numV, 0);
	size_t maxDeg = 0;
	for (auto i = __numV; i--;)
	{
		const auto deg = RR._FRsets[i].size();
		coverage[i] = deg;
		if (deg > maxDeg) maxDeg = deg;
	}
	RRsets degMap(maxDeg + 1); 
	for (auto i = __numV; i--;)
	{
		degMap[coverage[i]].push_back(i);
	}
	Nodelist sortedNode(__numV); 
	Nodelist nodePosition(__numV); 
	Nodelist degreePosition(maxDeg + 2); 
	uint32_t idxSort = 0;
	size_t idxDegree = 0;
	for (auto& nodes : degMap)
	{	
		degreePosition[idxDegree + 1] = degreePosition[idxDegree] + (uint32_t)nodes.size();
		idxDegree++;
		for (auto& node : nodes) 
		{
			nodePosition[node] = idxSort;
			sortedNode[idxSort++] = node;
		}
	}
	std::vector<bool> edgeMark(num_RRsets, false);
	Seeds.clear();
	size_t total_deg = 0;
	double expInf = 0.0;
	while(expInf<Q)
	{
		const auto seed = sortedNode.back(); // The last node in sortedNode
		sortedNode.pop_back();
		total_deg += coverage[seed];
		Seeds.push_back(seed);
		coverage[seed] = 0;
		for (auto edgeIdx : RR._FRsets[seed])
		{
			if (edgeMark[edgeIdx]) continue;
			edgeMark[edgeIdx] = true;
			for (auto nodeIdx : RR._RRsets[edgeIdx])
			{
				if (coverage[nodeIdx] == 0) continue; 
				const auto currPos = nodePosition[nodeIdx]; 
				const auto currDeg = coverage[nodeIdx]; 
				const auto startPos = degreePosition[currDeg]; 
				const auto startNode = sortedNode[startPos]; 
				std::swap(sortedNode[currPos], sortedNode[startPos]);
				nodePosition[nodeIdx] = startPos;
				nodePosition[startNode] = currPos;
				degreePosition[currDeg]++;
				coverage[nodeIdx]--;
			}
		}
		expInf = 1.0 * total_deg * Q / num_RRsets;
	}
	return expInf;
}


void set_cascade_model(const CascadeModel model)
{
	__model=model;
	RR.set_cascade_model(model);
}


void set_parameters(const CascadeModel model, vector<double> cost, uint8_t mode, double prob, Graph F_graph, double eta, double eps_MC, double delta_MC, double tau, double delta_1, double lambda, double delta_2, double eps_Inf, double delta_Inf, double veri_epsilon_Inf, string result_dir, uint E_PCG, double RR_ratio, uint32_t dist, double xi)
{
	RR.set_model_mode(model, mode);
	veriRR.set_model_mode(model,mode);
	__cost=cost;
	__mode=mode;
	__prob=prob;
	__eta=eta;
	__eps_MC=eps_MC;
	delta_MC=delta_MC;
	__tau=tau;
	__delta_1=delta_1;
	__lambda=lambda;
	__delta_2=delta_2;
	__F_graph=F_graph;
	__epsilon_Inf=eps_Inf;
	__delta_Inf=delta_Inf;
	__veri_epsilon_Inf=veri_epsilon_Inf;
	__result_dir=result_dir;
	__E_PCG=E_PCG;
	__RR_ratio=RR_ratio;
	__dist=dist;
	__xi=xi;
}

void set_F_graph(const Graph F_graph)
{
	__F_graph=F_graph;
}

double effic_inf_valid_algo(const double delta = 1e-3, const double eps = 0.01)
{
	return RR.effic_inf_valid_algo(Seeds, delta, eps);
}

vector<uint32_t> mine_ECG(const double delta, const double sigma, const double gamma)
{
	Timer mine_runtime("mine_runtime_ECG");
	double ksi=__xi*__eta; 
	num_RRsets=(2.0*(1+1.0*__epsilon_Inf/3)*log(1.0/__delta_Inf)*__numV/(__epsilon_Inf*__epsilon_Inf*ksi));
	if(__E_PCG==1) num_RRsets=ceil(num_RRsets*__RR_ratio);
	cout<<"The num of RR-sets needed is "<<num_RRsets<<endl;
	RR.build_n_RRsets(num_RRsets);
	__veri_num_RRsets=2.0*(1+1.0*__veri_epsilon_Inf/3)*log(1.0/__delta_Inf)*__numV/(__veri_epsilon_Inf*__veri_epsilon_Inf*__eta);
	if(__E_PCG==1) __veri_num_RRsets=ceil(num_RRsets*__RR_ratio);
	veriRR.build_n_RRsets(__veri_num_RRsets);
	max_ratio_lazy_ECG_mine((1+__veri_epsilon_Inf)*__eta);
	auto mine_ECG_time=mine_runtime.get_total_time();
	cout<<"The runtime of MINE-ECG is "<<mine_ECG_time<<endl;
	return Seeds;
}

vector<uint32_t> mine_PCG(const double delta, const double sigma, const double gamma)
{
	double ksi=__xi*__eta;
	num_RRsets=(2.0*(1+1.0*__epsilon_Inf/3)*log(1.0/__delta_Inf)*__numV/(__epsilon_Inf*__epsilon_Inf*ksi));
	RR.build_n_RRsets(num_RRsets);
	__kappa=ceil( log(1/__delta_2)/(2*__lambda*__lambda) );  // # realizations
	__theta=ceil( log(1/__delta_1)/(2*__tau*__tau) );  // The num of roots
	__vecCover.resize(__kappa, vector<bool>(__theta, false));
	__vecEta.resize(__kappa, false);
	__vecCount.resize(__kappa, 0.0);
	//cout<<"The num of realization is "<<__kappa<<", The num of roots in each realization is "<<__theta<<endl;
	Timer MINE_MFFcol("MINE_MFFcol");
	RR.generate_MFFcol_20231029(__kappa, __theta);
	auto MINE_MFFcol_time=MINE_MFFcol.get_total_time();
	cout<<"The time for MINE MFFcol generation is "<<MINE_MFFcol_time<<endl;
	Timer MINE_Greedy("MINE_Greedy");
	max_ratio_lazy(__eta);
	auto MINE_Greedy_MC_time=MINE_Greedy.get_total_time();
	RR.release_memory();
	return Seeds;
}

vector<uint32_t> mine_PCG_GRR(const double delta, const double sigma, const double gamma)
{
	double ksi=__xi*__eta;
	num_RRsets=(2.0*(1+1.0*__epsilon_Inf/3)*log(1.0/__delta_Inf)*__numV/(__epsilon_Inf*__epsilon_Inf*ksi));
	RR.build_n_RRsets(num_RRsets);
	__kappa=ceil( log(1/__delta_2)/(2*__lambda*__lambda) );  // # realizations
	__theta=ceil( log(1/__delta_1)/(2*__tau*__tau) );  // The num of roots
	__vecCover.resize(__kappa, vector<bool>(__theta, false));
	__vecEta.resize(__kappa, false);
	__vecCount.resize(__kappa, 0.0);
	__vec_prob_GRR.resize(__kappa, 0.0);
	RR.generate_GRR(__kappa, __theta);
	cout<<"GRRs have been generated."<<endl;
	max_ratio_lazy(__eta);
	RR.release_memory();
	return Seeds;
}

vector<uint32_t> kdd(const double delta, const double sigma, const double gamma) 
{
	double Lambda=__eta;
	double W1[3]={__eta, gamma, delta/2};
	double W2[3]={__eta, sigma, delta/2};
	double ut=2*__numV*(3+W1[1])* (log(1/W1[2])+8*log(__numV)) /(3*W1[1]*W1[1]*W1[0]);
	double lt=2*__numV*log(1/W2[2])/(W2[1]*W2[1]*W2[0]);
	auto T=ceil(max(ut,lt)); 
	Timer kdd_RR_MC("KDD-RRset");
	RR.build_n_RRsets(T);
	auto kdd_RR_time=kdd_RR_MC.get_total_time();
	cout<<"The time for KDD RR-sets is "<<kdd_RR_time<<endl;
	num_RRsets=RR.get_RR_sets_size();
	Timer kdd_greedy_MC("KDD");
	max_ratio_lazy(Lambda);
	auto kdd_time=kdd_greedy_MC.get_total_time();
	return Seeds;
}

/// @brief Estimate the probability with RR-sets
/// @return the estimated probability
double prob_est_RR(const uint32_t seed)
{
	for(uint32_t i=0;i<__kappa;i++)
	{

		if(__vecEta[i])	continue;  // If this realization has already reached eta, there is no need to update it again.
		for(uint32_t FRset_real:RR._MFFcol[i][seed])
		{
			if(__vecCover[i][FRset_real]) { continue; }
			else
			{
				__vecCover[i][FRset_real]=true;
				__vecCount[i]+=1;
			}
		}
		if(1.0*__vecCount[i]*__numV/__theta >= __eta)	
		{
			__freq+=1;
			__vecEta[i]=true;
		}
	}
	return 1.0*__freq/__kappa;
}

double prob_est_GRR(const uint32_t seed)
{
	__freq=0.0;
	for(uint32_t i=0;i<__kappa;i++)
	{
		if(__vecEta[i])	continue; 
		for(uint32_t FRset_real:RR._GFFcol[i][seed])
		{
			if(__vecCover[i][FRset_real]) { continue; }
			else
			{
				__vecCover[i][FRset_real]=true;
				__vecCount[i]+=1;
			}
		}
		for (uint32_t i = 0; i < __kappa; i++)
		{
			__vec_prob_GRR[i]=1.0*__vecCount[i]/__theta;
		}
	}
	for(auto prob:__vec_prob_GRR)
	{
		__freq+=prob;
	}
	return 1.0*__freq/__kappa;
}

double prob_est_MC(const vector<uint32_t>& vecSeed)
{
	ASSERTT(!__F_graph.empty(), "The input forward graph __F_graph for prob_est_MC() is empty.");
	Timer MC_time("MC_time");
	uint32_t nodeId, currProgress = 0;
	double prob=0.0;
	queue<uint32_t> Que;
	vector<uint32_t> vecActivated;
	uint32_t simulations=ceil( log(2/__delta_MC)/(2*__eps_MC*__eps_MC) );
	ASSERTT(simulations>0, "The value of MC simulation times is invalid!");
	bool* activated = (bool *)calloc(__numV, sizeof(bool));
	uint32_t* visited = (uint32_t *)calloc(__numV, sizeof(uint32_t)); 
	vector<double> vecThr(__numV);  // Threshold of LT
	vector<double> vecActivateWeight(__numV, 0.0);  // total weight of active neighbors in LT
	for (auto seedId : vecSeed) activated[seedId] = true;
	for (uint32_t i = 0; i < simulations; i++)
	{
		for (auto seed : vecSeed)
		{
			Que.push(seed);
		}

		// BFS traversal
		if (__model == IC)
		{
			while (!Que.empty())
			{
				nodeId = Que.front();
				Que.pop();
				for (auto& nbr : __F_graph[nodeId])
				{
					if (activated[get<0>(nbr)]) continue;
					if (dsfmt_gv_genrand_open_close() <= get<1>(nbr))
					{
						activated[get<0>(nbr)] = true;
						vecActivated.push_back(get<0>(nbr)); //Records which nodes are activated by the node.
						Que.push(get<0>(nbr));
					}
				}
			}
		}
		else if (__model == LT)
		{
			while (!Que.empty())
			{
				nodeId = Que.front();
				Que.pop();
				for (auto& nbr : __F_graph[nodeId])
				{
					if (activated[get<0>(nbr)]) continue;
					if (visited[get<0>(nbr)] < i + 1)
					{
						// First time visit this node
						visited[get<0>(nbr)] = i + 1;
						vecThr[get<0>(nbr)] = dsfmt_gv_genrand_open_close();
						vecActivateWeight[get<0>(nbr)] = 0.0;
					}
					vecActivateWeight[get<0>(nbr)] += get<1>(nbr);
					if (vecActivateWeight[get<0>(nbr)] >= vecThr[get<0>(nbr)])
					{
						activated[get<0>(nbr)] = true;
						vecActivated.push_back(get<0>(nbr));
						Que.push(get<0>(nbr));
					}
				}
			}
		}
		uint32_t num_activated=count(activated, activated+__numV, true);
		if(count(activated, activated+__numV, true) >=__eta) { prob=prob+1; }
		for (auto activatedNode : vecActivated) activated[activatedNode] = false;
		vecActivated.clear();
	}
	free(activated);
	free(visited);
	return 1.0*prob/simulations;
}

pair<double, double> actual_prob_eval(const vector<uint32_t> & vecSeed, uint32_t simulations)
{
	Timer prob_eval_time("prob_eval_time");
	uint32_t nodeId;
	uint32_t exceed=0; // Record the num of simulations that exceed Q
	queue<uint32_t> Que;
	vector<uint32_t> vecActivated;
	double spread=0.0;
	bool* activated = (bool *)calloc(__numV, sizeof(bool));
	uint32_t* visited = (uint32_t *)calloc(__numV, sizeof(uint32_t)); 
	vector<double> vecThr(__numV);  // Threshold of LT
	vector<double> vecActivateWeight(__numV, 0.0);  // total weight of active neighbors in LT
	for (auto seedId : vecSeed) activated[seedId] = true;
	for (uint32_t i = 0; i < simulations; i++)
	{
		for (auto seed : vecSeed)
		{
			Que.push(seed);
		}

		// BFS traversal
		if (__model == IC)
		{
			while (!Que.empty())
			{
				nodeId = Que.front();
				Que.pop();
				for (auto& nbr : __F_graph[nodeId])
				{
					if (activated[get<0>(nbr)]) continue;
					if (dsfmt_gv_genrand_open_close() <= get<1>(nbr))
					{
						activated[get<0>(nbr)] = true;
						vecActivated.push_back(get<0>(nbr)); //Records which nodes are activated by the node.
						Que.push(get<0>(nbr));
					}
				}
			}
		}
		else if (__model == LT)
		{
			while (!Que.empty())
			{
				nodeId = Que.front();
				Que.pop();
				for (auto& nbr : __F_graph[nodeId])
				{
					if (activated[get<0>(nbr)]) continue;
					if (visited[get<0>(nbr)] < i + 1)
					{
						visited[get<0>(nbr)] = i + 1;
						vecThr[get<0>(nbr)] = dsfmt_gv_genrand_open_close();
						vecActivateWeight[get<0>(nbr)] = 0.0;
					}
					vecActivateWeight[get<0>(nbr)] += get<1>(nbr);
					if (vecActivateWeight[get<0>(nbr)] >= vecThr[get<0>(nbr)])
					{
						activated[get<0>(nbr)] = true;
						vecActivated.push_back(get<0>(nbr));
						Que.push(get<0>(nbr));
					}
				}
			}
		}
		uint32_t active = count(activated, activated+__numV, true);
		if(active >=__eta) { exceed=exceed+1; }
		spread += active;
		for (auto activatedNode : vecActivated) activated[activatedNode] = false;
		vecActivated.clear();
	}
	free(activated);
	free(visited);
	auto simu_time=prob_eval_time.get_total_time();
	cout<<"The time for simulation is "<<simu_time<<endl;
	return make_pair(1.0*exceed/simulations, 1.0*spread/simulations);
}

	/// Efficiently estimate the influence spread with sampling error epsilon within probability 1-delta
	double effic_inf_valid_algo(const vector<uint32_t>& vecSeed, const double delta = 1e-3, const double eps = 0.01)
	{
		const double c = 2.0 * (exp(1.0) - 2.0);
		const double LambdaL = 1.0 + 2.0 * c * (1.0 + eps) * log(2.0 / delta) / (eps * eps);
		size_t numHyperEdge = 0;
		size_t numCoverd = 0;
		std::vector<bool> vecBoolSeed(__numV);
		for (auto seed : vecSeed) vecBoolSeed[seed] = true;
		std::vector<bool> __vecVisitBool(__numV, false);
		std::vector<uint32_t> __vecVisitNode(__numV);

		while (numCoverd < LambdaL)
		{
			numHyperEdge++;
			size_t numVisitNode = 0, currIdx = 0;
			const auto uStart = dsfmt_gv_genrand_uint32_range(__numV);
			if (vecBoolSeed[uStart])
			{
				numCoverd++;
				continue;
			}
			__vecVisitNode[numVisitNode++] = uStart;
			__vecVisitBool[uStart] = true;
			while (currIdx < numVisitNode)
			{
				const auto expand = __vecVisitNode[currIdx++];
				if (__model == IC)
				{
					for (auto& nbr : __R_graph[expand])
					{
						const auto nbrId = get<0>(nbr);
						if (__vecVisitBool[nbrId])
							continue;
						const auto randDouble = dsfmt_gv_genrand_open_close();
						if (randDouble > get<1>(nbr))
							continue;
						if (vecBoolSeed[nbrId])
						{
							numCoverd++;
							goto postProcess;
						}
						__vecVisitNode[numVisitNode++] = nbrId;
						__vecVisitBool[nbrId] = true;
					}
				}
			}
		postProcess:
			for (auto i = 0; i < numVisitNode; i++)
				__vecVisitBool[__vecVisitNode[i]] = false;
		}
		return 1.0 * numCoverd * __numV / numHyperEdge;
	}

};//cls


using TAlg = Algorithm;
using PAlg = std::shared_ptr<TAlg>;