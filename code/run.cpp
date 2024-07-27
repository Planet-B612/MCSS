#include "head.h"
#include "./dSFMT/dSFMT.h"
#include "CommonStruc.h"
#include "CommonFunc.h"
#include "graph.h"
#include <iostream>
#include "Memory.h"
#include "Algorithm.h"
#include <vector>
#include <algorithm>
#include "Timer.h"
#include <string>
#include <ctime>

using namespace std;

int main(int argn, char **argv)
{
	usint dataset_No=3;
	double prob=0.6; 
	double eta_0=0.08;
	usint setting=3;
	CascadeModel model=IC;
	double eps_MC=0.07; 		
	double delta_MC=0.01; 
	double tau=0.01;  // 0.01 by default.
	double sigma=0.05;
	double epsilon_Inf=0.1;   // eps of influence estimation
	double veri_epsilon_Inf=0.01; // eps of influence estimation during inf validation
	double delta_Inf=0.01;
	bool Rnd_cost=true;
	double alpha=0.02; 
	uint32_t simRnd=10000;
	uint E_PCG=0; // 0: PCG, 1: ECG + ratio output, 2: ECG
	uint32_t dist=10;  // 10 by default
	double xi=0.002;
	uint format_graph=0;

	bool isReverse=true;

	for (int i = 0; i < argn; i++)
    {
		if (argv[i] == string("-dataset_No"))
			dataset_No = stoi(argv[i + 1]);
		if (argv[i] == string("-setting"))
			setting = stoi(argv[i + 1]);
		if (argv[i] == string("-prob"))
			prob = atof(argv[i + 1]);
		if (argv[i] == string("-eta_0"))
			eta_0 = atof(argv[i + 1]);		
		if (argv[i] == string("-model"))
			model = argv[i + 1]=="LT" ? LT : IC;
		if (argv[i] == string("-tau"))
			tau = atof(argv[i + 1]);
		if (argv[i] == string("-sigma"))
			sigma = atof(argv[i + 1]);
		if (argv[i] == string("-eps_MC"))
			eps_MC = atof(argv[i + 1]);
		if (argv[i] == string("-delta_MC"))
			delta_MC = atof(argv[i + 1]);
		if (argv[i] == string("-alpha"))
			alpha = atof(argv[i + 1]);
		if (argv[i] == string("-epsilon_Inf"))
			epsilon_Inf = atof(argv[i + 1]);
		if (argv[i] == string("-simRnd"))
			simRnd = atoi(argv[i + 1]);
		if (argv[i] == string("-dist"))
			dist = atoi(argv[i + 1]);
		if (argv[i] == string("-E_PCG"))
			E_PCG = atoi(argv[i + 1]);
		if (argv[i] == string("-format_graph"))
			format_graph = atoi(argv[i + 1]);
    }

	double delta_1=delta_MC/3.0;
	double lambda=eps_MC; // probability estimation error of RR
	double delta_2=delta_MC-delta_1; 
	uint8_t mode=0;  // 0: no estimation, 1: MC, 2:MRR
	double gamma=sigma;

	vector<string> dataset={"facebook", "livejournal", "pokec", "dblp", "friendster"};
	vector<string> alg_arr={"KDD", "MINE-ECG", "GRR"};
	if(E_PCG==0) alg_arr[1]="MINE-PCG";
	vector<usint> algs={1};
	vector<usint> data={3};
	vector<double> RR_ratio={1.0};
	dsfmt_gv_init_gen_rand(static_cast<uint32_t>(time(nullptr)));
	cout<<"# The parameters are: Aiming prob="<<prob<<", eps_MC="<<eps_MC<<", delta_MC="<<delta_MC<<", tau="<<tau<<", delta_1="<<delta_1<<", delta_2="<<delta_2<<", lambda="<<lambda<<", gamma="<<gamma<<", sigma="<<sigma<<", alpha="<<alpha<<", dist="<<dist<<endl;
	cout<<"The num of simulations in each MC is "<<ceil( log(2/delta_MC)/(2*eps_MC*eps_MC) )<<". The num of realizations in MRRcol is "<<ceil( log(1/delta_2)/(2*lambda*lambda) )<<". The num of roots for each MFF is "<<ceil( log(1/delta_1)/(2*tau*tau) )<<endl;
	string result_dir="backup.txt";
	std::fstream result_bk(result_dir, ios::app);
	ASSERT(!result_bk.fail());
	result_bk<<"# The parameters are: Aiming prob="<<prob<<", eps_MC="<<eps_MC<<", delta_MC="<<delta_MC<<", tau="<<tau<<", delta_1="<<delta_1<<", delta_2="<<delta_2<<", lambda="<<lambda<<", gamma="<<gamma<<", sigma="<<sigma<<", alpha="<<alpha<<", dist="<<dist<<endl;
	for(usint k: data)
	{
		dataset_No=k;
		string dataset_dir = "./graphInfo/"+dataset[dataset_No];
		if(format_graph==1) GraphBase::format_graph(dataset_dir, !isReverse);
		if(format_graph==2)	GraphBase::format_graph(dataset_dir, isReverse);
		if(format_graph!=0) return 0;
		Graph R_graph = GraphBase::load_graph(dataset_dir, isReverse);
		uint32_t numV=R_graph.size();
		for(usint i=0;i<1;i++)
		{
			eta_0=eta_0+0.02;
			for(double ratio:RR_ratio)
			{
				for(usint alg:algs)
				{
					setting=alg;
					double eta=eta_0*numV;
					cout<<"Running alg: "<<alg_arr[setting]<<", at eta = "<<eta_0<<", dataset = "<<dataset[dataset_No]<<", # node = "<<numV<<", eta = "<<eta<<endl;
					Graph F_graph={};
					vector<double> cost_0(numV,1.0);
					double total_cost=0.0;
					uint8_t cost_aware=0;
					if(setting==0)
					{
						mode=1;
						cost_aware=0;
					}
					if(setting==1)
					{
						mode=2;
						cost_aware=1;
					}
					if(setting==2)
					{
						mode=3;
						cost_aware=1;
					}

					// load_cost
					vector<double> cost_ca(numV);
					string cost_file;
					if(Rnd_cost)
					{
						cost_file="./graphInfo/"+dataset[dataset_No]+"_cost_Rand.txt";
					}
					else
					{
						cost_file="./graphInfo/"+dataset[dataset_No]+"_cost_001DEG.txt";
					}
					std::ifstream inFile;
					inFile.open(cost_file);
					if(!inFile)
					{
						cout<<"cannot open the cost file at "<<cost_file<<endl;
						exit(1);
					}
					double nodeCost;
					double c_min=0.0;
					inFile.seekg(0, std::ios_base::beg);
					for(size_t i=0;i<numV;i++)
					{
						inFile>>nodeCost;
						cost_ca[i]=nodeCost;
						c_min=min(c_min, nodeCost);
					}
					//cout<<"The min cost is: "<<c_min<<endl;
					inFile.close();
					F_graph = GraphBase::load_graph(dataset_dir, !isReverse);  
					TAlg Alg(R_graph);
					if(cost_aware==0)	Alg.set_parameters(model, cost_0, mode, prob, F_graph, eta, eps_MC, delta_MC, tau, delta_1, lambda, delta_2, epsilon_Inf, delta_Inf, veri_epsilon_Inf, result_dir, E_PCG, ratio, dist, xi);
					else  				Alg.set_parameters(model, cost_ca, mode, prob, F_graph, eta, eps_MC, delta_MC, tau, delta_1, lambda, delta_2, epsilon_Inf, delta_Inf, veri_epsilon_Inf, result_dir, E_PCG, ratio, dist, xi);
					vector<uint32_t> seeds;
					seeds.clear();
					pair<double, double> ASM_res;
					Timer Alg_time("Alg_time");
					if(setting==0)
					// MinSeed
					seeds=Alg.kdd(delta_Inf, 1.0*alpha/3.0, 1.0*alpha/3.0);
					if(setting==1)
					//Mine
					{
						if(E_PCG==0)
							seeds=Alg.mine_PCG(delta_Inf, sigma, gamma); // compared with ASM, gamma seems should be arg.epsilon, to keep ratio the same.
						else
							seeds=Alg.mine_ECG(delta_Inf, sigma, gamma); // compared with ASM, gamma seems should be arg.epsilon, to keep ratio the same.
					}
					if(setting==2)
					// GRR
					seeds=Alg.mine_PCG_GRR(delta_Inf, sigma, gamma);
					auto run_time=Alg_time.get_total_time();
					auto memory=getProcMemory();
				
					for(auto node : seeds)	
					{
						total_cost+=cost_ca[node];
					}
					double actual_prob=-1.0;
					double actual_inf=0.0;

					Graph().swap(R_graph); 
					vector<double>().swap(cost_ca);  vector<double>().swap(cost_0);
					if(setting!=2) F_graph = GraphBase::load_graph(dataset_dir, !isReverse); 
					Alg.set_F_graph(F_graph);
					//cout<<"The seed num is "<<seeds.size()<<endl;
					pair<double,double> simul_Results=Alg.actual_prob_eval(seeds, simRnd);	
					actual_prob=simul_Results.first;
					actual_inf=simul_Results.second;

					result_bk<<"("<<dataset[dataset_No]<<", "<<eta_0<<", "<<alg_arr[setting]<<", "<<total_cost<<", "<<actual_prob<<", "<<run_time<<", "<<actual_inf<<", "<<memory<<")"<<endl;
					cout<<"("<<dataset[dataset_No]<<", "<<eta_0<<", "<<alg_arr[setting]<<", "<<total_cost<<", "<<actual_prob<<", "<<run_time<<", "<<actual_inf<<", "<<memory<<")"<<endl;

				}
			}
		}
		eta_0=0.08;
	}
	result_bk.close();
	return 0;
}