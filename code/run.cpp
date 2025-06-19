#include "../dSFMT/dSFMT.h"
#include "graph.h"
#include <iostream>
#include <vector>
#include "CommonStruc.h"
#include <cstring>
#include "Timer.h"
#include "Memory.h"
#include "MemoryUsage.h"
#include "Algorithm.h"
// #include "Argument.h"
#include <queue>
using namespace std;


int main(int argn, char **argv)
{
    
    Argument arg;  // claimed in Argument.h
    R_graph.clear(), O_graph.clear();  // global variables
    dsfmt_gv_init_gen_rand(static_cast<uint32_t>(time(nullptr)));  // the type must be uint32_t, to be accord with the function definition
    arg.arg_update(argn, argv);
    // std::fstream result_bk("../results/backup", ios::app);
    // assert(!result_bk.fail());
    double avg_cost = 0.0;
    double avg_time = 0.0;
    while(arg.times < arg.run_times)
    {
        arg.Initialization();
        auto k = arg.dataset_No;
        arg.load_cost_graph(k);
        #ifdef DEBUG
        if(arg.graph_sort_check(O_graph)==false)
        {
            cout << "O_Graph sort check failed!" << endl;
            exit(1);
        }
        if(arg.graph_sort_check(R_graph) == false)
        {
            cout << "R_Graph sort check failed!" << endl;
            exit(1);
        }
        #endif
        // if(k==0) arg.eta_0=0.95;
        // if(k==1) arg.eta_0=0.1;
        if (arg.eta_0 > 1) arg.eta_0 = arg.eta_0/(1.0*arg.numV);
        double total_cost=0.0;
        __eta_left = (arg.eta_0 ) * arg.numV;
        root_num=ceil(1.0/arg.eta_0);  // initial root number
        cout << " Running MINE_Alg at eta = " << arg.eta_0 << ", dataset = " << arg.dataset[k] << ", # node = " << arg.numV << ", eta = " << __eta_left <<", batch = "<<arg.batch<<", eps = "<<arg.eps<<", model = "<<arg.model <<", Rnd_cost = "<<arg.Rnd_cost<< ", real_time_pw = "<<arg.real_time_pw << endl;
        // result_bk<< "Running MINE_Alg at eta = " << arg.eta_0 + i * 0.01 << ", dataset = " << arg.dataset[k] << ", # node = " << arg.numV << ", eta = " << __eta_left <<", batch = "<<arg.batch<<", eps = "<<arg.eps<<", model = "<<arg.model <<", Rnd_cost = "<<arg.Rnd_cost<<", time = "<<arg.time<< ", real_time_pw = "<<arg.real_time_pw <<", q_ratio = "<<arg.q_ratio << endl;
        TAlg Alg(arg);

        vector<int> seeds;
        auto RR_info = Alg.AdaptiveSelect();
        seeds = seed_set;
        auto memory = getProcMemory();
        for (auto node : seeds)
        {
            total_cost += cost[node];
        }
        avg_cost += total_cost;
        avg_time += get<4>(RR_info);
        string results;
        results = "(" + arg.dataset[k] + ", eta = " + to_string(arg.eta_0) + ", Alg = " + "MINE, cost = " + to_string(total_cost) + ", prob = " + to_string(1.0) + ", time = " + to_string(get<4>(RR_info))  + ", memory = " + to_string(memory) + ", total_mRR = " + to_string(get<0>(RR_info)) +", mRR_update = " + to_string(get<1>(RR_info)) + +", mRR_add_back = " + to_string(get<2>(RR_info)) +", mRR_delete = " + to_string(get<3>(RR_info)) +  ")";
        // result_bk << results << endl;
        cout << results << endl;
        results = arg.dataset[k] + " " + to_string(arg.eta_0) + " MINE" + " " + to_string(total_cost) + " " + to_string(1.0) + " " + to_string(get<4>(RR_info)) + " " + to_string(get<5>(RR_info)) + " " + to_string(get<0>(RR_info)) + " " + to_string(get<1>(RR_info)) + " " + to_string(get<2>(RR_info)) + " " + to_string(get<3>(RR_info));
        cout << results << endl;
        if(arg.seed_out)
        {
            arg.seed_record(seeds,k,arg.times);
        }
        Alg.release_memory();
        arg.times++;
    }
    cout << "Average cost: " << avg_cost / arg.run_times << endl;
    cout << "Average time: " << avg_time / arg.run_times << endl;
    return 0;
}