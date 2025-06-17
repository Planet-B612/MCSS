
#if !defined(DSFMT_MEXP)
#ifdef __GNUC__
#define DSFMT_MEXP 19937
#endif
#endif
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
    __Activated.clear();
    activated_nodes.clear();
    activated_nodes.push_back({});
    seed_set.clear();
    cost.clear();
    dsfmt_gv_init_gen_rand(static_cast<uint32_t>(time(nullptr)));  // the type must be uint32_t, to be accord with the function definition
    arg.arg_update(argn, argv);
    // arg.Initialization();
    // std::fstream result_bk("../results/backup", ios::app);
    // assert(!result_bk.fail());
    // for (auto k : arg.data)
    // {
        auto k = arg.dataset_No;
        arg.load_cost_graph(k);
        #ifndef NDEBUG
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

        // Alg.RR.test_mRR();
        // Alg.RR.vecRoot_num.resize(1,0);
        // Alg.RR._mRRsets.resize(1);
        // Alg.RR.build_one_mRRset_tree(0,2,0.0);

        vector<int> seeds;
        auto RR_info = Alg.AdaptiveSelect();
        seeds = seed_set;
        auto memory = getProcMemory();
        for (auto node : seeds)
        {
            total_cost += cost[node];
        }
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
        // }
    // }
    return 0;
}