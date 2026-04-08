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
    double avg_cost = 0.0, avg_time = 0.0, avg_build_time = 0.0, avg_revise_time = 0.0, avg_memory = 0.0;
    while(arg.times < arg.run_times && do_verify==false)
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
        cout << " Running MINE_Alg at eta = " << __eta_left << ", dataset = " << arg.dataset[k] << ", # node = " << arg.numV << ", batch = "<<arg.batch<<", eps = "<<arg.eps<<", model = "<<arg.model <<", Rnd_cost = "<<arg.Rnd_cost<< ", real_time_pw = "<<arg.real_time_pw <<", delta_amp = "<<arg.delta_amp<<", adapt_IM = "<<adapt_IM << endl;
        // result_bk<< "Running MINE_Alg at eta = " << arg.eta_0 + i * 0.01 << ", dataset = " << arg.dataset[k] << ", # node = " << arg.numV << ", eta = " << __eta_left <<", batch = "<<arg.batch<<", eps = "<<arg.eps<<", model = "<<arg.model <<", Rnd_cost = "<<arg.Rnd_cost<<", time = "<<arg.time<< ", real_time_pw = "<<arg.real_time_pw <<", q_ratio = "<<arg.q_ratio << endl;
        TAlg Alg(arg);

        vector<int> seeds;

        if(mRR_time_test)
        {
            root_num=10;
            residual=0.0;
            cout<<"tesing time of building mRRsets..."<<endl;
            vector<double> time_mine_mRR, time_mRR_fresh;
            for(int num=2e6;num<=1e7;num+=2e6)
            {
                cout<<"generating "<<num<<" mRRsets..."<<endl;
                auto mine_start = std::chrono::high_resolution_clock::now();
                // for (auto i = 0; i < num; i++)
                // {
                //     cout<<"building mRR-set "<<i+1<<endl;
                //     Alg.RR.build_one_mRRset_tree(i, 100, 0);
                // }
                Alg.RR.build_n_mRRsets_tree(num, 0);
                auto mine_end = std::chrono::high_resolution_clock::now();
                double build_time_mine = std::chrono::duration<double>(mine_end - mine_start).count();
                time_mine_mRR.push_back(build_time_mine);
                Alg.RR.refresh_FRmRRsets(0);
                cout<<"my mRRsets built, now building fresh mRRsets..."<<endl;

                auto fresh_start = std::chrono::high_resolution_clock::now();
                for (auto i = 0; i < num; i++)
                {
                    Alg.RR.build_one_mRRset_fresh_vec(i, root_num, 0);
                }
                auto fresh_end = std::chrono::high_resolution_clock::now();
                double build_time_fresh = std::chrono::duration<double>(fresh_end - fresh_start).count();
                time_mRR_fresh.push_back(build_time_fresh);
                Alg.RR.refresh_FRmRRsets(0);
            }
            for(int i=0;i<time_mine_mRR.size();i++)
            {
                cout << "num_mRRsets: " << (i+1)*2e6 << ", time_mine_mRR: " << time_mine_mRR[i] << " s, time_fresh_mRR: " << time_mRR_fresh[i] << " s" << endl;
            }
            exit(0);
        }

        tuple<int, int, int, int, double, double> RR_info;
        double adaIM_time=0.0;
        if(adapt_IM)
        {
            auto adaIM_beg = std::chrono::high_resolution_clock::now();
            Alg.AdaptiveIM();
            auto adaIM_end = std::chrono::high_resolution_clock::now();
            adaIM_time=std::chrono::duration<double>(adaIM_end-adaIM_beg).count();
        }
        else
        {
            RR_info= Alg.AdaptiveSelect();
        }
        seeds = seed_set;
        cout<<"The number of seeds is "<<seeds.size()<<endl;
        auto memory = getProcMemory();
        for (auto node : seeds)
        {
            total_cost += cost[node];
        }
        avg_cost += total_cost;
        if(adapt_IM)
        {
            avg_time += adaIM_time;
        }
        else
        {
            avg_time += get<4>(RR_info);
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
        arg.times++;
    }
    cout << "Average cost: " << avg_cost / arg.run_times << endl;
    cout << "Average time: " << avg_time / arg.run_times << endl;
    cout << "Average memory: " << avg_memory / arg.run_times << " MB" << endl;
    cout << "Average build mRRset time: " << avg_build_time / arg.run_times << " s" << endl;
    cout << "Average revise mRRset time: " << avg_revise_time / arg.run_times << " s" << endl;

    if(do_verify)
    {
        arg.Initialization();
        auto k = arg.dataset_No;
        arg.load_cost_graph(k);
        TAlg Alg(arg);
        Alg.accuracy_verification();
    }
    return 0;
}