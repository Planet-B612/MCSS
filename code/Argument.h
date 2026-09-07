#pragma once
#include <string>
#include <iostream>
#include "CommonStruc.h"
#include "graph.h"
#include <immintrin.h>
#include <cassert>
using namespace std;

#define VAR_NAME(x) #x
// #define DEBUG

// global variables
Graph R_graph, O_graph;  
vector<float> Inv_inDeg;
int __numV_left;
double __eta_left;
int root_num;
double residual=0.0, decimal=1.0;
vector<int> __Activated;
// vector<vector<int>> activated_nodes;
vector<int> seed_set;
vector<float> cost;
int round_num=0;
bool use_UB = true;
int RR_thr=3000000;
double RR_step=100000;
int num_gen = 0;
int num_addback = 0;
int deg_amplifier=1.0; // amplify the in_deg_threshold to make diffusion easier
bool adapt_IM=false;
bool mRR_time_test = false;
bool ablation_update=false;
bool ablation_add_root=false;
const __m512i zeros512 = _mm512_setzero_epi32();
const __m512i ones512 = _mm512_set1_epi32(1);
const __m512i minus_one512 = _mm512_set1_epi32(-1);
__m512i numV512;
string mRR_output_path = "../mRR_output.txt";
string layer_output_path = "../layer_output.txt";
string FR_output_path = "../FR_output.txt";
int num_update=0, num_add_root=0;


bool do_verify = false;
int pol_node_num_for_accuracy_verification = 20000;
vint pol_node_num_for_accuracy_verification_list = {10000,20000,30000,40000,50000};
// vint pol_node_num_for_accuracy_verification_list = {50000};
int seed_num_for_accuracy_verification = 500;
int MC_round=10000;
double eps_for_verification=0.1;
int eta_for_verification=20000;

// #define debug
// #define debugroots
// vector<bool> revised;  // for debug
// double root_amplifier=1;  // for debug
// bool suppress_delete_root=false;  // for debug
// eps is also changed to be 0.2

class Argument{
public:
    float eta_0 = 5000;
    uint eta_start = 0;
    uint eta_end = 1;
    double eta_step = 0.01;
    string model = "IC";
    bool Rnd_cost = false;
    //int simRnd = 100;
    float eps = 0.9;
    double delta=0.01;
    //double delta_Inf = 0.01;  // 1/numV by default
    uint format_graph = 0; // 0: do not format graph, 1: form the forward graph, 2: form the reverse graph.
    vector<string> dataset = {"facebook", "dblp", "flickr","nethept","epinions", "youtube", "pokec", "orkut", "livejournal", "friendster","DBLP_sym","Youtube_sym","twitter","citeseer","Flickr_sym","wikitalk","wikitalkar","sample", "Twitter", "congress"};
    // vector<int> data={4};
    int dataset_No = 4;  // 17: sample, 18: Twitter, 10: DBLP_sym, 11: Youtube_sym, 4: epinions, 8: livejournal, 9: friendster
    int cur_data;//=data[0];
    string graph_path="/data/fc/graphInfo/";
    string pw_path="/data/fc/realization/";
    string result_dir = "../backup.txt";
    int run_times=1;
    int times=0;
    int batch=24;
    int linear_search_thr=0;  // recommended 50 for formal running
    bool seed_out=false;
    bool gene_ini_pw=false;
    bool real_time_pw=false;
    int numV=1;
    int eta_left_threshold = 5;
    float regen_threshold=0.15;
    int root_num_bound = 250; // try to delete unnecessary mRR-sets when the number of roots in an mRR exceeds this value. sample:0, facebook: 25, dblp: 250
    float left_num = 100.0;
    float over_pnodes = 10.0;
    std::ofstream result_bk;
    int delta_amp=1;
    
    Argument()
    {
        #ifdef xxx
        // dataset_No=5;
        real_time_pw = true;
        deg_amplifier=1.3;
        #endif
    }
    void arg_update(int argn, char** argv)
    {
        for (int i = 0; i < argn; i++)
        {
            if (argv[i] == string("-model"))
                model = argv[i + 1];
            // if (argv[i] == string("-format_graph"))
            //     format_graph = stoi(argv[i + 1]);
            // if (argv[i] == string("-simRnd"))
            //     simRnd = stoi(argv[i + 1]);
            if (argv[i] == string("-dataset_No"))
                dataset_No = stoi(argv[i + 1]);
            if (argv[i] == string("-Rand_cost"))
                Rnd_cost = stoi(argv[i + 1]);
            if (argv[i] == string("-eta_0"))
                eta_0 = stof(argv[i + 1]);
            if (argv[i] == string("-eps"))
                eps = stof(argv[i + 1]);
            if (argv[i] == string("-batch"))
                batch = stof(argv[i + 1]);
            if (argv[i] == string("-run_times"))
                run_times = stoi(argv[i + 1]);
            if (argv[i] == string("-times"))
                times = stoi(argv[i + 1]);
            if (argv[i] == string("-gene_ini_pw"))
                {
                    gene_ini_pw = stoi(argv[i + 1]);
                    cout<<"gene_ini_pw: "<<gene_ini_pw<<endl;
                }
            if (argv[i] == string("-real_time_pw"))
                real_time_pw = stoi(argv[i + 1]);
            if (argv[i] == string("-left_num"))
                left_num = stof(argv[i + 1]);
            if (argv[i] == string("-over_pnodes"))
                over_pnodes = stof(argv[i + 1]);
            if (argv[i] == string("-del_num"))
                pol_node_num_for_accuracy_verification= stoi(argv[i + 1]);
            if (argv[i] == string("-seed_num"))
                seed_num_for_accuracy_verification= stoi(argv[i + 1]);
            if (argv[i] == string("-do_verify"))
                do_verify = stoi(argv[i + 1]);
            if (argv[i] == string("-MC_round"))
                MC_round = stoi(argv[i + 1]);
            if (argv[i] == string("-eps_for_verification"))
                eps_for_verification = stof(argv[i + 1]);
            if (argv[i] == string("-eta_for_verification"))
                eta_for_verification = stof(argv[i + 1]);
            if (argv[i] == string("-regen"))
                regen_threshold = stof(argv[i + 1]);
            if (argv[i] == string("root_nbound"))
                root_num_bound = stoi(argv[i + 1]);
            if (argv[i] == string("-delta_amp"))
                delta_amp = stoi(argv[i + 1]);
            if (argv[i] == string("-adapt_IM"))
                adapt_IM = stoi(argv[i + 1]);
            if (argv[i] == string("-mRR_time_test"))
                mRR_time_test = stoi(argv[i + 1]);
            if (argv[i] == string("-ablation_update"))
                ablation_update = stoi(argv[i + 1]);
            if (argv[i] == string("-ablation_add_root"))
                ablation_add_root = stoi(argv[i + 1]);
        }      
    }   
    void Initialization()
    {
        __Activated.clear();
        // activated_nodes.clear();
        // activated_nodes.push_back({});
        seed_set.clear();
        cost.clear();
        if(batch<1)
        {
            cout << "The batch size should be larger than 0, please check the input." << endl;
            exit(1);
        }
        if(do_verify)
        {
            cout<< "The accuracy verification is enabled, the number of polluted nodes is " << pol_node_num_for_accuracy_verification << ", the number of seeds is " << seed_num_for_accuracy_verification << ", the MC_round is "<<MC_round<<", the eta_for_verification is "<<eta_for_verification<<", eps_for_verification is "<<eps_for_verification<<endl;
        }
        if(model=="IC")
        {
            cout << "The diffusion model is IC." << endl;
        }
        else
        {
            cout << "The diffusion model is LT." << endl;
        }
        // if (format_graph != 0)
        // {
        //     for(auto k:data)
        //     {
        //         string dataset_dir = graph_path + dataset[k];
        //         GraphBase::format_graph(dataset_dir, format_graph - 1);
        //     }
        //     exit(0);
        // }
        // string parameter_str = "# The parameters are: Rnd_cost=" + to_string(Rnd_cost) + ", model=" + model+", eps_Inf="+to_string(eps);
        // cout << parameter_str << endl;
        // result_bk.open(result_dir, ios::app);
        // assert(!result_bk.fail());
        // result_bk << parameter_str << endl;
        // result_bk.close();
    }
    void load_cost_graph(int k)
    {
        cur_data=k;
        string dataset_dir = graph_path + dataset[k];
        cout<<dataset_dir<<", "<<dataset[k]<<endl;
        mRR_output_path="../mRR_output_" + dataset[k] + "_" + model+ ".txt";
        layer_output_path="../layer_output_" + dataset[k] + "_" + model+ ".txt";
        FR_output_path="../FR_output_" + dataset[k] + "_" + model+ ".txt";
        GraphBase::load_graph_directly_nbr_sorted(dataset_dir, O_graph, R_graph);
        // R_graph = GraphBase::load_graph(dataset_dir, 1);
        // O_graph = GraphBase::load_graph(dataset_dir, 0);
        numV=static_cast<int>(R_graph.size());
        numV512 = _mm512_set1_epi32(numV);
        __numV_left=numV;
        round_num=0;

        Inv_inDeg.resize(numV);
        for(int i=0;i<numV;i++)
        {
            Inv_inDeg[i]=static_cast<float>(1.0/R_graph[i].size())*deg_amplifier;
        }

        cost.resize(numV);
        fill(cost.begin(), cost.end(), 1.0);
        __Activated.resize(numV);
        fill(__Activated.begin(), __Activated.end(), false);
        //eta=(eta_0 + 0.02 * k)*numV;
        string cost_file;
        if (Rnd_cost)
        {
            cost_file = graph_path + dataset[k] + "_cost_Rand.txt";
            cout << "Using the cost file at " << cost_file << endl;
        }
        else
        {
            cost_file = graph_path + dataset[k] + "_cost_001DEG.txt";
            cout << "Using the cost file at " << cost_file << endl;
        }
        std::ifstream inFile;
        inFile.open(cost_file);
        if (!inFile)
        {
            cout << "cannot open the cost file at " << cost_file << endl;
            exit(1);
        }
        inFile.seekg(0, std::ios_base::beg);
        for (int i = 0; i < numV; i++)
        {
            inFile >> cost[i];
            // assert(cost[i] >= 0 && cost[i] <= 1);
        }
        inFile.close();
        // cout<<"==================Using uniform costs!!!==============="<<endl;
        // for (size_t i = 0; i < numV; i++)
        // {
        //     cost[i]=1.0;
        // }
        if(gene_ini_pw)
        {
            for (int i = times; i < run_times; i++){
                vector<vector<int>> PO(numV);
                ofstream out_pw;
                if(model=="IC")
                {
                    out_pw.open(pw_path + dataset[dataset_No] + "_pw_ic" + to_string(i) + ".txt");
                    assert((!out_pw.fail()));
                    for(int i=0;i<(numV);i++)
                    {
                        auto nbrs=(O_graph)[i];
                        for(auto edge:nbrs)
                        {
                            if((dsfmt_gv_genrand_open_close())<Inv_inDeg[edge])
                            {
                                PO[i].push_back(edge);
                            }
                        }
                    }
                }
                else
                {
                    out_pw.open(pw_path + dataset[dataset_No] + "_pw_lt" + to_string(i) + ".txt");
                    assert((!out_pw.fail()));
                    for(int v=0;v<(numV);v++)
                    {
                        auto nbrs=(R_graph)[v];
                        auto nbrs_size=nbrs.size();
                        if(nbrs_size==0)
                        {
                            continue;
                        }
                        int index = dsfmt_gv_genrand_uint32_range(nbrs_size);
                        int u = nbrs[index];
                        PO[u].push_back(v);
                    }
                }
                for(int i=0;i<numV;i++)
                {
                    for(long unsigned int j=0;j<PO[i].size();j++)
                    {
                        out_pw<<i<<" "<<PO[i][j]<<endl;
                    }
                    // out_pw<<endl;
                }
                out_pw.close();
            }
            cout << "generate PO for"<<model<<"at "<<dataset[dataset_No]<<" Done" <<endl;
            exit(0);
        }
        // else
        // {
        //     string cost_file;
        //     if (Rnd_cost)
        //     {
        //         cost_file = graph_path + dataset[k] + "_cost_Rand.txt";
        //     }
        //     else
        //     {
        //         cost_file = graph_path + dataset[k] + "_cost_001DEG.txt";
        //     }
        //     cout << "Using the cost file at " << cost_file << endl;
        //     std::ifstream inFile;
        //     inFile.open(cost_file);
        //     if (!inFile)
        //     {
        //         cout << "cannot open the cost file at " << cost_file << endl;
        //         exit(1);
        //     }
        //     inFile.seekg(0, std::ios_base::beg);
        //     for (int i = 0; i < numV; i++)
        //     {
        //         inFile >> cost[i];
        //         // assert(cost[i] >= 0 && cost[i] <= 1);
        //     }
        //     inFile.close();
        // }
        return;
    }
    
    void seed_record(vector<int> seeds,int k,int i)
    {
        ofstream out_seeds("../results/seed/MINE_" + dataset[k] + "_" + to_string(eta_0*numV) + "_" + to_string(batch) + "_" + to_string(i) + ".txt", ios::out);
        assert((!out_seeds.fail()));
        for (auto node : seeds)
        {
            out_seeds << node << endl;
        }
        out_seeds.close();
    }

    bool graph_sort_check(const Graph &graph)
    {
        for (ulint i = 0; i < graph.size(); i++)
        {
            auto node = graph[i];
            if(node.size() <2)
            {
                continue;
            }
            for (ulint j = 0; j < node.size() - 1; j++)
            {
                if (node[j] > node[j + 1])
                {
                    cout << "Error: the graph is not sorted!" << endl;
                    return false;
                }
            }
        }
        return true;
    }
};