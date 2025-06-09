#pragma once
#include <string>
#include <iostream>
#include "CommonStruc.h"
#include "graph.h"
using namespace std;

// global variables
Graph R_graph, O_graph;  
int __numV_left;
double __eta_left;
int root_num;
vector<bool> __Activated;
vector<vector<int>> activated_nodes;
vector<int> seed_set;
vector<float> cost;
int round_num=0;
bool use_UB = true;
int RR_thr=3000000;
double RR_step=100000;
int num_gen = 0;
int num_addback = 0;
int deg_amplifier=1; // amplify the in_deg_threshold to make diffusion easier

// #define debug
// #define debugroots
// vector<bool> revised;  // for debug
// double root_amplifier=1;  // for debug
// bool suppress_delete_root=false;  // for debug
// eps is also changed to be 0.2

class Argument{
public:
    float eta_0 = 0.05;
    uint eta_start = 0;
    uint eta_end = 1;
    double eta_step = 0.01;
    string model = "IC";
    bool Rnd_cost = true;
    //int simRnd = 100;
    float eps = 0.5;
    double delta=0.01;
    //double delta_Inf = 0.01;  // 1/numV by default
    uint format_graph = 0; // 0: do not format graph, 1: form the forward graph, 2: form the reverse graph.
    vector<string> dataset = {"sample", "facebook", "dblp", "flickr","nethept","epinions", "youtube", "pokec", "orkut", "livejournal", "friendster","DBLP_sym","Youtube_sym","twitter","citeseer","Flickr_sym","wikitalk","wikitalkar"};
    // vector<int> data={4};
    int dataset_No = 5;
    int cur_data;//=data[0];
    string graph_path="/data/fc/graphInfo/";
    string pw_path="/data/fc/realization/";
    string result_dir = "../results/backup.txt";
    int run_times=1;
    int times=0;
    int batch=4;
    int linear_search_thr=0;  // recommended 50 for formal running
    bool seed_out=true;
    bool gene_ini_pw=false;
    bool real_time_pw=true;
    int numV=1;
    float left_num = 100.0;
    float over_pnodes = 10.0;
    vector<double> Inv_inDeg;
    std::ofstream result_bk;
    
    Argument()
    {
        #ifdef debug
            data={0};
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
            if (argv[i] == string("-Rnd_cost"))
                Rnd_cost = stoi(argv[i + 1]);
            if (argv[i] == string("-eta_0"))
                eta_0 = stof(argv[i + 1]);
            if (argv[i] == string("-batch"))
                batch = stof(argv[i + 1]);
            if (argv[i] == string("-run_times"))
                run_times = stoi(argv[i + 1]);
            if (argv[i] == string("-times"))
                times = stoi(argv[i + 1]);
            if (argv[i] == string("-gene_ini_pw"))
                gene_ini_pw = stoi(argv[i + 1]);
            if (argv[i] == string("-real_time_pw"))
                real_time_pw = stoi(argv[i + 1]);
            if (argv[i] == string("-left_num"))
                left_num = stof(argv[i + 1]);
            if (argv[i] == string("-over_pnodes"))
                over_pnodes = stof(argv[i + 1]);
        }      
    }
    void Initialization()
    {
        activated_nodes.clear();
        activated_nodes.push_back({});
        seed_set.clear();
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
        GraphBase::load_graph_directly_nbr_sorted(dataset_dir, O_graph, R_graph);
        // R_graph = GraphBase::load_graph(dataset_dir, 1);
        // O_graph = GraphBase::load_graph(dataset_dir, 0);
        numV=static_cast<int>(R_graph.size());
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
        }
        else
        {
            cost_file = graph_path + dataset[k] + "_cost_001DEG.txt";
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
            assert(cost[i] >= 0 && cost[i] <= 1);
        }
        inFile.close();
        // cout<<"==================Using uniform costs!!!==============="<<endl;
        // for (size_t i = 0; i < numV; i++)
        // {
        //     cost[i]=1.0;
        // }
        if(gene_ini_pw)
        {
            for (int i = 0; i < run_times; i++){
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
            cout << "generate PO Done" <<endl;
            exit(0);
        }
        return;
    }
    
    void seed_record(vector<int> seeds,int k,int i)
    {
        ofstream out_seeds("../results/newseeds/MINE_" + dataset[k] + "_" + to_string(eta_0) + +"_" + to_string(i) + ".txt", ios::out);
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