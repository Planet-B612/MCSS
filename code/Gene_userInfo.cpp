#include "graph.h"
#include <iostream>
#include "Memory.h"
#include <string>
#include "../dSFMT/dSFMT.h"
// #include "../../dSFMT/dSFMT.c"
using namespace std;

int main(int argn, char **argv)
{
    bool Rnd_cost = false;
    dsfmt_gv_init_gen_rand(static_cast<uint32_t>(time(nullptr)));
    vector<string> dataset = {"facebook", "epinions", "dblp","livejournal", "youtube", "Twitter"};
    vint data = {0};
    string graph_path="/data/fc/graphInfo/";

    for (int i = 0; i < argn; i++)
    {
        if (argv[i] == string("-Rnd_cost"))
        {
            if (argv[i + 1] == string("true"))
                Rnd_cost = true;
            else
                Rnd_cost = false;
        }
    }

    for (int dataset_No: data)
    {
        string dataset_dir = graph_path + dataset[dataset_No];
        Graph O_graph, R_graph;
        GraphBase::load_graph_directly_nbr_sorted(dataset_dir, O_graph, R_graph);
        int numV = O_graph.size();
        if (Rnd_cost)
        {
            ofstream outFile(graph_path + dataset[dataset_No] + "_cost_Rand.txt");
            for (int i = 0; i < numV; i++)
            {
                double cost = dsfmt_gv_genrand_open_close();
                outFile << cost << '\n';
            }
            outFile.close();
            cout<<"generate random cost file done"<<endl;
        }
        else
        {
            ofstream outFile(graph_path + dataset[dataset_No] + "_cost_001DEG.txt");
            for (int i = 0; i < numV; i++)
            {
                double cost = 0.01 + 0.01 * (O_graph[i].size());
                outFile << cost << '\n';
            }
            outFile.close();
            cout<<"generate 0.01*deg cost file done"<<endl;
        }
    }
}