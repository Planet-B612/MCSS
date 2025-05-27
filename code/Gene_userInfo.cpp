#include "CommonFunc.h"
#include "graph.h"
#include <iostream>
#include "Memory.h"
#include "Algorithm.h"
#include <string>
#include "../dSFMT/dSFMT.h"
// #include "../../dSFMT/dSFMT.c"
using namespace std;

int main()
{
    bool Rnd_cost = true;
    bool isReverse = true;
    bool Gene_prob = 0;
    float prob_thresh = 0.2;
    dsfmt_gv_init_gen_rand(static_cast<uint32_t>(time(nullptr)));
    vector<string> dataset = {"sample", "facebook", "dblp", "flickr", "pokec", "orkut", "livejournal", "friendster"};

    for (usint dataset_No = 0; dataset_No < 1; dataset_No++)
    {
        string dataset_dir = "/data/fc/graphInfo/new/" + dataset[dataset_No];
        Graph F_graph = GraphBase::load_graph(dataset_dir, !isReverse);
        int numV = F_graph.size();
        // int numV=65608366;
        if (Gene_prob == true)
        {
            ofstream outFile("/data/fc/graphInfo/new/" + dataset[dataset_No] + "_prob.txt");
            for (int i = 0; i < numV; i++)
            {
                double cost = prob_thresh * dsfmt_gv_genrand_open_close();
                outFile << cost << '\n';
            }
            outFile.close();
            continue;
        }
        else if (Rnd_cost)
        {
            ofstream outFile("/data/fc/graphInfo/new/" + dataset[dataset_No] + "_cost_Rand.txt");
            for (int i = 0; i < numV; i++)
            {
                double cost = dsfmt_gv_genrand_open_close();
                outFile << cost << '\n';
            }
            outFile.close();
        }
        else
        {
            ofstream outFile("/data/fc/graphInfo/new/" + dataset[dataset_No] + "_cost_001DEG.txt");
            for (int i = 0; i < numV; i++)
            {
                double cost = 0.01 + 0.01 * (F_graph[i].size());
                outFile << cost << '\n';
            }
            outFile.close();
        }
    }
}