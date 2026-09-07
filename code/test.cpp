#define DEBUG
#define PREFETCH

#include <atomic>
#include <mutex>
#include <omp.h>
#include <immintrin.h>
#include <bitset>
#include <iostream>
#include <vector>
#include <chrono>
// #include <algorithm>
#include "../dSFMT/dSFMT.h"
#include "Argument.h"
#include "CommonStruc.h"
#include "CommonFunc.h"
#include "../dSFMT/dSFMT.h"
#include "mRRcollection.h"
// #include "Algorithm.h"
// #include "graph.h"
#include <cstring>
#include "Timer.h"
#include "Memory.h"
// #include "MemoryUsage.h"
#include <queue>
#include <algorithm>   // std::shuffle
#include <random>     
#include <ctime>
#include <ratio>
#include <algorithm>
#include <immintrin.h> // SIMD intrinsics
#include <numeric>
#include <cmath>
#include <cassert>
// #include "test_ic.h"


using namespace std;
using namespace std::chrono;

#define Test_Time 0

#include <iostream>
#include <vector>

int main(int argn, char **argv)
{   
    vint seeds; 
    int RR_num=1000;
    Argument arg;
    arg.dataset_No=4;
    arg.model="IC";
    int root_num_increase=10, seed_num=500, round_num=20;
    root_num=10;
    for(int i=0;i<argn;i++)
    {
        if(argv[i]==string("-dataset_No"))
        {
            arg.dataset_No=stoi(argv[i+1]);
        }
        if(argv[i]==string("-model"))
        {
            arg.model=argv[i+1];
        }
        if(argv[i]==string("-RR_num"))
        {
            RR_num=stoi(argv[i+1]);
        }
        if(argv[i]==string("-root_num"))
        {
            root_num=stoi(argv[i+1]);
        }
        if(argv[i]==string("-root_num_increase"))
        {
            root_num_increase=stoi(argv[i+1]);
        }
        if(argv[i]==string("-seed_num"))
        {
            seed_num=stoi(argv[i+1]);
        }
        if(argv[i]==string("-round_num"))
        {
            round_num=stoi(argv[i+1]);
        }
    }
    R_graph.clear(), O_graph.clear();
    dsfmt_gv_init_gen_rand(static_cast<uint32_t>(time(nullptr)));
    // arg.arg_update(argn, argv);
    arg.Initialization();
    arg.load_cost_graph(arg.dataset_No);
    mRRcollection RR(arg);
    RR.vv_virtual_roots.resize(RR_num);
    RR.build_n_mRRsets_tree(RR_num,0);

    for(int round=0;round<round_num;round++)
    {
        root_num+=root_num_increase;
        for (int j = 0; j < seed_num; j++)
        {
            int seed = dsfmt_gv_genrand_uint32_range(arg.numV);
            while (__Activated[seed])
            {
                seed = dsfmt_gv_genrand_uint32_range(arg.numV);
            }
            __Activated[seed] = 1;
            seeds.push_back(seed);
        }
        RR.realization(seeds);
        cout<<"updating mRR-sets in round "<<round<<endl;
        RR.build_n_mRRsets_tree(RR_num,0);
        cout<<"FR_reverse_check begins."<<endl;
        if(RR.FR_reverse_check(""))
        {
            cout<<"\033[31m"<<__LINE__<<"Error"<<"\033[0m"<<": FR_reverse_check failed!"<<endl;
            exit(1);
        }
        cout<<"Finished updating and adding roots."<<endl;
        cout<<"num_update: "<<num_update<<", num_add_root: "<<num_add_root<<endl;
        seeds.clear();
        num_update=0;
        num_add_root=0;
    }
}

// int main() {
//     dsfmt_gv_init_gen_rand(static_cast<uint32_t>(time(nullptr)));

//     __m512i minus_one = _mm512_set1_epi32(-1);
//     __m512i zero = _mm512_setzero_si512();
//     vint_aligned arr;
//     int Vnum=1e6, copy_size=10000, RR_size=1000;
//     vint __vecNewTree(Vnum, -1), eraseNodes;
//     mRRset mRR_copy(copy_size);
//     for(int i=0;i<copy_size;i++)
//     {
//         for(int j=0;j<RR_size;j++)
//         {
//             mRR_copy[i].push_back(dsfmt_gv_genrand_uint32_range(Vnum));
//         }
//         arr.push_back(i);
//     }
//     cout<<"mRR_copy has been generated"<<endl;

//     high_resolution_clock::time_point startTime = high_resolution_clock::now();	
//     for (auto &RR : mRR_copy)
//     {
//         for (auto &node : RR)
//         {
//             if (__vecNewTree[node] < 0) // previously in mRR but now not in mRR
//             {
//                 auto it=std::lower_bound(RR.begin(), RR.end(), 400);
//                 if(it != RR.end() && *it == 400)
//                 {
//                     RR.erase(it);
//                 }
//             }
//             // __vecNewTree[node] = -1;
//         }
//         // auto it=std::lower_bound(RR.begin(), RR.end(), 400);
//         // if(it != RR.end() && *it == 400)
//         // {
//         //     RR.erase(it);
//         // }
//     }
//     high_resolution_clock::time_point vallina_time = high_resolution_clock::now();	
//     cout<< "The time for vallina is "<<std::chrono::duration<double>(vallina_time - startTime).count()<<" s"<<endl;

//     // for (auto &RR : mRR_copy)
//     // {
//     //     for (auto &node : RR)
//     //     {
//     //         if (__vecNewTree[node] < 0) // previously in mRR but now not in mRR
//     //         {
//     //             simd_ordered_erase(RR, 600);
//     //         }
//     //         __vecNewTree[node] = -1;
//     //     }
//     // }
//     // high_resolution_clock::time_point simd_time = high_resolution_clock::now();	
//     // cout<< "The time for simd is "<<std::chrono::duration<double>(simd_time-vallina_time).count()<<" s"<<endl;

//     for(int i=0;i<copy_size;i++)
//     {
//         auto &RR = mRR_copy[i];
//         int j=0;
//         // if(i%100==0)
//         // {
//         //     cout<<"processing "<<i<<"th mRR..."<<endl;
//         // }
//         auto size=RR.size();
//         for (; j + 16 <= size; j += 16)
//         {
//             __m512i idx = _mm512_load_si512(reinterpret_cast<const void*>(RR.data() + j));

//             __m512i vals = _mm512_i32gather_epi32(idx, __vecNewTree.data(), 4);

//             __mmask16 mask = _mm512_cmplt_epi32_mask(vals, zero);

//             // scatter: __vecNewTree[node] = -1
//             // _mm512_i32scatter_epi32(__vecNewTree.data(), idx, minus_one, 4);

//             // 把需要 erase 的 node 收集起来
//             alignas(64) int nodes[16];
//             _mm512_store_si512(reinterpret_cast<void*>(nodes), idx);

//             while (mask)
//             {
//                 int k = __builtin_ctz(mask);
//                 auto it=std::lower_bound(RR.begin(), RR.end(), 400);
//                 if(it != RR.end() && *it == 400)
//                 {
//                     RR.erase(it);
//                 }
//                 mask &= mask - 1;
//             }
//         }

//         for (; j < size; ++j)
//         {
//             // int node = RR[j];
//             // if (__vecNewTree[node] < 0)
//             // {
//             //     auto it=std::lower_bound(RR.begin(), RR.end(), 400);
//             //     if(it != RR.end() && *it == 400)
//             //     {
//             //         RR.erase(it);
//             //     }
//             // }
//             __vecNewTree[RR[j]] = -1;
//         }
//     }
//     // cout<<"erasing begins."<<endl; 
//     int k=0;
//     // for (auto &RR : mRR_copy)
//     // {
//         // for (auto &node : RR)
//         // {
//         //     if (__vecNewTree[node] < 0) // previously in mRR but now not in mRR
//         //     {
//         //         // auto it=std::lower_bound(RR.begin(), RR.end(), 400);
//         //         // if(it != RR.end() && *it == 400)
//         //         // {
//         //         //     RR.erase(it);
//         //         // }
//         //     }
//         //     // __vecNewTree[node] = -1;
//         // }
//     //     simd_ordered_erase(RR, 400);
//     // }
//     high_resolution_clock::time_point SIMD_time = high_resolution_clock::now();	

//     cout<< "The time for SIMD is "<<std::chrono::duration<double>(SIMD_time - vallina_time).count()<<" s"<<endl;
//     return 0;
// }



// int main(int argn, char **argv)
// {
//     int num=3;
//     Graph g(num);
//     for(int i=0;i<num;i++)
//     {
//         for(int j=0;j<num;j++)
//         {
//             g[i].push_back(j);
//         }
//     }
//     for(const auto &nbrs:g)
//     {
//         for(const auto &nbr:nbrs)
//         {
//             cout<<nbr<<", ";
//         }
//         cout<<endl;
//     }

       // 严格递增、无重复元素的数组
    // vint_aligned data = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160};
    // int n = data.size();

    // // cout<<simd_ordered_insert(data, 711)<<endl;
    // simd_ordered_insert(data, 711);

    // // cout<<simd_ordered_insert(data, 40)<<endl;
    // simd_ordered_insert(data, 40);

    // // cout<<simd_ordered_insert(data, 5)<<endl;
    // simd_ordered_insert(data, 5);

    // // cout<<simd_ordered_insert(data, 120)<<endl;
    // simd_ordered_insert(data, 120);

    // // cout<<simd_ordered_insert(data, 20)<<endl;
    // simd_ordered_insert(data, 20);
    // // cout<<simd_ordered_insert(data, 60)<<endl;
    // simd_ordered_insert(data, 60);

    // high_resolution_clock::time_point startTime = high_resolution_clock::now();	
    // std::random_device rd;                 // 硬件熵源（一次就够）
    // std::mt19937 g(rd());
    // int n=1e6;
    // vector<int> vec;
    // for(int i=0;i<n;i++)
    // {
    //     vec.push_back(i);
    // }
    // std::shuffle(vec.begin(), vec.end(), g);
    // sort(vec.begin(),vec.end());
    // auto now = std::chrono::high_resolution_clock::now();
    // cout<<"The time is "<<std::chrono::duration<double>(now - startTime).count()<<endl;
    // for(int i=2; i<20; i++)
    // {
    //     cout<<endl;
    // }
    // std::string filename = "/data/fc/graphInfo/Twitter";
    // std::vector<std::string> lines;
    // std::string line;
    
    // // 读取文件
    // std::ifstream inFile(filename);
    // while (std::getline(inFile, line)) {
    //     lines.push_back(line);
    // }
    // inFile.close();
    
    // // 修改第一行
    // int firstNum = std::stoi(lines[0]);
    // lines[0] = std::to_string(firstNum) + " " + std::to_string(lines.size());
    
    // // 写回文件
    // std::ofstream outFile(filename);
    // for (const auto& l : lines) {
    //     outFile << l << std::endl;
    // }
    // outFile.close();
    
    // std::cout << "修改完成! 新第一行: " << lines[0] << std::endl;

    // double eps = 0.7;
    // int __eta_left = 60;
    // int eta = 20000;
    // int __numV_left = 4846609 + __eta_left - eta;
    // double delta=eps/(100.0*(1-1/2.71828)*(1-eps)*__eta_left);
    // int batch_size = 8;
    // double approx=1.0-power((1-1.0/batch_size),batch_size);
    // double elta=eps/(100.0*(1-1/2.71828)*(1-eps)*__eta_left);
    // double eps_hat=99.0*eps/(100.0-eps);
    // const double alpha = sqrt(log(6.0 / delta));
    // const double beta = sqrt((logcnk(__numV_left, batch_size) + log(6.0 / delta)) / approx);
    // int theta = 2 * (alpha + beta)* (alpha + beta);
    // cout << "alpha: " << alpha << ", beta: " << beta << ", theta: " << theta << endl;
    // Argument arg;
    // arg.dataset_No = 0;
    // arg.real_time_pw = true;  // generate possible world in real time
    // string graph_path="/data/gongyao/graphInfo/sample";
    // R_graph.clear(), O_graph.clear();  // global variables
    // GraphBase::load_graph_directly_nbr_sorted(graph_path, O_graph, R_graph);
    // __Activated.clear();
    // activated_nodes.clear();
    // activated_nodes.push_back({});
    // seed_set.clear();
    // cost.clear();
    // dsfmt_gv_init_gen_rand(static_cast<uint32_t>(time(nullptr)));  // the type must be uint32_t, to be accord with the function definition
    // arg.arg_update(argn, argv);
    // arg.load_cost_graph(0);
    // mRRcollection RR(arg);
    // root_num = 3;
    // // RR[0] = {}
    // RR.build_n_mRRsets_tree(1);
    // // RR.output_info(0);
    // RR.add_root(0,2);
    // RR.add_root(0,2);
    // // RR.output_info(0);
    // // RR.delete_root(0,1);
    // // RR.output_info(0);
    // __Activated[1] = true;
    // __Activated[2] = true;
    // vector<int> del_nodes = {1, 2};
    // RR.mRR_update(0,del_nodes);
    // __Activated[3] = true;
    // __Activated[4] = true;
    // return 0;
// }

// int main()
// {
//     // int arr[10]={1,2,3,4,5,6,7,8,9};
//     // cout<<arr[20]<<endl;
//     // exit(0);
//     for (int i = 5; i--;)
//     {
//         cout<<i<<endl;        
//     }
//     // int num=1000000;
//     // vector<int> a={125,234,567};
//     // for(auto i=0;i<10;i++)
//     // {
//     //     cout<<a.capacity()<<endl;
//     //     a.push_back(i);
//     // }
//     // vector<vector<int>> b;
//     // for(auto i=0;i<num+20;i++)
//     // {
//     //     a.push_back(i);
//     // }
//     // auto start = std::chrono::high_resolution_clock::now();
//     // for(auto i=0;i<num;i++)
//     // {
//     //     b.emplace_back(a.begin()+i,a.begin()+i+10);
//     // }
//     // auto end = std::chrono::high_resolution_clock::now();
// 	// std::chrono::duration<double> elapsed = end - start;
//     // std::cout << elapsed.count() << " 秒, " << std::endl;

//     return 0;
// }

/*
int main_a()
{
    string graph_path="/data/fc/graphInfo/new/livejournal";
    Graph O_graph = GraphBase::load_graph(graph_path, 0);

    int num=2000000, map_num=1000;
    std::unordered_map<int,Nodelist> a;
    vector<std::unordered_map<int,Nodelist>> vec_map;
    vector<vector<int>> vec_key;
    int key_num=0;
    int cnt_vec=0, cnt=0;
    disp_mem_usage();
    
    for(auto i=0;i<num;i++)
    {
        // if(key_num<num/2)
        // {
        //     vec_key.insert(vec_key.end(), O_graph[i].begin(), O_graph[i].end());
        //     key_num+=O_graph[i].size();
        // }
        // a[i]={i+1, i+2, i+3};
        // for(auto j=0;j<num;j++)
        // {
        //     if(std::fmod(j,j/100)==0)
        //     {
        //         a[i].push_back(j);
        //     }
        // }
        // int key=0;
        // if(O_graph[i].empty())
        // {
        //     key=rand();
        // }
        // else
        // {
        //     key=O_graph[i][0];
        // }
        a[i]=O_graph[i];
        // for(auto j=0;j<10;j++)
        // {
        //     a[i].push_back(i+j);
        // }
        // a[i]={i-5,i-4,i-3,i-2,i-1,i,i+};
    }
        
    // vec_key.resize(num/2);
    
    for(auto i=0;i<map_num;i++)
    {
        std::unordered_map<int,Nodelist> a_i;
        for(auto j=(num/map_num)*(i);j<(num/map_num)*(i+1);j++)
        {
            // a_i[j]={j+1, j+2, j+3};
            //  for(auto k=0;k<num;k++)
            // {
            //     if(std::fmod(k,rand()/100)==0)
            //     {
            //         a_i[j].push_back(k);
            //     }
            // }
            // a_i[j].resize(50,j);
            // int key=0;
            // if(O_graph[j].empty())
            // {
            //     key=rand();
            // }
            // else
            // {
            //     key=O_graph[j][0];
            // }
            a_i[j]=O_graph[j];
            // for(auto k=0;k<10;k++)
            // {
            //     a_i[j].push_back(i+k);
            // }
        }
        vec_map.push_back(a_i);
    }
    for(auto &mp:vec_map)
    {
        vector<int> keys;
        for(auto entry:mp)
        {
            keys.push_back(entry.first);
        }
        vec_key.push_back(keys);
    }
    
    // for(auto i=0;i<num;i++)
    // {
    //     a.erase(i);
    // }
    
    auto vec_map_size=vec_map.size();
    for(auto i=0;i<vec_map_size;i++)
    {
        for(auto j:vec_key[i])
        {
            vec_map[i].erase(j);
        }
    }


    // for(auto key:vec_key)
    // {
    //     if(vec_map[50].find(key)!=vec_map[50].end())
    //     {
    //         cnt++;
    //     }
    // }

    disp_mem_usage();
    // for(auto mp:vec_map)
    // {
    //     for(auto entry:mp)
    //     {
    //         cnt_vec+=entry.second.size();
    //     }
    // }
    // for(auto entry:a)
    // {
    //     cnt+=entry.second.size();
    // }
    // cout<<cnt<<", "<<cnt_vec<<endl;
    
    // for(auto k=0;k<10;k++)
    // {
    //     for(auto i=0;i<num/map_num-1;i++)
    //     {
    //         swap(vec_map[i],vec_map[i+1]);
    //     }
    // }
    disp_mem_usage();


    // int*k;
    // k=new int (3);
    // cout<<*k<<endl;
    // delete k;
    // int k=3;
    // vector<int> vec={1,2,3,45,5,66,6,7,43};
    // vec.erase(vec.begin()+1,vec.begin()+k+1);
    // for(auto i:vec)
    // {
    //    cout<<i<<", ";
    // }
    //     cout<<endl;
    // vector<vector<int>> vec;
    // vec.push_back({1,2,3,45,5});
    // vec.back().push_back(999);
    // for(auto list:vec)
    // {
    //     for(auto n:list)
    //     {
    //         cout<<n<<", ";
    //     }
    //     cout<<endl;
    // }
    // int numV=1e7, cnt=0;
    // dsfmt_gv_init_gen_rand(static_cast<uint32_t>(time(nullptr)));
    // robin_hood::unordered_map<int,int> a;
    // vector<int> aha;
    // for(auto i=0;i<1000000;i++)
    // {
    //     a[dsfmt_gv_genrand_uint32_range(numV)]=i;
    //     aha.push_back(i);
    // }
    
    // auto start = std::chrono::high_resolution_clock::now();
    // vector<int> state(numV,-1);
    // vector<int> aha(numV);
    // for(auto i=0;i<1e7;i++)
    // {
    //     if(a.find(i)==a.end())
    //     {
    //         cnt++;
    //     }
    // }
    // for(auto entry:a)
    // {
    //     state[entry.first]=1;
    // }
    // for(auto entry:a)
    // {
    //     aha[entry.first]=2;
    // }
    // auto end = std::chrono::high_resolution_clock::now();
	// std::chrono::duration<double> elapsed = end - start;
    // std::cout << elapsed.count() << " 秒, " <<"cnt = "<<cnt<< std::endl;


    return 0;
}


int main_1() {
    const int ARRAY_SIZE = 100; // 每个数组的大小
    const int NUM_ARRAYS = 100000; // 数组的数量
 
    int** arrays = new int*[NUM_ARRAYS]; // 指向数组的指针数组
 
    for (int i = 0; i < NUM_ARRAYS; ++i) {
        arrays[i] = new int[ARRAY_SIZE](); // 使用()初始化为0
    }
 
    // 使用数组...
 
    // 清理内存
    disp_mem_usage();
    for (int i = 0; i < NUM_ARRAYS; ++i) {
        delete[] arrays[i];
    }
    delete[] arrays;
    disp_mem_usage();

    int** arrays_ = new int*[NUM_ARRAYS]; // 指向数组的指针数组
 
    for (int i = 0; i < NUM_ARRAYS; ++i) {
        arrays_[i] = new int[ARRAY_SIZE](); // 使用()初始化为0
    }
    disp_mem_usage();
 
    return 0;
}

int main_2(int argn, char **argv)
{


    // std::unordered_map<int,vector<int>> a; //std::vector<int>
    std::unordered_map<int,vector<int>*> *a= new std::unordered_map<int,vector<int>*>();
    vector<int> *p=new vector<int>();
    vector<vector<int>*> vec_ptr= vector<vector<int>*>();
    // (*vec_ptr).resize(1000001);
    disp_mem_usage();
    for(auto i=1;i<1000001;i++)
    {
        vector<int> *vec=new vector<int>();
        for(auto j=0;j<10;j++)
        {
            vec->push_back(i*j);
        }
        // (*vec_ptr)[i]=vec;
        vec_ptr.push_back(vec);
    }
    // cout<<(vec_ptr).size()<<endl;
    for(auto i=1;i<1000000;i++)
    {

        (*a)[i]=(vec_ptr)[i];
        // vec_ptr->push_back(vec);  
    }
    // for(auto j:*(*a)[2])
    // {
    //     cout<<j<<", ";
    // }
    // cout<<endl;


    // vector<int> *q=new vector<int>();
    // vector<int> *r=new vector<int>();
    // vector<int> *s=new vector<int>();
    // vector<int> *t=new vector<int>();
    
    // cout<<&p<<endl;
    // // p=
    // disp_mem_usage();
    // for(int i=0;i<1000006;i++)
    // {
    //     p->push_back(i);
    //     // q->push_back(i+1);
    //     // r->push_back(i+2);
    //     // s->push_back(i+3);
    //     // r->push_back(i+10);
    // }
    // for(int i=0;i<1000000;i++)
    // {
    //     vector<int> vec;
    //     vec={(*p)[i], (*p)[i+1], (*p)[i+2], (*p)[i+3], (*p)[+4], (*p)[i+5]};
    // }
    // for(int i=1000000;i<2000000;i++)
    // {
    //     q->push_back(i);
    // }
    //     for(int i=2000000;i<3000000;i++)
    // {
    //     r->push_back(i);
    // }
    //         for(int i=3000000;i<4000000;i++)
    // {
    //     s->push_back(i);
    // }
    //         for(int i=4000000;i<5000000;i++)
    // {
    //     t->push_back(i);
    // }
    // a[0]=*p;
    // a[1]=*q;
    // a[2]=*r;
    // a[3]=*s;
    // a[4]=*t;
    // pair<int,std::vector<int>> pr=make_pair(0,*p);
    // a[0]=*p;
    // p->push_back(11),p->push_back(22),p->push_back(33);
    // cout<<&p[0]<<endl;
    // for(auto i:*p)
    // {
    //     cout<<i<<endl;
    // }
    // int *ind=new int[1000000];
    // vector<pair<int, int>> *pr=new vector<pair<int, int>>();
    // for(int i=0;i<1000000;i++)
    // {
    //     a[i]=(*p)[i];
    //     // a[i]=i;
    // }
    // for(auto i:*pr)
    // {
    //     a.insert(i);
    // }
    // for(auto i=0; i<1000000;i++)
    // {
    //     a.erase(i);
    // }
    // vector<double> *vec=new vector<double>(10000000,10.2365);
    
    // a.erase(0);
    disp_mem_usage();
    // a.clear();
    // for(int i=0;i<1000001;i++) 
    // {
    //     // delete (*a)[i];
    //     (*a)[i]=NULL;
    //     // (*a).erase(i);
    // }
    delete a;
    delete p;
    for(auto &ptr:vec_ptr)
    {
        // vector<int>().swap(*ptr);
        delete ptr;
    }
    malloc_trim(0);
    // delete vec_ptr;
    // delete q;
    // delete r;
    // delete s;
    // delete t;
    // malloc_trim(0);
    // delete vec;
    // delete [] ind;
    // vec.clear();
    // cout<<vec.size();
    // delete pr;
    disp_mem_usage();
    // cout<<&p<<endl;
    // if(p==NULL)
    // {
    //     cout<<"is null"<<endl;
    // }
    // for(auto i:*p)
    // {
    //     cout<<i<<endl;
    // }


    cout<<"Printing the memory after creating a new map b: "<<endl;
    std::unordered_map<int,vector<int>*> *b= new std::unordered_map<int,vector<int>*>();
    // vector<int> *p=new vector<int>();
    vector<vector<int>*> vec_ptr_n= vector<vector<int>*>();
    // (*vec_ptr).resize(1000001);
    for(auto i=1000001;i<2000003;i++)
    {
        vector<int> *vec=new vector<int>();
        for(auto j=10;j<20;j++)
        {
            vec->push_back(i*j);
        }
        // (*vec_ptr)[i]=vec;
        vec_ptr_n.push_back(vec);
    }
    // cout<<(vec_ptr).size()<<endl;
    for(auto i=1;i<1000000;i++)
    {
        (*b)[i]=(vec_ptr_n)[i];
        // vec_ptr->push_back(vec);  
    }
    disp_mem_usage();
    for(auto &ptr:vec_ptr_n)
    {
        // vector<int>().swap(*ptr);
        delete ptr;
    }
    // (*b).clear();
    delete b;
    malloc_trim(0);
    disp_mem_usage();

    cout<<"Printing the memory after creating a new map c: "<<endl;
    std::unordered_map<int,vector<int>*> *c= new std::unordered_map<int,vector<int>*>();
    // vector<int> *p=new vector<int>();
    vector<vector<int>*> vec_ptr_m= vector<vector<int>*>();
    // (*vec_ptr).resize(1000001);
    for(auto i=2000001;i<2000004;i++)
    {
        vector<int> *vec=new vector<int>();
        for(auto j=20;j<30;j++)
        {
            vec->push_back(i*j);
        }
        // (*vec_ptr)[i]=vec;
        vec_ptr_m.push_back(vec);
    }
    // cout<<(vec_ptr).size()<<endl;
    for(auto i=1;i<1000000;i++)
    {
        (*c)[i]=(vec_ptr_m)[i];
        // vec_ptr->push_back(vec);  
    }
    disp_mem_usage();
    for(auto &ptr:vec_ptr_m)
    {
        // vector<int>().swap(*ptr);
        delete ptr;
    }
    // (*c).clear();
    delete c;
    malloc_trim(0);
    disp_mem_usage();

    cout<<"Printing the memory after creating a new map d: "<<endl;
    std::unordered_map<int,vector<int>*> *d= new std::unordered_map<int,vector<int>*>();
    // vector<int> *p=new vector<int>();
    vector<vector<int>*> vec_ptr_k= vector<vector<int>*>();
    // (*vec_ptr).resize(1000001);
    for(auto i=3000001;i<4000004;i++)
    {
        vector<int> *vec=new vector<int>();
        for(auto j=1;j<11;j++)
        {
            vec->push_back(i*j);
        }
        // (*vec_ptr)[i]=vec;
        vec_ptr_k.push_back(vec);
    }
    // cout<<(vec_ptr).size()<<endl;
    for(auto i=1;i<1000000;i++)
    {
        (*d)[i]=(vec_ptr_k)[i];
        // vec_ptr->push_back(vec);  
    }
    disp_mem_usage();
    for(auto &ptr:vec_ptr_k)
    {
        // vector<int>().swap(*ptr);
        delete ptr;
    }
    // d->clear();
    delete d;
    malloc_trim(0);
    // malloc_trim(0);
    disp_mem_usage();


    cout<<"Printing the memory after creating a new map e: "<<endl;
    std::unordered_map<int,vector<int>*> *e= new std::unordered_map<int,vector<int>*>();
    // vector<int> *p=new vector<int>();
    vector<vector<int>*> vec_ptr_l= vector<vector<int>*>();
    // (*vec_ptr).resize(1000001);
    for(auto i=4000001;i<5000004;i++)
    {
        vector<int> *vec=new vector<int>();
        for(auto j=22;j<32;j++)
        {
            vec->push_back(i*j);
        }
        // (*vec_ptr)[i]=vec;
        vec_ptr_l.push_back(vec);
    }
    // cout<<(vec_ptr).size()<<endl;
    for(auto i=1;i<1000000;i++)
    {
        (*e)[i]=(vec_ptr_l)[i];
        // vec_ptr->push_back(vec);  
    }
    disp_mem_usage();
    for(auto &ptr:vec_ptr_l)
    {
        // vector<int>().swap(*ptr);
        delete ptr;
    }
    // e->clear();
    delete e;
    malloc_trim(0);
    disp_mem_usage();





    // mRRset mRR;
    // robin_hood::unordered_map<int,std::vector<int>> a;
    // a[1]={1,2,3,4,5,6,7,8,9};
    // a[2]={1,2,3,4,5,6,7,8,9};
    // a[3]={1,2,3,4,5,6,7,8,9};
    // a[4]={1,2,3,4,5,6,7,8,9};
    // a[5]={1,2,3,4,5,6,7,8,9};
    // a[6]={1,2,3,4,5,6,7,8,9};
    // a[7]={1,2,3,4,5,6,7,8,9};
    // a[8]={1,2,3,4,5,6,7,8,9};
    // // std::cout<<", "<<&<<endl;
    // cout<<", "<<&a[8][0]<<", "<<&a[8]<<endl;
    // a[9]={1,2,3,4,5,6,7,8,9};
    // a[10]={1,2,3,4,5,6,7,8,9};
    // cout<<a.begin()->first<<endl;
    // a.erase(5);
    // a.erase(8);
    // a.erase(1);
    // a.erase(3);
    // a.clear();
    // a[15]={1,2,3,4,5,6,7,8,9};
    // a[9]={1,2,3,4,5,6,7,8,9};
    // a[10]={1,2,3,4,5,6,7,8,9};
    // a.clear();
    



 
    // auto start = std::chrono::high_resolution_clock::now();



    #if Test_Time!=0
    {    
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        std::cout << elapsed.count() << " 秒" << std::endl;
    }
    #endif
    return 0;

}
*/


