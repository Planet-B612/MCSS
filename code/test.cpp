#include <iostream>
#include <vector>
#include <chrono>
// #include <algorithm>
// #include "CommonStruc.h"
// #include "../dSFMT/dSFMT.h"
// #include "Memory.h"
// #include "MemoryUsage.h"
// #include "robin_hood.h"
// #include "graph.h"
// #include <malloc.h>

// #include <bits/stdc++.h>
// #include <unordered_map>

using namespace std;

#define Test_Time 0

// int main()
// {
//     string graph_path="/data/fc/graphInfo/new/livejournal";
//     GraphBase::format_graph(graph_path, 0);
//     return 0;
// }

int main()
{
    // int arr[10]={1,2,3,4,5,6,7,8,9};
    // cout<<arr[20]<<endl;
    // exit(0);
    for (int i = 5; i--;)
    {
        cout<<i<<endl;        
    }
    // int num=1000000;
    // vector<int> a={125,234,567};
    // for(auto i=0;i<10;i++)
    // {
    //     cout<<a.capacity()<<endl;
    //     a.push_back(i);
    // }
    // vector<vector<int>> b;
    // for(auto i=0;i<num+20;i++)
    // {
    //     a.push_back(i);
    // }
    // auto start = std::chrono::high_resolution_clock::now();
    // for(auto i=0;i<num;i++)
    // {
    //     b.emplace_back(a.begin()+i,a.begin()+i+10);
    // }
    // auto end = std::chrono::high_resolution_clock::now();
	// std::chrono::duration<double> elapsed = end - start;
    // std::cout << elapsed.count() << " 秒, " << std::endl;

    return 0;
}

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