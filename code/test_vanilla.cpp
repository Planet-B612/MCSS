#include <iostream>
#include <vector>
#include <chrono>
#include <algorithm>
#include "CommonStruc.h"
#include "../dSFMT/dSFMT.h"
#include "Memory.h"
#include "MemoryUsage.h"
// #include "CommonFunc.h"
#include <set>
// #include "test.h"
// #include "robin_hood.h"
// #include "patchmap.hpp"
#include <unordered_map>
#include <algorithm>
#include <random>

using namespace std;

#define Test_Time 0


void quick_sort(vector<int>  &q,int l,int r)
{
    if (l>=r) return; 
/*
l>=r意味着数组中只有一个元素或者是空数组，注意这边一定要用>=，而不能是==，否则会有边界问题
*/
    int x = q[l]; //这边建议直接选q[l]，是后续分界点
    int i = l-1,j = r+1; //注意这边要比边界再往外一位，后续要用do while循环，保持一致
    while(i<j) //大循环，只要两个指针没有相遇，就不断循环
    {
        do i++; while(q[i]<x); //注意这里两个都不能有等于号，要严格大于/小于
        do j--; while(q[j]>x);
        if(i<j) swap(q[i],q[j]); //这里注意要判断是否i<j，否则很可能会有额外交换
    }
    quick_sort(q,l,j); 
/*
递归处理左右两段，同样注意这里的边界问题，上面分界点用x = q[l]的话这里就必须用j和j+1，否则会有边界问题
*/
    quick_sort(q,j+1,r);
    return;
}

int factorial(int n)
{
    if(n==0) return 1;
    return n*factorial(n-1);
}

long double cnk(int n, int k)
{
    if(k>n || k<0 || n<0) {cout<<"value error"<<endl; return 0;};
    if(k==0 || k==n) return 1;
    if(k==1) return n;
    if(k==2) return n*(n-1)/2;
    long double res=1;
    for(int i=0;i<k;i++)
    {
        res=res*(n-i)/(i+1);
    }
    return res;
}

int compare(char i, char j, unordered_map<char, int> &mp) //i>j:1, i<j:0, i=j:-1
{
    if(mp[i]>mp[j])
    return 1;
    else if(mp[i]<mp[j])
    return 0;
    else 
    {
        if(i>j) return 1;
        else return 0;
    }
    // return -1;
}

static bool comp(pair<int, double> a, pair<int, double> b)
{
    return a.second > b.second;
} 

int main(int argn, char **argv)
{
    // vint vec;
    // auto it=lower_bound(vec.begin(), vec.end(), 2);
    // if(it!=vec.end())
    // {
    //     cout<<"it is 2"<<endl;
    // }
    // else
    // {
    //     cout<<"it is not 2"<<endl;
    // }
    // vector<tuple<int,int,int,int>> vec;
    // int a[4]={1,2,3,4};
    // vec.emplace_back(a[0], a[1], a[2], a[3]);
    // vint vec={0,1};//={1};
    // // vec.reserve(100);
    // int *p=&(*(vec.begin()));
    // int *p1=&(*(vec.end()));
    // auto it=lower_bound(vec.begin(), vec.begin()+1002, 2);
    // cout<<it-vec.begin()<<endl;
    // cout<<*it<<endl;
    // int *p2=&(*(it));
    // cout<<p2<<endl;
    // // cout<<p<<endl;
    // // cout<<p1<<endl;
    // vec.insert(it, 2);
    // // cout<<&vec[0]<<endl;
    // // cout<<&vec[1]<<endl;
    // // auto &frset = vec;
    // // auto it=lower_bound(frset.begin(), frset.end(), 45155);
    // // frset.insert(it, 45155);
    // for(int i=2;i<1002;i++)
    // {
    //     if(vec[i]!=0)
    //     {
    //         cout<<"not 0"<<i<<endl;
    //         cout<<vec[i]<<endl;
    //         break;
    //     }
    // }
    // for(int i=0;i<vec.size();i++)
    // {
    //     cout<<vec[i]<<", ";
    // }
    // cout<<endl;
    // int num=1e6;
    // vector<int> vec(num);
    // int a[num];
    // for(int i=0;i<num;i++)
    // {
    //     vec[i]=i;
    //     a[i]=i;
    // }
    // auto start = std::chrono::high_resolution_clock::now();
    // vec.erase(vec.begin());
    // auto now = std::chrono::high_resolution_clock::now();
    // std::chrono::duration<double> elapsed = now - start;
    // cout<<"vec erase time = "<<elapsed.count() << " 秒"<<endl;

    // auto start1 = std::chrono::high_resolution_clock::now();
    // std::move(a+1, a + num, a);
    // auto now1 = std::chrono::high_resolution_clock::now();
    // std::chrono::duration<double> elapsed1 = now1 - start1;
    // cout<<"array erase time = "<<elapsed1.count() << " 秒"<<endl;
    // int n=6,k=2,d=2;
    // for (int i = 0; i < argn; i++)
    // {
    //     if (argv[i] == string("-n"))
    //         n = stoi(argv[i + 1]);
    //     if (argv[i] == string("-k"))
    //         k = stoi(argv[i + 1]);
    //     if (argv[i] == string("-d"))
    //         d = stoi(argv[i + 1]);
    // }
    // if(k>d || d>=n)
    // {
    //     cout<<"value error"<<endl;
    //     exit(1);
    // }
    // long double fresh=0;
    // long double revise=0;
    // fresh=1.0/cnk(n-d,k);
    // cout<<"fresh = "<<fresh<<endl;
    // for(int i=0;i<k+1;i++)
    // {
    //     revise+=1.0*cnk(k,i)*cnk(d,k-i)/cnk(n-d-i,k-i);
    // }
    // revise/=cnk(n,k);
    // cout<<"revise = "<<revise<<endl;

    // vector<vector<int>> vec={
    //     {1,2,3},
    //     {4,5,6},
    //     {7,8,9}
    // };
    // vec.clear();
    // for(auto i:vec)
    // {
    //     for(auto j:i)
    //     {
    //         cout<<j<<", ";
    //     }
    //     cout<<endl;
    // }
    // unordered_map<int, vector<int>> mp;
    // mp.reserve(1);
    // vector<int> 
    // vector<int> vec_1={2,4,5,6,8,10,12,14,16,18,20,24,25,28,32,33,36,40,44,48}, &vec=vec_1; //, vec_1={1,3,5,7,9,11,13,15,17,19};  
    // vec.insert(vec.end(), vec_1.begin(), vec_1.end());
    // vec.insert(vec.end(), vec_1.begin(), vec_1.end());
    // Nodelist vec_1;//, &vec=vec_1;
    // int num=10000;
    // vec_1.reserve(num);
    // for(auto i=0;i<num;i++)
    // {
    //     vec_1.push_back(i);
    // }
    std::random_device rd;
    std::mt19937 gen(rd());
    std::binomial_distribution<int> dist(100, 0.5);
    double avg=0.0, var=0.0, rand=0.0, num=10000;
    vector<double> rands;
    for(int i=0;i<num;i++)
    {
        rands.push_back(dist(gen));
        avg+=rands[i];
        // cout<<rands[i]<<", ";
    }
    // cout<<endl;
    avg=avg/num;
    for(auto rd:rands)
    {
        var+=(rd-avg)*(rd-avg);
    }
    var/=(num-1);
    cout<<"The avg = "<< avg<<", the var = "<<var<<endl;
    // std::shuffle(vec_1.begin(), vec_1.end(), gen);
    // auto start = std::chrono::high_resolution_clock::now();
    // make_min_heap(vec_1);
    // auto now = std::chrono::high_resolution_clock::now();
    // std::chrono::duration<double> elapsed = now - start;
    // cout<<"heap time = "<<elapsed.count() << " 秒"<<endl;
    // std::shuffle(vec_1.begin(), vec_1.end(), gen);
    // auto start1 = std::chrono::high_resolution_clock::now();
    // std::sort(vec_1.begin(), vec_1.end(), [](int a, int b) { return a > b; });
    // auto now1 = std::chrono::high_resolution_clock::now();
    // std::chrono::duration<double> elapsed1 = now1 - start1;
    // cout<<"sort time = "<<elapsed1.count() << " 秒"<<endl;
    // vec.reserve(100); // 预分配 100 个元素的空间
    // for (int i = 0; i < 10; ++i) {
    //     vec.push_back(i);
    // }
    // std::vector<int> &vec_1=vec;

    // 删除 90 个元素
    // vec.erase(vec.begin(), vec.begin() + 3);
    // cout<<static_cast<uint32_t>(time(nullptr))<<endl;
    // dsfmt_t dsfmt; // 初始化生成器，使用固定种子 
    // dsfmt_gv_init_gen_rand(1743605959); // 固定种子 42 // 生成 5 个随机数 
    // for (int i = 0; i < 5; i++) 
    // { 
    //     printf("%f\n", dsfmt_gv_genrand_uint32_range(100)); // 生成随机数 
    // } 
    // return 0; 


    // srand(100);
    // for(auto i=0;i<10;i++)
    // {
    //     cout<<rand()<<", ";
    // }
    // // auto r=rand();
    // // cout<<r<<", "<<endl;
    // // for(int i=0;i<1;i++)
    // // {
    // //     cout<<rand()<<", ";
    // // }
    // cout<<endl;



    // int val=8;
    // auto itt=lower_bound(vec.begin(), vec.end(), val);
    // cout<<itt-vec.begin()<<endl;
    // vec.erase(itt);
    // if(itt!=vec.end())
    // vec.insert(itt, val);
    // else
    // vec.push_back(val);
    // for(auto ch:vec)
    // {
    //     cout<<ch<<", ";
    // }
    // cout<<endl;


    // int b=123;
    // int *p=&b;
    // int c=234;
    // auto p_1=&c;
    // cout<<p<<endl;
    // cout<<p_1<<endl;
    // if(p==p_1)
    // {
    //     cout<<"p==p_1"<<endl;
    // }
    // else
    // {
    //     cout<<"p!=p_1"<<endl;
    // }
    // dsfmt_gv_init_gen_rand(static_cast<uint32_t>(time(nullptr)));
    // vector<double> a;
    // int num=1000;
    // a.reserve(num);
    // auto start = std::chrono::high_resolution_clock::now();
    // for(int i=0;i<num;i++)
    // {
    //     a.push_back(dsfmt_gv_genrand_open_close());
    // }
    // auto now = std::chrono::high_resolution_clock::now();
    // std::chrono::duration<double> elapsed = now - start;
    // cout<<elapsed.count() << " 秒"<<endl;
    // vector<int> vec={2,4,6,8,10,12,14,16,18,20};  //, a={1,2,3,4,5,6};
    // int vec_size=vec.size(), pre_ind=9;
    // int l=0, r=vec_size-1, mid, val=1, res=0;
    // while(l<=r)
    // {
    //     mid=(l+r)/2;
    //     if(vec[mid]>val)
    //     {
    //         res=mid;
    //         r=mid-1;
    //     }
    //     else
    //     {
    //         l=mid+1;
    //     }
    // }
    // cout<<"mid = "<<mid<<"; l = "<<l<<"; res = "<<res<<"; r = "<<r<<endl;
    // vec[pre_ind]=val;
    // for(int i=pre_ind;i>l;i--)
    // {
    //     swap(vec[i], vec[i-1]);
    // }
    // exit(0);


    // vector<char> vec={'a','b','c','d','e','f','g','h','i','j','k','l'};  //size=12
    // unordered_map<char, int> mp;
    // int it=0;
    // for(auto i:vec)
    // {
    //     mp[i]=2*it++;
    //     cout<<i<<": "<<mp[i]<<", ";
    // }
    // cout<<endl;
    // int pre_ind=9;
    // int val=9;
    // mp[vec[pre_ind]]=val;
    // cout<<"The new value at vec["<<pre_ind<<"] is "<<val<<" "<<endl;
    // for(auto ch:vec)
    // {
    //     cout<<ch<<" = "<<mp[ch]<<", ";
    // }
    // cout<<endl;

    // int l=0, r=9, new_ind;
    // while(l<r)
    // {
    //     new_ind=(l+r)/2;
    //     // if(vec[new_ind]>val) r=new_ind;
    //     if(compare(vec[new_ind], vec[pre_ind], mp)>0) r=new_ind;
    //     else l=new_ind+1;
    // }
    // // new_ind++;
    // cout<<"new_ind: "<<new_ind<<endl;
    // for(int i=pre_ind;i>new_ind;i--)
    // {
    //     swap(vec[i], vec[i-1]);
    // }


    // for(auto ch:vec)
    // {
    //     cout<<ch<<" = "<<mp[ch]<<", ";
    // }
    // cout<<endl;

  
    
    


    // int target=10;
    // int l=0, r=vec.size()-1, mid;
    // while(l<r)
    // {
    //     mid=(l+r)/2;
    //     if(vec[mid]==target)
    //     {
    //         cout<<mid<<endl;
    //         break;
    //     }
    //     else if(vec[mid]<target)
    //     {
    //         l=mid+1;
    //     }
    //     else
    //     {
    //         r=mid-1;
    //     }
    // }
    // if(l==r)
    // {
    //     cout<<l<<endl;
    // }
    // else
    // {
    //     cout<<"Not found"<<endl;
    // }
    
    // int new_ind, pre_ind=8, new_val=7;
    // auto it=upper_bound(vec.begin(), vec.end(), new_val);
    // new_ind=it-vec.begin();
    // vec[pre_ind]=new_val;
    // for(int i=pre_ind;i>new_ind;i--)
    // {
    //     swap(vec[i], vec[i-1]);
    // }
    // cout<<it-vec.begin()<<endl;
    

    // int_vec_patchmap mRR;
    // // cout<<&vec[1]<<endl;
    // mRR[334], mRR[1];
    // vector<int> &vec_1=mRR[334];
    // // vec_1.reserve(100);
    // vector<int> &vec_2=mRR[1];
    // cout<<&vec_1<<", "<<&vec_2<<endl;
    // vec_1.insert(vec_1.begin(), vec.begin(), vec.end());
    // vec_2.insert(vec_2.begin(), vec.begin(), vec.end());
    // cout<<&vec_1<<", "<<&vec_2<<endl;
    // vec_nodes.push_back(vec);
    // cout<<&vec_nodes[0][1]<<endl;
    // vec.insert(vec.begin()+2, {100});
    // a.push_back(1);
    // cout<<a.capacity()<<endl;
    // a.push_back(1);
    // cout<<a.capacity()<<endl;
    // a.push_back(1);
    // cout<<a.capacity()<<endl;
    // a.push_back(1);
    // cout<<a.capacity()<<endl;
    // a.push_back(1);
    // cout<<a.capacity()<<endl;
    // cout<<vec.capacity()<<", "<<a.capacity()<<endl;
    // std::random_shuffle(vec.begin(), vec.end());
    // vec.pop_back();
    // vec.pop_back();
    // cout<<vec.back()<<endl;
    // vec.reserve(10);
    // vec.clear();
    // vec.resize(15, 1);
    // vec.erase(vec.begin()+2+1);
    // auto i = vec.size();
    // cout<<i<<", "<<vec.size()<<endl;
    // for (int i=0;i<10; i++)
    // {
    //     cout<<a[i]<<", ";
    // }
    // cout<<endl;
    // cout<<vec.size()<<endl;
    // quick_sort(vec,0, vec.size()-1);
    
    
    // for(auto i:vec)
    // {
    //     cout<<i<<", ";
    // }
    // cout<<endl;
    // cout<<endl;
    // cout<<vec.end()-vec.begin()<<endl;
    // auto lb=lower_bound(vec.begin(),vec.end(),10);
    // cout<<lb-vec.begin()<<endl;
    // cout<<*lb<<endl;
    // lb=lower_bound(vec.begin(),vec.end(),9);
    // cout<<lb-vec.begin()<<endl;
    // cout<<*lb<<endl;
    // vector<vector<int>> bi_vec= {{1,2,3},{4,5,6},{7,8,9}};
    // for(auto i=2;i>1;i--)
    // {
    //     vector<int>().swap(bi_vec[i]);
    // }
    // for(auto vec:bi_vec)
    // {
    //     for(auto i:vec)
    //     {
    //         cout<<i<<", ";
    //     }
    //     cout<<endl;
    // }
    // cout<<bi_vec.size()<<endl;

//     vector<int> a(1000000);
//     auto pre_time=std::chrono::high_resolution_clock::now();
    

// auto start = std::chrono::high_resolution_clock::now();
// 			std::chrono::duration<double> elapsed = start - pre_time;
//             cout<<elapsed.count() << " 秒"<<endl;
    return 0;
}