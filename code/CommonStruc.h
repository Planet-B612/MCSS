#pragma once
#include <set>
// #include "robin_hood.h"
// #include "patchmap.hpp"
#include "../hashmap_dir/sparsepp/spp.h"
using spp::sparse_hash_map;
// #include "../sparsehash11/sparse_hash_map"
// #include "../parallel_hashmap/phmap.h"
// using phmap::flat_hash_map;

#if !defined(DSFMT_MEXP)
#ifdef __GNUC__
#define DSFMT_MEXP 19937
#endif
#endif

#if 0  // STL
typedef std::unordered_map<int, std::vector<int>>  int_vec_patchmap;
typedef std::unordered_map<int, int>  int_int_patchmap;
typedef std::unordered_map<int, double>  int_dbl_patchmap;
typedef std::unordered_map<int, std::pair<double,int>>  int_dpair_patchmap;
typedef std::unordered_map<int, std::pair<int,int>>  Chain;
#endif

#if 0  // robin_hood
typedef robin_hood::unordered_map<int, std::vector<int>>  int_vec_patchmap;
typedef robin_hood::unordered_map<int, int>  int_int_patchmap;
typedef robin_hood::unordered_map<int, double>  int_dbl_patchmap;
typedef robin_hood::unordered_map<int, std::pair<double,int>>  int_dpair_patchmap;
typedef robin_hood::unordered_map<int, std::pair<int,int>>  Chain;
#endif

#if 1  // sparsepp
typedef sparse_hash_map<int, std::vector<int>>  int_vec_patchmap;
typedef sparse_hash_map<int, int>  int_int_patchmap;
typedef sparse_hash_map<int, double>  int_dbl_patchmap;
typedef sparse_hash_map<int, std::pair<double,int>>  int_dpair_patchmap;
typedef sparse_hash_map<int, std::tuple<int,int, int, int>>  int_location_map;
typedef sparse_hash_map<int, std::pair<int,int>>  Chain;
typedef std::set<std::pair<std::pair<double,double>,int>> set_timer;
#endif
#if 0 // parallel hash map
typedef flat_hash_map<int, std::vector<int>>  int_vec_patchmap;
typedef flat_hash_map<int, int>  int_int_patchmap;
typedef flat_hash_map<int, double>  int_dbl_patchmap;
typedef flat_hash_map<int, std::pair<double,int>>  int_dpair_patchmap;
typedef flat_hash_map<int, std::pair<int,int>>  Chain;
#endif

//=======================================================================================================

#if 0  // patchmap
typedef whash::patchmap<int, std::vector<int>>  int_vec_patchmap;
typedef whash::patchmap<int, int>  int_int_patchmap;
typedef whash::patchmap<int, double>  int_dbl_patchmap;
typedef whash::patchmap<int, std::pair<double,int>>  int_dpair_patchmap;
typedef whash::patchmap<int, std::pair<int,int>>  Chain;
#endif

#if 0
// google's sparse_hash_map C++11
#define GOOGLE_SPARSE_11
typedef google::sparse_hash_map<int, std::vector<int>>  int_vec_patchmap;
typedef google::sparse_hash_map<int, int>  int_int_patchmap;
typedef google::sparse_hash_map<int, double>  int_dbl_patchmap;
typedef google::sparse_hash_map<int, std::pair<double,int>>  int_dpair_patchmap;
typedef google::sparse_hash_map<int, std::pair<int,int>>  Chain;
#endif


/// A set of reverse reachable sets
// typedef int_vec_patchmap mRRset;
// typedef std::vector<mRRset> mRRsets;
// typedef std::vector<Chain> Mchain;
// typedef std::vector<Mchain> Mchains;

typedef std::vector<int> vint;
typedef std::vector<vint> vvint;
typedef std::vector<bool> vbool;

/// Node list
typedef std::vector<int> Nodelist;
/// Edge structure: neighbor id, the edge weight
// typedef std::pair<int, float> Edge;
/// Edgelist structure from one source/target node
typedef std::vector<int> Edgelist;
/// Graph structure
typedef std::vector<Edgelist> Graph;
/// One forward reachable set
// typedef spp::sparse_hash_set<int> int_unordered_set;
// typedef spp::sparse_hash_set<int> FRset;
typedef std::vector<int> FRset;
/// A set of forward reachable sets
typedef std::vector<FRset> FRsets;
typedef std::vector<FRsets> FRcollection;
/// One reverse reachable set
typedef std::vector<vint> mRRset;
// typedef std::vector<int_vec_patchmap> mRRset;
typedef std::vector<mRRset> mRRsets;
typedef std::vector<Chain> Mchain;
typedef std::vector<Mchain> Mchains;


// /// Define the way of storing results
// typedef std::tuple<double, double, double, double> Res;
// typedef std::vector<Res> vecRes;

using namespace std;
// size_t is "unsigned (long) int", 0 -- 4 294 967 295
// int is the same as size_t, while size_t is more suitable to sizeof operator
typedef unsigned int uint;
typedef unsigned short int usint;
typedef unsigned char uint8;
typedef char int8;
typedef long long int64;
typedef unsigned long long uint64;
typedef pair<int, int> ipair;
typedef pair<double, double> dpair;
typedef unsigned long int ulint;

/// Cascade models: IC, LT
// enum CascadeModel { IC, LT };