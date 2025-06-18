#pragma once
#include <vector>
#include <set>



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

typedef std::set<int> sint;
typedef std::vector<sint> vsint;

// /// Define the way of storing results
// typedef std::tuple<double, double, double, double> Res;
// typedef std::vector<Res> vecRes;

// size_t is "unsigned (long) int", 0 -- 4 294 967 295
// int is the same as size_t, while size_t is more suitable to sizeof operator
typedef unsigned int uint;
typedef unsigned short int usint;
typedef unsigned char uint8;
typedef char int8;
typedef long long int64;
typedef unsigned long long uint64;
typedef std::pair<int, int> ipair;
typedef std::pair<double, double> dpair;
typedef unsigned long int ulint;

/// Cascade models: IC, LT
// enum CascadeModel { IC, LT };