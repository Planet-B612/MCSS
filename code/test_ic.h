#pragma once
#include "Argument.h"


void test_mRR()
{
	gene_syn_mRR();
	vint del_nodes = {16,19};
	for (const auto &node : del_nodes)
	{
		__Activated[node] = true;
	}
	vv_virtual_roots[0] = {10};
	vecRoot_num[0] = 3;
	add_root_lt(0, 1);
	// delete_root(0, 1);
	mRR_update_lt(0, del_nodes);
	add_root_lt(0, 1);
	cout<<"test mRR finished."<<endl;
}

void gene_syn_mRR()
{
	// __Activated[17] = true;
	mRRset mRR;
	mRRset mRR_layer;
	if(model=="IC")
	{
		mRR.push_back({5, 4, 21, 0, 6, 22, 18, 15, 13});
		mRR_layer.push_back({0,1,3,6});
		mRR.push_back({20, 12, 19, 10, 14, 1, 16, 3, 7});
		mRR_layer.push_back({0,1,3,6,8});
		vec_mRR_layer.push_back(mRR_layer);
	}
	else
	{
		mRR.push_back({11,15,6,0,18,19,2});
		mRR.push_back({22,23,3,16,9,15,13,12,10});
	}
	for (const auto &RR : mRR)
	{
		for (const auto &node : RR)
		{
			_FRsets[node].push_back(0);
		}
	}
	_mRRsets.push_back(mRR);
	vv_virtual_roots.resize(1);
	vv_polluted_nodes.resize(1);
	vecRoot_num = {2};
}

bool synthetic_check(int mRRid, string str, int newtree, int tree, int visitBool, int seq, int fr, bool reverse = false, bool identical_ele=false)
{
	bool a = false, b = false, c = false, d = false, e = false, f = false, g=false;
	if (newtree)
		a = vec_value_check(__vecNewTree, -1, 1, str + " __vecNewTree includes non -1 values.");
	if (tree)
		b = vec_value_check(__vecTree, -1, 1, str + " __vecTree includes non -1 values.");
	if (visitBool)
		c = vec_value_check(__vecVisitBool, 0, 1, str + " __vecVisitBool includes TRUE values.");
	// if (seq)
	// 	d = vec_value_check(__vecSeq, -1, 1, str + " __vecSeq includes non -1 values.");
	if (fr)
		e = FR_check(mRRid, str + " FR_check");
	if(reverse) f=FR_reverse_check(str+" FR_reverse_check");
	if(identical_ele)
	{
		g=identical_element_check(mRRid, str);
	}
	if (a || b || c || d || e || f || g)
	{
		cout << "Error in synthetic_check of " << str << endl;
		return true;
	}
	else
	{
		return false;
	}
	// for(const auto &RR:_mRRsets[mRRid])
	// {
	// 	if(RR.empty())
	// 	{
	// 		cout<<"Error in synthetic_check of "<<str<<", mRR "<<mRRid<<" contains an empty RR."<<endl;
	// 		return true;
	// 	}
	// }
	// return false;
}
bool identical_element_check(int rid, const string& str)
{
    mRRset& mRR = _mRRsets[rid];
    mRRset& mRR_layer = vec_mRR_layer[rid];

    int pre_node = -1;

    for (const auto& RR : mRR)
    {
        const int* p = RR.data();
        size_t n = RR.size();

        if (n == 0) {
            continue;
        }

        /*
            原始逻辑中 pre_node 是跨 RR 保留的：

                for RR in mRR:
                    for node in RR:
                        if pre_node == node:
                            error
                        pre_node = node

            所以需要先检查：
                上一个 RR 的最后一个节点 == 当前 RR 的第一个节点
        */
        if (pre_node == p[0])
        {
            cout << "pre_node " << pre_node
                 << " == node " << p[0] << endl;
            return true;
        }

        int bad_pre = -1;
        int bad_cur = -1;

        if (check_adjacent_equal_avx512(p, n, bad_pre, bad_cur))
        {
            cout << "pre_node " << bad_pre
                 << " == node " << bad_cur << endl;
            return true;
        }

        pre_node = p[n - 1];
    }

    if (model == "IC")
    {
        for (const auto& layer : mRR_layer)
        {
            const int* p = layer.data();
            size_t n = layer.size();

            int bad_pre = -1;
            int bad_cur = -1;

            if (check_strictly_increasing_avx512(p, n, bad_pre, bad_cur))
            {
                cout << "pre_idx " << bad_pre
                     << " <= idx " << bad_cur << endl;
                return true;
            }
        }
    }

    return false;
}


bool v_roots_check(int mRRid, string str)
{
	for (const auto &root : vv_virtual_roots[mRRid])
	{
		if (__vecTree[root] < 0 && __vecNewTree[root] < 0 && __vecVisitBool[root] == false)
		{
			cout << "Error in v_roots_check of " << str << ", the root " << root << " is not in mRR " << mRRid << " or new mRR." << endl;
			return true;
		}
	}
	return false;
}

bool del_nodes_check(int mRRid, vint &del_nodes)
{
	int i = 0;
	for (const auto &node : del_nodes)
	{
		if (__vecTree[node] < 0)
		{
			i++;
		}
	}
	if (i > 0)
	{
		out_mRRset(_mRRsets[mRRid], mRRid);
		cout << "Error in del_nodes_check of mRR " << mRRid << ", at least one del_node is not in mRR or new mRR." << endl;
		vec_out(del_nodes, "del_nodes: ");
		return true;
	}
	return false;
}

template <typename T, typename T1>
bool vec_value_check(T &vec, T1 val, int equality, string str) // equality: 1: should be equal to val, -1: shoud not equal to val, 2: should be greater than, -2: should be smaller than
{
	bool flag = false;
	vint vec_ind={};
	if (equality == 1)
	{
		size_t n = vec.size(), i=0;
		__m512i v_val = _mm512_set1_epi32(val);
		for(;i+16<n;i+=16)
		{
			__m512i v_data = _mm512_load_epi32((const __m512i*)(vec.data() + i));
        	__mmask16 eq_mask = _mm512_cmpeq_epi32_mask(v_data, v_val);
        	__mmask16 neq_mask = _mm512_knot(eq_mask);
			if (neq_mask == 0) continue;
			flag = true;
			__m512i v_base = _mm512_set1_epi32((int)i);
        	__m512i v_indices = _mm512_add_epi32(
            _mm512_set_epi32(15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0), v_base);
        
        // 用掩码压缩存储：只存储不等于 val 的索引
        alignas(64) int temp_indices[16];
        _mm512_mask_compressstoreu_epi32(temp_indices, neq_mask, v_indices);
        
        // 计算实际存储的数量
        int count = __builtin_popcount((unsigned)neq_mask);
        vec_ind.insert(vec_ind.end(), temp_indices, temp_indices + count);
		}
		for (ulint i = 0; i < vec.size(); i++)
		{
			if (vec[i] != val)
			{
				flag = true;
				vec_ind.push_back(i);
			}
		}
		if (flag)
		{
			cout << str<<VAR_NAME(vec) << ", vec Value errors: " << endl;
			for (auto i : vec_ind)
			{
				cout << vec[i] << ", ";
			}
			vec_out(vec_ind, "vec_ind: ");
			return true;
		}
		return false;
	}
	else if (equality == -1)
	{
		for (auto i : vec)
		{
			if (i == val)
			{
				cout << str << " vec Value error: " << i << " = " << val << endl;
				vec_out(vec);
				// exit(0);
				return true;
			}
		}
		return false;
	}
	else if (equality == 2)
	{
		for (auto i : vec)
		{
			if (i <= val)
			{
				cout << str << " vec Value error: " << i << " != " << val << endl;
				vec_out(vec);
				// exit(0);
				return true;
			}
		}
		return false;
	}
	else if (equality == -2)
	{
		for (auto i : vec)
		{
			if (i >= val)
			{
				cout << str << " vec Value error: " << i << " != " << val << endl;
				vec_out(vec);
				// exit(0);
				return true;
			}
		}
		return false;
	}
	return false;
}

bool FR_reverse_check_hash(int rid, mRRset &mRR_copy, string str)
{
	// sint &mRR_hash = vec_hash_mRR[rid];
	// for (const auto &node : mRR_hash)
	// {
	// 	if (vec_hash_FR[node].find(rid) == vec_hash_FR[node].end())
	// 	{
	// 		cout << "FR_check error: " << rid << " is not in the FR of " << node << endl;
	// 		return true;
	// 	}
	// }
	// for (const auto &RR : mRR_copy)
	// {
	// 	for (const auto &node : RR)
	// 	{
	// 		if (vec_hash_FR[node].find(rid) == vec_hash_FR[node].end() && __vecNewTree[node] > -1)
	// 		{
	// 			cout << "FR_check error: " << rid << " is not in the FR of " << node << endl;
	// 			return true;
	// 		}
	// 	}
	// }
	return false;
}

bool FR_reverse_check(string str)
{
	for (int i = 0; i < __numV; i++)
	{
		bool find_it = false;
		if (__Activated[i])
			continue;
		auto &frset = _FRsets[i];
		for (const auto &rid : frset)
		{
			__m512i i_512 = _mm512_set1_epi32(i);
			for (const auto &RR : _mRRsets[rid])
			{
				ulint j=0, RR_size=RR.size();
				for(;j+16<=RR_size;j+=16)
				{
					__m512i v_data = _mm512_load_epi32((const void*)(RR.data() + j));
        			__mmask16 mask = _mm512_cmpeq_epi32_mask(v_data, i_512);
					if (mask != 0)
					{
						find_it=true;
						break;
					}
				}
				for (;j<RR_size;j++)
				{
					if (RR[j] == i)
					{
						find_it = true;
						break;
					}
				}
			}
			if (!find_it)
			{
				cout << str + " Error in FR_reverse_check, the node " << i << " is not in any RR of mRR " << rid << "; but the node's _FRsets contains the mRRid" << endl;
				out_vec(_FRsets[i]);
				out_layer(vec_mRR_layer[rid], rid);
				out_mRRset(_mRRsets[rid], rid);
				return true;
			}
		}
	}
	return false;
}

bool FR_check_hash(int rid, string str)
{
	sint &mRR_hash = vec_hash_mRR[rid];
	bool flag = false;
	for (const auto &node : mRR_hash)
	{
		if (vec_hash_FR[node].find(rid) == vec_hash_FR[node].end())
		{
			// cout<<"RRid in vec_hash_mRR of "<<node <<" includes: ";
			// for(const auto &id:vec_hash_FR[node])
			// {
			// 	cout<<id<<", ";
			// }
			// cout<<endl;
			flag = true;
			cout << "FR_check error: " << rid << " is not in the FR of " << node << endl;
			cout<<"RRid in vec_hash_mRR of "<<node <<" includes: ";
			for(const auto index:vec_hash_FR[node])
			{
				cout<<index<<", "<<endl;
			}
			return true;
		}
	}
	if(flag)
	{
		return true;
	}
	return false;
}

bool FR_check(int rid, const string& str)
{
    mRRset& mRR = _mRRsets[rid];

    std::atomic<bool> found_error{false};
    std::atomic<int> error_node{-1};

    size_t total_size = 0;
    for (const auto& RR : mRR) {
        total_size += RR.size();
    }

    std::vector<int> all_nodes;
    all_nodes.reserve(total_size);

    for (const auto& RR : mRR) {
        all_nodes.insert(all_nodes.end(), RR.begin(), RR.end());
    }

    #pragma omp parallel for
    for (size_t idx = 0; idx < all_nodes.size(); ++idx)
    {
        if (found_error.load(std::memory_order_relaxed)) {
            continue;
        }

        int node = all_nodes[idx];
        const auto& frset = _FRsets[node];

        auto it = std::lower_bound(frset.begin(), frset.end(), rid);

        if (it == frset.end() || *it != rid)
        {
            bool expected = false;

            if (found_error.compare_exchange_strong(
                    expected,
                    true,
                    std::memory_order_relaxed))
            {
                error_node.store(node, std::memory_order_relaxed);
            }
        }
    }

    if (found_error.load(std::memory_order_relaxed))
    {
        int node = error_node.load(std::memory_order_relaxed);
		out_mRRset(_mRRsets[rid], rid);
		out_vec(_FRsets[node]);
        cout << str << " Error in FR_full_check, mRR " << rid
             << " contains the node " << node
             << "; but the mRRid is not in node's _FRsets\n";

        return true;
    }

    return false;
}

	// for (auto &RR : mRR)
	// {
	// 	for (auto &node : RR)
	// 	{
	// 		auto &frset = _FRsets[node];
	// 		auto it = lower_bound(frset.begin(), frset.end(), rid);
	// 		if (it == frset.end())
	// 		{
	// 			set_out({node});
	// 			output_info(rid, false);
	// 			cout << str + " Error in FR_full_check, mRR " << rid << " contains the node " << node << "; but the mRRid is not in node's _FRsets" << endl;
	// 			// output_info(rid);
	// 			return true;
	// 		}
	// 	}
	// }
	// for(int i=0;i<__numV;i++)
	// {
	// 	// if(_FRsets[i].find(rid)!=_FRsets[i].end())
	// 	auto &frset= _FRsets[i];
	// 	auto it= lower_bound(frset.begin(), frset.end(), rid);
	// 	if( it != frset.end() )
	// 	{
	// 		bool find_it=false;
	// 		for(auto &adj_list:mRR)
	// 		{
	// 			// if(adj_list.find(i)==adj_list.end())
	// 			// {
	// 			// 	continue;
	// 			// }
	// 			// else
	// 			// {
	// 			// 	find_it=true;
	// 			// 	break;
	// 			// }
	// 		}
	// 		if(!find_it)
	// 		{
	// 			set_out({i});
	// 			// output_info(rid);
	// 			cout<<"Error in FR_full_check, _FRset[node] contains mRRid "<<rid<<", but the node "<<i<<" is not in any adj_list of the mRR "<<endl;
	// 			output_info(rid);
	// 			return i;
	// 		}
	// 	}
	// }
// 	return false;
// }

bool FR_sorted_check(string str)
{
	int flag=0;
	for (ulint j=0;j<_FRsets.size();j++)
	{
		auto &vec=_FRsets[j];
		size_t n = vec.size(), i=1;
    	if (n <= 1) return false;
		for (; i + 16 <= n; i += 16) 
		{
			__m512i v_prev = _mm512_load_epi32((const __m512i*)(vec.data() + i - 1));
			__m512i v_curr = _mm512_load_epi32((const __m512i*)(vec.data() + i));
			__mmask16 mask = _mm512_cmpgt_epi32_mask(v_prev, v_curr);
			if (mask != 0) 
			{
				flag=1;
				break;
			}
    	}
		for (; i < n; ++i) 
		{
			if (vec[i - 1] > vec[i]) 
			{
				flag=1;
				break;
			}
		}
		if (flag==1)
		{
			cout << "Error in FR_sorted_check of " + str + " in checking the " + to_string(j) + "-th FRset. The _FRsets is not sorted." << endl;
			exit(1);
		}
	}
	return false;
}

bool FR_insert_check(int rid, int node)
{
	// if(_FRsets[node].find(rid)!=_FRsets[node].end())
	auto &frset = _FRsets[node];
	auto it = lower_bound(frset.begin(), frset.end(), rid);
	if (it != frset.end() && *it == rid)
	{
		set_out({node});
		// __log_message("INFO", __FILE__, __LINE__, __func__);
		output_info(rid);
		// exit(1);
		return 1;
	}
	return 0;
}

template <typename T>
void vec_out(T &vec, string str = "")
{
	std::fstream result_bk(result, ios::app);
	assert(!result_bk.fail());
	result_bk << str + "vec: ";
	for (auto i : vec)
	{
		result_bk << i << ", ";
	}
	result_bk << endl;
	result_bk.close();
	// cout<<str+"vec: ";
	// for(auto i:vec)
	// {
	// 	cout<<i<<", ";
	// }
	// cout<<endl;
}

void mRR_out(int mRRid, string str = "")
{
	std::fstream result_bk(result, ios::app);
	assert(!result_bk.fail());
	result_bk << str + "mRRset: " << mRRid << endl;
	auto &mRR = _mRRsets[mRRid];
	for (ulint i = 0; i < mRR.size(); i++)
	{
		result_bk << i << "-th RR: ";
		for (auto &j : mRR[i])
		{
			result_bk << j << ", ";
		}
		result_bk << endl;
	}
	result_bk.close();
}

void set_out(vint p_nodes)
{
	// std::fstream result_bk(result, ios::out);
	(*__arg).result_bk.open(result);
	assert(!(*__arg).result_bk.fail());
	for (auto i : p_nodes)
	{
		(*__arg).result_bk << i << "'s hashset values: " << endl;
		for (auto j : _FRsets[i])
		{
			(*__arg).result_bk << j << ", ";
		}
		(*__arg).result_bk << endl;
	}
	(*__arg).result_bk.close();
}

bool output_info(int mRRid, bool erase = false, vint p_nodes = {})
{
	// cout<<vec_virtual_roots[mRRid].size()<<endl;
	// if(p_nodes.size()==0 && vec_virtual_roots[mRRid].size()==0)
	// {
	// 	return true;
	// }
	if (erase)
		(*__arg).result_bk.open(result);
	else
		(*__arg).result_bk.open(result, ios::app);
	assert(!(*__arg).result_bk.fail());
	(*__arg).result_bk << "===========================================================" << endl;
	(*__arg).result_bk << "The vec_virtual_roots is: " << endl;
	for (auto i : vv_virtual_roots[mRRid])
	{
		(*__arg).result_bk << i << ", ";
	}
	(*__arg).result_bk << endl;
	(*__arg).result_bk << "The p_nodes are: " << endl;
	for (auto node : p_nodes)
	{
		(*__arg).result_bk << node << ", ";
	}
	(*__arg).result_bk << endl;
	(*__arg).result_bk << "The mRRid is: " << mRRid << endl;
	for (auto i = 0; i < __numV; i++)
	{
		(*__arg).result_bk << __vecVisitBool[i] << ", ";
	}
	(*__arg).result_bk << endl;
	(*__arg).result_bk << "The trees are: " << endl;
	auto k = 0;
	for (auto &RR : _mRRsets[mRRid])
	{
		// if(k==0) 		continue;
		(*__arg).result_bk << k << "-th tree is: " << endl;
		for (auto &node : RR)
		{
			(*__arg).result_bk << node << ", ";
			(*__arg).result_bk << endl;
		}
		k++;
	}
	(*__arg).result_bk.close();
	return false;
}


static inline bool check_adjacent_equal_avx512(
    const int* p,
    size_t n,
    int& pre_val,
    int& cur_val
)
{
    if (n < 2) {
        return false;
    }

    size_t i = 0;

    /*
        比较：
        p[i + 0]  == p[i + 1]
        p[i + 1]  == p[i + 2]
        ...
        p[i + 15] == p[i + 16]

        所以每轮需要访问到 p[i + 16]。
        条件是 i + 16 < n。
    */
    for (; i + 16 < n; i += 16)
    {
        __m512i cur  = _mm512_load_si512((const void*)(p + i));
        __m512i next = _mm512_loadu_si512((const void*)(p + i + 1));

        __mmask16 mask = _mm512_cmpeq_epi32_mask(cur, next);

        if (mask != 0)
        {
            int lane = __builtin_ctz(mask);

            pre_val = p[i + lane];
            cur_val = p[i + lane + 1];

            return true;
        }
    }

    // tail scalar
    for (; i + 1 < n; ++i)
    {
        if (p[i] == p[i + 1])
        {
            pre_val = p[i];
            cur_val = p[i + 1];
            return true;
        }
    }

    return false;
}


static inline bool check_strictly_increasing_avx512(
    const int* p,
    size_t n,
    int& pre_val,
    int& cur_val
)
{
    if (n == 0) {
        return false;
    }

    /*
        原始代码里：

            int pre_idx = -1;
            for (idx : layer) {
                if (idx <= pre_idx) error;
                pre_idx = idx;
            }

        所以第一个元素也要和 -1 比较。
        如果 p[0] <= -1，需要报错。
    */
    if (p[0] <= -1)
    {
        pre_val = -1;
        cur_val = p[0];
        return true;
    }

    if (n < 2) {
        return false;
    }

    size_t i = 0;

    /*
        检查：
        p[i + 1]  <= p[i + 0]
        p[i + 2]  <= p[i + 1]
        ...
        p[i + 16] <= p[i + 15]

        即 next <= cur 时报错。
    */
    for (; i + 16 < n; i += 16)
    {
        __m512i cur  = _mm512_load_si512((const void*)(p + i));
        __m512i next = _mm512_loadu_si512((const void*)(p + i + 1));

        /*
            AVX-512 没有直接的 <= signed epi32 mask intrinsic，
            可以用：

                next <= cur

            等价于：

                cur >= next

            使用 _mm512_cmpge_epi32_mask(cur, next)
        */
        __mmask16 mask = _mm512_cmpge_epi32_mask(cur, next);

        if (mask != 0)
        {
            int lane = __builtin_ctz(mask);

            pre_val = p[i + lane];
            cur_val = p[i + lane + 1];

            return true;
        }
    }

    // tail scalar
    for (; i + 1 < n; ++i)
    {
        if (p[i + 1] <= p[i])
        {
            pre_val = p[i];
            cur_val = p[i + 1];
            return true;
        }
    }

    return false;
}
