#pragma once
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <queue>
#include "CommonStruc.h"
using namespace std;
#include "../dSFMT/dSFMT.h"
#include <tuple>
using std::tuple;
/// Assertion
template <typename _Ty>
static inline void ASSERT(_Ty v)
{
	if (!(v))
	{
		cerr << "ASSERT FAIL @ " << __FILE__ << ":" << __LINE__ << endl;
		exit(1);
	}
}

template <typename _Ty>
static inline void ASSERTT(_Ty v, string s)
{
	if (!(v))
	{
		cerr << s << ": ASSERT FAIL @ " << __FILE__ << ":" << __LINE__ << endl;
		exit(1);
	}
}

bool comp(tuple<int, int, int> a, tuple<int, int, int> b)
{
	return get<1>(a) < get<1>(b);
}; // sort in ascending order. the smaller, the earlier

// /// Max and Min
// template <typename _Ty>
// static inline _Ty MAXX(_Ty a, _Ty b) {return a>b ? a : b;}

// template <typename _Ty>
// static inline _Ty MINN(_Ty a, _Ty b) {return a>b ? b : a;}

/// Math, pow2
static inline double pow2(const double t)
{
	return t * t;
}

static inline double power(const double t, const uint8_t k)
{
	ASSERTT(k > 0 && k - (int)k == 0, "This power function is only suitable for positive integers");
	double val = t;
	for (uint8_t i = 0; i < k - 1; i++)
	{
		val = val * t;
	}
	return val;
}

/// Math, log2
static inline double log2(const size_t n)
{
	return log(n) / log(2);
}

/// Math, logcnk
static inline double logcnk(const size_t n, size_t k)
{
	k = k < n - k ? k : n - k;
	double res = 0;
	for (size_t i = 1; i <= k; i++)
		res += log(double(n - k + i) / i);
	return res;
}

/// Log information
template <typename _Ty>
static inline void loginfo(_Ty val)
{
	cout << val << endl;
}

/// Log information
template <typename _Ty>
static inline void loginfo(const string title, _Ty val)
{
	cout << title << ": " << val << endl;
}

/// Normalize the probabilities to a accumulative format, e.g., [0.2, 0.5, 0.3]->[0.2, 0.7, 1.0]
// static inline void to_normal_accum_prob(Graph &vecGraph)
// {
// 	for (auto &nbrs : vecGraph)
// 	{
// 		float accumVal = float(0.0);
// 		for (auto &nbr : nbrs)
// 		{
// 			accumVal += get<1>(nbr);
// 			get<1>(nbr) = accumVal;
// 		}
// 		// Normalization
// 		for (auto &nbr : nbrs)
// 		{
// 			get<1>(nbr) /= accumVal;
// 		}
// 	}
// }

/// Make the vector to a min-heap.
inline void make_min_heap(vint &vec)
{
	// Min heap
	const auto size = vec.size();
	if (2 <= size)
	{
		for (auto hole = (size + 1) / 2; hole--;)
		{
			const auto val = vec[hole];
			size_t i, child;
			for (i = hole; i * 2 + 1 < size; i = child)
			{
				// Find smaller child
				child = i * 2 + 2;
				if (child == size || vec[child - 1] < vec[child])
				{
					// One child only or the left child is smaller than the right one
					--child;
				}

				// Percolate one level
				if (vec[child] < val)
				{
					vec[i] = vec[child];
				}
				else
				{
					break;
				}
			}
			vec[i] = val;
		}
	}
}

/// Replace the value for the first element and down-heap this element.
inline void min_heap_replace_min_value(vint &vec, const int &val)
{
	// Increase the value of the first element
	const auto size = vec.size();
	size_t i, child;
	for (i = 0; i * 2 + 1 < size; i = child)
	{
		// Find smaller child
		child = i * 2 + 2;
		if (child == size || vec[child - 1] < vec[child])
		{
			// One child only or the left child is smaller than the right one
			--child;
		}

		// Percolate one level
		if (vec[child] < val)
		{
			vec[i] = vec[child];
		}
		else
		{
			break;
		}
	}
	vec[i] = val;
}

//=========================================================================================================
//=================-------------------- MAKE MAX HEAP --------------------=================================
//=========================================================================================================

/// Make the vector to a max-heap.
static inline void make_max_heap(vector<tuple<int, double, int>> &vec)
{
	// Max heap
	const auto size = vec.size();
	if (2 <= size)
	{
		for (auto hole = (size + 1) / 2; hole--;)
		{
			const auto val = vec[hole];
			size_t i, child;
			for (i = hole; i * 2 + 1 < size; i = child)
			{
				// Find smaller child
				child = i * 2 + 2;
				if (child == size || get<1>(vec[child - 1]) > get<1>(vec[child]))
				{
					// One child only or the left child is greater than the right one
					--child;
				}
				// Percolate one level
				if (get<1>(vec[child]) > get<1>(val))
				{
					vec[i] = vec[child];
				}
				else
				{
					break;
				}
			}
			vec[i] = val;
		}
	}
}

/// Replace the value for the first element and down-heap this element.
static inline void max_heap_replace_max_value(vector<tuple<int, double, int>> &vec, tuple<int, double, int> &val)
{
	// Increase the value of the first element
	const auto size = vec.size();
	size_t i, child;
	auto hole = vec[0];
	for (i = 0; i * 2 + 1 < size; i = child)
	{
		// Find smaller child
		child = i * 2 + 2;
		if (child == size || get<1>(vec[child - 1]) > get<1>(vec[child]))
		{
			// One child only or the left child is greater than the right one
			--child;
		}

		// Percolate one level
		if (get<1>(vec[child]) > get<1>(val))
		{
			vec[i] = vec[child];
		}
		else
		{
			break;
		}
	}
	// hole.first = val;
	hole = val;
	vec[i] = hole;
}

/// Make the vector to a max-heap.
static inline void make_max_heap(vector< tuple<int, int, double, int> >& vec)
{
	// Max heap
	const auto size = vec.size();
	if (2 <= size)
	{
		for (auto hole = (size + 1) / 2; hole--;)
		{
			const auto val = vec[hole];
			size_t i, child;
			for (i = hole; i * 2 + 1 < size; i = child)
			{
				// Find smaller child
				child = i * 2 + 2;
				if (child == size || get<2>(vec[child - 1]) > get<2>(vec[child]))
				{
					// One child only or the left child is greater than the right one
					--child;
				}
				// Percolate one level
				if ( get<2>(vec[child]) > get<2>(val) )
				{
					vec[i] = vec[child];
				}
				else
				{
					break;
				}
			}
			vec[i] = val;
		}
	}
}

/// Replace the value for the first element and down-heap this element.
static inline void max_heap_replace_max_value(vector< tuple<int, int, double, int> >& vec, tuple<int, int, double, int>& val)
{
	// Increase the value of the first element
	const auto size = vec.size();
	size_t i, child;
	auto hole = vec[0];
	for (i = 0; i * 2 + 1 < size; i = child)
	{
		// Find smaller child
		child = i * 2 + 2;
		if (child == size || get<2>(vec[child - 1]) > get<2>(vec[child]) )
		{
			// One child only or the left child is greater than the right one
			--child;
		}

		// Percolate one level
		if (get<2>(vec[child]) > get<2>(val) )
		{
			vec[i] = vec[child];
		}
		else
		{
			break;
		}
	}
	//hole.first = val;
	hole=val;
	vec[i] = hole;
}

static inline void make_max_heap(vector< pair<int, double> >& vec)
{
	// Max heap
	const auto size = vec.size();
	if (2 <= size)
	{
		for (auto hole = (size + 1) / 2; hole--;)
		{
			const auto val = vec[hole];
			size_t i, child;
			for (i = hole; i * 2 + 1 < size; i = child)
			{
				// Find smaller child
				child = i * 2 + 2;
				if (child == size || vec[child - 1].second > vec[child].second)
				{
					// One child only or the left child is greater than the right one
					--child;
				}
				// Percolate one level
				if ( vec[child].second > val.second )
				{
					vec[i] = vec[child];
				}
				else
				{
					break;
				}
			}
			vec[i] = val;
		}
	}
}

static inline void max_heap_replace_max_value(vector< pair<int, double> >& vec, pair<int, double>& val)
{
	// Increase the value of the first element
	const auto size = vec.size();
	size_t i, child;
	auto hole = vec[0];
	for (i = 0; i * 2 + 1 < size; i = child)
	{
		// Find smaller child
		child = i * 2 + 2;
		if (child == size || vec[child - 1].second > vec[child].second )
		{
			// One child only or the left child is greater than the right one
			--child;
		}

		// Percolate one level
		if (vec[child].second > val.second )
		{
			vec[i] = vec[child];
		}
		else
		{
			break;
		}
	}
	//hole.first = val;
	hole=val;
	vec[i] = hole;
}