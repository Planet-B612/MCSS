// Minimal repro for calling delete_root() mid mRR_update_and_add_roots.
//
// Build & run (ASan recommended):
//   make repro_delete_root && ./repro_delete_root
//
// Or through the full class path (same bug, needs no real graph file):
//   make test_asan && ./test -repro_delete_root

#define DEBUG

#include <atomic>
#include <mutex>
#include <omp.h>
#include <immintrin.h>
#include <iostream>
#include <vector>
#include <chrono>
#include <cstring>
#include <ctime>
#include <cassert>
#include "../dSFMT/dSFMT.h"
#include "Argument.h"
#include "mRRcollection.h"

using namespace std;

// Same structural mistake as delete_root mid-update: keep a reference into
// mRR[min_tree], swap later trees away, then resize mRR from the back until
// min_tree is destroyed, then push_back through the dangling reference.
static void repro_pattern_only()
{
	cout << "--- pattern-only repro (vector resize + dangling ref) ---" << endl;
	vector<vector<int>> mRR = {
		{0, 1, 2, 3},
		{10, 11, 12, 13, 14}, // min_tree
		{20, 21, 22},
		{30, 31, 32},
	};
	const int min_tree = 1;
	vector<vector<int>> mRR_copy(mRR.size() - min_tree);

	// Move suffix of min_tree into copy[0] (as in the update), then swap later trees.
	mRR_copy[0].assign(mRR[min_tree].begin() + 3, mRR[min_tree].end());
	mRR[min_tree].resize(3);
	for (int i = min_tree + 1; i < static_cast<int>(mRR.size()); ++i)
		mRR[i].swap(mRR_copy[i - min_tree]);

	auto &min_tree_RR = mRR[min_tree];

	// delete_root(mRRid, 3) with empty v_roots: drop last 3 trees, including min_tree.
	const int num_del = 3;
	mRR.resize(mRR.size() - num_del);

	// BFS would do this next — use-after-free if ASan is on.
	min_tree_RR.push_back(99);
	cout << "UNEXPECTED: pattern-only finished (no ASan?)" << endl;
}

int main(int argc, char **argv)
{
	dsfmt_gv_init_gen_rand(42);

	bool pattern_only = false;
	for (int i = 1; i < argc; ++i)
	{
		if (argv[i] == string("-pattern_only"))
			pattern_only = true;
	}

	if (pattern_only)
	{
		repro_pattern_only();
		return 0;
	}

	cout << "--- full-path repro via mRR_update_and_add_roots ---" << endl;
	const int n = 40;
	Argument arg;
	arg.model = "IC";
	arg.numV = n;
	arg.real_time_pw = true; // avoid loading a large PO into a tiny PO[]

	R_graph.assign(n, vint_aligned());
	O_graph.assign(n, vint_aligned());
	Inv_inDeg.assign(n, 1.0f);
	__Activated.assign(n, 0);
	numV512 = _mm512_set1_epi32(n);
	root_num = 1;

	mRRcollection RR(arg);
	RR.repro_delete_root_mid_update();
	return 0;
}
