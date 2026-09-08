all: run runprefetch runlog rundebug test

run: *.cpp *.h
		g++ run.cpp -Wall -g -mavx512f -std=c++17 -O3 ../dSFMT/dSFMT.c -o run

runprefetch: *.cpp *.h
		g++ run.cpp -Wall -g  -mavx512f -DPREFETCH -std=c++17 -O3 ../dSFMT/dSFMT.c -o runprefetch

runlog: *.cpp *.h
		g++ run.cpp -Wall -g -DLOG -std=c++17 -O3 ../dSFMT/dSFMT.c -o runlog

rundebug: *.cpp *.h
		g++ run.cpp -Wall -g -DDEBUG -DLOG -std=c++17 -O3 ../dSFMT/dSFMT.c -o rundebug

test: *.cpp *.h
		g++ test.cpp -Wall -g -mavx512f -fopenmp -std=c++17 -O3 ../dSFMT/dSFMT.c -o test

# ASan builds to catch UAF from delete_root mid mRR_update_and_add_roots
test_asan: test.cpp *.h
		g++ test.cpp -Wall -g -mavx512f -fopenmp -std=c++17 -O1 -fsanitize=address -fno-omit-frame-pointer ../dSFMT/dSFMT.c -o test_asan

repro_delete_root: repro_delete_root.cpp *.h
		g++ repro_delete_root.cpp -Wall -g -mavx512f -std=c++17 -O1 -fsanitize=address -fno-omit-frame-pointer ../dSFMT/dSFMT.c -o repro_delete_root

Gene_userInfo: *.cpp *.h
		g++ Gene_userInfo.cpp -Wall -std=c++17 -O3 ../dSFMT/dSFMT.c -o Gene_userInfo