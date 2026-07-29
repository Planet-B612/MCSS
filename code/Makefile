all: run runprefetch runlog rundebug test

run: *.cpp *.h
		g++ run.cpp -Wall -mavx512f -O3 -std=c++17 -O3 ../dSFMT/dSFMT.c -o run

runprefetch: *.cpp *.h
		g++ run.cpp -Wall  -mavx512f -O3 -DPREFETCH -std=c++17 -O3 ../dSFMT/dSFMT.c -o runprefetch

runlog: *.cpp *.h
		g++ run.cpp -Wall -g -DLOG -std=c++17 -O3 ../dSFMT/dSFMT.c -o runlog

rundebug: *.cpp *.h
		g++ run.cpp -Wall -g -DDEBUG -DLOG -std=c++17 -O3 ../dSFMT/dSFMT.c -o rundebug

test: *.cpp *.h
		g++ test_vanilla.cpp -Wall -g -std=c++17 -O3 ../dSFMT/dSFMT.c -o test

Gene_userInfo: *.cpp *.h
		g++ Gene_userInfo.cpp -Wall -std=c++17 -O3 ../dSFMT/dSFMT.c -o Gene_userInfo