all: run runprefetch runlog rundebug test

run: *.cpp *.h
		g++ run.cpp -Wall -mavx512f -std=c++17 -O3 ../dSFMT/dSFMT.c -o run

runprefetch: *.cpp *.h
		g++ run.cpp -Wall  -mavx512f -DPREFETCH -std=c++17 -O3 ../dSFMT/dSFMT.c -o runprefetch

runlog: *.cpp *.h
		g++ run.cpp -Wall -g -DLOG -std=c++17 -O3 ../dSFMT/dSFMT.c -o runlog

rundebug: *.cpp *.h
		g++ run.cpp -Wall -g -DDEBUG -DLOG -std=c++17 -O3 ../dSFMT/dSFMT.c -o rundebug

test: *.cpp *.h
		g++ test.cpp -Wall -g -mavx512f -std=c++17 -O3 ../dSFMT/dSFMT.c -o test

Gene_userInfo: *.cpp *.h
		g++ Gene_userInfo.cpp -Wall -std=c++17 -O3 ../dSFMT/dSFMT.c -o Gene_userInfo