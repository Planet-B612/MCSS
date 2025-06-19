all: run rundebug test

run: *.cpp *.h
		g++ run.cpp -Wall -pg -std=c++11 -O3 ../dSFMT/dSFMT.c -o run

rundebug: *.cpp *.h
		g++ run.cpp -Wall -g -DDEBUG -std=c++11 -O3 ../dSFMT/dSFMT.c -o rundebug

test: *.cpp *.h
		g++ test.cpp -Wall -g -std=c++11 -O3 ../dSFMT/dSFMT.c -o test

# Gene_userInfo: *.cpp *.h
# 		g++ Gene_userInfo.cpp -Wall -std=c++17 -O3 ../dSFMT/dSFMT.c -o Gene_userInfo