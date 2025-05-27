all: run test

run: *.cpp *.h
		g++ run.cpp -Wall -DNDEBUG -std=c++11 -O3 ../dSFMT/dSFMT.c -o run

# run: *.cpp *.h
# 		g++ run.cpp -Wall -std=c++11 -O3 ../dSFMT/dSFMT.c -o run

test: *.cpp *.h
		g++ test_vanilla.cpp -Wall -std=c++11 -O3 ../dSFMT/dSFMT.c -o test

# Gene_userInfo: *.cpp *.h
# 		g++ Gene_userInfo.cpp -Wall -std=c++17 -O3 ../dSFMT/dSFMT.c -o Gene_userInfo