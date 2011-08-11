debug:
	clang++ -g -Wall -Wextra -O3 -funroll-loops main.cpp -I$(HOME)/usr/include/eigen3

release:
	g++ -Wall -Wextra -O3 -funroll-loops -DNDEBUG main.cpp -I$(HOME)/usr/include/eigen3
