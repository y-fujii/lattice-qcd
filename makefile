debug:
	clang++ -g -Wall -Wextra -O3 -funroll-loops main.cpp

release:
	g++ -Wall -Wextra -O3 -funroll-loops -DNDEBUG main.cpp
