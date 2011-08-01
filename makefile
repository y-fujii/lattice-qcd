debug:
	g++ -g -Wall -Wextra -O3 -funroll-loops main.cpp

release:
	g++45 -Wall -Wextra -O3 -funroll-loops -DNDEBUG main.cpp
