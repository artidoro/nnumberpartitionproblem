all: kk

%: %.cc
	g++ -std=c++11 $< -o $@