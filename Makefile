COMPILER = g++
CFLAGS = -O3 -Wall -pedantic
MPFLAG = -fopenmp

all: rf_ace

rf_ace: src/rf_ace.cpp src/getopt_pp.cpp src/randomforest.cpp src/GBT.cpp src/rootnode.cpp src/node.cpp src/treedata.cpp src/mtrand.cpp src/datadefs.cpp
	$(COMPILER) $(CFLAGS) src/rf_ace.cpp src/getopt_pp.cpp src/randomforest.cpp src/GBT.cpp src/rootnode.cpp src/node.cpp src/treedata.cpp src/mtrand.cpp src/datadefs.cpp -o bin/rf_ace



benchmark: src/benchmark.cpp src/randomforest.cpp src/rootnode.cpp src/node.cpp src/treedata.cpp src/datadefs.cpp 
	$(COMPILER) $(CFLAGS) src/benchmark.cpp src/randomforest.cpp src/rootnode.cpp src/node.cpp src/treedata.cpp src/datadefs.cpp -o bin/benchmark

GBT_test: src/GBT_test.cpp src/getopt_pp.cpp src/GBT.cpp src/rootnode.cpp src/node.cpp src/treedata.cpp src/mtrand.cpp src/datadefs.cpp
	$(COMPILER) $(CFLAGS) src/GBT_test.cpp src/getopt_pp.cpp src/GBT.cpp src/rootnode.cpp src/node.cpp src/treedata.cpp src/mtrand.cpp src/datadefs.cpp -o bin/GBT_test

clean:
	rm -rf bin/rf_ace bin/benchmark bin/GBT_test src/*.o
