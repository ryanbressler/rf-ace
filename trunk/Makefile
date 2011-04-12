COMPILER = g++
CFLAGS = -g -pg -O3 -Wall -pedantic
RF_ACE_LIBS = -lboost_program_options

all: rf_ace benchmark GBT_test

rf_ace: src/rf_ace.cpp src/randomforest.cpp src/node.cpp src/treedata.cpp src/mtrand.cpp src/datadefs.cpp
	$(COMPILER) $(CFLAGS) src/rf_ace.cpp src/randomforest.cpp src/node.cpp src/treedata.cpp src/mtrand.cpp src/datadefs.cpp -o bin/rf_ace $(RF_ACE_LIBS)

benchmark: src/benchmark.cpp src/randomforest.cpp src/node.cpp src/treedata.cpp src/datadefs.cpp 
	$(COMPILER) $(CFLAGS) src/benchmark.cpp src/randomforest.cpp src/node.cpp src/treedata.cpp src/datadefs.cpp -o bin/benchmark

GBT_test: src/GBT_test.cpp src/GBT.cpp src/node.cpp src/treedata.cpp src/datadefs.cpp
	$(COMPILER) $(CFLAGS) src/GBT_test.cpp src/GBT.cpp src/node.cpp src/treedata.cpp src/datadefs.cpp -o bin/GBT_test

clean:
	rm -rf bin/rf_ace bin/benchmark bin/GBT_test src/*.o
