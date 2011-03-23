COMPILER = g++
CFLAGS = -g -pg -O3 -Wall -pedantic
RF_ACE_LIBS = -lboost_program_options

all: rf_ace benchmark

rf_ace: src/rf_ace.cpp src/randomforest.cpp src/node.cpp src/treedata.cpp src/datadefs.cpp
	$(COMPILER) $(CFLAGS) $(RF_ACE_LIBS) src/rf_ace.cpp src/randomforest.cpp src/node.cpp src/treedata.cpp src/datadefs.cpp -o bin/rf_ace

benchmark: src/benchmark.cpp 
	$(COMPILER) $(CFLAGS) src/benchmark.cpp -o bin/benchmark

clean:
	rm -rf bin/rf_ace bin/benchmark src/*.o
