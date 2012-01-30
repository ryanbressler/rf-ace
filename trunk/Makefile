COMPILER = g++
CFLAGS = -O3 -Wall -Wextra -pedantic -I/usr/lib64/glib-2.12/include -I/usr/include/glib-2.12 -I/usr/ -Isrc/
SOURCEFILES = src/progress.cpp src/statistics.cpp src/math.cpp src/stochasticforest.cpp src/rootnode.cpp src/node.cpp src/splitter.cpp src/treedata.cpp src/mtrand.cpp src/datadefs.cpp src/utils.cpp src/partitionsequence.cpp
STATICFLAGS = -static-libgcc -static
TESTFILES = test/argparse_test.hpp test/datadefs_test.hpp test/stochasticforest_test.hpp test/rootnode_test.hpp test/node_test.hpp test/treedata_test.hpp test/partitionsequence_test.hpp test/splitter_test.hpp
TESTFLAGS = -L/home/erkkila2/lib -lcppunit -ldl -pedantic -I/home/erkkila2/include -I/usr/lib64/glib-2.12/include -I/usr/include/glib-2.12 -I/usr/ -Isrc/
MPFLAG = -fopenmp
.PHONY: all test clean  # Squash directory checks for the usual suspects

all: rf_ace

rf_ace: $(SOURCEFILES)
	$(COMPILER) $(CFLAGS) src/rf_ace.cpp $(SOURCEFILES) -o bin/rf_ace

debug: $(SOURCEFILES)
	$(COMPILER) $(CFLAGS) src/rf_ace.cpp $(SOURCEFILES) -o bin/rf_ace -ggdb

static: $(SOURCEFILES)
	$(COMPILER) $(CFLAGS) src/rf_ace.cpp $(STATICFLAGS) $(SOURCEFILES) -o bin/rf_ace

#benchmark: src/benchmark.cpp $(SOURCEFILES)
#	$(COMPILER) $(CFLAGS) src/benchmark.cpp $(SOURCEFILES) -o bin/benchmark

GBT_benchmark: test/GBT_benchmark.cpp $(SOURCEFILES)
	$(COMPILER) $(CFLAGS) test/GBT_benchmark.cpp $(SOURCEFILES) -o bin/GBT_benchmark

test: $(SOURCEFILES) $(TESTFILES)
	rm -f bin/test; $(COMPILER) $(TESTFLAGS) test/run_tests.cpp $(SOURCEFILES) -o bin/test -ggdb; ./bin/test

clean:
	rm -rf bin/rf_ace bin/benchmark bin/GBT_benchmark bin/test bin/*.dSYM/ src/*.o
