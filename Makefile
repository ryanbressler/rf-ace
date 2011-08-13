COMPILER = g++
CFLAGS = -O3 -Wall -pedantic -I/usr/lib64/glib-2.12/include -I/usr/include/glib-2.12 -I/usr/ -Isrc/
SOURCEFILES = src/randomforest.cpp src/GBT.cpp src/rootnode.cpp src/node.cpp src/treedata.cpp src/mtrand.cpp src/datadefs.cpp
STATICFLAGS = -static-libgcc -static
TESTFILES = test/argparse_test.hpp test/datadefs_test.hpp
TESTFLAGS = -lcppunit -pedantic -I/usr/lib64/glib-2.12/include -I/usr/include/glib-2.12 -I/usr/ -Isrc/
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
