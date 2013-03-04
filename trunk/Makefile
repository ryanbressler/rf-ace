COMPILER = g++
CFLAGS = -O3 -std=c++0x -Wall -Wextra -pedantic -Isrc/ -lz
TFLAGS = -pthread
SOURCEFILES = src/murmurhash3.cpp src/datadefs.cpp src/progress.cpp src/statistics.cpp src/math.cpp src/stochasticforest.cpp src/rootnode.cpp src/node.cpp src/treedata.cpp src/utils.cpp src/distributions.cpp src/reader.cpp src/feature.cpp
STATICFLAGS = -static-libgcc -static
TESTFILES = test/rface_test.hpp test/distributions_test.hpp test/argparse_test.hpp test/datadefs_test.hpp test/stochasticforest_test.hpp test/utils_test.hpp test/math_test.hpp test/rootnode_test.hpp test/node_test.hpp test/treedata_test.hpp
NEWTESTFILES = test/reader_newtest.hpp
TESTFLAGS = -std=c++0x -L${HOME}/lib/ -L/usr/local/lib -lcppunit -ldl -pedantic -I${HOME}/include/ -I/usr/local/include -Itest/ -Isrc/
.PHONY: all test clean  # Squash directory checks for the usual suspects

all: rf-ace

rf-ace: $(SOURCEFILES)
	$(COMPILER) $(CFLAGS) src/rf_ace.cpp $(SOURCEFILES) $(TFLAGS) -o bin/rf-ace

rf-ace-i386: $(SOURCEFILES)
	$(COMPILER) $(CFLAGS) -m32 src/rf_ace.cpp $(SOURCEFILES) $(TFLAGS) -o bin/rf-ace-i386

rf-ace-amd64: $(SOURCEFILES)
	$(COMPILER) $(CFLAGS) -m64 src/rf_ace.cpp $(SOURCEFILES) $(TFLAGS) -o bin/rf-ace-amd64

no-threads: $(SOURCEFILES)
	$(COMPILER) $(CFLAGS) -DNOTHREADS $(SOURCEFILES) src/rf_ace.cpp -o bin/rf-ace

debug: $(SOURCEFILES)
	$(COMPILER) $(CFLAGS) src/rf_ace.cpp $(SOURCEFILES) $(TFLAGS) -o bin/rf-ace -g -ggdb -pg

static: $(SOURCEFILES)
	$(COMPILER) $(CFLAGS) src/rf_ace.cpp $(STATICFLAGS) $(SOURCEFILES) $(TFLAGS) -o bin/rf-ace

GBT_benchmark: test/GBT_benchmark.cpp $(SOURCEFILES)
	$(COMPILER) $(CFLAGS) test/GBT_benchmark.cpp $(SOURCEFILES) $(TFLAGS) -o bin/GBT_benchmark

test: $(SOURCEFILES) $(TESTFILES)
	rm -f bin/test; $(COMPILER) $(TESTFLAGS) test/run_tests.cpp $(SOURCEFILES) $(TFLAGS) -o bin/test -ggdb; ./bin/test

newtest: $(SOURCEFILES) $(NEWTESTFILES)
	rm -f bin/newtest; $(COMPILER) $(CFLAGS) test/run_newtests.cpp $(SOURCEFILES) $(TFLAGS) -o bin/newtest; ./bin/newtest

newtest-no-threads: $(SOURCEFILES) $(NEWTESTFILES)
	rm -f bin/newtest; $(COMPILER) $(CFLAGS) -DNOTHREADS test/run_newtests.cpp $(SOURCEFILES) -o bin/newtest; ./bin/newtest

test-no-threads: $(SOURCEFILES)
	rm -f bin/test; $(COMPILER) $(TESTFLAGS) -DNOTHREADS test/run_tests.cpp $(SOURCEFILES) -o bin/test -ggdb; ./bin/test  

clean:
	rm -rf bin/rf-ace bin/benchmark bin/GBT_benchmark bin/test bin/*.dSYM/ src/*.o
