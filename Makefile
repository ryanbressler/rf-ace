COMPILER = g++
CFLAGS = -O3 -std=c++0x -Wall -Wextra -pedantic -Isrc/
LFLAGS = -lz
TFLAGS = -pthread

COMMON_OBJ = obj/densetreedata.o obj/murmurhash3.o obj/datadefs.o obj/progress.o obj/statistics.o obj/math.o obj/CatPredictThread.o obj/NumPredictThread.o obj/TreeGrowerThread.o obj/stochasticforest.o obj/rootnode.o obj/node.o obj/utils.o obj/distributions.o obj/reader.o obj/feature.o

BIN_OBJ = $(COMMON_OBJ) obj/rf_ace.o

TESTOBJ = $(COMMON_OBJ) obj/run_newtests.o

SOURCEFILES = src/densetreedata.cpp src/murmurhash3.cpp src/datadefs.cpp src/progress.cpp src/statistics.cpp src/math.cpp src/CatPredictThread.cpp src/NumPredictThread.cpp src/TreeGrowerThread.cpp src/stochasticforest.cpp src/rootnode.cpp src/node.cpp src/utils.cpp src/distributions.cpp src/reader.cpp src/feature.cpp

STATICFLAGS = -static-libgcc -static

TESTDEPS = test/rface_test.hpp test/distributions_test.hpp test/argparse_test.hpp test/datadefs_test.hpp test/stochasticforest_test.hpp test/utils_test.hpp test/math_test.hpp test/rootnode_test.hpp test/node_test.hpp test/densetreedata_test.hpp

TESTFLAGS = -std=c++0x -L${HOME}/lib/ -L/usr/local/lib -lcppunit -ldl -pedantic -I${HOME}/include/ -I/usr/local/include -Itest/ -Isrc/
.PHONY: all test clean  # Squash directory checks for the usual suspects

# RF-ACE targets
# ==============
all: rf-ace

obj/datadefs.o: src/datadefs.hpp src/datadefs.cpp
	$(COMPILER) $(CFLAGS) src/datadefs.cpp -c -o obj/datadefs.o

obj/reader.o: src/reader.hpp src/reader.cpp
	$(COMPILER) $(CFLAGS) src/reader.cpp -c -o obj/reader.o

obj/utils.o: src/utils.hpp src/utils.cpp
	$(COMPILER) $(CFLAGS) src/utils.cpp -c -o obj/utils.o

obj/murmurhash3.o: src/murmurhash3.hpp src/murmurhash3.cpp
	$(COMPILER) $(CFLAGS) src/murmurhash3.cpp -c -o obj/murmurhash3.o

obj/progress.o: src/progress.hpp src/progress.cpp
	$(COMPILER) $(CFLAGS) src/progress.cpp -c -o obj/progress.o

obj/feature.o: src/feature.hpp src/feature.cpp
	$(COMPILER) $(CFLAGS) src/feature.cpp -c -o obj/feature.o

obj/node.o: src/node.hpp src/node.cpp
	$(COMPILER) $(CFLAGS) src/node.cpp -c -o obj/node.o

obj/rootnode.o: src/rootnode.hpp src/rootnode.cpp
	$(COMPILER) $(CFLAGS) src/rootnode.cpp -c -o obj/rootnode.o

obj/statistics.o: src/statistics.hpp src/statistics.cpp
	$(COMPILER) $(CFLAGS) src/statistics.cpp -c -o obj/statistics.o

obj/distributions.o: src/distributions.hpp src/distributions.cpp
	$(COMPILER) $(CFLAGS) src/distributions.cpp -c -o obj/distributions.o

obj/math.o: src/math.hpp src/math.cpp
	$(COMPILER) $(CFLAGS) src/math.cpp -c -o obj/math.o

obj/densetreedata.o: src/densetreedata.hpp src/densetreedata.cpp
	$(COMPILER) $(CFLAGS) src/densetreedata.cpp -c -o obj/densetreedata.o

obj/CatPredictThread.o: src/CatPredictThread.hpp src/CatPredictThread.cpp obj/rootnode.o
	$(COMPILER) $(CFLAGS) src/CatPredictThread.cpp -c -o obj/CatPredictThread.o

obj/NumPredictThread.o: src/NumPredictThread.hpp src/NumPredictThread.cpp
	$(COMPILER) $(CFLAGS) src/NumPredictThread.cpp -c -o obj/NumPredictThread.o

obj/TreeGrowerThread.o: src/TreeGrowerThread.hpp src/TreeGrowerThread.cpp
	$(COMPILER) $(CFLAGS) src/TreeGrowerThread.cpp -c -o obj/TreeGrowerThread.o

obj/stochasticforest.o: src/stochasticforest.hpp src/stochasticforest.cpp
	$(COMPILER) $(CFLAGS) src/stochasticforest.cpp -c -o obj/stochasticforest.o

obj/rf_ace.o: src/rf_ace.hpp src/rf_ace.cpp
	$(COMPILER) $(CFLAGS) src/rf_ace.cpp -c -o obj/rf_ace.o

rf-ace: $(BIN_OBJ)
	$(COMPILER) $(CFLAGS) $(LFLAGS) $(TFLAGS) $(BIN_OBJ) -o bin/rf-ace

rf-ace-i386: $(BIN_OBJ)
	$(COMPILER) $(CFLAGS) $(LFLAGS) $(TFLAGS) -m32 $(BIN_OBJ) -o bin/rf-ace-i386

rf-ace-amd64: $(BIN_OBJ)
	$(COMPILER) $(CFLAGS) $(LFLAGS) $(TFLAGS) -m64 $(BIN_OBJ) -o bin/rf-ace-amd64

no-threads: $(BIN_OBJ)
	$(COMPILER) $(CFLAGS) -DNOTHREADS $(LFLAGS) $(BIN_OBJ) -o bin/rf-ace

debug: $(BIN_OBJ)
	$(COMPILER) $(CFLAGS) $(LFLAGS) $(TFLAGS) $(BIN_OBJ) -o bin/rf-ace -g -ggdb -pg

static: $(SOURCEFILES)
	$(COMPILER) $(CFLAGS) src/rf_ace.cpp $(STATICFLAGS) $(SOURCEFILES) $(TFLAGS) -o bin/rf-ace

static-no-threads: $(SOURCEFILES)
	$(COMPILER) $(CFLAGS) -DNOTHREADS $(STATICFLAGS) $(SOURCEFILES) src/rf_ace.cpp -o bin/rf-ace

GBT_benchmark: test/GBT_benchmark.cpp $(SOURCEFILES)
	$(COMPILER) $(CFLAGS) $(TFLAGS) test/GBT_benchmark.cpp $(SOURCEFILES) -o bin/GBT_benchmark

# Test targets
# ============
obj/run_newtests.o: test/run_newtests.cpp
	$(COMPILER) $(CFLAGS) test/run_newtests.cpp -c -o obj/run_newtests.o

test: $(TESTOBJ)
	rm -f bin/newtest; $(COMPILER) $(CFLAGS) $(TFLAGS) $(TESTOBJ) -o bin/newtest -ggdb; ./bin/newtest

test-no-threads: $(TESTOBJ)
	rm -f bin/newtest; $(COMPILER) $(CFLAGS) -DNOTHREADS $(TESTOBJ) -o bin/newtest -ggdb; ./bin/newtest

# Misc
# ====
clean:
	rm -f bin/rf-ace bin/rf-ace-i386 bin/rf-ace-amd64 bin/benchmark bin/GBT_benchmark bin/test bin/*.dSYM/ src/*.o obj/*.o
