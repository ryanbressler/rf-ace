COMPILER = g++
CFLAGS = -pg -O3 -Wall -pedantic

all: rf_ace

rf_ace: src/rf_ace.cpp src/node.cpp
	$(COMPILER) $(CFLAGS) src/rf_ace.cpp src/node.cpp -o bin/rf_ace

benchmark: src/benchmark.cpp src/node.cpp
	$(COMPILER) $(CFLAGS) src/benchmark.cpp src/node.cpp -o bin/benchmark

clean:
	rm -rf bin/rf_ace bin/benchmark src/*.o