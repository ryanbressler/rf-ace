#!/bin/bash

## Prepare C++ compiler flags
export PKG_CPPFLAGS="`Rscript -e 'Rcpp:::CxxFlags()'` -std=c++0x -Wall -Wextra -pedantic"

## Prepare library flags
export PKG_LIBS=`Rscript -e "Rcpp:::LdFlags()"`

## Make shared library
R CMD SHLIB -o lib/rf_ace_R.so src/rf_ace_R.cpp src/progress.cpp src/statistics.cpp src/math.cpp src/stochasticforest.cpp src/rootnode.cpp src/node.cpp src/treedata.cpp src/datadefs.cpp src/utils.cpp src/distributions.cpp