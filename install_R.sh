#!/bin/bash

#PKG_CPPFLAGS=-std=c++0x
PKG_CPPFLAGS=`Rscript -e "Rcpp:::CxxFlags()"`
PKG_CPPFLAGS=`echo -n $PKG_CPPFLAGS " -std=c++0x"`
export PKG_CPPFLAGS #=`Rscript -e 'Rcpp:::CxxFlags()'`
export PKG_LIBS=`Rscript -e 'Rcpp:::LdFlags()'`

R CMD SHLIB -o rf_ace_R.so src/rf_ace_R.cpp src/treedata.cpp