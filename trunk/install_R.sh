#!/bin/bash

## Prepare C++ compiler flags
PKG_CPPFLAGS=`Rscript -e "Rcpp:::CxxFlags()"`
PKG_CPPFLAGS=`echo -n $PKG_CPPFLAGS " -std=c++0x"`
export PKG_CPPFLAGS

## Prepare library flags
export PKG_LIBS=`Rscript -e "Rcpp:::LdFlags()"`

## Make shared library
R CMD SHLIB -o lib/rf_ace_R.so src/rf_ace_R.cpp src/treedata.cpp