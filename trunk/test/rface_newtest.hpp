#ifndef RFACE_NEWTEST_HPP
#define RFACE_NEWTEST_HPP

#include <cstdlib>
#include "options.hpp"
#include "treedata.hpp"
#include "rf_ace.hpp"
#include "newtest.hpp"

void rface_newtest_train_test_classification();

void rface_newtest() {

  newtest( "Making an RF classification experiment", &rface_newtest_train_test_classification );

}


void rface_newtest_train_test_classification() {

  string fileName = "test_103by300_mixed_nan_matrix.afm";

  Treedata trainData(fileName,'\t',':',false);

  size_t targetIdx = trainData.getFeatureIdx("N:output");

  vector<num_t> weights = trainData.getFeatureWeights();
  weights[targetIdx] = 0;

  RFACE rface;

  ForestOptions forestOptions;

  forestOptions.setRFDefaults();
  forestOptions.mTry = 30;

  rface.train(&trainData,targetIdx,weights,&forestOptions);

}

#endif
