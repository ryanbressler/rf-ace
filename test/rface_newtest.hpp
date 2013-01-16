#ifndef RFACE_NEWTEST_HPP
#define RFACE_NEWTEST_HPP

#include <cstdlib>
#include <cmath>
#include "options.hpp"
#include "treedata.hpp"
#include "rf_ace.hpp"
#include "newtest.hpp"

using namespace std;
using datadefs::num_t;

void rface_newtest_RF_train_test_classification();
void rface_newtest_RF_train_test_regression();
void rface_newtest_GBT_train_test_classification();
void rface_newtest_GBT_train_test_regression();

void rface_newtest() {
  
  newtest( "Making an RF classification experiment", &rface_newtest_RF_train_test_classification );
  newtest( "Making an RF regression experiment", &rface_newtest_RF_train_test_regression );
  newtest( "Making a GBT classification experiment", &rface_newtest_GBT_train_test_classification );
  newtest( "Making a GBT regrssion experiment", &rface_newtest_GBT_train_test_regression );
  
}

RFACE::TestOutput make_predictions(ForestOptions& forestOptions,const string& targetStr) {

  string fileName = "test_103by300_mixed_nan_matrix.afm";
  Treedata trainData(fileName,'\t',':',false);
  size_t targetIdx = trainData.getFeatureIdx(targetStr);
  vector<num_t> weights = trainData.getFeatureWeights();
  weights[targetIdx] = 0;

  RFACE rface;

  rface.train(&trainData,targetIdx,weights,&forestOptions);
  
  return( rface.test(&trainData) );
  
}

num_t classification_error(RFACE::TestOutput& testOutput) { 
 
  return(0);

}

void rface_newtest_RF_train_test_classification() {
  
  // Set parameters
  RFACE rface;
  ForestOptions forestOptions;
  forestOptions.setRFDefaults();
  forestOptions.mTry = 50;

  RFACE::TestOutput predictions = make_predictions(forestOptions,"C:class");

  // Calculate error
  num_t pError = 0.0;
  num_t n = static_cast<num_t>(predictions.catPredictions.size());
  for ( size_t i = 0; i < predictions.catPredictions.size(); ++i ) {
    pError += (predictions.catPredictions[i] != predictions.catTrueData[i]) / n;
  }
  pError = sqrt(pError);
  
  newassert(pError < 0.2);

}

void rface_newtest_RF_train_test_regression() {
  
  // Set parameters
  RFACE rface;
  ForestOptions forestOptions;
  forestOptions.setRFDefaults();
  forestOptions.mTry = 50;

  // Train & Test
  RFACE::TestOutput predictions = make_predictions(forestOptions,"N:output");
  
  // Calculate error
  num_t RMSE = 0.0;
  num_t n = static_cast<num_t>(predictions.numPredictions.size());
  for ( size_t i = 0; i < predictions.numPredictions.size(); ++i ) {
    num_t e = predictions.numPredictions[i] - predictions.numTrueData[i];
    RMSE += powf(e,2)/n;
  }
  RMSE = sqrt(RMSE);

  newassert(RMSE < 1.0);

}

void rface_newtest_GBT_train_test_classification() {

  

}

void rface_newtest_GBT_train_test_regression() { 

}

#endif
