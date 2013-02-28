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
void rface_newtest_QRF_train_test_regression();
void rface_newtest_GBT_train_test_classification();
void rface_newtest_GBT_train_test_regression();
void rface_newtest_RF_save_load_classification();
void rface_newtest_RF_save_load_regression();
void rface_newtest_QRF_save_load_regression();
void rface_newtest_GBT_save_load_classification();
void rface_newtest_GBT_save_load_regression();

void rface_newtest() {
  
  newtest( "Testing RF for classification", &rface_newtest_RF_train_test_classification );
  newtest( "Testing RF for regression", &rface_newtest_RF_train_test_regression );
  newtest( "Testing Quantile Regression Random Forest", &rface_newtest_QRF_train_test_regression );
  //newtest( "Testing GBT for classification", &rface_newtest_GBT_train_test_classification );
  //newtest( "Testing GBT for regression", &rface_newtest_GBT_train_test_regression );
  newtest( "Testing save/load RF for classification", &rface_newtest_RF_save_load_classification );
  newtest( "Testing save/load RF for regression", &rface_newtest_RF_save_load_regression );
  newtest( "Testing save/load Quantile Regression Random Forest", &rface_newtest_QRF_save_load_regression );
  //newtest( "Testing save/load GBT for classification", &rface_newtest_GBT_save_load_classification );
  //newtest( "Testing save/load GBT for regression", &rface_newtest_GBT_save_load_regression );

}

RFACE::TestOutput make_predictions(ForestOptions& forestOptions, const string& targetStr) {

  string fileName = "test_103by300_mixed_nan_matrix.afm";
  Treedata trainData(fileName,'\t',':',false);
  size_t targetIdx = trainData.getFeatureIdx(targetStr);
  vector<num_t> weights = trainData.getFeatureWeights();
  weights[targetIdx] = 0;

  RFACE rface;

  rface.train(&trainData,targetIdx,weights,&forestOptions);

  return( rface.test(&trainData) );

}

RFACE::QuantilePredictionOutput make_quantile_predictions(ForestOptions& forestOptions, const string& targetStr) {

  string fileName = "test_103by300_mixed_nan_matrix.afm";
  Treedata trainData(fileName,'\t',':',false);
  size_t targetIdx = trainData.getFeatureIdx(targetStr);
  vector<num_t> weights = trainData.getFeatureWeights();
  weights[targetIdx] = 0;

  RFACE rface;

  rface.train(&trainData,targetIdx,weights,&forestOptions);

  return( rface.predictQuantiles(&trainData,forestOptions.quantiles,forestOptions.nSamplesForQuantiles) );

}



RFACE::TestOutput make_save_load_predictions(ForestOptions& forestOptions, const string& targetStr) {

  string fileName = "test_103by300_mixed_nan_matrix.afm";
  Treedata trainData(fileName,'\t',':',false);
  size_t targetIdx = trainData.getFeatureIdx(targetStr);
  vector<num_t> weights = trainData.getFeatureWeights();
  weights[targetIdx] = 0;

  RFACE rface;

  rface.train(&trainData,targetIdx,weights,&forestOptions);
  
  rface.save("foo.sf");

  RFACE rface2;
  
  rface2.load("foo.sf");

  return( rface2.test(&trainData) );
  
}

RFACE::QuantilePredictionOutput make_save_load_quantile_predictions(ForestOptions& forestOptions, const string& targetStr) {
  
  string fileName = "test_103by300_mixed_nan_matrix.afm";
  Treedata trainData(fileName,'\t',':',false);
  size_t targetIdx = trainData.getFeatureIdx(targetStr);
  vector<num_t> weights = trainData.getFeatureWeights();
  weights[targetIdx] = 0;
  
  RFACE rface;
  
  rface.train(&trainData,targetIdx,weights,&forestOptions);
  
  rface.save("foo.sf");
  
  RFACE rface2;
  
  rface2.load("foo.sf");
  
  return( rface2.predictQuantiles(&trainData,forestOptions.quantiles,forestOptions.nSamplesForQuantiles) );

}


num_t classification_error(const RFACE::TestOutput& predictions) { 
 
  num_t pError = 0.0;
  num_t n = static_cast<num_t>(predictions.catPredictions.size());
  for ( size_t i = 0; i < predictions.catPredictions.size(); ++i ) {
    pError += (predictions.catPredictions[i] != predictions.catTrueData[i]) / n;
  }
  
  return(pError);

}

num_t regression_error(const RFACE::TestOutput& predictions) {

  num_t RMSE = 0.0;
  num_t n = static_cast<num_t>(predictions.numPredictions.size());
  for ( size_t i = 0; i < predictions.numPredictions.size(); ++i ) {
    num_t e = predictions.numPredictions[i] - predictions.numTrueData[i];
    RMSE += powf(e,2)/n;
  }
  
  return( sqrt(RMSE) );

}

vector<num_t> quantile_regression_error(const RFACE::QuantilePredictionOutput& qPredOut) {

  vector<num_t> QDEV(qPredOut.quantiles.size(),0.0);
  num_t n = static_cast<num_t>(qPredOut.predictions.size());

  for ( size_t q = 0; q < qPredOut.quantiles.size(); ++q ) {
    for ( size_t i = 0; i < qPredOut.predictions.size(); ++i ) {
      bool b = qPredOut.trueData[i] < qPredOut.predictions[i][q];
      QDEV[q] += b/n;
    }
    QDEV[q] = fabs(QDEV[q] - qPredOut.quantiles[q]);
  }

  return(QDEV);

}

void rface_newtest_RF_train_test_classification() {
  
  ForestOptions forestOptions(forest_t::RF);
  forestOptions.mTry = 30;

  num_t pError = classification_error( make_predictions(forestOptions,"C:class") );

  newassert(pError < 0.2);

}

void rface_newtest_RF_train_test_regression() {
  
  ForestOptions forestOptions(forest_t::RF);
  forestOptions.mTry = 30;

  num_t RMSE = regression_error( make_predictions(forestOptions,"N:output") );
  
  newassert(RMSE < 1.0);

}

void rface_newtest_QRF_train_test_regression() {

  ForestOptions forestOptions(forest_t::QRF);
  forestOptions.mTry = 30;

  vector<num_t> QDEV = quantile_regression_error( make_quantile_predictions(forestOptions,"N:output") );

  utils::write(cout,QDEV.begin(),QDEV.end());

  newassert( math::mean(QDEV) < 0.05 );
  
}

void rface_newtest_GBT_train_test_classification() {

  ForestOptions forestOptions(forest_t::GBT);

  num_t pError = classification_error( make_predictions(forestOptions,"C:class") );

  newassert( pError < 0.2 );

}

void rface_newtest_GBT_train_test_regression() { 

  ForestOptions forestOptions(forest_t::GBT);

  num_t RMSE = regression_error( make_predictions(forestOptions,"N:output") );

  newassert(RMSE < 1.0);

}

void rface_newtest_RF_save_load_classification() {

  ForestOptions forestOptions(forest_t::RF);
  forestOptions.mTry = 30;

  num_t pError1 = classification_error( make_predictions(forestOptions,"C:class") );
  num_t pError2 = classification_error( make_save_load_predictions(forestOptions,"C:class") );

  newassert( fabs(pError1 - pError2) < 1e-1 );

}

void rface_newtest_RF_save_load_regression() {

  ForestOptions forestOptions(forest_t::RF);
  forestOptions.mTry = 30;

  num_t RMSE1 = classification_error( make_predictions(forestOptions,"N:output") );
  num_t RMSE2 = classification_error( make_save_load_predictions(forestOptions,"N:output") );

  newassert( fabs(RMSE1 - RMSE2) < 1e-1 );

}

void rface_newtest_QRF_save_load_regression() {

  ForestOptions forestOptions(forest_t::QRF);
  forestOptions.mTry = 30;
  
  vector<num_t> QRMSE1 = quantile_regression_error( make_quantile_predictions(forestOptions,"N:output") );
  vector<num_t> QRMSE2 = quantile_regression_error( make_save_load_quantile_predictions(forestOptions,"N:output") );
  
  newassert( fabs(QRMSE1[2] - QRMSE2[2]) < 1e-1 );


}

void rface_newtest_GBT_save_load_classification() {

  ForestOptions forestOptions(forest_t::GBT);

  num_t pError1 = classification_error( make_predictions(forestOptions,"C:class") );
  num_t pError2 = classification_error( make_save_load_predictions(forestOptions,"C:class") );

  newassert( fabs(pError1 - pError2) < 1e-1 );
  
}

void rface_newtest_GBT_save_load_regression() {

  ForestOptions forestOptions(forest_t::GBT);

  num_t RMSE1 = classification_error( make_predictions(forestOptions,"N:output") );
  num_t RMSE2 = classification_error( make_save_load_predictions(forestOptions,"N:output") );

  newassert( fabs(RMSE1 - RMSE2) < 1e-1 );

}

#endif
