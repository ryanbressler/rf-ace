#include <Rcpp.h>
#include <cstdlib>
#include <vector>

#include "rf_ace.hpp"
#include "treedata.hpp"
#include "datadefs.hpp"
#include "options.hpp"

using namespace std;
using datadefs::num_t;

void parseDataFrame(SEXP dataFrameObj, vector<Feature>& dataMatrix, vector<string>& sampleHeaders) {

  Rcpp::DataFrame df(dataFrameObj);

  //Rcpp::CharacterVector colNames = df.attr("names");
  //Rcpp::CharacterVector rowNames = df.attr("row.names");

  vector<string> featureHeaders = df.attr("names");
  vector<string> foo = df.attr("row.names");
  sampleHeaders = foo;

  dataMatrix.resize( 0 );

  //cout << "nf = " << featureHeaders.size() << endl;
  //cout << "ns = " << sampleHeaders.size() << endl;

  // Read one column of information, which in this case is assumed to be one sample
  for ( size_t i = 0; i < featureHeaders.size(); ++i ) {
    Rcpp::List vec = df[i];
    assert(vec.length() == sampleHeaders.size() );
    //cout << " " << foo[0] << flush;
    //cout << " df[" << i << "].length() = " << vec.length() << endl;
    if ( featureHeaders[i].substr(0,2) != "N:" ) {
      vector<string> sVec(sampleHeaders.size());
      for ( size_t j = 0; j < sampleHeaders.size(); ++j ) {
        //cout << Rcpp::as<string>(vec[j]) << endl;
        sVec[j] = Rcpp::as<string>(vec[j]);
      }
      //cout << "sVec = ";
      //utils::write(cout,sVec.begin(),sVec.end());
      //cout << endl;
      dataMatrix.push_back( Feature(sVec,featureHeaders[i]) );
    } else {
      vector<num_t> sVec(sampleHeaders.size());
      for ( size_t j = 0; j < sampleHeaders.size(); ++j ) {
        sVec[j] = Rcpp::as<num_t>(vec[j]);
      }
      dataMatrix.push_back( Feature(sVec,featureHeaders[i]) );
    }

    //  cout << "df[" << j << "," << i << "] = " << Rcpp::as<num_t>(vec[j]) << endl;
    // }
  }

  assert( dataMatrix.size() == featureHeaders.size() );

}

RcppExport void rfaceSave(SEXP rfaceObj, SEXP fileName) {

  Rcpp::XPtr<RFACE> rface(rfaceObj);

  rface->save(Rcpp::as<string>(fileName));

}

RcppExport SEXP rfaceLoad(SEXP rfaceFile) {

  
  Rcpp::XPtr<RFACE> rface( new RFACE, true);

  rface->load(Rcpp::as<string>(rfaceFile));

  return(rface);

}

RcppExport SEXP rfaceTrain(SEXP trainDataFrameObj, SEXP targetStrR, SEXP featureWeightsR, SEXP nTreesR, SEXP mTryR, SEXP nodeSizeR, SEXP nMaxLeavesR, SEXP nThreadsR) {

  ForestOptions forestOptions;

  string targetStr = Rcpp::as<string>(targetStrR);
  forestOptions.nTrees     = Rcpp::as<size_t>(nTreesR);
  forestOptions.mTry       = Rcpp::as<size_t>(mTryR);
  forestOptions.nodeSize   = Rcpp::as<size_t>(nodeSizeR);
  forestOptions.nMaxLeaves = Rcpp::as<size_t>(nMaxLeavesR);
  size_t nThreads = Rcpp::as<size_t>(nThreadsR);

  vector<Feature> dataMatrix;
  vector<string> sampleHeaders;

  parseDataFrame(trainDataFrameObj,dataMatrix,sampleHeaders);

  bool useContrasts = false;
  Treedata trainData(dataMatrix,useContrasts,sampleHeaders);

  size_t targetIdx = trainData.getFeatureIdx(targetStr);

  if ( targetIdx == trainData.end() ) {
    int integer;
    if ( datadefs::isInteger(targetStr,integer) && integer >= 0 && integer < static_cast<int>(trainData.nFeatures()) ) {
      targetIdx = static_cast<size_t>(integer);
    } else {
      cerr << "Invalid target: " << targetStr << endl;
      exit(1);
    }
  }

  Rcpp::XPtr<RFACE> rface( new RFACE, true);

  Rcpp::NumericVector foo(featureWeightsR);
  vector<num_t> featureWeights(foo.size());
  for ( size_t i = 0; i < featureWeights.size(); ++i ) {
    featureWeights[i] = foo[i];
  }

  //vector<num_t> featureWeights(trainData.nFeatures(),1.0);
  featureWeights[targetIdx] = 0.0;

  int seed = 0;

  rface->train(&trainData,targetIdx,featureWeights,&forestOptions,seed,nThreads);

  return(rface);

}

RcppExport SEXP rfacePredict(SEXP rfaceObj, SEXP testDataFrameObj, SEXP nThreadsR) {

  size_t nThreads = Rcpp::as<size_t>(nThreadsR);

  Rcpp::XPtr<RFACE> rface(rfaceObj);

  vector<Feature> testDataMatrix;
  vector<string> sampleHeaders;

  parseDataFrame(testDataFrameObj,testDataMatrix,sampleHeaders);

  bool useContrasts = false;

  Treedata testData(testDataMatrix,useContrasts,sampleHeaders);

  RFACE::TestOutput testOutput = rface->test(&testData,nThreads);

  Rcpp::List predictions;

  if ( testOutput.isTargetNumerical ) {
    
    predictions = Rcpp::List::create(Rcpp::Named("targetName")=testOutput.targetName,
				     Rcpp::Named("sampleNames")=testOutput.sampleNames,
				     Rcpp::Named("trueData")=testOutput.numTrueData,
				     Rcpp::Named("predData")=testOutput.numPredictions,
				     Rcpp::Named("predError")=testOutput.confidence);
    
  } else {

    predictions = Rcpp::List::create(Rcpp::Named("targetName")=testOutput.targetName,
				     Rcpp::Named("sampleNames")=testOutput.sampleNames,
				     Rcpp::Named("trueData")=testOutput.catTrueData,
				     Rcpp::Named("predData")=testOutput.catPredictions,
				     Rcpp::Named("predError")=testOutput.confidence);
    
  }
  
  return(predictions);

}

RcppExport SEXP rfaceFilter(SEXP filterDataFrameObj,  SEXP targetStrR, SEXP featureWeightsR, SEXP nTreesR, SEXP mTryR, SEXP nodeSizeR, SEXP nMaxLeavesR, SEXP nThreadsR) {

  string targetStr = Rcpp::as<string>(targetStrR);

  ForestOptions forestOptions;
  forestOptions.nTrees = Rcpp::as<size_t>(nTreesR);
  forestOptions.mTry = Rcpp::as<size_t>(mTryR);
  forestOptions.nodeSize = Rcpp::as<size_t>(nodeSizeR);
  forestOptions.nMaxLeaves = Rcpp::as<size_t>(nMaxLeavesR);

  size_t nThreads = Rcpp::as<size_t>(nThreadsR);

  FilterOptions filterOptions;

  vector<Feature> dataMatrix;
  vector<string> sampleHeaders;

  parseDataFrame(filterDataFrameObj,dataMatrix,sampleHeaders);

  bool useContrasts = true;

  Treedata filterData(dataMatrix,useContrasts,sampleHeaders);
 
  size_t seed = 0;

  size_t targetIdx = filterData.getFeatureIdx(targetStr);

  if ( targetIdx == filterData.end() ) {
    int integer;
    if ( datadefs::isInteger(targetStr,integer) && integer >= 0 && integer < static_cast<int>(filterData.nFeatures()) ) {
      targetIdx = static_cast<size_t>(integer);
    } else {
      cerr << "Invalid target: " << targetStr << endl;
      exit(1);
    }
  }

  Rcpp::NumericVector foo(featureWeightsR);
  vector<num_t> featureWeights(foo.size());
  for ( size_t i = 0; i < featureWeights.size(); ++i ) {
    featureWeights[i] = foo[i];
  }
  featureWeights[targetIdx] = 0.0;

  RFACE rface;

  RFACE::FilterOutput filterOutput = rface.filter(&filterData,targetIdx,featureWeights,&forestOptions,&filterOptions,seed,nThreads);

  Rcpp::List filterOutputR = Rcpp::List::create(Rcpp::Named("featureNames")=filterOutput.featureNames,
                                                Rcpp::Named("pValues")=filterOutput.pValues,
                                                Rcpp::Named("importances")=filterOutput.importances,
                                                Rcpp::Named("correlations")=filterOutput.correlations,
                                                Rcpp::Named("sampleCounts")=filterOutput.sampleCounts);
  

  return(filterOutputR);
  
}

