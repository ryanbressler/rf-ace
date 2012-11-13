#include <Rcpp.h>
#include <cstdlib>
#include <vector>

#include "rf_ace.hpp"
#include "treedata.hpp"
#include "utils.hpp"
#include "datadefs.hpp"

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

  rface->saveForest(Rcpp::as<string>(fileName));

}

RcppExport SEXP rfaceLoad(SEXP rfaceFile, SEXP nThreads) {

  
  Rcpp::XPtr<RFACE> rface( new RFACE, true);

  rface->loadForest(Rcpp::as<string>(rfaceFile));

  return(rface);

}

RcppExport SEXP rfaceTrain(SEXP trainDataFrameObj, SEXP targetStrR, SEXP nTreesR, SEXP mTryR, SEXP nodeSizeR, SEXP nMaxLeavesR, SEXP nThreadsR) {

  ForestOptions forestOptions;

  targetStr = Rcpp::as<string>(targetStrR);
  forestOptions.nTrees     = Rcpp::as<size_t>(nTreesR);
  forestOptions.mTry       = Rcpp::as<size_t>(mTryR);
  forestOptions.nodeSize   = Rcpp::as<size_t>(nodeSizeR);
  forestOptions.nMaxLeaves = Rcpp::as<size_t>(nMaxLeavesR);
  nThreads   = Rcpp::as<size_t>(nThreadsR);

  vector<Feature> dataMatrix;

  parseDataFrame(trainDataFrameObj,dataMatrix,sampleHeaders);

  bool useContrasts = false;
  Treedata trainData(dataMatrix,useContrasts,sampleHeaders);

  size_t targetIdx = trainData.getFeatureIdx(targetStr);

  if ( targetIdx == trainData.end() ) {
    int integer;
    if ( utils::isInteger(targetStr,integer) && integer >= 0 && integer < trainData.nFeatures() ) {
      targetIdx = static_cast<size_t>(integer);
    } else {
      cerr << "Invalid target: " << targetStr << endl;
      exit(1);
    }
  }

  Rcpp::XPtr<RFACE> rface( new RFACE, true);

  vector<num_t> featureWeights(trainData.nFeatures(),1.0);
  featureWeights[targetIdx] = 0.0;

  int seed = 0;

  rface->train(trainData,targetIdx,featureWeights,forestOptions,seed,nThreads);

  return(rface);

}

