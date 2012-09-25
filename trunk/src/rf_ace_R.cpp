#include <Rcpp.h>

#include <cstdlib>
#include <cassert>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>
#include <cmath>
#include <stdio.h>
#include <iomanip>
#include <algorithm>
#include <map>
#include <unordered_map>
#include <set>
//#include <tuple>

#include "rf_ace.hpp"
#include "stochasticforest.hpp"
#include "treedata.hpp"
#include "datadefs.hpp"
#include "utils.hpp"
#include "options.hpp"
#include "statistics.hpp"
#include "progress.hpp"
#include "math.hpp"
#include "distributions.hpp"
#include "timer.hpp"

using namespace std;
using datadefs::num_t;

Timer* TIMER_G = NULL;

RcppExport SEXP rfaceTrain(SEXP trainDataFile, SEXP targetStr, SEXP nTrees, SEXP mTry, SEXP nodeSize, SEXP nMaxLeaves) {

  rface::printHeader(cout);

  TIMER_G = new Timer();

  TIMER_G->tic("TOTAL");

  options::General_options params;

  params.targetStr  = Rcpp::as<string>(targetStr);
  params.nTrees     = Rcpp::as<size_t>(nTrees);
  params.mTry       = Rcpp::as<size_t>(mTry);
  params.nodeSize   = Rcpp::as<size_t>(nodeSize);
  params.nMaxLeaves = Rcpp::as<size_t>(nMaxLeaves);

  Treedata trainData(Rcpp::as<string>(trainDataFile),&params);

  rface::updateTargetStr(trainData,params);

  rface::pruneFeatureSpace(trainData,params);

  //rface::setEnforcedForestParameters(trainData,params);

  rface::printGeneralSetup(trainData,params);

  params.print();

  params.validateParameters();

  if ( params.modelType == options::RF ) {
    cout << "===> Growing RF predictor... " << flush;
  } else if ( params.modelType == options::GBT ) {
    cout << "===> Growing GBT predictor... " << flush;
  } else if ( params.modelType == options::CART ) {
    cout << "===> Growing CART predictor... " << flush;
  } else {
    cerr << "Unknown forest type!" << endl;
    exit(1);
  }

  Rcpp::XPtr<StochasticForest> predictor( new StochasticForest(&trainData,&params) );
  cout << "DONE" << endl << endl;

  if ( params.modelType == options::GBT ) {
    cout << "GBT diagnostics disabled temporarily" << endl << endl;
    return predictor;
  }

  size_t targetIdx = trainData.getFeatureIdx(params.targetStr);
  vector<num_t> data = utils::removeNANs(trainData.getFeatureData(targetIdx));

  num_t oobError = predictor->getOobError();
  num_t ibOobError =  predictor->getError();

  cout << "RF training error measures (NULL == no model):" << endl;
  if ( trainData.isFeatureNumerical(targetIdx) ) {
    num_t nullError = math::var(data);
    cout << "              NULL std = " << sqrt( nullError ) << endl;
    cout << "               OOB std = " << sqrt( oobError ) << endl;
    cout << "            IB+OOB std = " << sqrt( ibOobError ) << endl;
    cout << "  % explained by model = " << 1 - oobError / nullError << " = 1 - (OOB var) / (NULL var)" << endl;
  } else {
    num_t nullError = math::nMismatches( data, math::mode(data) );
    cout << "       NULL % mispred. = " << 1.0 * nullError / data.size() << endl;
    cout << "        OOB % mispred. = " << oobError << endl;
    cout << "     IB+OOB % mispred. = " << ibOobError << endl;
    cout << "  % explained by model = " << 1 - oobError / nullError << " ( 1 - (OOB # mispred.) / (NULL # mispred.) )" << endl;
  }
  cout << endl;

  TIMER_G->toc("TOTAL");
  TIMER_G->print();

  delete TIMER_G;

  return predictor; 

}

RcppExport SEXP rfacePredict(SEXP predictorObj, SEXP testDataFile) {

  Rcpp::XPtr<StochasticForest> predictor(predictorObj);
  
  vector<num_t> prediction,confidence;

  options::General_options params;

  Treedata testData(Rcpp::as<string>(testDataFile),&params);

  //predictor->predict(&testData,prediction,confidence);

  return Rcpp::wrap(NULL);

}
