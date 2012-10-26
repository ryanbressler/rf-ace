#include <cstdlib>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <iomanip>

#include "rf_ace.hpp"
#include "datadefs.hpp"
#include "options.hpp"

using namespace std;
using datadefs::num_t;
using datadefs::ftable_t;

int main(const int argc, char* const argv[]) {

  RFACE rface;

  rface.printHeader(cout);
  
  // Structs that store all the user-specified command-line arguments
  options::General_options gen_op(argc,argv);
  
  // With no input arguments the help is printed
  if ( argc == 1 || gen_op.printHelp ) {
    gen_op.help();
    return(EXIT_SUCCESS);
  }

  if ( gen_op.isFilter ) {

    if ( gen_op.nPerms > 1 ) {
      gen_op.useContrasts = true;
    }

    rface.filter(gen_op);

  } else if ( gen_op.isSet(gen_op.recombinePerms_s,gen_op.recombinePerms_l) ) {
    
    cout << " *(EXPERIMENTAL) RF-ACE RECOMBINER (" << gen_op.recombinePerms << " permutations) ACTIVATED* " << endl;
        
    if ( gen_op.recombinePerms == 0 ) {
      cerr << "Currently the number of permutations to be recombined ( -" 
	   << gen_op.recombinePerms_s << " / --" << gen_op.recombinePerms_l << endl
	   << " ) needs to be explicitly specified." << endl;
      exit(1);
    }

    rface.recombine(gen_op);

  } else {

    rface.trainAndTest(gen_op);

  }


}


#ifdef RFACER
#include <Rcpp.h>

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

RcppExport void rfaceSave(SEXP predictorObj, SEXP fileName) {

  Rcpp::XPtr<StochasticForest> predictor(predictorObj);

  predictor->printToFile(Rcpp::as<string>(fileName));

}

RcppExport SEXP rfaceLoad(SEXP predictorFile, SEXP nThreads) {

  options::General_options params;

  params.forestInput = Rcpp::as<string>(predictorFile);
  params.nThreads = Rcpp::as<size_t>(nThreads);

  params.initRandIntGens();

  Rcpp::XPtr<StochasticForest> predictor( new StochasticForest(params), true);

  return(predictor);

}

RcppExport SEXP rfaceTrain(SEXP trainDataFrameObj, SEXP targetStr, SEXP nTrees, SEXP mTry, SEXP nodeSize, SEXP nMaxLeaves, SEXP nThreads) {

  rface.printHeader(cout);

  TIMER_G = new Timer();

  TIMER_G->tic("TOTAL");

  options::General_options params;

  params.targetStr  = Rcpp::as<string>(targetStr);
  params.nTrees     = Rcpp::as<size_t>(nTrees);
  params.mTry       = Rcpp::as<size_t>(mTry);
  params.nodeSize   = Rcpp::as<size_t>(nodeSize);
  params.nMaxLeaves = Rcpp::as<size_t>(nMaxLeaves);
  params.nThreads   = Rcpp::as<size_t>(nThreads);

  params.initRandIntGens();

  vector<Feature> dataMatrix;
  vector<string> sampleHeaders;

  TIMER_G->tic("READ");
  parseDataFrame(trainDataFrameObj,dataMatrix,sampleHeaders);
  TIMER_G->toc("READ");

  //return(Rcpp::wrap(NULL));

  Treedata trainData(dataMatrix,&params,sampleHeaders);

  //StochasticForest predictor = rface.buildPredictor(trainData,params);

  //Rcpp::XPtr<StochasticForest> predictorObj( &predictor, true );

  rface.updateTargetStr(trainData,params);
  
  rface.pruneFeatureSpace(trainData,params);
  
  //rface.setEnforcedForestParameters(trainData,params);
  
  rface.printGeneralSetup(trainData,params);
  
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
  
  Rcpp::XPtr<StochasticForest> predictor( new StochasticForest(&trainData,params), true );
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

RcppExport SEXP rfacePredict(SEXP predictorObj, SEXP testDataFrameObj) {

  TIMER_G = new Timer();

  TIMER_G->tic("PREDICTION");

  Rcpp::XPtr<StochasticForest> predictor(predictorObj);
  //Rcpp::XPtr<options::General_options> params(paramsObj);

  vector<Feature> testDataMatrix;
  vector<string> sampleHeaders;

  parseDataFrame(testDataFrameObj,testDataMatrix,sampleHeaders);

  Treedata testData(testDataMatrix,predictor->params(),sampleHeaders);

  string targetName = predictor->params()->targetStr;
  size_t targetIdx = testData.getFeatureIdx(targetName);
  

  vector<string> trueData(testData.nSamples(),datadefs::STR_NAN);
  if ( targetIdx != testData.end() ) {
    trueData = testData.getRawFeatureData(targetIdx);
  }

  vector<string> prediction;
  vector<num_t> confidence;
  predictor->predict(&testData,prediction,confidence);

  Rcpp::List predictions;

  for(size_t i = 0; i < prediction.size(); ++i) {
    predictions.push_back( Rcpp::List::create(Rcpp::Named("target")=targetName,
					      Rcpp::Named("sample")=testData.getSampleName(i),
					      Rcpp::Named("true")=trueData[i],
					      Rcpp::Named("prediction")=prediction[i],
					      Rcpp::Named("error")=confidence[i]));
  }

  TIMER_G->toc("PREDICTION");
  TIMER_G->print();

  delete TIMER_G;


  return predictions;

}

#endif
