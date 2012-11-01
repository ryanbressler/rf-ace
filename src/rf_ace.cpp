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

struct ProgramLogic {

  bool filter;
  bool recombine;
  bool trainModel;
  bool testModel;
  bool saveModel;
  bool loadModel;

  void parse(options::General_options* params) {
    
    bool isTargetSet     = params->isSet(params->targetStr_s,params->targetStr_l);
    bool isFilterSet     = params->isFilter;
    bool isRecombineSet  = params->isSet(params->recombinePerms_s,params->recombinePerms_l);
    bool isInputSet      = params->isSet(params->input_s,params->input_l);
    bool isPredDataSet   = params->isSet(params->predictionData_s,params->predictionData_l);
    bool isOutputSet     = params->isSet(params->output_s,params->output_l);
    bool isForestFileSet = params->isSet(params->forestInput_s,params->forestInput_l);

    bool getAssociations = isFilterSet || isRecombineSet;

    if ( isFilterSet && isRecombineSet ) {
      cerr << "ERROR: cannot have both --filter and --recombine set, choose either one!" << endl;
      exit(1);
    }

    this->filter     =  isFilterSet     &&  isInputSet    && isTargetSet      && isOutputSet;
    this->recombine  =  isRecombineSet  &&  isInputSet    && isOutputSet;
    this->trainModel = !getAssociations &&  isInputSet    && isTargetSet;
    this->testModel  = !getAssociations &&  isPredDataSet && isOutputSet;
    this->saveModel  = !getAssociations &&  isInputSet    && isTargetSet      && !isPredDataSet && isOutputSet;
    this->loadModel  = !getAssociations && !isInputSet    && isForestFileSet  && isOutputSet;

    if ( this->recombine && params->recombinePerms == 0 ) {
      cerr << "Currently the number of permutations to be recombined ( -"
           << params->recombinePerms_s << " / --" << params->recombinePerms_l << endl
           << " ) needs to be explicitly specified." << endl;
      exit(1);
    }

    if ( isFilterSet && !isTargetSet ) {
      cerr << "ERROR: unable to apply filter without a target!" << endl;
      exit(1);
    }

    if ( this->trainModel && !isTargetSet ) {
      cerr << "ERROR: unable to train a model without a target!" << endl;
      exit(1);
    }

    if ( ! ( this->trainModel || this->loadModel ) && ( this->saveModel ) ) {
      cerr << "ERROR: unable to load/train model for saving!" << endl;
      exit(1);
    }

    if ( ! ( this->trainModel || this->loadModel ) && ( this->testModel ) ) {
      cerr << "ERROR: unable to load/train model for testing!" << endl;
      exit(1);
    }

    if ( ! ( this->filter || this->recombine || this->trainModel || this->testModel || this->saveModel || this->loadModel ) ) {
      cerr << "ERROR: unable to resolve execution logic!" << endl;
      params->helpHint();
      exit(1);
    }

  }

  void print() {
    cout << "Parsed program logic:" << endl
	 << "  - filter     == " << this->filter << endl
	 << "  - recombine  == " << this->recombine << endl
	 << "  - trainModel == " << this->trainModel << endl
	 << "  - testModel  == " << this->testModel << endl
	 << "  - saveModel  == " << this->saveModel << endl
	 << "  - loadModel  == " << this->loadModel << endl << endl;
  }

} programLogic;

vector<num_t> readFeatureWeights(Treedata& treeData, const size_t targetIdx, const string& fileName, const num_t featureWeight);

int main(const int argc, char* const argv[]) {

  options::General_options params(argc,argv);

  // With no input arguments the help is printed
  if ( argc == 1 || params.printHelp ) {
    params.help();
    return(EXIT_SUCCESS);
  }

  programLogic.parse(&params);
  programLogic.print();

  RFACE rface(params);

  if ( programLogic.filter ) {

    cout << "===> Reading file '" << params.input << "', please wait... " << flush;
    Treedata treeData(params.input,&params);
    cout << "DONE" << endl;

    size_t targetIdx = treeData.getFeatureIdx(params.targetStr);

    assert( targetIdx != treeData.end() );

    vector<num_t> featureWeights = readFeatureWeights(treeData,targetIdx,params.featureWeightFile,params.defaultFeatureWeight);
    
    rface.filter(treeData,featureWeights);

  } else if ( programLogic.recombine ) {
    
    cout << " *(EXPERIMENTAL) RF-ACE RECOMBINER (" << params.recombinePerms << " permutations) ACTIVATED* " << endl;
        
    rface.recombine();

  } else {

    if ( programLogic.loadModel ) {

      rface.load(params.forestInput);
      
    }

    if ( programLogic.trainModel ) {

      // Read train data into Treedata object
      cout << "===> Reading train file '" << params.input << "', please wait... " << flush;
      Treedata trainData(params.input,&params);
      cout << "DONE" << endl;

      size_t targetIdx = trainData.getFeatureIdx(params.targetStr);

      assert( targetIdx != trainData.end() );
      
      vector<num_t> featureWeights = readFeatureWeights(trainData,targetIdx,params.featureWeightFile,params.defaultFeatureWeight);

      rface.train(trainData,featureWeights);
      
    }

    if ( programLogic.testModel ) {

      cout << "===> Reading test file '" << params.predictionData << "', please wait..." << flush;
      Treedata testData(params.predictionData,&params);
      cout << "DONE" << endl;
      
      cout << "===> Making predictions with test data... " << flush;
      rface.test(testData);
      cout << "DONE" << endl;

      cout << "Prediction file '" << params.output << "' created. Format:" << endl
	   << "TARGET   SAMPLE_ID  TRUE_DATA(*)  PREDICTION    CONFIDENCE(**)" << endl
	   << endl
	   << "  (*): should target variable have true data for test samples, write them," << endl
	   << "       otherwise write NA" << endl
	   << " (**): confidence is the st.dev for regression and % of mispred. for classification" << endl
	   << endl
	   << "RF-ACE completed successfully." << endl
	   << endl;


    }

    if ( programLogic.saveModel ) {

      cout << "===> Writing predictor to file... " << flush;
      rface.save( params.output );
      cout << "DONE" << endl
	   << endl
	   << "RF-ACE predictor built and saved to a file '" << params.output << "'" << endl
	   << endl;

    }
    

  }

  return( EXIT_SUCCESS );

}

vector<num_t> readFeatureWeights(Treedata& treeData, const size_t targetIdx, const string& fileName, const num_t defaulFeatureWeight) {

  vector<num_t> weights(0);

  if ( fileName == "" ) {
    weights.resize(treeData.nFeatures(),1);
  } else {
    
    weights.resize(treeData.nFeatures(),defaulFeatureWeight);
    
    vector<string> weightStrings = utils::readListFromFile(fileName,'\n');
    
    for ( size_t i = 0; i < weightStrings.size(); ++i ) {
      vector<string> weightPair = utils::split(weightStrings[i],'\t');
      string featureName = weightPair[0];
      size_t featureIdx = treeData.getFeatureIdx(featureName);
      
      if ( featureIdx == treeData.end() ) {
	cerr << "Unknown feature name in feature weights: " << featureName << endl;
	exit(1);
      }
      
      weights[featureIdx] = utils::str2<num_t>(weightPair[1]);
      
      cout << "read " << featureName << " => " << weights[featureIdx] << endl;

    }
    
  }

  weights[targetIdx] = 0;

  return(weights);

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
