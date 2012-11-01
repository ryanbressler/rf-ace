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

