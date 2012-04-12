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

#include "rf_ace.hpp"
#include "stochasticforest.hpp"
#include "treedata.hpp"
#include "datadefs.hpp"
#include "options.hpp"
#include "statistics.hpp"
#include "progress.hpp"
#include "math.hpp"
#include "utils.hpp"

using namespace std;
using datadefs::num_t;

int main(const int argc, char* const argv[]) {

  // Structs that store all the user-specified command-line parameters
  options::General_options gen_op; 
  options::PredictorBuilder_options PB_op; 
  options::StochasticForest_options SF_op;

  // Load general parameters
  gen_op.loadUserParams(argc,argv);

  // With no input parameters or with help flag raised help is printed
  if(argc == 1 || gen_op.printHelp ) {

    options::printPredictorBuilderOverview();
    gen_op.help();
    PB_op.help();
    SF_op.help();
    options::printPredictorBuilderExamples();

    return(EXIT_SUCCESS);
  }

  rface::validateRequiredParameters(gen_op);

  // Print the intro header
  options::printHeader(cout);

  // Read train data into Treedata object
  cout << "Reading file '" << gen_op.input << "', please wait... " << flush;
  Treedata treeData(gen_op.input,gen_op.dataDelimiter,gen_op.headerDelimiter,gen_op.seed);
  cout << "DONE" << endl;
  
  rface::updateTargetStr(treeData,gen_op);

  // Load predictor builder parameters 
  PB_op.loadUserParams(argc,argv);

  // Initialize parameters struct for the stochastic forest and load defaults
  StochasticForest::Parameters parameters;
  if ( PB_op.isGBT ) {
    parameters.model = StochasticForest::GBT;
    parameters.inBoxFraction = 0.5;
    parameters.sampleWithReplacement = false;
    parameters.isRandomSplit = false;
    SF_op.setGBTDefaults();
  } else if ( PB_op.isRF ) {
    parameters.model = StochasticForest::RF;
    parameters.inBoxFraction = 1.0;
    parameters.sampleWithReplacement = true;
    parameters.isRandomSplit = true;
    SF_op.setRFDefaults();
  } else {
    cerr << "Model needs to be specified explicitly" << endl;
    return( EXIT_FAILURE );
  }

  // These are to override the default parameter settings
  SF_op.loadUserParams(argc,argv);

  rface::updateMTry(treeData,SF_op);

  // Copy command line parameters to parameters struct for the stochastic forest
  parameters.nTrees       = SF_op.nTrees;
  parameters.mTry         = SF_op.mTry;
  parameters.nMaxLeaves   = SF_op.nMaxLeaves;
  parameters.nodeSize     = SF_op.nodeSize;
  parameters.useContrasts = false;
  parameters.shrinkage    = SF_op.shrinkage;
    
  rface::pruneFeatureSpace(treeData,gen_op);
  rface::printGeneralSetup(treeData,gen_op);
  rface::printStochasticForestSetup(SF_op);

  // Store the start time (in clock cycles) just before the analysis
  clock_t clockStart( clock() );
      
  if ( PB_op.isGBT ) {
    cout << "===> Growing GBT predictor... " << flush;
  } else {
    cout << "===> Growing RF predictor... " << flush;
  }
  
  StochasticForest SF(&treeData,gen_op.targetStr,parameters);
  cout << "DONE" << endl;

  cout << "===> Writing predictor to file... " << flush;
  SF.printToFile( gen_op.output );
  cout << "DONE" << endl << endl;

  cout << "OOB error = " << SF.getOobError();

  size_t targetIdx = treeData.getFeatureIdx(gen_op.targetStr);
  vector<num_t> data = utils::removeNANs(treeData.getFeatureData(targetIdx));
  if ( treeData.isFeatureNumerical(targetIdx) ) {
    cout << " TOTAL error = " << math::squaredError(data) / data.size() << endl;
  } else {
    map<num_t,size_t> freq = math::frequency(data);
    cout << " TOTAL error = " << 1.0 * ( data.size() - freq[ math::mode<num_t>(data) ] ) / data.size() << endl;
  }
  
  cout << endl;
  cout << 1.0 * ( clock() - clockStart ) / CLOCKS_PER_SEC << " seconds elapsed." << endl << endl;
  
  cout << "RF-ACE-BUILD-PREDICTOR completed successfully." << endl;
  cout << endl;
  
  return(EXIT_SUCCESS);
}



