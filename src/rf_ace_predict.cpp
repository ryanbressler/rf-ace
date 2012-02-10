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
//#include "argparse.hpp"
#include "stochasticforest.hpp"
#include "treedata.hpp"
#include "datadefs.hpp"
#include "options.hpp"
#include "statistics.hpp"
#include "progress.hpp"
#include "math.hpp"
#include "utils.hpp"

using namespace std;
//using namespace options;
//using namespace statistics;
using datadefs::num_t;


void printPredictionToFile(StochasticForest& SF, Treedata& treeData, const string& targetName, const string& fileName);

vector<string> readFeatureMask(const string& fileName);

int main(const int argc, char* const argv[]) {

  // Structs that store all the user-specified command-line arguments
  options::General_options gen_op(argc,argv);
  options::Predictor_options pred_op(argc,argv);

  // Print the intro header
  options::printHeader(cout);

  // With no input arguments the help is printed
  if(argc == 1 || gen_op.printHelp ) {
    options::printPredictorOverview();
    gen_op.help();
    pred_op.help();
    return(EXIT_SUCCESS);
  }

  // Read data into Treedata object
  cout << "Reading file '" << gen_op.input << "', please wait... " << flush;
  Treedata treedata(gen_op.input,gen_op.dataDelimiter,gen_op.headerDelimiter);
  cout << "DONE" << endl;

  cout << "===> Loading GBT predictor... " << flush;
  StochasticForest SF(&treedata,pred_op.forest);
  gen_op.targetStr = SF.getTargetName();
  cout << "DONE" << endl;

  options::validateOptions(gen_op);


  rface::printGeneralSetup(treedata,gen_op);
    
  // Store the start time (in clock cycles) just before the analysis
  clock_t clockStart( clock() );
  
  cout << "===> Making predictions with test data... " << flush;
    
  printPredictionToFile(SF,treedata,gen_op.targetStr,gen_op.output);
    
  cout << "DONE" << endl;
  
  cout << endl;
    
  cout << 1.0 * ( clock() - clockStart ) / CLOCKS_PER_SEC << " seconds elapsed." << endl << endl;
  
  
  cout << "Prediction file '" << gen_op.output << "' created. Format:" << endl;
  cout << "TARGET   SAMPLE_ID   DATA      PREDICTION   CONFIDENCE" << endl; 
  cout << endl;
  
  
  cout << "RF-ACE completed successfully." << endl;
  cout << endl;
  
  return(EXIT_SUCCESS);
}

void printPredictionToFile(StochasticForest& SF, Treedata& treeData, const string& targetName, const string& fileName) {

  ofstream toPredictionFile(fileName.c_str());
  
  size_t targetIdx = treeData.getFeatureIdx(targetName);
  
  if ( treeData.isFeatureNumerical(targetIdx)) {
    
    vector<num_t> prediction;
    vector<num_t> confidence;
    SF.predict(prediction,confidence);
    
    for(size_t i = 0; i < prediction.size(); ++i) {
      toPredictionFile << targetName << "\t" << treeData.getSampleName(i) << "\t" << treeData.getRawFeatureData(targetIdx,i)
		       << "\t" << prediction[i] << "\t" << setprecision(3) << confidence[i] << endl;
    }
    
  } else {

    vector<string> prediction;
    vector<num_t> confidence;
    SF.predict(prediction,confidence);
    
    for(size_t i = 0; i < prediction.size(); ++i) {
      toPredictionFile << targetName << "\t" << treeData.getSampleName(i) << "\t" << treeData.getRawFeatureData(targetIdx,i)
                       << "\t" << prediction[i] << "\t" << setprecision(3) << confidence[i] << endl;
    }
    
  }

  toPredictionFile.close();
   
}


