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
  options::GBT_options GBT_op(argc,argv);

  // Print the intro header
  options::printHeader(cout);

  // With no input arguments the help is printed
  if(argc == 1 || gen_op.printHelp ) {
    gen_op.help();
    GBT_op.help();
    return(EXIT_SUCCESS);
  }

  options::validateOptions(gen_op);

  // Read train data into Treedata object
  cout << "Reading file '" << gen_op.input << "', please wait... " << flush;
  Treedata treedata(gen_op.input,gen_op.dataDelimiter,gen_op.headerDelimiter);
  cout << "DONE" << endl;

  // Check if the target is specified as an index
  int integer;
  if ( datadefs::isInteger(gen_op.targetStr,integer) ) {

    if ( integer < 0 || integer >= static_cast<int>( treedata.nFeatures() ) ) {
      cerr << "Feature index (" << integer << ") must be within bounds 0 ... " << treedata.nFeatures() - 1 << endl;
      return EXIT_FAILURE;
    }

    // Extract the name of the feature, as upon masking the indices will become rearranged
    gen_op.targetStr = treedata.getFeatureName(static_cast<size_t>(integer));

  }

  rface::pruneFeatureSpace(treedata,gen_op);

  rface::printGeneralSetup(treedata,gen_op);

  rface::printGBTSetup(GBT_op);

  // After masking, it's safe to refer to features as indices 
  // TODO: rf_ace.cpp: this should be made obsolete; instead of indices, use the feature headers
  //size_t targetIdx = treedata.getFeatureIdx(gen_op.targetStr);
    
  // Store the start time (in clock cycles) just before the analysis
  clock_t clockStart( clock() );
      
  cout << "===> Growing GBT predictor... " << flush;
  
  StochasticForest SF(&treedata,gen_op.targetStr,GBT_op.nTrees);
  SF.learnGBT(GBT_op.nMaxLeaves, GBT_op.shrinkage, GBT_op.subSampleSize);
  
  cout << "DONE" << endl;
  
  if( gen_op.output != "" ) {
    cout << "===> Writing predictor to file... " << flush;
    SF.printToFile( gen_op.output );
    cout << "DONE" << endl;
  }
  
  
  // THIS CONTENT WILL BE MOVED TO THE PREDICTOR PROGRAM
  if ( false ) {
    
    cout << "===> Making predictions with test data... " << flush;
    
    Treedata treedata_test(gen_op.input,gen_op.dataDelimiter,gen_op.headerDelimiter);
    
    printPredictionToFile(SF,treedata_test,gen_op.targetStr,gen_op.output);
    
  } 

  if ( false ) {
    
    cout << "===> Making predictions with train data... " << flush;
    
    printPredictionToFile(SF,treedata,gen_op.targetStr,gen_op.output);
    
  }
  
  cout << "DONE" << endl;
  
  
  
  
  cout << endl;
  
  /*
    if ( gen_op.log != "" ) {
    
    ofstream toLogFile(gen_op.log.c_str());
    printHeader(toLogFile);
    RF_stat.print(toLogFile);
    toLogFile.close();
    
    toLogFile.open("contrasts.tsv");
    RF_stat.printContrastImportance(toLogFile);
    toLogFile.close();
    
    }
  */
  
  cout << 1.0 * ( clock() - clockStart ) / CLOCKS_PER_SEC << " seconds elapsed." << endl << endl;
  
  if ( false ) {
    cout << "Prediction file '" << gen_op.output << "' created. Format:" << endl;
    cout << "TARGET   SAMPLE_ID   DATA      PREDICTION   CONFIDENCE" << endl; 
    cout << endl;
  }
  
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
    SF.predict(&treeData,prediction,confidence);
    
    for(size_t i = 0; i < prediction.size(); ++i) {
      toPredictionFile << targetName << "\t" << treeData.getSampleName(i) << "\t" << treeData.getRawFeatureData(targetIdx,i)
		       << "\t" << prediction[i] << "\t" << setprecision(3) << confidence[i] << endl;
    }
    
  } else {

    vector<string> prediction;
    vector<num_t> confidence;
    SF.predict(&treeData,prediction,confidence);
    
    for(size_t i = 0; i < prediction.size(); ++i) {
      toPredictionFile << targetName << "\t" << treeData.getSampleName(i) << "\t" << treeData.getRawFeatureData(targetIdx,i)
                       << "\t" << prediction[i] << "\t" << setprecision(3) << confidence[i] << endl;
    }
    
  }

  toPredictionFile.close();
   
}


