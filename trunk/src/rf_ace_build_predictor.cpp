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
using namespace options;
using namespace statistics;
using datadefs::num_t;


void printPredictionToFile(StochasticForest& SF, Treedata& treeData, const string& targetName, const string& fileName);

vector<string> readFeatureMask(const string& fileName);

int main(const int argc, char* const argv[]) {

  cout << endl << " * RF-ACE PREDICTOR BUILDER * " << endl;

  // Structs that store all the user-specified command-line arguments
  General_options gen_op(argc,argv);
  GBT_options GBT_op(argc,argv);

  // Print the intro header
  printHeader(cout);

  // With no input arguments the help is printed
  if(argc == 1 || gen_op.printHelp ) {
    gen_op.help();
    GBT_op.help();
    return(EXIT_SUCCESS);
  }

  validateOptions(gen_op);

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

  if ( treedata.nFeatures() == 0 ) {
    cerr << "All features were removed!" << endl;
    exit(1);
  }

  // After masking, it's safe to refer to features as indices 
  // TODO: rf_ace.cpp: this should be made obsolete; instead of indices, use the feature headers
  size_t targetIdx = treedata.getFeatureIdx(gen_op.targetStr);

  size_t nAllFeatures = treedata.nFeatures();
  size_t nRealSamples = treedata.nRealSamples(targetIdx);
  num_t realFraction = 1.0*nRealSamples / treedata.nSamples();

  //Before number crunching, print values of parameters of RF-ACE
  int maxwidth = 17;
  cout << "General configuration:" << endl;
  cout << "    nfeatures" << setw(8) << "" << "= " << nAllFeatures << endl;
  cout << "    nsamples"  << setw(9) << "" << "= " << treedata.nRealSamples(targetIdx) << " / " << treedata.nSamples() << " ( " << 100.0 * ( 1 - realFraction ) << " % missing )" << endl;
  cout << "    tree type" << setw(8) << "" << "= ";
  if(treedata.isFeatureNumerical(targetIdx)) { cout << "Regression CART" << endl; } else { cout << treedata.nCategories(targetIdx) << "-class CART" << endl; }
  cout << "  --" << gen_op.dataDelimiter_l << setw( maxwidth - gen_op.dataDelimiter_l.size() ) << ""
       << "= '" << gen_op.dataDelimiter << "'" << endl;
  cout << "  --" << gen_op.headerDelimiter_l << setw( maxwidth - gen_op.headerDelimiter_l.size() ) << ""
       << "= '" << gen_op.headerDelimiter << "'" << endl;
  cout << "  --" << gen_op.input_l << setw( maxwidth - gen_op.input_l.size() ) << ""
       << "= " << gen_op.input << endl;
  cout << "  --" << gen_op.targetStr_l << setw( maxwidth - gen_op.targetStr_l.size() ) << ""
       << "= " << gen_op.targetStr << " ( index " << targetIdx << " )" << endl;
  cout << "  --" << gen_op.output_l << setw( maxwidth - gen_op.output_l.size() ) << ""
       << "= "; if ( gen_op.output != "" ) { cout << gen_op.output << endl; } else { cout << "NOT SET" << endl; }
  cout << "  --" << gen_op.log_l << setw( maxwidth - gen_op.log_l.size() ) << ""
       << "= "; if( gen_op.log != "" ) { cout << gen_op.log << endl; } else { cout << "NOT SET" << endl; }
  cout << endl;

  cout << "Gradient boosting tree configuration for prediction:" << endl;
  cout << "  --" << GBT_op.nTrees_l << setw( maxwidth - GBT_op.nTrees_l.size() ) << ""
       << "= " << GBT_op.nTrees << endl;
  cout << "  --" << GBT_op.nMaxLeaves_l << setw( maxwidth - GBT_op.nMaxLeaves_l.size() ) << ""
       << "= " << GBT_op.nMaxLeaves << endl;
  cout << "  --" << GBT_op.shrinkage_l << setw( maxwidth - GBT_op.shrinkage_l.size() ) << ""
       << "= " << GBT_op.shrinkage << endl;
  cout << "  --" << GBT_op.subSampleSize_l << setw( maxwidth - GBT_op.subSampleSize_l.size() ) << ""
       << "= " << GBT_op.subSampleSize << endl;
  cout << endl;
    
  //If the target has no real samples, the program will just exit
  if(nRealSamples == 0) {
    cout << "Target has no real samples. Quitting." << endl;
    return EXIT_SUCCESS;
  }

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
  
  
  
  if ( false ) {
    
    cout << "===> Making predictions with test data... " << flush;
    
    Treedata treedata_test(gen_op.input,gen_op.dataDelimiter,gen_op.headerDelimiter);
    
    printPredictionToFile(SF,treedata_test,gen_op.targetStr,gen_op.output);
    
  } else {
    
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


