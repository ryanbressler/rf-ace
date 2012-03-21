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

  // Structs that store all the user-specified command-line arguments
  options::General_options gen_op(argc,argv);
  options::GBT_options GBT_op(argc,argv);
  options::RF_options RF_op(argc,argv);
  options::PredictorBuilder_options PB_op(argc,argv);

  // Print the intro header
  options::printHeader(cout);

  // With no input arguments the help is printed
  if(argc == 1 || gen_op.printHelp ) {
    options::printPredictorBuilderOverview();
    gen_op.help();
    PB_op.help();
    GBT_op.help();
    RF_op.help();
    PB_op.help();

    return(EXIT_SUCCESS);
  }

  gen_op.validate();
  PB_op.validate();
  
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
    
  // Store the start time (in clock cycles) just before the analysis
  clock_t clockStart( clock() );
        
  if ( PB_op.isGBT ) {
    rface::printGBTSetup(GBT_op);
    cout << "===> Growing GBT predictor... " << flush;  
    StochasticForest SF(&treedata,gen_op.targetStr,GBT_op.nTrees);
    SF.learnGBT(GBT_op.nMaxLeaves, GBT_op.shrinkage, GBT_op.subSampleSize);
    cout << "DONE" << endl;
    cout << "===> Writing predictor to file... " << flush;
    SF.printToFile( gen_op.output );
    cout << "DONE" << endl << endl;
    cout << "OOB error = " << SF.oobError() << endl;
  }

  if ( PB_op.isRF ) {
    rface::printRFSetup(RF_op);
    cout << "===> Growing RF predictor... " << flush;
    bool useContrasts = false;
    StochasticForest SF(&treedata,gen_op.targetStr,RF_op.nTrees);
    SF.learnRF(RF_op.mTryFraction,RF_op.nMaxLeaves,RF_op.nodeSize,useContrasts);
    cout << "DONE" << endl;
    cout << "===> Writing predictor to file... " << flush;
    SF.printToFile( gen_op.output );
    cout << "DONE" << endl << endl;
    cout << "OOB error = " << SF.oobError() << endl;
  }

  cout << endl;
  cout << 1.0 * ( clock() - clockStart ) / CLOCKS_PER_SEC << " seconds elapsed." << endl << endl;
    
  cout << "RF-ACE-BUILD-PREDICTOR completed successfully." << endl;
  cout << endl;
  
  return(EXIT_SUCCESS);
}



