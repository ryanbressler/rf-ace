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
#include "utils.hpp"
#include "options.hpp"
#include "statistics.hpp"
#include "progress.hpp"
#include "math.hpp"

using namespace std;
//using namespace statistics;
//using namespace options;
using datadefs::num_t;

statistics::RF_statistics executeRandomForest(Treedata& treedata,
					      const options::General_options& gen_op,
					      const options::RF_options& RF_op,
					      vector<num_t>& pValues,
					      vector<num_t>& importanceValues);

int main(const int argc, char* const argv[]) {
  
  // Structs that store all the user-specified command-line arguments
  options::General_options gen_op(argc,argv);
  options::RF_options RF_op(argc,argv); 
  statistics::RF_statistics RF_stat;
  
  // Print the intro header
  options::printHeader(cout);
  
  // With no input arguments the help is printed
  if(argc == 1 || gen_op.printHelp ) {
    options::printFilterOverview();
    gen_op.help();
    RF_op.help();
    return(EXIT_SUCCESS);
  }
  
  gen_op.validate();
  RF_op.validate();

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
      
  // After masking, it's safe to refer to features as indices 
  // TODO: rf_ace.cpp: this should be made obsolete; instead of indices, use the feature headers
  size_t targetIdx = treedata.getFeatureIdx(gen_op.targetStr);
    
  if(treedata.nSamples() < 2 * RF_op.nodeSize) {
    cerr << "Not enough samples (" << treedata.nSamples() << ") to perform a single split" << endl;
    return EXIT_FAILURE;
  }
  
  rface::printGeneralSetup(treedata,gen_op);
  rface::printRFSetup(RF_op);
      
  // Store the start time (in clock cycles) just before the analysis
  clock_t clockStart( clock() );
  
  ////////////////////////////////////////////////////////////////////////
  //  STEP 1 -- MULTIVARIATE ASSOCIATIONS WITH RANDOM FOREST ENSEMBLES  //
  ////////////////////////////////////////////////////////////////////////     
  vector<num_t> pValues; 
  vector<num_t> importanceValues; 
  set<string> featureNames;
  
  cout << "===> Uncovering associations... " << flush;
  RF_stat = executeRandomForest(treedata,gen_op,RF_op,pValues,importanceValues);
  cout << "DONE" << endl;
  
  /////////////////////////////////////////////////
  //  STEP 2 -- FEATURE FILTERING WITH P-VALUES  //
  /////////////////////////////////////////////////
  
  cout << "===> Filtering features... " << flush;
  
  size_t nFeatures = treedata.nFeatures();
    
  size_t nSignificantFeatures = 0;

  if( gen_op.output != "" ) {
    
    ofstream toAssociationFile(gen_op.output.c_str());
    toAssociationFile.precision(8);
    
    vector<size_t> refIcs( treedata.nFeatures() );
    bool isIncreasingOrder = true;
    datadefs::sortDataAndMakeRef(isIncreasingOrder,pValues,refIcs);
    assert( pValues.size() == refIcs.size() );
    assert( treedata.nFeatures() == refIcs.size() );
    datadefs::sortFromRef<num_t>(importanceValues,refIcs);
    
    assert( gen_op.targetStr == treedata.getFeatureName(targetIdx) );
    
    for ( size_t i = 0; i < refIcs.size(); ++i ) {
      size_t featureIdx = refIcs[i];
      
      if ( pValues[i] > RF_op.pValueThreshold ) {
	continue;
      }
      
      if ( featureIdx == targetIdx ) {
	continue;
      }

      ++nSignificantFeatures;
      
      toAssociationFile << fixed << gen_op.targetStr.c_str() << "\t" << treedata.getFeatureName(featureIdx).c_str()
			<< "\t" << log10(pValues[i]) << "\t" << importanceValues[i] << "\t"
			<< treedata.pearsonCorrelation(targetIdx,featureIdx) << "\t" << treedata.nRealSamples(targetIdx,featureIdx) << endl;
    }
    
    toAssociationFile.close();
  }
  
   
  // Print some statistics
  // NOTE: we're subtracting the target from the total head count, that's why we need to subtract by 1
  cout << "DONE, " << nSignificantFeatures << " / " << nFeatures  - 1 << " features ( "
       << 100.0 * nSignificantFeatures / ( nFeatures - 1 ) << " % ) left " << endl;
  
  if ( gen_op.log != "" ) {
    
    ofstream toLogFile(gen_op.log.c_str());
    options::printHeader(toLogFile);
    RF_stat.print(toLogFile);
    toLogFile.close();
    
    toLogFile.open("contrasts.tsv");
    RF_stat.printContrastImportance(toLogFile);
    toLogFile.close();
    
  }
  
  cout << 1.0 * ( clock() - clockStart ) / CLOCKS_PER_SEC << " seconds elapsed." << endl << endl;
  
  cout << "Association file created, format:" << endl;
  cout << "TARGET   PREDICTOR   LOG10(P-VALUE)   IMPORTANCE   CORRELATION   NSAMPLES" << endl;
  cout << endl;
  cout << "RF-ACE completed successfully." << endl;
  cout << endl;
  
  return(EXIT_SUCCESS);
}

statistics::RF_statistics executeRandomForest(Treedata& treedata,
					      const options::General_options& gen_op,
					      const options::RF_options& RF_op,
					      vector<num_t>& pValues,
					      vector<num_t>& importanceValues) {
  
  //RF_statistics RF_stat;
  vector<vector<size_t> > nodeMat(RF_op.nPerms,vector<size_t>(RF_op.nTrees));
  
  vector<vector<num_t> >         importanceMat( RF_op.nPerms, vector<num_t>(treedata.nFeatures()) );
  vector<vector<num_t> > contrastImportanceMat( RF_op.nPerms, vector<num_t>(treedata.nFeatures()) );
  
  size_t nFeatures = treedata.nFeatures();
  
  pValues.clear();
  pValues.resize(nFeatures,1.0);
  importanceValues.resize(2*nFeatures);
  
  Progress progress;
  clock_t clockStart( clock() );
  vector<num_t> cSample(RF_op.nPerms);
  
  for(int permIdx = 0; permIdx < static_cast<int>(RF_op.nPerms); ++permIdx) {
    
    progress.update(1.0*permIdx/RF_op.nPerms);
    
    bool useContrasts;
    if(RF_op.nPerms > 1) {
      useContrasts = true;
    } else {
      useContrasts = false;
    }
    
    // Initialize the Random Forest object
    StochasticForest SF(&treedata,gen_op.targetStr,RF_op.nTrees);
    
    // Grow the Random Forest
    SF.learnRF(RF_op.mTryFraction,RF_op.nMaxLeaves,RF_op.nodeSize,useContrasts);
    
    // Get the number of nodes in each tree in the forest
    nodeMat[permIdx] = SF.nNodes();
        
    importanceMat[permIdx] = SF.importanceValues();
    contrastImportanceMat[permIdx] = SF.contrastImportanceValues();

    // Store the new percentile value in the vector cSample
    cSample[permIdx] = math::percentile( utils::removeNANs( contrastImportanceMat[permIdx] ) , 0.5);
    
  }
  
  // Remove possible NANs from the contrast sample
  cSample = utils::removeNANs(cSample);

  // Notify if the sample size of the null distribution is very low
  if ( cSample.size() < 5 ) {
    cerr << " Too few samples drawn ( " << cSample.size() << " < 5 ) from the null distribution. Consider adding more permutations. Quitting..." << endl;
    exit(0);
  }

  // Loop through each feature and calculate p-value for each
  for(size_t featureIdx = 0; featureIdx < treedata.nFeatures(); ++featureIdx) {
    
    size_t nRealSamples;
    vector<num_t> fSample(RF_op.nPerms);
    
    // Extract the sample for the real feature
    for(size_t permIdx = 0; permIdx < RF_op.nPerms; ++permIdx) {
      fSample[permIdx] = importanceMat[permIdx][featureIdx];
    }

    // Remove missing values
    fSample = utils::removeNANs(fSample);

    // If sample size is too small, assign p-value to 1.0
    if ( fSample.size() < 5 ) {

      pValues[featureIdx] = 1.0;

    } else {
      
      // Perform t-test against the contrast sample
      pValues[featureIdx] = math::ttest(fSample,cSample);
      
    }

    // If for some reason the t-test returns NAN, turn that into 1.0
    // NOTE: 1.0 is better number than NAN when sorting
    if ( datadefs::isNAN( pValues[featureIdx] ) ) {
      pValues[featureIdx] = 1.0;
    }
    
    // Calculate mean importace score from the sample
    datadefs::mean(fSample,importanceValues[featureIdx],nRealSamples);
    
  }
  
  // Store statistics of the run in an object
  statistics::RF_statistics RF_stat(importanceMat,contrastImportanceMat,nodeMat, 1.0 * ( clock() - clockStart ) / CLOCKS_PER_SEC );
  
  // Resize importance value container to proper dimensions
  importanceValues.resize( treedata.nFeatures() );
  
  // Return statistics
  return( RF_stat );
  
}

