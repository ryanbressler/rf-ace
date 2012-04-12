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

statistics::RF_statistics executeRandomForest(Treedata& treeData,
					      const options::General_options& gen_op,
					      const options::StochasticForest_options& SF_op,
					      const options::StatisticalTest_options& ST_op,
					      vector<num_t>& pValues,
					      vector<num_t>& importanceValues);

int main(const int argc, char* const argv[]) {
  
  // Structs that store all the user-specified command-line arguments
  options::General_options gen_op; gen_op.loadUserParams(argc,argv);
  options::StochasticForest_options SF_op; SF_op.setRFDefaults(); SF_op.loadUserParams(argc,argv); 
  options::StatisticalTest_options ST_op; ST_op.loadUserParams(argc,argv);
  statistics::RF_statistics RF_stat;
  
  // Print the intro header
  options::printHeader(cout);
  
  // With no input arguments the help is printed
  if(argc == 1 || gen_op.printHelp ) {
    options::printFilterOverview();
    gen_op.help();
    SF_op.help();
    ST_op.help();
    options::printFilterExamples();
    return(EXIT_SUCCESS);
  }

  rface::validateRequiredParameters(gen_op);

  // Read train data into Treedata object
  cout << "Reading file '" << gen_op.input << "', please wait... " << flush;
  Treedata treeData(gen_op.input,gen_op.dataDelimiter,gen_op.headerDelimiter,gen_op.seed);
  cout << "DONE" << endl;

  rface::updateMTry(treeData,SF_op);
  rface::updateTargetStr(treeData,gen_op);
  rface::pruneFeatureSpace(treeData,gen_op);
      
  // After masking, it's safe to refer to features as indices 
  // TODO: rf_ace.cpp: this should be made obsolete; instead of indices, use the feature headers
  size_t targetIdx = treeData.getFeatureIdx(gen_op.targetStr);
    
  if(treeData.nSamples() < 2 * SF_op.nodeSize) {
    cerr << "Not enough samples (" << treeData.nSamples() << ") to perform a single split" << endl;
    return EXIT_FAILURE;
  }
  
  rface::printGeneralSetup(treeData,gen_op);
  rface::printStochasticForestSetup(SF_op);
      
  // Store the start time (in clock cycles) just before the analysis
  clock_t clockStart( clock() );
  
  ////////////////////////////////////////////////////////////////////////
  //  STEP 1 -- MULTIVARIATE ASSOCIATIONS WITH RANDOM FOREST ENSEMBLES  //
  ////////////////////////////////////////////////////////////////////////     
  vector<num_t> pValues; 
  vector<num_t> importanceValues; 
  set<string> featureNames;
  
  cout << "===> Uncovering associations... " << flush;
  RF_stat = executeRandomForest(treeData,gen_op,SF_op,ST_op,pValues,importanceValues);
  cout << "DONE" << endl;
  
  /////////////////////////////////////////////////
  //  STEP 2 -- FEATURE FILTERING WITH P-VALUES  //
  /////////////////////////////////////////////////
  
  cout << "===> Filtering features... " << flush;
  
  size_t nFeatures = treeData.nFeatures();
    
  size_t nSignificantFeatures = 0;

  if( gen_op.output != "" ) {
    
    ofstream toAssociationFile(gen_op.output.c_str());
    toAssociationFile.precision(8);
    
    vector<size_t> refIcs( treeData.nFeatures() );
    bool isIncreasingOrder = true;
    datadefs::sortDataAndMakeRef(isIncreasingOrder,pValues,refIcs);
    assert( pValues.size() == refIcs.size() );
    assert( treeData.nFeatures() == refIcs.size() );
    datadefs::sortFromRef<num_t>(importanceValues,refIcs);
    
    assert( gen_op.targetStr == treeData.getFeatureName(targetIdx) );
    
    for ( size_t i = 0; i < refIcs.size(); ++i ) {
      size_t featureIdx = refIcs[i];
      
      if ( pValues[i] > ST_op.pValueThreshold ) {
	continue;
      }
      
      if ( featureIdx == targetIdx ) {
	continue;
      }

      ++nSignificantFeatures;
      
      toAssociationFile << fixed << gen_op.targetStr.c_str() << "\t" << treeData.getFeatureName(featureIdx).c_str()
			<< "\t" << log10(pValues[i]) << "\t" << importanceValues[i] << "\t"
			<< treeData.pearsonCorrelation(targetIdx,featureIdx) << "\t" << treeData.nRealSamples(targetIdx,featureIdx) << endl;
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

statistics::RF_statistics executeRandomForest(Treedata& treeData,
					      const options::General_options& gen_op,
					      const options::StochasticForest_options& SF_op,
					      const options::StatisticalTest_options& ST_op,
					      vector<num_t>& pValues,
					      vector<num_t>& importanceValues) {
  
  vector<vector<size_t> > nodeMat(ST_op.nPerms,vector<size_t>(SF_op.nTrees));
  
  vector<vector<num_t> >         importanceMat( ST_op.nPerms, vector<num_t>(treeData.nFeatures()) );
  vector<vector<num_t> > contrastImportanceMat( ST_op.nPerms, vector<num_t>(treeData.nFeatures()) );
  
  size_t nFeatures = treeData.nFeatures();
  
  pValues.clear();
  pValues.resize(nFeatures,1.0);
  importanceValues.resize(2*nFeatures);
  
  StochasticForest::Parameters parameters;
  parameters.model = StochasticForest::RF;
  parameters.inBoxFraction = 1.0;
  parameters.sampleWithReplacement = true;
  parameters.isRandomSplit = true;
  parameters.nTrees       = SF_op.nTrees;
  parameters.mTry         = SF_op.mTry;
  parameters.nMaxLeaves   = SF_op.nMaxLeaves;
  parameters.nodeSize     = SF_op.nodeSize;
  parameters.useContrasts = true;
  parameters.shrinkage    = SF_op.shrinkage;

  if(ST_op.nPerms > 1) {
    parameters.useContrasts = true;
  } else {
    parameters.useContrasts = false;
  }


  Progress progress;
  clock_t clockStart( clock() );
  vector<num_t> cSample(ST_op.nPerms);
  
  for(int permIdx = 0; permIdx < static_cast<int>(ST_op.nPerms); ++permIdx) {
    
    progress.update(1.0*permIdx/ST_op.nPerms);
    
    // Initialize the Random Forest object
    StochasticForest SF(&treeData,gen_op.targetStr,parameters);
        
    // Get the number of nodes in each tree in the forest
    nodeMat[permIdx] = SF.nNodes();
        
    SF.getImportanceValues(importanceMat[permIdx],contrastImportanceMat[permIdx]);

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
  for(size_t featureIdx = 0; featureIdx < treeData.nFeatures(); ++featureIdx) {
    
    vector<num_t> fSample(ST_op.nPerms);
    
    // Extract the sample for the real feature
    for(size_t permIdx = 0; permIdx < ST_op.nPerms; ++permIdx) {
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
    importanceValues[featureIdx] = math::mean(fSample);
    
  }
  
  // Store statistics of the run in an object
  statistics::RF_statistics RF_stat(importanceMat,contrastImportanceMat,nodeMat, 1.0 * ( clock() - clockStart ) / CLOCKS_PER_SEC );
  
  // Resize importance value container to proper dimensions
  importanceValues.resize( treeData.nFeatures() );
  
  // Return statistics
  return( RF_stat );
  
}
