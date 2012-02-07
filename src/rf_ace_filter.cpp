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

#include "argparse.hpp"
#include "stochasticforest.hpp"
#include "treedata.hpp"
#include "datadefs.hpp"
#include "utils.hpp"
#include "options.hpp"
#include "statistics.hpp"
#include "progress.hpp"
#include "math.hpp"

using namespace std;
using namespace statistics;
using namespace options;
using datadefs::num_t;

RF_statistics executeRandomForest(Treedata& treedata,
				  const General_options& gen_op,
				  const RF_options& RF_op,
				  vector<num_t>& pValues,
				  vector<num_t>& importanceValues);

int main(const int argc, char* const argv[]) {
  
  cout << endl << " * RF-ACE FILTER * " << endl;
  
  // Structs that store all the user-specified command-line arguments
  General_options gen_op(argc,argv);
  RF_options RF_op(argc,argv); 
  RF_statistics RF_stat;
  
  // Print the intro header
  printHeader(cout);
  
  // With no input arguments the help is printed
  if(argc == 1 || gen_op.printHelp ) {
    gen_op.help();
    RF_op.help();
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
  
  // Perform masking, if requested
  if ( gen_op.whiteList != "" ) {
    
    cout << "Reading whitelist '" << gen_op.whiteList << "', please wait... " << flush;
    set<string> whiteFeatureNames = utils::readFeatureMask(treedata,gen_op.whiteList);
    cout << "DONE" << endl;
    cout << "Applying feature mask, removing " << treedata.nFeatures() - whiteFeatureNames.size() 
	 << " / " << treedata.nFeatures() << " features, please wait... " << flush;
    treedata.whiteList(whiteFeatureNames);
    cout << "DONE" << endl;
  } else if ( gen_op.blackList != "" ) {
    
    cout << "Reading blacklist '" << gen_op.blackList << "', please wait... " << flush;
    set<string> blackFeatureNames = utils::readFeatureMask(treedata,gen_op.blackList);
    cout << "DONE" << endl;
    cout << "Applying blacklist, keeping " << treedata.nFeatures() - blackFeatureNames.size() 
	 << " / " << treedata.nFeatures() << " features, please wait... " << flush;
    treedata.blackList(blackFeatureNames);
    cout << "DONE" << endl;
  }
  
  if ( gen_op.pruneFeatures ) {
    
    cout << "Pruning features with less than " << gen_op.pruneFeatures << " real samples... " << flush;
    size_t nFeaturesOld = treedata.nFeatures();
    utils::pruneFeatures(treedata,gen_op.targetStr,gen_op.pruneFeatures);
    cout << "DONE, " << nFeaturesOld - treedata.nFeatures() << " features ( "
         << ( 100.0*(nFeaturesOld - treedata.nFeatures()) / nFeaturesOld ) << "% ) pruned" << endl;
    
  }
  
  if ( treedata.nFeatures() == 0 ) {
    cerr << "All features were removed!" << endl;
    exit(1);
  }
  
  // After masking, it's safe to refer to features as indices 
  // TODO: rf_ace.cpp: this should be made obsolete; instead of indices, use the feature headers
  size_t targetIdx = treedata.getFeatureIdx(gen_op.targetStr);
  
  //If default mTry is to be used...
  if ( RF_op.mTry == RF_DEFAULT_M_TRY ) {
    RF_op.mTry = static_cast<size_t>( 0.1*static_cast<num_t>(treedata.nFeatures()));
    if ( RF_op.mTry == 0 ) {
      RF_op.mTry = 2;
    }
  }
  
  if(treedata.nFeatures() < RF_op.mTry) {
    cerr << "Not enough features (" << treedata.nFeatures()-1 << ") to test with mtry = "
         << RF_op.mTry << " features per split" << endl;
    return EXIT_FAILURE;
  }
  
  if(treedata.nSamples() < 2 * RF_op.nodeSize) {
    cerr << "Not enough samples (" << treedata.nSamples() << ") to perform a single split" << endl;
    return EXIT_FAILURE;
  }
  
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
  
  cout << "Random Forest configuration:" << endl;
  cout << "  --" << RF_op.nTrees_l << setw( maxwidth - RF_op.nTrees_l.size() ) << ""
       << "= "; if(RF_op.nTrees == 0) { cout << "DEFAULT" << endl; } else { cout << RF_op.nTrees << endl; }
  cout << "  --" << RF_op.mTry_l << setw( maxwidth - RF_op.mTry_l.size() ) << ""
       << "= "; if(RF_op.mTry == 0) { cout << "DEFAULT" << endl; } else { cout << RF_op.mTry << endl; }
  cout << "  --" << RF_op.nMaxLeaves_l << setw( maxwidth - RF_op.nMaxLeaves_l.size() ) << ""
       << "= " << RF_op.nMaxLeaves << endl;
  cout << "  --" << RF_op.nodeSize_l << setw( maxwidth - RF_op.nodeSize_l.size() ) << ""
       << "= "; if(RF_op.nodeSize == 0) { cout << "DEFAULT" << endl; } else { cout << RF_op.nodeSize << endl; }
  cout << endl;
  
  cout << "Significance analysis configuration:" << endl;
  cout << "  --" << RF_op.nPerms_l << setw( maxwidth - RF_op.nPerms_l.size() ) << ""
       << "= " << RF_op.nPerms << endl;
  cout << "    test type" << setw(8) << "" << "= T-test" << endl;
  cout << "  --pthresold" << setw(8) << "" << "= " << RF_op.pValueThreshold << endl;
  cout << endl;
  
  //If the target has no real samples, the program will just exit
  if(nRealSamples == 0) {
    cout << "Target has no real samples. Quitting." << endl;
    return EXIT_SUCCESS;
  }
  
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
  
  // Go through each feature, and keep those having p-value higher than the threshold. 
  // Save the kept and removed features, and remember to accumulate the counter
  for ( size_t featureIdx = 0; featureIdx < nFeatures; ++featureIdx ) {
    
    if ( featureIdx == targetIdx || pValues[featureIdx] <= RF_op.pValueThreshold ) {
      featureNames.insert(treedata.getFeatureName(featureIdx));
      pValues[nSignificantFeatures] = pValues[featureIdx];
      importanceValues[nSignificantFeatures] = importanceValues[featureIdx];
      ++nSignificantFeatures;
    } 
  }
  
  // Resize containers
  treedata.whiteList( featureNames );
  pValues.resize( nSignificantFeatures );
  importanceValues.resize ( nSignificantFeatures );
  
  targetIdx = treedata.getFeatureIdx(gen_op.targetStr);
  assert( gen_op.targetStr == treedata.getFeatureName(targetIdx) );
  
  if( gen_op.output != "" ) {
    
    ofstream toAssociationFile(gen_op.output.c_str());
    toAssociationFile.precision(8);
    
    vector<size_t> refIcs( treedata.nFeatures() );
    bool isIncreasingOrder = true;
    datadefs::sortDataAndMakeRef(isIncreasingOrder,pValues,refIcs);
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
      
      num_t log10p = log10(pValues[i]);
      if ( log10p < -30.0 ) {
	log10p = -30.0;
      }
      
      toAssociationFile << fixed << gen_op.targetStr.c_str() << "\t" << treedata.getFeatureName(featureIdx).c_str()
			<< "\t" << log10p << "\t" << importanceValues[i] << "\t"
			<< treedata.pearsonCorrelation(targetIdx,featureIdx) << "\t" << treedata.nRealSamples(targetIdx,featureIdx) << endl;
    }
    
    toAssociationFile.close();
  }
  
  
  // Print some statistics
  // NOTE: we're subtracting the target from the total head count, that's why we need to subtract by 1
  cout << "DONE, " << treedata.nFeatures() - 1 << " / " << nAllFeatures  - 1 << " features ( "
       << 100.0 * ( treedata.nFeatures() - 1 ) / ( nAllFeatures - 1 ) << " % ) left " << endl;
  
  if ( gen_op.log != "" ) {
    
    ofstream toLogFile(gen_op.log.c_str());
    printHeader(toLogFile);
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

RF_statistics executeRandomForest(Treedata& treedata,
				  const General_options& gen_op,
				  const RF_options& RF_op,
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
    SF.learnRF(RF_op.mTry,RF_op.nMaxLeaves,RF_op.nodeSize,useContrasts);
    
    // Get the number of nodes in each tree in the forest
    nodeMat[permIdx] = SF.nNodes();
    
    // Compute importance scores for real and contrast features
    importanceValues = SF.featureImportance();
    
    // Get importance scores
    copy(importanceValues.begin(), // Origin START
	 importanceValues.begin() + nFeatures, // Origin END
	 importanceMat[permIdx].begin()); // Destination START
    
    // Get contrast importance scores
    copy(importanceValues.begin() + nFeatures, // Origin START
	 importanceValues.begin() + 2*nFeatures, // Origin END
	 contrastImportanceMat[permIdx].begin()); // Destination START
    
    cSample[permIdx] = math::percentile( utils::removeNANs( contrastImportanceMat[permIdx] ) , 0.95);
    
  }
  
  for(size_t featureIdx = 0; featureIdx < treedata.nFeatures(); ++featureIdx) {
    
    size_t nRealSamples;
    vector<num_t> fSample(RF_op.nPerms);
    
    for(size_t permIdx = 0; permIdx < RF_op.nPerms; ++permIdx) {
      fSample[permIdx] = importanceMat[permIdx][featureIdx];
    }
    
    pValues[featureIdx] = datadefs::ttest(fSample,cSample);
    
    if ( datadefs::isNAN( pValues[featureIdx] ) ) {
      pValues[featureIdx] = 1.0;
    }
    
    datadefs::mean(fSample,importanceValues[featureIdx],nRealSamples);
    
  }
  
  RF_statistics RF_stat(importanceMat,contrastImportanceMat,nodeMat, 1.0 * ( clock() - clockStart ) / CLOCKS_PER_SEC );
  
  importanceValues.resize( treedata.nFeatures() );
  
  return( RF_stat );
  
}

