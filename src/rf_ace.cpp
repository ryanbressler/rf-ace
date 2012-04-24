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
using datadefs::num_t;

statistics::RF_statistics executeRandomForest(Treedata& treeData,
					      const options::General_options& gen_op,
					      vector<num_t>& pValues,
					      vector<num_t>& importanceValues);

void rf_ace_filter(options::General_options& gen_op);

void rf_ace(options::General_options& gen_op);

void printPredictionToFile(StochasticForest& SF, Treedata& treeData, const string& targetName, const string& fileName);

int main(const int argc, char* const argv[]) {

  options::printHeader(cout);

  ArgParse parser(argc,argv);

  // Structs that store all the user-specified command-line arguments
  options::General_options gen_op;
  gen_op.loadUserParams(parser);
  
  // With no input arguments the help is printed
  if ( argc == 1 || gen_op.printHelp ) {
    //options::printFilterOverview();
    gen_op.help();
    return(EXIT_SUCCESS);

  }

  rface::validateRequiredParameters(gen_op);

  gen_op.isFilter ? rf_ace_filter(gen_op) : rf_ace(gen_op) ;
  
}

void rf_ace_filter(options::General_options& gen_op) {
  
  statistics::RF_statistics RF_stat;

  // Read train data into Treedata object
  cout << "Reading file '" << gen_op.input << "', please wait... " << flush;
  Treedata treeData(gen_op.input,gen_op.dataDelimiter,gen_op.headerDelimiter,gen_op.seed);
  cout << "DONE" << endl;

  rface::updateTargetStr(treeData,gen_op);
  rface::pruneFeatureSpace(treeData,gen_op);
  rface::updateMTry(treeData,gen_op);
      
  // After masking, it's safe to refer to features as indices 
  // TODO: rf_ace.cpp: this should be made obsolete; instead of indices, use the feature headers
  size_t targetIdx = treeData.getFeatureIdx(gen_op.targetStr);
    
  if(treeData.nSamples() < 2 * gen_op.nodeSize) {
    cerr << "Not enough samples (" << treeData.nSamples() << ") to perform a single split" << endl;
    exit(1);
  }
  
  rface::printGeneralSetup(treeData,gen_op);
  //rface::printStochasticForestSetup(SF_op);
      
  // Store the start time (in clock cycles) just before the analysis
  clock_t clockStart( clock() );
  
  ////////////////////////////////////////////////////////////////////////
  //  STEP 1 -- MULTIVARIATE ASSOCIATIONS WITH RANDOM FOREST ENSEMBLES  //
  ////////////////////////////////////////////////////////////////////////     
  vector<num_t> pValues; 
  vector<num_t> importanceValues; 
  set<string> featureNames;
  
  cout << "===> Uncovering associations... " << flush;
  RF_stat = executeRandomForest(treeData,gen_op,pValues,importanceValues);
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
     
    // Initialize mapping vector to identity mapping: range 0,1,...,N-1
    // NOTE: when not sorting we use the identity map, otherwise the sorted map
    vector<size_t> featureIcs = utils::range( treeData.nFeatures() );

    // If we prefer sorting the outputs wrt. significance (either p-value or importance)
    if ( !gen_op.noSort ) {
      // If there are more than one permutation, we can compute the p-values and thus sort wrt. them
      if ( gen_op.nPerms > 1 ) {
	bool isIncreasingOrder = true;
	datadefs::sortDataAndMakeRef(isIncreasingOrder,pValues,featureIcs);
	datadefs::sortFromRef<num_t>(importanceValues,featureIcs);
      } else { // ... otherwise we can sort wrt. importance scores 
	bool isIncreasingOrder = false;
	datadefs::sortDataAndMakeRef(isIncreasingOrder,importanceValues,featureIcs);
	datadefs::sortFromRef<num_t>(pValues,featureIcs);
      }
    }
    
    assert( gen_op.targetStr == treeData.getFeatureName(targetIdx) );
    
    for ( size_t i = 0; i < featureIcs.size(); ++i ) {
      size_t featureIdx = featureIcs[i];
      
      if ( gen_op.nPerms > 1 && pValues[i] > gen_op.pValueThreshold ) {
	continue;
      }
      
      if ( featureIdx == targetIdx ) {
	continue;
      }

      ++nSignificantFeatures;
      
      // cout << " " << i << ":" << featureIdx;

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
  
  exit(0);
}

statistics::RF_statistics executeRandomForest(Treedata& treeData,
					      const options::General_options& gen_op,
					      vector<num_t>& pValues,
					      vector<num_t>& importanceValues) {
  
  vector<vector<size_t> > nodeMat(gen_op.nPerms,vector<size_t>(gen_op.nTrees));
  
  vector<vector<num_t> >         importanceMat( gen_op.nPerms, vector<num_t>(treeData.nFeatures()) );
  vector<vector<num_t> > contrastImportanceMat( gen_op.nPerms, vector<num_t>(treeData.nFeatures()) );
  
  size_t nFeatures = treeData.nFeatures();
  
  pValues.clear();
  pValues.resize(nFeatures,1.0);
  importanceValues.resize(2*nFeatures);
  
  StochasticForest::Parameters parameters;
  parameters.model = StochasticForest::RF;
  parameters.inBoxFraction = 1.0;
  parameters.sampleWithReplacement = true;
  parameters.isRandomSplit = true;
  parameters.nTrees       = gen_op.nTrees;
  parameters.mTry         = gen_op.mTry;
  parameters.nMaxLeaves   = gen_op.nMaxLeaves;
  parameters.nodeSize     = gen_op.nodeSize;
  parameters.useContrasts = true;
  parameters.shrinkage    = gen_op.shrinkage;

  if(gen_op.nPerms > 1) {
    parameters.useContrasts = true;
  } else {
    parameters.useContrasts = false;
  }
  
  
  Progress progress;
  clock_t clockStart( clock() );
  vector<num_t> cSample(gen_op.nPerms);
  
  for(int permIdx = 0; permIdx < static_cast<int>(gen_op.nPerms); ++permIdx) {
    
    progress.update(1.0*permIdx/gen_op.nPerms);
    
    // Initialize the Random Forest object
    StochasticForest SF(&treeData,gen_op.targetStr,parameters);
    
    // Get the number of nodes in each tree in the forest
    nodeMat[permIdx] = SF.nNodes();
    
    SF.getImportanceValues(importanceMat[permIdx],contrastImportanceMat[permIdx]);
    
    // Store the new percentile value in the vector cSample
    cSample[permIdx] = math::mean( utils::removeNANs( contrastImportanceMat[permIdx] ) );
    
  }

  if ( gen_op.contrastOutput != "" ) {
    
    utils::printToFile(cSample.begin(),cSample.end(),'\t',gen_op.contrastOutput);

  }
  
  // Remove possible NANs from the contrast sample
  cSample = utils::removeNANs(cSample);
  
  // Notify if the sample size of the null distribution is very low
  if ( gen_op.nPerms > 1 && cSample.size() < 5 ) {
    cerr << " Too few samples drawn ( " << cSample.size() << " < 5 ) from the null distribution. Consider adding more permutations. Quitting..." << endl;
    exit(0);
  }
  
  // Loop through each feature and calculate p-value for each
  for(size_t featureIdx = 0; featureIdx < treeData.nFeatures(); ++featureIdx) {
    
    vector<num_t> fSample(gen_op.nPerms);
    
    // Extract the sample for the real feature
    for(size_t permIdx = 0; permIdx < gen_op.nPerms; ++permIdx) {
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


void rf_ace(options::General_options& gen_op) {
  
  // Read train data into Treedata object
  cout << "Reading file '" << gen_op.input << "', please wait... " << flush;
  Treedata treeData(gen_op.input,gen_op.dataDelimiter,gen_op.headerDelimiter,gen_op.seed);
  cout << "DONE" << endl;
  
  rface::updateTargetStr(treeData,gen_op);

  // Initialize parameters struct for the stochastic forest and load defaults
  StochasticForest::Parameters parameters;
  //if ( PB_op.isGBT ) {
  // parameters.model = StochasticForest::GBT;
  // parameters.inBoxFraction = 0.5;
  // parameters.sampleWithReplacement = false;
  // parameters.isRandomSplit = false;
  // gen_op.setGBTDefaults();
  // else if ( PB_op.isRF ) {
  parameters.model = StochasticForest::RF;
  parameters.inBoxFraction = 1.0;
  parameters.sampleWithReplacement = true;
  parameters.isRandomSplit = true;
  // gen_op.setRFDefaults();
  //} else {
  //   cerr << "Model needs to be specified explicitly" << endl;
  //   exit(1);
  // }

  // These are to override the default parameter settings
  //gen_op.loadUserParams(parser);

  rface::pruneFeatureSpace(treeData,gen_op);
  rface::updateMTry(treeData,gen_op);

  // Copy command line parameters to parameters struct for the stochastic forest
  parameters.nTrees       = gen_op.nTrees;
  parameters.mTry         = gen_op.mTry;
  parameters.nMaxLeaves   = gen_op.nMaxLeaves;
  parameters.nodeSize     = gen_op.nodeSize;
  parameters.useContrasts = false;
  parameters.shrinkage    = gen_op.shrinkage;
    
  rface::printGeneralSetup(treeData,gen_op);
  //rface::printStochasticForestSetup(_op);

  // Store the start time (in clock cycles) just before the analysis
  clock_t clockStart( clock() );
      
  //if ( PB_op.isGBT ) {
  //   cout << "===> Growing GBT predictor... " << flush;
  //} else {
  cout << "===> Growing RF predictor... " << flush;
  // }
  
  StochasticForest SF(&treeData,gen_op.targetStr,parameters);
  cout << "DONE" << endl << endl;

  size_t targetIdx = treeData.getFeatureIdx(gen_op.targetStr);
  vector<num_t> data = utils::removeNANs(treeData.getFeatureData(targetIdx));

  num_t oobError = SF.getOobError();
  num_t ibOobError =  SF.getError();
  num_t nullError;
  if ( treeData.isFeatureNumerical(targetIdx) ) {
    nullError = math::squaredError(data) / data.size();
  } else {
    map<num_t,size_t> freq = math::frequency(data);
    nullError = 1.0 * ( data.size() - freq[ math::mode<num_t>(data) ] ) / data.size();
  }

  cout << "RF training error measures:" << endl;
  cout << "  data variance = " << nullError << endl;
  cout << "      OOB error = " << oobError << endl;
  cout << "   IB+OOB error = " << ibOobError << endl;
  cout << "    1 - OOB/var = " << 1 - oobError / nullError << endl;
  cout << endl;

  if ( gen_op.predictionData == "" ) {

    cout << "===> Writing predictor to file... " << flush;
    SF.printToFile( gen_op.output );
    cout << "DONE" << endl << endl;
    
    cout << endl;
    cout << 1.0 * ( clock() - clockStart ) / CLOCKS_PER_SEC << " seconds elapsed." << endl << endl;
    
    cout << "RF-ACE predictor built and saved to a file." << endl;
    cout << endl;

  } else {

    cout << "===> Making predictions with test data... " << flush;

    printPredictionToFile(SF,treeData,gen_op.targetStr,gen_op.output);

    cout << "DONE" << endl;

    cout << endl;

    cout << 1.0 * ( clock() - clockStart ) / CLOCKS_PER_SEC << " seconds elapsed." << endl << endl;


    cout << "Prediction file '" << gen_op.output << "' created. Format:" << endl;
    cout << "TARGET   SAMPLE_ID     PREDICTION    CONFIDENCE" << endl;
    cout << endl;


    cout << "RF-ACE completed successfully." << endl;
    cout << endl;

  }
  
  exit(0);
}

vector<string> readFeatureMask(const string& fileName);

void printPredictionToFile(StochasticForest& SF, Treedata& treeData, const string& targetName, const string& fileName) {

  ofstream toPredictionFile(fileName.c_str());
  
  //size_t targetIdx = treeData.getFeatureIdx(targetName);
  
  if ( SF.isTargetNumerical() ) {
    
    vector<num_t> prediction;
    vector<num_t> confidence;
    SF.predict(prediction,confidence);
    
    for(size_t i = 0; i < prediction.size(); ++i) {
      toPredictionFile << targetName << "\t" << treeData.getSampleName(i) << "\t" << prediction[i] << "\t" << setprecision(3) << confidence[i] << endl;
    }
    
  } else {

    vector<string> prediction;
    vector<num_t> confidence;
    SF.predict(prediction,confidence);
    
    for(size_t i = 0; i < prediction.size(); ++i) {
      toPredictionFile << targetName << "\t" << treeData.getSampleName(i) << "\t" << prediction[i] << "\t" << setprecision(3) << confidence[i] << endl;
    }
    
  }

  toPredictionFile.close();
   
}


