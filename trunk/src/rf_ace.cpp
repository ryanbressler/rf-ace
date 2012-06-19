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
#include <map>

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
using datadefs::num_t;

statistics::RF_statistics executeRandomForest(Treedata& treeData,
					      options::General_options& gen_op,
					      vector<num_t>& pValues,
					      vector<num_t>& importanceValues,
					      vector<num_t>& contrastImportanceSample,
					      set<size_t>& featuresInAllForests);

void rf_ace_filter(options::General_options& gen_op);

void rf_ace(options::General_options& gen_op);

void rf_ace_recombine(options::General_options& gen_op);

void setEnforcedForestParameters(Treedata& treeData, options::General_options& gen_op); 

void printGeneralSetup(Treedata& treeData, const options::General_options& gen_op); 

void printAssociationsToFile(options::General_options& gen_op, 
			     vector<string>& featureNames,
			     vector<num_t>& pValues,
			     vector<num_t>& importanceValues,
			     vector<num_t>& contrastImportanceSample,
			     vector<num_t>& correlations,
			     vector<size_t>& sampleCounts);


void printPredictionToFile(StochasticForest& SF, Treedata& treeDataTest, const string& targetName, const string& fileName);

void printHeader(ostream& out) {
  out << endl
      << "-----------------------------------------------------------" << endl
      << "|  RF-ACE version:  1.0.5, May 29 2012                    |" << endl
      << "|    Compile date:  " << __DATE__ << ", " << __TIME__ << "                 |" << endl 
      << "|   Report issues:  code.google.com/p/rf-ace/issues/list  |" << endl
      << "-----------------------------------------------------------" << endl
      << endl;
}

int main(const int argc, char* const argv[]) {

  // Store the start time (in clock cycles) just before the analysis
  clock_t clockStart( clock() );
  
  printHeader(cout);
  
  // Structs that store all the user-specified command-line arguments
  options::General_options gen_op(argc,argv);
  gen_op.loadUserParams();
  
  // With no input arguments the help is printed
  if ( argc == 1 || gen_op.printHelp ) {
    gen_op.help();
    return(EXIT_SUCCESS);
  }

  if ( gen_op.isFilter ) {

    rf_ace_filter(gen_op);

  } else if ( gen_op.isSet(gen_op.recombinePerms_s,gen_op.recombinePerms_l) ) {
    
    cout << " *(EXPERIMENTAL) RF-ACE RECOMBINER (" << gen_op.recombinePerms << " permutations) ACTIVATED* " << endl;
        
    if ( gen_op.recombinePerms == 0 ) {
      cerr << "Currently the number of permutations to be recombined ( -" 
	   << gen_op.recombinePerms_s << " / --" << gen_op.recombinePerms_l << endl
	   << " ) needs to be explicitly specified." << endl;
      exit(1);
    }

    rf_ace_recombine(gen_op);

  } else {

    rf_ace(gen_op);

  }

  cout << endl;
  cout << 1.0 * ( clock() - clockStart ) / CLOCKS_PER_SEC << " seconds elapsed." << endl << endl;

}

void rf_ace_filter(options::General_options& gen_op) {
  
  statistics::RF_statistics RF_stat;

  // Read train data into Treedata object
  cout << "===> Reading file '" << gen_op.input << "', please wait... " << flush;
  Treedata treeData(gen_op.input,gen_op.dataDelimiter,gen_op.headerDelimiter,gen_op.seed);
  cout << "DONE" << endl;

  rface::updateTargetStr(treeData,gen_op);
  rface::pruneFeatureSpace(treeData,gen_op);

  // Some default and enforced parameter settings for RF, CART, and GBT
  setEnforcedForestParameters(treeData,gen_op);

  if(treeData.nSamples() < 2 * gen_op.nodeSize) {
    cerr << "Not enough samples (" << treeData.nSamples() << ") to perform a single split" << endl;
    exit(1);
  }
  
  printGeneralSetup(treeData,gen_op);

  gen_op.printParameters();

  gen_op.validateParameters();
      
  vector<num_t> pValues; 
  vector<num_t> importanceValues; 
  vector<num_t> contrastImportanceSample;
  set<size_t> featuresInAllForests;
  
  cout << "===> Uncovering associations... " << flush;
  RF_stat = executeRandomForest(treeData,gen_op,pValues,importanceValues,contrastImportanceSample,featuresInAllForests);
  cout << "DONE" << endl;
  
  cout << "===> Filtering features... " << flush;

  size_t targetIdx = treeData.getFeatureIdx(gen_op.targetStr);
  
  vector<string> featureNames( treeData.nFeatures() );
  vector<num_t> correlations(  treeData.nFeatures() );
  vector<size_t> sampleCounts( treeData.nFeatures() );

  size_t nIncludedFeatures = 0;

  for ( size_t i = 0; i < treeData.nFeatures(); ++i ) {

    // If we don't want to report all features in the association file...
    if ( !gen_op.reportAllFeatures ) {
     
      bool featureNotInForests = featuresInAllForests.find(i) == featuresInAllForests.end();
      bool tooHighPValue = gen_op.nPerms > 1 && pValues[i] > gen_op.pValueThreshold;
      bool tooLowImportance = importanceValues[i] < gen_op.importanceThreshold;
     
      if ( featureNotInForests || tooLowImportance || tooHighPValue ) {
	continue;
      }
    }

    // In any case, we will omit target feature from the association list
    if ( gen_op.targetStr == treeData.getFeatureName(i) ) {
      continue;
    }
    
    pValues[nIncludedFeatures] = pValues[i];
    importanceValues[nIncludedFeatures] = importanceValues[i];
    featureNames[nIncludedFeatures] = treeData.getFeatureName(i);
    correlations[nIncludedFeatures] = treeData.pearsonCorrelation(targetIdx,i);
    sampleCounts[nIncludedFeatures] = treeData.nRealSamples(targetIdx,i);
    ++nIncludedFeatures;
  }
  
  pValues.resize(nIncludedFeatures);
  importanceValues.resize(nIncludedFeatures);
  featureNames.resize(nIncludedFeatures);
  correlations.resize(nIncludedFeatures);
  sampleCounts.resize(nIncludedFeatures);
  
  printAssociationsToFile(gen_op,
			  featureNames,
			  pValues,
			  importanceValues,
			  contrastImportanceSample,
			  correlations,
			  sampleCounts);
  
  if ( gen_op.log != "" ) {
    
    ofstream toLogFile(gen_op.log.c_str());
    printHeader(toLogFile);
    RF_stat.print(toLogFile);
    toLogFile.close();
    
    toLogFile.open("contrasts.tsv");
    RF_stat.printContrastImportance(toLogFile);
    toLogFile.close();
    
  }
  
  cout << "Association file created, format:" << endl;
  cout << "TARGET   PREDICTOR   P-VALUE   IMPORTANCE   CORRELATION   NSAMPLES" << endl;
  cout << endl;
  cout << "RF-ACE completed successfully." << endl;
  cout << endl;
  
}

void printGeneralSetup(Treedata& treeData, const options::General_options& gen_op) {
  
  // After masking, it's safe to refer to features as indices
  // TODO: rf_ace.cpp: this should be made obsolete; instead of indices, use the feature headers
  size_t targetIdx = treeData.getFeatureIdx(gen_op.targetStr);
  
  size_t nAllFeatures = treeData.nFeatures();
  size_t nRealSamples = treeData.nRealSamples(targetIdx);
  num_t realFraction = 1.0*nRealSamples / treeData.nSamples();
  
  //Before number crunching, print values of parameters of RF-ACE
  cout << endl;
  cout << "Input data:" << endl;
  cout << " - " << nAllFeatures << " features" << endl;
  cout << " - " << treeData.nRealSamples(targetIdx) << " samples / " << treeData.nSamples() << " ( " << 100.0 * ( 1 - realFraction ) << " % missing )" << endl;
  
}


statistics::RF_statistics executeRandomForest(Treedata& treeData,
					      options::General_options& gen_op,
					      vector<num_t>& pValues,
					      vector<num_t>& importanceValues,
					      vector<num_t>& contrastImportanceSample,
					      set<size_t>& featuresInAllForests) {
  
  vector<vector<size_t> > nodeMat(gen_op.nPerms,vector<size_t>(gen_op.nTrees));
  
  vector<vector<num_t> >         importanceMat( gen_op.nPerms, vector<num_t>(treeData.nFeatures()) );
  vector<vector<num_t> > contrastImportanceMat( gen_op.nPerms, vector<num_t>(treeData.nFeatures()) );
  
  size_t nFeatures = treeData.nFeatures();
  
  pValues.clear();
  pValues.resize(nFeatures,1.0);
  importanceValues.resize(2*nFeatures);
  featuresInAllForests.clear();
  
  Progress progress;
  clock_t clockStart( clock() );
  contrastImportanceSample.resize(gen_op.nPerms);
  
  size_t targetIdx = treeData.getFeatureIdx(gen_op.targetStr);
  assert( targetIdx != treeData.end() );

  for(size_t permIdx = 0; permIdx < gen_op.nPerms; ++permIdx) {
    
    progress.update(1.0*permIdx/gen_op.nPerms);
    
    StochasticForest SF(&treeData,&gen_op);

    // Get the number of nodes in each tree in the forest
    nodeMat[permIdx] = SF.nNodes();
    
    SF.getImportanceValues(importanceMat[permIdx],contrastImportanceMat[permIdx]);
    
    // Will update featuresInAllForests to contain all features in the current forest
    math::setUnion(featuresInAllForests,SF.getFeaturesInForest());
    
    // Store the new percentile value in the vector contrastImportanceSample
    contrastImportanceSample[permIdx] = math::mean( contrastImportanceMat[permIdx] );
    
  }
  
  assert( !datadefs::containsNAN(contrastImportanceSample) );
  
  // Notify if the sample size of the null distribution is very low
  if ( gen_op.nPerms > 1 && contrastImportanceSample.size() < 5 ) {
    cerr << " Too few samples drawn ( " << contrastImportanceSample.size() 
	 << " < 5 ) from the null distribution. Consider adding more permutations. Quitting..." 
	 << endl;
    exit(0);
  }
  
  // Loop through each feature and calculate p-value for each
  for(size_t featureIdx = 0; featureIdx < treeData.nFeatures(); ++featureIdx) {
    
    vector<num_t> featureImportanceSample(gen_op.nPerms);
    
    // Extract the sample for the real feature
    for(size_t permIdx = 0; permIdx < gen_op.nPerms; ++permIdx) {
      featureImportanceSample[permIdx] = importanceMat[permIdx][featureIdx];
    }
    
    assert( !datadefs::containsNAN(featureImportanceSample) );

    // If sample size is too small, assign p-value to 1.0
    if ( featureImportanceSample.size() < 5 ) {
      
      pValues[featureIdx] = 1.0;

    } else {
      
      // Perform t-test against the contrast sample
      pValues[featureIdx] = math::ttest(featureImportanceSample,contrastImportanceSample);
      
    }

    // If for some reason the t-test returns NAN, turn that into 1.0
    // NOTE: 1.0 is better number than NAN when sorting
    if ( datadefs::isNAN( pValues[featureIdx] ) ) {
      pValues[featureIdx] = 1.0;
    }
    
    // Calculate mean importace score from the sample
    importanceValues[featureIdx] = math::mean(featureImportanceSample);

  }
  
  // Store statistics of the run in an object
  statistics::RF_statistics RF_stat(importanceMat,contrastImportanceMat,nodeMat, 1.0 * ( clock() - clockStart ) / CLOCKS_PER_SEC );
  
  // Resize importance value container to proper dimensions
  importanceValues.resize( treeData.nFeatures() );
  
  // Return statistics
  return( RF_stat );
  
}

StochasticForest rf_ace_build_predictor(Treedata& trainData, options::General_options& gen_op) {

  rface::updateTargetStr(trainData,gen_op);

  rface::pruneFeatureSpace(trainData,gen_op);

  setEnforcedForestParameters(trainData,gen_op);

  // We never want to use contrasts when we are building a predictor
  gen_op.useContrasts = false;

  printGeneralSetup(trainData,gen_op);

  gen_op.printParameters();

  gen_op.validateParameters();

  if ( gen_op.modelType == options::RF ) {
    cout << "===> Growing RF predictor... " << flush;
  } else if ( gen_op.modelType == options::GBT ) {
    cout << "===> Growing GBT predictor... " << flush;
  } else if ( gen_op.modelType == options::CART ) {
    cout << "===> Growing CART predictor... " << flush;
  } else {
    cerr << "Unknown forest type!" << endl;
    exit(1);
  }

  StochasticForest SF(&trainData,&gen_op);
  cout << "DONE" << endl << endl;

  size_t targetIdx = trainData.getFeatureIdx(gen_op.targetStr);
  vector<num_t> data = utils::removeNANs(trainData.getFeatureData(targetIdx));

  num_t oobError = SF.getOobError();
  num_t ibOobError =  SF.getError();

  cout << "RF training error measures (NULL == no model):" << endl;
  if ( trainData.isFeatureNumerical(targetIdx) ) {
    num_t nullError = math::var(data);
    cout << "              NULL std = " << sqrt( nullError ) << endl;
    cout << "               OOB std = " << sqrt( oobError ) << endl;
    cout << "            IB+OOB std = " << sqrt( ibOobError ) << endl;
    cout << "  % explained by model = " << 1 - oobError / nullError << " = 1 - (OOB var) / (NULL var)" << endl;
  } else {
    num_t nullError = math::nMismatches( data, math::mode(data) );
    cout << "       NULL % mispred. = " << 1.0 * nullError / data.size() << endl;
    cout << "        OOB % mispred. = " << oobError << endl;
    cout << "     IB+OOB % mispred. = " << ibOobError << endl;
    cout << "  % explained by model = " << 1 - oobError / nullError << " ( 1 - (OOB # mispred.) / (NULL # mispred.) )" << endl;
  }
  cout << endl;

  return( SF );

}

void rf_ace(options::General_options& gen_op) {
  
  if ( gen_op.isSet(gen_op.input_s,gen_op.input_l) ) {

    // Read train data into Treedata object
    cout << "===> Reading file '" << gen_op.input << "', please wait... " << flush;
    Treedata trainData(gen_op.input,gen_op.dataDelimiter,gen_op.headerDelimiter,gen_op.seed);
    cout << "DONE" << endl;
    
    StochasticForest SF = rf_ace_build_predictor(trainData,gen_op);

    if ( gen_op.isSet(gen_op.predictionData_s,gen_op.predictionData_l) ) {
      
      cout << "===> Making predictions with test data... " << flush;

      Treedata treeDataTest(gen_op.predictionData,gen_op.dataDelimiter,gen_op.headerDelimiter,gen_op.seed);

      printPredictionToFile(SF,treeDataTest,gen_op.targetStr,gen_op.output);

      cout << "DONE" << endl << endl;

      cout << "Prediction file '" << gen_op.output << "' created. Format:" << endl
	   << "TARGET   SAMPLE_ID  TRUE_DATA(*)  PREDICTION    CONFIDENCE(**)" << endl
	   << endl
	   << "  (*): should target variable have true data for test samples, write them," << endl
	   << "       otherwise write NA" << endl
	   << " (**): confidence is the st.dev for regression and % of mispred. for classification" << endl
	   << endl
	   << "RF-ACE completed successfully." << endl
	   << endl;
      
    } else {
      
      cout << "===> Writing predictor to file... " << flush;
      SF.printToFile( gen_op.output );
      cout << "DONE" << endl 
	   << endl
	   << "RF-ACE predictor built and saved to a file '" << gen_op.output << "'" << endl
	   << endl;

    }

    return;

  } 

  if ( !gen_op.isSet(gen_op.forestInput_s,gen_op.forestInput_l) ) {
    cerr << "If no input data is specified, forest file for prediction is assumed" << endl;
    exit(1);
  }

  if ( !gen_op.isSet(gen_op.predictionData_s,gen_op.predictionData_l) ) {
    cerr << "If forest predictor is given as input, prediction data is assumed" << endl;
    exit(1);
  }

  cout << "Yay, selected the right branch!" << endl;

  StochasticForest SF(&gen_op);

  cout << "===> Making predictions with test data... " << flush;

  Treedata treeDataTest(gen_op.predictionData,gen_op.dataDelimiter,gen_op.headerDelimiter,gen_op.seed);

  printPredictionToFile(SF,treeDataTest,gen_op.targetStr,gen_op.output);

  cout << "DONE" << endl << endl;

  cout << "Prediction file '" << gen_op.output << "' created. Format:" << endl;
  cout << "TARGET   SAMPLE_ID     PREDICTION    CONFIDENCE" << endl;
  cout << endl;

  cout << "RF-ACE completed successfully." << endl;
  cout << endl;

}

void setEnforcedForestParameters(Treedata& treeData, options::General_options& gen_op) {

  size_t targetIdx = treeData.getFeatureIdx(gen_op.targetStr);

  // Allow trees to grow to maximal depth, if not told otherwise
  gen_op.setIfNotSet(gen_op.nMaxLeaves_s,gen_op.nMaxLeaves_l,gen_op.nMaxLeaves,treeData.nRealSamples(targetIdx));
  
  if ( gen_op.modelType == options::RF ) {

    // RF mTry is by default set to 10% of features
    gen_op.setIfNotSet(gen_op.mTry_s,gen_op.mTry_l,gen_op.mTry,static_cast<size_t>(0.1*treeData.nFeatures()));

    // Minimum mTry is 1
    if ( gen_op.mTry < 1 ) {
      gen_op.mTry = 1;
    }

  } else if ( gen_op.modelType == options::CART ) {

    // In CART mode only one tree is grown
    gen_op.nTrees = 1;
  }

}

vector<string> readFeatureMask(const string& fileName);


void printAssociationsToFile(options::General_options& gen_op, 
			     vector<string>& featureNames, 
			     vector<num_t>& pValues,
			     vector<num_t>& importanceValues,
			     vector<num_t>& contrastImportanceSample,
			     vector<num_t>& correlations,
			     vector<size_t>& sampleCounts) {

  assert( featureNames.size() == pValues.size() );
  assert( featureNames.size() == importanceValues.size() );
  assert( featureNames.size() == correlations.size() );
  assert( featureNames.size() == sampleCounts.size() );

  ofstream toAssociationFile(gen_op.output.c_str());

  // Initialize mapping vector to identity mapping: range 0,1,...,N-1
  // NOTE: when not sorting we use the identity map, otherwise the sorted map
  vector<size_t> featureIcs = utils::range( featureNames.size() );

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
  datadefs::sortFromRef<string>(featureNames,featureIcs);
  datadefs::sortFromRef<num_t>(correlations,featureIcs);
  datadefs::sortFromRef<size_t>(sampleCounts,featureIcs);

  size_t nFeatures = featureIcs.size();

  for ( size_t i = 0; i < nFeatures; ++i ) {

    num_t pValue = gen_op.isAdjustedPValue ? pValues[i] * nFeatures / ( i + 1 ) : pValues[i];
    pValue = pValue > 1.0 ? 1.0 : pValue;

    toAssociationFile << gen_op.targetStr.c_str() << "\t" << featureNames[i]
                      << "\t" << scientific << pValue << "\t" << scientific << importanceValues[i] << "\t"
		      << scientific << correlations[i] << "\t" << sampleCounts[i] << endl;

  }

  if ( gen_op.reportContrasts ) {
    for ( size_t i = 0; i < contrastImportanceSample.size(); ++i ) {
      toAssociationFile << gen_op.targetStr.c_str() << "\t" << datadefs::CONTRAST
			<< "\t" << scientific << 1.0 << "\t" << scientific << contrastImportanceSample[i] << "\t" 
			<< scientific << 0.0 << "\t" << 0 << endl;
    }
  }

  toAssociationFile.close();


  // Print some statistics
  // NOTE: we're subtracting the target from the total head count, that's why we need to subtract by 1
  cout << "DONE, " << nFeatures << " features printed to file" << endl;

}

void rf_ace_recombine(options::General_options& gen_op) {

  // Read all lines from file
  vector<string> associations = utils::readListFromFile(gen_op.input,'\n');

  // Initialize containers for storing the maps from associated features to 
  // variables
  map<string,vector<num_t> > associationMap;
  map<string,num_t> correlationMap;
  map<string,size_t> sampleCountMap;

  // For reference extract the first association ...
  vector<string> association = utils::split(associations[0],'\t');

  // ... and from the first association extract the target feature
  gen_op.targetStr = association[0];

  // Go through all associations in the list and update the map containers
  for ( size_t i = 0; i < associations.size(); ++i ) {
    
    // Extract the association from line "i"
    association = utils::split(associations[i],'\t');

    // Extract the feature name from line "i"
    string featureName = association[1];

    // Make sure the target feature is listed as first in each entry in the file
    // that is to be recombined (thus, we know that the entries are related to 
    // the same variable)
    assert( gen_op.targetStr == association[0] );
    
    // Extract the importance value ...
    num_t importanceValue = utils::str2<num_t>(association[3]);
    
    // ... and push it back to the container 
    associationMap[featureName].push_back(importanceValue);

    // 
    correlationMap[featureName] = utils::str2<num_t>(association[4]);
    sampleCountMap[featureName] = utils::str2<size_t>(association[5]);
  }

  assert( associationMap.find( datadefs::CONTRAST ) != associationMap.end() );

  // Subtract one because contrast is included
  size_t nRealFeatures = associationMap.size() - 1;

  vector<string> featureNames(nRealFeatures);
  vector<num_t>  pValues(nRealFeatures);
  vector<num_t>  importanceValues(nRealFeatures);
  vector<num_t>  correlations(nRealFeatures);
  vector<size_t> sampleCounts(nRealFeatures);

  vector<num_t> contrastImportanceSample = associationMap[datadefs::CONTRAST];
  contrastImportanceSample.resize(gen_op.recombinePerms,0.0);

  // Notify if the sample size of the null distribution is very low
  if ( contrastImportanceSample.size() < 5 ) {
    cerr << " Too few samples drawn ( " << contrastImportanceSample.size() << " < 5 ) from the null distribution. Consider adding more permutations. Quitting..." << endl;
    exit(0);
  }

  // Keep count of the total number of features for which there are 
  // enough ( >= 5 ) importance values
  size_t nIncludedFeatures = 0;
  
  // Go through all features in the container
  for ( map<string,vector<num_t> >::const_iterator it(associationMap.begin() ); it != associationMap.end(); ++it ) {
    
    // For clarity map the iterator into more representative variable names
    string featureName = it->first;
    vector<num_t> importanceSample = it->second; 
    importanceSample.resize(gen_op.recombinePerms,0.0);
    
    num_t pValue = math::ttest(importanceSample,contrastImportanceSample);
    num_t importanceValue = math::mean(importanceSample);

    bool featureIsContrast = featureName == datadefs::CONTRAST; 
    bool tooHighPValue = gen_op.nPerms > 1 && pValue > gen_op.pValueThreshold;
    bool tooLowImportance = importanceValue < gen_op.importanceThreshold;
    
    if ( tooLowImportance || 
	 tooHighPValue || 
	 featureIsContrast ) {
      continue;
    }
    
    featureNames[nIncludedFeatures] = featureName;
    pValues[nIncludedFeatures] = pValue;
    importanceValues[nIncludedFeatures] = importanceValue;
    correlations[nIncludedFeatures] = correlationMap[featureName];
    sampleCounts[nIncludedFeatures] = sampleCountMap[featureName];
    
    ++nIncludedFeatures;
    
  }
  
  featureNames.resize(nIncludedFeatures);
  pValues.resize(nIncludedFeatures);
  importanceValues.resize(nIncludedFeatures);
  correlations.resize(nIncludedFeatures);
  sampleCounts.resize(nIncludedFeatures);
  
  gen_op.reportContrasts = false;
  
  printAssociationsToFile(gen_op,
			  featureNames,
			  pValues,
			  importanceValues,
			  contrastImportanceSample,
			  correlations,
			  sampleCounts);
  
}

void printPredictionToFile(StochasticForest& SF, Treedata& treeDataTest, const string& targetName, const string& fileName) {

  ofstream toPredictionFile(fileName.c_str());
  
  size_t targetIdx = treeDataTest.getFeatureIdx(targetName);
  
  if ( SF.isTargetNumerical() ) {
        
    vector<num_t> trueData(treeDataTest.nSamples(),datadefs::NUM_NAN);
    if ( targetIdx != treeDataTest.end() ) {
      trueData = treeDataTest.getFeatureData(targetIdx);
    }
    
    vector<num_t> prediction;
    vector<num_t> confidence;
    SF.predictWithTestData(&treeDataTest,prediction,confidence);
    
    for(size_t i = 0; i < prediction.size(); ++i) {
      toPredictionFile << targetName << "\t" << treeDataTest.getSampleName(i) << "\t" << trueData[i] << "\t" << prediction[i] << "\t" << setprecision(3) << confidence[i] << endl;
    }
    
  } else {
    
    vector<string> trueData(treeDataTest.nSamples(),datadefs::STR_NAN);
    if ( targetIdx != treeDataTest.end() ) {
      trueData = treeDataTest.getRawFeatureData(targetIdx);
    }
    
    vector<string> prediction;
    vector<num_t> confidence;
    SF.predictWithTestData(&treeDataTest,prediction,confidence);
    
    for(size_t i = 0; i < prediction.size(); ++i) {
      toPredictionFile << targetName << "\t" << treeDataTest.getSampleName(i) << "\t" << trueData[i] << "\t" << prediction[i] << "\t" << setprecision(3) << confidence[i] << endl;
    }
    
  }
  
  toPredictionFile.close();
  
}


