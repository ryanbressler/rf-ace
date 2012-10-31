#ifndef RF_ACE_HPP
#define RF_ACE_HPP

#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>
#include <cmath>
#include <stdio.h>
#include <iomanip>
#include <algorithm>
#include <map>
#include <unordered_map>
#include <set>
#include <cstdlib>
#include <string>
#include <vector>

#include "stochasticforest.hpp"
#include "treedata.hpp"
#include "options.hpp"
#include "utils.hpp"
#include "math.hpp"
#include "datadefs.hpp"
#include "statistics.hpp"
#include "progress.hpp"
#include "distributions.hpp"
#include "timer.hpp"
// #include "log.h"


using namespace std;
using datadefs::num_t;
using datadefs::ftable_t;

using namespace std;

class RFACE { 

public:

  RFACE(options::General_options& params): trainedModel_(NULL) {
    
    // Structs that store all the user-specified command-line arguments
    params_ = new options::General_options(params);

    timer_ = new Timer();
    
  }

  ~RFACE() {

    timer_->print();
    delete timer_;
    delete params_;

    if ( trainedModel_ ) {
      delete trainedModel_;
      trainedModel_ = NULL;
    }
    
  }

  void printHeader(ostream& out) {
    out << endl
	<< "-----------------------------------------------------------" << endl
	<< "|  RF-ACE version:  1.0.7, Aug 28 2012                    |" << endl
	<< "|    Compile date:  " << __DATE__ << ", " << __TIME__ << "                 |" << endl
	<< "|   Report issues:  code.google.com/p/rf-ace/issues/list  |" << endl
	<< "-----------------------------------------------------------" << endl
	<< endl;
  }

  void pruneFeatureSpace(Treedata& treeData) {

    size_t targetIdx = treeData.getFeatureIdx(params_->targetStr);

    if ( treeData.nRealSamples(targetIdx) == 0 ) {
      cerr << "Target feature '" << params_->targetStr << "' does not have any real samples!" << endl;
      exit(1);
    }
    
    if ( params_->whiteList != "" ) {
      
      cout << "===> Reading whitelist '" << params_->whiteList << "', please wait... " << flush;
      set<string> whiteFeatureNames = utils::readFeatureMask(treeData,params_->whiteList);
      cout << "DONE" << endl;
      cout << "===> Applying feature mask, removing " << treeData.nFeatures() - whiteFeatureNames.size()
	   << " / " << treeData.nFeatures() << " features, please wait... " << flush;
      
      // Add the target feature into the white list, otherwise it may get removed
      whiteFeatureNames.insert(params_->targetStr);

      treeData.whiteList(whiteFeatureNames);
      cout << "DONE" << endl;
    } 

    if ( params_->blackList != "" ) {
      
      cout << "===> Reading blacklist '" << params_->blackList << "', please wait... " << flush;
      set<string> blackFeatureNames = utils::readFeatureMask(treeData,params_->blackList);
      cout << "DONE" << endl;
      cout << "===> Applying blacklist, keeping " << treeData.nFeatures() - blackFeatureNames.size()
	   << " / " << treeData.nFeatures() << " features, please wait... " << flush;

      // Remove the target feature from the black list, otherwise it will get removed
      if ( blackFeatureNames.find(params_->targetStr) != blackFeatureNames.end() ) {
	cout << " Target found in the blacklist -- omitting... " << flush;
	blackFeatureNames.erase(params_->targetStr);
      }
      
      treeData.blackList(blackFeatureNames);
      cout << "DONE" << endl;
    }
    
    if ( params_->pruneFeatures ) {
      
      cout << "===> Pruning features with less than " << params_->pruneFeatures << " real samples... " << flush;
      size_t nFeaturesOld = treeData.nFeatures();
      utils::pruneFeatures(treeData,params_->targetStr,params_->pruneFeatures);
      cout << "DONE, " << nFeaturesOld - treeData.nFeatures() << " features ( "
	   << ( 100.0*(nFeaturesOld - treeData.nFeatures()) / nFeaturesOld ) << "% ) pruned" << endl;
      
    }

    if ( treeData.nFeatures() == 0 ) {
      cout << "All features were removed!" << endl;

      ofstream toLogFile(params_->log.c_str());

      toLogFile << "All features were removed!" << endl;

      toLogFile.close();

      exit(0);

    }
    
  }

  void updateTargetStr(Treedata& treeData) {

    // Check if the target is specified as an index
    int integer;
    if ( datadefs::isInteger(params_->targetStr,integer) ) {
      
      if ( integer < 0 || integer >= static_cast<int>( treeData.nFeatures() ) ) {
	cerr << "Feature index (" << integer << ") must be within bounds 0 ... " << treeData.nFeatures() - 1 << endl;
	exit(1);
      }
      
      // Extract the name of the feature, as upon masking the indices will become rearranged
      params_->targetStr = treeData.getFeatureName(static_cast<size_t>(integer));
      
    }

  }

  void updateMTry(Treedata& treeData) {
    if ( params_->mTry == options::RF_DEFAULT_M_TRY ) {
      params_->mTry = static_cast<size_t>( floor(0.1*treeData.nFeatures()) );
    }
  }

  void printGeneralSetup(Treedata& treeData) {

    // After masking, it's safe to refer to features as indices
    // TODO: rf_ace.cpp: this should be made obsolete; instead of indices, use the feature headers
    size_t targetIdx = treeData.getFeatureIdx(params_->targetStr);

    size_t nAllFeatures = treeData.nFeatures();
    size_t nRealSamples = treeData.nRealSamples(targetIdx);
    num_t realFraction = 1.0*nRealSamples / treeData.nSamples();

    //Before number crunching, print values of parameters of RF-ACE
    cout << endl;
    cout << "Input data:" << endl;
    cout << " - " << nAllFeatures << " features" << endl;
    cout << " - " << treeData.nRealSamples(targetIdx) << " samples / " << treeData.nSamples() << " ( " << 100.0 * ( 1 - realFraction ) << " % missing )" << endl;

  }

  void setEnforcedForestParameters(Treedata& treeData) {

    if ( params_->modelType == options::RF ) {

      size_t targetIdx = treeData.getFeatureIdx(params_->targetStr);

      // Allow trees to grow to maximal depth, if not told otherwise
      params_->setIfNotSet(params_->nMaxLeaves_s,params_->nMaxLeaves_l,params_->nMaxLeaves,treeData.nRealSamples(targetIdx));

      // RF mTry is by default set to 10% of features
      params_->setIfNotSet(params_->mTry_s,params_->mTry_l,params_->mTry,static_cast<size_t>(0.1*treeData.nFeatures()));

      // Minimum mTry is 1
      if ( params_->mTry < 1 ) {
	params_->mTry = 1;
      }

    } else if ( params_->modelType == options::CART ) {

      // In CART mode only one tree is grown
      params_->nTrees = 1;
    }

  }

  void train(Treedata& trainData) {
    
    if ( trainedModel_ ) {
      delete trainedModel_;
      trainedModel_ = NULL;
    }

    updateTargetStr(trainData);
    
    pruneFeatureSpace(trainData);
    
    setEnforcedForestParameters(trainData);
    
    printGeneralSetup(trainData);
    
    params_->print();
    
    params_->validateParameters();
    
    if ( params_->modelType == options::RF ) {
      cout << "===> Growing RF predictor... " << flush;
    } else if ( params_->modelType == options::GBT ) {
      cout << "===> Growing GBT predictor... " << flush;
    } else if ( params_->modelType == options::CART ) {
      cout << "===> Growing CART predictor... " << flush;
    } else {
      cerr << "Unknown forest type!" << endl;
      exit(1);
    }
    
    trainedModel_ = new StochasticForest(&trainData,params_);
    cout << "DONE" << endl << endl;
    
    if ( params_->modelType == options::GBT ) {
      cout << "GBT diagnostics disabled temporarily" << endl << endl;
      return;
    }
    
    size_t targetIdx = trainData.getFeatureIdx(params_->targetStr);
    vector<num_t> data = utils::removeNANs(trainData.getFeatureData(targetIdx));

    num_t oobError = trainedModel_->getOobError();
    num_t ibOobError =  trainedModel_->getError();
    
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

  }

  void filter(Treedata& treeData, const vector<num_t>& weights) {

    statistics::RF_statistics RF_stat;

    this->updateTargetStr(treeData);
    this->pruneFeatureSpace(treeData);

    // Some default and enforced parameter settings for RF, CART, and GBT
    this->setEnforcedForestParameters(treeData);

    if(treeData.nSamples() < 2 * params_->nodeSize) {
      cerr << "Not enough samples (" << treeData.nSamples() << ") to perform a single split" << endl;
      exit(1);
    }

    this->printGeneralSetup(treeData);

    params_->print();

    params_->validateParameters();

    vector<num_t> pValues;
    vector<num_t> importanceValues;
    vector<num_t> contrastImportanceSample;
    set<size_t> featuresInAllForests;

    cout << "===> Uncovering associations... " << flush;
    RF_stat = executeRandomForest(treeData,pValues,importanceValues,contrastImportanceSample,featuresInAllForests);
    cout << "DONE" << endl;

    // Initialize mapping vector to identity mapping: range 0,1,...,N-1
    // NOTE: when not sorting we use the identity map, otherwise the sorted map
    vector<size_t> featureIcs = utils::range( treeData.nFeatures() );

    // If there are more than one permutation, we can compute the p-values and thus sort wrt. them
    if ( params_->nPerms > 1 ) {
      bool isIncreasingOrder = true;
      utils::sortDataAndMakeRef(isIncreasingOrder,pValues,featureIcs);
      utils::sortFromRef<num_t>(importanceValues,featureIcs);

      // Apply p-value adjustment with the Benjamini-Hochberg method if needed
      if ( params_->isAdjustedPValue ) {
	size_t nTests = treeData.nFeatures() - 1;
	math::adjustPValues(pValues,nTests);
      }

    } else { // ... otherwise we can sort wrt. importance scores
      bool isIncreasingOrder = false;
      utils::sortDataAndMakeRef(isIncreasingOrder,importanceValues,featureIcs);
      utils::sortFromRef<num_t>(pValues,featureIcs);
    }

    size_t nIncludedFeatures = 0;

    for ( size_t i = 0; i < treeData.nFeatures(); ++i ) {

      // If we don't want to report all features in the association file...
      if ( !params_->reportAllFeatures ) {

	bool featureNotInForests = featuresInAllForests.find(featureIcs[i]) == featuresInAllForests.end();
	bool tooHighPValue = params_->nPerms > 1 && pValues[i] > params_->pValueThreshold;
	bool tooLowImportance = importanceValues[i] < params_->importanceThreshold;

	if ( featureNotInForests || tooLowImportance || tooHighPValue ) {
	  continue;
	}
      }

      // In any case, we will omit target feature from the association list
      if ( params_->targetStr == treeData.getFeatureName(featureIcs[i]) ) {
	continue;
      }

      pValues[nIncludedFeatures] = pValues[i];
      importanceValues[nIncludedFeatures] = importanceValues[i];
      featureIcs[nIncludedFeatures] = featureIcs[i];
      ++nIncludedFeatures;
    }

    pValues.resize(nIncludedFeatures);
    importanceValues.resize(nIncludedFeatures);
    featureIcs.resize(nIncludedFeatures);

    printAssociationsToFile(treeData,
			    featureIcs,
			    pValues,
			    importanceValues,
			    contrastImportanceSample);

    if ( params_->log != "" ) {

      ofstream toLogFile(params_->log.c_str());
      this->printHeader(toLogFile);
      RF_stat.print(toLogFile);
      toLogFile.close();

      toLogFile.open("contrasts.tsv");
      RF_stat.printContrastImportance(toLogFile);
      toLogFile.close();

    }

    size_t nFeatures = treeData.nFeatures();

    cout << endl
	 << "Significant associations (" << nIncludedFeatures << "/" << nFeatures - 1 << ") written to file '" << params_->output << "'. Format:" << endl
	 << "TARGET   PREDICTOR   P-VALUE   IMPORTANCE   CORRELATION   NSAMPLES" << endl
	 << endl
	 << "RF-ACE completed successfully." << endl
	 << endl;

  }

  statistics::RF_statistics executeRandomForest(Treedata& treeData,
						vector<num_t>& pValues,
						vector<num_t>& importanceValues,
						vector<num_t>& contrastImportanceSample,
						set<size_t>& featuresInAllForests) {

    vector<vector<size_t> > nodeMat(params_->nPerms,vector<size_t>(params_->nTrees));

    vector<vector<num_t> >         importanceMat( params_->nPerms, vector<num_t>(treeData.nFeatures()) );
    vector<vector<num_t> > contrastImportanceMat( params_->nPerms, vector<num_t>(treeData.nFeatures()) );

    size_t nFeatures = treeData.nFeatures();

    pValues.clear();
    pValues.resize(nFeatures,1.0);
    importanceValues.resize(2*nFeatures);
    featuresInAllForests.clear();

    Progress progress;
    clock_t timeStart( time(0) );
    contrastImportanceSample.resize(params_->nPerms);

    //size_t targetIdx = treeData.getFeatureIdx(params_->targetStr);
    //assert( targetIdx != treeData.end() );

    ftable_t frequency;

    timer_->tic("MODEL_BUILD");

    for(size_t permIdx = 0; permIdx < params_->nPerms; ++permIdx) {

      progress.update(1.0*permIdx/params_->nPerms);

      StochasticForest SF(&treeData,params_);

      // Get the number of nodes in each tree in the forest
      for ( size_t treeIdx = 0; treeIdx < SF.nTrees(); ++treeIdx ) {
	nodeMat[permIdx][treeIdx] = SF.nNodes(treeIdx);
      }

      SF.getImportanceValues(importanceMat[permIdx],contrastImportanceMat[permIdx]);

      // Will update featuresInAllForests to contain all features in the current forest
      math::setUnion(featuresInAllForests,SF.getFeaturesInForest());

      // Store the new percentile value in the vector contrastImportanceSample
      contrastImportanceSample[permIdx] = math::mean( utils::removeNANs( contrastImportanceMat[permIdx] ) );

      if ( params_->isSet(params_->pairInteractionOutput_s,params_->pairInteractionOutput_l) ) {
	updateFeatureFrequency(frequency,&SF);
      }

    }

    timer_->toc("MODEL_BUILD");
    timer_->tic("MODEL_TEST");

    assert( !datadefs::containsNAN(contrastImportanceSample) );

    // Notify if the sample size of the null distribution is very low
    if ( params_->nPerms > 1 && contrastImportanceSample.size() < 5 ) {
      cerr << " Too few samples drawn ( " << contrastImportanceSample.size()
	   << " < 5 ) from the null distribution. Consider adding more permutations. Quitting..."
	   << endl;
      exit(0);
    }

    // Loop through each feature and calculate p-value for each
    for(size_t featureIdx = 0; featureIdx < treeData.nFeatures(); ++featureIdx) {

      vector<num_t> featureImportanceSample(params_->nPerms);

      // Extract the sample for the real feature
      for ( size_t permIdx = 0; permIdx < params_->nPerms; ++permIdx ) {
	featureImportanceSample[permIdx] = importanceMat[permIdx][featureIdx];
      }

      featureImportanceSample = utils::removeNANs(featureImportanceSample);

      assert( !datadefs::containsNAN(featureImportanceSample) );

      // If sample size is too small, assign p-value to 1.0
      if ( featureImportanceSample.size() < 5 ) {

	pValues[featureIdx] = 1.0;

      } else {

	// Perform WS-approximated t-test against the contrast sample
	bool WS = true;
	pValues[featureIdx] = math::ttest(contrastImportanceSample,featureImportanceSample,WS);

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
    statistics::RF_statistics RF_stat(importanceMat,contrastImportanceMat,nodeMat, time(0) - timeStart );

    // Resize importance value container to proper dimensions
    importanceValues.resize( treeData.nFeatures() );

    if ( params_->isSet(params_->pairInteractionOutput_s,params_->pairInteractionOutput_l) ) {
      printPairInteractionsToFile(&treeData,frequency,params_->pairInteractionOutput);
    }

    timer_->toc("MODEL_TEST");

    // Return statistics
    return( RF_stat );

  }

  void test(Treedata& testData) {

    assert(trainedModel_);

    printPredictionToFile(*trainedModel_,testData,params_->targetStr,params_->output);
	
  }
  
  void load(const string& file) {
    
    if ( trainedModel_ ) {
      delete trainedModel_;
      trainedModel_ = NULL;
    }
    
    trainedModel_ = new StochasticForest(params_);
  }

  void save(const string& file) {

    assert(trainedModel_);
    
    trainedModel_->printToFile( file );

  }
  
  void printAssociationsToFile(Treedata& treeData,
			       vector<size_t>& featureIcs,
			       vector<num_t>& pValues,
			       vector<num_t>& importanceValues,
			       vector<num_t>& contrastImportanceSample) {



    assert( featureIcs.size() == pValues.size() );
    assert( featureIcs.size() == importanceValues.size() );

    ofstream toAssociationFile(params_->output.c_str());

    size_t targetIdx = treeData.getFeatureIdx(params_->targetStr);

    for ( size_t i = 0; i < featureIcs.size(); ++i ) {

      assert(featureIcs[i] != targetIdx);

      string featureName = treeData.getFeatureName(featureIcs[i]);
      num_t correlation = treeData.pearsonCorrelation(targetIdx,featureIcs[i]);
      size_t sampleCount = treeData.nRealSamples(targetIdx,featureIcs[i]);

      toAssociationFile << params_->targetStr.c_str() << "\t" << featureName
			<< "\t" << scientific << pValues[i] << "\t" << scientific << importanceValues[i] << "\t"
			<< scientific << correlation << "\t" << sampleCount << endl;

    }

    if ( params_->reportContrasts ) {
      for ( size_t i = 0; i < contrastImportanceSample.size(); ++i ) {
	toAssociationFile << params_->targetStr << "\t" << datadefs::CONTRAST
			  << "\t" << scientific << 1.0 << "\t" << scientific << contrastImportanceSample[i] << "\t"
			  << scientific << 0.0 << "\t" << 0 << endl;
      }
    }

    toAssociationFile.close();

  }

  void recombine() {

    // Read all lines from file
    vector<string> associations = utils::readListFromFile(params_->input,'\n');

    // Initialize containers for storing the maps from associated features to
    // variables
    map<string,vector<num_t> > associationMap;
    map<string,num_t> correlationMap;
    map<string,size_t> sampleCountMap;

    // For reference extract the first association ...
    vector<string> association = utils::split(associations[0],'\t');

    // ... and from the first association extract the target feature
    params_->targetStr = association[0];

    // Go through all associations in the list and update the map containers
    for ( size_t i = 0; i < associations.size(); ++i ) {

      // Extract the association from line "i"
      association = utils::split(associations[i],'\t');

      // Extract the feature name from line "i"
      string featureName = association[1];

      // Make sure the target feature is listed as first in each entry in the file
      // that is to be recombined (thus, we know that the entries are related to
      // the same variable)
      assert( params_->targetStr == association[0] );

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
    contrastImportanceSample.resize(params_->recombinePerms,0.0);

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
      importanceSample.resize(params_->recombinePerms,0.0);

      num_t pValue = math::ttest(importanceSample,contrastImportanceSample);
      num_t importanceValue = math::mean(importanceSample);

      bool featureIsContrast = featureName == datadefs::CONTRAST;
      bool tooHighPValue = params_->nPerms > 1 && pValue > params_->pValueThreshold;
      bool tooLowImportance = importanceValue < params_->importanceThreshold;

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

    ofstream toAssociationFile(params_->output.c_str());

    for ( size_t i = 0; i < featureNames.size(); ++i ) {

      toAssociationFile << params_->targetStr.c_str() << "\t" << featureNames[i]
			<< "\t" << scientific << pValues[i] << "\t" << scientific << importanceValues[i] << "\t"
			<< scientific << correlations[i] << "\t" << sampleCounts[i] << endl;

    }

    toAssociationFile.close();


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
      SF.predict(&treeDataTest,prediction,confidence);

      for(size_t i = 0; i < prediction.size(); ++i) {
	toPredictionFile << targetName << "\t" << treeDataTest.getSampleName(i) << "\t"
			 << trueData[i] << "\t" << prediction[i] << "\t"
			 << setprecision(3) << confidence[i] << endl;
      }

    } else {

      vector<string> trueData(treeDataTest.nSamples(),datadefs::STR_NAN);
      if ( targetIdx != treeDataTest.end() ) {
	trueData = treeDataTest.getRawFeatureData(targetIdx);
      }

      vector<string> prediction;
      vector<num_t> confidence;
      SF.predict(&treeDataTest,prediction,confidence);

      for(size_t i = 0; i < prediction.size(); ++i) {
	toPredictionFile << targetName << "\t" << treeDataTest.getSampleName(i) << "\t"
			 << trueData[i] << "\t" << prediction[i] << "\t"
			 << setprecision(3) << confidence[i] << endl;
      }

    }

    toPredictionFile.close();

  }

  void updateFeatureFrequency(ftable_t& frequency, StochasticForest* SF) {

    // We loop through all the trees
    for ( size_t treeIdx = 0; treeIdx < SF->nTrees(); ++treeIdx ) {

      set<size_t> featuresInTree = SF->tree(treeIdx)->getFeaturesInTree();

      for ( set<size_t>::const_iterator it1(featuresInTree.begin() ); it1 != featuresInTree.end(); ++it1 ) {
	set<size_t>::const_iterator it2( it1 );
	++it2;
	while ( it2 != featuresInTree.end() ) {

	  size_t lokey,hikey;
	  if ( *it1 <= *it2 ) {
	    lokey = *it1;
	    hikey = *it2;
	  } else {
	    lokey = *it2;
	    hikey = *it1;
	  }

	  if ( frequency.find(lokey) == frequency.end() || frequency[lokey].find(hikey) == frequency[lokey].end() ) {
	    frequency[lokey][hikey] = 1;
	  } else {
	    ++frequency[lokey][hikey];
	  }

	  ++it2;

	}
      }

    }

  }

  template<typename T, size_t pos, bool isAscending>
  struct SortBy {

    bool operator()(const vector<T>& a, const vector<T>& b) {
      return( isAscending ? a[pos] < b[pos] : a[pos] > b[pos] );
    }

  };

  void printPairInteractionsToFile(Treedata* trainData, ftable_t& frequency, const string& fileName) {

    ofstream toFile(fileName.c_str());

    vector<vector<size_t> > fTuples;

    for ( ftable_t::const_iterator it1( frequency.begin() ); it1 != frequency.end(); ++it1 ) {
      for ( unordered_map<size_t,size_t>::const_iterator it2( it1->second.begin()); it2 != it1->second.end(); ++it2 ) {
	fTuples.push_back( {it1->first,it2->first,it2->second} );
	//toFile << trainData->getFeatureName( it1->first ) << "\t" << trainData->getFeatureName( it2->first ) << "\t" << it2->second << endl;
      }
    }

    SortBy<size_t,2,false> sorter;

    //cout << sorter(fTuples[0],fTuples[1]) << endl;

    sort(fTuples.begin(),fTuples.end(),sorter);

    for ( size_t i = 0; i < fTuples.size(); ++i ) {
      toFile << trainData->getFeatureName( fTuples[i][0] ) << "\t" << trainData->getFeatureName( fTuples[i][1] ) << "\t" << 100.0 * fTuples[i][2] / ( params_->nPerms * params_->nTrees ) << endl;
    }

    toFile.close();

  }


private:

  StochasticForest* trainedModel_;

  Timer* timer_;
  options::General_options* params_;
  
};

#endif
