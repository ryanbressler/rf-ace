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
#include "progress.hpp"
#include "distributions.hpp"

using namespace std;
using datadefs::num_t;
using datadefs::forest_t;
using datadefs::ftable_t;

using namespace std;

class RFACE { 

public:

  RFACE():
    trainedModel_(NULL) {
  }

  ~RFACE() {

    if ( trainedModel_ ) {
      delete trainedModel_;
      trainedModel_ = NULL;
    }
    
  }

  struct FilterOutput {
    size_t nAllFeatures;
    size_t nSignificantFeatures;
    string targetName;
    vector<string> featureNames;
    vector<num_t> pValues;
    vector<num_t> importances;
    vector<num_t> correlations;
    vector<num_t> sampleCounts;
  };

  struct TestOutput {
    string targetName;
    bool isTargetNumerical;
    vector<num_t> numPredictions;
    vector<num_t> numTrueData;
    vector<string> catPredictions;
    vector<string> catTrueData;
    vector<num_t> confidence;
    vector<string> sampleNames;
  };
  
  void train(Treedata* trainData, 
	     const size_t targetIdx, 
	     const vector<num_t>& featureWeights, 
	     ForestOptions* forestOptions, 
	     const int seed = 0,
	     const size_t nThreads = 1) {
    
    forestOptions->useContrasts = false;

    forestOptions->validate();

    assert( ! forestOptions->useContrasts );

    if ( trainedModel_ ) {
      delete trainedModel_;
      trainedModel_ = NULL;
    }
    
    vector<distributions::Random> randoms = this->makeRandomNumberGenerators(nThreads,seed);
    
    trainedModel_ = new StochasticForest();

    if ( forestOptions->forestType == forest_t::RF ) {
      trainedModel_->learnRF(trainData,targetIdx,forestOptions,featureWeights,randoms);
    } else if ( forestOptions->forestType == forest_t::GBT ) {
      trainedModel_->learnGBT(trainData,targetIdx,forestOptions,featureWeights,randoms);
    } else if ( forestOptions->forestType == forest_t::CART ) {
      trainedModel_->learnRF(trainData,targetIdx,forestOptions,featureWeights,randoms);
    } else {
      cerr << "Unknown forest type!" << endl;
      exit(1);
    }
  }

  FilterOutput filter(Treedata* filterData, 
		      const size_t targetIdx, 
		      const vector<num_t>& featureWeights, 
		      ForestOptions* forestOptions,
		      FilterOptions* filterOptions,
		      const int seed = 0,
		      const size_t nThreads = 1,
		      const string& forestFile = "" ) {

    forestOptions->useContrasts = true;

    forestOptions->validate();
    filterOptions->validate();

    assert( forestOptions->useContrasts );

    if ( filterData->nSamples() < 2 * forestOptions->nodeSize ) {
      cerr << "Not enough samples (" << filterData->nSamples() << ") to perform a single split" << endl;
      exit(1);
    }

    FilterOutput filterOutput;

    vector<num_t> contrastImportanceSample;
    set<size_t> featuresInAllForests;

    vector<distributions::Random> randoms = makeRandomNumberGenerators(nThreads,seed);

    cout << endl << "Uncovering associations... " << flush;
    executeRandomForest(filterData,targetIdx,featureWeights,forestOptions,filterOptions,filterOutput,randoms,forestFile);
    cout << "DONE" << endl;

    return( filterOutput );

  }

  void executeRandomForest(Treedata* filterData,
			   const size_t targetIdx,
			   const vector<num_t>& featureWeights,
			   ForestOptions* forestOptions,
			   FilterOptions* filterOptions,
			   FilterOutput& filterOutput,
			   vector<distributions::Random>& randoms,
			   const string& forestFile) {
    
    //vector<vector<size_t> > nodeMat(filterOptions->nPerms,vector<size_t>(forestOptions->nTrees));

    vector<vector<num_t> >         importanceMat( filterOptions->nPerms, vector<num_t>(filterData->nFeatures()) );
    vector<vector<num_t> > contrastImportanceMat( filterOptions->nPerms, vector<num_t>(filterData->nFeatures()) );

    size_t nFeatures = filterData->nFeatures();

    filterOutput.targetName = filterData->getFeatureName(targetIdx);
    filterOutput.pValues.resize(nFeatures,1.0);
    filterOutput.importances.resize(nFeatures);
    filterOutput.sampleCounts.resize(nFeatures);
    filterOutput.correlations.resize(nFeatures);
    filterOutput.featureNames.resize(nFeatures);

    Progress progress;
    vector<num_t> contrastImportanceSample(filterOptions->nPerms);

    ftable_t frequency;

    ofstream toFile;
    if ( forestFile != "" ) {
      toFile.open(forestFile.c_str());
      toFile.close();
    }

    for(size_t permIdx = 0; permIdx < filterOptions->nPerms; ++permIdx) {

      filterData->permuteContrasts(&randoms[0]);

      progress.update(1.0*permIdx/filterOptions->nPerms);

      StochasticForest SF;

      SF.learnRF(filterData,targetIdx,forestOptions,featureWeights,randoms);

      if ( forestFile != "" ) {
	//ofstream toFile;
	toFile.open(forestFile.c_str(),ios::app);
	SF.saveForest(toFile);
	toFile.close();
      }

      SF.getMeanMinimalDepthValues(filterData,importanceMat[permIdx],contrastImportanceMat[permIdx]);

      // Store the new percentile value in the vector contrastImportanceSample
      contrastImportanceSample[permIdx] = math::mean( utils::removeNANs( contrastImportanceMat[permIdx] ) );

    }

    contrastImportanceSample = utils::removeNANs( contrastImportanceSample );

    //assert( !datadefs::containsNAN(contrastImportanceSample) );

    // Notify if the sample size of the null distribution is very low
    if ( filterOptions->nPerms > 1 && contrastImportanceSample.size() < 5 ) {
      cerr << " WARNING: Too few samples drawn ( " << contrastImportanceSample.size()
	   << " < 5 ) from the null distribution. Consider adding more permutations. Quitting..."
	   << endl;
      exit(1);
    }
    
    // Loop through each feature and calculate p-value for each
    for ( size_t featureIdx = 0; featureIdx < filterData->nFeatures(); ++featureIdx ) {

      vector<num_t> featureImportanceSample(filterOptions->nPerms);
      
      // Extract the sample for the real feature
      for ( size_t permIdx = 0; permIdx < filterOptions->nPerms; ++permIdx ) {
	featureImportanceSample[permIdx] = importanceMat[permIdx][featureIdx];
      }
      
      featureImportanceSample = utils::removeNANs(featureImportanceSample);
      
      assert( !datadefs::containsNAN(featureImportanceSample) );
      
      // If sample size is too small, assign p-value to 1.0
      if ( featureImportanceSample.size() < 5 ) {
	
	filterOutput.pValues[featureIdx] = 1.0;
	
      } else {
	
	// Perform WS-approximated t-test against the contrast sample
	bool WS = true;
	filterOutput.pValues[featureIdx] = math::ttest(contrastImportanceSample,featureImportanceSample,WS);
	
	// If for some reason the t-test returns NAN, turn that into 1.0
	// NOTE: 1.0 is better number than NAN when sorting
	if ( datadefs::isNAN( filterOutput.pValues[featureIdx] ) ) {
	  filterOutput.pValues[featureIdx] = 1.0;
	}
      }
      
      // Calculate mean importace score from the sample
      filterOutput.importances[featureIdx] = math::mean(featureImportanceSample);
      filterOutput.correlations[featureIdx] = filterData->isFeatureNumerical(featureIdx) ? filterData->pearsonCorrelation(targetIdx,featureIdx) : datadefs::NUM_NAN;
      filterOutput.sampleCounts[featureIdx] = filterData->isFeatureNumerical(featureIdx) ? filterData->nRealSamples(targetIdx,featureIdx) : filterData->nRealSamples(targetIdx);
      filterOutput.featureNames[featureIdx] = filterData->getFeatureName(featureIdx);
    }
    
    sortFilterOutput(&filterOutput);
    
    size_t nSelectedFeatures = 0;

    for ( size_t i = 0; i < nFeatures; ++i ) {

      if ( filterOutput.pValues[i] > filterOptions->pValueThreshold ) {
	continue;
      }

      filterOutput.pValues[nSelectedFeatures] = filterOutput.pValues[i];
      filterOutput.importances[nSelectedFeatures] = filterOutput.importances[i];
      filterOutput.correlations[nSelectedFeatures] = filterOutput.correlations[i];
      filterOutput.sampleCounts[nSelectedFeatures] = filterOutput.sampleCounts[i];
      filterOutput.featureNames[nSelectedFeatures] = filterOutput.featureNames[i];
      ++nSelectedFeatures;
    }

    filterOutput.nAllFeatures = nFeatures - 1;
    filterOutput.nSignificantFeatures = nSelectedFeatures;
    filterOutput.pValues.resize(nSelectedFeatures);
    filterOutput.importances.resize(nSelectedFeatures);
    filterOutput.correlations.resize(nSelectedFeatures);
    filterOutput.sampleCounts.resize(nSelectedFeatures);
    filterOutput.featureNames.resize(nSelectedFeatures);

  }

  void sortFilterOutput(FilterOutput* filterOutput) {

    vector<size_t> sortIcs = utils::range(filterOutput->pValues.size());

    bool isIncreasingOrder = true;

    utils::sortDataAndMakeRef(isIncreasingOrder,filterOutput->pValues,sortIcs);
    utils::sortFromRef(filterOutput->importances,sortIcs);
    utils::sortFromRef(filterOutput->correlations,sortIcs);
    utils::sortFromRef(filterOutput->sampleCounts,sortIcs);
    utils::sortFromRef(filterOutput->featureNames,sortIcs);

  }

  TestOutput test(Treedata* testData,const size_t nThreads = 1) {

    assert(trainedModel_);

    TestOutput testOutput;

    size_t targetIdx = testData->getFeatureIdx(trainedModel_->getTargetName());

    testOutput.targetName = trainedModel_->getTargetName();
    vector<num_t> confidence;
    if ( trainedModel_->isTargetNumerical() ) {
      testOutput.isTargetNumerical = true;
      if ( targetIdx != testData->end() ) {
	testOutput.numTrueData = testData->getFeatureData(targetIdx);
      } else {
	testOutput.numTrueData = vector<num_t>(testData->nSamples(),datadefs::NUM_NAN);
      }
      vector<num_t> predictions;
      trainedModel_->predict(testData,predictions,confidence,nThreads);
      testOutput.numPredictions = predictions;
    } else {
      testOutput.isTargetNumerical = false;
      if ( targetIdx != testData->end() ) {
	testOutput.catTrueData = testData->getRawFeatureData(targetIdx);
      } else {
	testOutput.catTrueData = vector<string>(testData->nSamples(),datadefs::STR_NAN);
      }
      vector<string> predictions;
      trainedModel_->predict(testData,predictions,confidence,nThreads);
      testOutput.catPredictions = predictions;
    }
    
    testOutput.sampleNames.resize( testData->nSamples() );
    for ( size_t i = 0; i < testData->nSamples(); ++i ) {
      testOutput.sampleNames[i] = testData->getSampleName(i);
    }
    testOutput.confidence = confidence;
    return( testOutput );
  }
  
  void load(const string& fileName) {
    
    if ( trainedModel_ ) {
      delete trainedModel_;
      trainedModel_ = NULL;
    }
    
    trainedModel_ = new StochasticForest();
    trainedModel_->loadForest(fileName);
  }

  void save(const string& fileName) {

    assert(trainedModel_);
    
    ofstream toFile(fileName);

    trainedModel_->saveForest( toFile );

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

  /*
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
    toFile << trainData->getFeatureName( fTuples[i][0] ) << "\t" << trainData->getFeatureName( fTuples[i][1] ) << "\t" << 100.0 * fTuples[i][2] / ( filterOptions->nPerms * forestOptions->nTrees ) << endl;
    }
    
    toFile.close();
    
    }
  */

  vector<distributions::Random> makeRandomNumberGenerators(const size_t nThreads, int seed) {

    assert( nThreads >= 1 );

    vector<distributions::Random> randoms(nThreads);

    if ( seed < 0 ) {
      cerr << "Invalid random seed (" << seed << ")" << endl;
    }

    for ( size_t threadIdx = 0; threadIdx < nThreads; ++threadIdx ) {
      randoms[threadIdx].seed(seed + threadIdx);
    }

    return(randoms);

  }


private:

  StochasticForest* trainedModel_;

};

#endif
