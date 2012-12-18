#include <cstdlib>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <iomanip>

#include "rf_ace.hpp"
#include "datadefs.hpp"
#include "options.hpp"
#include "timer.hpp"

using namespace std;
using datadefs::num_t;

void printHeader(ostream& out) {
  out << endl
      << "-----------------------------------------------------------" << endl
      << "|  RF-ACE version:  1.1.0, Dec 5th 2012                   |" << endl
      << "|    Compile date:  " << __DATE__ << ", " << __TIME__ << "                 |" << endl
      << "|   Report issues:  code.google.com/p/rf-ace/issues/list  |" << endl
      << "-----------------------------------------------------------" << endl
      << endl;
}

size_t getTargetIdx(Treedata& treeData, const string& targetAsStr);

vector<num_t> readFeatureWeights(const Treedata& treeData, const size_t targetIdx, const Options& options);

void printDataStatistics(Treedata& treeData, const size_t targetIdx);

void writeFilterOutputToFile(RFACE::FilterOutput& filterOutput, const string& fileName);

void printPredictionsToFile(RFACE::TestOutput& testOutput, const string& fileName);

int main(const int argc, char* const argv[]) {

  printHeader(cout);

  Options options;
  options.load(argc,argv);

  Timer timer;

  options.io.validate();

  // With no input arguments the help is printed
  if ( argc == 1 || options.generalOptions.printHelp ) {
    options.help();
    return(EXIT_SUCCESS);
  }

  RFACE rface;

  RFACE::FilterOutput filterOutput;
  RFACE::TestOutput testOutput;

  timer.tic("Total time elapsed");

  if ( options.io.filterDataFile != "" ) {

    options.forestOptions.print();

    options.filterOptions.print();

    bool useContrasts = true;
    cout << "Reading file '" << options.io.filterDataFile << "' for filtering, please wait... " << flush;
    Treedata filterData(options.io.filterDataFile,options.generalOptions.dataDelimiter,options.generalOptions.headerDelimiter,useContrasts);
    cout << "DONE" << endl;

    size_t targetIdx = getTargetIdx(filterData,options.generalOptions.targetStr);

    assert( targetIdx != filterData.end() );

    printDataStatistics(filterData,targetIdx);

    vector<num_t> featureWeights = readFeatureWeights(filterData,targetIdx,options);
    
    if ( options.generalOptions.seed < 0 ) {
      options.generalOptions.seed = distributions::generateSeed();
    }

    filterOutput = rface.filter(&filterData,targetIdx,featureWeights,&options.forestOptions,&options.filterOptions,options.generalOptions.seed,options.generalOptions.nThreads,options.io.saveForestFile);

    options.io.saveForestFile = "";

  } 

  if ( options.io.associationsFile != "" ) {
    
    writeFilterOutputToFile(filterOutput,options.io.associationsFile);
    
    cout << endl
         << "Significant associations (" << filterOutput.nSignificantFeatures << "/" << filterOutput.nAllFeatures << ") written to file '" << options.io.associationsFile << "'. Format:" << endl
         << "TARGET   PREDICTOR   P-VALUE   MEAN.MIN.DEPTH   CORRELATION   NSAMPLES" << endl
         << endl
         << "RF-ACE completed successfully." << endl
         << endl;
  } 

  if ( options.io.loadForestFile != "" ) {
    
    rface.load(options.io.loadForestFile);
    
  }


  
  if ( options.io.trainDataFile != "" ) {

    options.forestOptions.print();
    
    // Read train data into Treedata object
    cout << "Reading train file '" << options.io.trainDataFile << "', please wait... " << flush;
    Treedata trainData(options.io.trainDataFile,options.generalOptions.dataDelimiter,options.generalOptions.headerDelimiter);
    cout << "DONE" << endl;
    
    size_t targetIdx = getTargetIdx(trainData,options.generalOptions.targetStr);
    
    assert( targetIdx != trainData.end() );

    printDataStatistics(trainData,targetIdx);
    
    vector<num_t> featureWeights = readFeatureWeights(trainData,targetIdx,options);
    
    if ( options.generalOptions.seed < 0 ) {
      options.generalOptions.seed = distributions::generateSeed();
    }

    rface.train(&trainData,targetIdx,featureWeights,&options.forestOptions,options.generalOptions.seed,options.generalOptions.nThreads);
    
    vector<num_t> data = utils::removeNANs(trainData.getFeatureData(targetIdx));
    
    num_t oobError = 0; //trainedModel->getOobError();
    num_t ibOobError = 0;// trainedModel->getError();
    
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
  
  if ( options.io.testDataFile != "" ) {
    
    cout << "Reading test file '" << options.io.testDataFile << "', please wait..." << flush;
    Treedata testData(options.io.testDataFile,options.generalOptions.dataDelimiter,options.generalOptions.headerDelimiter);
    cout << "DONE" << endl;
    
    cout << "Making predictions with test data... " << flush;
    testOutput = rface.test(&testData,options.generalOptions.nThreads);
    cout << "DONE" << endl;
    
  }

  if ( options.io.predictionsFile != "" ) {
    printPredictionsToFile(testOutput,options.io.predictionsFile);
    
    cout << "Prediction file '" << options.io.predictionsFile << "' created. Format:" << endl
	 << "TARGET   SAMPLE_ID  TRUE_DATA(*)  PREDICTION    CONFIDENCE(**)" << endl
	 << endl
	 << "  (*): should target variable have true data for test samples, write them," << endl
	 << "       otherwise write NA" << endl
	 << " (**): confidence is the st.dev for regression and % of mispred. for classification" << endl
	 << endl
	 << "RF-ACE completed successfully." << endl
	 << endl; 
  }
  
  if ( options.io.saveForestFile != "" ) {
    
    cout << "Writing predictor to file... " << flush;
    rface.save( options.io.saveForestFile );
    cout << "DONE" << endl
	 << endl
	 << "RF-ACE predictor built and saved to a file '" << options.io.saveForestFile << "'" << endl
	 << endl;
    
  }

  timer.toc("Total time elapsed");

  timer.print();
  
  return( EXIT_SUCCESS );
  
}



vector<num_t> readFeatureWeights(const Treedata& treeData, const size_t targetIdx, const Options& options) {

  size_t nFeatures = treeData.nFeatures();
  vector<num_t> weights(0);

  if ( options.io.whiteListFile == "" && options.io.blackListFile == "" && options.io.featureWeightsFile == "" ) {
    weights.resize(nFeatures,1.0);
    for ( size_t i = 0; i < treeData.nFeatures(); ++i ) {
      if ( treeData.isFeatureTextual(i) ) {
	num_t entropy = treeData.getFeatureEntropy(i);
	cout << "Feature '" << treeData.getFeatureName(i) << "' is textual and has sqrt(entropy) " << sqrtf(entropy) << endl;
	weights[i] = sqrtf(entropy);
      }
    }
  }

  if ( options.io.whiteListFile != "" ) {
    weights.resize(nFeatures,0.0);
    vector<string> featureNames = utils::readListFromFile(options.io.whiteListFile,'\n');
    for ( size_t i = 0; i < featureNames.size(); ++i ) {
      size_t featureIdx = treeData.getFeatureIdx(featureNames[i]);
      if ( featureIdx == treeData.end() ) {
	cout << "WARNING: could not locate feature '" << featureNames[i] << "'" << endl;
      } else { 
	weights[featureIdx] = 1.0;
      }
    }
  }

  if ( options.io.blackListFile != "" ) {
    weights.resize(nFeatures,1.0);
    vector<string> featureNames = utils::readListFromFile(options.io.blackListFile,'\n');
    for ( size_t i = 0; i < featureNames.size(); ++i ) {
      size_t featureIdx = treeData.getFeatureIdx(featureNames[i]);
      if ( featureIdx == treeData.end() ) {
        cout << "WARNING: could not locate feature '" << featureNames[i] << "'" << endl;
      } else {
        weights[featureIdx] = 0.0;
      }
    }
  }


  if ( options.io.featureWeightsFile != "" ) {
    
    weights.resize(nFeatures,options.generalOptions.defaultFeatureWeight);
    
    vector<string> weightStrings = utils::readListFromFile(options.io.featureWeightsFile,'\n');
    
    for ( size_t i = 0; i < weightStrings.size(); ++i ) {
      vector<string> weightPair = utils::split(weightStrings[i],'\t');
      string featureName = weightPair[0];
      size_t featureIdx = treeData.getFeatureIdx(featureName);
      if ( featureIdx == treeData.end() ) {
	cout << "Unknown feature name in feature weights: " << featureName << endl;
      } else {
	weights[featureIdx] = utils::str2<num_t>(weightPair[1]);
      }
    }
    
  }

  weights[targetIdx] = 0;

  return(weights);

}

void printDataStatistics(Treedata& treeData, const size_t targetIdx) {

  //cout << endl << "Data statistics:" << endl;

  size_t nSamples = treeData.nSamples();
  size_t nFeatures = treeData.nFeatures();
  size_t nRealSamples = treeData.nRealSamples(targetIdx);
  num_t realFraction = 1.0*nRealSamples / treeData.nSamples();

  cout << endl << "Feature '" << treeData.getFeatureName(targetIdx) << "' chosen as target with " << nRealSamples << " / " << nSamples << " samples ( " << 100.0 * ( 1 - realFraction ) << " % missing ) among " << nFeatures << " features" << endl;

}

size_t getTargetIdx(Treedata& treeData, const string& targetAsStr) {

  // Check if the target is specified as an index
  int integer;
  if ( datadefs::isInteger(targetAsStr,integer) ) {

    if ( integer < 0 || integer >= static_cast<int>( treeData.nFeatures() ) ) {
      cerr << "Feature index (" << integer << ") must be within bounds 0 ... " << treeData.nFeatures() - 1 << endl;
      exit(1);
    }

    // Extract the name of the feature, as upon masking the indices will become rearranged
    return( static_cast<size_t>(integer) );

  } else {
    
    size_t targetIdx = treeData.getFeatureIdx(targetAsStr);

    if ( targetIdx == treeData.end() ) {
      cerr << "Target " << targetAsStr << " not found in data!" << endl;
      exit(1);
    }

    return( targetIdx );

  }

}

void writeFilterOutputToFile(RFACE::FilterOutput& filterOutput, const string& fileName) {

  size_t nFeatures = filterOutput.pValues.size();

  ofstream toAssociationFile(fileName.c_str());

  for ( size_t i = 0; i < nFeatures; ++i ) {

    toAssociationFile.setf(ios::scientific);
    toAssociationFile << filterOutput.targetName << "\t" << filterOutput.featureNames[i] << "\t" 
		      << filterOutput.pValues[i] << "\t" << filterOutput.importances[i] << "\t"
		      << filterOutput.correlations[i] << flush;
    toAssociationFile.unsetf(ios::scientific);
    toAssociationFile << "\t" << filterOutput.sampleCounts[i] << endl;
    
  }
  
  toAssociationFile.close();
  
}


void printPredictionsToFile(RFACE::TestOutput& testOutput, const string& fileName) {

  ofstream toPredictionFile(fileName.c_str());

  if ( testOutput.isTargetNumerical ) {

    for(size_t i = 0; i < testOutput.numPredictions.size(); ++i) {
      toPredictionFile << testOutput.targetName << "\t" << testOutput.sampleNames[i] << "\t"
		       << testOutput.numTrueData[i] << "\t" << testOutput.numPredictions[i] << "\t"
		       << setprecision(3) << testOutput.confidence[i] << endl;
    }

  } else {

    for(size_t i = 0; i < testOutput.catPredictions.size(); ++i) {
      toPredictionFile << testOutput.targetName << "\t" << testOutput.sampleNames[i] << "\t"
		       << testOutput.catTrueData[i] << "\t" << testOutput.catPredictions[i] << "\t"
		       << setprecision(3) << testOutput.confidence[i] << endl;
    }

  }

  toPredictionFile.close();

}
