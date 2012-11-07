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

using namespace std;
using datadefs::num_t;

void printHeader(ostream& out) {
  out << endl
      << "-----------------------------------------------------------" << endl
      << "|  RF-ACE version:  1.0.7, Aug 28 2012                    |" << endl
      << "|    Compile date:  " << __DATE__ << ", " << __TIME__ << "                 |" << endl
      << "|   Report issues:  code.google.com/p/rf-ace/issues/list  |" << endl
      << "-----------------------------------------------------------" << endl
      << endl;
}

size_t getTargetIdx(Treedata& treeData, const string& targetAsStr);

vector<num_t> readFeatureWeights(Treedata& treeData, const size_t targetIdx, const string& fileName, const num_t featureWeight);

void printDataStatistics(Treedata& treeData, const size_t targetIdx);

void writeFilterOutputToFile(RFACE::FilterOutput& filterOutput, const string& fileName);

void printPredictionsToFile(RFACE::TestOutput& testOutput, const string& fileName);

int main(const int argc, char* const argv[]) {

  printHeader(cout);

  Options options;
  options.load(argc,argv);

  options.print();

  // With no input arguments the help is printed
  if ( argc == 1 || options.generalOptions.printHelp ) {
    options.help();
    return(EXIT_SUCCESS);
  }

  RFACE rface;

  RFACE::FilterOutput filterOutput;
  RFACE::TestOutput testOutput;

  if ( options.io.filterDataFile != "" ) {

    bool useContrasts = true;
    cout << "===> Reading file '" << options.io.filterDataFile << "' for filtering, please wait... " << flush;
    Treedata filterData(options.io.filterDataFile,options.generalOptions.dataDelimiter,options.generalOptions.headerDelimiter,useContrasts);
    cout << "DONE" << endl;

    size_t targetIdx = getTargetIdx(filterData,options.generalOptions.targetStr);

    assert( targetIdx != filterData.end() );

    cout << "Data statistics:" << endl;
    printDataStatistics(filterData,targetIdx);

    vector<num_t> featureWeights = readFeatureWeights(filterData,targetIdx,options.io.featureWeightsFile,options.generalOptions.defaultFeatureWeight);
    
    filterOutput = rface.filter(&filterData,targetIdx,featureWeights,&options.forestOptions,&options.filterOptions,&options.generalOptions);

  } 



  if ( options.io.associationsFile != "" ) {
    
    writeFilterOutputToFile(filterOutput,options.io.associationsFile);
    
    cout << endl
         << "Significant associations (" << filterOutput.nSignificantFeatures << "/" << filterOutput.nAllFeatures << ") written to file '" << options.io.associationsFile << "'. Format:" << endl
         << "TARGET   PREDICTOR   P-VALUE   IMPORTANCE   CORRELATION   NSAMPLES" << endl
         << endl
         << "RF-ACE completed successfully." << endl
         << endl;
  } 



  if ( false ) {
    
    //cout << " *(EXPERIMENTAL) RF-ACE RECOMBINER (" << params.recombinePerms << " permutations) ACTIVATED* " << endl;
        
    //rface.recombine();
  }


  
  if ( options.io.loadForestFile != "" ) {
    
    rface.load(options.io.loadForestFile);
    
  }


  
  if ( options.io.trainDataFile != "" ) {
    
    // Read train data into Treedata object
    cout << "===> Reading train file '" << options.io.trainDataFile << "', please wait... " << flush;
    Treedata trainData(options.io.trainDataFile,options.generalOptions.dataDelimiter,options.generalOptions.headerDelimiter);
    cout << "DONE" << endl;
    
    size_t targetIdx = getTargetIdx(trainData,options.generalOptions.targetStr);
    
    assert( targetIdx != trainData.end() );
    
    vector<num_t> featureWeights = readFeatureWeights(trainData,targetIdx,options.io.featureWeightsFile,options.generalOptions.defaultFeatureWeight);
    
    rface.train(trainData,targetIdx,featureWeights,&options.forestOptions,&options.generalOptions);
    
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
    
    cout << "===> Reading test file '" << options.io.testDataFile << "', please wait..." << flush;
    Treedata testData(options.io.testDataFile,options.generalOptions.dataDelimiter,options.generalOptions.headerDelimiter);
    cout << "DONE" << endl;
    
    cout << "===> Making predictions with test data... " << flush;
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
    
    cout << "===> Writing predictor to file... " << flush;
    rface.save( options.io.saveForestFile );
    cout << "DONE" << endl
	 << endl
	 << "RF-ACE predictor built and saved to a file '" << options.io.saveForestFile << "'" << endl
	 << endl;
    
  }

  return( EXIT_SUCCESS );

}

vector<num_t> readFeatureWeights(Treedata& treeData, const size_t targetIdx, const string& fileName, const num_t defaulFeatureWeight) {

  vector<num_t> weights(0);

  if ( fileName == "" ) {
    weights.resize(treeData.nFeatures(),1);
  } else {
    
    weights.resize(treeData.nFeatures(),defaulFeatureWeight);
    
    vector<string> weightStrings = utils::readListFromFile(fileName,'\n');
    
    for ( size_t i = 0; i < weightStrings.size(); ++i ) {
      vector<string> weightPair = utils::split(weightStrings[i],'\t');
      string featureName = weightPair[0];
      size_t featureIdx = treeData.getFeatureIdx(featureName);
      
      if ( featureIdx == treeData.end() ) {
	cerr << "Unknown feature name in feature weights: " << featureName << endl;
	exit(1);
      }
      
      weights[featureIdx] = utils::str2<num_t>(weightPair[1]);
      
      cout << "read " << featureName << " => " << weights[featureIdx] << endl;

    }
    
  }

  weights[targetIdx] = 0;

  return(weights);

}

void printDataStatistics(Treedata& treeData, const size_t targetIdx) {

  cout << "here" << endl;

  size_t nSamples = treeData.nSamples();
  size_t nFeatures = treeData.nFeatures();
  size_t nRealSamples = treeData.nRealSamples(targetIdx);
  num_t realFraction = 1.0*nRealSamples / treeData.nSamples();

  cout << " - " << nFeatures << " features" << endl;
  cout << " - " << nRealSamples << " samples / " << nSamples << " ( " << 100.0 * ( 1 - realFraction ) << " % missing )" << endl;

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

    cout << "returning..." << endl;

    return( targetIdx );

  }

}

void writeFilterOutputToFile(RFACE::FilterOutput& filterOutput, const string& fileName) {

  size_t nFeatures = filterOutput.pValues.size();

  ofstream toAssociationFile(fileName.c_str());

  for ( size_t i = 0; i < nFeatures; ++i ) {

    toAssociationFile << filterOutput.targetName << "\t" << filterOutput.featureNames[i] << "\t" 
		      << scientific << filterOutput.pValues[i] << "\t" << scientific << filterOutput.importances[i] << "\t"
		      << scientific << filterOutput.correlations[i] << "\t" << filterOutput.sampleCounts[i] << endl;
    
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
