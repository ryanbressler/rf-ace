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
//#include "exceptions.hpp"

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

void printQRFPredictionsToFile(RFACE::QRFPredictionOutput& qPredOut, const bool printDistributions, const string& fileName);

int main(const int argc, char* const argv[]) {

  printHeader(cout);

  Options options(forest_t::QRF);

  try {
    options.load(argc,argv);
  } catch(exception& e) {
    cerr << "Could not parse command-line arguments" << endl;
    return(EXIT_FAILURE);
  }

  Timer timer;

  options.io.validate();

  // With no input arguments the help is printed
  if ( argc == 1 || options.generalOptions.printHelp ) {
    options.help();
    return(EXIT_SUCCESS);
  }

  RFACE rface(options.generalOptions.nThreads,options.generalOptions.seed);

  RFACE::FilterOutput filterOutput;
  RFACE::QRFPredictionOutput qPredOut;

  timer.tic("Total time elapsed");

  if ( options.io.filterDataFile != "" ) {

    options.forestOptions.print();

    options.filterOptions.print();

    bool useContrasts = true;
    cout << "-Reading file '" << options.io.filterDataFile << "' for filtering" << endl;
    Treedata filterData(options.io.filterDataFile,options.generalOptions.dataDelimiter,options.generalOptions.headerDelimiter,useContrasts);

    size_t targetIdx = getTargetIdx(filterData,options.generalOptions.targetStr);

    assert( targetIdx != filterData.end() );

    printDataStatistics(filterData,targetIdx);

    vector<num_t> featureWeights = readFeatureWeights(filterData,targetIdx,options);
    
    if ( options.generalOptions.seed < 0 ) {
      options.generalOptions.seed = distributions::generateSeed();
    }

    filterOutput = rface.filter(&filterData,targetIdx,featureWeights,&options.forestOptions,&options.filterOptions,options.io.saveForestFile);

    options.io.saveForestFile = "";

  } 

  if ( options.io.associationsFile != "" ) {
    cout << "-Writing associations to file '" << options.io.associationsFile << "'" << endl;
    writeFilterOutputToFile(filterOutput,options.io.associationsFile);
  } 

  if ( options.io.loadForestFile != "" && 
       options.io.testDataFile != "" && 
       options.forestOptions.forestType == forest_t::QRF && 
       options.io.predictionsFile != "" ) {

    cout << "-Loading model '" << options.io.loadForestFile << "', making on-the-fly predictions and saving to file '" << options.io.predictionsFile << "'" << endl;
    Treedata testData(options.io.testDataFile,options.generalOptions.dataDelimiter,options.generalOptions.headerDelimiter);
    qPredOut = rface.loadForestAndPredictQRF(options.io.loadForestFile,&testData,options.forestOptions);
    printQRFPredictionsToFile(qPredOut,options.forestOptions.distributions,options.io.predictionsFile);
    return(EXIT_SUCCESS);
  } 

  if ( options.io.loadForestFile != "" ) {
    cout << "-Loading model '" << options.io.loadForestFile << "'" << endl;
    rface.load(options.io.loadForestFile);
  }

  if ( options.io.trainDataFile != "" ) {

    options.forestOptions.print();
    
    // Read train data into Treedata object
    cout << "-Reading train file '" << options.io.trainDataFile << "'" << endl;
    Treedata trainData(options.io.trainDataFile,options.generalOptions.dataDelimiter,options.generalOptions.headerDelimiter);
    
    size_t targetIdx = getTargetIdx(trainData,options.generalOptions.targetStr);
    
    assert( targetIdx != trainData.end() );

    printDataStatistics(trainData,targetIdx);
    
    vector<num_t> featureWeights = readFeatureWeights(trainData,targetIdx,options);
    
    cout << "-Training the model" << endl;
    rface.train(&trainData,targetIdx,featureWeights,&options.forestOptions);

    if ( trainData.feature(targetIdx)->isNumerical() ) {
      vector<num_t> data = utils::removeNANs(trainData.feature(targetIdx)->getNumData());
    
      num_t oobError = 0; //trainedModel->getOobError();
      num_t ibOobError = 0; // trainedModel->getError();
      
      cout << "RF training error measures (NULL == no model):" << endl;
      num_t nullError = math::var(data);
      cout << "              NULL std = " << sqrt( nullError ) << endl;
      cout << "               OOB std = " << sqrt( oobError ) << endl;
      cout << "            IB+OOB std = " << sqrt( ibOobError ) << endl;
      cout << "  % explained by model = " << 1 - oobError / nullError << " = 1 - (OOB var) / (NULL var)" << endl;
    } else {

      vector<cat_t> data = utils::removeNANs(trainData.feature(targetIdx)->getCatData());

      num_t oobError = 0; //trainedModel->getOobError();
      num_t ibOobError = 0; // trainedModel->getError();

      num_t nullError = math::nMismatches( data, math::mode(data) );
      cout << "       NULL % mispred. = " << 1.0 * nullError / data.size() << endl;
      cout << "        OOB % mispred. = " << oobError << endl;
      cout << "     IB+OOB % mispred. = " << ibOobError << endl;
      cout << "  % explained by model = " << 1 - oobError / nullError << " ( 1 - (OOB # mispred.) / (NULL # mispred.) )" << endl;
    }
    cout << endl;
    
  }
  
  if ( options.io.testDataFile != "" ) {  
    cout << "-Reading test file '" << options.io.testDataFile << "'" << endl;
    Treedata testData(options.io.testDataFile,options.generalOptions.dataDelimiter,options.generalOptions.headerDelimiter);
    cout << "-Making predictions" << endl;
    qPredOut = rface.predictQRF(&testData,options.forestOptions);
  }

  if ( options.io.predictionsFile != "" ) {
    cout << "-Writing predictions to file '" << options.io.predictionsFile << "'" << endl; 
    printQRFPredictionsToFile(qPredOut,options.forestOptions.distributions,options.io.predictionsFile);
  }
    

  if ( options.io.saveForestFile != "" ) {
    cout << "-Writing predictor to file '" << options.io.saveForestFile << "'" << endl;
    rface.save( options.io.saveForestFile );    
  }

  timer.toc("Total time elapsed");
  timer.print();
  
  return( EXIT_SUCCESS );
  
}



vector<num_t> readFeatureWeights(const Treedata& treeData, const size_t targetIdx, const Options& options) {

  size_t nFeatures = treeData.nFeatures();
  vector<num_t> weights(0);

  if ( options.io.whiteListFile == "" && options.io.blackListFile == "" && options.io.featureWeightsFile == "" ) {

    weights = treeData.getFeatureWeights();

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

  cout << endl << "Feature '" << treeData.feature(targetIdx)->name() << "' chosen as target with " << nRealSamples << " / " << nSamples << " samples ( " << 100.0 * ( 1 - realFraction ) << " % missing ) among " << nFeatures << " features" << endl;

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

  toAssociationFile << "TARGET\tPREDICTOR\tP_VALUE\tIMPORTANCE\tCORRELATION\tN_SAMPLES" << endl;

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
    toPredictionFile << "SAMPLE_ID\t" << testOutput.targetName << "_TRUE\t" << testOutput.targetName << "_PRED\t" << testOutput.targetName << "_STD" << endl; 
    for(size_t i = 0; i < testOutput.numPredictions.size(); ++i) {
      toPredictionFile << testOutput.sampleNames[i] << "\t"
		       << testOutput.numTrueData[i] << "\t" << testOutput.numPredictions[i] << "\t"
		       << setprecision(3) << testOutput.confidence[i] << endl;
    }

  } else {
    toPredictionFile << "SAMPLE_ID\t" << testOutput.targetName << "_TRUE\t" << testOutput.targetName << "_PRED\t" << testOutput.targetName << "_ERR%" << endl;
    for(size_t i = 0; i < testOutput.catPredictions.size(); ++i) {
      toPredictionFile << testOutput.sampleNames[i] << "\t"
		       << testOutput.catTrueData[i] << "\t" << testOutput.catPredictions[i] << "\t"
		       << setprecision(3) << testOutput.confidence[i] << endl;
    }

  }

  toPredictionFile.close();

}

void printQRFPredictionsToFile(RFACE::QRFPredictionOutput& qPredOut, const bool printDistributions, const string& fileName) {

  ofstream toPredictionFile(fileName.c_str());

  size_t nSamples = qPredOut.sampleNames.size();

  toPredictionFile << "SAMPLE_ID\t" << qPredOut.targetName << "_TRUE";

  if ( qPredOut.isTargetNumerical ) {
    
    // Print quantile headers
    for ( size_t q = 0; q < qPredOut.quantiles.size(); ++q ) {
      toPredictionFile << "\t" << qPredOut.targetName << "_Q" << qPredOut.quantiles[q];
    }

    // If we print also the distributions, print the header for that
    if ( printDistributions ) {
      toPredictionFile << "\t" << qPredOut.targetName << "_DISTRIBUTION";
    }
    toPredictionFile << endl;

    // Loop through all samples
    for ( size_t i = 0; i < nSamples; ++i ) {

      // Print sample name and true numerical data
      toPredictionFile << qPredOut.sampleNames[i] << "\t" << qPredOut.trueNumData[i] << "\t";

      // Print all requested quantiles
      utils::write(toPredictionFile,qPredOut.numPredictions[i].begin(),qPredOut.numPredictions[i].end(),'\t');

      // Print distribution if requested
      if ( printDistributions ) {
	toPredictionFile << "\t";
	utils::write(toPredictionFile,qPredOut.numDistributions[i].begin(),qPredOut.numDistributions[i].end(),',');
      }
      toPredictionFile << endl;
    }

  } else {
    
    // Print category headers
    for ( size_t c = 0; c < qPredOut.categories.size(); ++c ) {
      toPredictionFile << "\t" << qPredOut.targetName << "_" << qPredOut.categories[c];
    }
    toPredictionFile << endl;
    
    // Write category probabilities for all samples
    for ( size_t i = 0; i < nSamples; ++i ) {
      toPredictionFile << qPredOut.sampleNames[i] << "\t" << qPredOut.trueCatData[i] << "\t";
      utils::write(toPredictionFile,qPredOut.catPredictions[i].begin(),qPredOut.catPredictions[i].end(),'\t');
      toPredictionFile << endl;
    }
    
  }
  
}
