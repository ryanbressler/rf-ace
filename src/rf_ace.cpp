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
using namespace options;
using namespace statistics;
using datadefs::num_t;

void printHeader(ostream& output) {

  output << endl;
  output << " ------------------------------------------------------- " << endl;
  output << "|  RF-ACE version:  0.9.9, January 20th, 2012           |" << endl;
  output << "|    Project page:  http://code.google.com/p/rf-ace     |" << endl;
  output << "|     Report bugs:  timo.p.erkkila@tut.fi               |" << endl;                     
  output << " ------------------------------------------------------- " << endl;
  output << endl;
  
}

void printHelpHint() {
  cout << endl;
  cout << "To get started, type \"-h\" or \"--help\"" << endl;
}

void pruneFeatures(Treedata& treeData, const General_options& gen_op);

RF_statistics executeRandomForest(Treedata& treedata,
				  const RF_options& RF_op,
				  const General_options& gen_op,
				  vector<num_t>& pValues,
				  vector<num_t>& importanceValues);

void printPredictionToFile(StochasticForest& SF, Treedata& treeData, const string& targetName, const string& fileName);

vector<string> readFeatureMask(const string& fileName);


int main(const int argc, char* const argv[]) {

  // Structs that store all the user-specified command-line arguments
  General_options gen_op;
  RF_options RF_op; 
  RF_statistics RF_stat;

  GBT_options GBT_op;

  // Print the intro header
  printHeader(cout);

  // With no input arguments the help is printed
  if(argc == 1) {
    printHelp(gen_op,RF_op,GBT_op);
    return(EXIT_SUCCESS);
  }

  // Read the user parameters ... 
  ArgParse parser(argc,argv);
  
  // First read general options
  parser.getFlag(gen_op.printHelp_s, gen_op.printHelp_l, gen_op.printHelp);
  parser.getArgument<string>(gen_op.trainInput_s, gen_op.trainInput_l, gen_op.trainInput); 
  parser.getArgument<string>(gen_op.targetStr_s, gen_op.targetStr_l, gen_op.targetStr); 
  parser.getArgument<string>(gen_op.associationOutput_s, gen_op.associationOutput_l, gen_op.associationOutput);
  parser.getArgument<string>(gen_op.whiteListInput_s, gen_op.whiteListInput_l, gen_op.whiteListInput);
  parser.getArgument<string>(gen_op.blackListInput_s, gen_op.blackListInput_l, gen_op.blackListInput);
  parser.getArgument<string>(gen_op.testInput_s, gen_op.testInput_l, gen_op.testInput);
  parser.getArgument<string>(gen_op.predictionOutput_s, gen_op.predictionOutput_l, gen_op.predictionOutput);
  parser.getArgument<string>(gen_op.logOutput_s,gen_op.logOutput_l,gen_op.logOutput);
  parser.getArgument<string>(gen_op.forestOutput_s,gen_op.forestOutput_l,gen_op.forestOutput);
  parser.getArgument<num_t>(gen_op.pValueThreshold_s, gen_op.pValueThreshold_l, gen_op.pValueThreshold);
  string dataDelimiter,headerDelimiter;
  parser.getArgument<string>(gen_op.dataDelimiter_s, gen_op.dataDelimiter_l, dataDelimiter);
  parser.getArgument<string>(gen_op.headerDelimiter_s, gen_op.headerDelimiter_l, headerDelimiter);
  parser.getArgument<size_t>(gen_op.pruneFeatures_s, gen_op.pruneFeatures_l, gen_op.pruneFeatures);
  parser.getFlag(gen_op.noFilter_s, gen_op.noFilter_l, gen_op.noFilter);
  stringstream ss(dataDelimiter);
  ss >> gen_op.dataDelimiter;
  //cout << gen_op.dataDelimiter;

  ss.clear();
  ss.str("");
  ss << headerDelimiter;
  ss >> gen_op.headerDelimiter;

  // Then read Random Forest specific options
  parser.getArgument<size_t>(RF_op.nTrees_s,RF_op.nTrees_l,RF_op.nTrees);
  parser.getArgument<size_t>(RF_op.mTry_s, RF_op.mTry_l, RF_op.mTry); 
  parser.getArgument<size_t>(RF_op.nMaxLeaves_s, RF_op.nMaxLeaves_l, RF_op.nMaxLeaves);
  parser.getArgument<size_t>(RF_op.nodeSize_s, RF_op.nodeSize_l, RF_op.nodeSize); 
  parser.getArgument<size_t>(RF_op.nPerms_s, RF_op.nPerms_l, RF_op.nPerms); 

  // And last read Gradient Boosting Trees options
  parser.getArgument<size_t>(GBT_op.nTrees_s, GBT_op.nTrees_l, GBT_op.nTrees);
  parser.getArgument<size_t>(GBT_op.nMaxLeaves_s, GBT_op.nMaxLeaves_l, GBT_op.nMaxLeaves);
  parser.getArgument<num_t>(GBT_op.shrinkage_s, GBT_op.shrinkage_l, GBT_op.shrinkage);
  parser.getArgument<num_t>(GBT_op.subSampleSize_s, GBT_op.subSampleSize_l, GBT_op.subSampleSize);

  // See if the help flag was raised
  if(gen_op.printHelp) {
    printHelp(gen_op,RF_op,GBT_op);
    return(EXIT_SUCCESS);
  }
  
  // Extract some handy boolean flags from the options
  bool whiteListInputExists = gen_op.whiteListInput != "";
  bool blackListInputExists = gen_op.blackListInput != "";
  //bool maskInputExists = whiteListInputExists || blackListInputExists;
  bool trainInputExists = gen_op.trainInput != "";
  bool testInputExists = gen_op.testInput != "";
  bool targetExists = gen_op.targetStr != "";
  bool predictionOutputExists = gen_op.predictionOutput != "";
  bool associationOutputExists = gen_op.associationOutput != "";
  bool logOutputExists = gen_op.logOutput != "";
  bool forestOutputExists = gen_op.forestOutput != "";

  // Print help and exit if input file is not specified
  if ( !trainInputExists ) {
    cerr << "Input file not specified" << endl;
    printHelpHint();
    return EXIT_FAILURE;
  }

  // Print help and exit if target index is not specified
  if ( !targetExists ) {
    cerr << "target(s) (-i/--target) not specified" << endl;
    printHelpHint();
    return(EXIT_FAILURE);
  }

  if ( !associationOutputExists && !predictionOutputExists && !forestOutputExists ) {
    cerr << "Cannot do anything!" << endl;
    printHelpHint();
    return(EXIT_FAILURE);
  }

  if ( associationOutputExists && gen_op.noFilter ) {
    cerr << "Cannot generate associations ( -O / --associations ) when filtering is turned OFF ( --noFilter flag raised )" << endl;
    exit(1);
  }

  // Read train data into Treedata object
  cout << "Reading file '" << gen_op.trainInput << "', please wait... " << flush;
  Treedata treedata(gen_op.trainInput,gen_op.dataDelimiter,gen_op.headerDelimiter);
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
  if ( whiteListInputExists ) {
        
    cout << "Reading whitelist '" << gen_op.whiteListInput << "', please wait... " << flush;
    vector<string> whiteFeatureNames = readFeatureMask(gen_op.whiteListInput);
    cout << "DONE" << endl;
    cout << "Applying feature mask, removing " << treedata.nFeatures() - whiteFeatureNames.size() 
	 << " / " << treedata.nFeatures() << " features, please wait... " << flush;
    treedata.whiteList(whiteFeatureNames);
    cout << "DONE" << endl;
  } else if ( blackListInputExists ) {

    cout << "Reading blacklist '" << gen_op.blackListInput << "', please wait... " << flush;
    vector<string> blackFeatureNames = readFeatureMask(gen_op.blackListInput);
    cout << "DONE" << endl;
    cout << "Applying blacklist, keeping " << treedata.nFeatures() - blackFeatureNames.size() 
	 << " / " << treedata.nFeatures() << " features, please wait... " << flush;
    treedata.blackList(blackFeatureNames);
    cout << "DONE" << endl;
  }

  if ( gen_op.pruneFeatures ) {

    cout << "Pruning features with less than " << gen_op.pruneFeatures << " real samples... " << flush;
    size_t nFeaturesOld = treedata.nFeatures();
    pruneFeatures(treedata,gen_op);
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
  cout << "  --" << gen_op.trainInput_l << setw( maxwidth - gen_op.trainInput_l.size() ) << ""
       << "= " << gen_op.trainInput << endl;
  cout << "  --" << gen_op.targetStr_l << setw( maxwidth - gen_op.targetStr_l.size() ) << ""
       << "= " << gen_op.targetStr << " ( index " << targetIdx << " )" << endl;
  cout << "  --" << gen_op.associationOutput_l << setw( maxwidth - gen_op.associationOutput_l.size() ) << ""
       << "= "; if ( associationOutputExists ) { cout << gen_op.associationOutput << endl; } else { cout << "NOT SET" << endl; }
  cout << "  --" << gen_op.testInput_l << setw( maxwidth - gen_op.testInput_l.size() ) << ""
       << "= "; if( testInputExists ) { cout << gen_op.testInput << endl; } else { cout << "NOT SET" << endl; }
  cout << "  --" << gen_op.predictionOutput_l << setw( maxwidth - gen_op.predictionOutput_l.size() ) << ""
       << "= "; if( predictionOutputExists ) { cout << gen_op.predictionOutput << endl; } else { cout << "NOT SET" << endl; }
  cout << "  --" << gen_op.logOutput_l << setw( maxwidth - gen_op.logOutput_l.size() ) << ""
       << "= "; if( logOutputExists ) { cout << gen_op.logOutput << endl; } else { cout << "NOT SET" << endl; }
  cout << "  --" << gen_op.forestOutput_l << setw( maxwidth - gen_op.forestOutput_l.size() ) << ""
       << "= "; if( forestOutputExists ) { cout << gen_op.forestOutput << endl; } else { cout << "NOT SET" << endl; }
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
  cout << "  --pthresold" << setw(8) << "" << "= " << gen_op.pValueThreshold << endl;
  cout << endl;

  bool growGBTPredictor = predictionOutputExists || forestOutputExists ;
  if ( growGBTPredictor ) {
    cout << "Gradient boosting tree configuration for prediction:" << endl;
    cout << "  --" << GBT_op.nTrees_l << setw( maxwidth - GBT_op.nTrees_l.size() ) << ""
         << "= " << GBT_op.nTrees << endl;
    cout << "  --" << GBT_op.nMaxLeaves_l << setw( maxwidth - GBT_op.nMaxLeaves_l.size() ) << ""
         << "= " << GBT_op.nMaxLeaves << endl;
    cout << "  --" << GBT_op.shrinkage_l << setw( maxwidth - GBT_op.shrinkage_l.size() ) << ""
         << "= " << GBT_op.shrinkage << endl;
    cout << "  --" << GBT_op.subSampleSize_l << setw( maxwidth - GBT_op.subSampleSize_l.size() ) << ""
         << "= " << GBT_op.subSampleSize << endl;
    cout << endl;
  }
    
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
  vector<string> featureNames;
  
  if ( !gen_op.noFilter ) {

    cout << "===> Uncovering associations... " << flush;
    RF_stat = executeRandomForest(treedata,RF_op,gen_op,pValues,importanceValues);
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
      
      if ( featureIdx == targetIdx || pValues[featureIdx] <= gen_op.pValueThreshold ) {
	featureNames.push_back(treedata.getFeatureName(featureIdx));
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
  
    if( associationOutputExists ) {

      ofstream toAssociationFile(gen_op.associationOutput.c_str());
      toAssociationFile.precision(8);

      vector<size_t> refIcs( treedata.nFeatures() );
      bool isIncreasingOrder = true;
      datadefs::sortDataAndMakeRef(isIncreasingOrder,pValues,refIcs);
      datadefs::sortFromRef<num_t>(importanceValues,refIcs);
      //datadefs::sortFromRef<string>(featureNames,refIcs);

      assert( gen_op.targetStr == treedata.getFeatureName(targetIdx) );

      for ( size_t i = 0; i < refIcs.size(); ++i ) {
	size_t featureIdx = refIcs[i];

	if ( pValues[i] > gen_op.pValueThreshold ) {
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
  }  
  

  ///////////////////////////////////////////////////////////////////////////
  //  STEP 3 ( OPTIONAL ) -- DATA PREDICTION WITH GRADIENT BOOSTING TREES  //
  ///////////////////////////////////////////////////////////////////////////
  if ( growGBTPredictor ) {      

    cout << "===> Growing GBT predictor... " << flush;
    
    StochasticForest SF(&treedata,gen_op.targetStr,GBT_op.nTrees);
    SF.learnGBT(GBT_op.nMaxLeaves, GBT_op.shrinkage, GBT_op.subSampleSize);
    
    cout << "DONE" << endl;

    if( forestOutputExists ) {
      cout << "===> Writing predictor to file... " << flush;
      SF.printToFile( gen_op.forestOutput );
      cout << "DONE" << endl;
    }

    if ( predictionOutputExists ) {
      
      // TODO: rf_ace.cpp: test input handling is still malfunctioning!
      if ( testInputExists ) {

	cout << "===> Making predictions with test data... " << flush;
	
	Treedata treedata_test(gen_op.testInput,gen_op.dataDelimiter,gen_op.headerDelimiter);
	
	printPredictionToFile(SF,treedata_test,gen_op.targetStr,gen_op.predictionOutput);
	
      } else {
	
	cout << "===> Making predictions with train data... " << flush;

	printPredictionToFile(SF,treedata,gen_op.targetStr,gen_op.predictionOutput);
	
      }
      
      cout << "DONE" << endl;
      
    }
    
  }
  cout << endl;
  
  if ( logOutputExists ) {
    
    ofstream toLogFile(gen_op.logOutput.c_str());
    printHeader(toLogFile);
    RF_stat.print(toLogFile);
    toLogFile.close();

    toLogFile.open("contrasts.tsv");
    RF_stat.printContrastImportance(toLogFile);
    toLogFile.close();
    
  }
  
  cout << 1.0 * ( clock() - clockStart ) / CLOCKS_PER_SEC << " seconds elapsed." << endl << endl;
  
  if ( associationOutputExists ) {
    cout << "Association file '" << gen_op.associationOutput << "' created. Format:" << endl;
    cout << "TARGET   PREDICTOR   LOG10(P-VALUE)   IMPORTANCE   CORRELATION   NSAMPLES" << endl;
    cout << endl;
  }
  
  if ( predictionOutputExists ) {
    cout << "Prediction file '" << gen_op.predictionOutput << "' created. Format:" << endl;
    cout << "TARGET   SAMPLE_ID   DATA      PREDICTION   CONFIDENCE" << endl; 
    cout << endl;
  }
  
  cout << "RF-ACE completed successfully." << endl;
  cout << endl;
      
  return(EXIT_SUCCESS);
}

void pruneFeatures(Treedata& treeData, const General_options& gen_op) {

  vector<string> featureNames;

  size_t targetIdx = treeData.getFeatureIdx(gen_op.targetStr);

  for ( size_t featureIdx = 0; featureIdx < treeData.nFeatures(); ++featureIdx ) {
    
    if ( treeData.nRealSamples(targetIdx,featureIdx) < gen_op.pruneFeatures ) {
      featureNames.push_back(treeData.getFeatureName(featureIdx));
    }

  }
  
  treeData.blackList(featureNames);

}

RF_statistics executeRandomForest(Treedata& treedata,
				  const RF_options& RF_op,
				  const General_options& gen_op,
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



void printPredictionToFile(StochasticForest& SF, Treedata& treeData, const string& targetName, const string& fileName) {
  
  ofstream toPredictionFile(fileName.c_str());

  size_t targetIdx = treeData.getFeatureIdx(targetName);

  if ( treeData.isFeatureNumerical(targetIdx)) {
    
    vector<num_t> prediction;
    vector<num_t> confidence;
    SF.predict(&treeData,prediction,confidence);
    
    for(size_t i = 0; i < prediction.size(); ++i) {
      toPredictionFile << targetName << "\t" << treeData.getSampleName(i) << "\t" << treeData.getRawFeatureData(targetIdx,i)
		       << "\t" << prediction[i] << "\t" << setprecision(3) << confidence[i] << endl;
    }
    
  } else {

    vector<string> prediction;
    vector<num_t> confidence;
    SF.predict(&treeData,prediction,confidence);
    
    for(size_t i = 0; i < prediction.size(); ++i) {
      toPredictionFile << targetName << "\t" << treeData.getSampleName(i) << "\t" << treeData.getRawFeatureData(targetIdx,i)
                       << "\t" << prediction[i] << "\t" << setprecision(3) << confidence[i] << endl;
    }
    
  }

  toPredictionFile.close();
   
}

vector<string> readFeatureMask(const string& fileName) {

  ifstream featurestream;
  featurestream.open(fileName.c_str());
  assert(featurestream.good());
  
  string newFeature;

  set<string> featureMaskSet;
  vector<string> featureMaskVec;
  
  while ( getline(featurestream,newFeature) ) {
    if ( featureMaskSet.find(newFeature) == featureMaskSet.end() ) {
      featureMaskSet.insert(newFeature);
      featureMaskVec.push_back(newFeature);
    }
  }

  return( featureMaskVec );
}

