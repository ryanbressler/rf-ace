#ifndef OPTIONS_HPP
#define OPTIONS_HPP

#include <cstdlib>
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <thread>
#include "argparse.hpp"
#include "datadefs.hpp"
#include "treedata.hpp"
#include "utils.hpp"
#include "distributions.hpp"

using namespace std;
using datadefs::num_t;

namespace options {

  enum ForestType {GBT,RF,CART,UNKNOW};

  // Default general configuration
  const bool       GENERAL_DEFAULT_PRINT_HELP = false;
  const bool       GENERAL_DEFAULT_IS_FILTER = false;
  const size_t     GENERAL_DEFAULT_RECOMBINE_PERMS = 0;
  const bool       GENERAL_DEFAULT_REPORT_CONTRASTS = false;
  const char       GENERAL_DEFAULT_DATA_DELIMITER = '\t';
  const char       GENERAL_DEFAULT_HEADER_DELIMITER = ':';
  const size_t     GENERAL_DEFAULT_MIN_SAMPLES = 5;
  const int        GENERAL_DEFAULT_SEED = -1;
  const ForestType GENERAL_DEFAULT_MODEL_TYPE = RF;
  const size_t     GENERAL_DEFAULT_N_THREADS = 1;
  const bool       GENERAL_DEFAULT_IS_MAX_THREADS = false;

  // Statistical test default configuration
  const size_t     ST_DEFAULT_N_PERMS = 20;
  const num_t      ST_DEFAULT_P_VALUE_THRESHOLD = 0.05;
  const bool       ST_DEFAULT_IS_ADJUSTED_P_VALUE = false;
  const num_t      ST_DEFAULT_IMPORTANCE_THRESHOLD = 0.00001;
  const bool       ST_NORMALIZE_IMPORTANCE_VALUES = false;
  const bool       ST_DEFAULT_REPORT_NONEXISTENT_FEATURES = false;

  // Random Forest default configuration
  const size_t     RF_DEFAULT_N_TREES = 100;
  const size_t     RF_DEFAULT_M_TRY = 0;
  const size_t     RF_DEFAULT_N_MAX_LEAVES = 0;
  const size_t     RF_DEFAULT_NODE_SIZE = 3;
  const num_t      RF_DEFAULT_IN_BOX_FRACTION = 1.0;
  const num_t      RF_DEFAULT_SAMPLE_WITH_REPLACEMENT = true;
  const bool       RF_DEFAULT_USE_CONTRASTS = true;
  const bool       RF_DEFAULT_IS_RANDOM_SPLIT = true;
  const num_t      RF_DEFAULT_SHRINKAGE = 0.0;

  // Gradient Boosting Trees default configuration
  const size_t     GBT_DEFAULT_N_TREES = 100;
  const size_t     GBT_DEFAULT_M_TRY = 0;
  const size_t     GBT_DEFAULT_N_MAX_LEAVES = 6;
  const size_t     GBT_DEFAULT_NODE_SIZE = 3;
  const num_t      GBT_DEFAULT_IN_BOX_FRACTION = 0.5;
  const num_t      GBT_DEFAULT_SAMPLE_WITH_REPLACEMENT = false;
  const bool       GBT_DEFAULT_USE_CONTRASTS = false;
  const bool       GBT_DEFAULT_IS_RANDOM_SPLIT = false;
  const num_t      GBT_DEFAULT_SHRINKAGE = 0.1;

  // CART default configuration
  const size_t     CART_DEFAULT_N_TREES = 1;
  const size_t     CART_DEFAULT_M_TRY = 0;
  const size_t     CART_DEFAULT_N_MAX_LEAVES = 0;
  const size_t     CART_DEFAULT_NODE_SIZE = 3;
  const num_t      CART_DEFAULT_IN_BOX_FRACTION = 1.0;
  const num_t      CART_DEFAULT_SAMPLE_WITH_REPLACEMENT = false;
  const bool       CART_DEFAULT_USE_CONTRASTS = false;
  const bool       CART_DEFAULT_IS_RANDOM_SPLIT = false;
  const num_t      CART_DEFAULT_SHRINKAGE = 0;

  // Determines the amount of indentation in help print-outs
  const size_t maxWidth = 17;
  
  class Base_options {
    
    //private:
    //ArgParse* parser_;
    
  public:
    // EXPERIMENTAL
    size_t recombinePerms; const string recombinePerms_s; const string recombinePerms_l;

    // I/O related and general paramters
    bool printHelp; const string printHelp_s; const string printHelp_l;
    bool isFilter; const string isFilter_s; const string isFilter_l;
    bool reportContrasts; const string reportContrasts_s; const string reportContrasts_l;
    string input; const string input_s; const string input_l;
    string output; const string output_s; const string output_l;
    string forestInput; const string forestInput_s; const string forestInput_l;
    string targetStr; const string targetStr_s; const string targetStr_l;
    string whiteList; const string whiteList_s; const string whiteList_l;
    string blackList; const string blackList_s; const string blackList_l;
    string predictionData; const string predictionData_s; const string predictionData_l;
    string log; const string log_s; const string log_l;
    char dataDelimiter; const string dataDelimiter_s; const string dataDelimiter_l;
    char headerDelimiter; const string headerDelimiter_s; const string headerDelimiter_l;
    size_t pruneFeatures; const string pruneFeatures_s; const string pruneFeatures_l;
    int seed; string seed_s; string seed_l;
    size_t nThreads; const string nThreads_s; const string nThreads_l;
    bool isMaxThreads; const string isMaxThreads_s; const string isMaxThreads_l;

    // Forest Type
    ForestType modelType; const string modelType_s; const string modelType_l;

    // Random Forest related parameters
    size_t nTrees; const string nTrees_s; const string nTrees_l;
    size_t mTry; const string mTry_s; const string mTry_l;
    size_t nMaxLeaves; const string nMaxLeaves_s; const string nMaxLeaves_l;
    size_t nodeSize; const string nodeSize_s; const string nodeSize_l;
    num_t shrinkage; const string shrinkage_s; const string shrinkage_l;

    // Parameters that cannot be controlled through the interface
    num_t inBoxFraction;
    bool sampleWithReplacement;
    bool isRandomSplit;
    bool useContrasts;
    vector<distributions::RandInt> randIntGens;

    // Statistical test related parameters
    size_t nPerms; const string nPerms_s; const string nPerms_l;
    num_t pValueThreshold; const string pValueThreshold_s; const string pValueThreshold_l;
    bool isAdjustedPValue; const string isAdjustedPValue_s; const string isAdjustedPValue_l;
    bool normalizeImportanceValues; const string normalizeImportanceValues_s; const string normalizeImportanceValues_l;
    num_t importanceThreshold; const string importanceThreshold_s; const string importanceThreshold_l;
    bool reportAllFeatures; const string reportAllFeatures_s; const string reportAllFeatures_l;

    Base_options(): 

      // EXPERIMENTAL PARAMETERS
      recombinePerms(GENERAL_DEFAULT_RECOMBINE_PERMS),recombinePerms_s("R"),recombinePerms_l("recombine"),
      // I/O and general parameters
      printHelp(GENERAL_DEFAULT_PRINT_HELP),printHelp_s("h"),printHelp_l("help"),
      isFilter(GENERAL_DEFAULT_IS_FILTER),isFilter_s("f"),isFilter_l("filter"),
      reportContrasts(GENERAL_DEFAULT_REPORT_CONTRASTS),reportContrasts_s("C"),reportContrasts_l("listAllContrasts"),
      input(""),input_s("I"),input_l("input"),
      output(""),output_s("O"),output_l("output"),
      forestInput(""),forestInput_s("F"),forestInput_l("forestInput"),
      targetStr(""),targetStr_s("i"),targetStr_l("target"),
      whiteList(""),whiteList_s("W"),whiteList_l("whiteList"),
      blackList(""),blackList_s("B"),blackList_l("blackList"),
      predictionData(""),predictionData_s("T"),predictionData_l("test"),
      log(""),log_s("L"),log_l("log"),
      dataDelimiter(GENERAL_DEFAULT_DATA_DELIMITER),dataDelimiter_s("D"),dataDelimiter_l("dataDelim"),
      headerDelimiter(GENERAL_DEFAULT_HEADER_DELIMITER),headerDelimiter_s("H"),headerDelimiter_l("headDelim"),
      pruneFeatures(GENERAL_DEFAULT_MIN_SAMPLES),pruneFeatures_s("X"),pruneFeatures_l("pruneFeatures"),
      seed(GENERAL_DEFAULT_SEED),seed_s("S"),seed_l("seed"),
      nThreads(GENERAL_DEFAULT_N_THREADS),nThreads_s("e"),nThreads_l("nThreads"),
      isMaxThreads(GENERAL_DEFAULT_IS_MAX_THREADS),isMaxThreads_s("R"),isMaxThreads_l("maxThreads"),
      // Forest Type
      modelType(GENERAL_DEFAULT_MODEL_TYPE),modelType_s("M"),modelType_l("modelType"),
      // Random Forest related parameters
      // NOTE: Defaults will be loaded inside the constructor
      nTrees_s("n"),nTrees_l("nTrees"),
      mTry_s("m"),mTry_l("mTry"),
      nMaxLeaves_s("a"),nMaxLeaves_l("nMaxLeaves"),
      nodeSize_s("s"),nodeSize_l("nodeSize"),
      shrinkage_s("k"),shrinkage_l("shrinkage"),
      // Statistical test related parameters
      nPerms(ST_DEFAULT_N_PERMS),nPerms_s("p"),nPerms_l("nPerms"),
      pValueThreshold(ST_DEFAULT_P_VALUE_THRESHOLD),pValueThreshold_s("t"),pValueThreshold_l("pValueTh"),
      isAdjustedPValue(ST_DEFAULT_IS_ADJUSTED_P_VALUE),isAdjustedPValue_s("d"),isAdjustedPValue_l("adjustP"),
      normalizeImportanceValues(ST_NORMALIZE_IMPORTANCE_VALUES),normalizeImportanceValues_s("r"),normalizeImportanceValues_l("normImportance"),
      importanceThreshold(ST_DEFAULT_IMPORTANCE_THRESHOLD),importanceThreshold_s("o"),importanceThreshold_l("importanceTh"),
      reportAllFeatures(ST_DEFAULT_REPORT_NONEXISTENT_FEATURES),reportAllFeatures_s("A"),reportAllFeatures_l("listAllFeatures")

    { /* EMPTY CONSTRUCTOR */}
  };

  class General_options : public Base_options {
    
  private: 
    ArgParse parser_;

  public:

    General_options():
      Base_options() { 

      if ( modelType == RF ) {
        setRFDefaults();
      } else if ( modelType == GBT ) {
        setGBTDefaults();
      } else if ( modelType == CART ) {
        setCARTDefaults();
      }

      this->loadUserParams(); 
    }
    
    General_options(const int argc, char* const argv[]):
      Base_options(),
      parser_(argc,argv) { 
      
      if ( modelType == RF ) {
	setRFDefaults(); 
      } else if ( modelType == GBT ) {
	setGBTDefaults();
      } else if ( modelType == CART ) {
	setCARTDefaults();
      } 

      this->loadUserParams();
    }
    
    ~General_options() {
      //if ( parser_ ) {
      //delete parser_;
      //}
    }

    void loadUserParams() {

      // If forest type is explicitly specified, update it
      if ( this->isSet(modelType_s,modelType_l) ) {
        string modelTypeAsStr("");
        parser_.getArgument<string>(modelType_s, modelType_l, modelTypeAsStr);
        if ( modelTypeAsStr == "RF" ) {
          modelType = RF;
        } else if ( modelTypeAsStr == "GBT" ) {
          modelType = GBT;
        } else if ( modelTypeAsStr == "CART" ) {
          modelType = CART;
        } else {
          cerr << "Invalid Forest Type!" << endl;
          exit(1);
        }
      }

      if ( modelType == RF ) {
        setRFDefaults();
      } else if ( modelType == GBT ) {
        setGBTDefaults();
      } else if ( modelType == CART ) {
	setCARTDefaults();
      } else {
	cerr << "Unknown forest type!" << endl;
	exit(1);
      }
            
      // EXPERIMENTAL PARAMETERS
      parser_.getArgument<size_t>(recombinePerms_s,recombinePerms_l,recombinePerms);

      // I/O related and general parameters
      parser_.getFlag(printHelp_s, printHelp_l, printHelp);
      parser_.getFlag(isFilter_s,isFilter_l,isFilter);
      parser_.getArgument<size_t>(recombinePerms_s,recombinePerms_l,recombinePerms);
      parser_.getFlag(reportContrasts_s,reportContrasts_l,reportContrasts);
      parser_.getArgument<string>(input_s, input_l, input);
      parser_.getArgument<string>(targetStr_s, targetStr_l, targetStr);
      parser_.getArgument<string>(output_s, output_l, output);
      parser_.getArgument<string>(forestInput_s, forestInput_l, forestInput);
      parser_.getArgument<string>(whiteList_s, whiteList_l, whiteList);
      parser_.getArgument<string>(blackList_s, blackList_l, blackList);
      parser_.getArgument<string>(predictionData_s, predictionData_l, predictionData);
      parser_.getArgument<string>(log_s,log_l,log);
      string dataDelimiter,headerDelimiter;
      parser_.getArgument<string>(dataDelimiter_s, dataDelimiter_l, dataDelimiter);
      parser_.getArgument<string>(headerDelimiter_s, headerDelimiter_l, headerDelimiter);
      parser_.getArgument<size_t>(pruneFeatures_s, pruneFeatures_l, pruneFeatures);
      stringstream ss(dataDelimiter);
      ss >> dataDelimiter;

      ss.clear();
      ss.str("");
      ss << headerDelimiter;
      ss >> headerDelimiter;

      parser_.getArgument<int>(seed_s, seed_l, seed);
      
      // If no seed was provided, generate one
      if ( seed == GENERAL_DEFAULT_SEED ) {
	seed = distributions::generateSeed();
      }

      parser_.getArgument<size_t>(nThreads_s, nThreads_l, nThreads);
      parser_.getFlag(isMaxThreads_s, isMaxThreads_l, isMaxThreads);

      if ( isMaxThreads ) {
	if ( datadefs::MAX_THREADS == 0 ) {
	  cout << "Cannot extract information about maximum available threads. Using " << nThreads << " threads instead." << endl;
	} else {
	  nThreads = datadefs::MAX_THREADS;
	}
      }
      
      randIntGens.resize(nThreads);
      for ( size_t threadIdx = 0; threadIdx < nThreads; ++threadIdx ) {
	randIntGens[threadIdx].seed(seed + threadIdx);
      }
    
      // Random Forest related parameters
      parser_.getArgument<size_t>( nTrees_s, nTrees_l, nTrees );
      parser_.getArgument<size_t>( mTry_s, mTry_l, mTry );
      parser_.getArgument<size_t>( nMaxLeaves_s, nMaxLeaves_l, nMaxLeaves );
      parser_.getArgument<size_t>( nodeSize_s, nodeSize_l, nodeSize );
      parser_.getArgument<num_t>(  shrinkage_s, shrinkage_l, shrinkage );

      // Statistical test parameters
      parser_.getArgument<size_t>(nPerms_s, nPerms_l, nPerms);
      parser_.getArgument<num_t>(pValueThreshold_s,  pValueThreshold_l,  pValueThreshold);
      parser_.getFlag(isAdjustedPValue_s, isAdjustedPValue_l, isAdjustedPValue);
      parser_.getFlag(normalizeImportanceValues_s, normalizeImportanceValues_l, normalizeImportanceValues);
      parser_.getArgument<num_t>(importanceThreshold_s, importanceThreshold_l, importanceThreshold);
      parser_.getFlag(reportAllFeatures_s, reportAllFeatures_l, reportAllFeatures);

    }

    void validateParameters() {
      
      if ( nPerms > 1 && nPerms < 6 ) {
        cerr << "Use more than 5 permutations in statistical test!" << endl;
        exit(1);
      }
      
      if ( nThreads == 0 ) {
	cerr << "The number of threads must be greater than zero!" << endl;
	exit(1);
      }

      if ( pValueThreshold < 0.0 || pValueThreshold > 1.0 ) {
        cerr << "P-value threshold in statistical test must be within 0...1" << endl;
        exit(1);
      }

      // Print help and exit if input file is not specified
      /*
	if ( input == "" ) {
	cerr << "Input file not specified" << endl;
	this->helpHint();
	exit(1);
	}
	
	// Print help and exit if target index is not specified
	if ( !this->isSet(recombinePerms_s,recombinePerms_l) && targetStr == "" ) {
	cerr << "target not specified" << endl;
	this->helpHint();
	exit(1);
	}
	
	if ( output == "" ) {
	cerr << "You forgot to specify an output file!" << endl;
	this->helpHint();
	exit(1);
	}
      */

      assert( nTrees > 0 );

      if ( isRandomSplit && mTry == 0 ) {
        cerr << "if RF-ACE splits with random features, mTry must be greater than 0" << endl;
	exit(1);
      }

      assert( nMaxLeaves > 0 );
      assert( nodeSize > 0 );

      if ( sampleWithReplacement ) {
	assert( 0.0 < inBoxFraction );
      } else {
	assert( 0.0 < inBoxFraction && inBoxFraction <= 1.0 );
      }

    }

    string printOpt(const string& shortOpt, const string& longOpt) {
      stringstream out;
      out << "  -" << shortOpt << " / --" << longOpt << setw( options::maxWidth - longOpt.size() ) << "";
      return( out.str() );
    }
    
    void printParameters() {

      cout << endl;
      cout << "General configuration:" << endl;
      cout << printOpt(dataDelimiter_s,dataDelimiter_l) << "= '" << dataDelimiter << "' [AFM only]" << endl;
      cout << printOpt(headerDelimiter_s,headerDelimiter_l) << "= '" << headerDelimiter << "' [AFM only]" << endl;
      cout << printOpt(input_s,input_l) << "= " << input << endl;
      cout << printOpt(targetStr_s,targetStr_l) << "= " << targetStr << endl;
      cout << printOpt(output_s,output_l) << "= "; if ( output != "" ) { cout << output << endl; } else { cout << "NOT SET" << endl; }
      cout << printOpt(log_s,log_l) << "= "; if( log != "" ) { cout << log << endl; } else { cout << "NOT SET" << endl; }
      cout << printOpt(seed_s,seed_l) << "= " << seed << endl;
      cout << printOpt(nThreads_s,nThreads_l) << "= " << nThreads << endl;
      cout << endl;

      if(modelType == options::RF ) { 
	cout << "Random Forest configuration:" << endl;
	cout << printOpt(nTrees_s,nTrees_l) << "= " << nTrees << endl;
	cout << printOpt(mTry_s,mTry_l) << "= " << mTry << endl;
      } else if(modelType == options::GBT ) { 
	cout << "Gradient Boosting Trees configuration:" << endl;
	cout << printOpt(nTrees_s,nTrees_l) << "= " << nTrees << endl;
	cout << printOpt(shrinkage_s,shrinkage_l) << "= " << shrinkage << endl;
      } else if(modelType == options::CART ) { 
	cout << "CART configuration:" << endl;
      }      
      cout << printOpt(nMaxLeaves_s,nMaxLeaves_l) << "= " << nMaxLeaves << endl;
      cout << printOpt(nodeSize_s,nodeSize_l) << "= " << nodeSize << endl; 
      cout << endl;

      if ( isFilter ) {
	cout << "Filter configuration:" << endl;
	cout << printOpt(nPerms_s,nPerms_l) << "= " << nPerms << endl;
	cout << printOpt(pValueThreshold_s,pValueThreshold_l) << "= " << pValueThreshold << " (upper limit)" <<endl;
	cout << printOpt(normalizeImportanceValues_s,normalizeImportanceValues_l) << "= ";
	if ( normalizeImportanceValues ) { cout << "YES" << endl; } else { cout << "NO" << endl; }
	cout << printOpt(importanceThreshold_s,importanceThreshold_l) << "= " << importanceThreshold << " (lower limit)" << endl;
	cout << endl;
      }
      
    }

    bool isSet(const string& shortOpt, const string& longOpt) {

      bool ret = false;

      parser_.getFlag(shortOpt,longOpt,ret);
      return( ret );

    }

    template <typename T>
    void setIfNotSet(const string& shortOpt, const string& longOpt, T& oldVal, const T& newVal) {
      if ( !isSet(shortOpt,longOpt) ) {
	oldVal = newVal;
      }
    }

    void helpHint() {
      cout << endl;
      cout << "To get started, type -" << printHelp_s << " / --" << printHelp_l << endl;
    }


    void help() {

      cout << "GENERAL ARGUMENTS:" << endl;
      cout << " -" << input_s << " / --" << input_l << setw( maxWidth - input_l.size() )
           << " " << "Input data file (.afm or .arff)" << endl;
      cout << " -" << targetStr_s << " / --" << targetStr_l << setw( maxWidth - targetStr_l.size() )
           << " " << "Target, specified as integer or string that is to be matched with the content of input" << endl;
      cout << " -" << output_s << " / --" << output_l << setw( maxWidth - output_l.size() )
           << " " << "Output file" << endl;
      cout << " -" << forestInput_s << " / --" << forestInput_l << setw( maxWidth - forestInput_l.size() )
	   << " " << "[Prediction only] Forest input file (.sf), for making predictions" << endl;
      cout << " -" << predictionData_s << " / --" << predictionData_l << setw( maxWidth - predictionData_l.size() )
	   << " " << "[Prediction only] Test data file (.afm or .arff) for prediction" << endl;
      cout << " -" << log_s << " / --" << log_l << setw( maxWidth - log_l.size() )
           << " " << "Log output file" << endl;
      cout << " -" << whiteList_s << " / --" << whiteList_l << setw( maxWidth - whiteList_l.size() )
           << " " << "White list of features to be included in the input file(s)." << endl;
      cout << " -" << blackList_s << " / --" << blackList_l << setw( maxWidth - blackList_l.size() )
           << " " << "Black list of features to be excluded from the input file(s)." << endl;
      cout << " -" << dataDelimiter_s << " / --" << dataDelimiter_l << setw( maxWidth - dataDelimiter_l.size() )
           << " " << "[AFM only] Data delimiter (default \\t)" << endl;
      cout << " -" << headerDelimiter_s << " / --" << headerDelimiter_l << setw( maxWidth - headerDelimiter_l.size() )
           << " " << "[AFM only] Header delimiter that separates the N and C symbols in feature headers from the rest (default " << GENERAL_DEFAULT_HEADER_DELIMITER << ")" << endl;
      cout << " -" << pruneFeatures_s << " / --" << pruneFeatures_l << setw( maxWidth - pruneFeatures_l.size() )
           << " " << "Features with less than n ( default " << GENERAL_DEFAULT_MIN_SAMPLES << " ) samples will be removed" << endl;
      cout << " -" << seed_s << " / --" << seed_l << setw( maxWidth - seed_l.size() ) 
	   << " " << "Seed (positive integer) for the Mersenne Twister random number generator" << endl;
      cout << " -" << modelType_s << " / --" << modelType_l << setw( maxWidth - modelType_l.size() )
	   << " " << "Model Type: RF (default), GBT, or CART" << endl;
      cout << endl;

      cout << "STOCHASTIC FOREST ARGUMENTS:" << endl;
      cout << " -" << nTrees_s << " / --" << nTrees_l << setw( maxWidth - nTrees_l.size() )
	   << " " << "[RF and GBT only] Number of trees in the forest" << endl;
      cout << " -" << mTry_s << " / --" << mTry_l << setw( maxWidth - mTry_l.size() )
	   << " " << "[RF only] Fraction of randomly drawn features per node split" << endl;
      cout << " -" << nMaxLeaves_s << " / --" << nMaxLeaves_l << setw( maxWidth - nMaxLeaves_l.size() )
	   << " " << "Maximum number of leaves per tree" << endl;
      cout << " -" << nodeSize_s << " / --" << nodeSize_l << setw( maxWidth - nodeSize_l.size() )
	   << " " << "Smallest number of train samples per leaf node" << endl;
      cout << " -" << shrinkage_s << " / --" << shrinkage_l << setw( maxWidth - shrinkage_l.size() )
	   << " " << "[GBT only] Shrinkage applied to evolving the residual" << endl;
      cout << endl;

      cout << "FILTER ARGUMENTS:" << endl;
      cout << " -" << isFilter_s << " / --" << isFilter_l << setw( maxWidth - isFilter_l.size() )
           << " " << "Set this flag if you want to perform feature selection (default OFF)" << endl;
      cout << " -" << nPerms_s << " / --" << nPerms_l << setw( maxWidth - nPerms_l.size() )
           << " " << "[Filter only] Number of permutations in statistical test (default " << ST_DEFAULT_N_PERMS << ")" << endl;
      cout << " -" << pValueThreshold_s << " / --" << pValueThreshold_l << setw( maxWidth - pValueThreshold_l.size() )
           << " " << "[Filter only] P-value threshold in statistical test (default " << ST_DEFAULT_P_VALUE_THRESHOLD << ")" << endl;
      cout << " -" << normalizeImportanceValues_s << " / --" << normalizeImportanceValues_l << setw( maxWidth - normalizeImportanceValues_l.size() )
           << " " << "[Filter only] set this flag if wou want importance values normalized (default OFF)" << endl;
      cout << " -" << importanceThreshold_s << " / --" << importanceThreshold_l << setw( maxWidth - importanceThreshold_l.size() )
	   << " " << "[Filter only] Importance threshold in statistical test (default " << ST_DEFAULT_IMPORTANCE_THRESHOLD << ")" << endl;
      cout << " -" << isAdjustedPValue_s << " / --" << isAdjustedPValue_l << setw( maxWidth - isAdjustedPValue_l.size() ) 
	   << " " << "[Filter only] set this flag if wou want Benjamini-Hochberg multiple testing correction (default OFF)" << endl;
      cout << " -" << reportAllFeatures_s << " / --" << reportAllFeatures_l << setw( maxWidth - reportAllFeatures_l.size() )
           << " " << "[Filter only] Set this flag if you want the association file to list all features regardless of statistical significance (default OFF)" << endl; 
      cout << " -" << reportContrasts_s << " / --" << reportContrasts_l << setw( maxWidth - reportContrasts_l.size() )
	   << " " << "[Filter only] Set this flag if you want to report contrasts in the output association file (default OFF)" << endl;
      cout << endl;

      cout << "EXAMPLES:" << endl << endl;
      
      cout << "Performing feature selection using 'target' as the target variable:" << endl
	   << "bin/rf-ace --filter -I data.arff -i target -O associations.tsv" << endl << endl;

      cout << "Performing feature selection with 50 permutations and p-value threshold of 0.001:" << endl
	   << "bin/rf-ace --filter -I data.arff -i 5 -p 50 -t 0.001 -O associations.tsv" << endl << endl;
      
      cout << "Building Random Forest with 1000 trees and mTry of 10 and predicting with test data:" << endl
	   << "bin/rf-ace -I data.arff -i target -T testdata.arff -n 1000 -m 10 -O predictions.tsv" << endl << endl;

      cout << "Building Random Forest predictor and saving it to a file for later use:" << endl
	   << "bin/rf-ace -I data.arff -i target -O rf_predictor.sf" << endl << endl;

      cout << "Loading Random Forest predictor from file and predicting with test data:" << endl
	   << "bin/rf-ace -F rf_predictor.sf -T testdata.arff -O predictions.tsv" << endl << endl;
      
    }
    
    void setRFDefaults() {

      modelType             = RF;
      inBoxFraction         = RF_DEFAULT_IN_BOX_FRACTION;
      sampleWithReplacement = RF_DEFAULT_SAMPLE_WITH_REPLACEMENT;
      isRandomSplit         = RF_DEFAULT_IS_RANDOM_SPLIT;
      useContrasts          = RF_DEFAULT_USE_CONTRASTS;
      nTrees                = RF_DEFAULT_N_TREES;
      mTry                  = RF_DEFAULT_M_TRY;
      nMaxLeaves            = RF_DEFAULT_N_MAX_LEAVES;
      nodeSize              = RF_DEFAULT_NODE_SIZE;
      shrinkage             = RF_DEFAULT_SHRINKAGE;
    }

    void setGBTDefaults() {

      modelType             = GBT;
      inBoxFraction         = GBT_DEFAULT_IN_BOX_FRACTION;
      sampleWithReplacement = GBT_DEFAULT_SAMPLE_WITH_REPLACEMENT;
      isRandomSplit         = GBT_DEFAULT_IS_RANDOM_SPLIT;
      useContrasts          = GBT_DEFAULT_USE_CONTRASTS;
      nTrees                = GBT_DEFAULT_N_TREES;
      mTry                  = GBT_DEFAULT_M_TRY;
      nMaxLeaves            = GBT_DEFAULT_N_MAX_LEAVES;
      nodeSize              = GBT_DEFAULT_NODE_SIZE;
      shrinkage             = GBT_DEFAULT_SHRINKAGE;

    }

    void setCARTDefaults() {

      modelType             = CART;
      inBoxFraction         = CART_DEFAULT_IN_BOX_FRACTION;
      sampleWithReplacement = CART_DEFAULT_SAMPLE_WITH_REPLACEMENT;
      isRandomSplit         = CART_DEFAULT_IS_RANDOM_SPLIT;
      useContrasts          = CART_DEFAULT_USE_CONTRASTS;
      nTrees                = CART_DEFAULT_N_TREES;
      mTry                  = CART_DEFAULT_M_TRY;
      nMaxLeaves            = CART_DEFAULT_N_MAX_LEAVES;
      nodeSize              = CART_DEFAULT_NODE_SIZE;
      shrinkage             = CART_DEFAULT_SHRINKAGE;

    }
    
  };
  
}

#endif
