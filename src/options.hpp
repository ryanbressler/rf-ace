#ifndef OPTIONS_HPP
#define OPTIONS_HPP

#include <cstdlib>
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <algorithm>

#include "utils.hpp"
#include "argparse.hpp"
#include "datadefs.hpp"
#include "distributions.hpp"

using namespace std;
using datadefs::num_t;
using datadefs::forest_t;

class HelpStyler {
public:
  HelpStyler(): maxWidth_(17) {}
  ~HelpStyler() {}
  
  void printHelpLine(const string& shortOpt, const string& longOpt, const string& description) const {
    cout << " -" << shortOpt << " / --" << longOpt << setw( maxWidth_ - longOpt.length() ) << " " << description << endl;
  }

  template<typename T>
  inline void printOption(const string& shortOpt, const string& longOpt, const T& value) const {
    cout << " -" << shortOpt << " / --" << longOpt << setw( maxWidth_ - longOpt.length() ) << " = " << value << endl;
  }
 
  template<typename StartIter, typename StopIter> 
  inline void printOption(const string& shortOpt, const string& longOpt, StartIter startIt, StopIter stopIt) const {
    cout << " -" << shortOpt << " / --" << longOpt << setw( maxWidth_ - longOpt.length() ) << " = " << flush;
    if ( startIt != stopIt ) {
      cout << "[ " << flush;
      utils::write<StartIter,StopIter>(cout,startIt,stopIt,' ');
      cout << " ]" << endl;
    } else {
      cout << "NOT SET" << endl;
    }
  }
 
  void printHelpHint(const string& shortOpt, const string& longOpt) const {
    cout << endl;
    cout << "To get started, type -" << shortOpt << " / --" << longOpt << endl;
  }
  
private:
  size_t maxWidth_;
  
};

template<> inline void HelpStyler::printOption<bool>(const string& shortOpt, const string& longOpt, const bool& value) const {
  cout << " -" << shortOpt << " / --" << longOpt << setw( maxWidth_ - longOpt.length() ) << " = " << flush;
  if ( value ) { cout << "SET" << flush; } else { cout << "NOT SET" << flush; }
  cout << endl;
}

class ForestOptions : public HelpStyler {
public:

  forest_t forestType; const string forestType_s; const string forestType_l;
  size_t nTrees; const string nTrees_s; const string nTrees_l;
  size_t mTry; const string mTry_s; const string mTry_l;
  size_t nMaxLeaves; const string nMaxLeaves_s; const string nMaxLeaves_l;
  size_t nodeSize; const string nodeSize_s; const string nodeSize_l;
  num_t shrinkage; const string shrinkage_s; const string shrinkage_l;
  num_t contrastFraction; const string contrastFraction_s; const string contrastFraction_l;
  bool noNABranching; const string noNABranching_s; const string noNABranching_l;
  vector<num_t> quantiles; const string quantiles_s; const string quantiles_l;
  size_t nSamplesForQuantiles; const string nSamplesForQuantiles_s; const string nSamplesForQuantiles_l;

  num_t inBoxFraction;
  bool sampleWithReplacement;
  bool isRandomSplit;
  bool useContrasts;
  
  ForestOptions(const forest_t ft):
    forestType(ft), forestType_s("f"), forestType_l("forestType"),
    nTrees_s("n"),nTrees_l("nTrees"),
    mTry_s("m"),mTry_l("mTry"),
    nMaxLeaves_s("a"),nMaxLeaves_l("nMaxLeaves"),
    nodeSize_s("s"),nodeSize_l("nodeSize"),
    shrinkage_s("k"),shrinkage_l("shrinkage"),
    contrastFraction_s("c"), contrastFraction_l("contrastFraction"),
    noNABranching(datadefs::SF_DEFAULT_NO_NA_BRANCHING),noNABranching_s("N"), noNABranching_l("noNABranching"),
    quantiles_s("q"), quantiles_l("quantiles"),
    nSamplesForQuantiles_s("r"), nSamplesForQuantiles_l("qSamples") {

    forestType = forest_t::QRF;

    if ( forestType == forest_t::RF ) {
      this->setRFDefaults();
    } else if ( forestType == forest_t::QRF ) {
      this->setQRFDefaults();
    } else if ( forestType == forest_t::GBT ) {
      this->setGBTDefaults();
    } else {
      cerr << "ERROR: wrong forest type!" << endl;
      exit(1);
    }

  }

  ~ForestOptions() {}

  void load(const int argc, char* const argv[]) {
    ArgParse parser(argc,argv);

    bool isSet = false;
    parser.getFlag(forestType_s, forestType_l, isSet);
    if ( isSet ) {
      string forestTypeAsStr = "";
      parser.getArgument<string>(forestType_s, forestType_l, forestTypeAsStr);
      //forestType = datadefs::forestTypeAssign.at(forestTypeAsStr);
      if ( forestType == forest_t::RF ) {
	this->setRFDefaults();
      } else if ( forestType == forest_t::QRF ) {
	this->setQRFDefaults();
      } else if ( forestType == forest_t::GBT ) {
	this->setGBTDefaults();
      } else {
	cerr << "GeneralOptions::load() -- unknown forest type: " << forestTypeAsStr << endl;
	exit(1);
      }
    }
  
    parser.getArgument<size_t>( nTrees_s,           nTrees_l,           nTrees );
    parser.getArgument<size_t>( mTry_s,             mTry_l,             mTry );
    parser.getArgument<size_t>( nMaxLeaves_s,       nMaxLeaves_l,       nMaxLeaves );
    parser.getArgument<size_t>( nodeSize_s,         nodeSize_l,         nodeSize );
    parser.getArgument<num_t>(  shrinkage_s,        shrinkage_l,        shrinkage );
    parser.getArgument<num_t>(  contrastFraction_s, contrastFraction_l, contrastFraction );
    parser.getFlag(             noNABranching_s,    noNABranching_l,    noNABranching );

    string quantilesAsStr;
    parser.getArgument<string>( quantiles_s,        quantiles_l,        quantilesAsStr );
    vector<string> quantileList = utils::split(quantilesAsStr,',');
    quantiles.resize(quantileList.size());
    transform(quantileList.begin(),quantileList.end(),quantiles.begin(),utils::str2<num_t>);

    parser.getArgument<size_t>( nSamplesForQuantiles_s, nSamplesForQuantiles_l, nSamplesForQuantiles );

  }

  void setRFDefaults() {
    forestType            = forest_t::RF;
    inBoxFraction         = datadefs::RF_DEFAULT_IN_BOX_FRACTION;
    sampleWithReplacement = datadefs::RF_DEFAULT_SAMPLE_WITH_REPLACEMENT;
    isRandomSplit         = datadefs::RF_DEFAULT_IS_RANDOM_SPLIT;
    useContrasts          = datadefs::RF_DEFAULT_USE_CONTRASTS;
    contrastFraction      = datadefs::RF_DEFAULT_CONTRAST_FRACTION;
    nTrees                = datadefs::RF_DEFAULT_N_TREES;
    mTry                  = datadefs::RF_DEFAULT_M_TRY;
    nMaxLeaves            = datadefs::RF_DEFAULT_N_MAX_LEAVES;
    nodeSize              = datadefs::RF_DEFAULT_NODE_SIZE;
    shrinkage             = datadefs::RF_DEFAULT_SHRINKAGE;
    quantiles             = datadefs::RF_DEFAULT_QUANTILES;
    nSamplesForQuantiles  = datadefs::RF_DEFAULT_N_SAMPLES_FOR_QUANTILES;

  }

  void setQRFDefaults() {
    forestType            = forest_t::QRF;
    inBoxFraction         = datadefs::QRF_DEFAULT_IN_BOX_FRACTION;
    sampleWithReplacement = datadefs::QRF_DEFAULT_SAMPLE_WITH_REPLACEMENT;
    isRandomSplit         = datadefs::QRF_DEFAULT_IS_RANDOM_SPLIT;
    useContrasts          = datadefs::QRF_DEFAULT_USE_CONTRASTS;
    contrastFraction      = datadefs::QRF_DEFAULT_CONTRAST_FRACTION;
    nTrees                = datadefs::QRF_DEFAULT_N_TREES;
    mTry                  = datadefs::QRF_DEFAULT_M_TRY;
    nMaxLeaves            = datadefs::QRF_DEFAULT_N_MAX_LEAVES;
    nodeSize              = datadefs::QRF_DEFAULT_NODE_SIZE;
    shrinkage             = datadefs::QRF_DEFAULT_SHRINKAGE;
    quantiles             = datadefs::QRF_DEFAULT_QUANTILES;
    nSamplesForQuantiles  = datadefs::QRF_DEFAULT_N_SAMPLES_FOR_QUANTILES;
  }

  void setGBTDefaults() {
    forestType            = forest_t::GBT;
    inBoxFraction         = datadefs::GBT_DEFAULT_IN_BOX_FRACTION;
    sampleWithReplacement = datadefs::GBT_DEFAULT_SAMPLE_WITH_REPLACEMENT;
    isRandomSplit         = datadefs::GBT_DEFAULT_IS_RANDOM_SPLIT;
    useContrasts          = datadefs::GBT_DEFAULT_USE_CONTRASTS;
    contrastFraction      = datadefs::GBT_DEFAULT_CONTRAST_FRACTION;
    nTrees                = datadefs::GBT_DEFAULT_N_TREES;
    mTry                  = datadefs::GBT_DEFAULT_M_TRY;
    nMaxLeaves            = datadefs::GBT_DEFAULT_N_MAX_LEAVES;
    nodeSize              = datadefs::GBT_DEFAULT_NODE_SIZE;
    shrinkage             = datadefs::GBT_DEFAULT_SHRINKAGE;
    quantiles             = datadefs::GBT_DEFAULT_QUANTILES;
    nSamplesForQuantiles  = datadefs::GBT_DEFAULT_N_SAMPLES_FOR_QUANTILES;
  }

  void validate() {
    
    if ( forestType == forest_t::UNKNOWN ) {
      cerr << "ERROR: forest type must be set!" << endl;
      exit(1);
    }
    
    if ( nTrees == 0 ) {
      cerr << "ERROR: the number of trees must be set!" << endl;
      exit(1);
    } 
    
    if ( isRandomSplit && mTry == 0 ) {
      cerr << "ERROR: For random split mTry must be set!" << endl;
      exit(1);
    }
    
    if ( nMaxLeaves == 0 ) {
      cerr << "ERROR: nMaxLeaves must be set!" << endl;
      exit(1);
    }
    
    if ( nodeSize == 0 ) {
      cerr << "ERROR: nodeSize must be set!" << endl;
      exit(1);
    }
    
    if ( useContrasts && ( contrastFraction <= 0.0 || contrastFraction >= 1.0 ) ) {
      cerr << "ERROR: contrastFraction must be between (0,1)" << endl;
      exit(1);
    }
    
    if ( inBoxFraction <= 0.0 || inBoxFraction > 1.0 ) {
      cerr << "ERROR: inBoxFraction must be between (0,1]" << endl;
      exit(1);
    }

    if ( forestType == forest_t::QRF ) {
      for ( size_t i = 0; i < quantiles.size(); ++i ) {
	if ( 0.0 >= quantiles[i] || quantiles[i] > 1.0 ) {
	  cerr << "ERROR: quantiles must be real numbers between (0,1]" << endl;
	  exit(1);
	}
      }
      if (nSamplesForQuantiles == 0) {
	cerr << "ERROR: using 0 samples for quantile predictions" << endl;
	exit(1);
      }
    }

  }

  void help() {
    cout << "Forest Options:" << endl;
    this->printHelpLine(forestType_s,forestType_l,"Forest type: RF (default), QRF, or GBT");
    this->printHelpLine(nTrees_s,nTrees_l,"Number of trees in the forest");
    this->printHelpLine(mTry_s,mTry_l,"[RF+QRF] Fraction of randomly drawn features per node split");
    this->printHelpLine(nMaxLeaves_s,nMaxLeaves_l,"Maximum number of leaves per tree");
    this->printHelpLine(nodeSize_s,nodeSize_l,"Smallest number of train samples per leaf node");
    this->printHelpLine(shrinkage_s,shrinkage_l,"[GBT] Shrinkage applied to evolving the residual");
    this->printHelpLine(contrastFraction_s,contrastFraction_l,"[Filter only] the fraction of contrast features sampled to approximate the null distribution");
    this->printHelpLine(noNABranching_s,noNABranching_l,"If set, the splitter will NOT create a third branch for cases where the splitter is NA");
    this->printHelpLine(quantiles_s,quantiles_l,"[QRF] comma-separated list of quantiles to build a Quantile Random Forest from");
    this->printHelpLine(nSamplesForQuantiles_s,nSamplesForQuantiles_l,"[QRF] specify the number of samples per tree for calculating the quantiles");
  }

  void print() {
    if ( forestType == forest_t::RF ) { 
      cout << "Random Forest (RF) configuration:" << endl;
      this->printOption(nTrees_s,nTrees_l,nTrees);
      this->printOption(mTry_s,mTry_l,mTry);
      this->printOption(nodeSize_s,nodeSize_l,nodeSize);
      this->printOption(nMaxLeaves_s,nMaxLeaves_l,nMaxLeaves);
    } else if ( forestType == forest_t::QRF ) {
      cout << "Quantile Random Forest (QRF) configuration:" << endl;
      this->printOption(nTrees_s,nTrees_l,nTrees);
      this->printOption(mTry_s,mTry_l,mTry);
      this->printOption(nodeSize_s,nodeSize_l,nodeSize);
      this->printOption(nMaxLeaves_s,nMaxLeaves_l,nMaxLeaves);
      this->printOption(quantiles_s,quantiles_l,quantiles.begin(),quantiles.end());
      this->printOption(nSamplesForQuantiles_s,nSamplesForQuantiles_l,nSamplesForQuantiles);
    } else if ( forestType == forest_t::GBT ) { 
      cout << "Gradient Boosting Tree (GBT) configuration:" << endl;
      this->printOption(nTrees_s,nTrees_l,nTrees);
      this->printOption(mTry_s,mTry_l,mTry);
      this->printOption(nodeSize_s,nodeSize_l,nodeSize);
      this->printOption(nMaxLeaves_s,nMaxLeaves_l,nMaxLeaves);
      this->printOption(shrinkage_s,shrinkage_l,shrinkage);
    } else {
      cerr << "ERROR: unknown model type to print parameters for!" << endl;
      exit(1);
    }
    this->printOption(noNABranching_s,noNABranching_l,noNABranching);
    cout << endl;
  }
   
};

class FilterOptions : public HelpStyler {
public:

  // Statistical test related parameters
  size_t nPerms; const string nPerms_s; const string nPerms_l;
  num_t pValueThreshold; const string pValueThreshold_s; const string pValueThreshold_l;
  //bool isAdjustedPValue; const string isAdjustedPValue_s; const string isAdjustedPValue_l;
  //bool normalizeImportanceValues; const string normalizeImportanceValues_s; const string normalizeImportanceValues_l;
  //num_t importanceThreshold; const string importanceThreshold_s; const string importanceThreshold_l;
  //bool reportAllFeatures; const string reportAllFeatures_s; const string reportAllFeatures_l;

  FilterOptions():
    nPerms(datadefs::FILTER_DEFAULT_N_PERMS),nPerms_s("p"),nPerms_l("nPerms"),
    pValueThreshold(datadefs::FILTER_DEFAULT_P_VALUE_THRESHOLD),pValueThreshold_s("t"),pValueThreshold_l("pValueTh") {}
    //isAdjustedPValue(datadefs::FILTER_DEFAULT_IS_ADJUSTED_P_VALUE),isAdjustedPValue_s("d"),isAdjustedPValue_l("adjustP"),
    //normalizeImportanceValues(datadefs::FILTER_NORMALIZE_IMPORTANCE_VALUES),normalizeImportanceValues_s("r"),normalizeImportanceValues_l("normImportance"),
    //importanceThreshold(datadefs::FILTER_DEFAULT_IMPORTANCE_THRESHOLD),importanceThreshold_s("o"),importanceThreshold_l("importanceTh") {}
    //reportAllFeatures(datadefs::FILTER_DEFAULT_REPORT_NONEXISTENT_FEATURES),reportAllFeatures_s("A"),reportAllFeatures_l("listAllFeatures") {}

  ~FilterOptions() {}

  void load(const int argc, char* const argv[]) {
    ArgParse parser(argc,argv);
    parser.getArgument<size_t>(nPerms_s,nPerms_l,nPerms);
    parser.getArgument<num_t>(pValueThreshold_s,pValueThreshold_l,pValueThreshold);
    //parser.getArgument<bool>(isAdjustedPValue_s,isAdjustedPValue_l,isAdjustedPValue);
    //parser.getArgument<num_t>(importanceThreshold_s,importanceThreshold_l,importanceThreshold);
    //parser.getArgument<bool>(normalizeImportanceValues_s,normalizeImportanceValues_l,normalizeImportanceValues);
    //parser.getArgument<bool>(reportAllFeatures_s,reportAllFeatures_l,reportAllFeatures);
  }

  void validate() {
    assert( nPerms >= 5 );
    assert( 0 <= pValueThreshold && pValueThreshold <= 1.0 );
  }

  void help() {
    cout << "Filter Options:" << endl;
    this->printHelpLine(nPerms_s,nPerms_l,"Number of permutations in statistical test");
    this->printHelpLine(pValueThreshold_s,pValueThreshold_l,"P-value threshold in statistical test");
    //this->printHelpLine(isAdjustedPValue_s,isAdjustedPValue_l,"Flag to turn ON Benjamini-Hochberg multiple testing correction");
    //this->printHelpLine(importanceThreshold_s,importanceThreshold_l,"Importance threshold");
    //this->printHelpLine(normalizeImportanceValues_s,normalizeImportanceValues_l,"Flag to turn ON normalization of importance scores");
    //this->printHelpLine(reportAllFeatures_s,reportAllFeatures_l,"Flag to turn ON reporting of all features, regardless of statistical significance");
  }

  void print() {
    cout << "Filter options:" << endl;
    this->printOption(nPerms_s,nPerms_l,nPerms);
    this->printOption(pValueThreshold_s,pValueThreshold_l,pValueThreshold);
    cout << endl;
    //cout << "isAdjustedPValue = " << isAdjustedPValue << endl;
    //cout << "normalizeImportanceValues = " << normalizeImportanceValues << endl;
    //cout << "importanceThreshold = " << importanceThreshold << endl;
    //cout << "reportAllFeatures = " << reportAllFeatures << endl;
  }

};

class IO : public HelpStyler {
public:
  
  string filterDataFile; const string filterDataFile_s; const string filterDataFile_l;
  string trainDataFile; const string trainDataFile_s; const string trainDataFile_l;
  string testDataFile; const string testDataFile_s; const string testDataFile_l;
  string loadForestFile; const string loadForestFile_s; const string loadForestFile_l;
  string saveForestFile; const string saveForestFile_s; const string saveForestFile_l;
  string associationsFile; const string associationsFile_s; const string associationsFile_l;
  string predictionsFile; const string predictionsFile_s; const string predictionsFile_l;
  string pairInteractionsFile; const string pairInteractionsFile_s; const string pairInteractionsFile_l;
  string logFile; const string logFile_s; const string logFile_l;
  string featureWeightsFile; const string featureWeightsFile_s; const string featureWeightsFile_l;
  string whiteListFile; const string whiteListFile_s; const string whiteListFile_l;
  string blackListFile; const string blackListFile_s; const string blackListFile_l;
  
  
  IO():
    filterDataFile_s("F"), filterDataFile_l("filterData"),
    trainDataFile_s("I"), trainDataFile_l("trainData"),
    testDataFile_s("T"), testDataFile_l("testData"),
    loadForestFile_s("L"), loadForestFile_l("loadForest"),
    saveForestFile_s("V"), saveForestFile_l("saveForest"),
    associationsFile_s("A"), associationsFile_l("associations"),
    predictionsFile_s("P"), predictionsFile_l("predictions"),
    pairInteractionsFile_s("R"), pairInteractionsFile_l("pairInteractions"),
    logFile_s("G"), logFile_l("log"),
    featureWeightsFile_s("w"), featureWeightsFile_l("featureWeights"),
    whiteListFile_s("W"), whiteListFile_l("whiteList"),
    blackListFile_s("B"), blackListFile_l("blackList") {}

  ~IO() {}

  void load(const int argc, char* const argv[]) {
    ArgParse parser(argc,argv);
    parser.getArgument<string>(filterDataFile_s,filterDataFile_l,filterDataFile);
    parser.getArgument<string>(trainDataFile_s,trainDataFile_l,trainDataFile);
    parser.getArgument<string>(testDataFile_s,testDataFile_l,testDataFile);
    parser.getArgument<string>(loadForestFile_s,loadForestFile_l,loadForestFile);
    parser.getArgument<string>(saveForestFile_s,saveForestFile_l,saveForestFile);
    parser.getArgument<string>(associationsFile_s,associationsFile_l,associationsFile);
    parser.getArgument<string>(predictionsFile_s,predictionsFile_l,predictionsFile);
    parser.getArgument<string>(pairInteractionsFile_s,pairInteractionsFile_l,pairInteractionsFile);
    parser.getArgument<string>(logFile_s,logFile_l,logFile);
    parser.getArgument<string>(featureWeightsFile_s,featureWeightsFile_l,featureWeightsFile);
    parser.getArgument<string>(whiteListFile_s,whiteListFile_l,whiteListFile);
    parser.getArgument<string>(blackListFile_s,blackListFile_l,blackListFile);
  }

  void help() {
    cout << "File Options:" << endl;
    this->printHelpLine(filterDataFile_s,filterDataFile_l,"Load data file (.afm or .arff) for feature selection");
    this->printHelpLine(trainDataFile_s,trainDataFile_l,"Load data file (.afm or .arff) for training a model");
    this->printHelpLine(featureWeightsFile_s,featureWeightsFile_l,"Load feature weights from file");
    this->printHelpLine(whiteListFile_s,whiteListFile_l,"Load white list from file");
    this->printHelpLine(blackListFile_s,blackListFile_l,"Load black list from file");
    this->printHelpLine(testDataFile_s,testDataFile_l,"Load data file (.afm or .arff) for testing a model");
    this->printHelpLine(loadForestFile_s,loadForestFile_l,"Load model from file (.sf)");
    this->printHelpLine(saveForestFile_s,saveForestFile_l,"Save model to file (.sf)");
    this->printHelpLine(associationsFile_s,associationsFile_l,"Save associations to file");
    this->printHelpLine(predictionsFile_s,predictionsFile_l,"Save predictions to file");
    this->printHelpLine(pairInteractionsFile_s,pairInteractionsFile_l,"Save pair interactions to file");
    this->printHelpLine(logFile_s,logFile_l,"Save log to file");
  }

  void print() {
    cout << "filterDataFile = " << filterDataFile << endl;
    cout << "trainDataFile = " << trainDataFile << endl;
    cout << "testDataFile = " << testDataFile << endl;
    cout << "loadForestFile = " << loadForestFile << endl;
    cout << "saveForestFile = " << saveForestFile << endl;
    cout << "associationsFile = " << associationsFile << endl;
    cout << "predictionsFile = " << predictionsFile << endl;
    cout << "pairInteractionsFile = " << pairInteractionsFile << endl;
    cout << "logFile = " << logFile << endl;
    cout << "featureWeightsFile = " << featureWeightsFile << endl;
    cout << "whiteListFile = " << whiteListFile << endl;
    cout << "blackListFile = " << blackListFile << endl;
  }
  
  void validate() {
    
    size_t iter = 0;

    iter += featureWeightsFile != "" ? 1 : 0 ;
    iter += whiteListFile != "" ? 1 : 0 ;
    iter += blackListFile != "" ? 1 : 0 ;

    if ( iter > 1 ) {
      cerr << "ERROR: Specify only one of the following: feature weights, whitelist, or blacklist" << endl;
      exit(1);
    }

  }

};

class GeneralOptions : public HelpStyler {
public:

  bool printHelp; const string printHelp_s; const string printHelp_l;
  string targetStr; const string targetStr_s; const string targetStr_l;
  char dataDelimiter; const string dataDelimiter_s; const string dataDelimiter_l;
  char headerDelimiter; const string headerDelimiter_s; const string headerDelimiter_l;
  size_t pruneFeatures; const string pruneFeatures_s; const string pruneFeatures_l;
  int seed; string seed_s; string seed_l;
  size_t nThreads; const string nThreads_s; const string nThreads_l;
  bool isMaxThreads; const string isMaxThreads_s; const string isMaxThreads_l;
  num_t defaultFeatureWeight; const string defaultFeatureWeight_s; const string defaultFeatureWeight_l;  

  GeneralOptions():
    printHelp(datadefs::GENERAL_DEFAULT_PRINT_HELP),printHelp_s("h"),printHelp_l("help"),
    targetStr(""),targetStr_s("i"),targetStr_l("target"),
    dataDelimiter(datadefs::GENERAL_DEFAULT_DATA_DELIMITER),dataDelimiter_s("D"),dataDelimiter_l("dataDelim"),
    headerDelimiter(datadefs::GENERAL_DEFAULT_HEADER_DELIMITER),headerDelimiter_s("H"),headerDelimiter_l("headDelim"),
    pruneFeatures(datadefs::GENERAL_DEFAULT_MIN_SAMPLES),pruneFeatures_s("X"),pruneFeatures_l("pruneFeatures"),
    seed(datadefs::GENERAL_DEFAULT_SEED),seed_s("S"),seed_l("seed"),
    nThreads(datadefs::GENERAL_DEFAULT_N_THREADS),nThreads_s("e"),nThreads_l("nThreads"),
    isMaxThreads(datadefs::GENERAL_DEFAULT_IS_MAX_THREADS),isMaxThreads_s("R"),isMaxThreads_l("maxThreads"),
    defaultFeatureWeight(datadefs::GENERAL_DEFAULT_FEATURE_WEIGHT),defaultFeatureWeight_s("d"),defaultFeatureWeight_l("defaultWeight") {}
  ~GeneralOptions() {}

  void load(const int argc, char* const argv[]) {
    ArgParse parser(argc,argv);
    parser.getFlag(printHelp_s, printHelp_l, printHelp);
    parser.getArgument<string>(targetStr_s, targetStr_l, targetStr);
    string dataDelimiterAsStr, headerDelimiterAsStr;
    parser.getArgument<string>(dataDelimiter_s, dataDelimiter_l, dataDelimiterAsStr);
    parser.getArgument<string>(headerDelimiter_s, headerDelimiter_l, headerDelimiterAsStr);
    parser.getArgument<size_t>(pruneFeatures_s, pruneFeatures_l, pruneFeatures);
    parser.getArgument<num_t>(defaultFeatureWeight_s, defaultFeatureWeight_l, defaultFeatureWeight);

    stringstream ss(dataDelimiterAsStr);
    ss >> dataDelimiter;

    ss.clear();
    ss.str("");
    ss << headerDelimiterAsStr;
    ss >> headerDelimiter;

    parser.getArgument<int>(seed_s, seed_l, seed);
    parser.getArgument<size_t>(nThreads_s, nThreads_l, nThreads);
    parser.getFlag(isMaxThreads_s, isMaxThreads_l, isMaxThreads);
  }

  void validate() {
    assert( targetStr != "" );
    //assert( dataDelimiter != '' );
    //assert( headerDelimiter != '' );
    assert( defaultFeatureWeight >= 0.0 );
    assert( nThreads >= 1 );
    //assert( seed >= 0 );
  }

  void help() {
    cout << "General Options:" << endl;
    this->printHelpLine(targetStr_s,targetStr_l,"Name or index for the variable to make modeling for");
    this->printHelpLine(dataDelimiter_s,dataDelimiter_l,"[AFM only] Field delimiter");
    this->printHelpLine(headerDelimiter_s,headerDelimiter_l,"[AFM only] Feature type and name delimiter");
    this->printHelpLine(pruneFeatures_s,pruneFeatures_l,"Prune Features");
    this->printHelpLine(seed_s,seed_l,"Seed for random number generator");
    this->printHelpLine(nThreads_s,nThreads_l,"Number of threads if using multithreading");
    this->printHelpLine(isMaxThreads_s,isMaxThreads_l,"Flag to make use of all available threads");
    this->printHelpLine(defaultFeatureWeight_s,defaultFeatureWeight_l,"Default feature weight, if using feature weighting");
  }

  void print() {
    cout << "printHelp = " << printHelp << endl;
    cout << "targetStr = " << targetStr << endl;
    cout << "dataDelimiter = " << dataDelimiter << endl;
    cout << "headerDelimiter = " << headerDelimiter << endl;
    cout << "pruneFeatures = " << pruneFeatures << endl;
    cout << "seed = " << seed << endl;
    cout << "nThreads = " << nThreads << endl;
    cout << "isMaxThreads = " << isMaxThreads << endl;
    cout << "defaultFeatureWeight = " << defaultFeatureWeight << endl;
  }

};

class Options : public HelpStyler {

public:
  
  ForestOptions forestOptions;
  IO io;
  GeneralOptions generalOptions;
  FilterOptions filterOptions;
  
  Options(forest_t ft): 
    forestOptions(ft) {}

  ~Options() {}
  
  void load(const int argc, char* const argv[]) {
    forestOptions.load(argc,argv);
    io.load(argc,argv);
    generalOptions.load(argc,argv);
    filterOptions.load(argc,argv);
  }

  void help() {
  
    forestOptions.help();
    io.help();
    filterOptions.help();
    generalOptions.help();
    this->printExamples();
  
  }

  void print() {

    forestOptions.print();
    io.print();
    filterOptions.print();
    generalOptions.print();
  }

  void printExamples() {
    
    cout << endl;
    cout << "Feature selection examples:" << endl
	 << "bin/rf-ace -" << io.filterDataFile_s << " data.afm -" << generalOptions.targetStr_s << " target -" << forestOptions.nTrees_s << " 100 -" << forestOptions.mTry_s << " 30 -" << io.associationsFile_s << " associations.tsv" << endl
	 << "bin/rf-ace -" << io.filterDataFile_s << " data.afm -" << generalOptions.targetStr_s << " 5      -" << forestOptions.nTrees_s << " 100 -" << forestOptions.mTry_s << " 30 -" << filterOptions.nPerms_s << " 50 -" << filterOptions.pValueThreshold_s << " 0.001 -" << io.associationsFile_s << " associations.tsv" << endl << endl;
    
    cout << "Model training & prediction examples:" << endl
	 << "bin/rf-ace -" << io.trainDataFile_s << " data.afm -" << generalOptions.targetStr_s << " target -" << io.testDataFile_s << " testdata.afm -" << forestOptions.nTrees_s << " 100 -" << forestOptions.mTry_s << " 30 -" << io.predictionsFile_s << " predictions.tsv" << endl
	 << "bin/rf-ace -" << io.trainDataFile_s << " data.afm -" << generalOptions.targetStr_s << " target -" << io.saveForestFile_s << " model.sf" << endl
	 << "bin/rf-ace -" << io.loadForestFile_s << " model.sf -" << io.testDataFile_s << " testdata.afm -" << io.predictionsFile_s << " predictions.tsv" << endl << endl;
    
  }
  
};

#endif
