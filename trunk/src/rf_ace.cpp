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
//#include <omp.h>

#include "argparse.hpp"
#include "stochasticforest.hpp"
#include "treedata.hpp"
#include "datadefs.hpp"

using namespace std;
using datadefs::num_t;

const bool   GENERAL_DEFAULT_PRINT_HELP = false;
const bool   GENERAL_DEFAULT_NO_FILTER = false;
const bool   GENERAL_DEFAULT_NO_PREDICTION = false; // TEMPORARY VARIABLE
const num_t  GENERAL_DEFAULT_P_VALUE_THRESHOLD = 0.10;

const bool   RF_IS_OPTIMIZED_NODE_SPLIT = false;
const size_t RF_DEFAULT_N_TREES = 100; // zero means it will be estimated from the data by default
const size_t RF_DEFAULT_M_TRY = 0; // same here ...
const size_t RF_DEFAULT_N_MAX_LEAVES = 10;
const size_t RF_DEFAULT_NODE_SIZE = 5; // ... and here
const size_t RF_DEFAULT_N_PERMS = 50;

const bool   GBT_IS_OPTIMIZED_NODE_SPLIT = false;
const size_t GBT_DEFAULT_N_TREES = 100;
const size_t GBT_DEFAULT_N_MAX_LEAVES = 6;
const num_t  GBT_DEFAULT_SHRINKAGE = 0.1;
const num_t  GBT_DEFAULT_SUB_SAMPLE_SIZE = 0.5;

struct General_options {

  bool   printHelp;
  string printHelp_s;
  string printHelp_l;

  string trainInput;
  string trainInput_s;
  string trainInput_l;

  string associationOutput;
  string associationOutput_s;
  string associationOutput_l;

  string targetStr;
  string targetStr_s;
  string targetStr_l;

  string featureMaskInput;
  string featureMaskInput_s;
  string featureMaskInput_l;

  string testInput;
  string testInput_s;
  string testInput_l;

  string predictionOutput;
  string predictionOutput_s;
  string predictionOutput_l;

  bool   noFilter;
  string noFilter_s;
  string noFilter_l;

  num_t  pValueThreshold;
  string pValueThreshold_s;
  string pValueThreshold_l;

  General_options():
    printHelp(GENERAL_DEFAULT_PRINT_HELP),
    printHelp_s("h"),
    printHelp_l("help"),    
    
    trainInput(""),
    trainInput_s("I"),
    trainInput_l("traindata"),
    
    associationOutput(""),
    associationOutput_s("O"),
    associationOutput_l("associations"),
    
    targetStr(""),
    targetStr_s("i"),
    targetStr_l("target"),

    featureMaskInput(""),
    featureMaskInput_s("M"),
    featureMaskInput_l("fmask"),
    
    testInput(""),
    testInput_s("T"),
    testInput_l("testdata"),

    predictionOutput(""),
    predictionOutput_s("P"),
    predictionOutput_l("predictions"),
  
    noFilter(GENERAL_DEFAULT_NO_FILTER),
    noFilter_s("f"),
    noFilter_l("noFilter"),

    pValueThreshold(GENERAL_DEFAULT_P_VALUE_THRESHOLD),
    pValueThreshold_s("t"),
    pValueThreshold_l("pthreshold") {}
};

struct RF_options {
  
  bool isOptimizedNodeSplit;
  string isOptimizedNodeSplit_s;
  string isOptimizedNodeSplit_l;
  
  size_t nTrees;
  string nTrees_s;
  string nTrees_l;
  
  size_t mTry;
  string mTry_s;
  string mTry_l;

  size_t nMaxLeaves;
  string nMaxLeaves_s;
  string nMaxLeaves_l;

  size_t nodeSize;
  string nodeSize_s;
  string nodeSize_l;

  size_t nPerms;
  string nPerms_s;
  string nPerms_l;

  RF_options():
    isOptimizedNodeSplit(RF_IS_OPTIMIZED_NODE_SPLIT),
    isOptimizedNodeSplit_s("o"),
    isOptimizedNodeSplit_l("RF_optimize"),

    nTrees(RF_DEFAULT_N_TREES),
    nTrees_s("n"),
    nTrees_l("RF_ntrees"),
    
    mTry(RF_DEFAULT_M_TRY),
    mTry_s("m"),
    mTry_l("RF_mtry"),

    nMaxLeaves(RF_DEFAULT_N_MAX_LEAVES),
    nMaxLeaves_s("a"),
    nMaxLeaves_l("RF_maxleaves"),

    nodeSize(RF_DEFAULT_NODE_SIZE),
    nodeSize_s("s"),
    nodeSize_l("RF_nodesize"),

    nPerms(RF_DEFAULT_N_PERMS),
    nPerms_s("p"),
    nPerms_l("RF_nperms") {}
};

struct GBT_options {

  bool isOptimizedNodeSplit;
  string isOptimizedNodeSplit_s;
  string isOptimizedNodeSplit_l;

  size_t nTrees;
  string nTrees_s;
  string nTrees_l;
  
  size_t nMaxLeaves;
  string nMaxLeaves_s;
  string nMaxLeaves_l;
  
  num_t shrinkage;
  string shrinkage_s;
  string shrinkage_l;
  
  num_t subSampleSize;
  string subSampleSize_s;
  string subSampleSize_l;

  GBT_options():
    isOptimizedNodeSplit(GBT_IS_OPTIMIZED_NODE_SPLIT),
    isOptimizedNodeSplit_s("b"),
    isOptimizedNodeSplit_l("GBT_optimize"),
    
    nTrees(GBT_DEFAULT_N_TREES),
    nTrees_s("r"),
    nTrees_l("GBT_ntrees"),
    
    nMaxLeaves(GBT_DEFAULT_N_MAX_LEAVES),
    nMaxLeaves_s("l"),
    nMaxLeaves_l("GBT_maxleaves"),
    
    shrinkage(GBT_DEFAULT_SHRINKAGE),
    shrinkage_s("z"),
    shrinkage_l("GBT_shrinkage"),

    subSampleSize(GBT_DEFAULT_SUB_SAMPLE_SIZE),
    subSampleSize_s("u"),
    subSampleSize_l("GBT_samplesize") {}
};

class Progress {
public:
  Progress(): width_(3) { cout << setw(width_) << "0" << "%" << flush; }
  ~Progress() { reset(); }

  void update(const num_t fraction) { reset(); cout << setw(width_) << static_cast<size_t>(fraction*100) << "%" << flush; }

private:
  
  void reset() { for(size_t i = 0; i <= width_; ++i) { cout << "\b"; } }

  size_t width_;

};


void printHeader() {
  cout << endl;
  cout << " ------------------------------------------------------- " << endl;
  cout << "|         RF-ACE v0.8.2, September 16th 2011            |" << endl;
  cout << "|                                                       |" << endl;
  cout << "|    Project page: http://code.google.com/p/rf-ace      |" << endl;
  cout << " ------------------------------------------------------- " << endl;
  cout << endl;
}

void printHelp(const General_options& geno, const RF_options& rfo, const GBT_options& gbto) {

  size_t maxwidth = 20;
  cout << endl;
  
  cout << "REQUIRED ARGUMENTS:" << endl;
  cout << " -" << geno.trainInput_s << " / --" << geno.trainInput_l << setw( maxwidth - geno.trainInput_l.size() )
       << " " << "Train data input file, associations will be sought from this data. Supported formats: AFM and ARFF" << endl;
  cout << " -" << geno.targetStr_s << " / --" << geno.targetStr_l << setw( maxwidth - geno.targetStr_l.size() )
       << " " << "Target, specified as integer or string that is to be matched with the content of input" << endl;
  cout << endl;

  cout << "OPTIONAL ARGUMENTS:" << endl;
  cout << " -" << geno.associationOutput_s << " / --" << geno.associationOutput_l << setw( maxwidth - geno.associationOutput_l.size() )
       << " " << "Association output file" << endl;
  cout << " -" << geno.testInput_s << " / --" << geno.testInput_l << setw( maxwidth - geno.testInput_l.size() )
       << " " << "Test data input file, predictions will be made from this data" << endl;
  cout << " -" << geno.predictionOutput_s << " / --" << geno.predictionOutput_l << setw( maxwidth - geno.predictionOutput_l.size() )
       << " " << "Prediction output file" << endl;  
  cout << " -" << geno.featureMaskInput_s << " / --" << geno.featureMaskInput_l << setw( maxwidth - geno.featureMaskInput_l.size() )
       << " " << "Feature mask input file. String of ones and zeroes, zeroes indicating removal of the feature in the matrix" << endl;
  cout << endl;
  
  cout << "OPTIONAL ARGUMENTS -- RANDOM FOREST:" << endl;
  cout << " -" << rfo.isOptimizedNodeSplit_s << " / --" << rfo.isOptimizedNodeSplit_l
       << setw( maxwidth - rfo.isOptimizedNodeSplit_l.size() )
       << " " << "Perform optimized node splitting with Random Forests (currently ENFORCED)" << endl;
  cout << " -" << rfo.nTrees_s << " / --" << rfo.nTrees_l << setw( maxwidth - rfo.nTrees_l.size() )
       << " " << "Number of trees per RF (default nSamples/realSampleFraction)" << endl;
  cout << " -" << rfo.mTry_s << " / --" << rfo.mTry_l << setw( maxwidth - rfo.mTry_l.size() )
       << " " << "Number of randomly drawn features per node split (default sqrt(nFeatures))" << endl;
  cout << " -" << rfo.nMaxLeaves_s << " / --" << rfo.nMaxLeaves_l << setw( maxwidth - rfo.nMaxLeaves_l.size() )
       << " " << "Maximum number of leaves per tree (default " << RF_DEFAULT_N_MAX_LEAVES << ")" << endl;
  cout << " -" << rfo.nodeSize_s << " / --" << rfo.nodeSize_l << setw( maxwidth - rfo.nodeSize_l.size() )
       << " " << "Minimum number of train samples per node, affects tree depth (default max{5,nSamples/100})" << endl;
  cout << " -" << rfo.nPerms_s << " / --" << rfo.nPerms_l << setw( maxwidth - rfo.nPerms_l.size() ) 
       << " " << "Number of Random Forests (default " << RF_DEFAULT_N_PERMS << ")" << endl;
  cout << endl;

  cout << "OPTIONAL ARGUMENTS -- FEATURE FILTER:" << endl;
  cout << " -" << geno.noFilter_s << " / --" << geno.noFilter_l << setw( maxwidth - geno.noFilter_l.size() )
       << " " << "Set this flag to turn RF filtering OFF (default ON)" << endl;
  cout << " -" << geno.pValueThreshold_s << " / --" << geno.pValueThreshold_l << setw( maxwidth - geno.pValueThreshold_l.size() )
       << " " << "p-value threshold below which associations are listed (default "
       << GENERAL_DEFAULT_P_VALUE_THRESHOLD << ")" << endl;
  cout << endl;
  
  cout << "OPTIONAL ARGUMENTS -- GRADIENT BOOSTING TREES:" << endl;
  cout << " -" << gbto.isOptimizedNodeSplit_s << " / --" << gbto.isOptimizedNodeSplit_l
       << setw( maxwidth - gbto.isOptimizedNodeSplit_l.size() ) 
       << " " << "Perform optimized node splitting with Gradient Boosting Trees (currently ENFORCED)" << endl;
  cout << " -" << gbto.nTrees_s << " / --" << gbto.nTrees_l << setw( maxwidth - gbto.nTrees_l.size() ) 
       << " " << "Number of trees in the GBT (default " << GBT_DEFAULT_N_TREES << ")" << endl; 
  cout << " -" << gbto.nMaxLeaves_s << " / --" << gbto.nMaxLeaves_l << setw( maxwidth - gbto.nMaxLeaves_l.size() ) 
       << " " << "Maximum number of leaves per tree (default " << GBT_DEFAULT_N_MAX_LEAVES << ")" << endl;
  cout << " -" << gbto.shrinkage_s << " / --" << gbto.shrinkage_l << setw( maxwidth - gbto.shrinkage_l.size() ) 
       << " " << "Shrinkage applied to evolving the residual (default " << GBT_DEFAULT_SHRINKAGE << ")" << endl;
  cout << " -" << gbto.subSampleSize_s << " / --" << gbto.subSampleSize_l << setw( maxwidth - gbto.subSampleSize_l.size() ) 
       << " " << "Sample size fraction for training the trees (default " << GBT_DEFAULT_SUB_SAMPLE_SIZE << ")" << endl;
  cout << endl;
  
  cout << "EXAMPLES:" << endl;
  cout << endl;

  cout << "List features associated with feature called \"AFFECTION\":" << endl;
  cout << "bin/rf_ace --traindata featurematrix.afm --target AFFECTION --associations associations.tsv" << endl << endl;

  cout << "Train the model for \"AFFECTION\" with \"original.afm\", and predict new data for \"AFFECTION\" with \"newdata.afm\":" << endl;
  cout << "bin/rf_ace --traindata original.afm --target AFFECTION --testdata newdata.afm --predictions.tsv" << endl << endl;
}

void printHelpHint() {
  cout << endl;
  cout << "To get started, type \"-h\" or \"--help\"" << endl;
}

void executeRandomForestFilter(Treedata& treedata,
                               const size_t targetIdx,
                               const RF_options& op,
                               vector<num_t>& pValues,
                               vector<num_t>& importanceValues);

void readFeatureMask(const string& fileName, const size_t nFeatures, vector<size_t>& keepFeatureIcs);


int main(const int argc, char* const argv[]) {

  General_options gen_op;
  RF_options RF_op_copy; // We need a copy (backup) since the options will differ between RF permutations, if enabled ...
  GBT_options GBT_op;

  //Print the intro header
  printHeader();

  //With no input arguments the help is printed
  if(argc == 1) {
    printHelp(gen_op,RF_op_copy,GBT_op);
    return(EXIT_SUCCESS);
  }

  //Read the user parameters ... 
  ArgParse parser(argc,argv);

  
  // first General Options
  parser.getFlag(gen_op.printHelp_s, gen_op.printHelp_l, gen_op.printHelp);
  parser.getArgument<string>(gen_op.trainInput_s, gen_op.trainInput_l, gen_op.trainInput); 
  parser.getArgument<string>(gen_op.targetStr_s, gen_op.targetStr_l, gen_op.targetStr); 
  parser.getArgument<string>(gen_op.associationOutput_s, gen_op.associationOutput_l, gen_op.associationOutput);
  parser.getArgument<string>(gen_op.featureMaskInput_s, gen_op.featureMaskInput_l, gen_op.featureMaskInput);
  parser.getArgument<string>(gen_op.testInput_s, gen_op.testInput_l, gen_op.testInput);
  parser.getArgument<string>(gen_op.predictionOutput_s, gen_op.predictionOutput_l, gen_op.predictionOutput);
  parser.getFlag(gen_op.noFilter_s, gen_op.noFilter_l, gen_op.noFilter);
  parser.getArgument<num_t>(gen_op.pValueThreshold_s, gen_op.pValueThreshold_l, gen_op.pValueThreshold);

  parser.getFlag(RF_op_copy.isOptimizedNodeSplit_s, RF_op_copy.isOptimizedNodeSplit_l, RF_op_copy.isOptimizedNodeSplit);
  parser.getArgument<size_t>(RF_op_copy.nTrees_s,RF_op_copy.nTrees_l,RF_op_copy.nTrees);
  parser.getArgument<size_t>(RF_op_copy.mTry_s, RF_op_copy.mTry_l, RF_op_copy.mTry); 
  parser.getArgument<size_t>(RF_op_copy.nMaxLeaves_s, RF_op_copy.nMaxLeaves_l, RF_op_copy.nMaxLeaves);
  parser.getArgument<size_t>(RF_op_copy.nodeSize_s, RF_op_copy.nodeSize_l, RF_op_copy.nodeSize); 
  parser.getArgument<size_t>(RF_op_copy.nPerms_s, RF_op_copy.nPerms_l, RF_op_copy.nPerms); 
  //parser.getArgument<num_t>(RF_op_copy.pValueThreshold_s, RF_op_copy.pValueThreshold_l, RF_op_copy.pValueThreshold); 

  parser.getFlag(GBT_op.isOptimizedNodeSplit_s, GBT_op.isOptimizedNodeSplit_l, GBT_op.isOptimizedNodeSplit);
  parser.getArgument<size_t>(GBT_op.nTrees_s, GBT_op.nTrees_l, GBT_op.nTrees);
  parser.getArgument<size_t>(GBT_op.nMaxLeaves_s, GBT_op.nMaxLeaves_l, GBT_op.nMaxLeaves);
  parser.getArgument<num_t>(GBT_op.shrinkage_s, GBT_op.shrinkage_l, GBT_op.shrinkage);
  parser.getArgument<num_t>(GBT_op.subSampleSize_s, GBT_op.subSampleSize_l, GBT_op.subSampleSize);
  
  bool makePrediction = ( gen_op.testInput != "" ) && (gen_op.predictionOutput != "" );
  
  if(gen_op.printHelp) {
    printHelp(gen_op,RF_op_copy,GBT_op);
    return(EXIT_SUCCESS);
  }

  bool writeAssociationsToFile = gen_op.associationOutput != "";
  bool writePredictionsToFile = gen_op.predictionOutput != "";

  //Print help and exit if input file is not specified
  if ( gen_op.trainInput == "" ) {
    cerr << "Input file not specified" << endl;
    printHelpHint();
    return EXIT_FAILURE;
  }

  //Print help and exit if target index is not specified
  if ( gen_op.targetStr == "" ) {
    cerr << "target(s) (-i/--target) not specified" << endl;
    printHelpHint();
    return(EXIT_FAILURE);
  }

  if ( gen_op.associationOutput == "" && gen_op.predictionOutput == "" ) {
    cerr << "No output files specified" << endl;
    printHelpHint();
    return(EXIT_FAILURE);
  }

  //Read data into Treedata object
  cout << "Reading file '" << gen_op.trainInput << "', please wait... " << flush;
  Treedata treedata_copy(gen_op.trainInput);
  cout << "DONE" << endl;

  int integer;
  if ( datadefs::isInteger(gen_op.targetStr,integer) ) {

    if ( integer < 0 || integer >= static_cast<int>( treedata_copy.nFeatures() ) ) {
      cerr << "Feature index (" << integer << ") must be within bounds 0 ... " << treedata_copy.nFeatures() - 1 << endl;
      return EXIT_FAILURE;
    }

    gen_op.targetStr = treedata_copy.getFeatureName( static_cast<size_t>( integer ) );
  }

  vector<size_t> keepFeatureIcs;
  if(gen_op.featureMaskInput != "") {
    cout << "Reading masking file '" << gen_op.featureMaskInput << "', please wait... " << flush;
    readFeatureMask(gen_op.featureMaskInput,treedata_copy.nFeatures(),keepFeatureIcs);
    treedata_copy.keepFeatures(keepFeatureIcs);
    cout << "DONE" << endl;
  }
  cout << endl;

  if ( RF_op_copy.mTry == RF_DEFAULT_M_TRY ) {
    RF_op_copy.mTry = static_cast<size_t>( 0.1 * treedata_copy.nFeatures() );
  }

  //Check which feature names match with the specified target identifier
  set<size_t> targetIcs;
  treedata_copy.getMatchingTargetIcs(gen_op.targetStr,targetIcs);
  if(targetIcs.size() == 0) {
    cerr << "No features match the specified target identifier '" << gen_op.targetStr << "'" << endl;
    if ( gen_op.featureMaskInput != "" ) {
      cerr << "NOTE: feature mask, being set, may have caused the target(s) to be erased from the feature matrix" << endl;
    }
    return EXIT_FAILURE;
  }

  if(gen_op.featureMaskInput != "" && targetIcs.size() > 1) {
    cout << "WARNING: feature mask is specified in the presence of multiple targets. All targets will be analyzed with the same mask set." << endl;
  }

  // NOTE: this is to override node split optimization settings
  bool noOptimization = false;
  parser.getFlag("q","optimization_off",noOptimization);
  if ( noOptimization) {
    cout << "TEST MODE: optimized node split turned OFF" << endl << endl;
    RF_op_copy.isOptimizedNodeSplit = false;
    GBT_op.isOptimizedNodeSplit = false;
  } else {
    RF_op_copy.isOptimizedNodeSplit = true;
    GBT_op.isOptimizedNodeSplit = true;
  }

  ofstream toAssociationFile(gen_op.associationOutput.c_str());
  toAssociationFile.precision(8);

  ofstream toPredictionFile(gen_op.predictionOutput.c_str());

  //Before starting number crunching, print values of parameters of RF-ACE
  int maxwidth = 19;
  cout << endl;
  cout << "RF-ACE parameter configuration:" << endl;
  cout << endl;
  cout << "General configuration:" << endl;
  cout << "  --" << gen_op.trainInput_l << setw( maxwidth - gen_op.trainInput_l.size() ) << ""
       << "= " << gen_op.trainInput << endl;
  cout << "  --" << gen_op.targetStr_l << setw( maxwidth - gen_op.targetStr_l.size() ) << ""
       << "= " << gen_op.targetStr << endl;
  cout << "  --" << gen_op.associationOutput_l << setw( maxwidth - gen_op.associationOutput_l.size() ) << ""
       << "= "; if ( writeAssociationsToFile ) { cout << gen_op.associationOutput << endl; } else { cout << "NOT SET" << endl; }
  cout << "  --" << gen_op.testInput_l << setw( maxwidth - gen_op.testInput_l.size() ) << ""
       << "= "; if( makePrediction ) { cout << gen_op.testInput << endl; } else { cout << "NOT SET" << endl; }
  cout << "  --" << gen_op.predictionOutput_l << setw( maxwidth - gen_op.predictionOutput_l.size() ) << ""
       << "= "; if( writePredictionsToFile ) { cout << gen_op.predictionOutput << endl; } else { cout << "NOT SET" << endl; }
  cout << endl;

  cout << "Random Forest configuration:" << endl;
  cout << "  --" << RF_op_copy.nTrees_l << setw( maxwidth - RF_op_copy.nTrees_l.size() ) << ""
       << "= "; if(RF_op_copy.nTrees == 0) { cout << "DEFAULT" << endl; } else { cout << RF_op_copy.nTrees << endl; }
  cout << "  --" << RF_op_copy.mTry_l << setw( maxwidth - RF_op_copy.mTry_l.size() ) << ""
       << "= "; if(RF_op_copy.mTry == 0) { cout << "DEFAULT" << endl; } else { cout << RF_op_copy.mTry << endl; }
  cout << "  --" << RF_op_copy.nMaxLeaves_l << setw( maxwidth - RF_op_copy.nMaxLeaves_l.size() ) << ""
       << "= " << RF_op_copy.nMaxLeaves << endl;
  cout << "  --" << RF_op_copy.nodeSize_l << setw( maxwidth - RF_op_copy.nodeSize_l.size() ) << ""
       << "= "; if(RF_op_copy.nodeSize == 0) { cout << "DEFAULT" << endl; } else { cout << RF_op_copy.nodeSize << endl; }
  cout << "  --" << RF_op_copy.nPerms_l << setw( maxwidth - RF_op_copy.nPerms_l.size() ) << ""
       << "= " << RF_op_copy.nPerms << endl;
  cout << endl;

  if ( !gen_op.noFilter ) {
    cout << "Feature filter ENABLED. Configuration:" << endl;
    cout << "  --pthresold          = " << gen_op.pValueThreshold << endl;
    cout << endl;
  } else {
    cout << "Feature filter DISABLED" << endl;
    cout << endl;
  }

  if ( makePrediction ) {
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

  //The program starts a loop in which an RF-ACE model will be built for each spcified target feature
  size_t iter = 1;
  for(set<size_t>::const_iterator it(targetIcs.begin()); it != targetIcs.end(); ++it, ++iter) {

    //Copy the data and options into new objects, which allows the program to alter the other copy without losing data
    Treedata treedata = treedata_copy;
    RF_options RF_op = RF_op_copy;

    //Extract the target index from the pointer and the number of real samples from the treedata object
    size_t targetIdx = *it;
    string targetName = treedata.getFeatureName(targetIdx);
    size_t nRealSamples = treedata.nRealSamples(targetIdx);
    num_t realFraction = 1.0*nRealSamples / treedata.nSamples();

    size_t maxwidth = 1 + static_cast<int>(targetIcs.size()) / 10;

    cout << "== " << setw(maxwidth) << iter << "/" << setw(maxwidth) << targetIcs.size() 
	 << " target " << treedata.getFeatureName(targetIdx) << ", " << flush;
    
    //If the target has no real samples, the program will just exit
    if(nRealSamples == 0) {
      cout << "Omitting: it has no real samples." << endl;
      continue;
    }

    if(treedata.isFeatureNumerical(targetIdx)) {
      cout << "regression ";
    } else {
      cout << treedata.nCategories(targetIdx) << "-class ";
    }
    cout << "CARTs. " << nRealSamples << " / " << treedata.nSamples() << " samples ( " << 100 * ( 1 - realFraction ) << " % missing )" << endl;
            
    //If default mTry is to be used...
    if(RF_op.mTry == RF_DEFAULT_M_TRY) {
      RF_op.mTry = static_cast<size_t>( floor( sqrt( 1.0*treedata.nFeatures() ) ) );   
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

    //cout << "== " << RF_op.nPerms << " RFs; " << RF_op.nTrees << " trees per RF; " 
    //	 << RF_op.mTry << " features tested per split; minimum node size of " << RF_op.nodeSize << endl;
    
    ////////////////////////////////////////////////////////////////////////
    //  STEP 1 -- MULTIVARIATE ASSOCIATIONS WITH RANDOM FOREST ENSEMBLES  //
    ////////////////////////////////////////////////////////////////////////     
    vector<num_t> pValues; //(treedata.nFeatures());
    vector<num_t> importanceValues; //(treedata.nFeatures());
    cout << "    => Uncovering associations... " << flush;
    executeRandomForestFilter(treedata,targetIdx,RF_op,pValues,importanceValues);
    cout << "DONE" << endl;

    //////////////////////////////////////////////////////////////
    //  STEP 2 ( OPTIONAL ) -- FEATURE FILTERING WITH P-VALUES  //
    //////////////////////////////////////////////////////////////
    

    vector<size_t> keepFeatureIcs;    
    // If filtering is to be applied ...
    if(!gen_op.noFilter) {
      cout << "    => Filtering features... " << flush;
      size_t nFeatures = treedata.nFeatures();

      //Store removed features and their indices here, just in case
      vector<string> removedFeatures;
      vector<size_t> removedFeatureIcs;

      // So we've already kept the target, hence start counting from 1
      size_t nKeepFeatures = 0;

      // Go through each feature, and keep those having p-value higher than the threshold. 
      // Save the kept and removed features, and remember to accumulate the counter
      for ( size_t featureIdx = 0; featureIdx < nFeatures; ++featureIdx ) {
	
	if (featureIdx == targetIdx || pValues[featureIdx] <= gen_op.pValueThreshold ) {
	  keepFeatureIcs.push_back(featureIdx);
	  pValues[nKeepFeatures] = pValues[featureIdx];
	  importanceValues[nKeepFeatures] = importanceValues[featureIdx];
	  ++nKeepFeatures;
	} else {
	  removedFeatureIcs.push_back(featureIdx);
	  removedFeatures.push_back(treedata.getFeatureName(featureIdx));
	}
      }

      // Resize containers
      treedata.keepFeatures( keepFeatureIcs );
      pValues.resize( keepFeatureIcs.size() );
      importanceValues.resize ( keepFeatureIcs.size() );
      
      // Reassign target index variable point to the target, which is the last feature in the list
      targetIdx = keepFeatureIcs.size() - 1;

      // Print some statistics
      cout << "DONE, " << treedata.nFeatures() << " / " << treedata_copy.nFeatures() << " features ( "
	   << 100.0 * treedata.nFeatures() / treedata_copy.nFeatures() << " % ) left " << endl;
    }

    ///////////////////////////////////////////////////////////////////////////
    //  STEP 3 ( OPTIONAL ) -- DATA PREDICTION WITH GRADIENT BOOSTING TREES  //
    ///////////////////////////////////////////////////////////////////////////
    if( makePrediction && writePredictionsToFile ) {      
      cout << "    => Predicting... " << flush;
      
      StochasticForest SF(&treedata,targetIdx,GBT_op.nTrees);
      SF.learnGBT(GBT_op.nMaxLeaves, GBT_op.shrinkage, GBT_op.subSampleSize);
      
      //StochasticForest SF(&treedata,targetIdx,10000);
      //SF.learnRF(static_cast<size_t>(sqrt(1.0*treedata.nFeatures())),5,false,true);

      Treedata treedata_test(gen_op.testInput);
      if ( !gen_op.noFilter ) {
	treedata_test.keepFeatures(keepFeatureIcs);
      }

      // An assertion to check that the train and test matrices have features in the same order
      assert(treedata.nFeatures() == treedata_test.nFeatures());
      for ( size_t i = 0; i < treedata.nFeatures(); ++i ) {
	assert(treedata.getFeatureName(i) == treedata_test.getFeatureName(i));
      }

      vector<num_t> confidence;
      if(treedata_test.isFeatureNumerical(targetIdx)) {

	vector<num_t> prediction;
	SF.predict(&treedata_test,prediction,confidence);

	for(size_t i = 0; i < prediction.size(); ++i) {
	  toPredictionFile << targetName.c_str() << "\t" << "sampleID" << "\t" << treedata_test.getRawFeatureData(targetIdx,i) 
			   << "\t" << prediction[i] << "\t" << setprecision(3) << confidence[i] << endl;
	  
	}
       
      } else {
	
	vector<string> prediction;
	SF.predict(&treedata_test,prediction,confidence);
		
	for(size_t i = 0; i < prediction.size(); ++i) {
	  toPredictionFile << targetName.c_str() << "\t" << "sampleID" << "\t" << treedata_test.getRawFeatureData(targetIdx,i) 
			   << "\t" << prediction[i] << "\t" << setprecision(3) << confidence[i] << endl;
	}
	
      }
     
      cout << "DONE" << endl;
    }
    cout << endl;
            
    vector<size_t> refIcs( treedata.nFeatures() );
    bool isIncreasingOrder = true;
    datadefs::sortDataAndMakeRef(isIncreasingOrder,pValues,refIcs); // BUG
    datadefs::sortFromRef<num_t>(importanceValues,refIcs);

    if( writeAssociationsToFile ) {
    
      for ( size_t featureIdx = 0; featureIdx < treedata.nFeatures(); ++featureIdx ) {
        
	if ( pValues[featureIdx] > gen_op.pValueThreshold ) {
          continue;
	}
        
        if ( refIcs[featureIdx] == targetIdx ) {
          continue;
        }
        
        if ( RF_op.nPerms > 1 ) {
	  toAssociationFile << fixed << targetName.c_str() << "\t" << treedata.getFeatureName(refIcs[featureIdx]).c_str() 
			    << "\t" << pValues[featureIdx] << "\t" << importanceValues[featureIdx] << "\t"
                            << treedata.pearsonCorrelation(targetIdx,refIcs[featureIdx]) << endl;
	} else {
          toAssociationFile << fixed << targetName.c_str() << "\t" << treedata.getFeatureName(refIcs[featureIdx]).c_str() 
			    << "\tNA\t" << importanceValues[featureIdx] << "\t"
			    << treedata.pearsonCorrelation(targetIdx,refIcs[featureIdx]) << endl;
        }    
      }
    }
  }
  
  toAssociationFile.close();
  toPredictionFile.close();

  if ( writeAssociationsToFile ) {
    cout << "Association file '" << gen_op.associationOutput << "' created. Format:" << endl;
    cout << "TARGET   PREDICTOR   P-VALUE   IMPORTANCE   CORRELATION" << endl;
    cout << endl;
  }

  if ( writePredictionsToFile ) {
    cout << "Prediction file '" << gen_op.predictionOutput << "' created. Format:" << endl;
    cout << "TARGET   SAMPLE_ID   DATA      PREDICTION   CONFIDENCE" << endl; 
    cout << endl;
  }
  
  cout << "RF-ACE completed successfully." << endl;
  cout << endl;
      
  return(EXIT_SUCCESS);
}


void executeRandomForestFilter(Treedata& treedata,
                               const size_t targetIdx,
                               const RF_options& op,
                               vector<num_t>& pValues,
                               vector<num_t>& importanceValues) {

  vector<vector<num_t> > importanceMat(op.nPerms);
  pValues.resize(treedata.nFeatures());
  importanceValues.resize(treedata.nFeatures());
  size_t nNodesInAllForests = 0;

  Progress progress;
  for(int permIdx = 0; permIdx < static_cast<int>(op.nPerms); ++permIdx) {
    //cout << "  RF " << permIdx + 1 << ": ";
    //Treedata td_thread = treedata;

    progress.update(1.0*permIdx/op.nPerms);
  
    bool useContrasts;
    if(op.nPerms > 1) {
      useContrasts = true;
    } else {
      useContrasts = false;
    }

    StochasticForest SF(&treedata,targetIdx,op.nTrees);
    SF.learnRF(op.mTry,op.nMaxLeaves,op.nodeSize,useContrasts,op.isOptimizedNodeSplit);
    size_t nNodesInForest = SF.nNodes();
    nNodesInAllForests += nNodesInForest;
    importanceMat[permIdx] = SF.featureImportance();
    //printf("  RF %i: %i nodes (avg. %6.3f nodes/tree)\n",permIdx+1,static_cast<int>(nNodesInForest),1.0*nNodesInForest/op.nTrees);
    progress.update( 1.0 * ( 1 + permIdx ) / op.nPerms );
  }

  if(op.nPerms > 1) {
    for(size_t featureIdx = 0; featureIdx < treedata.nFeatures(); ++featureIdx) {

      size_t nRealSamples;
      vector<num_t> fSample(op.nPerms);
      vector<num_t> cSample(op.nPerms);
      for(size_t permIdx = 0; permIdx < op.nPerms; ++permIdx) {
        fSample[permIdx] = importanceMat[permIdx][featureIdx];
        cSample[permIdx] = importanceMat[permIdx][featureIdx + treedata.nFeatures()];
      }
      
      datadefs::utest(fSample,cSample,pValues[featureIdx]);
      datadefs::mean(fSample,importanceValues[featureIdx],nRealSamples);
    }
  } else {
    importanceValues = importanceMat[0];
  }

  importanceValues.resize(treedata.nFeatures());
  
}

void readFeatureMask(const string& fileName, const size_t nFeatures, vector<size_t>& keepFeatureIcs) {

  ifstream featurestream;
  featurestream.open(fileName.c_str());
  assert(featurestream.good());
  
  string maskAsStr;
  getline(featurestream,maskAsStr);
  
  assert(maskAsStr.size() == nFeatures);

  keepFeatureIcs.clear();

  for(size_t featureIdx = 0; featureIdx < nFeatures; ++featureIdx) {

    if ( maskAsStr[featureIdx] == '1' ) {
      keepFeatureIcs.push_back(featureIdx);
    } else if ( maskAsStr[featureIdx] != '0' ) { 
      cerr << "Mask file formatted incorrectly, must contain only 0's and 1's" << endl;
      assert(false);
    }
    
  }

}

