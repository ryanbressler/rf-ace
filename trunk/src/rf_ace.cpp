#include <cstdlib>
#include <cassert>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>
#include <cmath>
#include <stdio.h>
//#include <omp.h>

#include "argparse.hpp"
#include "randomforest.hpp"
#include "GBT.hpp"
#include "treedata.hpp"
#include "datadefs.hpp"

using namespace std;
using datadefs::num_t;

const bool DEFAULT_PRINT_HELP = false;
const size_t DEFAULT_NTREES = 0;
const size_t DEFAULT_RF_MTRY = 0;
const size_t DEFAULT_RF_NODESIZE = 0;
const size_t DEFAULT_RF_NPERMS = 1;
const num_t DEFAULT_RF_PVALUETHRESHOLD = 0.10;
const bool DEFAULT_APPLY_RF_FILTER = false;
const bool DEFAULT_OPTIMIZED_NODE_SPLIT = false;
const bool DEFAULT_ENABLE_GBT = false; // TEMPORARY VARIABLE
const size_t DEFAULT_GBT_N_MAX_LEAVES = 6;
const num_t DEFAULT_GBT_SHRINKAGE = 0.5;
const num_t DEFAULT_GBT_SUB_SAMPLE_SIZE = 0.5;


void printHeader()
{
  cout << endl;
  cout << " --------------------------------------------------------------- " << endl;
  cout << "| RF-ACE -- efficient feature selection with heterogeneous data |" << endl;
  cout << "|                                                               |" << endl;
  cout << "|  Version:      RF-ACE v0.7.0, July 24th, 2011                 |" << endl;
  cout << "|  Project page: http://code.google.com/p/rf-ace                |" << endl;
  cout << "|  Contact:      timo.p.erkkila@tut.fi                          |" << endl;
  cout << "|                kari.torkkola@gmail.com                        |" << endl;
  cout << "|                                                               |" << endl;
  cout << "|              DEVELOPMENT VERSION, BUGS EXIST!                 |" << endl;
  cout << " --------------------------------------------------------------- " << endl;
}

void printHelp()
{
  cout << endl;
  cout << "REQUIRED ARGUMENTS:" << endl;
  cout << " -I / --input         Input feature file (AFM or ARFF)" << endl;
  cout << " -i / --target        Target, specified as integer or string that is to be matched with the content of input" << endl;
  cout << " -O / --output        Output association file" << endl;
  cout << endl;
  cout << "OPTIONAL ARGUMENTS:" << endl;
  cout << " -n / --ntrees        (RF+GBT) Number of trees per RF (default 10*nsamples/nrealsamples)" << endl;
  cout << " -m / --mtry          (RF) Number of randomly drawn features per node split (default sqrt(nfeatures))" << endl;
  cout << " -s / --nodesize      (RF) Minimum number of train samples per node, affects tree depth (default max{5,nsamples/100})" << endl;
  cout << " -p / --nperms        (RF) Number of Random Forests (default " << DEFAULT_RF_NPERMS << ")" << endl;
  cout << " -t / --pthreshold    (RF) p-value threshold below which associations are listed (default " << DEFAULT_RF_PVALUETHRESHOLD << ")" << endl;
  cout << " -l / --maxleaves     (GBT) Maximum number of leaves per tree (default " << DEFAULT_GBT_N_MAX_LEAVES << ")" << endl;
  cout << " -z / --shrinkage     (GBT) Shrinkage applied to evolving the residual (default " << DEFAULT_GBT_SHRINKAGE << ")" << endl;
  cout << " -u / --subsamplesize (GBT) Sample size fraction for training the trees (default " << DEFAULT_GBT_SUB_SAMPLE_SIZE << ")" << endl;
  cout << endl;
  cout << "OPTIONAL HANDLES, SET TO *FALSE* BY DEFAULT:" << endl;
  cout << " -f / --filterRF      Apply feature filtering with Random Forests" << endl;
  cout << " -o / --optimizedRF   Perform optimized node splitting with Random Forests" << endl;
  cout << " -g / --enableGBT     Enable Gradient Boosting Trees" << endl;
  cout << endl;
}

void printHelpHint()
{
  cout << endl;
  cout << "To get started, type \"-h\" or \"--help\"" << endl;
}

struct Options {

  // IO
  bool printHelpHandle;
  string input;
  string output;
  
  // GENERAL
  string targetStr;
  size_t nTrees;
  bool applyFilter;
  bool isOptimizedNodeSplit;
  
  // RF
  size_t mTry;
  size_t nodeSize;
  size_t nPerms;
  num_t pValueThreshold;
  
  // GBT
  size_t nMaxLeaves;
  num_t shrinkage;
  num_t subSampleSize;

  Options():   
    printHelpHandle(DEFAULT_PRINT_HELP),
    input(""),
    output(""),
    targetStr(""),
    nTrees(DEFAULT_NTREES),
    applyFilter(DEFAULT_APPLY_RF_FILTER),
    isOptimizedNodeSplit(DEFAULT_OPTIMIZED_NODE_SPLIT), 
    mTry(DEFAULT_RF_MTRY),
    nodeSize(DEFAULT_RF_NODESIZE),
    nPerms(DEFAULT_RF_NPERMS),
    pValueThreshold(DEFAULT_RF_PVALUETHRESHOLD),
    nMaxLeaves(DEFAULT_GBT_N_MAX_LEAVES),
    shrinkage(DEFAULT_GBT_SHRINKAGE),
    subSampleSize(DEFAULT_GBT_SUB_SAMPLE_SIZE) { /* EMPTY CONSTRUCTOR */ }
};

void executeRandomForestFilter(Treedata& treedata,
                               const size_t targetIdx,
                               const Options& op,
                               vector<num_t>& pValues,
                               vector<num_t>& importanceValues);


int main(int argc, char* argv[])
{

  Options op_copy;

  //Print the intro tab
  printHeader();

  //With no input arguments the help is printed
  if(argc == 1)
    {
      printHelp();
      return(EXIT_SUCCESS);
    }

  //Set user parameters to default values
  op_copy.printHelpHandle = DEFAULT_PRINT_HELP;
  op_copy.input = "";
  op_copy.targetStr = "";
  op_copy.nTrees = DEFAULT_NTREES;
  op_copy.mTry = DEFAULT_RF_MTRY;
  op_copy.nodeSize = DEFAULT_RF_NODESIZE;
  op_copy.nPerms = DEFAULT_RF_NPERMS;
  op_copy.pValueThreshold = DEFAULT_RF_PVALUETHRESHOLD;
  op_copy.applyFilter = DEFAULT_APPLY_RF_FILTER;
  op_copy.isOptimizedNodeSplit = DEFAULT_OPTIMIZED_NODE_SPLIT;
  op_copy.output = "";
  bool enableGBT = DEFAULT_ENABLE_GBT;
  op_copy.nMaxLeaves = DEFAULT_GBT_N_MAX_LEAVES;
  op_copy.shrinkage = DEFAULT_GBT_SHRINKAGE;
  op_copy.subSampleSize = DEFAULT_GBT_SUB_SAMPLE_SIZE;

  //Read the user parameters 
  ArgParse parser(argc,argv);
  parser.getArgument<bool>("h","help",op_copy.printHelpHandle);
  parser.getArgument<string>("I","input",op_copy.input); 
  parser.getArgument<string>("i","target",op_copy.targetStr); 
  parser.getArgument<size_t>("n","ntrees",op_copy.nTrees); 
  parser.getArgument<size_t>("m","mtry",op_copy.mTry); 
  parser.getArgument<size_t>("s","nodesize",op_copy.nodeSize); 
  parser.getArgument<size_t>("p","nperms",op_copy.nPerms); 
  parser.getArgument<num_t>("t","pthreshold",op_copy.pValueThreshold); 
  parser.getArgument<bool>("f","filterRF",op_copy.applyFilter); 
  parser.getArgument<bool>("o","optimizedRF",op_copy.isOptimizedNodeSplit);
  parser.getArgument<string>("O","output",op_copy.output);
  parser.getArgument<bool>("g","enableGBT",enableGBT);
  parser.getArgument<size_t>("l","maxleaves",op_copy.nMaxLeaves);
  parser.getArgument<num_t>("z","shrinkage",op_copy.shrinkage);
  parser.getArgument<num_t>("u","subsamplesize",op_copy.subSampleSize);

  if(op_copy.printHelpHandle)
    {
      printHelp();
      return(EXIT_SUCCESS);
    }

  //Print help and exit if input file is not specified
  if(op_copy.input == "")
    {
      cerr << "Input file not specified" << endl;
      printHelpHint();
      return EXIT_FAILURE;
    }

  //Print help and exit if target index is not specified
  if(op_copy.targetStr == "")
    {
      cerr << "target(s) (-i/--target) not specified" << endl;
      printHelpHint();
      return EXIT_FAILURE;
    }

  //Print help and exit if output file is not specified
  if(op_copy.output == "")
    {
      cerr << "Output file not specified" << endl;
      printHelpHint();
      return EXIT_FAILURE;
    }

  //Read data into Treedata object
  cout << endl << "Reading file '" << op_copy.input << "'" << endl;
  Treedata treedata_copy(op_copy.input);

  //Check which feature names match with the specified target identifier
  set<size_t> targetIcs;
  treedata_copy.getMatchingTargetIcs(op_copy.targetStr,targetIcs);
  if(targetIcs.size() == 0)
    {
      cerr << "No features match the specified target identifier '" << op_copy.targetStr << "'" << endl;
      return EXIT_FAILURE;
    }

  //The program starts a loop in which an RF-ACE model will be built for each spcified target feature
  size_t iter = 1;
  for(set<size_t>::const_iterator it(targetIcs.begin()); it != targetIcs.end(); ++it, ++iter)
    {
      //Copy the data and options into new objects, which allows the program to alter the other copy without losing data
      Treedata treedata = treedata_copy;
      Options op = op_copy;

      //Extract the target index from the pointer and the number of real samples from the treedata object
      size_t targetIdx = *it;
      size_t nRealSamples = treedata.nRealSamples(targetIdx);
      
      //If the target has no real samples, the program will just exit
      if(nRealSamples == 0)
	{
	  cout << "Omitting target '" <<  treedata.getFeatureName(targetIdx) << "', for it has no real samples..." << endl;
	  continue;
	}
      
      num_t realFraction = 1.0*nRealSamples / treedata.nSamples();
      
      //If default nTrees is to be used...
      if(op.nTrees == DEFAULT_NTREES)
	{
	  op.nTrees = static_cast<size_t>( 10.0*treedata.nSamples() / realFraction );
	}
      
      //If default mTry is to be used...
      if(op.mTry == DEFAULT_RF_MTRY)
	{
	  op.mTry = static_cast<size_t>( floor( sqrt( 1.0*treedata.nFeatures() ) ) );   
	}
      
      //If default nodeSize is to be used...
      if(op.nodeSize == DEFAULT_RF_NODESIZE)
	{
	  op.nodeSize = 5;
	  size_t altNodeSize = static_cast<size_t>( ceil( 1.0*nRealSamples/100 ) );
	  if(altNodeSize > op.nodeSize)
	    {
	      op.nodeSize = altNodeSize;
	    }
	}
      
      //Print values of parameters of RF-ACE
      cout << endl;
      cout << "RF-ACE parameter configuration:" << endl;
      cout << "  --input              = " << op.input << endl;
      cout << "  --nsamples           = " << nRealSamples << " / " << treedata.nSamples() << " (" << 100 * ( 1 - realFraction )<< "% missing)" << endl;
      cout << "  --nfeatures          = " << treedata.nFeatures() << endl;
      cout << "  --targetidx          = " << targetIdx << ", header '" << treedata.getFeatureName(targetIdx) << "'" << endl;
      cout << "  --ntrees             = " << op.nTrees << endl;
      cout << "  --(RF)mtry           = " << op.mTry << endl;
      cout << "  --(RF)nodesize       = " << op.nodeSize << endl;
      cout << "  --(RF)nperms         = " << op.nPerms << endl;
      cout << "  --(RF)pthresold      = " << op.pValueThreshold << endl;
      cout << "  --filterRF           = "; if(op.applyFilter) { cout << "YES" << endl; } else { cout << "NO" << endl; }
      cout << "  --optimizedRF        = "; if(op.isOptimizedNodeSplit) { cout << "YES" << endl; } else { cout << "NO" << endl; }
      cout << "  --enableGBT          = "; if(enableGBT) { cout << "YES" << endl;} else { cout << "NO" << endl; }
      cout << "  --(GBT)maxleaves     = " << op.nMaxLeaves << endl;
      cout << "  --(GBT)shrinkage     = " << op.shrinkage << endl;
      cout << "  --(GBT)subsamplesize = " << op.subSampleSize << endl;
      cout << "  --output             = " << op.output << endl << endl;
      
      if(!op.applyFilter)
	{
	  // NOTE: THIS WILL TAKE EFFECT AFTER GBT IS IN PLACE
	  //if(optimizedRF)
	  //  {
	  //    cout << "NOTE: setting filter optimizations ON does not help when the filter is OFF" << endl;
	  //  }
	}
      
      if(treedata.nFeatures() < op.mTry)
	{
	  cerr << "Not enough features (" << treedata.nFeatures()-1 << ") to test with mtry = " << op.mTry << " features per split" << endl;
	  return EXIT_FAILURE;
	}
      
      if(treedata.nSamples() < 2 * op.nodeSize)
	{
	  cerr << "Not enough samples (" << treedata.nSamples() << ") to perform a single split" << endl;
	  return EXIT_FAILURE;
	}

      
      //////////////////////////////////////////////////////////////////////
      //  ANALYSIS 1 -- Random Forest ensemble with artificial contrasts  //
      //////////////////////////////////////////////////////////////////////
      
      vector<num_t> pValues; //(treedata.nFeatures());
      vector<num_t> importanceValues; //(treedata.nFeatures());
      
      if(op.applyFilter)
	{
	  cout << "Applying Random Forest filter..." << endl; 
	  executeRandomForestFilter(treedata,targetIdx,op,pValues,importanceValues);
	   
	  size_t nFeatures = treedata.nFeatures();
	  vector<size_t> keepFeatureIcs(1);
	  keepFeatureIcs[0] = targetIdx;
	  vector<string> removedFeatures;
	  vector<size_t> removedFeatureIcs;
	  for(size_t featureIdx = 0; featureIdx < nFeatures; ++featureIdx)
	    {
	      if(featureIdx != targetIdx && importanceValues[featureIdx] > datadefs::EPS)
		{
		  keepFeatureIcs.push_back(featureIdx);
		}
	      else
		{
		  removedFeatureIcs.push_back(featureIdx);
		  removedFeatures.push_back(treedata.getFeatureName(featureIdx));
		}
	    }
	  
	  treedata.keepFeatures(keepFeatureIcs);
	  targetIdx = 0;
	  cout << endl;
	  cout << "Features filtered: " << removedFeatures.size() << " (" << treedata.nFeatures() << " left)" << endl;
	  cout << "TEST: target is '" << treedata.getFeatureName(targetIdx) << "'" << endl;
	}
      
      // THIS WILL BE REPLACED BY GBT
      executeRandomForestFilter(treedata,targetIdx,op,pValues,importanceValues);
    	  
      
      /////////////////////////////////////////////
      //  ANALYSIS 2 -- Gradient Boosting Trees  //
      /////////////////////////////////////////////
      
      if(enableGBT)
	{
	  cout << "Gradient Boosting Trees *ENABLED* (NON-FUNCTIONAL, WILL NOT AFFECT OUTPUT)" << endl;
	  
	  
	  GBT gbt(&treedata, targetIdx, op.nTrees, op.nMaxLeaves, op.shrinkage, op.subSampleSize);
	  
	  
	  cout <<endl<< "PREDICTION:" << endl;
	  vector<num_t> prediction(treedata.nSamples());
	  gbt.predictForest(&treedata, prediction);
	}
      else
      	{
	  cout << "Gradient Boosting Trees *DISABLED*" << endl;
	}
      
      
      
      
      ///////////////////////
      //  GENERATE OUTPUT  //
      ///////////////////////  
      vector<size_t> refIcs(treedata.nFeatures());
      //vector<string> fnames = treedata.featureheaders();
      bool isIncreasingOrder = false;
      datadefs::sortDataAndMakeRef(isIncreasingOrder,importanceValues,refIcs); // BUG
      datadefs::sortFromRef<num_t>(pValues,refIcs);
      //targetIdx = refIcs[targetIdx];
      
      string targetName = treedata.getFeatureName(targetIdx);
      
      //MODIFICATION: ALL ASSOCIATIONS WILL BE PRINTED
      if(true)
	{
	  //cout << "Writing associations to file '" << output << "'" << endl;
	  //ofstream os(output.c_str());
	  FILE* file;
	  if(iter == 1)
	    {
	      file = fopen(op.output.c_str(),"w");
	    }
	  else
	    {
	      file = fopen(op.output.c_str(),"a");
	    }
	  
	  for(size_t featureIdx = 0; featureIdx < treedata.nFeatures(); ++featureIdx)
	    {
	      
	      //MODIFICATION: ALL ASSOCIATIONS WILL BE PRINTED
	      /*
		if(importanceValues[featureIdx] < datadefs::EPS)
		{
		break;
		}
	      */
	      
	      if(refIcs[featureIdx] == targetIdx)
		{
		  //cout << refIcs[featureIdx] << " == " << targetIdx << " (" << targetHeader << ")" << endl;
		  continue;
		}
	      
	      if(op.nPerms > 1)
		{
		  fprintf(file,"%s\t%s\t%9.8f\t%9.8f\t%9.8f\n",targetName.c_str(),treedata.getFeatureName(refIcs[featureIdx]).c_str(),pValues[featureIdx],importanceValues[featureIdx],treedata.pearsonCorrelation(targetIdx,refIcs[featureIdx]));
		}
	      else
		{
		  fprintf(file,"%s\t%s\tNaN\t%9.8f\t%9.8f\n",targetName.c_str(),treedata.getFeatureName(refIcs[featureIdx]).c_str(),importanceValues[featureIdx],treedata.pearsonCorrelation(targetIdx,refIcs[featureIdx]));
		}
	      //os << target_str << "\t" << treedata.get_featureheader(ref_ics[i]) << "\t" 
	      //   << pvalues[i] << "\t" << ivalues[i] << "\t" << treedata.corr(targetidx,ref_ics[i]) << endl;
	    }
	  
	  //CURRENTLY NOT POSSIBLE TO REPORT FILTERED FEATURES
	  /*
	    if(reportFiltered)
	    {
	    for(size_t featureIdx = 0; featureIdx < removedFeatureIcs.size(); ++featureIdx)
	    {
	    fprintf(po,"%s\t%s\tNaN\tNaN\t%9.8f\n",targetName.c_str(),removedFeatures[featureIdx].c_str(),treedata.pearsonCorrelation(oldTargetIdx,removedFeatureIcs[featureIdx]));
	    }
	    }
	  */


	  fclose(file);
	  cout << endl << "Association file created. Format:" << endl;
	  cout << "TARGET   PREDICTOR   P-VALUE   IMPORTANCE   CORRELATION" << endl << endl;
	  cout << "Done." << endl;
	}
      else
	{
	  cout << endl << "No significant associations found..." << endl;
	}
    }
      
      
  return(EXIT_SUCCESS);
}


void executeRandomForestFilter(Treedata& treedata,
                               const size_t targetIdx,
                               const Options& op,
                               vector<num_t>& pValues,
                               vector<num_t>& importanceValues)
{

  vector<vector<num_t> > importanceMat(op.nPerms);
  pValues.resize(treedata.nFeatures());
  importanceValues.resize(treedata.nFeatures());
  size_t nNodesInAllForests = 0;

  clock_t time_start(clock());

  cout << "Growing " << op.nPerms << " Random Forests (RFs), please wait..." << endl;
  //#pragma omp parallel for
  for(int permIdx = 0; permIdx < static_cast<int>(op.nPerms); ++permIdx)
    {
      //cout << "  RF " << permIdx + 1 << ": ";
      //Treedata td_thread = treedata;
      bool useContrasts = false;
      Randomforest RF(&treedata,targetIdx,op.nTrees,op.mTry,op.nodeSize,useContrasts,op.isOptimizedNodeSplit);
      size_t nNodesInForest = RF.nNodes();
      nNodesInAllForests += nNodesInForest;
      importanceMat[permIdx] = RF.featureImportance();
      printf("  RF %i: %i nodes (avg. %6.3f nodes/tree)\n",permIdx+1,static_cast<int>(nNodesInForest),1.0*nNodesInForest/op.nTrees);
    }

  num_t time_diff = 1.0*(clock() - time_start) / CLOCKS_PER_SEC;
  cout << op.nPerms << " RFs, " << op.nPerms * op.nTrees << " trees, and " << nNodesInAllForests
       << " nodes generated in " << time_diff << " seconds (" << 1.0*nNodesInAllForests / time_diff
       << " nodes per second)" << endl;
  
  if(op.nPerms > 1)
    {
      for(size_t featureIdx = 0; featureIdx < treedata.nFeatures(); ++featureIdx)
	{
	  
	  size_t nRealSamples;
	  vector<num_t> fSample(op.nPerms);
	  vector<num_t> cSample(op.nPerms);
	  for(size_t permIdx = 0; permIdx < op.nPerms; ++permIdx)
	    {
	      fSample[permIdx] = importanceMat[permIdx][featureIdx];
	      cSample[permIdx] = importanceMat[permIdx][featureIdx + treedata.nFeatures()];
	    }
	  datadefs::utest(fSample,cSample,pValues[featureIdx]);
	  datadefs::meanVals(fSample,importanceValues[featureIdx],nRealSamples);
	}
    }
  else
    {
      importanceValues = importanceMat[0];
    }

  importanceValues.resize(treedata.nFeatures());
  
}
