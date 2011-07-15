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
const size_t DEFAULT_MTRY = 0;
const size_t DEFAULT_NODESIZE = 0;
const size_t DEFAULT_NPERMS = 1;
const num_t DEFAULT_PVALUETHRESHOLD = 0.10;
const bool DEFAULT_APPLY_FILTER = false;

void printHeader()
{
  cout << endl;
  cout << " --------------------------------------------------------------- " << endl;
  cout << "| RF-ACE -- efficient feature selection with heterogeneous data |" << endl;
  cout << "|                                                               |" << endl;
  cout << "|  Version:      RF-ACE v0.6.2, July 15th, 2011                 |" << endl;
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
  cout << " -I / --input        input feature file (AFM or ARFF)" << endl;
  cout << " -i / --target       target, specified as integer or string that is to be matched with the content of input" << endl;
  cout << " -O / --output       output association file" << endl;
  cout << endl;
  cout << "OPTIONAL ARGUMENTS:" << endl;
  cout << " -n / --ntrees       number of trees per RF (default 10*nsamples/nrealsamples)" << endl;
  cout << " -m / --mtry         number of randomly drawn features per node split (default sqrt(nfeatures))" << endl;
  cout << " -s / --nodesize     minimum number of train samples per node, affects tree depth (default max{5,nsamples/100})" << endl;
  cout << " -p / --nperms       number of Random Forests (default " << DEFAULT_NPERMS << ")" << endl;
  cout << " -t / --pthreshold   p-value threshold below which associations are listed (default " << DEFAULT_PVALUETHRESHOLD << ")" << endl;
  cout << endl;
  cout << "OPTIONAL HANDLES:" << endl;
  cout << " -f / --filter       Apply feature filtering with Random Forests" << endl;
}

void printHelpHint()
{
  cout << endl;
  cout << "To get started, type \"-h\" or \"--help\"" << endl;
}

void executeRandomForestFilter(Treedata& treedata,
			       size_t targetIdx,
			       size_t nTrees,
			       size_t mTry,
			       size_t nodeSize,
			       size_t nPerms,
			       vector<num_t>& pValues,
			       vector<num_t>& importanceValues);

int main(int argc, char* argv[])
{

  //Print the intro tab
  printHeader();

  //With no input arguments the help is printed
  if(argc == 1)
    {
      printHelp();
      return(EXIT_SUCCESS);
    }

  //Set user parameters to default values
  bool printHelpHandle = DEFAULT_PRINT_HELP;
  string input = "";
  string targetStr = "";
  size_t nTrees = DEFAULT_NTREES;
  size_t mTry = DEFAULT_MTRY;
  size_t nodeSize = DEFAULT_NODESIZE;
  size_t nPerms = DEFAULT_NPERMS;
  num_t pValueThreshold = DEFAULT_PVALUETHRESHOLD;
  bool applyFilter = DEFAULT_APPLY_FILTER;
  string output = "";

  //Read the user parameters 
  ArgParse parser(argc,argv);
  parser.getArgument<bool>("h","help",printHelpHandle);
  parser.getArgument<string>("I","input",input); 
  parser.getArgument<string>("i","target",targetStr); 
  parser.getArgument<size_t>("n","ntrees",nTrees); 
  parser.getArgument<size_t>("m","mtry",mTry); 
  parser.getArgument<size_t>("s","nodesize",nodeSize); 
  parser.getArgument<size_t>("p","nperms",nPerms); 
  parser.getArgument<num_t>("t","pthreshold",pValueThreshold); 
  parser.getArgument<bool>("f","filter",applyFilter); 
  parser.getArgument<string>("O","output",output); 

  if(printHelpHandle)
    {
      printHelp();
      return(EXIT_SUCCESS);
    }

  //Print help and exit if input file is not specified
  if(input == "")
    {
      cerr << "Input file not specified" << endl;
      printHelpHint();
      return EXIT_FAILURE;
    }

  //Print help and exit if target index is not specified
  if(targetStr == "")
    {
      cerr << "target(s) (-i/--target) not specified" << endl;
      printHelpHint();
      return EXIT_FAILURE;
    }

  //Print help and exit if output file is not specified
  if(output == "")
    {
      cerr << "Output file not specified" << endl;
      printHelpHint();
      return EXIT_FAILURE;
    }

  //Read data into Treedata object
  cout << endl << "Reading file '" << input << "'" << endl;
  Treedata treedata(input);

  //Check which feature names match with the specified target identifier
  set<size_t> targetIcs;
  treedata.getMatchingTargetIcs(targetStr,targetIcs);
  if(targetIcs.size() == 0)
    {
      cerr << "No features match the specified target identifier '" << targetStr << "'" << endl;
      return EXIT_FAILURE;
    }

  //The program starts a loop in which an RF-ACE model will be built for each spcified target feature
  size_t iter = 1;
  for(set<size_t>::const_iterator it(targetIcs.begin()); it != targetIcs.end(); ++it, ++iter)
    {
      //Extract the target index from the pointer and the number of real samples from the treedata object
      size_t targetIdx = *it;
      size_t nRealSamples = treedata.nRealSamples(targetIdx);
      
      //If the target has no real samples, the program will just exit
      if(nRealSamples == 0)
	{
	  cout << "Target has no real samples, quitting..." << endl;
	  return EXIT_SUCCESS;
	}
      
      num_t realFraction = 1.0*nRealSamples / treedata.nSamples();
      
      //If default nTrees is to be used...
      if(nTrees == DEFAULT_NTREES)
	{
	  nTrees = static_cast<size_t>( 10.0*treedata.nSamples() / realFraction );
	}
      
      //If default mTry is to be used...
      if(mTry == DEFAULT_MTRY)
	{
	  mTry = static_cast<size_t>( floor( sqrt( 1.0*treedata.nFeatures() ) ) );   
	}
      
      //If default nodeSize is to be used...
      if(nodeSize == DEFAULT_NODESIZE)
	{
	  nodeSize = 5;
	  size_t altNodeSize = static_cast<size_t>( ceil( 1.0*nRealSamples/100 ) );
	  if(altNodeSize > nodeSize)
	    {
	      nodeSize = altNodeSize;
	    }
	}
      
      //Print values of parameters of RF-ACE
      cout << endl;
      cout << "RF-ACE parameter configuration:" << endl;
      cout << "  --input      = " << input << endl;
      cout << "  --nsamples   = " << nRealSamples << " / " << treedata.nSamples() << " (" << 100 * ( 1 - realFraction )<< "% missing)" << endl;
      cout << "  --nfeatures  = " << treedata.nFeatures() << endl;
      cout << "  --targetidx  = " << targetIdx << ", header '" << treedata.getFeatureName(targetIdx) << "'" << endl;
      cout << "  --ntrees     = " << nTrees << endl;
      cout << "  --mtry       = " << mTry << endl;
      cout << "  --nodesize   = " << nodeSize << endl;
      cout << "  --nperms     = " << nPerms << endl;
      cout << "  --pthresold  = " << pValueThreshold << endl;
      cout << "  --output     = " << output << endl << endl;
      
      //Some assertions to check whether the parameter values and data dimensions allow construction of decision trees in the first place
      //NOTE: consider writing a more sophisticated test and provide better print-outs than just assertion failures
      assert(treedata.nFeatures() >= mTry);
      assert(treedata.nSamples() > 2*nodeSize);
      
      //////////////////////////////////////////////////////////////////////
      //  ANALYSIS 1 -- Random Forest ensemble with artificial contrasts  //
      //////////////////////////////////////////////////////////////////////
      
      vector<num_t> pValues; //(treedata.nFeatures());
      vector<num_t> importanceValues; //(treedata.nFeatures());
      
      if(applyFilter)
	{
	  cout << "Applying Random Forest filter..." << endl; 
	  executeRandomForestFilter(treedata,targetIdx,nTrees,mTry,nodeSize,nPerms,pValues,importanceValues);
	   
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
      executeRandomForestFilter(treedata,targetIdx,nTrees,mTry,nodeSize,nPerms,pValues,importanceValues);
    	  
      
      /////////////////////////////////////////////
      //  ANALYSIS 2 -- Gradient Boosting Trees  //
      /////////////////////////////////////////////
      
      if(false)
	{
	  cout << "Gradient Boosting Trees *ENABLED*" << endl;
	  
	  size_t nMaxLeaves = 6;
	  num_t shrinkage = 0.5;
	  num_t subSampleSize = 0.5;
	  GBT gbt(&treedata, targetIdx, nTrees, nMaxLeaves, shrinkage, subSampleSize);
	  
	  
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
	  FILE* po;
	  if(iter == 1)
	    {
	      po = fopen(output.c_str(),"w");
	    }
	  else
	    {
	      po = fopen(output.c_str(),"a");
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
	      
	      if(nPerms > 1)
		{
		  fprintf(po,"%s\t%s\t%9.8f\t%9.8f\t%9.8f\n",targetName.c_str(),treedata.getFeatureName(refIcs[featureIdx]).c_str(),pValues[featureIdx],importanceValues[featureIdx],treedata.pearsonCorrelation(targetIdx,refIcs[featureIdx]));
		}
	      else
		{
		  fprintf(po,"%s\t%s\tNaN\t%9.8f\t%9.8f\n",targetName.c_str(),treedata.getFeatureName(refIcs[featureIdx]).c_str(),importanceValues[featureIdx],treedata.pearsonCorrelation(targetIdx,refIcs[featureIdx]));
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


	  fclose(po);
	  cout << endl << "Association file created. Format:" << endl;
	  cout << "TARGET   PREDICTOR   P-VALUE   IMPORTANCE   CORRELATION" << endl << endl;
	  cout << "Done." << endl;
	}
      else
	{
	  cout << endl << "No significant associations found, quitting..." << endl;
	}
    }
      
      
  return(EXIT_SUCCESS);
}


void executeRandomForestFilter(Treedata& treedata,
                               size_t targetIdx,
                               size_t nTrees,
                               size_t mTry,
                               size_t nodeSize,
			       size_t nPerms,
                               vector<num_t>& pValues,
                               vector<num_t>& importanceValues)
{

  vector<vector<num_t> > importanceMat(nPerms);
  pValues.resize(treedata.nFeatures());
  importanceValues.resize(treedata.nFeatures());
  size_t nNodesInAllForests = 0;

  clock_t time_start(clock());

  cout << "Growing " << nPerms << " Random Forests (RFs), please wait..." << endl;
  //#pragma omp parallel for
  for(int permIdx = 0; permIdx < static_cast<int>(nPerms); ++permIdx)
    {
      //cout << "  RF " << permIdx + 1 << ": ";
      //Treedata td_thread = treedata;
      bool useContrasts = false;
      Randomforest RF(&treedata,targetIdx,nTrees,mTry,nodeSize,useContrasts);
      size_t nNodesInForest = RF.nNodes();
      nNodesInAllForests += nNodesInForest;
      importanceMat[permIdx] = RF.featureImportance();
      printf("  RF %i: %i nodes (avg. %6.3f nodes/tree)\n",permIdx+1,static_cast<int>(nNodesInForest),1.0*nNodesInForest/nTrees);
    }

  num_t time_diff = 1.0*(clock() - time_start) / CLOCKS_PER_SEC;
  cout << nPerms << " RFs, " << nPerms*nTrees << " trees, and " << nNodesInAllForests
       << " nodes generated in " << time_diff << " seconds (" << 1.0*nNodesInAllForests / time_diff
       << " nodes per second)" << endl;
  
  if(nPerms > 1)
    {
      for(size_t featureIdx = 0; featureIdx < treedata.nFeatures(); ++featureIdx)
	{
	  
	  size_t nRealSamples;
	  vector<num_t> fSample(nPerms);
	  vector<num_t> cSample(nPerms);
	  for(size_t permIdx = 0; permIdx < nPerms; ++permIdx)
	    {
	      fSample[permIdx] = importanceMat[permIdx][featureIdx];
	      cSample[permIdx] = importanceMat[permIdx][featureIdx + treedata.nFeatures()];
	    }
	  datadefs::utest(fSample,cSample,pValues[featureIdx]);
	  datadefs::mean(fSample,importanceValues[featureIdx],nRealSamples);
	}
    }
  else
    {
      importanceValues = importanceMat[0];
    }

  importanceValues.resize(treedata.nFeatures());
  
}
