#include <cstdlib>
#include <cassert>
#include <string>
#include <iostream>
#include <fstream>
#include <ctime>
#include <cmath>
#include <stdio.h>

#include "getopt_pp.h"
#include "randomforest.hpp"
//#include "GBT.hpp"
#include "treedata.hpp"
#include "datadefs.hpp"

using namespace std;
using datadefs::num_t;
using namespace GetOpt;

enum MatrixFormat { FEATURE_ROWS, FEATURE_COLUMNS };

const size_t DEFAULT_TARGETIDX = 0;
const size_t DEFAULT_NTREES = 0;
const size_t DEFAULT_MTRY = 0;
const size_t DEFAULT_NODESIZE = 5;
const size_t DEFAULT_NPERMS = 20;
const num_t DEFAULT_PTHRESHOLD = 0.10;
const num_t DEFAULT_ALPHA = 0.95;

int main(int argc, char* argv[])
{

  if(argc == 1 || argc == 2)
    {

      if(argc == 2)
	{
	  string helphandle(argv[1]);
	  if(helphandle != "-h" && helphandle != "--help")
	    {
	      cerr << "use -h or --help to get started" << endl;
	      return EXIT_FAILURE;
	    }
	}

      cout << endl;
      cout << "REQUIRED ARGUMENTS:" << endl;
      cout << "-I / --input        input feature matrix" << endl;
      cout << "-O / --output       output association file" << endl;
      cout << endl;
      cout << "OPTIONAL ARGUMENTS:" << endl;
      cout << "-i / --targetidx    target index, ref. to feature matrix (default " << DEFAULT_TARGETIDX << ")" << endl;
      cout << "-n / --ntrees       number of trees per RF (default nsamples)" << endl;
      cout << "-m / --mtry         number of randomly drawn features per node split (default sqrt(nfeatures))" << endl;
      cout << "-s / --nodesize     minimum number of train samples per node, affects tree depth (default " << DEFAULT_NODESIZE << ")" << endl;
      cout << "-p / --nperms       number of Random Forests (default " << DEFAULT_NPERMS << ")" << endl;
      cout << "-t / --pthreshold   p-value threshold below which associations are listed (default " << DEFAULT_PTHRESHOLD << ")" << endl;
      cout << "-a / --alpha        percentile of contrast importances, defines stringency of the t-test (default " << DEFAULT_ALPHA << ")" << endl;
      cout << endl;
      return EXIT_SUCCESS;
    }

  cout << endl;
  cout << "  ----------------------------------" << endl;
  cout << "  ---    RF-ACE version 0.3.0    ---" << endl;
  cout << "  ----------------------------------" << endl;

  //using namespace GetOpt;
  string input = "";
  size_t targetidx = DEFAULT_TARGETIDX;
  size_t ntrees = DEFAULT_NTREES;
  size_t mtry = DEFAULT_MTRY;
  size_t nodesize = DEFAULT_NODESIZE;
  size_t nperms = DEFAULT_NPERMS;
  num_t pthreshold = DEFAULT_PTHRESHOLD;
  num_t alpha = DEFAULT_ALPHA;
  //bool is_featurerows = true;
  string output = "";

  GetOpt_pp ops(argc, argv);
        
  ops >> Option('I',"input",input);
  ops >> Option('i',"targetidx",targetidx);
  ops >> Option('n',"ntrees",ntrees);
  ops >> Option('m',"mtry",mtry);
  ops >> Option('p',"nperms",nperms);
  ops >> Option('t',"pthreshold",pthreshold);
  ops >> Option('O',"output",output);

  if(input == "")
    {
      cerr << "input file not specified" << endl;
      return EXIT_FAILURE;
    }

  if(output == "")
    {
      cerr << "output file not specified" << endl;
      return EXIT_FAILURE;
    }

  //Read data into Treedata object
  cout << endl << "Reading file '" << input << "'" << endl;
  Treedata treedata(input);

 
  if(targetidx >= treedata.nfeatures())
    {
      cerr << "targetidx must point to a valid feature (0.." << treedata.nfeatures()-1 << ")" << endl;
      return EXIT_FAILURE;
    }

  size_t nrealsamples = treedata.nrealsamples(targetidx);

  if(nrealsamples == 0)
    {
      cout << "Target has no real samples, quitting..." << endl;
      return EXIT_SUCCESS;
    }

  num_t realfraction = static_cast<num_t>(nrealsamples) / static_cast<num_t>(treedata.nsamples());

  if(ntrees == DEFAULT_NTREES)
    {
      ntrees = static_cast<size_t>(2.0 * static_cast<num_t>(treedata.nsamples()) / realfraction );
    }
  
  if(mtry == DEFAULT_MTRY)
    {
      mtry = static_cast<size_t>(floor(sqrt(static_cast<num_t>(treedata.nfeatures()))));   
    }

  //treedata.print(targetidx);
  
  cout << endl;
  cout << "RF-ACE parameter configuration:" << endl;
  cout << "  --input      = " << input << endl;
  cout << "  --nsamples   = " << nrealsamples << " / " << treedata.nsamples() << " (" << 100 * ( 1 - realfraction )<< "% missing)" << endl;
  cout << "  --nfeatures  = " << treedata.nfeatures() << endl;
  cout << "  --targetidx  = " << targetidx << ", header '" << treedata.get_featureheader(targetidx) << "'" << endl;
  cout << "  --ntrees     = " << ntrees << endl;
  cout << "  --mtry       = " << mtry << endl;
  cout << "  --nodesize   = " << nodesize << endl;
  cout << "  --nperms     = " << nperms << endl;
  cout << "  --pthresold  = " << pthreshold << endl;
  cout << "  --alpha      = " << alpha << endl;
  cout << "  --output     = " << output << endl << endl;

  assert(treedata.nfeatures() >= mtry);
  assert(treedata.nsamples() > 2*nodesize);

  Randomforest RF(&treedata,ntrees,mtry,nodesize);
  RF.select_target(targetidx);
  //RF.blacklist_and_kill(0.8,blistheaders,blistcorrelations);

  //size_t nperms = 9;
  //num_t alpha = 0.50;
  vector<num_t> pvalues(treedata.nfeatures());
  vector<num_t> ivalues(treedata.nfeatures());
  //vector<num_t> ivalues(treedata.nfeatures());
  
  cout << "Growing " << nperms << " Random Forests (RFs), please wait..." << endl;
  RF.grow_forest(nperms,alpha,pvalues,ivalues);
  //cout << "Time elapsed: " << float(clock() - time_start)/CLOCKS_PER_SEC << " seconds" << endl;
  
  vector<size_t> ref_ics(treedata.nfeatures());
  //vector<string> fnames = treedata.featureheaders();
  datadefs::sort_and_make_ref<num_t>(pvalues,ref_ics);
  datadefs::sort_from_ref<num_t>(ivalues,ref_ics);
  //datadefs::sort_from_ref<string>(fnames,ref_ics);
  
  string target_str = treedata.get_featureheader(targetidx);

  if(pvalues[0] <= pthreshold)
    {
      //cout << "Writing associations to file '" << output << "'" << endl;
      //ofstream os(output.c_str());
      FILE* po;
      po = fopen(output.c_str(),"w");
      
      for(size_t i = 0; i < treedata.nfeatures(); ++i)
	{
	  if(pvalues[i] > pthreshold)
	    {
	      break;
	    }

	  fprintf(po,"%s\t%s\t%9.8f\t%9.8f\t%9.8f\n",target_str.c_str(),treedata.get_featureheader(ref_ics[i]).c_str(),pvalues[i],ivalues[i],treedata.corr(targetidx,ref_ics[i]));

	  //os << target_str << "\t" << treedata.get_featureheader(ref_ics[i]) << "\t" 
	  //   << pvalues[i] << "\t" << ivalues[i] << "\t" << treedata.corr(targetidx,ref_ics[i]) << endl;
	}
      //os.close();
      fclose(po);
      cout << endl << "Association file created. Format:" << endl;
      cout << "TARGET   PREDICTOR   P-VALUE   IMPORTANCE   CORRELATION" << endl << endl;
      cout << "Done." << endl;
    }
  else
    {
      cout << endl << "No significant associations found, quitting..." << endl;
    }

  
 
  return(EXIT_SUCCESS);
}
