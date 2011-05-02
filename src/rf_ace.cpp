#include <cstdlib>
#include <cassert>
#include <string>
#include <iostream>
#include <fstream>
#include <ctime>
#include <cmath>


//#include <boost/algorithm/string.hpp>
//#include <boost/program_options.hpp>

#include "getopt_pp.h"
#include "randomforest.hpp"
//#include "GBT.hpp"
#include "treedata.hpp"
#include "datadefs.hpp"

using namespace std;
using datadefs::num_t;
using namespace GetOpt;

//namespace po = boost::program_options;

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

  if(argc == 1)
    {
      cout << "REQUIRED:" << endl;
      cout << "-I / --input        input feature matrix" << endl;
      cout << "-O / --output       output association file" << endl << endl;
      cout << "OPTIONAL:" << endl;
      cout << "-i / --targetidx    target index, ref. to feature matrix (default " << DEFAULT_TARGETIDX << ")" << endl;
      cout << "-n / --ntrees       number of trees per RF (default nsamples)" << endl;
      cout << "-m / --mtry         number of randomly drawn features per node split (default sqrt(nfeatures))" << endl;
      cout << "-s / --nodesize     minimum number of train samples per node, affects tree depth (default " << DEFAULT_NODESIZE << ")" << endl;
      cout << "-p / --nperms       number of Random Forests (default " << DEFAULT_NPERMS << ")" << endl;
      cout << "-t / --pthreshold   p-value threshold below which associations are listed (default " << DEFAULT_PTHRESHOLD << ")" << endl;
      cout << "-a / --alpha        percentile of contrast importances, defines stringency of the t-test (default " << DEFAULT_ALPHA << ")" << endl;
      return EXIT_SUCCESS;
    }

  //using namespace GetOpt;
  string input = "";
  size_t targetidx = DEFAULT_TARGETIDX;
  size_t ntrees = DEFAULT_NTREES;
  size_t mtry = DEFAULT_MTRY;
  size_t nodesize = DEFAULT_NODESIZE;
  size_t nperms = DEFAULT_NPERMS;
  num_t pthreshold = DEFAULT_PTHRESHOLD;
  num_t alpha = DEFAULT_ALPHA;
  bool is_featurerows = true;
  string output = "";

  GetOpt_pp ops(argc, argv);
        
  ops >> Option('I',"input",input);
  ops >> Option('i',"targetidx",targetidx);
  ops >> Option('n',"ntrees",ntrees);
  ops >> Option('m',"mtry",mtry);
  ops >> Option('p',"nperms",nperms);
  ops >> Option('O',"output",output);

  //Read data into Treedata object
  Treedata treedata(input, is_featurerows);
  
  if(ntrees == DEFAULT_NTREES)
    {
      ntrees = treedata.nsamples();
    }
  
  if(mtry == DEFAULT_MTRY)
    {
      mtry = static_cast<size_t>(floor(sqrt(static_cast<num_t>(treedata.nfeatures()))));   
    }

  //treedata.print(targetidx);
  
  cout << "Selected options:" << endl;
  cout << "--targetidx = " << targetidx << " (" << treedata.get_featureheader(targetidx) << ")" << endl;
  cout << "--input     = " << input << endl;
  cout << "--ntrees    = " << ntrees << endl;
  cout << "--mtry      = " << mtry << endl;
  cout << "--nodesize  = " << nodesize << endl;
  cout << "--nperms    = " << nperms << endl;
  cout << "--pthresold = " << pthreshold << endl;
  cout << "--alpha     = " << alpha << endl;
  cout << "--output    = " << output << endl;

  assert(treedata.nfeatures() >= mtry);
  assert(treedata.nsamples() > 2*nodesize);

  Randomforest RF(&treedata,ntrees,mtry,nodesize);
  RF.select_target(targetidx);

  //size_t nperms = 9;
  //num_t alpha = 0.50;
  vector<num_t> pvalues(treedata.nfeatures());
  vector<num_t> ivalues(treedata.nfeatures());
  //vector<num_t> ivalues(treedata.nfeatures());
  
  //clock_t time_start(clock());
  RF.grow_forest(nperms,alpha,pvalues,ivalues);
  //cout << "Time elapsed: " << float(clock() - time_start)/CLOCKS_PER_SEC << " seconds" << endl;
  
  vector<size_t> ref_ics(treedata.nfeatures());
  //vector<string> fnames = treedata.featureheaders();
  datadefs::sort_and_make_ref<num_t>(pvalues,ref_ics);
  datadefs::sort_from_ref<num_t>(ivalues,ref_ics);
  //datadefs::sort_from_ref<string>(fnames,ref_ics);
  
  string target_str = treedata.get_featureheader(targetidx);
  
  ofstream os(output.c_str());
  for(size_t i = 0; i < treedata.nfeatures(); ++i)
    {
      if(pvalues[i] > pthreshold)
	{
	  break;
	}
      os << target_str << "\t" << treedata.get_featureheader(ref_ics[i]) << "\t" 
	 << pvalues[i] << "\t" << ivalues[i] << "\t" << treedata.corr(targetidx,ref_ics[i]) << endl;
    }
  os.close();
 
  return(EXIT_SUCCESS);
}
