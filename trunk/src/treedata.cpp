#include "treedata.hpp"
#include<cstdlib>
#include<fstream>
#include<cassert>
#include<iostream>
#include<sstream>
#include<utility>
#include<algorithm>

using namespace std;

Treedata::Treedata(string fname, bool is_featurerows):
  istarget_(false),
  targetidx_(0),
  internaltargetidx_(0),
  isnumtarget_(false),
  catmatrix_(0),
  nummatrix_(0),
  nsamples_(0),
  nfeatures_(0),
  ncatfeatures_(0),
  nnumfeatures_(0),
  catfeatureheaders_(0),
  numfeatureheaders_(0),
  catfeatureics_(0),
  numfeatureics_(0)
{

  cout << "Treedata: reading matrix from file '" << fname << "'" << endl;

  //Initialize stream to read from file
  ifstream featurestream;
  featurestream.open(fname.c_str());
  assert(featurestream.good());

  string field;
  string row;

  //Remove upper left element from the matrix as useless
  getline(featurestream,field,'\t');

  //Count the number of columns
  getline(featurestream,row);
  stringstream ss(row);
  size_t ncols = 0;
  while(getline(ss,field,'\t'))
    {
      ++ncols;
    }

  //Count the number of rows
  size_t nrows = 0;
  while(getline(featurestream,row))
    {
      ++nrows;
    }

  //Reset streams and remove upper left element from the matrix as useless
  featurestream.clear();
  featurestream.seekg(0, ios::beg);
  getline(featurestream,field,'\t');
  ss.clear();
  ss.str("");

  //These are temporary containers
  vector<string> colheaders(ncols);
  vector<string> rowheaders(nrows);
  vector<vector<string> > datamatrix(nrows);
  for(size_t i = 0; i < nrows; ++i)
    {
      datamatrix[i] = colheaders;
    }

  cout << "read " << datamatrix.size() << " rows and " << datamatrix[0].size() << " columns." << endl;

  //Read first row into the stream
  getline(featurestream,row);
  ss << row;

  //Read column headers from the stream
  for(size_t i = 0; i < ncols; ++i)
    {
      getline(ss,colheaders[i],'\t');
      cout << '\t' << colheaders[i];
    }
  cout << endl;

  //Go through the rest of the rows
  for(size_t i = 0; i < nrows; ++i)
    {
      //Read row from the stream
      getline(featurestream,row);
      ss.clear();
      ss.str("");

      //Read the string back to a stream (REDUNDANCY)
      ss << row;

      //Read one element from the row stream
      getline(ss,rowheaders[i],'\t');
      cout << rowheaders[i];
      for(size_t j = 0; j < ncols; ++j)
	{
	  getline(ss,datamatrix[i][j],'\t');
	  cout << '\t' << datamatrix[i][j];
	}
      cout << endl;
    }
  cout << endl;
  
  //If the data is row-formatted...
  if(is_featurerows)
    {
      //We can extract the number of features and samples from the row and column counts, respectively
      nfeatures_ = nrows;
      nsamples_ = ncols;

      //Thus, sample headers are column headers
      sampleheaders_ = colheaders;
      for(size_t i = 0; i < nfeatures_; ++i)
	{
	  //First letters in the row headers determine whether the feature is numerical or categorical
	  if(rowheaders[i][0] == 'N')
	    {
	      numfeatureheaders_.push_back(rowheaders[i]);
	      numfeatureics_.push_back(i);
	      ++nnumfeatures_;
	    }
	  else if(rowheaders[i][0] == 'C')
	    {
	      catfeatureheaders_.push_back(rowheaders[i]);
	      catfeatureics_.push_back(i);
	      ++ncatfeatures_;
	    }
	  else
	    {
	      cerr << "Data type must be either N or C!" << endl;
	      assert(false);
	    }
	}
    }
  else
    {
      cerr << "samples as rows not yet supported!" << endl;
      assert(false);
    }

  //Transform raw data to the internal format. First categorical data...
  for(size_t i = 0; i < ncatfeatures_; ++i)
    {
      vector<cat_t> foo(nsamples_);
      datadefs::strv2catv(datamatrix[catfeatureics_[i]],foo);
      catmatrix_.push_back(foo);
    }
  
  
  //... and next numerical data.
  for(size_t i = 0; i < nnumfeatures_; ++i)
    {
      vector<num_t> foo(nsamples_);
      datadefs::strv2numv(datamatrix[numfeatureics_[i]],foo);
      nummatrix_.push_back(foo);
    }
}

Treedata::~Treedata()
{
}

size_t Treedata::nfeatures()
{
  return(nfeatures_);
}

size_t Treedata::nsamples()
{
  return(nsamples_);
}

void Treedata::print()
{
  cout << "Printing categorical data (missing values encoded to " << datadefs::cat_nan << "):" << endl;
  for(size_t j = 0; j < nsamples_; ++j)
    {
      cout << '\t' << sampleheaders_[j];
    }
  cout << endl;
  for(size_t i = 0; i < ncatfeatures_; ++i)
    {
      cout << catfeatureics_[i] << ':' << catfeatureheaders_[i] << ':';
      for(size_t j = 0; j < nsamples_; ++j)
        {
          cout << '\t' << catmatrix_[i][j];
        }
      cout << endl;
    }
  cout << "Printing numerical data (missing values encoded to " << datadefs::num_nan << "):" << endl;
  for(size_t j = 0; j < nsamples_; ++j)
    {
      cout << '\t' << sampleheaders_[j];
    }
  cout << endl;
  for(size_t i = 0; i < nnumfeatures_; ++i)
    {
      cout << numfeatureics_[i] << ':' << numfeatureheaders_[i] << ':';
      for(size_t j = 0; j < nsamples_; ++j)
        {
          cout << '\t' << nummatrix_[i][j];
        }
      cout << endl;
    }
}

void Treedata::range(vector<size_t>& ics)
{
  assert(nsamples_ == ics.size());
  for(size_t i = 0; i < nsamples_; ++i)
    {
      ics[i] = i;
    }
}

template <typename T1,typename T2> 
void Treedata::join_pairedv(vector<T1>& v1, vector<T2>& v2, vector<pair<T1,T2> >& p)
{
  assert(v1.size() == v2.size() && v2.size() == p.size() && p.size() == nsamples_);
  for(size_t i = 0; i < nsamples_; ++i)
    {
      p[i] = make_pair(v1[i],v2[i]);
    }
}

template <typename T1,typename T2> 
void Treedata::separate_pairedv(vector<pair<T1,T2> >& p, vector<T1>& v1, vector<T2>& v2)
{
  assert(v1.size() == v2.size() && v2.size() == p.size() && p.size() == nsamples_);
  for(size_t i = 0; i < nsamples_; ++i)
    {
      v1[i] = p[i].first;
      v2[i] = p[i].second;
    }
}

template <typename T>
void Treedata::sort_from_ref(vector<T>& in, vector<size_t> const& reference)
{
  vector<T> foo = in;
  
  int n = in.size();
  for (int i = 0; i < n; ++i)
    {
      in[i] = foo[reference[i]];
    }
}

void Treedata::select_target(size_t targetidx)
{
  targetidx_ = targetidx; 
  istarget_ = true;
  bool isfound(false);
  for(size_t i = 0; i < ncatfeatures_; ++i)
    {
      if(catfeatureics_[i] == targetidx_)
	{
	  isnumtarget_ = false;
	  internaltargetidx_ = i;
	  isfound = true;
	}
    }
  if(!isfound)
    {
      for(size_t i = 0; i < nnumfeatures_; ++i)
	{
	  if(numfeatureics_[i] == targetidx_)
	    {
	      isnumtarget_ = true;
	      internaltargetidx_ = i;
	      isfound = true;
	    } 
	}
    }
  assert(isfound);

  Treedata::sort_all_wrt_target();
}

void Treedata::sort_all_wrt_feature(size_t featureidx)
{
  assert(istarget_);
}

void Treedata::sort_all_wrt_target()
{
  //Check that a target has been set
  assert(istarget_);

  //Generate and index vector that'll define the new order
  vector<size_t> neworder_ics(nsamples_);

  if(isnumtarget_)
    {
      //Generate a paired vector with which sorting will be performed
      vector<pair<num_t,size_t> > pairedv(nsamples_);
      
      //Generate indices from 0,1,...,(nsamples-1)
      Treedata::range(neworder_ics);
      
      //Join the target and index vector
      Treedata::join_pairedv<num_t,size_t>(nummatrix_[internaltargetidx_],neworder_ics,pairedv);
      
      //Sort the paired vector (indices will now define the new order)
      sort(pairedv.begin(),pairedv.end(),datadefs::ordering<size_t>());
      
      vector<num_t> foo(nsamples_);
      //Separate the target vector and new order
      Treedata::separate_pairedv<num_t,size_t>(pairedv,foo,neworder_ics);
    }
  else
    {
      //Generate a paired vector with which sorting will be performed
      vector<pair<cat_t,size_t> > pairedv(nsamples_);

      //Generate indices from 0,1,...,(nsamples-1)
      Treedata::range(neworder_ics);

      //Join the target and index vector
      Treedata::join_pairedv<cat_t,size_t>(catmatrix_[internaltargetidx_],neworder_ics,pairedv);

      //Sort the paired vector (indices will now define the new order)
      sort(pairedv.begin(),pairedv.end(),datadefs::ordering<size_t>());

      vector<cat_t> foo(nsamples_);
      //Separate the target vector and new order
      Treedata::separate_pairedv<cat_t,size_t>(pairedv,foo,neworder_ics);
    }      

  //Use the new order to sort the sample headers
  Treedata::sort_from_ref<string>(sampleheaders_,neworder_ics);
  
  //Sort categorical and numerical data as well
  for(size_t i = 0; i < ncatfeatures_; ++i)
    {
      Treedata::sort_from_ref<cat_t>(catmatrix_[i],neworder_ics);
    }
  for(size_t i = 0; i < nnumfeatures_; ++i)
    {
      Treedata::sort_from_ref<num_t>(nummatrix_[i],neworder_ics);
    }

}


void Treedata::transpose()
{
  
}



