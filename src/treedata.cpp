#include "treedata.hpp"
#include<cstdlib>
#include<fstream>
#include<cassert>
#include<iostream>
#include<sstream>
#include<utility>
#include<algorithm>
#include<ctime>

using namespace std;

Treedata::Treedata(string fname, bool is_featurerows):
  targetidx_(0),
  featurematrix_(0),
  isfeaturenum_(0),
  nsamples_(0),
  nfeatures_(0),
  sampleics_(0),
  featureics_(0),
  featureheaders_(0)
{

  //Initialize random number rgenerator
  time_t now;
  time(&now);
  srand((unsigned int)now);

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
  vector<vector<string> > rawmatrix(nrows);
  for(size_t i = 0; i < nrows; ++i)
    {
      rawmatrix[i] = colheaders;
    }

  cout << "read " << rawmatrix.size() << " rows and " << rawmatrix[0].size() << " columns." << endl;

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
	  getline(ss,rawmatrix[i][j],'\t');
	  cout << '\t' << rawmatrix[i][j];
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
      featureheaders_ = rowheaders;
      for(size_t i = 0; i < nfeatures_; ++i)
	{
	  //First letters in the row headers determine whether the feature is numerical or categorical
	  if(rowheaders[i][0] == 'N')
	    {
	      isfeaturenum_.push_back(true);
	    }
	  else if(rowheaders[i][0] == 'C')
	    {
	      isfeaturenum_.push_back(false);
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

  //Store indices 0,1,...,(nsamples_-1) to sampleics_ 
  vector<size_t> sampleics(nsamples_);
  Treedata::range(sampleics);
  sampleics_ = sampleics;

  //Store indices 0,1,...,(nfeatures_-1) to featureics_;
  vector<size_t> featureics(nfeatures_);
  Treedata::range(featureics);
  featureics_ = featureics;

  //Transform raw data to the internal format.
  for(size_t i = 0; i < nfeatures_; ++i)
    {
      vector<num_t> featurev(nsamples_);
      map<string,size_t> str2valmap;
      if(isfeaturenum_[i])
	{
	  datadefs::strv2numv(rawmatrix[i],featurev);
	  //featurematrix_.push_back(featurev);
	}
      else
	{
	  datadefs::strv2catv(rawmatrix[i],featurev,str2valmap);
	}
      featurematrix_.push_back(featurev);
    } 

  //By default, make the first feature the target
  targetidx_ = 0;
  Treedata::select_target(targetidx_);
  Treedata::sort_all_wrt_target();
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
  cout << "Printing feature matrix (missing values encoded to " << datadefs::num_nan << "):" << endl;
  for(size_t j = 0; j < nsamples_; ++j)
    {
      cout << '\t' << sampleheaders_[j];
    }
  cout << endl;
  for(size_t i = 0; i < nfeatures_; ++i)
    {
      cout << i << ':' << featureheaders_[i] << ':';
      for(size_t j = 0; j < nsamples_; ++j)
        {
          cout << '\t' << featurematrix_[i][j];
        }
      cout << endl;
    }
}

template <typename T1,typename T2> 
void Treedata::make_pairedv(vector<T1> const& v1, vector<T2> const& v2, vector<pair<T1,T2> >& p)
{
  assert(v1.size() == v2.size() && v2.size() == p.size() && p.size() == nsamples_);
  for(size_t i = 0; i < nsamples_; ++i)
    {
      p[i] = make_pair(v1[i],v2[i]);
    }
}

template <typename T1,typename T2> 
void Treedata::separate_pairedv(vector<pair<T1,T2> > const& p, vector<T1>& v1, vector<T2>& v2)
{
  assert(v1.size() == v2.size() && v2.size() == p.size() && p.size() == nsamples_);
  for(size_t i = 0; i < nsamples_; ++i)
    {
      v1[i] = p[i].first;
      v2[i] = p[i].second;
    }
}

template <typename T> void Treedata::sort_and_make_ref(vector<T>& v, vector<size_t>& ref_ics)
{
  
  //Make a paired vector of v and ref_ics
  vector<pair<T,size_t> > pairedv(ref_ics.size());
  ref_ics = sampleics_;
  Treedata::make_pairedv<T,size_t>(v,ref_ics,pairedv);
  
  //Sort the paired vector with respect to the first element
  sort(pairedv.begin(),pairedv.end(),datadefs::ordering<T>());
  
  //Separate the paired vector
  Treedata::separate_pairedv<T,size_t>(pairedv,v,ref_ics);
}


template <typename T>
void Treedata::sort_from_ref(vector<T>& in, vector<size_t> const& ref_ics)
{
  vector<T> foo = in;
  
  int n = in.size();
  for (int i = 0; i < n; ++i)
    {
      in[i] = foo[ref_ics[i]];
    }
}

void Treedata::select_target(size_t targetidx)
{
  targetidx_ = targetidx; 
  Treedata::sort_all_wrt_target();
}

void Treedata::sort_all_wrt_feature(size_t featureidx)
{
  //Generate and index vector that'll define the new order
  vector<size_t> ref_ics(nsamples_);
  
  //Make a dummy variable so that the original feature matrix will stay intact 
  vector<num_t> dummy = featurematrix_[featureidx];

  //Sort specified feature and make reference indices
  Treedata::sort_and_make_ref<num_t>(dummy,ref_ics);
  
  //Use the new order to sort the sample headers
  Treedata::sort_from_ref<string>(sampleheaders_,ref_ics);

  //Sort the other features as well
  for(size_t i = 0; i < nfeatures_; ++i)
    {
      Treedata::sort_from_ref<num_t>(featurematrix_[i],ref_ics);
    }
}

void Treedata::sort_all_wrt_target()
{
  Treedata::sort_all_wrt_feature(targetidx_);
}

void Treedata::range(vector<size_t>& ics)
{
  for(size_t i = 0; i < ics.size(); ++i)
    {
      ics[i] = i;
    }
}


void Treedata::permute(vector<size_t>& ics)
{
  size_t n = ics.size();
  for (size_t i = 0; i < n; ++i)
    {
      size_t j = rand() % (i + 1);
      ics[i] = ics[j];
      ics[j] = i;
    }
}

void Treedata::permute(vector<num_t>& data)
{
  size_t n = data.size();
  vector<size_t> ics(n);
  vector<num_t> foo = data;
  Treedata::permute(ics);
  for(size_t i = 0; i < n; ++i)
    {
      data[i] = foo[ics[i]];
    }
}

void Treedata::bootstrap(vector<size_t>& ics, vector<size_t>& oob_ics, size_t& noob)
{

  for(size_t i = 0; i < nsamples_; ++i)
    {
      ics[i] = rand() % nsamples_;
    }
  sort(ics.begin(),ics.end());
  vector<size_t> allics(nsamples_);
  Treedata::range(allics);
  vector<size_t>::iterator it = set_difference(allics.begin(),allics.end(),ics.begin(),ics.end(),oob_ics.begin());
  noob = distance(oob_ics.begin(),it);
}

//void Treedata::find_split(size_t featureidx,
//                         size_t& split_pos,
//                         num_t& impurity_left,
//                         num_t& impurity_right)
//{
// Treedata::find_split(featureidx,sampleics_,split_pos,impurity_left,impurity_right);
//}

void Treedata::find_split(size_t featureidx,
			  vector<size_t>& sampleics,
			  vector<size_t>& sampleics_left,
			  vector<size_t>& sampleics_right,
			  size_t& n_left,
			  size_t& n_right,
			  num_t& impurity,
			  num_t& impurity_left,
			  num_t& impurity_right)
{


  size_t n(sampleics.size());
  size_t n_nonnan(0);
  num_t mu,se;

  
  datadefs::sqerr(featurematrix_[featureidx],sampleics,n_nonnan,mu,se);
  cout << "feature idx " << featureidx << " squared error: mu = " << mu << "\tsq.err. = " << se << "\t " << n-n_nonnan << " missing values" << endl;

  vector<num_t> fv(n);
  vector<num_t> tv(n);
  

}

void Treedata::split_at_pos(size_t featureidx,
			    vector<size_t>& sampleics,
			    num_t& impurity_left,
			    num_t& impurity_right)
{

}






