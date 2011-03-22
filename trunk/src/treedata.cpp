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
  nrealvalues_(0),
  ncatvalues_(0),
  nsamples_(0),
  nfeatures_(0),
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
      //cout << '\t' << colheaders[i];
    }
  //cout << endl;

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
      //cout << rowheaders[i];
      for(size_t j = 0; j < ncols; ++j)
        {
          getline(ss,rawmatrix[i][j],'\t');
          //cout << '\t' << rawmatrix[i][j];
        }
      //cout << endl;
    }
  //cout << endl;
  
  //If the data is row-formatted...
  if(is_featurerows)
    {
      //We can extract the number of features and samples from the row and column counts, respectively
      nfeatures_ = nrows; 
      nsamples_ = ncols;

      //Thus, sample headers are column headers
      sampleheaders_ = colheaders;

      //... and feature headers are row headers
      featureheaders_ = rowheaders;

      for(size_t i = 0; i < nfeatures_; ++i)
        {
          //First letters in the row headers determine whether the feature is numerical or categorical
          if(rowheaders[i][0] == 'N')
            {
              isfeaturenum_.push_back(true);
            }
          else if(rowheaders[i][0] == 'C' or rowheaders[i][0] == 'B')
            {
              isfeaturenum_.push_back(false);
            }
          else
            {
              cerr << "Data type must be either N (numerical), C (categorical), or B (binary)!" << endl;
              assert(false);
            }
        }
    }
  else
    {
      cerr << "samples as rows not yet supported!" << endl;
      assert(false);
    }

  vector<size_t> ncatvalues(nfeatures_);
  ncatvalues_ = ncatvalues;

  //Transform raw data to the internal format.
  for(size_t i = 0; i < nfeatures_; ++i)
    {
      vector<num_t> featurev(nsamples_);
      map<string,size_t> str2valmap;
      if(isfeaturenum_[i])
        {
          datadefs::strv2numv(rawmatrix[i],featurev);
          ncatvalues_[i] = 0;
        }
      else
        {
          datadefs::strv2catv(rawmatrix[i],featurev,str2valmap);
	  ncatvalues_[i] = str2valmap.size();
        }
      featurematrix_.push_back(featurev);
    } 
  
  vector<size_t> nrealvalues(nfeatures_);
  nrealvalues_ = nrealvalues;
  
  for(size_t featureidx = 0; featureidx < nfeatures_; ++featureidx)
    {
      Treedata::count_real_values(featureidx,nrealvalues_[featureidx]);
    }
  
  //By default, make the first feature the target
  targetidx_ = 0;
  Treedata::select_target(targetidx_); //This will also sort the data wrt targetidx_
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

void Treedata::count_real_values(size_t featureidx, size_t& nreal)
{
  nreal = 0;
  for(size_t i = 0; i < nsamples_; ++i)
    {
      if(!datadefs::is_nan(featurematrix_[featureidx][i]))
        {
          ++nreal;
        }
    }
}

template <typename T1,typename T2> 
void Treedata::make_pairedv(vector<T1> const& v1, vector<T2> const& v2, vector<pair<T1,T2> >& p)
{
  assert(v1.size() == v2.size() && v2.size() == p.size());
  for(size_t i = 0; i < p.size(); ++i)
    {
      p[i] = make_pair(v1[i],v2[i]);
    }
}

template <typename T1,typename T2> 
void Treedata::separate_pairedv(vector<pair<T1,T2> > const& p, vector<T1>& v1, vector<T2>& v2)
{
  assert(v1.size() == v2.size() && v2.size() == p.size());
  for(size_t i = 0; i < p.size(); ++i)
    {
      v1[i] = p[i].first;
      v2[i] = p[i].second;
    }
}

template <typename T> void Treedata::sort_and_make_ref(vector<T>& v, vector<size_t>& ref_ics)
{
  assert(v.size() == ref_ics.size());
  
  //Make a paired vector of v and ref_ics
  vector<pair<T,size_t> > pairedv(ref_ics.size());
  
  Treedata::range(ref_ics);
  Treedata::make_pairedv<T,size_t>(v,ref_ics,pairedv);
  
  //Sort the paired vector with respect to the first element
  sort(pairedv.begin(),pairedv.end(),datadefs::ordering<T>());
  
  //Separate the paired vector
  Treedata::separate_pairedv<T,size_t>(pairedv,v,ref_ics);
}


template <typename T>
void Treedata::sort_from_ref(vector<T>& in, vector<size_t> const& ref_ics)
{
  assert(in.size() == ref_ics.size());
  
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

size_t Treedata::get_target()
{
  return(targetidx_);
}

size_t Treedata::nrealvalues()
{
  return(nrealvalues_[targetidx_]);
}

size_t Treedata::nrealvalues(size_t featureidx)
{
  return(nrealvalues_[featureidx]);
}

void Treedata::sort_all_wrt_feature(size_t featureidx)
{
  assert(featureidx < nfeatures_);

  //Generate and index vector that'll define the new order
  vector<size_t> ref_ics(nsamples_);
  
  //Make a dummy variable so that the original feature matrix will stay intact 
  vector<num_t> dummy(nsamples_);

  dummy = featurematrix_[featureidx];

  //Sort specified feature and make reference indices
  Treedata::sort_and_make_ref<num_t>(dummy,ref_ics);
  
  //Use the new order to sort the sample headers
  Treedata::sort_from_ref<string>(sampleheaders_,ref_ics);

  //Sort features
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
  for (size_t i = 0; i < ics.size(); ++i)
    {
      size_t j = rand() % (i + 1);
      ics[i] = ics[j];
      ics[j] = i;
    }
}

void Treedata::permute(vector<num_t>& x)
{
  size_t n(x.size());
  vector<size_t> ics(n);
  vector<num_t> foo = x;
  Treedata::permute(ics);
  for(size_t i = 0; i < n; ++i)
    {
      x[i] = foo[ics[i]];
    }
}

void Treedata::bootstrap(vector<size_t>& ics, vector<size_t>& oob_ics, size_t& noob)
{
  assert(ics.size() == oob_ics.size());
  size_t n(ics.size());
  for(size_t i = 0; i < n; ++i)
    {
      ics[i] = rand() % n;
    }
  sort(ics.begin(),ics.end());
  vector<size_t> allics(n);
  Treedata::range(allics);
  vector<size_t>::iterator it = set_difference(allics.begin(),allics.end(),ics.begin(),ics.end(),oob_ics.begin());
  noob = distance(oob_ics.begin(),it);
}

void Treedata::split_target_wrt_feature(size_t featureidx,
					const size_t min_split,
					vector<size_t>& sampleics,
					vector<size_t>& sampleics_left,
					vector<size_t>& sampleics_right)
{
  
  //We allow the possibility to have target feature as splitter candidate, but we don't allow one to split with it 
  if(targetidx_ == featureidx)
    {
      sampleics_left.clear();
      sampleics_right.clear();
      return;
    }

  //Number of samples
  size_t n_tot(sampleics.size());

  if(isfeaturenum_[featureidx])
    {
    
      //Feature values
      vector<num_t> fv(n_tot);
      
      //Target values
      vector<num_t> tv(n_tot);
  
      //Collect data for the feature and target into the vectors
      for(size_t i = 0; i < n_tot; ++i)
	{
	  fv[i] = featurematrix_[featureidx][sampleics[i]];
	  tv[i] = featurematrix_[targetidx_][sampleics[i]];
	}
      
      //Make reference indices that define the sorting wrt. feature
      vector<size_t> ref_ics(n_tot);
      
      //Sort feature vector and collect reference indices
      Treedata::sort_and_make_ref<num_t>(fv,ref_ics);
      
      //Use the reference indices to sort the target vector and sample indices
      Treedata::sort_from_ref<num_t>(tv,ref_ics);
      Treedata::sort_from_ref<size_t>(sampleics,ref_ics);
      
    }
  else
    {
      
    }  
}

void Treedata::split_target(const size_t min_split,
			    vector<size_t>& sampleics,
			    vector<size_t>& sampleics_left,
			    vector<size_t>& sampleics_right)
{

  sort(sampleics.begin(),sampleics.end());

  if(isfeaturenum_[targetidx_])
    {
      Treedata::incremental_split(min_split,sampleics,sampleics_left,sampleics_right);
    }
  else
    {

      size_t n_tot(sampleics.size());
      size_t n_left(0);
      size_t n_right(n_tot);

      vector<num_t> tv(n_tot);
      for(size_t i = 0; i < n_tot; ++i)
	{
	  tv[i] = featurematrix_[targetidx_][sampleics[i]];
	}

      num_t impurity_left(0.0);
      num_t impurity_right(0.0);
      
      map<num_t,size_t> freq_right;
      map<num_t,size_t> freq_left;
      datadefs::gini(tv,freq_right,impurity_right);
      num_t impurity_tot(impurity_right);
      //size_t ncats_tot(freq_right.size());
      
      set<num_t>::iterator bestcatit;
      num_t bestimpurity(impurity_tot);
      
      set<num_t> categories_right;
      for(map<num_t,size_t>::const_iterator it(freq_right.begin()); it != freq_right.end(); ++it)
        {
          categories_right.insert(it->first);
        }
      
      while(true)
	{
	  
	  size_t ncats_right(categories_right.size());
	  if(ncats_right == 0)
	    {
	      break;
	    }
	  
	  bestcatit = categories_right.end();
	  for(set<num_t>::const_iterator catit(categories_right.begin()); catit != categories_right.end(); ++catit)
	    {
	      
	      size_t n_diff(freq_right[*catit]);
	      
	      n_left += n_diff;
	      n_right -= n_diff;
	      
	      freq_left.insert(pair<num_t,size_t>(*catit,n_diff));
	      freq_right.erase(*catit);
	      
	      datadefs::gini(n_left,freq_left,impurity_left);
	      datadefs::gini(n_right,freq_right,impurity_right);
	      
	      if((n_left*impurity_left+n_right*impurity_right) < n_tot*bestimpurity)
		{
		  bestcatit = catit;
		  bestimpurity = (impurity_left + impurity_right);
		}
	      
	      n_left -= n_diff;
	      n_right += n_diff;
	      
	      freq_left.erase(*catit);
	      freq_right.insert(pair<num_t,size_t>(*catit,n_diff));
	      
	    }
	  if(bestcatit == categories_right.end())
	    {
	      break;
	    }
	  
	  //num_t bestcat(categories_right[bestcatidx]);
	  size_t n_diff(freq_right[*bestcatit]);
	  n_left += n_diff;
	  n_right -= n_diff;
	  
	  freq_left.insert(pair<num_t,size_t>(*bestcatit,n_diff));
	  freq_right.erase(*bestcatit);
	  categories_right.erase(bestcatit);

	  //cout << "catleft=[";
	  //for(map<num_t,size_t>::iterator it(freq_left.begin()); it != freq_left.end(); ++it)
	  //  {
	  //    cout << " " << it->first << ":" << it->second;
	  //  }
	  //cout << " ](" << n_left << ")  catright=[";
	  //for(map<num_t,size_t>::iterator it(freq_right.begin()); it != freq_right.end(); ++it)
	  //  {
	  //    cout << " " << it->first << ":" << it->second;
	  //  }
	  //cout << " ](" << n_right << ")  impurity=" << impurity_left+impurity_right << "  bestimpurity=" << bestimpurity << endl;
	}
      
      set<num_t> categories_left;
      for(map<num_t,size_t>::const_iterator it(freq_left.begin()); it != freq_left.end(); ++it)
        {
          categories_left.insert(it->first);
        }

      Treedata::split_samples(targetidx_,sampleics,categories_left,sampleics_left,sampleics_right);
      
      //cout << "target splitted: left=[";
      //for(size_t i = 0; i < sampleics_left.size(); ++i)
      //  {
      //    cout << " " << sampleics_left[i];
      //  }
      //cout << " ] right=[";
      //for(size_t i = 0; i < sampleics_right.size(); ++i)
      //  {
      //    cout << " " << sampleics_right[i];
      //  }
      //cout << " ]" << endl;
      
    }
}

void Treedata::incremental_split(const size_t min_split, vector<size_t>& sampleics, vector<size_t>& sampleics_left, vector<size_t>& sampleics_right)
{

  //Number of samples
  size_t n_tot(sampleics.size());
  size_t n_right(n_tot);
  size_t n_left(0);

  vector<num_t> tv(n_tot);
  for(size_t i = 0; i < n_tot; ++i)
    {
      tv[i] = featurematrix_[targetidx_][sampleics[i]];
    }
  
  //num_t mu_old(0.0);
  num_t mu_right(0.0);
  num_t mu_left(0.0);
  num_t impurity_left(0.0);
  num_t impurity_right(0.0);
  
  size_t nreal;
  datadefs::sqerr(tv,mu_right,impurity_right,nreal);
  num_t impurity_tot(impurity_right);
  num_t bestimpurity(impurity_tot);
  int bestsplitidx = -1;
  
  while(n_left < n_tot - min_split)
    {
      int idx(n_left);
      ++n_left;
      --n_right;
      datadefs::update_sqerr(tv[idx],n_left,mu_left,impurity_left,n_right,mu_right,impurity_right);
      if((impurity_left+impurity_right) < bestimpurity && n_left >= min_split)
	{
	  bestsplitidx = idx;
	  bestimpurity = (impurity_left + impurity_right);
	}
    }
  
  Treedata::split_samples(sampleics,bestsplitidx,sampleics_left,sampleics_right);
  
}

void Treedata::categorical_split(vector<size_t>& sampleics, vector<size_t>& sampleics_left, vector<size_t>& sampleics_right)
{

}




void Treedata::impurity(size_t featureidx, vector<size_t>& sampleics, num_t& impurity)
{

  size_t n(sampleics.size());
  vector<num_t> data(n);

  for(size_t i = 0; i < n; ++i)
    {
      data[i] = featurematrix_[featureidx][sampleics[i]];
    }

  if(isfeaturenum_[featureidx])
    {
      num_t mu,se;
      size_t nreal;
      datadefs::sqerr(data,mu,se,nreal);
      impurity = se/nreal;
    }
  else
    {
      map<num_t,size_t> freq;
      num_t impurity;
      datadefs::gini(data,freq,impurity);
    }

}

void Treedata::split_samples(vector<size_t>& sampleics, 
			      size_t splitidx, 
			      vector<size_t>& sampleics_left, 
			      vector<size_t>& sampleics_right)
{
  //Split the sample indices
  size_t n_tot(sampleics.size());
  size_t n_left(splitidx + 1);
  size_t n_right(n_tot - n_left);

  if(n_left == 0)
    {
      sampleics_right = sampleics;
      vector<size_t> sampleics_left_copy(0);
      sampleics_left = sampleics_left_copy;
      return;
    }

  vector<size_t> sampleics_left_copy(n_left);
  vector<size_t> sampleics_right_copy(n_right);
  size_t i(0);
  while(i < n_left)
    {
      sampleics_left_copy[i] = sampleics[i];
      ++i;
    }
  while(i < n_tot)
    {
      sampleics_right_copy[i-n_left] = sampleics[i];
      ++i;
    }

  sampleics_left = sampleics_left_copy;
  sampleics_right = sampleics_right_copy;

}

void Treedata::split_samples(size_t featureidx,
			     vector<size_t>& sampleics,
			     set<num_t>& categories_left,
			     vector<size_t>& sampleics_left,
			     vector<size_t>& sampleics_right)
{

  sampleics_left.clear();
  sampleics_right.clear();

  for(size_t i = 0; i < sampleics.size(); ++i)
    {
      if(categories_left.find(featurematrix_[featureidx][sampleics[i]]) != categories_left.end())
	{
	  sampleics_left.push_back(sampleics[i]);
	}
      else
	{
	  sampleics_right.push_back(sampleics[i]);
	}
    }

}

