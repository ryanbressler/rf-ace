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
  contrastmatrix_(0),
  isfeaturenum_(0),
  nsamples_(0),
  nfeatures_(0),
  featureheaders_(0),
  sampleheaders_(0)
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

  //Transform raw data to the internal format.
  featurematrix_.resize(nfeatures_);
  contrastmatrix_.resize(nfeatures_);
  for(size_t i = 0; i < nfeatures_; ++i)
    {
      vector<num_t> featurev(nsamples_);
      //map<string,size_t> str2valmap;
      if(isfeaturenum_[i])
        {
          datadefs::strv2numv(rawmatrix[i],featurev);
          //ncatvalues_[i] = 0;
        }
      else
        {
          datadefs::strv2catv(rawmatrix[i],featurev);
	  //ncatvalues_[i] = str2valmap.size();
        }
      featurematrix_[i] = featurev;
    } 
  
  /*
    featureheaders_.resize(2*nfeatures_);
    isfeaturenum_.resize(2*nfeatures_);
    for(size_t i = nfeatures_; i < 2*nfeatures_; ++i)
    {
    featurematrix_[i] = featurematrix_[i-nfeatures_];
    featureheaders_[i] = "CONTRAST";
    isfeaturenum_[i] = isfeaturenum_[i-nfeatures_];
    }
  */
    
  //Treedata::permute_contrasts();

  targetidx_ = 0;
  Treedata::select_target(targetidx_);

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

string Treedata::get_featureheader(size_t featureidx)
{
  return(featureheaders_[featureidx]);
}

string Treedata::get_targetheader()
{
  return(featureheaders_[targetidx_]);
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

void Treedata::generate_contrasts()
{
  size_t nrealvalues;
  datadefs::count_real_values(featurematrix_[targetidx_],nrealvalues);
  for(size_t i = 0; i < nfeatures_; ++i)
    {
      vector<num_t> x(nrealvalues);
      contrastmatrix_[i].resize(nrealvalues);
      for(size_t j = 0; j < nrealvalues; ++j)
        {
          contrastmatrix_[i][j] = featurematrix_[i][j];
        }
    }

  Treedata::permute_contrasts();
}

void Treedata::permute_contrasts()
{
  for(size_t i = 0; i < nfeatures_; ++i)
    {
      Treedata::permute(contrastmatrix_[i]);
    }
}

void Treedata::select_target(size_t targetidx)
{
  targetidx_ = targetidx; 
  Treedata::sort_all_wrt_target();
  Treedata::generate_contrasts();
}

size_t Treedata::get_target()
{
  return(targetidx_);
}

bool Treedata::isfeaturenum(size_t featureidx)
{
  return(isfeaturenum_[featureidx]);
}

size_t Treedata::nrealvalues()
{
  size_t nreal;
  datadefs::count_real_values(featurematrix_[targetidx_],nreal);
  return(nreal);
}

size_t Treedata::nrealvalues(size_t featureidx)
{
  size_t nreal;
  datadefs::count_real_values(featurematrix_[featureidx],nreal);
  return(nreal);
}

void Treedata::remove_nans(size_t featureidx, 
			   vector<size_t>& sampleics, 
			   size_t& nreal)
{
  
  nreal = sampleics.size();
  int maxidx(nreal-1);
  for(int i = maxidx; i >= 0; --i)
    {
      if(datadefs::is_nan(featurematrix_[featureidx][sampleics[i]]))
	{
	  --nreal;
	  sampleics.erase(sampleics.begin()+i);
	}
    }
}

void Treedata::sort_all_wrt_feature(size_t featureidx)
{
  //assert(featureidx < nfeatures_);

  //Generate and index vector that'll define the new order
  vector<size_t> ref_ics(nsamples_);
  
  //Make a dummy variable so that the original feature matrix will stay intact 
  vector<num_t> dummy(nsamples_);

  dummy = featurematrix_[featureidx];

  //Sort specified feature and make reference indices
  datadefs::sort_and_make_ref<num_t>(dummy,ref_ics);

  //Use the new order to sort the sample headers
  datadefs::sort_from_ref<string>(sampleheaders_,ref_ics);

  //Sort features
  for(size_t i = 0; i < nfeatures_; ++i)
    {
      datadefs::sort_from_ref<num_t>(featurematrix_[i],ref_ics);
    }
}

void Treedata::sort_all_wrt_target()
{
  Treedata::sort_all_wrt_feature(targetidx_);
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
  //vector<num_t> foo = x;
  Treedata::permute(ics);
  for(size_t i = 0; i < n; ++i)
    {
      num_t temp(x[i]);
      x[i] = x[ics[i]];
      x[ics[i]] = temp;
    }
}

/*
  void Treedata::generate_contrasts()
  {
  size_t nrealvalues;
  datadefs::count_real_values(featurematrix_[targetidx_],nrealvalues);
  for(size_t i = 0; i < nfeatures_; ++i)
  {
  vector<num_t> x(nrealvalues);
  contrastmatrix_[i].resize(nrealvalues);
  for(size_t j = 0; j < nrealvalues; ++j)
  {
  contrastmatrix_[i][j] = featurematrix_[i][j];
  }
  }
  
  Treedata::permute_contrasts();
  }
*/


void Treedata::bootstrap(vector<size_t>& ics, vector<size_t>& oob_ics)
{
  //assert(ics.size() == oob_ics.size());
  size_t n(ics.size());
  for(size_t i = 0; i < n; ++i)
    {
      ics[i] = rand() % n;
    }
  sort(ics.begin(),ics.end());
  vector<size_t> allics(n);
  datadefs::range(allics);

  vector<size_t> foo(n);
  vector<size_t>::iterator it = set_difference(allics.begin(),allics.end(),ics.begin(),ics.end(),foo.begin());
  size_t noob = distance(foo.begin(),it);

  foo.resize(noob);
  oob_ics = foo;
}

void Treedata::split_target_with_num_feature(size_t featureidx,
					     const size_t min_split,
					     vector<size_t>& sampleics,
					     vector<size_t>& sampleics_left,
					     vector<size_t>& sampleics_right,
					     num_t& splitvalue)
{
  
  assert(isfeaturenum_[featureidx]);
  assert(featureidx != targetidx_);
  
  //Number of samples
  size_t n_tot(sampleics.size());
  
  //Feature values
  vector<num_t> fv(n_tot);
  
  //Collect data for the feature and target into the vectors
  for(size_t i = 0; i < n_tot; ++i)
    {
      fv[i] = featurematrix_[featureidx][sampleics[i]];
    }
  
  //Make reference indices that define the sorting wrt. feature
  vector<size_t> ref_ics(n_tot);
  
  //Sort feature vector and collect reference indices
  datadefs::sort_and_make_ref<num_t>(fv,ref_ics);
  
  //Use the reference indices to sort sample indices
  datadefs::sort_from_ref<size_t>(sampleics,ref_ics);
  
  Treedata::incremental_target_split(featureidx,min_split,sampleics,sampleics_left,sampleics_right,splitvalue);
  
}

void Treedata::split_target_with_cat_feature(size_t featureidx,
					     const size_t min_split,
					     vector<size_t>& sampleics,
					     vector<size_t>& sampleics_left,
					     vector<size_t>& sampleics_right,
					     set<num_t>& values_left)
{
  assert(!isfeaturenum_[featureidx]);
  assert(featureidx != targetidx_);
  
  Treedata::categorical_target_split(featureidx,sampleics,sampleics_left,sampleics_right,values_left);

}

void Treedata::split_target(const size_t min_split,
			    vector<size_t>& sampleics,
			    vector<size_t>& sampleics_left,
			    vector<size_t>& sampleics_right)
{

  num_t splitvalue;

  size_t n_tot(sampleics.size());

  if(isfeaturenum_[targetidx_])
    {
      Treedata::incremental_target_split(targetidx_,min_split,sampleics,sampleics_left,sampleics_right,splitvalue);
    }
  else
    {
      set<num_t> values_left;
      Treedata::categorical_target_split(targetidx_,sampleics,sampleics_left,sampleics_right,values_left); 
    }

  assert(sampleics.size() == n_tot);
  assert(sampleics_left.size() + sampleics_right.size() == n_tot);
}

void Treedata::incremental_target_split(size_t featureidx,
					const size_t min_split, 
					vector<size_t>& sampleics, 
					vector<size_t>& sampleics_left, 
					vector<size_t>& sampleics_right,
					num_t& splitvalue)
{

  //The splitter feature needs to be numerical
  assert(isfeaturenum_[featureidx]);

  //Number of samples
  size_t n_tot(sampleics.size());
  size_t n_right(n_tot);
  size_t n_left(0);

  //Check that there are enough samples to make the split in the first place
  assert(n_tot >= 2*min_split);

  //Feature values
  vector<num_t> fv(n_tot);
  
  //Collect data for the feature and target into the vectors
  for(size_t i = 0; i < n_tot; ++i)
    {
      fv[i] = featurematrix_[featureidx][sampleics[i]];
    }
  
  //Make reference indices that define the sorting wrt. feature
  vector<size_t> ref_ics(n_tot);
  
  //Sort feature vector and collect reference indices
  datadefs::sort_and_make_ref<num_t>(fv,ref_ics);
  
  //Use the reference indices to sort sample indices
  datadefs::sort_from_ref<size_t>(sampleics,ref_ics);
        
  //Next we collect the target data in the order specified by the reordered sampleics
  vector<num_t> tv(n_tot);
  for(size_t i = 0; i < n_tot; ++i)
    {
      tv[i] = featurematrix_[targetidx_][sampleics[i]];
    }

  //Count how many real values the feature and target has
  size_t nreal_f,nreal_t;
  datadefs::count_real_values(fv,nreal_f);
  datadefs::count_real_values(tv,nreal_t);
  assert(nreal_t == n_tot && nreal_f == n_tot);

  int bestsplitidx = -1;
 
  //If the target is numerical, we use the iterative squared error formula to update impurity scores while we traverse "right"
  if(isfeaturenum_[targetidx_])
    {
      num_t mu_right(0.0);
      num_t se_right(0.0);
      num_t mu_left(0.0);
      num_t se_left(0.0);
      num_t se_best(0.0);
      size_t nreal_right(0);
      
      datadefs::sqerr(tv,mu_right,se_right,nreal_right);
      assert(n_tot == nreal_right);
      se_best = se_right;

      size_t idx(0);
      while(n_left < n_tot - min_split)
	{
	  datadefs::update_sqerr(tv[idx],n_left,mu_left,se_left,n_right,mu_right,se_right);
	  if((se_left+se_right) < se_best && n_left >= min_split)
	    {
	      bestsplitidx = idx;
	      se_best = se_left + se_right;
	    }
	  ++idx;
	}
    }
  else //Otherwise we use the iterative gini index formula to update impurity scores while we traverse "right"
    {
      map<num_t,size_t> freq_left;
      map<num_t,size_t> freq_right;
      num_t sf_left(0.0);
      num_t sf_right(0.0);
      //num_t sf_best(0.0);
      size_t nreal_right(0);
 
      datadefs::sqfreq(tv,freq_right,sf_right,nreal_right);
      num_t nsf_best(sf_right/nreal_right);
      assert(n_tot == nreal_right);
      
      size_t idx(0);
      while(n_left < n_tot - min_split)
        {
	  datadefs::update_sqfreq(tv[idx],n_left,freq_left,sf_left,n_right,freq_right,sf_right);
          if(sf_left/n_left + sf_right/n_right > nsf_best && n_left >= min_split) 
            {
              bestsplitidx = idx;
              nsf_best = sf_left/n_left + sf_right/n_right;
            }
	  ++idx;
        }      
    }

  //Finally, make the split at the best point
  Treedata::split_samples(sampleics,bestsplitidx,sampleics_left,sampleics_right);
  
  assert(sampleics.size() == n_tot);
  assert(sampleics_left.size() + sampleics_right.size() == n_tot);

  //Return the split value
  splitvalue = fv[bestsplitidx];

  cout << "Feature " << featureidx << " splits target " << targetidx_ << " [";
  for(size_t i = 0; i < sampleics_left.size(); ++i)
    {
      cout << " " << featurematrix_[targetidx_][sampleics_left[i]];
    }
  cout << " ] <==> [";
  for(size_t i = 0; i < sampleics_right.size(); ++i)
    {
      cout << " " << featurematrix_[targetidx_][sampleics_right[i]];
    }
  cout << " ]" << endl;

}

//THIS NEEDS REWORKING AND OPTIMIZING
void Treedata::categorical_target_split(size_t featureidx, 
					vector<size_t>& sampleics, 
					vector<size_t>& sampleics_left, 
					vector<size_t>& sampleics_right,
					set<num_t>& categories_left)
{

  assert(!isfeaturenum_[featureidx]);

  //Sample size
  size_t n_tot(sampleics.size());
  
  //Check that sample size is positive
  assert(n_tot > 0);
  //size_t n_left(0);

  //Feature data
  vector<num_t> fv(n_tot);
  for(size_t i = 0; i < n_tot; ++i)
    {
      fv[i] = featurematrix_[featureidx][sampleics[i]];
    }
  
  //Make reference indices that define the sorting wrt. feature
  vector<size_t> ref_ics(n_tot);

  //Sort feature vector and collect reference indices
  datadefs::sort_and_make_ref<num_t>(fv,ref_ics);

  //Use the reference indices to sort sample indices
  datadefs::sort_from_ref<size_t>(sampleics,ref_ics);
  
  //Target data
  vector<num_t> tv(n_tot);
  for(size_t i = 0; i < n_tot; ++i)
    {
      tv[i] = featurematrix_[targetidx_][sampleics[i]];
    }

  //Count how many real values the feature and target has
  size_t nreal_f,nreal_t;
  datadefs::count_real_values(fv,nreal_f);
  datadefs::count_real_values(tv,nreal_t);
  assert(nreal_t == n_tot && nreal_f == n_tot);

  size_t nreal;
  num_t impurity_tot;

  if(isfeaturenum_[targetidx_])
    {
      num_t mu;
      //size_t nreal;
      datadefs::sqerr(tv,mu,impurity_tot,nreal);
      assert(nreal == n_tot);
      impurity_tot /= n_tot;
    }
  else
    {
      datadefs::gini(tv,impurity_tot,nreal);
      assert(nreal == n_tot);
    }
  
  //Impurity scores for the left and right branches
  num_t impurity_left(0.0);
  num_t impurity_right(0.0);
  
  map<num_t,vector<size_t> > fmap;
  datadefs::map_data(fv,fmap,n_tot);

  for(map<num_t,vector<size_t> >::iterator it(fmap.begin()); it != fmap.end(); ++it)
    {
      //cout << "fmap(" << it->first << ") [";
      for(size_t i = 0; i < it->second.size(); ++i)
	{
	  it->second[i] = sampleics[it->second[i]];
	  //cout << " " << it->second[i];
	}
      //cout << " ]" << endl;
    }

  
  map<num_t,vector<size_t> > fmap_left;
  map<num_t,vector<size_t> > fmap_right = fmap;

  num_t bestimpurity(impurity_tot);

  while(fmap_right.size() > 1)
    {

      map<num_t,vector<size_t> > fmap_right_copy = fmap_right;
      
      //num_t bestkey(datadefs::num_nan);
      map<num_t,vector<size_t> >::iterator bestit(fmap_right.end());

      for(map<num_t,vector<size_t> >::iterator it(fmap_right.begin()); it != fmap_right.end(); ++it)
	{

	  //num_t key(it->first);
	  fmap_left.insert(*it);
	  fmap_right_copy.erase(it->first);
	  //cout << "[";
	  vector<num_t> data_left;
	  //for(map<num_t,vector<size_t> >::const_iterator it2(fmap_left.begin()); it2 != fmap_left.end(); ++it2)
	  //  {
	      //cout << " " << it2->first << "(";
	  //    for(size_t i = 0; i < it2->second.size(); ++i)
	  //	{
	  //	  data_left.push_back(featurematrix_[targetidx_][it2->second[i]]);
	  //	  cout << " " << featurematrix_[targetidx_][it2->second[i]];
	  //	}
	  //    //cout << ")";
	  //  }
	  //cout << " ] <==> [";
	  vector<num_t> data_right;
          //for(map<num_t,vector<size_t> >::const_iterator it2(fmap_right_copy.begin()); it2 != fmap_right_copy.end(); ++it2)
          //  {
	  //    //cout << " " << it2->first << "(";
          //    for(size_t i = 0; i < it2->second.size(); ++i)
          //      {
          //        data_right.push_back(featurematrix_[targetidx_][it2->second[i]]);
	  //		  cout << " " << featurematrix_[targetidx_][it2->second[i]];
	  //	}
	      //cout << ")";
	  // }
	  //cout << " ] impurity_left=";

	  //num_t impurity_left, impurity_right;
	  
	  size_t n_left,n_right;
	  
	  if(isfeaturenum_[targetidx_])
	    {
	      num_t mu;
	      //size_t nreal;
	      datadefs::sqerr(data_left,mu,impurity_left,n_left);
	      //assert(nreal == );
	      impurity_left /= n_left;
	      datadefs::sqerr(data_right,mu,impurity_right,n_right);
	      //assert(nreal == data_right.size());
	      impurity_right /= n_right;
	    }
	  else
	    {
	      datadefs::gini(data_left,impurity_left,n_left);
	      datadefs::gini(data_right,impurity_right,n_right);
	    }	    

	  //size_t n_left(data_left.size());
          //size_t n_right(data_right.size());
	  num_t impurity_new = (n_left*impurity_left+n_right*impurity_right) / n_tot;

	  //cout << impurity_left << "  impurity_right=" << impurity_right << " (total=" << impurity_new << "\tcurr.best=" << bestimpurity << ")" << endl;
	  
	  if(impurity_new < bestimpurity)
	    {

	      bestit = it;
	      bestimpurity = impurity_new;
	    }
	  
	  fmap_left.erase(it->first);
	  fmap_right_copy.insert(*it);
	 	  
	}

      if(bestit == fmap_right.end())
	{
	  break;
	}

      //cout << "it's now...";
      fmap_left.insert(*bestit);
      fmap_right.erase(bestit->first);
      //cout << " or never." << endl;
      
      

    }
 
  sampleics_left.clear();
  categories_left.clear();
  for(map<num_t,vector<size_t> >::const_iterator it(fmap_left.begin()); it != fmap_left.end(); ++it)
    {
      categories_left.insert(it->first);
      //cout << " " << it2->first << "(";
      for(size_t i = 0; i < it->second.size(); ++i)
	{
	  sampleics_left.push_back(it->second[i]);
	  //cout << " " << featurematrix_[targetidx_][it2->second[i]];
	}
      //cout << ")";
    }

  //cout << "right" << endl;

  sampleics_right.clear();
  for(map<num_t,vector<size_t> >::const_iterator it(fmap_right.begin()); it != fmap_right.end(); ++it)
    {
      //cout << " " << it2->first << "(";
      for(size_t i = 0; i < it->second.size(); ++i)
        {
          sampleics_right.push_back(it->second[i]);
          //cout << " " << featurematrix_[targetidx_][it2->second[i]];
        }
      //cout << ")";
    }


  //Treedata::split_samples(targetidx_,sampleics,categories_left,sampleics_left,sampleics_right);
 
  cout << "Feature " << featureidx << " splits target " << targetidx_ << " [";
  for(size_t i = 0; i < sampleics_left.size(); ++i)
    {
      cout << " " << featurematrix_[targetidx_][sampleics_left[i]];
    }
  cout << " ] <==> [";
  for(size_t i = 0; i < sampleics_right.size(); ++i)
    {
      cout << " " << featurematrix_[targetidx_][sampleics_right[i]];
    }
  cout << " ]" << endl;
  
}

num_t Treedata::at(size_t featureidx, size_t sampleidx)
{
  return(featurematrix_[featureidx][sampleidx]);
}


num_t Treedata::atp(size_t featureidx, size_t sampleidx)
{
  return(contrastmatrix_[featureidx][sampleidx]);
}


void Treedata::impurity(size_t featureidx, vector<size_t> const& sampleics, num_t& impurity, size_t& nreal)
{

  size_t n(sampleics.size());
  vector<num_t> data(n);

  for(size_t i = 0; i < n; ++i)
    {
      data[i] = featurematrix_[featureidx][sampleics[i]];
    }

  //cout << "before impurity -- ";
  if(isfeaturenum_[featureidx])
    {
      num_t mu,se;
      //size_t nreal;
      datadefs::sqerr(data,mu,se,nreal);
      impurity = se/nreal;
    }
  else
    {
      map<num_t,size_t> freq;
      num_t sf;
      //datadefs::count_freq(data,freq);
      datadefs::sqfreq(data,freq,sf,nreal);
      impurity = 1-sf/pow(nreal,2);
    }
  //cout << "after impurity" << endl;
  
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



