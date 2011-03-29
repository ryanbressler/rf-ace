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
  contrastmatrix_.resize(nfeatures_);
  Treedata::select_target(targetidx_); //This will sort the data wrt. selected target and generate the contrasts
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

/*
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
*/

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


/*
  void Treedata::range(vector<size_t>& ics)
  {
  for(size_t i = 0; i < ics.size(); ++i)
  {
  ics[i] = i;
  }
  }
*/


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

void Treedata::generate_contrasts()
{
  size_t nrealvalues(nrealvalues_[targetidx_]);
  for(size_t i = 0; i < nfeatures_; ++i)
    {
      vector<num_t> x(nrealvalues);
      contrastmatrix_[i].resize(nrealvalues);
      for(size_t j = 0; j < nrealvalues; ++j)
        {
          contrastmatrix_[i][j] = featurematrix_[i][j];
        }
      Treedata::permute(contrastmatrix_[i]);
    }
}

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

  if(isfeaturenum_[targetidx_])
    {
      Treedata::incremental_target_split(targetidx_,min_split,sampleics,sampleics_left,sampleics_right,splitvalue);
    }
  else
    {
      set<num_t> values_left;
      Treedata::categorical_target_split(targetidx_,sampleics,sampleics_left,sampleics_right,values_left); 
    }
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

  //If target is split wrt. feature, then we first need to sort the target data (indices) wrt. feature data (indices)
  //if(featureidx != targetidx_)
  //  {
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
  //  }
  //else //Otherwise we just sort the sample indices
  //  {
  //    sort(sampleics.begin(),sampleics.end());
  //  }
      
  //Next we collect the target data in the order specified by the reordered sampleics
  vector<num_t> tv(n_tot);
  for(size_t i = 0; i < n_tot; ++i)
    {
      tv[i] = featurematrix_[targetidx_][sampleics[i]];
    }
 
  num_t impurity_left(0.0);
  num_t impurity_right(0.0);

  int bestsplitidx = -1;
 
  //If the target is numerical, we use the iterative squared error formula to update impurity scores while we traverse "right"
  if(isfeaturenum_[targetidx_])
    {
      num_t mu_right(0.0);
      num_t mu_left(0.0);
      
      size_t nreal;
      datadefs::sqerr(tv,mu_right,impurity_right,nreal);
      num_t impurity_tot(impurity_right);
      num_t bestimpurity(impurity_tot);

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
    }
  else //Otherwise we use the iterative gini index formula to update impurity scores while we traverse "right"
    {

      //map<num_t,size_t> freq_tot;
      map<num_t,size_t> freq_left;
      map<num_t,size_t> freq_right;

      datadefs::count_freq(tv,freq_right);
      datadefs::gini(freq_right,impurity_right);
      num_t impurity_tot(impurity_right);
      num_t bestimpurity(impurity_tot);

      while(n_left < n_tot - min_split)
        {
          int idx(n_left);
          ++n_left;
          --n_right;
	  datadefs::update_gini(tv[idx],n_left,freq_left,impurity_left,n_right,freq_right,impurity_right); //INCREMENTAL GINI INDEX COMPUTATION IS INEFFICIENT AT THE MOMENT
          if((n_left*impurity_left+n_right*impurity_right) < n_tot*bestimpurity && n_left >= min_split) //THIS NEEDS TO BE FIXED
            {
              bestsplitidx = idx;
              bestimpurity = (n_left*impurity_left+n_right*impurity_right) / n_tot;
            }
        }      
    }

  //Finally, make the split at the best point
  Treedata::split_samples(sampleics,bestsplitidx,sampleics_left,sampleics_right);
  
  //Return the split value
  splitvalue = fv[bestsplitidx];

  cout << "[";
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
  
  num_t impurity_tot;

  //Target data
  vector<num_t> tv(n_tot);
  for(size_t i = 0; i < n_tot; ++i)
    {
      tv[i] = featurematrix_[targetidx_][sampleics[i]];
    }

  if(isfeaturenum_[targetidx_])
    {
      num_t mu;
      size_t nreal;
      datadefs::sqerr(tv,mu,impurity_tot,nreal);
      assert(nreal == n_tot);
      impurity_tot /= n_tot;
    }
  else
    {
      datadefs::gini(tv,impurity_tot);
    }
  
  //Feature data
  vector<num_t> fv(n_tot);
  for(size_t i = 0; i < n_tot; ++i)
    {
      fv[i] = featurematrix_[featureidx][sampleics[i]];
    }

  //Impurity scores for the left and right branches
  num_t impurity_left(0.0);
  num_t impurity_right(0.0);
  
  map<num_t,vector<size_t> > fmap;
  datadefs::map_data(fv,fmap);

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
	  if(isfeaturenum_[targetidx_])
	    {
	      num_t mu;
	      size_t nreal;
	      datadefs::sqerr(data_left,mu,impurity_left,nreal);
	      assert(nreal == data_left.size());
	      impurity_left /= nreal;
	      datadefs::sqerr(data_right,mu,impurity_right,nreal);
	      assert(nreal == data_right.size());
	      impurity_right /= nreal;
	    }
	  else
	    {
	      datadefs::gini(data_left,impurity_left);
	      datadefs::gini(data_right,impurity_right);
	    }	    

	  size_t n_left(data_left.size());
          size_t n_right(data_right.size());
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
 
  cout << "[";
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


void Treedata::percolate_sampleidx(size_t sampleidx, Node** nodep)
{
  //cout << nodep->has_children() << endl;
  while((*nodep)->has_children())
    {
      size_t featureidx((*nodep)->get_splitter());
      num_t value(featurematrix_[featureidx][sampleidx]);
      //cout << featureidx << ":" << value << endl;
      *nodep = (*nodep)->percolate(value);
      //cout << &*nodep << endl;
    }
}


void Treedata::impurity(size_t featureidx, vector<size_t>& sampleics, num_t& impurity)
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
      size_t nreal;
      datadefs::sqerr(data,mu,se,nreal);
      impurity = se/nreal;
    }
  else
    {
      map<num_t,size_t> freq;
      num_t impurity;
      //datadefs::count_freq(data,freq);
      datadefs::gini(data,impurity);
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



