#include "treedata.hpp"
#include <cstdlib>
#include <fstream>
#include <cassert>
#include <iostream>
#include <sstream>
#include <utility>
#include <algorithm>
#include <ctime>

using namespace std;

Treedata::Treedata(string filename):
  targetidx_(0),
  featurematrix_(0),
  isfeaturenum_(0),
  nrealsamples_(0),
  nsamples_(0),
  nfeatures_(0),
  featureheaders_(0),
  sampleheaders_(0)
{

  //Initialize random number rgenerator
  time_t now;
  time(&now);
  //srand((unsigned int)now);
  
  //MTRand_int32 irand((unsigned int)now);
  irand_.seed((unsigned int)now);
  //datadefs::seedMT((size_t)now);

  //cout << "Treedata: reading matrix from file '" << filename << "'" << endl;

  //Initialize stream to read from file
  ifstream featurestream;
  featurestream.open(filename.c_str());
  assert(featurestream.good());

  //TODO: a sniffer function that find out the type of the input file based on its content.
  //*************************************************************************************** 
  Filetype filetype = UNKNOWN;
  Treedata::read_filetype(filename,filetype);

  vector<vector<string> > rawmatrix;
  if(filetype == AFM)
    {
      cout << "File type interpreted as Annotated Feature Matrix (AFM)" << endl;
      Treedata::read_afm(featurestream,rawmatrix);
    }
  else if(filetype == ARFF)
    {
      cout << "File type interpreted as Attribute-Relation File Format (ARFF)" << endl;
      Treedata::read_arff(featurestream,rawmatrix);
    }
  else
    {
      cout << "File type is unknown -- defaulting to Annotated Feature Matrix (AFM)" << endl;
      Treedata::read_afm(featurestream,rawmatrix);
    }      

  //TODO: move the following part, generation of artificial contrasts, to a separate function.
  //***************************************************** 

  //Transform raw data to the internal format.
  featurematrix_.resize(2*nfeatures_);
  //contrastmatrix_.resize(nfeatures_);
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
  
  featureheaders_.resize(2*nfeatures_);
  isfeaturenum_.resize(2*nfeatures_);
  for(size_t i = nfeatures_; i < 2*nfeatures_; ++i)
    {
      featurematrix_[i] = featurematrix_[i-nfeatures_];
      featureheaders_[i] = "CONTRAST";
      isfeaturenum_[i] = isfeaturenum_[i-nfeatures_];
    }
      
  //Treedata::permute_contrasts();

  targetidx_ = 0;
  Treedata::select_target(targetidx_);

}

Treedata::~Treedata()
{
}

void Treedata::read_filetype(string& filename, Filetype& filetype)
{

  stringstream ss(filename);
  string suffix = "";
  while(getline(ss,suffix,'.')) {}
  //datadefs::toupper(suffix);

  if(suffix == "AFM" || suffix == "afm")
    {
      filetype = AFM;
    }
  else if(suffix == "ARFF" || suffix == "arff")
    {
      filetype = ARFF;
    }
  else
    {
      filetype = UNKNOWN;
    }

}

void Treedata::read_afm(ifstream& featurestream, vector<vector<string> >& rawmatrix)
{

  string field;
  string row;

  //Remove upper left element from the matrix as useless
  getline(featurestream,field,'\t');

  //Count the number of columns
  getline(featurestream,row);
  stringstream ss(row);
  size_t ncols = 0;
  bool is_rowformatted = true;
  while(getline(ss,field,'\t'))
    {
      if(is_rowformatted && is_featureheader(field))
	{
	  is_rowformatted = false;
	}
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
  //vector<vector<string> > rawmatrix(nrows);
  rawmatrix.resize(nrows);
  for(size_t i = 0; i < nrows; ++i)
    {
      rawmatrix[i] = colheaders;
    }

  //cout << "read " << rawmatrix.size() << " rows and " << rawmatrix[0].size() << " columns." << endl;

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
  if(is_rowformatted)
    {
      cout << "AFM orientation: features as rows" << endl;

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
          else if(rowheaders[i][0] == 'C' || rowheaders[i][0] == 'B')
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

      cout << "AFM orientation: features as columns" << endl;
      //cout << "Transposing...";
      
      Treedata::transpose<string>(rawmatrix);
      //cout << "done" << endl;

      //We can extract the number of features and samples from the row and column counts, respectively
      nfeatures_ = ncols;
      nsamples_ = nrows;

      //Thus, sample headers are column headers
      sampleheaders_ = rowheaders;

      //... and feature headers are row headers
      featureheaders_ = colheaders;

      for(size_t i = 0; i < nfeatures_; ++i)
        {
          //First letters in the row headers determine whether the feature is numerical or categorical
          if(colheaders[i][0] == 'N')
            {
              isfeaturenum_.push_back(true);
            }
          else if(colheaders[i][0] == 'C' || colheaders[i][0] == 'B')
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
}

void Treedata::read_arff(ifstream& featurestream, vector<vector<string> >& rawmatrix)
{

  //string field;
  string row;

  bool hasrelation = false;
  size_t nattributes = 0;
  bool hasdata = false;

  while(getline(featurestream,row))
    {

      //cout << row << endl;

      if(hasdata && hasrelation)
	{
	  if(nattributes == 0)
	    {
	      cerr << "no attributes found from the ARFF file" << endl;
	      assert(false);
	    }

	  rawmatrix.resize(nattributes);
	  
	  while(getline(featurestream,row))
	    {
	      string field;
	      stringstream ss(row);
	      for(size_t attribute = 0; attribute < nattributes; ++attribute)
		{
		  getline(ss,field,',');
		  rawmatrix[attribute].push_back(field);
		}
	    }

	  break;
	}

      if(row[0] == '%' || row == "")
	{
	  continue;
	}

      if(!hasrelation && row.compare(0,9,"@relation") == 0)
	{
	  hasrelation = true;
	  cout << "found relation header: " << row << endl;
	}
      else if(row.compare(0,10,"@attribute") == 0)
	{
	  ++nattributes;
	  cout << "found attribute header: " << row << endl;
	}
      else if(!hasdata && row.compare(0,5,"@data") == 0)
	{
	  hasdata = true;
	  cout << "found data header:" << row << endl;
	}
      else
	{
	  cout << "problem" << endl;
	  assert(false);
	}
    }
}

void Treedata::parse_arff_attribute(const string& str, vector<string>& fields)
{

  fields.clear();

  stringstream ss(str);
  string token = "";

  while(getline(ss,token,' '))
    {
      fields.push_back(token);
      cout << "'" << token << "'";
    }
  cout << endl;
}

bool Treedata::is_featureheader(const string& str)
{
  return((str[0] == 'N' || str[0] == 'C' || str[0] == 'B') && str[1] == ':');
}

size_t Treedata::nfeatures()
{
  return(nfeatures_);
}

size_t Treedata::nsamples()
{
  return(nsamples_);
}

num_t Treedata::corr(size_t featureidx1, size_t featureidx2)
{
  num_t r;
  datadefs::pearson_correlation(featurematrix_[featureidx1],featurematrix_[featureidx2],r);
  return(r);
}

void Treedata::kill(const size_t featureidx)
{
  assert(featureidx != targetidx_);
  assert(featureidx < nfeatures_);

  featurematrix_.erase(featurematrix_.begin() + nfeatures_ + featureidx);
  featureheaders_.erase(featureheaders_.begin() + nfeatures_ + featureidx);
  
  featurematrix_.erase(featurematrix_.begin() + featureidx);
  featureheaders_.erase(featureheaders_.begin() + featureidx);
  
  --nfeatures_;
  if(featureidx < targetidx_)
    {
      --targetidx_;
    }
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


void Treedata::print(const size_t featureidx)
{
  cout << "Print " << featureheaders_[featureidx] << ":";
  for(size_t i = 0; i < nsamples_; ++i)
    {
      cout << " " << featurematrix_[featureidx][i];
    }
  cout << endl;
}


void Treedata::permute_contrasts()
{
  for(size_t i = nfeatures_; i < 2*nfeatures_; ++i)
    {
      Treedata::permute(featurematrix_[i]);
    }
}


void Treedata::select_target(size_t targetidx)
{
  targetidx_ = targetidx; 
  Treedata::sort_all_wrt_target();
  datadefs::count_real_values(featurematrix_[targetidx_],nrealsamples_);
  Treedata::permute_contrasts();
}


size_t Treedata::get_target()
{
  return(targetidx_);
}


bool Treedata::isfeaturenum(size_t featureidx)
{
  return(isfeaturenum_[featureidx]);
}


size_t Treedata::nrealsamples()
{
  //size_t nrealsamples;
  //datadefs::count_real_values(featurematrix_[targetidx_],nreal);
  return(nrealsamples_);
}


size_t Treedata::nrealsamples(size_t featureidx)
{ 
  size_t nrealsamples;
  datadefs::count_real_values(featurematrix_[featureidx],nrealsamples);
  return(nrealsamples);
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

template <typename T> void Treedata::transpose(vector<vector<T> >& mat)
{

  vector<vector<T> > foo = mat;

  size_t ncols = mat.size();
  size_t nrows = mat[0].size();

  mat.resize(nrows);
  for(size_t i = 0; i < nrows; ++i)
    {
      mat[i].resize(ncols);
    }

  for(size_t i = 0; i < nrows; ++i)
    {
      for(size_t j = 0; j < ncols; ++j)
	{
	  mat[i][j] = foo[j][i];
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
      size_t j = irand_() % (i + 1);
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
  size_t n = ics.size();
  for(size_t i = 0; i < n; ++i)
    {
      ics[i] = irand_() % n;
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

void Treedata::split_target(size_t featureidx,
			    const size_t min_split,
			    vector<size_t>& sampleics,
			    vector<size_t>& sampleics_left,
			    vector<size_t>& sampleics_right,
			    num_t& splitvalue,
			    set<num_t>& values_left)
{  
  if(isfeaturenum_[featureidx])
    {
      Treedata::incremental_target_split(featureidx,min_split,sampleics,sampleics_left,sampleics_right,splitvalue);
    }
  else
    {
      Treedata::categorical_target_split(featureidx,sampleics,sampleics_left,sampleics_right,values_left);
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
  //cout << nreal_f << " " << nreal_t << " " << n_tot << endl;
  assert(nreal_t == n_tot && nreal_f == n_tot);

  int bestsplitidx = -1;
 
  //If the target is numerical, we use the iterative squared error formula to update impurity scores while we traverse "right"
  if(isfeaturenum_[targetidx_])
    {
      num_t mu_right = 0.0;
      num_t se_right = 0.0;
      num_t mu_left = 0.0;
      num_t se_left = 0.0;
      num_t se_best = 0.0;
      size_t nreal_right = 0;
      
      datadefs::sqerr(tv,mu_right,se_right,nreal_right);
      assert(n_tot == nreal_right);
      se_best = se_right;

      size_t idx = 0;
      while(n_left < n_tot - min_split)
	{
	  datadefs::forward_backward_sqerr(tv[idx],n_left,mu_left,se_left,n_right,mu_right,se_right);
	  if( se_left + se_right < se_best && n_left >= min_split)
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
      size_t sf_left = 0;
      size_t sf_right = 0;
      size_t nreal_right = 0;
 
      datadefs::sqfreq(tv,freq_right,sf_right,nreal_right);
      num_t nsf_best = 1.0 * sf_right / nreal_right;
      assert(n_tot == nreal_right);
      
      size_t idx = 0;
      while(n_left < n_tot - min_split)
        {
	  datadefs::forward_backward_sqfreq(tv[idx],n_left,freq_left,sf_left,n_right,freq_right,sf_right);
          if(1.0 * n_right * sf_left + 1.0 * n_left * sf_right > n_left * n_right * nsf_best && n_left >= min_split) 
            {
              bestsplitidx = idx;
              nsf_best = 1.0 * sf_left / n_left + 1.0 * sf_right / n_right;
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

  if(false)
    {
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

  //Feature data
  vector<num_t> fv(n_tot);
  for(size_t i = 0; i < n_tot; ++i)
    {
      fv[i] = featurematrix_[featureidx][sampleics[i]];
    }

  //Target data
  vector<num_t> tv(n_tot);
  for(size_t i = 0; i < n_tot; ++i)
    {
      tv[i] = featurematrix_[targetidx_][sampleics[i]];
    }

  //Count how many real values the feature and target has (this is just to make sure there are no mising values thrown in)
  size_t nreal_f,nreal_t;
  datadefs::count_real_values(fv,nreal_f);
  datadefs::count_real_values(tv,nreal_t);
  assert(nreal_t == n_tot && nreal_f == n_tot);

  //Map all feature categories to the corresponding samples and represent it as map. The map is used to assign samples to left and right branches
  map<num_t,vector<size_t> > fmap_right;
  datadefs::map_data(fv,fmap_right,nreal_f);
  assert(n_tot == nreal_f);

  categories_left.clear();
  sampleics_left.clear();
  size_t n_left(0);
  size_t n_right(n_tot);

  if(isfeaturenum_[targetidx_])
    {
      
      num_t mu_right;
      num_t mu_left(0.0);
      num_t se_right;
      num_t se_left(0.0);
      datadefs::sqerr(tv,mu_right,se_right,n_right);
      assert(n_tot == n_right);
      num_t se_best(se_right);
      
      while(fmap_right.size() > 0)
	{
	  
	  map<num_t,vector<size_t> >::iterator it_best(fmap_right.end());

	  for(map<num_t,vector<size_t> >::iterator it(fmap_right.begin()); it != fmap_right.end(); ++it)
	    {
	      
	      //Take samples from right and put them left
	      for(size_t i = 0; i < it->second.size(); ++i)
		{
		  datadefs::forward_backward_sqerr(tv[it->second[i]],n_left,mu_left,se_left,n_right,mu_right,se_right);
		  //cout << n_left << "\t" << n_right << "\t" << se_left << "\t" << se_right << endl;
		}

	      if(se_left+se_right < se_best)
		{
		  it_best = it;
		  se_best = se_left+se_right;
		}

	      //Take samples from left back to right
	      for(size_t i = 0; i < it->second.size(); ++i)
                {
		  //cout << tv[it->second[i]] << ": ";
		  datadefs::forward_backward_sqerr(tv[it->second[i]],n_right,mu_right,se_right,n_left,mu_left,se_left);
		  //cout << n_left << "\t" << n_right << "\t" << se_left << "\t" << se_right << endl;
                }
	    }
	  
	  //If we found a categorical split that reduces impurity...
	  if(it_best == fmap_right.end())
	    {
	      break;
	    }


	  //Take samples from right and put them left
	  for(size_t i = 0; i < it_best->second.size(); ++i)
	    {
	      categories_left.insert(it_best->first);
	      //cout << "removing index " << it_best->second[i] << " (value " << tv[it_best->second[i]] << ") from right: ";
	      datadefs::forward_backward_sqerr(tv[it_best->second[i]],n_left,mu_left,se_left,n_right,mu_right,se_right);
	      sampleics_left.push_back(sampleics[it_best->second[i]]);
	      //cout << n_left << "\t" << n_right << "\t" << se_left << "\t" << se_right << endl;
	    }
	  fmap_right.erase(it_best->first);
	  
	}
    }  
  else
    {
      map<num_t,size_t> freq_left,freq_right;
      size_t sf_left = 0;
      size_t sf_right;
      datadefs::sqfreq(tv,freq_right,sf_right,n_right);
      assert(n_tot == n_right);
      //datadefs::count_freq(fv,freq_right,n_right);
      //for(map<num_t,size_t>::const_iterator it(freq_right.begin()); it != freq_right.end(); ++it)
      //	{
      //	  cout << " " << it->first << ":" << it->second;
      //	}
      //cout << endl;
      num_t nsf_best = 1.0 * sf_right / n_right;

      while(fmap_right.size() > 1)
        {

          map<num_t,vector<size_t> >::iterator it_best(fmap_right.end());

          for(map<num_t,vector<size_t> >::iterator it(fmap_right.begin()); it != fmap_right.end(); ++it)
            {

	      size_t n_left_c = n_left;
	      size_t n_right_c = n_right;
	      map<num_t,size_t> freq_left_c = freq_left;
	      map<num_t,size_t> freq_right_c = freq_right;

              //Take samples from right and put them left
	      //cout << "Moving " << it->second.size() << " samples corresponding to feature category " << it->first << " from right to left" << endl;
              for(size_t i = 0; i < it->second.size(); ++i)
                {
		  //cout << " " << tv[it->second[i]];
		  datadefs::forward_backward_sqfreq(tv[it->second[i]],n_left,freq_left,sf_left,n_right,freq_right,sf_right);
		  //cout << n_left << "\t" << n_right << "\t" << sf_left << "\t" << sf_right << endl;
                }
	      //cout << endl;

              if(1.0*n_right*sf_left + n_left*sf_right > n_left*n_right*nsf_best)
                {
                  it_best = it;
                  nsf_best = 1.0*sf_left/n_left + 1.0*sf_right/n_right;
                }

	      //cout << "Moving " << it->second.size() << " samples corresponding to feature category " << it->first << " from left to right" << endl;
              //Take samples from left back to right
              for(size_t i = 0; i < it->second.size(); ++i)
                {
		  //cout << " " << tv[it->second[i]];
		  datadefs::forward_backward_sqfreq(tv[it->second[i]],n_right,freq_right,sf_right,n_left,freq_left,sf_left);
		  //cout << n_left << "\t" << n_right << "\t" << sf_left << "\t" << sf_right << endl;
                }
	      //cout << endl;

	      assert(n_left_c == n_left);
	      assert(n_right_c == n_right);
	      assert(freq_left_c.size() == freq_left.size());
	      assert(freq_right_c.size() == freq_right.size());

            }

          //If we found a categorical split that reduces impurity...
          if(it_best == fmap_right.end())
            {
              break;
            }


          //Take samples from right and put them left
          for(size_t i = 0; i < it_best->second.size(); ++i)
            {
	      categories_left.insert(it_best->first);
	      //cout << "removing index " << it_best->second[i] << " (value " << tv[it_best->second[i]] << ") from right: ";
	      datadefs::forward_backward_sqfreq(tv[it_best->second[i]],n_left,freq_left,sf_left,n_right,freq_right,sf_right);
              sampleics_left.push_back(sampleics[it_best->second[i]]);
	      //cout << n_left << "\t" << n_right << "\t" << sf_left << "\t" << sf_right << endl;
            }
          fmap_right.erase(it_best->first);

        }
    }

  //Next we assign samples remaining on right
  sampleics_right.clear();
  for(map<num_t,vector<size_t> >::const_iterator it(fmap_right.begin()); it != fmap_right.end(); ++it)
    {
      for(size_t i = 0; i < it->second.size(); ++i)
	{
	  sampleics_right.push_back(sampleics[it->second[i]]);
	}
    }  

  if(false)
    {
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
}

num_t Treedata::at(size_t featureidx, size_t sampleidx)
{
  return(featurematrix_[featureidx][sampleidx]);
}



//This is slow if the feature has lots of missing values. Consider re-implementing for speed-up
num_t Treedata::randf(size_t featureidx)
{
  assert(nrealsamples_ > 0);
  num_t value(datadefs::num_nan);
  while(datadefs::is_nan(value))
    {
      value = featurematrix_[featureidx][irand_() % nrealsamples_];
    }
  return(value);
}



void Treedata::impurity(size_t featureidx, vector<size_t> const& sampleics, num_t& impurity, size_t& nreal)
{
  
  size_t n(sampleics.size());
  //vector<num_t> data(n);
  
  //for(size_t i = 0; i < n; ++i)
  //  {
  //    data[i] = featurematrix_[featureidx][sampleics[i]];
  //  }
  
  nreal = 0;
  if(isfeaturenum_[featureidx])
    {
      num_t mu = 0.0;
      num_t se = 0.0;
      for(size_t i = 0; i < n; ++i)
	{
	  datadefs::forward_sqerr(featurematrix_[featureidx][sampleics[i]],nreal,mu,se);  
	}
      impurity = se / nreal;
    }
  else
    {
      map<num_t,size_t> freq;
      size_t sf = 0;
      for(size_t i = 0; i < n; ++i)
	{
	  datadefs::forward_sqfreq(featurematrix_[featureidx][sampleics[i]],nreal,freq,sf);
	}
      impurity = 1.0 - 1.0 * sf / (nreal * nreal);
    }
}


void Treedata::split_samples(vector<size_t>& sampleics, 
			     size_t splitidx, 
			     vector<size_t>& sampleics_left, 
			     vector<size_t>& sampleics_right)
{
  //Split the sample indices
  size_t n_tot = sampleics.size();
  size_t n_left = splitidx + 1;
  size_t n_right = n_tot - n_left;

  if(n_left == 0)
    {
      sampleics_right = sampleics;
      sampleics_left.clear();
      return;
    }
  else if(n_right == 0)
    {
      sampleics_right.clear();
      sampleics_left = sampleics;
      return;
    }

  //vector<size_t> sampleics_left_copy(n_left);
  //vector<size_t> sampleics_right_copy(n_right);
  sampleics_left.resize(n_left);
  sampleics_right.resize(n_right);
  size_t i = 0;
  while(i < n_left)
    {
      sampleics_left[i] = sampleics[i];
      ++i;
    }
  while(i < n_tot)
    {
      sampleics_right[i-n_left] = sampleics[i];
      ++i;
    }
}

num_t Treedata::split_fitness(const size_t featureidx,
			      const size_t min_split,
			      vector<size_t> const& sampleics,
			      vector<size_t> const& sampleics_left,
			      vector<size_t> const& sampleics_right)
{
  
  assert(sampleics.size() == sampleics_left.size() + sampleics_right.size());

  //size_t n_tot = sampleics.size();
  size_t n_left = 0;
  size_t n_right = 0;
  if(isfeaturenum_[featureidx])
    {
      num_t mu_left = 0.0;
      num_t se_left = 0.0;
      num_t mu_right = 0.0;
      num_t se_right = 0.0;
 
      for(size_t i = 0; i < sampleics_left.size(); ++i)
	{
	  datadefs::forward_sqerr(featurematrix_[featureidx][sampleics_left[i]],n_right,mu_right,se_right);
	  //cout << "forward sqerr: " << featurematrix_[featureidx][sampleics_left[i]] << " " << n_right << " " << mu_right << " " << se_right << endl; 
	}

      for(size_t i = 0; i < sampleics_right.size(); ++i)
        {
	  datadefs::forward_sqerr(featurematrix_[featureidx][sampleics_right[i]],n_right,mu_right,se_right);
	  //cout << "forward sqerr: " << featurematrix_[featureidx][sampleics_right[i]] << " " << n_right << " " << mu_right << " " << se_right << endl;
        }

      if(n_right < 2*min_split)
        {
          return(0.0);
        }

      num_t se_tot = se_right;
      
      for(size_t i = 0; i < sampleics_left.size(); ++i)
	{
	  datadefs::forward_backward_sqerr(featurematrix_[featureidx][sampleics_left[i]],n_left,mu_left,se_left,n_right,mu_right,se_right);
	  //cout << "fw bw sqerr: " << featurematrix_[featureidx][sampleics_left[i]] << " " << n_left << " " << mu_left << " " << se_left << " " << n_right << " " << mu_right << " " << se_right << endl;
	}

      if(n_left < min_split || n_right < min_split)
	{
      	  return(0.0);
	}
      
      return(( se_tot - se_left - se_right ) / se_tot);
      
    }
  else
    {
      map<num_t,size_t> freq_left,freq_right;
      size_t sf_left = 0;
      size_t sf_right = 0;

      for(size_t i = 0; i < sampleics_left.size(); ++i)
        {
	  datadefs::forward_sqfreq(featurematrix_[featureidx][sampleics_left[i]],n_right,freq_right,sf_right);
	  //cout << "forward sqfreq: " << featurematrix_[featureidx][sampleics_left[i]] << " " << n_right << " " << sf_right << endl;
        }

      for(size_t i = 0; i < sampleics_right.size(); ++i)
        {
	  datadefs::forward_sqfreq(featurematrix_[featureidx][sampleics_right[i]],n_right,freq_right,sf_right);
	  //cout << "forward sqfreq: " << featurematrix_[featureidx][sampleics_right[i]] << " " << n_right << " " << sf_right << endl;
        }

      if(n_right < 2*min_split)
       {
         return(0.0);
       }

      size_t n_tot = n_right;
      size_t sf_tot = sf_right;

      for(size_t i = 0; i < sampleics_left.size(); ++i)
        {
	  datadefs::forward_backward_sqfreq(featurematrix_[featureidx][sampleics_left[i]],n_left,freq_left,sf_left,n_right,freq_right,sf_right);
	  //cout << "fw bw sqfreq: " << featurematrix_[featureidx][sampleics_left[i]] << " " << n_left << " "<< sf_left << " " << n_right << " " << sf_right << endl;
        }

      if(n_left < min_split || n_right < min_split)
        {
          return(0.0);
        }

      //cout << n_left << " " << n_right << " " << sf_tot << " " << n_tot << " " << n_right << " " << sf_left << " " << n_tot*n_left*sf_right << " " << n_left*n_right << " " << pow(n_tot,2) - sf_tot << endl;

      //num_t fitness = (-1.0*(n_left*n_right*sf_tot) + n_tot*n_right*sf_left + n_tot*n_left*sf_right) / (n_left*n_right*(pow(n_tot,2) - sf_tot));
      //cout << "Fitness " << fitness << endl;

      return( ( -1.0 * n_left*n_right*sf_tot + n_tot*n_right*sf_left + n_tot*n_left*sf_right ) / ( n_left*n_right * (1.0*n_tot*n_tot - sf_tot) ) ) ;
      
    }

}


