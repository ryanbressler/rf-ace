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

Treedata::Treedata(string fileName):
  features_(0),
  nSamples_(0),
  nFeatures_(0)
{

  //Initialize random number rgenerator
  time_t now;
  time(&now);
  //srand((unsigned int)now);
  
  //MTRand_int32 irand((unsigned int)now);
  randomInteger_.seed((unsigned int)now);
  //datadefs::seedMT((size_t)now);

  //cout << "Treedata: reading matrix from file '" << filename << "'" << endl;

  //Initialize stream to read from file
  ifstream featurestream;
  featurestream.open(fileName.c_str());
  assert(featurestream.good());

  //TODO: a sniffer function that find out the type of the input file based on its content.
  //*************************************************************************************** 
  FileType fileType = UNKNOWN;
  Treedata::readFileType(fileName,fileType);

  vector<vector<string> > rawMatrix;
  vector<string> featureHeaders;
  vector<bool> isFeatureNumerical;
  if(fileType == AFM)
    {
      cout << "File type interpreted as Annotated Feature Matrix (AFM)" << endl;
      try {
	Treedata::readAFM(featurestream,rawMatrix,featureHeaders,isFeatureNumerical);
      }
      catch(...)
	{
	  cerr << "The file could not be read. Is the file type really AFM? Quitting..." << endl;
	}
    }
  else if(fileType == ARFF)
    {
      cout << "File type interpreted as Attribute-Relation File Format (ARFF)" << endl;
      try {
	Treedata::readARFF(featurestream,rawMatrix,featureHeaders,isFeatureNumerical);
      }
      catch(...)
	{
	  cerr << "The file could not be read. Is the file type really ARFF? Quitting..." << endl;
	}
    }
  else
    {
      cout << "File type is unknown -- defaulting to Annotated Feature Matrix (AFM)" << endl;
      try
	{
	  Treedata::readAFM(featurestream,rawMatrix,featureHeaders,isFeatureNumerical);
	}
      catch(...)
	{
	  cerr << "The file could not be read, quitting..." << endl;
	}
    }      

  nSamples_ = rawMatrix[0].size();
  nFeatures_ = featureHeaders.size();
  features_.resize(2*nFeatures_);

  for(size_t i = 0; i < nFeatures_; ++i)
    {
      vector<num_t> featureData(nSamples_);
      features_[i].name = featureHeaders[i];
      features_[i].isNumerical = isFeatureNumerical[i];
      if(features_[i].isNumerical)
	{
          datadefs::strv2numv(rawMatrix[i],featureData);
	  features_[i].nCategories = 0;
        }
      else
        {
          datadefs::strv2catv(rawMatrix[i],featureData);
	  map<num_t,size_t> freq;
	  size_t nReal;
	  datadefs::count_freq(featureData,freq,nReal);
	  features_[i].nCategories = freq.size();
        }
      features_[i].data = featureData;
      //Treedata::permute(featureData);
      //features_[i].contrast = featureData;
    } 
  
  for(size_t i = nFeatures_; i < 2*nFeatures_; ++i)
    {
      //vector<num_t> contrastData(nSamples_);
      features_[i].name = "CONTRAST";
      features_[i].isNumerical = features_[ i-nFeatures_ ].isNumerical;
      features_[i].data = features_[ i-nFeatures_ ].data;
    }

  Treedata::permuteContrasts();

  //targetidx_ = 0;
  //Treedata::selectTarget(targetidx_);

}

Treedata::~Treedata()
{
}

void Treedata::readFileType(string& fileName, FileType& fileType)
{

  stringstream ss(fileName);
  string suffix = "";
  while(getline(ss,suffix,'.')) {}
  //datadefs::toupper(suffix);

  if(suffix == "AFM" || suffix == "afm")
    {
      fileType = AFM;
    }
  else if(suffix == "ARFF" || suffix == "arff")
    {
      fileType = ARFF;
    }
  else
    {
      fileType = UNKNOWN;
    }

}

void Treedata::readAFM(ifstream& featurestream, vector<vector<string> >& rawMatrix, vector<string>& featureHeaders, vector<bool>& isFeatureNumerical)
{

  string field;
  string row;

  //TODO: add Treedata::clearData(...)
  rawMatrix.clear();
  featureHeaders.clear();
  isFeatureNumerical.clear();

  //Remove upper left element from the matrix as useless
  getline(featurestream,field,'\t');

  //Count the number of columns (NEEDS REWORKING)
  getline(featurestream,row);
  stringstream ss(row);
  size_t nColumns = 0;
  bool isFeaturesAsRows = true;
  while(getline(ss,field,'\t'))
    {
      if(isFeaturesAsRows && isValidFeatureHeader(field))
	{
	  isFeaturesAsRows = false;
	}
      ++nColumns;
    }

  //Count the number of rows
  size_t nRows = 0;
  while(getline(featurestream,row))
    {
      ++nRows;
    }

  //Reset streams and remove upper left element from the matrix as useless
  featurestream.clear();
  featurestream.seekg(0, ios::beg);
  getline(featurestream,field,'\t');
  ss.clear();
  ss.str("");

  //These are temporary containers
  vector<string> columnHeaders(nColumns);
  vector<string> rowHeaders(nRows);
  //vector<vector<string> > rawmatrix(nrows);
  rawMatrix.resize(nRows);
  for(size_t i = 0; i < nRows; ++i)
    {
      rawMatrix[i] = columnHeaders;
    }

  //cout << "read " << rawmatrix.size() << " rows and " << rawmatrix[0].size() << " columns." << endl;

  //Read first row into the stream
  getline(featurestream,row);
  ss << row;

  //Read column headers from the stream
  for(size_t i = 0; i < nColumns; ++i)
    {
      getline(ss,columnHeaders[i],'\t');
      //cout << '\t' << colheaders[i];
    }
  //cout << endl;

  //Go through the rest of the rows
  for(size_t i = 0; i < nRows; ++i)
    {
      //Read row from the stream
      getline(featurestream,row);
      ss.clear();
      ss.str("");

      //Read the string back to a stream (REDUNDANCY)
      ss << row;

      //Read one element from the row stream
      getline(ss,rowHeaders[i],'\t');
      //cout << rowheaders[i];
      for(size_t j = 0; j < nColumns; ++j)
        {
          getline(ss,rawMatrix[i][j],'\t');
          //cout << '\t' << rawmatrix[i][j];
        }
      //cout << endl;
    }
  //cout << endl;

  //If the data is row-formatted...
  if(isFeaturesAsRows)
    {
      cout << "AFM orientation: features as rows" << endl;

      //... and feature headers are row headers
      featureHeaders = rowHeaders;

    }
  else
    {

      cout << "AFM orientation: features as columns" << endl;
      
      Treedata::transpose<string>(rawMatrix);
      
      //... and feature headers are row headers
      featureHeaders = columnHeaders;
      
    }

  size_t nFeatures = featureHeaders.size();
  isFeatureNumerical.resize(nFeatures);
  for(size_t i = 0; i < nFeatures; ++i)
    {
      if(Treedata::isValidNumericalHeader(featureHeaders[i]))
	{
	  isFeatureNumerical[i] = true;
	}
      else
	{
	  isFeatureNumerical[i] = false;
	}
    }
}

void Treedata::readARFF(ifstream& featurestream, vector<vector<string> >& rawMatrix, vector<string>& featureHeaders, vector<bool>& isFeatureNumerical)
{

  string row;

  bool hasRelation = false;
  bool hasData = false;

  size_t nFeatures = 0;
  //TODO: add Treedata::clearData(...)
  rawMatrix.clear();
  featureHeaders.clear();
  isFeatureNumerical.clear();
  
  //Read one line from the ARFF file
  while(getline(featurestream,row))
    {

      //This is the final branch: once relation and attributes are read, and we find data header, we'll start reading the data in 
      if(hasData && hasRelation)
	{
	  //There must be at least two attributes, otherwise the ARFF file makes no sense
	  if(nFeatures < 2)
	    {
	      cerr << "too few attributes ( < 2 ) found from the ARFF file" << endl;
	      assert(false);
	    }

	  rawMatrix.resize(nFeatures);

	  //Read data row-by-row
	  while(true)
	    {
	      //++nsamples_;
	      string field;
	      stringstream ss(row);
	      
	      for(size_t attributeIdx = 0; attributeIdx < nFeatures; ++attributeIdx)
		{
		  getline(ss,field,',');
		  //cout << " " << field;
		  rawMatrix[attributeIdx].push_back(field);
		}
	      //cout << endl;

	      if(!getline(featurestream,row))
		{
		  break;
		}
	    }

	  //sampleheaders_.resize(nsamples_);
	  //for(size_t i = 0; i < nsamples_; ++i)
	  //  {
	  //    sampleheaders_[i] = "foo";
	  //  }
	  //We're done, exit
	  break;
	}

      //Comment lines and empty lines are omitted
      if(row[0] == '%' || row == "")
	{
	  continue;	
	}

      //Read relation
      if(!hasRelation && row.compare(0,9,"@relation") == 0)
	{
	  hasRelation = true;
	  //cout << "found relation header: " << row << endl;
	}
      //Read attribute
      else if(row.compare(0,10,"@attribute") == 0)
	{
	  string attributeName = "";
	  bool isNumerical;
	  ++nFeatures;
	  //cout << "found attribute header: " << row << endl;
	  Treedata::parseARFFattribute(row,attributeName,isNumerical);
	  featureHeaders.push_back(attributeName);
	  isFeatureNumerical.push_back(isNumerical);
	  //isfeaturenum_.push_back(isNumeric);
	  //cout << "interpreted as: " << attributeName << " (";
	  //if(isNumeric)
	  // {
	  //   cout << "numeric)" << endl; 
	  //  }
	  //else
	  //  {
	  //    cout << "categorical)" << endl;
	  //  }
	}
      //Read data header
      else if(!hasData && row.compare(0,5,"@data") == 0)
	{
	  hasData = true;
	  //cout << "found data header:" << row << endl;
	}
      //If none of the earlier branches matched, we have a problem
      else
	{
	  cerr << "incorrectly formatted ARFF file" << endl;
	  assert(false);
	}
    }
}

void Treedata::parseARFFattribute(const string& str, string& attributeName, bool& isFeatureNumerical)
{

  stringstream ss(str);
  string attributeHeader = "";
  attributeName = "";
  string attributeType = "";

  getline(ss,attributeHeader,' ');
  getline(ss,attributeName,' ');
  getline(ss,attributeType);

  //string prefix;
  if(attributeType == "real")
    {
      isFeatureNumerical = true;
    }
  else
    {
      isFeatureNumerical = false;
    }
  //prefix.append(attributeName);
  //attributeName = prefix;
}

bool Treedata::isValidNumericalHeader(const string& str)
{
  return( str.compare(0,2,"N:") == 0 );
}

bool Treedata::isValidCategoricalHeader(const string& str)
{
  return( str.compare(0,2,"C:") == 0 || str.compare(0,2,"B:") == 0 );
}

bool Treedata::isValidFeatureHeader(const string& str)
{
  return( isValidNumericalHeader(str) || isValidCategoricalHeader(str) );
}

size_t Treedata::nFeatures()
{
  return(nFeatures_);
}

size_t Treedata::nSamples()
{
  return(nSamples_);
}

//WILL BECOME DEPRECATED
num_t Treedata::pearsonCorrelation(size_t featureidx1, size_t featureidx2)
{
  num_t r;
  datadefs::pearson_correlation(features_[featureidx1].data,features_[featureidx2].data,r);
  return(r);
}

/*
  void Treedata::killFeature(const size_t featureIdx)
  {
  assert(featureIdx < nFeatures_);
  
  features_.erase(features_.begin() + featureIdx);
  
  --nFeatures_;
  
  }
*/

string Treedata::getFeatureName(const size_t featureIdx)
{
  return(features_[featureIdx].name);
}


void Treedata::print()
{
  cout << "Printing feature matrix (missing values encoded to " << datadefs::NUM_NAN << "):" << endl;
  for(size_t j = 0; j < nSamples_; ++j)
    {
      cout << '\t' << "foo";
    }
  cout << endl;
  for(size_t i = 0; i < nFeatures_; ++i)
    {
      cout << i << ':' << features_[i].name << ':';
      for(size_t j = 0; j < nSamples_; ++j)
        {
          cout << '\t' << features_[i].data[j];
        }
      cout << endl;
    }
}


void Treedata::print(const size_t featureIdx)
{
  cout << "Print " << features_[featureIdx].name << ":";
  for(size_t i = 0; i < nSamples_; ++i)
    {
      cout << " " << features_[featureIdx].data[i];
    }
  cout << endl;
}


void Treedata::permuteContrasts()
{
  for(size_t i = nFeatures_; i < 2*nFeatures_; ++i)
    {
      Treedata::permute(features_[i].data);
    }
}

bool Treedata::isFeatureNumerical(size_t featureIdx)
{
  return(features_[featureIdx].isNumerical);
}


size_t Treedata::nRealSamples(const size_t featureIdx)
{ 
  size_t nRealSamples;
  datadefs::countRealValues(features_[featureIdx].data,nRealSamples);
  return(nRealSamples);
}

size_t Treedata::nCategories(const size_t featureIdx)
{
  return(features_[featureIdx].nCategories);
}


//WILL BE DEPRECATED
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


void Treedata::permute(vector<size_t>& ics)
{
  for (size_t i = 0; i < ics.size(); ++i)
    {
      size_t j = randomInteger_() % (i + 1);
      ics[i] = ics[j];
      ics[j] = i;
    }
}

void Treedata::permute(vector<num_t>& data)
{
  size_t n = data.size();
  vector<size_t> ics(n);

  Treedata::permute(ics);

  for(size_t i = 0; i < n; ++i)
    {
      num_t temp = data[i];
      data[i] = data[ics[i]];
      data[ics[i]] = temp;
    }
}


void Treedata::bootstrapFromRealSamples(const size_t featureIdx, vector<size_t>& ics, vector<size_t>& oobIcs)
{
    
  ics.clear();
  for(size_t i = 0; i < nSamples_; ++i)
    {
      if(!datadefs::isNAN(features_[featureIdx].data[i]))
	{
	  ics.push_back(i);
	}
    }
  
  size_t n = ics.size();
  vector<size_t> allIcs = ics;
  for(size_t i = 0; i < n; ++i)
    {
      ics[i] = allIcs[randomInteger_() % n];
    }

  vector<size_t> foo(n);
  vector<size_t>::iterator it = set_difference(allIcs.begin(),allIcs.end(),ics.begin(),ics.end(),foo.begin());
  size_t nOob = distance(foo.begin(),it);

  foo.resize(nOob);
  oobIcs = foo;
}


void Treedata::numericalFeatureSplit(vector<num_t>& tv,
				     const bool isTargetNumerical,
				     vector<num_t>& fv,
				     const size_t min_split, 
				     vector<size_t>& sampleIcs_left, 
				     vector<size_t>& sampleIcs_right,
				     num_t& splitValue)
{

  assert(tv.size() == fv.size());

  size_t n_tot = tv.size();
  size_t n_right = n_tot;
  size_t n_left = 0;

  sampleIcs_left.clear();
  sampleIcs_right.clear();
  
  //Check that there are enough samples to make the split in the first place
  assert(n_tot >= 2*min_split);
  
  //Make reference indices that define the sorting wrt. feature
  vector<size_t> refIcs;
  
  //Sort feature vector and collect reference indices
  datadefs::sortDataAndMakeRef(fv,sampleIcs_right);
  
  //Use the reference indices to sort sample indices
  //datadefs::sortFromRef<size_t>(sampleics,refIcs);
  datadefs::sortFromRef<num_t>(tv,sampleIcs_right);
  
  //Count how many real values the feature and target has
  size_t nreal_f,nreal_t;
  datadefs::countRealValues(fv,nreal_f);
  datadefs::countRealValues(tv,nreal_t);
  //cout << nreal_f << " " << nreal_t << " " << n_tot << endl;
  assert(nreal_t == n_tot && nreal_f == n_tot);

  int bestSplitIdx = -1;
 
  //If the target is numerical, we use the iterative squared error formula to update impurity scores while we traverse "right"
  if(isTargetNumerical)
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
	      bestSplitIdx = idx;
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
              bestSplitIdx = idx;
              nsf_best = 1.0 * sf_left / n_left + 1.0 * sf_right / n_right;
            }
	  ++idx;
        }      
    }
  
  splitValue = fv[bestSplitIdx];
  n_left = bestSplitIdx + 1;
  sampleIcs_left.resize(n_left);

  for(size_t i = 0; i < n_left; ++i)
    {
      sampleIcs_left[i] = sampleIcs_right[i];
    }
  sampleIcs_right.erase(sampleIcs_right.begin(),sampleIcs_right.begin() + n_left);

  //cout << sampleIcs_left.size() << " " << sampleIcs_right.size() << " == " << n_tot << endl; 
  assert(sampleIcs_left.size() + sampleIcs_right.size() == n_tot);

  if(false)
    {
      cout << "Feature splits target [";
      for(size_t i = 0; i < sampleIcs_left.size(); ++i)
	{
	  cout << " " << tv[sampleIcs_left[i]];
	}
      cout << " ] <==> [";
      for(size_t i = 0; i < sampleIcs_right.size(); ++i)
	{
	  cout << " " << tv[sampleIcs_right[i]];
	}
      cout << " ]" << endl;
    }
}

void Treedata::categoricalFeatureSplit(vector<num_t>& tv,
				       const bool isTargetNumerical,
				       vector<num_t>& fv,
				       vector<size_t>& sampleIcs_left,
				       vector<size_t>& sampleIcs_right,
				       set<num_t>& categories_left)
{

  assert(tv.size() == fv.size());

  size_t n_tot = tv.size();
  size_t n_right = n_tot;
  size_t n_left = 0;

  sampleIcs_left.clear();
  sampleIcs_right.clear();
  categories_left.clear();
  
  //Check that sample size is positive
  assert(n_tot > 0);

  //Count how many real values the feature and target has (this is just to make sure there are no missing values thrown in)
  size_t nreal_f,nreal_t;
  datadefs::countRealValues(fv,nreal_f);
  datadefs::countRealValues(tv,nreal_t);
  assert(nreal_t == n_tot && nreal_f == n_tot);

  //Map all feature categories to the corresponding samples and represent it as map. The map is used to assign samples to left and right branches
  map<num_t,vector<size_t> > fmap_right;
  datadefs::map_data(fv,fmap_right,nreal_f);
  assert(n_tot == nreal_f);

  if(isTargetNumerical)
    {
      
      num_t mu_right;
      num_t mu_left = 0.0;
      num_t se_right;
      num_t se_left = 0.0;
      datadefs::sqerr(tv,mu_right,se_right,n_right);
      assert(n_tot == n_right);
      num_t se_best = se_right;
      
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
		  se_best = se_left + se_right;
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
	      sampleIcs_left.push_back(it_best->second[i]);
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
              sampleIcs_left.push_back(it_best->second[i]);
	      //cout << n_left << "\t" << n_right << "\t" << sf_left << "\t" << sf_right << endl;
            }
          fmap_right.erase(it_best->first);

        }
    }

  //Next we assign samples remaining on right
  sampleIcs_right.clear();
  for(map<num_t,vector<size_t> >::const_iterator it(fmap_right.begin()); it != fmap_right.end(); ++it)
    {
      for(size_t i = 0; i < it->second.size(); ++i)
	{
	  sampleIcs_right.push_back(it->second[i]);
	}
    }  

  if(false)
    {
      cout << "Feature splits target [";
      for(size_t i = 0; i < sampleIcs_left.size(); ++i)
	{
	  cout << " " << tv[sampleIcs_left[i]];
	}
      cout << " ] <==> [";
      for(size_t i = 0; i < sampleIcs_right.size(); ++i)
	{
	  cout << " " << tv[sampleIcs_right[i]];
	}
      cout << " ]" << endl;
    }
}

void Treedata::getFeatureData(size_t featureIdx, vector<num_t>& data)
{
  data.resize(nSamples_);
  for(size_t i = 0; i < nSamples_; ++i)
    {
      data[i] = features_[featureIdx].data[i];
    }
}

void Treedata::getFeatureData(size_t featureIdx, const size_t sampleIdx, num_t& data)
{
  data = features_[featureIdx].data[sampleIdx];
}

void Treedata::getFeatureData(size_t featureIdx, const vector<size_t>& sampleIcs, vector<num_t>& data)
{
  data.resize(sampleIcs.size());
  for(size_t i = 0; i < sampleIcs.size(); ++i)
    {
      data[i] = features_[featureIdx].data[sampleIcs[i]];
    }

}


/*
  void Treedata::getContrastData(size_t featureIdx, vector<num_t>& data)
  {
  data.resize(nSamples_);
  for(size_t i = 0; i < nSamples_; ++i)
  {
  data[i] = features_[featureIdx].contrast[i];
  }
  }
  
  void Treedata::getContrastData(size_t featureIdx, const size_t sampleIdx, num_t& data)
  {
  data = features_[featureIdx].contrast[sampleIdx];
  }
  
  void Treedata::getContrastData(size_t featureIdx, const vector<size_t>& sampleIcs, vector<num_t>& data)
  {
  data.resize(sampleIcs.size());
  for(size_t i = 0; i < sampleIcs.size(); ++i)
  {
  data[i] = features_[featureIdx].contrast[sampleIcs[i]];
  }
  
  }
*/



/*
  void Treedata::getRandomData(const size_t featureIdx, num_t& data)
  {
  data = features_[featureIdx].data[irand_() % nSamples_];
  }
*/

/*
  void Treedata::getRandomData(const size_t featureIdx, num_t& data)
  {
  
  num_t featureSample = datadefs::NUM_NAN;
  num_t targetSample = datadefs::NUM_NAN;
  while(datadefs::isNAN(featureSample) || datadefs::isNAN(targetSample))
  {
  size_t randomIdx = irand_() % nsamples_;
  featureSample = features_[featureIdx].data[randomIdx];
  targetSample = features_[targetIdx].data[randomIdx];
  }
  data = featureSample;
  }
*/

/*
  This is slow if the feature has lots of missing values. Consider re-implementing for speed-up
  num_t Treedata::randf(size_t featureidx)
  {
  assert(nrealsamples_ > 0);
  num_t value = datadefs::NUM_NAN;
  while(datadefs::isNAN(value))
  {
  size_t randomIdx = irand_() % nsamples_;
  
  //value = features_[featureidx].data[irand_() % nrealsamples_];
  }
  return(value);
  }
*/

void Treedata::impurity(vector<num_t>& data, bool isFeatureNumerical, num_t& impurity, size_t& nreal)
{
  
  size_t n = data.size();
  
  
  nreal = 0;
  if(isFeatureNumerical)
    {
      num_t mu = 0.0;
      num_t se = 0.0;
      for(size_t i = 0; i < n; ++i)
	{
	  datadefs::forward_sqerr(data[i],nreal,mu,se);  
	}
      impurity = se / nreal;
    }
  else
    {
      map<num_t,size_t> freq;
      size_t sf = 0;
      for(size_t i = 0; i < n; ++i)
	{
	  datadefs::forward_sqfreq(data[i],nreal,freq,sf);
	}
      impurity = 1.0 - 1.0 * sf / (nreal * nreal);
    }
}

/*
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
*/

num_t Treedata::splitFitness(vector<num_t> const& data,
			     bool const& isFeatureNumerical,
			     size_t const& minSplit,
			     vector<size_t> const& sampleIcs_left,
			     vector<size_t> const& sampleIcs_right)
{
  
  assert(data.size() == sampleIcs_left.size() + sampleIcs_right.size());

  size_t n_left = 0;
  size_t n_right = 0;
  if(isFeatureNumerical)
    {
      num_t mu_left = 0.0;
      num_t se_left = 0.0;
      num_t mu_right = 0.0;
      num_t se_right = 0.0;
 
      for(size_t i = 0; i < sampleIcs_left.size(); ++i)
	{
	  datadefs::forward_sqerr(data[sampleIcs_left[i]],n_right,mu_right,se_right);
	  //cout << "forward sqerr: " << featurematrix_[featureidx][sampleics_left[i]] << " " << n_right << " " << mu_right << " " << se_right << endl; 
	}

      for(size_t i = 0; i < sampleIcs_right.size(); ++i)
        {
	  datadefs::forward_sqerr(data[sampleIcs_right[i]],n_right,mu_right,se_right);
	  //cout << "forward sqerr: " << featurematrix_[featureidx][sampleics_right[i]] << " " << n_right << " " << mu_right << " " << se_right << endl;
        }

      if(n_right < 2*minSplit)
        {
          return(0.0);
        }

      num_t se_tot = se_right;
      
      for(size_t i = 0; i < sampleIcs_left.size(); ++i)
	{
	  datadefs::forward_backward_sqerr(data[sampleIcs_left[i]],n_left,mu_left,se_left,n_right,mu_right,se_right);
	  //cout << "fw bw sqerr: " << featurematrix_[featureidx][sampleics_left[i]] << " " << n_left << " " << mu_left << " " << se_left << " " << n_right << " " << mu_right << " " << se_right << endl;
	}

      if(n_left < minSplit || n_right < minSplit)
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

      for(size_t i = 0; i < sampleIcs_left.size(); ++i)
        {
	  datadefs::forward_sqfreq(data[sampleIcs_left[i]],n_right,freq_right,sf_right);
	  //cout << "forward sqfreq: " << featurematrix_[featureidx][sampleics_left[i]] << " " << n_right << " " << sf_right << endl;
        }

      for(size_t i = 0; i < sampleIcs_right.size(); ++i)
        {
	  datadefs::forward_sqfreq(data[sampleIcs_right[i]],n_right,freq_right,sf_right);
	  //cout << "forward sqfreq: " << featurematrix_[featureidx][sampleics_right[i]] << " " << n_right << " " << sf_right << endl;
        }

      if(n_right < 2*minSplit)
       {
         return(0.0);
       }

      size_t n_tot = n_right;
      size_t sf_tot = sf_right;

      for(size_t i = 0; i < sampleIcs_left.size(); ++i)
        {
	  datadefs::forward_backward_sqfreq(data[sampleIcs_left[i]],n_left,freq_left,sf_left,n_right,freq_right,sf_right);
	  //cout << "fw bw sqfreq: " << featurematrix_[featureidx][sampleics_left[i]] << " " << n_left << " "<< sf_left << " " << n_right << " " << sf_right << endl;
        }

      if(n_left < minSplit || n_right < minSplit)
        {
          return(0.0);
        }

      //cout << n_left << " " << n_right << " " << sf_tot << " " << n_tot << " " << n_right << " " << sf_left << " " << n_tot*n_left*sf_right << " " << n_left*n_right << " " << pow(n_tot,2) - sf_tot << endl;

      //num_t fitness = (-1.0*(n_left*n_right*sf_tot) + n_tot*n_right*sf_left + n_tot*n_left*sf_right) / (n_left*n_right*(pow(n_tot,2) - sf_tot));
      //cout << "Fitness " << fitness << endl;

      return( ( -1.0 * n_left*n_right*sf_tot + n_tot*n_right*sf_left + n_tot*n_left*sf_right ) / ( n_left*n_right * (1.0*n_tot*n_tot - sf_tot) ) ) ;
      
    }

}


