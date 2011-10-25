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
  nFeatures_(0) {

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
  if ( !featurestream.good() ) {
    cerr << "Failed to open file '" << fileName << "' for reading" << endl;
    assert(false);
  }

  //TODO: a sniffer function that find out the type of the input file based on its content.
  //*************************************************************************************** 
  FileType fileType = UNKNOWN;
  Treedata::readFileType(fileName,fileType);

  vector<vector<string> > rawMatrix;
  vector<string> featureHeaders;
  vector<bool> isFeatureNumerical;
  if(fileType == AFM) {
    //cout << "File type interpreted as Annotated Feature Matrix (AFM)" << endl;
    try {
      Treedata::readAFM(featurestream,rawMatrix,featureHeaders,isFeatureNumerical);
    }
    catch(...) {
      cerr << "The file could not be read. Is the file type really AFM? Quitting..." << endl;
    }
  } else if(fileType == ARFF) {
    //cout << "File type interpreted as Attribute-Relation File Format (ARFF)" << endl;
    try {
      Treedata::readARFF(featurestream,rawMatrix,featureHeaders,isFeatureNumerical);
    }
    catch(...) {
      cerr << "The file could not be read. Is the file type really ARFF? Quitting..." << endl;
    }
  } else {
    //cout << "File type is unknown -- defaulting to Annotated Feature Matrix (AFM)" << endl;
    try {
      Treedata::readAFM(featurestream,rawMatrix,featureHeaders,isFeatureNumerical);
    }
    catch(...) {
      cerr << "The file could not be read, quitting..." << endl;
    }
  }      

  if ( !datadefs::is_unique(featureHeaders) ) {
    cerr << "Feature headers are not unique!" << endl;
    assert(false);
  }

  nSamples_ = rawMatrix[0].size();
  nFeatures_ = featureHeaders.size();
  features_.resize(2*nFeatures_);

  for(size_t i = 0; i < nFeatures_; ++i) {
    vector<num_t> featureData(nSamples_);
    features_[i].name = featureHeaders[i];
    features_[i].isNumerical = isFeatureNumerical[i];
    if(features_[i].isNumerical) {
      datadefs::strv2numv(rawMatrix[i],featureData);
      features_[i].nCategories = 0;
    } else {
      map<string,num_t> mapping;
      map<num_t,string> backMapping;
      datadefs::strv2catv(rawMatrix[i], featureData, mapping, backMapping);
      features_[i].mapping = mapping;
      features_[i].backMapping = backMapping;
      map<num_t,size_t> freq;
      size_t nReal;
      datadefs::count_freq(featureData, freq, nReal);
      features_[i].nCategories = freq.size();
    }
    features_[i].data = featureData;
    //features_[i].rawData = rawMatrix[i];
    //features_[i].contrast = featureData;
  } 
  
  for(size_t i = nFeatures_; i < 2*nFeatures_; ++i) {
    features_[i] = features_[ i-nFeatures_ ];
    features_[i].name = "CONTRAST";
    //features_[i].isNumerical = features_[ i-nFeatures_ ].isNumerical;
    //features_[i].data = features_[ i-nFeatures_ ].data;
  }

  Treedata::permuteContrasts();

  //targetidx_ = 0;
  //Treedata::selectTarget(targetidx_);

}



Treedata::Treedata(Treedata& treedata):
  features_(0),
  nSamples_(0),
  nFeatures_(0) {

  time_t now;
  time(&now);
  randomInteger_.seed((unsigned int)now);

  features_ = treedata.features_;
  nFeatures_ = treedata.nFeatures_;
  nSamples_ = treedata.nSamples_;

}


Treedata::Treedata(Treedata& treedata, const vector<size_t>& featureIcs):
  features_(0),
  nSamples_(0),
  nFeatures_(0) {

  time_t now;
  time(&now);
  randomInteger_.seed((unsigned int)now);

  size_t nFeaturesNew = featureIcs.size();

  //We'll leave room for the contrasts
  features_.resize(2*nFeaturesNew);
  for(size_t i = 0; i < nFeaturesNew; ++i) {
    features_[i] = treedata.features_[ featureIcs[i] ];
    features_[i + nFeaturesNew] = treedata.features_[ featureIcs[i] + treedata.nFeatures() ];
  }
  
  nSamples_ = treedata.nSamples();
  nFeatures_ = nFeaturesNew;

  //cout << nSamples_ << " " << nFeatures_ << endl;
  
}

Treedata::~Treedata() {
}

void Treedata::keepFeatures(const vector<size_t>& featureIcs) {
  
  size_t nFeaturesNew = featureIcs.size();
  
  vector<Treedata::Feature> featureCopy = features_;
  features_.resize(2*nFeaturesNew);
  for(size_t i = 0; i < nFeaturesNew; ++i) {
    features_[i] = featureCopy[featureIcs[i]];
    features_[ nFeaturesNew + i ] = featureCopy[ nFeatures_ + featureIcs[i] ];
  }
  nFeatures_ = nFeaturesNew;
}

void Treedata::readFileType(string& fileName, FileType& fileType) {

  stringstream ss(fileName);
  string suffix = "";
  while(getline(ss,suffix,'.')) {}
  //datadefs::toupper(suffix);

  if(suffix == "AFM" || suffix == "afm") {
    fileType = AFM;
  } else if(suffix == "ARFF" || suffix == "arff") {
    fileType = ARFF;
  } else {
    fileType = UNKNOWN;
  }

}

void Treedata::readAFM(ifstream& featurestream, 
		       vector<vector<string> >& rawMatrix, 
		       vector<string>& featureHeaders, 
		       //vector<string>& sampleHeaders,
		       vector<bool>& isFeatureNumerical) {

  string field;
  string row;

  rawMatrix.clear();
  featureHeaders.clear();
  isFeatureNumerical.clear();

  //Remove upper left element from the matrix as useless
  getline(featurestream,field,'\t');

  //Next read the first row, which should contain the column headers
  getline(featurestream,row);
  stringstream ss( datadefs::chomp(row) );
  bool isFeaturesAsRows = true;
  vector<string> columnHeaders;
  while ( getline(ss,field,'\t') ) {

    // If at least one of the column headers is a valid feature header, we assume features are stored as columns
    if ( isFeaturesAsRows && isValidFeatureHeader(field) ) {
      isFeaturesAsRows = false;
    }
    columnHeaders.push_back(field);
  }

  // We should have reached the end of file. NOTE: failbit is set since the last element read did not end at '\t'
  assert( ss.eof() );
  //assert( !ss.fail() );

  size_t nColumns = columnHeaders.size();

  vector<string> rowHeaders;
  vector<string> sampleHeaders; // THIS WILL BE DEFINED AS ONE OF THE INPUT ARGUMENTS

  //Go through the rest of the rows
  while ( getline(featurestream,row) ) {

    row = datadefs::chomp(row);

    //Read row from the stream
    ss.clear();
    ss.str("");

    //Read the string back to a stream
    ss << row;

    //Read the next row header from the stream
    getline(ss,field,'\t');
    rowHeaders.push_back(field);

    vector<string> rawVector(nColumns);
    for(size_t i = 0; i < nColumns; ++i) {
      getline(ss,rawVector[i],'\t');
    }
    assert(!ss.fail());
    assert(ss.eof());
    rawMatrix.push_back(rawVector);
  }

  //If the data is row-formatted...
  if(isFeaturesAsRows) {
    //cout << "AFM orientation: features as rows" << endl;

    //... and feature headers are row headers
    featureHeaders = rowHeaders;
    sampleHeaders = columnHeaders;

  } else {

    //cout << "AFM orientation: features as columns" << endl;
      
    Treedata::transpose<string>(rawMatrix);
      
    //... and feature headers are row headers
    featureHeaders = columnHeaders;
    sampleHeaders = rowHeaders;
      
  }

  size_t nFeatures = featureHeaders.size();
  isFeatureNumerical.resize(nFeatures);
  for(size_t i = 0; i < nFeatures; ++i) {
    if(Treedata::isValidNumericalHeader(featureHeaders[i])) {
      isFeatureNumerical[i] = true;
    } else {
      isFeatureNumerical[i] = false;
    }
  }
}

void Treedata::readARFF(ifstream& featurestream, vector<vector<string> >& rawMatrix, vector<string>& featureHeaders, vector<bool>& isFeatureNumerical) {

  string row;

  bool hasRelation = false;
  bool hasData = false;

  size_t nFeatures = 0;
  //TODO: add Treedata::clearData(...)
  rawMatrix.clear();
  featureHeaders.clear();
  isFeatureNumerical.clear();
  
  //Read one line from the ARFF file
  while ( getline(featurestream,row) ) {

    //This is the final branch: once relation and attributes are read, and we find data header, we'll start reading the data in 
    if ( hasData && hasRelation ) {
      //There must be at least two attributes, otherwise the ARFF file makes no sense
      if ( nFeatures < 2 ) {
        cerr << "too few attributes ( < 2 ) found from the ARFF file" << endl;
        assert(false);
      }

      rawMatrix.resize(nFeatures);

      //Read data row-by-row
      while(true) {
        //++nsamples_;
        string field;
        stringstream ss(row);
        
        for(size_t attributeIdx = 0; attributeIdx < nFeatures; ++attributeIdx) {
          getline(ss,field,',');
          //cout << " " << field;
          rawMatrix[attributeIdx].push_back(field);
        }
        //cout << endl;

        if(!getline(featurestream,row)) {
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
    if(row[0] == '%' || row == "") {
      continue;  
    }

    //Read relation
    if(!hasRelation && row.compare(0,9,"@relation") == 0) {
      hasRelation = true;
      //cout << "found relation header: " << row << endl;
    } else if(row.compare(0,10,"@attribute") == 0) {    //Read attribute 
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
      // else
      //  {
      //    cout << "categorical)" << endl;
      //  }
    } else if(!hasData && row.compare(0,5,"@data") == 0) {    //Read data header 

      hasData = true;
      //cout << "found data header:" << row << endl;
    } else {      //If none of the earlier branches matched, we have a problem
      cerr << "incorrectly formatted ARFF file" << endl;
      assert(false);
    }
  }
}

void Treedata::parseARFFattribute(const string& str, string& attributeName, bool& isFeatureNumerical) {

  stringstream ss(str);
  string attributeHeader = "";
  attributeName = "";
  string attributeType = "";

  getline(ss,attributeHeader,' ');
  getline(ss,attributeName,' ');
  getline(ss,attributeType);

  //string prefix;
  if(datadefs::toUpperCase(attributeType) == "NUMERIC" || datadefs::toUpperCase(attributeType) == "REAL" ) {
    isFeatureNumerical = true;
  } else {
    isFeatureNumerical = false;
  }
  //prefix.append(attributeName);
  //attributeName = prefix;
}

bool Treedata::isValidNumericalHeader(const string& str) {
  return( str.compare(0,2,"N:") == 0 );
}

bool Treedata::isValidCategoricalHeader(const string& str) {
  return( str.compare(0,2,"C:") == 0 || str.compare(0,2,"B:") == 0 );
}

bool Treedata::isValidFeatureHeader(const string& str) {
  return( isValidNumericalHeader(str) || isValidCategoricalHeader(str) );
}

size_t Treedata::nFeatures() {
  return(nFeatures_);
}

size_t Treedata::nSamples() {
  return(nSamples_);
}

//WILL BECOME DEPRECATED
num_t Treedata::pearsonCorrelation(size_t featureidx1, size_t featureidx2) {
  num_t r;
  datadefs::pearson_correlation(features_[featureidx1].data,features_[featureidx2].data,r);
  return(r);
}

void Treedata::getMatchingTargetIdx(const string& targetStr, size_t& targetIdx) {
  
  bool isFoundAlready = false;

  for ( size_t featureIdx = 0; featureIdx < nFeatures_; ++featureIdx ) {
    
    if ( features_[featureIdx].name == targetStr ) {
    
      if ( isFoundAlready ) {
	cerr << "Multiple instances of the same target found in the data!" << endl;
	assert(false);
      }
      
      isFoundAlready = true;
      targetIdx = featureIdx;
    
    }
  }

  if ( !isFoundAlready ) {
    cerr << "Feature '" << targetStr << "' not found" << endl;
    assert(false);
  }

}

string Treedata::getFeatureName(const size_t featureIdx) {
  return(features_[featureIdx].name);
}


void Treedata::print() {
  cout << "Printing feature matrix (missing values encoded to " << datadefs::NUM_NAN << "):" << endl;
  for(size_t j = 0; j < nSamples_; ++j) {
    cout << '\t' << "foo";
  }
  cout << endl;
  for(size_t i = 0; i < nFeatures_; ++i) {
    cout << i << ':' << features_[i].name << ':';
    for(size_t j = 0; j < nSamples_; ++j) {
      cout << '\t' << features_[i].data[j];
    }
    cout << endl;
  }
}


void Treedata::print(const size_t featureIdx) {
  cout << "Print " << features_[featureIdx].name << ":";
  for(size_t i = 0; i < nSamples_; ++i) {
    cout << " " << features_[featureIdx].data[i];
  }
  cout << endl;
}


void Treedata::permuteContrasts() {

  for(size_t i = nFeatures_; i < 2*nFeatures_; ++i) {
    Treedata::permute(features_[i].data);
  }

}

bool Treedata::isFeatureNumerical(size_t featureIdx) {

  return(features_[featureIdx].isNumerical);

}


size_t Treedata::nRealSamples(const size_t featureIdx) { 
  
  size_t nRealSamples;
  datadefs::countRealValues(features_[featureIdx].data,nRealSamples);
  return(nRealSamples);

}

size_t Treedata::nRealSamples(const size_t featureIdx1, const size_t featureIdx2) {

  size_t nRealSamples = 0;
  for(size_t i = 0; i < nSamples_; ++i ) {
    if( !datadefs::isNAN(features_[featureIdx1].data[i]) && !datadefs::isNAN(features_[featureIdx2].data[i])) {
      ++nRealSamples;
    }
  }
  return(nRealSamples);
}

size_t Treedata::nCategories(const size_t featureIdx) {

  return(features_[featureIdx].nCategories);

}


//WILL BE DEPRECATED
template <typename T> void Treedata::transpose(vector<vector<T> >& mat) {

  vector<vector<T> > foo = mat;

  size_t ncols = mat.size();
  size_t nrows = mat[0].size();

  mat.resize(nrows);
  for(size_t i = 0; i < nrows; ++i) {
    mat[i].resize(ncols);
  }

  for(size_t i = 0; i < nrows; ++i) {
    for(size_t j = 0; j < ncols; ++j) {
      mat[i][j] = foo[j][i];
    }
  }
}


void Treedata::permute(vector<size_t>& ics) {
  for (size_t i = 0; i < ics.size(); ++i) {
    size_t j = randomInteger_() % (i + 1);
    ics[i] = ics[j];
    ics[j] = i;
  }
}

void Treedata::permute(vector<num_t>& data) {
  size_t n = data.size();
  vector<size_t> ics(n);

  Treedata::permute(ics);

  for(size_t i = 0; i < n; ++i) {
    num_t temp = data[i];
    data[i] = data[ics[i]];
    data[ics[i]] = temp;
  }
}


void Treedata::bootstrapFromRealSamples(const bool withReplacement, 
                                        const num_t sampleSize, 
                                        const size_t featureIdx, 
                                        vector<size_t>& ics, 
                                        vector<size_t>& oobIcs) {
    
  //Check that the sampling parameters are appropriate
  assert(sampleSize > 0.0);
  if(!withReplacement && sampleSize > 1.0) {
    cerr << "when sampling without replacement, sample size must be less or equal to 100% (sampleSize <= 1.0)" << endl;
    assert(false);
  }

  //First we collect all indices that correspond to real samples
  vector<size_t> allIcs;
  for(size_t i = 0; i < nSamples_; ++i) {
    if(!datadefs::isNAN(features_[featureIdx].data[i])) {
      allIcs.push_back(i);
    }
  }
  
  //Extract the number of real samples, and see how many samples do we have to collect
  size_t nRealSamples = allIcs.size();
  size_t nSamples = static_cast<size_t>( floor( sampleSize * nRealSamples ) );
  ics.resize(nSamples);
  
  //If sampled with replacement...
  if(withReplacement) {
    //Draw nSamples random integers from range of allIcs
    for(size_t sampleIdx = 0; sampleIdx < nSamples; ++sampleIdx) {
      ics[sampleIdx] = allIcs[randomInteger_() % nRealSamples];
    }
  } else {  //If sampled without replacement...
    vector<size_t> foo(nRealSamples);
    Treedata::permute(foo);
    for(size_t i = 0; i < nSamples; ++i) {
      ics[i] = allIcs[foo[i]];
    }
  }

  sort(ics.begin(),ics.end());

  if(nSamples < nRealSamples) {
    oobIcs.resize(nRealSamples);
  } else {
    oobIcs.resize(nSamples);
  }

  //Then, as we now have the sample stored in ics, we'll check which of the samples, from allIcs, are not contained in ics and store them in oobIcs instead
  vector<size_t>::iterator it = set_difference(allIcs.begin(),allIcs.end(),ics.begin(),ics.end(),oobIcs.begin());
  size_t nOob = distance(oobIcs.begin(),it);
  oobIcs.resize(nOob);
  //cout << "nOob=" << nOob << endl;
}

void Treedata::getFeatureData(size_t featureIdx, vector<num_t>& data) {
  data.resize(nSamples_);
  for(size_t i = 0; i < nSamples_; ++i) {
    data[i] = features_[featureIdx].data[i];
  }
}

void Treedata::getFeatureData(size_t featureIdx, const size_t sampleIdx, num_t& data) {
  data = features_[featureIdx].data[sampleIdx];
}

void Treedata::getFeatureData(size_t featureIdx, const vector<size_t>& sampleIcs, vector<num_t>& data) {
  data.resize(sampleIcs.size());
  for(size_t i = 0; i < sampleIcs.size(); ++i) {
    data[i] = features_[featureIdx].data[sampleIcs[i]];
  }

}

string Treedata::getRawFeatureData(const size_t featureIdx, const size_t sampleIdx) {

  num_t value = features_[featureIdx].data[sampleIdx];

  if(datadefs::isNAN(value)) {
    return(datadefs::STR_NAN);
  } else {
    if(features_[featureIdx].isNumerical) {
      stringstream ss;
      ss << value;
      return(ss.str());
    } else {
      return(features_[featureIdx].backMapping[ value ]);
    }
  }
    
}

// DEPRECATED ??
void Treedata::impurity(vector<num_t>& data, bool isFeatureNumerical, num_t& impurity, size_t& nreal) {
  
  size_t n = data.size();
  
  
  nreal = 0;
  if(isFeatureNumerical) {
    num_t mu = 0.0;
    num_t se = 0.0;
    for(size_t i = 0; i < n; ++i) {
      datadefs::forward_sqerr(data[i],nreal,mu,se);  
    }
    impurity = se / nreal;
  } else {
    map<num_t,size_t> freq;
    size_t sf = 0;
    for(size_t i = 0; i < n; ++i) {
      datadefs::forward_sqfreq(data[i],nreal,freq,sf);
    }
    impurity = 1.0 - 1.0 * sf / (nreal * nreal);
  }
}
