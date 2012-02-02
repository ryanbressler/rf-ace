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

Treedata::Treedata(string fileName, char dataDelimiter, char headerDelimiter):
  dataDelimiter_(dataDelimiter),
  headerDelimiter_(headerDelimiter),
  features_(0),
  sampleHeaders_(0) {

  //Initialize stream to read from file
  ifstream featurestream;
  featurestream.open(fileName.c_str());
  if ( !featurestream.good() ) {
    cerr << "Failed to open file '" << fileName << "' for reading. Make sure the file exists. Quitting..." << endl;
    exit(1);
  }

  FileType fileType = UNKNOWN;
  Treedata::readFileType(fileName,fileType);

  vector<vector<string> > rawMatrix;
  vector<string> featureHeaders;
  vector<bool> isFeatureNumerical;
  if(fileType == AFM) {
    
    Treedata::readAFM(featurestream,rawMatrix,featureHeaders,sampleHeaders_,isFeatureNumerical);
    
  } else if(fileType == ARFF) {
    
    Treedata::readARFF(featurestream,rawMatrix,featureHeaders,isFeatureNumerical);
    sampleHeaders_.clear();
    sampleHeaders_.resize(rawMatrix[0].size(),"NO_SAMPLE_ID");
    
  } else {
    
    Treedata::readAFM(featurestream,rawMatrix,featureHeaders,sampleHeaders_,isFeatureNumerical);
    
  }      

  size_t nFeatures = featureHeaders.size();
  features_.resize(2*nFeatures);

  for(size_t i = 0; i < nFeatures; ++i) {
    
    if( name2idx_.find(featureHeaders[i]) != name2idx_.end() ) {
      cerr << "Duplicate feature header '" << featureHeaders[i] << "' found!" << endl;
      exit(1);
    }

    //name2idxHashTest_[featureHeaders[i]] = i;
    //cout << " " << name2idxHashTest_[featureHeaders[i]];
    name2idx_[featureHeaders[i]] = i;
    features_[i].name = featureHeaders[i];
    features_[i].isNumerical = isFeatureNumerical[i];

    if(features_[i].isNumerical) {

      datadefs::strv2numv(rawMatrix[i],features_[i].data);

      this->updateSortOrder(i);

    } else {

      features_[i].sortOrder.clear();

      datadefs::strv2catv(rawMatrix[i], 
			  features_[i].data, 
			  features_[i].mapping, 
			  features_[i].backMapping);

    }

  } 
  
  //cout << endl;

  for(size_t i = nFeatures; i < 2*nFeatures; ++i) {
    features_[i] = features_[ i - nFeatures ];
    string contrastName = features_[ i - nFeatures ].name;
    features_[i].name = contrastName.append("_CONTRAST");
    name2idx_[contrastName] = i;
  }

  // Initialize the Mersenne Twister RNG with the CPU cycle count as the seed
  time_t now;
  time(&now);
  unsigned int seed = clock() + now;
  randomInteger_.seed(seed);

  //cout << "permuting contrasts..." << endl;

  Treedata::permuteContrasts();

  //cout << "done permuting" << endl;

}

Treedata::~Treedata() {
  /* Empty destructor */
}

void Treedata::whiteList(const set<string>& featureNames ) {
  
  vector<bool> keepFeatureIcs(this->nFeatures(),false);

  for ( set<string>::const_iterator it( featureNames.begin() ); it != featureNames.end(); ++it ) {
    keepFeatureIcs[this->getFeatureIdx(*it)] = true;
  }

  this->whiteList(keepFeatureIcs);

}

void Treedata::blackList(const set<string>& featureNames ) {
  
  vector<bool> keepFeatureIcs(this->nFeatures(),true);
  
  for ( set<string>::const_iterator it( featureNames.begin() ); it != featureNames.end(); ++it ) {
    keepFeatureIcs[this->getFeatureIdx(*it)] = false;
  }
  
  this->whiteList(keepFeatureIcs);

}

void Treedata::whiteList(const vector<bool>& keepFeatureIcs) {
  
  assert( keepFeatureIcs.size() == this->nFeatures() );
  
  size_t nFeaturesNew = 0;
  for ( size_t i = 0; i < keepFeatureIcs.size(); ++i ) {
    if ( keepFeatureIcs[i] ) {
      ++nFeaturesNew;
    }
  }

  vector<Feature> featureCopy = features_;
  features_.resize(2*nFeaturesNew);

  map<string,size_t> name2idxCopy = name2idx_;

  name2idx_.clear();

  size_t nFeatures = keepFeatureIcs.size();
  size_t iter = 0;
  for ( size_t i = 0; i < nFeatures; ++i ) {

    if ( !keepFeatureIcs[i] ) {
      continue;
    }
    
    string featureName = featureCopy[i].name;

    if ( name2idxCopy.find(featureName) == name2idxCopy.end() ) {
      cerr << "Treedata::keepFeatures() -- feature '" << featureName << "' does not exist" << endl;
      exit(1);
    }

    features_[ iter ] = featureCopy[ name2idxCopy[featureName] ];
    name2idx_[ featureName ] = iter;

    featureName.append("_CONTRAST");

    if ( name2idxCopy.find(featureName) == name2idxCopy.end() ) {
      cerr << "Treedata::keepFeatures() -- feature '" << featureName << "' does not exist" << endl;
      exit(1);
    }

    features_[ nFeaturesNew + iter ] = featureCopy[ name2idxCopy[featureName] ];
    name2idx_[ featureName ] = nFeaturesNew + iter; 

    ++iter;

  }

  assert( iter == nFeaturesNew );
  assert( 2*iter == features_.size() );
  
}

void Treedata::updateSortOrder(const size_t featureIdx) {

  features_[featureIdx].sortOrder.resize( sampleHeaders_.size() );

  vector<size_t> refIcs( sampleHeaders_.size() );

  vector<num_t> foo = features_[featureIdx].data;
  bool isIncreasingOrder = true;
  datadefs::sortDataAndMakeRef(isIncreasingOrder,foo,refIcs);

  for( size_t i = 0; i < refIcs.size(); ++i ) {
    features_[featureIdx].sortOrder[refIcs[i]] = i;
  }

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
		       vector<string>& sampleHeaders,
		       vector<bool>& isFeatureNumerical) {

  string field;
  string row;

  rawMatrix.clear();
  featureHeaders.clear();
  sampleHeaders.clear();
  isFeatureNumerical.clear();

  //Remove upper left element from the matrix as useless
  getline(featurestream,field,dataDelimiter_);

  //Next read the first row, which should contain the column headers
  getline(featurestream,row);
  stringstream ss( datadefs::chomp(row) );
  bool isFeaturesAsRows = true;
  vector<string> columnHeaders;
  while ( getline(ss,field,dataDelimiter_) ) {

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
  //vector<string> sampleHeaders; // THIS WILL BE DEFINED AS ONE OF THE INPUT ARGUMENTS

  //Go through the rest of the rows
  while ( getline(featurestream,row) ) {

    row = datadefs::chomp(row);

    //Read row from the stream
    ss.clear();
    ss.str("");

    //Read the string back to a stream
    ss << row;

    //Read the next row header from the stream
    getline(ss,field,dataDelimiter_);
    rowHeaders.push_back(field);

    vector<string> rawVector(nColumns);
    for(size_t i = 0; i < nColumns; ++i) {
      getline(ss,rawVector[i],dataDelimiter_);
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

      break;
    }

    //Comment lines and empty lines are omitted
    if(row[0] == '%' || row == "") {
      continue;  
    }

    string rowU = datadefs::toUpperCase(row);

    //Read relation
    if(!hasRelation && rowU.compare(0,9,"@RELATION") == 0) {
      hasRelation = true;
      //cout << "found relation header: " << row << endl;
    } else if ( rowU.compare(0,10,"@ATTRIBUTE") == 0) {    //Read attribute 
      string attributeName = "";
      bool isNumerical;
      ++nFeatures;
      //cout << "found attribute header: " << row << endl;
      Treedata::parseARFFattribute(row,attributeName,isNumerical);
      featureHeaders.push_back(attributeName);
      isFeatureNumerical.push_back(isNumerical);

    } else if(!hasData && rowU.compare(0,5,"@DATA") == 0) {    //Read data header 

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
  
  stringstream ss(str);
  string typeStr;
  getline(ss,typeStr,headerDelimiter_);

  return(  typeStr == "N" );
}

bool Treedata::isValidCategoricalHeader(const string& str) {

  stringstream ss(str);
  string typeStr;
  getline(ss,typeStr,headerDelimiter_);

  return( typeStr == "C" || typeStr == "B" );
}

bool Treedata::isValidFeatureHeader(const string& str) {
  return( isValidNumericalHeader(str) || isValidCategoricalHeader(str) );
}

size_t Treedata::nFeatures() {
  return( features_.size() / 2 );
}

size_t Treedata::nSamples() {
  return( sampleHeaders_.size() );
}

// WILL BECOME DEPRECATED
num_t Treedata::pearsonCorrelation(size_t featureidx1, size_t featureidx2) {
  num_t r;
  datadefs::pearson_correlation(features_[featureidx1].data,features_[featureidx2].data,r);
  return(r);
}

size_t Treedata::getFeatureIdx(const string& featureName) {
  if ( name2idx_.find(featureName) == name2idx_.end() ) {
    cerr << "Treedata::getFeatureIdx() -- feature '" << featureName << "' does not exist" << endl;
    exit(1);
  }
  return( name2idx_[featureName] );
}

string Treedata::getFeatureName(const size_t featureIdx) {
  return(features_[featureIdx].name);
}

string Treedata::getSampleName(const size_t sampleIdx) {
  return(sampleHeaders_[sampleIdx]);
}


void Treedata::print() {
  cout << "Printing feature matrix (missing values encoded to " << datadefs::NUM_NAN << "):" << endl;
  for(size_t j = 0; j < Treedata::nSamples(); ++j) {
    cout << '\t' << "foo";
  }
  cout << endl;
  for(size_t i = 0; i < Treedata::nFeatures(); ++i) {
    cout << i << ':' << features_[i].name << ':';
    for(size_t j = 0; j < Treedata::nSamples(); ++j) {
      cout << '\t' << features_[i].data[j];
    }
    cout << endl;
  }
}


void Treedata::print(const size_t featureIdx) {
  cout << "Print " << features_[featureIdx].name << ":";
  for(size_t i = 0; i < Treedata::nSamples(); ++i) {
    cout << " " << features_[featureIdx].data[i];
  }
  cout << endl;
}


void Treedata::permuteContrasts() {

  size_t nFeatures = this->nFeatures();
  size_t nSamples = this->nSamples();

  for ( size_t i = nFeatures; i < 2*nFeatures; ++i ) {
    

    vector<size_t> sampleIcs(nSamples);
    datadefs::range(sampleIcs);

    vector<num_t> filteredData = this->getFilteredFeatureData(i,sampleIcs);
    
    Treedata::permute(filteredData);

    //datadefs::print(features_[i].data);

    for ( size_t j = 0; j < sampleIcs.size(); ++j ) {
      features_[i].data[ sampleIcs[j] ] = filteredData[j];
    }

    //datadefs::print(features_[i].data);
  }

}

bool Treedata::isFeatureNumerical(const size_t featureIdx) {
  return(features_[featureIdx].isNumerical);
}

size_t Treedata::nRealSamples(const size_t featureIdx) { 
  
  size_t nRealSamples;
  datadefs::countRealValues( features_[featureIdx].data, nRealSamples );
  return( nRealSamples );

}

size_t Treedata::nRealSamples(const size_t featureIdx1, const size_t featureIdx2) {

  size_t nRealSamples = 0;
  for( size_t i = 0; i < Treedata::nSamples(); ++i ) {
    if( !datadefs::isNAN( features_[featureIdx1].data[i] ) && !datadefs::isNAN( features_[featureIdx2].data[i] ) ) {
      ++nRealSamples;
    }
  }
  return( nRealSamples );
}

size_t Treedata::nCategories(const size_t featureIdx) {
  return( features_[featureIdx].mapping.size() );
}

size_t Treedata::nMaxCategories() {

  size_t ret = 0;
  for( size_t i = 0; i < Treedata::nFeatures(); ++i ) {
    if( ret < features_[i].mapping.size() ) {
      ret = features_[i].mapping.size();
    }
  }
  
  return( ret ); 

}

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
    cerr << "Treedata: when sampling without replacement, sample size must be less or equal to 100% (sampleSize <= 1.0)" << endl;
    exit(1);
  }

  //First we collect all indices that correspond to real samples
  vector<size_t> allIcs;
  for(size_t i = 0; i < Treedata::nSamples(); ++i) {
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


vector<num_t> Treedata::getFeatureData(size_t featureIdx) {
  
  vector<num_t> data( features_[featureIdx].data.size() );

  for(size_t i = 0; i < Treedata::nSamples(); ++i) {
    data[i] = features_[featureIdx].data[i];
  }

  return( data );
}


num_t Treedata::getFeatureData(size_t featureIdx, const size_t sampleIdx) {

  num_t data = features_[featureIdx].data[sampleIdx];

  return( data ); 
}

vector<num_t> Treedata::getFeatureData(size_t featureIdx, const vector<size_t>& sampleIcs) {
  
  vector<num_t> data(sampleIcs.size());
  
  for(size_t i = 0; i < sampleIcs.size(); ++i) {
    data[i] = features_[featureIdx].data[sampleIcs[i]];
  }

  return( data );

}

vector<num_t> Treedata::getFilteredFeatureData(const size_t featureIdx,
					       vector<size_t>& sampleIcs) {

  size_t n = sampleIcs.size();

  vector<num_t> featureData(n);

  size_t nReal = 0;

  for ( size_t i = 0; i < n; ++i ) {
    size_t idx = sampleIcs[i];
    num_t value = features_[featureIdx].data[idx];
    if ( !datadefs::isNAN(value) ) {
      featureData[nReal] = value;
      sampleIcs[nReal] = idx;
      ++nReal;
    }
  }
  sampleIcs.resize(nReal);
  featureData.resize(nReal);

  return(featureData);

}

void Treedata::getFilteredFeatureDataPair(const size_t featureIdx1, 
					  const size_t featureIdx2, 
					  vector<size_t>& sampleIcs, 
					  vector<num_t>& featureData1, 
					  vector<num_t>& featureData2) {

  size_t n = sampleIcs.size();
  featureData1.resize(n);
  featureData2.resize(n);
  size_t nReal = 0;
  for(size_t i = 0; i < n; ++i) {

    num_t v1 = features_[featureIdx1].data[sampleIcs[i]];
    num_t v2 = features_[featureIdx2].data[sampleIcs[i]];
    
    if(!datadefs::isNAN(v1) && !datadefs::isNAN(v2)) {
      sampleIcs[nReal] = sampleIcs[i];
      featureData1[nReal] = v1;
      featureData2[nReal] = v2;
      ++nReal;
    }
  }
  featureData1.resize(nReal);
  featureData2.resize(nReal);
  sampleIcs.resize(nReal);

}

void Treedata::getFilteredAndSortedFeatureDataPair(const size_t targetIdx, 
						   const size_t featureIdx, 
						   vector<size_t>& sampleIcs, 
						   vector<num_t>& targetData, 
						   vector<num_t>& featureData) {

  if ( !features_[featureIdx].isNumerical ) {
    cerr << "Treedata::getFilteredAndSortedDataPair() -- cannot perform for CATEGORICAL features" << endl;
    exit(1);
  }

  targetData.clear();
  //targetData.resize( sampleHeaders_.size(), datadefs::NUM_NAN );
  featureData.clear();
  //featureData.resize( sampleHeaders_.size(), datadefs::NUM_NAN );

  //vector<size_t> sampleIcsCopy( sampleHeaders_.size() );
  //size_t maxPos = 0;

  // A map: sortOrderKey -> (sampleIdx,multiplicity)
  map<size_t,pair<size_t,size_t> > mapOrder;
  
  // Count the number of real samples
  size_t nReal = 0;

  // Go through all sample indices
  for ( vector<size_t>::const_iterator it(sampleIcs.begin()); it != sampleIcs.end(); ++it ) {
    
    // Extract the target and feature values for the index
    num_t tVal = features_[targetIdx].data[*it];
    num_t fVal = features_[featureIdx].data[*it];

    // If the data are non-NA...
    if ( !datadefs::isNAN(fVal) && !datadefs::isNAN(tVal) ) {
    
      // Accumulate real data counter
      ++nReal;

      // Extract the ordered position of the sample
      size_t pos = features_[featureIdx].sortOrder[*it];
    
      // If the position is unused in the map...
      if ( mapOrder.find(pos) == mapOrder.end() ) {
	 
	// Add the ordered position, the original sample index, 
	// and initialize the sample counter to 1
	pair<size_t,size_t> foo(*it,1);
	mapOrder.insert(pair<size_t,pair<size_t,size_t> >(pos,foo));
      } else {

	// Otherwise accumulate multiplicity by one
	++mapOrder[pos].second;
      }
    }
  }

  targetData.resize(nReal);
  featureData.resize(nReal);
  sampleIcs.resize(nReal);

  size_t i = 0;
  
  for ( map<size_t,pair<size_t,size_t> >::const_iterator it(mapOrder.begin()); it != mapOrder.end(); ++it ) {
    
    for ( size_t j = 0; j < it->second.second; ++j ) {
      sampleIcs[i] = it->second.first;
      targetData[i] = features_[targetIdx].data[it->second.first];
      featureData[i] = features_[featureIdx].data[it->second.first];
      ++i;
    }
  }

  assert(i == nReal);

}

void Treedata::getFilteredAndSortedFeatureDataPair2(const size_t targetIdx,
						    const size_t featureIdx,
						    vector<size_t>& sampleIcs,
						    vector<num_t>& targetData,
						    vector<num_t>& featureData) {

  if ( !features_[featureIdx].isNumerical ) {
    cerr << "Treedata::getFilteredAndSortedDataPair() -- cannot perform for CATEGORICAL features" << endl;
    exit(1);
  }

  size_t n = sampleHeaders_.size();
  size_t s = sampleIcs.size();

  vector<num_t> targetDataCopy(n);
  vector<num_t> featureDataCopy(n);
  vector<size_t> sampleIcsCopy(n);
  vector<size_t> multiplicity(n, 0);

  //vector<size_t> sampleIcsCopy(  );
  size_t minPos = n;
  size_t maxPos = 0;

  // Count the number of real samples
  size_t nReal = 0;

  // Go through all sample indices
  for ( size_t i = 0; i < s; ++i ) {

    size_t ii = sampleIcs[i];

    // Extract the target and feature values for the index
    num_t tVal = features_[targetIdx].data[ii];
    num_t fVal = features_[featureIdx].data[ii];

    // If the data are non-NA...
    if ( !datadefs::isNAN(tVal) && !datadefs::isNAN(fVal) ) {

      // Accumulate real data counter
      ++nReal;

      // Extract the ordered position of the sample
      size_t pos = features_[featureIdx].sortOrder[ii];
      ++multiplicity[pos];

      if ( multiplicity[pos] == 1 ) {
	featureDataCopy[pos] = fVal;
	targetDataCopy[pos] = tVal;
	sampleIcsCopy[pos] = ii;

	if ( pos > maxPos ) {
	  maxPos = pos;
	}

	if ( pos < minPos ) {
	  minPos = pos;
	}

      }
      
    }
  }

  featureData.resize(nReal);
  targetData.resize(nReal);
  sampleIcs.resize(nReal);

  size_t iter = 0;
  for ( size_t i = minPos; i <= maxPos; ++i ) {
    for ( size_t j = 0; j < multiplicity[i]; ++j ) {
      featureData[iter] = featureDataCopy[i];
      targetData[iter] = targetDataCopy[i];
      sampleIcs[iter] = sampleIcsCopy[i];
      ++iter;
    }
  }

  assert(nReal == iter);
 
}

void Treedata::getFilteredAndSortedFeatureDataPair3(const size_t targetIdx,
						    const size_t featureIdx,
						    vector<size_t>& sampleIcs,
						    vector<num_t>& targetData,
						    vector<num_t>& featureData) {


  featureData = this->getFeatureData(featureIdx,sampleIcs);
  //targetData = this->getFeatureData(targetIdx,sampleIcs);

  bool isIncreasingOrder = true;
  vector<size_t> refIcs;

  datadefs::sortDataAndMakeRef(isIncreasingOrder,featureData,refIcs);

  vector<size_t> sampleIcsCopy = sampleIcs;

  for ( size_t i = 0; i < refIcs.size(); ++i ) {
    sampleIcs[i] = sampleIcsCopy[refIcs[i]];
  }
  sampleIcs.resize(refIcs.size());

  targetData = this->getFeatureData(targetIdx,sampleIcs);

}

string Treedata::getRawFeatureData(const size_t featureIdx, const size_t sampleIdx) {

  num_t data = features_[featureIdx].data[sampleIdx];

  return( this->getRawFeatureData(featureIdx,data) );
    
}

string Treedata::getRawFeatureData(const size_t featureIdx, const num_t data) {

  if ( datadefs::isNAN(data) ) {
    return( datadefs::STR_NAN );
  } else {
    if ( features_[featureIdx].isNumerical ) {
      stringstream ss;
      ss << data;
      return( ss.str() );
    } else {
      return( features_[featureIdx].backMapping[data] );
    }
  }
}

vector<string> Treedata::getRawFeatureData(const size_t featureIdx) {
  
  vector<string> rawData( sampleHeaders_.size() );

  for ( size_t i = 0; i < rawData.size(); ++i ) {
    rawData[i] = this->getRawFeatureData(featureIdx,i);
  }

  return( rawData );

}

void Treedata::replaceFeatureData(const size_t featureIdx, const vector<num_t>& featureData) {

  if(featureData.size() != features_[featureIdx].data.size() ) {
    cerr << "Treedata::replaceFeatureData(num_t) -- data dimension mismatch" << endl;
    exit(1);
  }

  // Since the data that was passed is numerical, we set isNumerical to true
  features_[featureIdx].isNumerical = true;

  // Data that is stored is directly the input data
  features_[featureIdx].data = featureData;

  // Update sort indices for fast lookup
  this->updateSortOrder(featureIdx);

  // Since the data is not categorical, there's no need to provide mappings
  features_[featureIdx].mapping.clear();
  features_[featureIdx].backMapping.clear();

}

void Treedata::replaceFeatureData(const size_t featureIdx, const vector<string>& rawFeatureData) {

  if(rawFeatureData.size() != features_[featureIdx].data.size() ) {
    cerr << "Treedata::replaceFeatureData(string) -- data dimension mismatch" << endl;
    exit(1);
  }

  // Since the data that was passed are string literals, we set isNumerical to false
  features_[featureIdx].isNumerical = false;

  // Categorical data does not need sorting, thus, it doesn't benefit from the sort indices either
  features_[featureIdx].sortOrder.clear();

  // The string literal data needs some processing 
  datadefs::strv2catv(rawFeatureData,
		      features_[featureIdx].data,
		      features_[featureIdx].mapping,
		      features_[featureIdx].backMapping);
}



