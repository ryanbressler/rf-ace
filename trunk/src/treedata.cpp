#include "treedata.hpp"
#include <cstdlib>
#include <fstream>
#include <cassert>
#include <iostream>
#include <sstream>
#include <utility>
#include <algorithm>
#include <ctime>

#include "math.hpp"
#include "utils.hpp"
//#include "timer.hpp"

//extern Timer* TIMER_G;

using namespace std;

Feature::Feature() {

  
}

Feature::Feature(const vector<num_t>& newData, const string& newName) {

  data = newData;
  name = newName;
  isNumerical = true;
  mapping.clear();
  backMapping.clear();

}

Feature::Feature(const vector<string>& newStringData, const string& newName) {

  utils::strv2catv(newStringData,data,mapping,backMapping);
  
  name = newName;
  isNumerical = false;

}


Treedata::Treedata(const vector<Feature>& features, options::General_options* parameters, const vector<string>& sampleHeaders):
  parameters_(parameters),
  features_(features),
  sampleHeaders_(sampleHeaders) {

  size_t nFeatures = features_.size();

  assert( nFeatures > 0 );

  // If we have contrasts, there would be 2*nFeatures, in which case
  // 4*nFeatures results in a reasonable max load factor of 0.5
  name2idx_.rehash(4*nFeatures);

  size_t nSamples = features_[0].data.size();

  for ( size_t featureIdx = 0; featureIdx < nFeatures; ++featureIdx ) {

    assert( features_[featureIdx].data.size() == nSamples );
    name2idx_[ features_[featureIdx].name ] = featureIdx;
  }

  assert( nSamples > 0 );

  if ( sampleHeaders_.size() == 0 ) {
    sampleHeaders_.resize(nSamples,"NO_SAMPLE_ID");
  } 

  assert( sampleHeaders_.size() == nSamples );

  if ( parameters_->useContrasts ) {
    this->createContrasts(); // Doubles matrix size
    this->permuteContrasts(parameters_->randIntGens[0]);
  }
  
}


/**
   Reads a data file into a Treedata object. The data file can be either AFM or ARFF
   NOTE: dataDelimiter and headerDelimiter are used only when the format is AFM, for 
   ARFF default delimiter (comma) is used 
*/
Treedata::Treedata(string fileName, options::General_options* parameters):
  parameters_(parameters) {

  //TIMER_G->tic("READ");

  //Initialize stream to read from file
  ifstream featurestream;
  featurestream.open(fileName.c_str());
  if ( !featurestream.good() ) {
    cerr << "Failed to open file '" << fileName << "' for reading. Make sure the file exists. Quitting..." << endl;
    exit(1);
  }

  // Interprets file type from the content of the file
  FileType fileType = UNKNOWN;
  Treedata::readFileType(fileName,fileType);

  // Reads raw data matrix from the input file
  // NOTE: should be optimized to scale for large data sets
  vector<vector<string> > rawMatrix;
  vector<string> featureHeaders;
  vector<bool> isFeatureNumerical;
  if(fileType == AFM) {
    
    // Reads from the AFM stream and inserts 
    // data to the following arguments
    Treedata::readAFM(featurestream,
		      rawMatrix,
		      featureHeaders,
		      sampleHeaders_,
		      isFeatureNumerical);
    
  } else if(fileType == ARFF) {
    
    // Reads from the ARFF stream and inserts 
    // data to the following arguments
    Treedata::readARFF(featurestream,
		       rawMatrix,
		       featureHeaders,
		       isFeatureNumerical);

    // ARFF doesn't contain sample headers 
    sampleHeaders_.clear();
    sampleHeaders_.resize(rawMatrix[0].size(),"NO_SAMPLE_ID");
    
  } else {
    
    // By default, reads from the AFM stream and inserts
    // data to the following arguments
    Treedata::readAFM(featurestream,
		      rawMatrix,
		      featureHeaders,
		      sampleHeaders_,
		      isFeatureNumerical);
    
  }      

  // Extract the number of features 
  size_t nFeatures = featureHeaders.size();

  features_.resize(nFeatures);

  // If we have contrasts, there would be 2*nFeatures, in which case
  // 4*nFeatures results in a reasonable max load factor of 0.5
  name2idx_.rehash(4*nFeatures);

  // Start reading data to the final container "features_"
  for(size_t i = 0; i < nFeatures; ++i) {
    
    // We require that no two features have identical header
    if( name2idx_.find(featureHeaders[i]) != name2idx_.end() ) {
      cerr << "Duplicate feature header '" << featureHeaders[i] << "' found!" << endl;
      exit(1);
    }

    // Map the i'th feature header to integer i
    // NOTE: could be replaced with a hash table
    name2idx_[featureHeaders[i]] = i;

    if ( isFeatureNumerical[i] ) {

      vector<num_t> data;

      // If type is numerical, read the raw data as numbers
      utils::strv2numv(rawMatrix[i],data);

      features_[i] = Feature(data,featureHeaders[i]);

    } else {

      features_[i] = Feature(rawMatrix[i],featureHeaders[i]);

    }

  } 
 
  if ( parameters_->useContrasts ) {
    this->createContrasts(); // Doubles matrix size
    this->permuteContrasts(parameters_->randIntGens[0]);
  }
  
  //TIMER_G->toc("READ");

}

Treedata::~Treedata() {
  /* Empty destructor */
}

void Treedata::createContrasts() {

  // Resize the feature data container to fit the
  // original AND contrast features ( so 2*nFeatures )
  size_t nFeatures = features_.size();
  features_.resize(2*nFeatures);

  // Generate contrast features
  for(size_t i = nFeatures; i < 2*nFeatures; ++i) {
    features_[i] = features_[ i - nFeatures ];
    features_[i].name = features_[ i - nFeatures ].name;
    features_[i].name.append("_CONTRAST");
    name2idx_[ features_[i].name ] = i;
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
  getline(featurestream,field,parameters_->dataDelimiter);

  //Next read the first row, which should contain the column headers
  getline(featurestream,row);
  stringstream ss( utils::chomp(row) );
  bool isFeaturesAsRows = true;
  vector<string> columnHeaders;
  while ( getline(ss,field,parameters_->dataDelimiter) ) {

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

    row = utils::chomp(row);

    //Read row from the stream
    ss.clear();
    ss.str("");

    //Read the string back to a stream
    ss << row;

    //Read the next row header from the stream
    getline(ss,field,parameters_->dataDelimiter);
    rowHeaders.push_back(field);

    vector<string> rawVector(nColumns);
    for(size_t i = 0; i < nColumns; ++i) {
      getline(ss,rawVector[i],parameters_->dataDelimiter);
      rawVector[i] = utils::trim(rawVector[i]);
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

    row = utils::chomp(row);

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
      break;
      //cout << "found data header:" << row << endl;
    } else {      //If none of the earlier branches matched, we have a problem
      cerr << "incorrectly formatted ARFF row '" << row << "'" << endl;
      assert(false);
    }
    
  }

  if ( !hasData ) {
    cerr << "Treedata::readARFF() -- could not find @data/@DATA identifier" << endl;
    exit(1);
  }

  if ( !hasRelation ) {
    cerr << "Treedata::readARFF() -- could not find @relation/@RELATION identifier" << endl;
    exit(1);
  }
    
  //Read data row-by-row
  while ( getline(featurestream,row) ) {
    
    row = utils::chomp(row);

    //Comment lines and empty lines are omitted
    if ( row == "" ) {
      continue;
    }
    
    // One sample is stored as row in the matrix
    rawMatrix.push_back( utils::split(row,',') );
    
    if ( rawMatrix.back().size() != nFeatures ) {
      cerr << "Treedata::readARFF() -- sample contains incorrect number of features" << endl;
      exit(1);
    }
    
  }
  
  this->transpose<string>(rawMatrix);
  
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
  if(datadefs::toUpperCase(attributeType) == "NUMERIC" ||
     datadefs::toUpperCase(attributeType) == "REAL" ) {
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
  getline(ss,typeStr,parameters_->headerDelimiter);

  return(  typeStr == "N" );
}

bool Treedata::isValidCategoricalHeader(const string& str) {

  stringstream ss(str);
  string typeStr;
  getline(ss,typeStr,parameters_->headerDelimiter);

  return( typeStr == "C" || typeStr == "B" );
}

bool Treedata::isValidFeatureHeader(const string& str) {
  return( isValidNumericalHeader(str) || isValidCategoricalHeader(str) );
}

size_t Treedata::nFeatures() {
  return( parameters_->useContrasts ? features_.size() / 2 : features_.size() );
}

size_t Treedata::nSamples() {
  return( sampleHeaders_.size() );
}

// WILL BECOME DEPRECATED
num_t Treedata::pearsonCorrelation(size_t featureIdx1, size_t featureIdx2) {
  
  vector<num_t> featureData1,featureData2;

  vector<size_t> sampleIcs = utils::range( this->nSamples() );

  this->getFilteredFeatureDataPair(featureIdx1,featureIdx2,sampleIcs,featureData1,featureData2);

  return( math::pearsonCorrelation(featureData1,featureData2) );

}

size_t Treedata::getFeatureIdx(const string& featureName) {
  
  unordered_map<string,size_t>::const_iterator it( name2idx_.find(featureName) );
  
  // If the feature does not exist, return "end", which is a value that 
  // points to outside the range of valid indices
  if ( it == name2idx_.end() ) {
    return( this->end() );
  }
  return( it->second );
}

string Treedata::getFeatureName(const size_t featureIdx) {
  return( features_.at(featureIdx).name );
}

string Treedata::getSampleName(const size_t sampleIdx) {
  return( sampleHeaders_.at(sampleIdx) );
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


void Treedata::permuteContrasts(distributions::RandInt& randInt) {

  size_t nFeatures = this->nFeatures();
  size_t nSamples = this->nSamples();

  for ( size_t i = nFeatures; i < 2*nFeatures; ++i ) {
    
    vector<size_t> sampleIcs = utils::range( nSamples );

    vector<num_t> filteredData = this->getFilteredFeatureData(i,sampleIcs);
    
    utils::permute(filteredData,randInt);
    //this->permute<num_t>(filteredData);

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

vector<string> Treedata::categories(const size_t featureIdx) {
  
  vector<string> categories;

  if( this->isFeatureNumerical(featureIdx) ) {
    return( categories );
  }
 
  for ( map<num_t,string>::const_iterator it( features_[featureIdx].backMapping.begin() ) ; it != features_[featureIdx].backMapping.end(); ++it ) {
    categories.push_back(it->second);
  }

  return( categories );

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

/*
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
*/


void Treedata::bootstrapFromRealSamples(distributions::RandInt& randInt,
					const bool withReplacement, 
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
      ics[sampleIdx] = allIcs[ randInt() % nRealSamples ];
    }
  } else {  //If sampled without replacement...
    vector<size_t> foo = utils::range(nRealSamples);
    utils::permute(foo,randInt);
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

// !! Correctness, Inadequate Abstraction: kill this method with fire. Refactor, REFACTOR, _*REFACTOR*_.
num_t Treedata::numericalFeatureSplit(const size_t targetIdx,
				      const size_t featureIdx,
				      const size_t minSamples,
				      vector<size_t>& sampleIcs_left,
				      vector<size_t>& sampleIcs_right,
				      num_t& splitValue) {

  num_t DI_best = 0.0;

  vector<num_t> tv,fv;

  sampleIcs_left.clear();

  this->getFilteredAndSortedFeatureDataPair3(targetIdx,featureIdx,sampleIcs_right,tv,fv);

  size_t n_tot = fv.size();
  size_t n_right = n_tot;
  size_t n_left = 0;

  if(n_tot < 2 * minSamples) {
    DI_best = 0.0;
    return( DI_best );
  }

  int bestSplitIdx = -1;

  //If the target is numerical, we use the incremental squared error formula
  if ( this->isFeatureNumerical(targetIdx) ) {

    // We start with all samples on the left branch
    num_t mu_tot = math::mean(tv);
    num_t mu_left = 0.0;
    num_t mu_right = mu_tot;

    // Make sure the squared error didn't become corrupted by NANs
    assert( !datadefs::isNAN(mu_tot) );

    // Add samples one by one from left to right until we hit the
    // minimum allowed size of the branch
    for( size_t i = 0; i < n_tot - minSamples; ++i ) {

      // Add n'th sample tv[i] from left to right and update
      // mean and squared error
      ++n_left;
      --n_right;
      mu_left  += ( tv[i] - mu_left  ) / n_left;
      mu_right -= ( tv[i] - mu_right ) / n_right;

      // If the sample is repeated and we can continue, continue
      if ( n_left < minSamples || (n_left < n_tot - minSamples && fv[ i + 1 ] == fv[ i ]) ) {
        continue;
      }
      
      // If the current split point yields a better split than the best,
      // update DI_best and bestSplitIdx
      num_t DI = math::deltaImpurity_regr(mu_tot,n_tot,mu_left,n_left,mu_right,n_right);
      //cout << "(" << n_left << "," << tv[i] << "," << DI << ")"; 
      if (  DI > DI_best ) {

        bestSplitIdx = i;
        DI_best = DI;

      }

    }

    // Make sure there were no NANs to corrupt the results
    assert( !datadefs::isNAN(DI_best) );


  } else { // Otherwise we use the iterative gini index formula to update impurity scores while we traverse "right"

    map<num_t,size_t> freq_right;
    size_t sf_right = 0;
    
    for ( size_t i = 0; i < n_tot; ++i ) {
      math::incrementSquaredFrequency(tv[i],freq_right,sf_right);
    }

    size_t sf_tot = sf_right;

    map<num_t,size_t> freq_left;
    size_t sf_left = 0;

    // Add samples one by one from right to left until we hit the
    // minimum allowed size of the branch
    for( size_t i = 0; i < n_tot - minSamples; ++i ) {

      // Add n'th sample tv[i] from right to left and update
      // mean and squared frequency
      math::incrementSquaredFrequency(tv[i],freq_left,sf_left);
      ++n_left;

      math::decrementSquaredFrequency(tv[i],freq_right,sf_right);
      --n_right;

      // If we have repeated samples and can continue, continue
      if ( n_left < minSamples || (n_left < n_tot - minSamples && fv[ i + 1 ] == fv[ i ]) ) {
        continue;
      }

      // If the split point "i-1" yields a better split than the previous one,
      // update se_best and bestSplitIdx
      num_t DI = math::deltaImpurity_class(sf_tot,n_tot,sf_left,n_left,sf_right,n_right);
      //cout << " tv=" << features_[targetIdx].backMapping[tv[i]] << " fv=" << fv[i] << " nl=" << n_left << " sfl=" << sf_left << " nr=" << n_right << " sfr=" << sf_right << " nt=" << n_tot << " sft=" << sf_tot << " DI=" << DI << endl;
      
      if ( DI > DI_best ) {
        bestSplitIdx = i;
	DI_best = DI;
      }

    }
    
    assert( !datadefs::isNAN(DI_best) );

  }


  if ( bestSplitIdx == -1 ) {
    DI_best = 0.0;
    return( DI_best );
  }

  splitValue = fv[bestSplitIdx];
  n_left = bestSplitIdx + 1;
  sampleIcs_left.resize(n_left);

  for(size_t i = 0; i < n_left; ++i) {
    sampleIcs_left[i] = sampleIcs_right[i];
  }
  sampleIcs_right.erase(sampleIcs_right.begin(),sampleIcs_right.begin() + n_left);
  n_right = sampleIcs_right.size();

  assert(n_left + n_right == n_tot);

  //cout << "N : " << n_left << " <-> " << n_right << " : fitness " << splitFitness << endl;

  return( DI_best );
  
}

// !! Inadequate Abstraction: Refactor me.
num_t Treedata::categoricalFeatureSplit(const size_t targetIdx,
					const size_t featureIdx,
					const size_t minSamples,
					vector<size_t>& sampleIcs_left,
					vector<size_t>& sampleIcs_right,
					set<num_t>& splitValues_left,
					set<num_t>& splitValues_right) {

  num_t DI_best = 0.0;

  vector<num_t> tv,fv;

  //cout << " -- sampleIcs_right.size() = " << sampleIcs_right.size();

  sampleIcs_left.clear();
  this->getFilteredFeatureDataPair(targetIdx,featureIdx,sampleIcs_right,tv,fv);

  //cout << " => " << sampleIcs_right.size() << endl;

  // Map all feature categories to the corresponding samples and represent it as map. The map is used to assign samples to left and right branches
  map<num_t,vector<size_t> > fmap_right;
  map<num_t,vector<size_t> > fmap_left;

  size_t n_tot = 0;
  datadefs::map_data(fv,fmap_right,n_tot);
  size_t n_right = n_tot;
  size_t n_left = 0;

  if(n_tot < 2 * minSamples) {
    DI_best = 0.0;
    return( DI_best );
  }

  if ( this->isFeatureNumerical(targetIdx) ) {

    num_t mu_tot = math::mean(tv);
    num_t mu_right = mu_tot;
    num_t mu_left = 0.0;

    while ( fmap_right.size() > 1 ) {

      map<num_t,vector<size_t> >::iterator it_best( fmap_right.end() );

      // We test each category one by one and see if the fitness becomes improved
      for ( map<num_t,vector<size_t> >::iterator it( fmap_right.begin() ); it != fmap_right.end() ; ++it ) {
	
        //cout << "Testing to split with feature '" << treedata->getRawFeatureData(featureIdx,it->first) << "'" << endl;
	
        // Take samples from right and put them left
        //cout << "from right to left: [";
	size_t n_left_c = n_left;
	size_t n_right_c = n_right;
	num_t mu_left_c = mu_left;
	num_t mu_right_c = mu_right;
	
        for(size_t i = 0; i < it->second.size(); ++i) {
          //cout << " " << it->second[i];
	  
	  ++n_left_c;
	  --n_right_c;
	  mu_left_c  += ( tv[ it->second[i] ] - mu_left_c  ) / n_left_c;
	  mu_right_c -= ( tv[ it->second[i] ] - mu_right_c ) / n_right_c;
	  
        }
        //cout << " ]" << endl;
	
        //If the split reduces impurity even further, save the point
	num_t DI = math::deltaImpurity_regr(mu_tot,n_tot,mu_left_c,n_left_c,mu_right_c,n_right_c);
        if ( DI > DI_best ) { //&& n_left >= minSamples && n_right >= minSamples )
	  
          it_best = it;
          DI_best = DI;
        }

      }

      // After testing all categories,
      // if we couldn't find any split that would reduce impurity,
      // we'll exit the loop
      if ( it_best == fmap_right.end() ) {
        //cout << " -- STOP --" << endl;
        break;
      }
      
      // Otherwise move samples from right to left
      for(size_t i = 0; i < it_best->second.size(); ++i) {
        //cout << " " << it->second[i];

	++n_left;
	--n_right;
	mu_left  += ( tv[ it_best->second[i] ] - mu_left  ) / n_left;
	mu_right -= ( tv[ it_best->second[i] ] - mu_right ) / n_right;

      }

      // Update the maps
      fmap_left.insert( *it_best );
      fmap_right.erase( it_best->first );

    }

    // Calculate the final split fitness
    //splitFitness = this->getNumericalSplitFitness(se_tot,se_best);

  } else {

    map<num_t,size_t> freq_left,freq_right;
    size_t sf_left = 0;
    size_t sf_right = 0;

    for( size_t i = 0; i < n_tot; ++i ) {
      math::incrementSquaredFrequency(tv[i], freq_right, sf_right);
    }

    size_t sf_tot = sf_right;

    while ( fmap_right.size() > 1 ) {

      map<num_t,vector<size_t> >::iterator it_best( fmap_right.end() );
      //cout << "There are " << fmap_right.size() << " categories on right" << endl;

      // We test each category one by one and see if the fitness becomes improved
      for ( map<num_t,vector<size_t> >::iterator it( fmap_right.begin() ); it != fmap_right.end() ; ++it ) {

        //cout << "Testing to split with feature '" << treedata->getRawFeatureData(featureIdx,it->first) << "'" << endl;

        // Take samples from right and put them left
        //cout << "from right to left: [";
        for(size_t i = 0; i < it->second.size(); ++i) {

          // Add sample to left
          ++n_left;
	  math::incrementSquaredFrequency(tv[ it->second[i] ], freq_left, sf_left);

          // Remove sample from right
          --n_right;
	  math::decrementSquaredFrequency(tv[ it->second[i] ], freq_right, sf_right);

        }
        //cout << " ]" << endl;

        //If the impurity becomes reduced even further, save the point
	num_t DI = math::deltaImpurity_class(sf_tot,n_tot,sf_left,n_left,sf_right,n_right);
        if ( DI > DI_best ) { //&& n_left >= minSamples && n_right >= minSamples )

          it_best = it;
	  DI_best = DI;
        }

        // Take samples from left and put them right
        //cout << "From left to right: [";
        for(size_t i = 0; i < it->second.size(); ++i) {

          // Add sample to right
          ++n_right;
	  math::incrementSquaredFrequency(tv[ it->second[i] ], freq_right, sf_right);

          // Remove sample from left
          --n_left;
	  math::decrementSquaredFrequency(tv[ it->second[i] ], freq_left, sf_left);

        }
        //cout << " ]" << endl;

      }

      // After testing all categories,
      // if we couldn't find any split that would reduce impurity,
      // we'll exit the loop
      if ( it_best == fmap_right.end() ) {
        //cout << " -- STOP --" << endl;
        break;
      }

      // Take samples from right and put them left
      for(size_t i = 0; i < it_best->second.size(); ++i) {

        // Add sample to left
        ++n_left;
	math::incrementSquaredFrequency(tv[ it_best->second[i] ], freq_left, sf_left);

        // Remove sample from right
        --n_right;
	math::decrementSquaredFrequency(tv[ it_best->second[i] ], freq_right, sf_right);

      }

      // Update the maps
      fmap_left.insert( *it_best );
      fmap_right.erase( it_best->first );

    }

  }

  if( n_left < minSamples || n_right < minSamples ) {
    DI_best = 0.0;
    return( DI_best );
  }

  // Assign samples and categories on the left. First store the original sample indices
  vector<size_t> sampleIcs = sampleIcs_right;

  assert( n_left + n_right == n_tot );

  // Then populate the left side (sample indices and split values)
  sampleIcs_left.resize(n_left);
  splitValues_left.clear();
  size_t iter = 0;
  for ( map<num_t,vector<size_t> >::const_iterator it(fmap_left.begin()); it != fmap_left.end(); ++it ) {
    for ( size_t i = 0; i < it->second.size(); ++i ) {
      sampleIcs_left[iter] = sampleIcs[it->second[i]];
      ++iter;
    }
    splitValues_left.insert( it->first );
  }
  assert( iter == n_left);
  assert( splitValues_left.size() == fmap_left.size() );

  // Last populate the right side (sample indices and split values)
  sampleIcs_right.resize(n_right);
  splitValues_right.clear();
  iter = 0;
  for ( map<num_t,vector<size_t> >::const_iterator it(fmap_right.begin()); it != fmap_right.end(); ++it ) {
    for ( size_t i = 0; i < it->second.size(); ++i ) {
      sampleIcs_right[iter] = sampleIcs[it->second[i]];
      ++iter;
    }
    splitValues_right.insert( it->first );
  }

  assert( iter == n_right );
  assert( splitValues_right.size() == fmap_right.size() );

  return( DI_best );

}


/*
  void Treedata::getFilteredAndSortedFeatureDataPair(const size_t targetIdx, 
  const size_t featureIdx, 
  vector<size_t>& sampleIcs, 
  vector<num_t>& targetData, 
  vector<num_t>& featureData) {
  
  if ( !features_[featureIdx].isNumerical ) {
  cerr << "Treedata::getFilteredAndSortedDataPair() -- cannot perform for CATEGORICAL features" << endl;
  exit(1);
  }
  
  //targetData.clear();
  //targetData.resize( sampleHeaders_.size(), datadefs::NUM_NAN );
  //featureData.clear();
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
  //num_t tVal = features_[targetIdx].data[*it];
  //num_t fVal = features_[featureIdx].data[*it];
  
  // If the data are non-NA...
  if ( !datadefs::isNAN(features_[featureIdx].data[*it]) && 
  !datadefs::isNAN(features_[targetIdx].data[*it]) ) {
  
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
*/


/*
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
  
  //vector<num_t> targetDataCopy(n);
  //vector<num_t> featureDataCopy(n);
  //vector<size_t> sampleIcsCopy(n);
  
  fill(temp_.multiplicity.begin(),temp_.multiplicity.end(),0);
  //vector<size_t> multiplicity(n, 0);
  
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
  if ( !datadefs::isNAN(fVal) && 
  !datadefs::isNAN(tVal) ) {
  
  // Accumulate real data counter
  ++nReal;
  
  // Extract the ordered position of the sample
  size_t pos = features_[featureIdx].sortOrder[ii];
  ++temp_.multiplicity[pos];
  
  if ( temp_.multiplicity[pos] == 1 ) {
  temp_.featureDataCopy[pos] = fVal;
  temp_.targetDataCopy[pos] = tVal;
  temp_.sampleIcsCopy[pos] = ii;
  
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
  for ( size_t j = 0; j < temp_.multiplicity[i]; ++j ) {
  featureData[iter] = temp_.featureDataCopy[i];
  targetData[iter] = temp_.targetDataCopy[i];
  sampleIcs[iter] = temp_.sampleIcsCopy[i];
  ++iter;
  }
  }
  
  assert(nReal == iter);
  
  }
*/


void Treedata::getFilteredAndSortedFeatureDataPair3(const size_t targetIdx,
						    const size_t featureIdx,
						    vector<size_t>& sampleIcs,
						    vector<num_t>& targetData,
						    vector<num_t>& featureData) {


  featureData = this->getFeatureData(featureIdx,sampleIcs);
  //targetData = this->getFeatureData(targetIdx,sampleIcs);

  bool isIncreasingOrder = true;
  vector<size_t> refIcs;

  utils::filterSort(isIncreasingOrder,featureData,refIcs);
  //datadefs::sortFromRef<num_t>(targetData,refIcs);
  //datadefs::sortFromRef<size_t>(sampleIcs,refIcs);
  
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

  // If the input data is NaN, we return NaN as string 
  if ( datadefs::isNAN(data) ) {
    return( datadefs::STR_NAN );
  }
  
  // If input feature is numerical, we just represent the numeric value as string
  if ( features_[featureIdx].isNumerical ) {
    //stringstream ss;
    //ss << data;
    //return( ss.str() );
    return( utils::num2str(data) );
  } else {
    
    if ( features_[featureIdx].backMapping.find(data) == features_[featureIdx].backMapping.end() ) {
      cerr << "Treedata::getRawFeatureData() -- unknown value to get" << endl;
      exit(1);
    }
    
    return( features_[featureIdx].backMapping[data] );
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
  //this->updateSortOrder(featureIdx);

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
  //features_[featureIdx].sortOrder.clear();

  // The string literal data needs some processing 
  utils::strv2catv(rawFeatureData,
		   features_[featureIdx].data,
		   features_[featureIdx].mapping,
		   features_[featureIdx].backMapping);
}



