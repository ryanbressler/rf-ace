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

using namespace std;

Treedata::Treedata(const vector<Feature>& features, const bool useContrasts, const vector<string>& sampleHeaders):
  useContrasts_(useContrasts),
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
    name2idx_[ features_[featureIdx].name() ] = featureIdx;
  }

  assert( nSamples > 0 );

  if ( sampleHeaders_.size() == 0 ) {
    sampleHeaders_.resize(nSamples,"NO_SAMPLE_ID");
  } 

  assert( sampleHeaders_.size() == nSamples );

  if ( useContrasts_ ) {
    this->createContrasts(); // Doubles matrix size
    // this->permuteContrasts(random);
  }
  
}


/**
   Reads a data file into a Treedata object. The data file can be either AFM or ARFF
   NOTE: dataDelimiter and headerDelimiter are used only when the format is AFM, for 
   ARFF default delimiter (comma) is used 
*/
Treedata::Treedata(string fileName, const char dataDelimiter, const char headerDelimiter, const bool useContrasts):
  useContrasts_(useContrasts) {

  FileType fileType = this->getFileType(fileName);

  if ( fileType == AFM ) {
    this->readAFM(fileName,dataDelimiter,headerDelimiter);
  } else if ( fileType == ARFF ) {
    this->readARFF(fileName);
  } else {
    cerr << "ERROR: unknown file type in '" << fileName << "'" << endl;
    exit(1);
  }

  if ( useContrasts_ ) {
    this->createContrasts();
  }
  
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
    features_[i].setName( features_[i].name().append("_CONTRAST") );
    name2idx_[ features_[i].name() ] = i;
  }

}

Treedata::FileType Treedata::getFileType(const string& fileName) {
  
  FileType fileType;
  
  string suffix = utils::tolower( utils::suffix(fileName) );
  
  if ( suffix == "afm" ) {
    fileType = AFM;
  } else if ( suffix == "arff" ) {
    fileType = ARFF;
  } else {
    fileType = AFM;
  }

  return(fileType);
}

bool Treedata::isValidNumericalHeader(const string& str, const char headerDelimiter) {
  if ( str.size() > 1 ) {
    return( str[0] == 'N' && str[1] == headerDelimiter );
  } else {
    return( false );
  }
}

bool Treedata::isValidCategoricalHeader(const string& str, const char headerDelimiter) {
  if ( str.size() > 1 ) {
    return( ( str[0] == 'C' || str[0] == 'B' ) && str[1] == headerDelimiter );
  } else {
    return(false);
  }
}

bool Treedata::isValidTextHeader(const string& str, const char headerDelimiter) {
  if ( str.size() > 1 ) {
    return( str[0] == 'T' && str[1] == headerDelimiter );
  } else {
    return(false);
  }
}

bool Treedata::isValidFeatureHeader(const string& str, const char headerDelimiter) {
  return( isValidNumericalHeader(str,headerDelimiter) || isValidCategoricalHeader(str,headerDelimiter) || isValidTextHeader(str,headerDelimiter) );
}

bool Treedata::isRowsAsSamplesInAFM(Reader& reader, const char headerDelimiter) {

  reader.rewind();
  reader.nextLine();
  reader.skipField();

  bool foundSomeFeatures = false;
  bool foundSomeSamples = false;

  while ( !reader.endOfLine() ) {
    string str;
    reader >> str;
    if ( this->isValidFeatureHeader(str,headerDelimiter) ) {
      foundSomeFeatures = true;
    } else {
      foundSomeSamples = true;
    }
  }

  reader.rewind();

  if ( foundSomeFeatures && foundSomeSamples ) {
    cerr << "ERROR: First for of AFM contains both features and sample names!" << endl;
    exit(1);
  }

  return( foundSomeFeatures ? true : false );
  
}

void Treedata::readAFM(const string& fileName, const char dataDelimiter, const char headerDelimiter) {

  Reader reader(fileName,dataDelimiter);

  string numPrefix = string("N") + headerDelimiter;
  string catPrefix = string("C") + headerDelimiter;
  string txtPrefix = string("T") + headerDelimiter;

  if ( this->isRowsAsSamplesInAFM(reader,headerDelimiter) ) { 
    
    size_t nSamples = reader.nLines() - 1;
    
    // Load first line into linefeed and skip the first field (top-left corner) being empty
    reader.nextLine();
    reader.skipField();
    
    // Prepare feature containers and name2idx mapping
    features_.resize(0);
    name2idx_.clear();
    for ( size_t i = 0; ! reader.endOfLine(); ++i ) {
      string featureName; reader >> featureName;
      if ( featureName.substr(0,2) == numPrefix ) {
	features_.push_back( Feature(Feature::Type::NUM,featureName,nSamples) );
      } else if ( featureName.substr(0,2) == catPrefix ) {
	features_.push_back( Feature(Feature::Type::CAT,featureName,nSamples) );
      } else if ( featureName.substr(0,2) == txtPrefix ) {
	features_.push_back( Feature(Feature::Type::TXT,featureName,nSamples) );
      } else {
	cerr << "ERROR reading AFM: unknown feature type for '" << featureName << "'. Are you sure you didn't mean TAFM (Transposed AFM)?" << endl;
	exit(1);
      }
      if ( name2idx_.find(featureName) == name2idx_.end() ) {
	name2idx_[featureName] = i;
      } else {
	cerr << "ERROR reading AFM: duplicate feature name found" << endl;
	exit(1);
      }
    }
    
    assert( reader.endOfLine() );
    
    size_t nFeatures = features_.size();
    
    // Read sample names and data
    sampleHeaders_.resize(nSamples);
    for ( size_t i = 0; i < nSamples; ++i ) {
      reader.nextLine();
      reader >> sampleHeaders_[i];
      for ( size_t j = 0; j < nFeatures; ++j ) {
	if ( features_[j].isNumerical() ) {
	  num_t val; reader >> val;
	  features_[j].setNumSampleValue(i,val);
	} else if ( features_[j].isCategorical() ) {
	  string str; reader >> str;
	  features_[j].setCatSampleValue(i,str);
	} else if ( features_[j].isTextual() ) {
	  string str; reader >> str;
	  features_[j].setTxtSampleValue(i,str);
	}
      }
      assert( reader.endOfLine() );
    }

  } else { 

    size_t nFeatures = reader.nLines() - 1;

    // Load first line into linefeed and skip the first field (top-left corner) being empty
    reader.nextLine();
    reader.skipField();

    sampleHeaders_.clear();
    while ( ! reader.endOfLine() ) {
      string sampleName; reader >> sampleName;
      sampleHeaders_.push_back( sampleName );
    }

    assert( reader.endOfLine() );

    size_t nSamples = sampleHeaders_.size();

    name2idx_.clear();
    for ( size_t i = 0; i < nFeatures; ++i ) {
      reader.nextLine();
      string featureName; reader >> featureName;
      if ( featureName.substr(0,2) == numPrefix ) {
	features_.push_back( Feature(Feature::Type::NUM,featureName,nSamples) );
	for ( size_t j = 0; j < nSamples; ++j ) {
	  num_t val; reader >> val;
	  features_[i].setNumSampleValue(j,val);
	}
      } else if ( featureName.substr(0,2) == catPrefix ) {
	features_.push_back( Feature(Feature::Type::CAT,featureName,nSamples) );
	for ( size_t j = 0; j < nSamples; ++j ) {
	  string str; reader >> str;
	  features_[i].setCatSampleValue(j,str);
	}
      } else if ( featureName.substr(0,2) == txtPrefix ) {
	features_.push_back( Feature(Feature::Type::TXT,featureName,nSamples) );
	for ( size_t j = 0; j < nSamples; ++j ) {
	  string str; reader >> str;
	  features_[i].setTxtSampleValue(j,str);
	}
      } else {
	cerr << "ERROR reading TAFM: unknown feature type for '" << featureName << "'. Are you sure you didn't mean AFM?" << endl;
	exit(1);
      }
      if ( name2idx_.find(featureName) == name2idx_.end() ) {
	name2idx_[featureName] = i;
      } else {
	cerr << "ERROR reading TAFM: duplicate feature name found" << endl;
	exit(1);
      }
    }

    assert( nFeatures == features_.size() );

  }

}

void Treedata::readARFF(const string& fileName) {

  Reader reader(fileName,',');

  cerr << "ERROR: ARFF implementation is missing!" << endl;
  exit(1);
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

size_t Treedata::nFeatures() const {
  return( useContrasts_ ? features_.size() / 2 : features_.size() );
}

size_t Treedata::nSamples() const {
  return( sampleHeaders_.size() );
}

// WILL BECOME DEPRECATED
num_t Treedata::pearsonCorrelation(size_t featureIdx1, size_t featureIdx2) {
  
  if ( ! this->feature(featureIdx1)->isNumerical() || ! this->feature(featureIdx2)->isNumerical() ) {
    return( datadefs::NUM_NAN );
  }

  vector<size_t> sampleIcs = utils::range( this->nSamples() );

  vector<size_t> missingIcs;

  this->separateMissingSamples(featureIdx1,sampleIcs,missingIcs);
  this->separateMissingSamples(featureIdx2,sampleIcs,missingIcs);

  vector<num_t> featureData1 = this->getFeatureData(featureIdx1,sampleIcs);
  vector<num_t> featureData2 = this->getFeatureData(featureIdx2,sampleIcs);

  return( math::pearsonCorrelation(featureData1,featureData2) );

}

size_t Treedata::getFeatureIdx(const string& featureName) const {
  
  unordered_map<string,size_t>::const_iterator it( name2idx_.find(featureName) );
  
  // If the feature does not exist, return "end", which is a value that 
  // points to outside the range of valid indices
  if ( it == name2idx_.end() ) {
    return( this->end() );
  }
  return( it->second );
}

string Treedata::getSampleName(const size_t sampleIdx) {
  return( sampleHeaders_.at(sampleIdx) );
}

vector<num_t> Treedata::getFeatureWeights() const {

  vector<num_t> weights(this->nFeatures(),1.0);

  for ( size_t i = 0; i < this->nFeatures(); ++i ) {
    const Feature* feature = this->feature(i);
    if ( feature->isTextual() ) {
      num_t entropy = feature->entropy();
      cout << "Feature '" << feature->name() << "' is textual and has sqrt(entropy) " << sqrtf(entropy) << endl;
      weights[i] = sqrtf(entropy);
    }
  }

  return(weights);

}

void Treedata::permuteContrasts(distributions::Random* random) {

  size_t nFeatures = this->nFeatures();
  size_t nSamples = this->nSamples();

  for ( size_t i = nFeatures; i < 2*nFeatures; ++i ) {
    
    if ( this->feature(i)->isTextual() ) { continue; }

    vector<size_t> sampleIcs = utils::range( nSamples );
    vector<size_t> missingIcs;
    
    this->separateMissingSamples(i,sampleIcs,missingIcs);

    vector<num_t> filteredData = this->getFeatureData(i,sampleIcs);
    
    utils::permute(filteredData,random);

    for ( size_t j = 0; j < sampleIcs.size(); ++j ) {
      features_[i].data[ sampleIcs[j] ] = filteredData[j];
    }

  }

}

size_t Treedata::nRealSamples(const size_t featureIdx) { 
  
  size_t nRealSamples;
  datadefs::countRealValues( features_[featureIdx].data, nRealSamples );
  return( nRealSamples );

}

size_t Treedata::nRealSamples(const size_t featureIdx1, const size_t featureIdx2) {

  size_t nRealSamples = 0;
  for( size_t i = 0; i < Treedata::nSamples(); ++i ) {
    if( !this->feature(featureIdx1)->isMissing(i) && !this->feature(featureIdx2)->isMissing(i) ) {
      ++nRealSamples;
    }
  }
  return( nRealSamples );
}

void Treedata::bootstrapFromRealSamples(distributions::Random* random,
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
  for(size_t i = 0; i < this->nSamples(); ++i) {
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
      ics[sampleIdx] = allIcs[ random->integer() % nRealSamples ];
    }
  } else {  //If sampled without replacement...
    vector<size_t> foo = utils::range(nRealSamples);
    utils::permute(foo,random);
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
  
void Treedata::separateMissingSamples(const size_t featureIdx,
				      vector<size_t>& sampleIcs,
				      vector<size_t>& missingIcs) {
  
  size_t nReal = 0;
  size_t nMissing = 0;

  missingIcs.resize(sampleIcs.size());

  vector<size_t>::const_iterator it(sampleIcs.begin());

  if ( this->feature(featureIdx)->isTextual() ) {
    for ( ; it != sampleIcs.end(); ++it ) {
      if ( this->feature(featureIdx)->hashSet[*it].size() > 0 ) {
	sampleIcs[nReal++] = *it;
      } else {
	missingIcs[nMissing++] = *it;
      }
    } 
  } else {
    for ( ; it != sampleIcs.end(); ++it ) {
      if ( ! datadefs::isNAN( this->feature(featureIdx)->data[*it] ) ) {
	sampleIcs[nReal++] = *it;
      } else {
	missingIcs[nMissing++] = *it;
      }
    }
  }

  sampleIcs.resize(nReal);
  missingIcs.resize(nMissing);

}


// !! Correctness, Inadequate Abstraction: kill this method with fire. Refactor, REFACTOR, _*REFACTOR*_.
num_t Treedata::numericalFeatureSplit(const size_t targetIdx,
				      const size_t featureIdx,
				      const size_t minSamples,
				      vector<size_t>& sampleIcs_left,
				      vector<size_t>& sampleIcs_right,
				      num_t& splitValue) {

  num_t DI_best = 0.0;

  sampleIcs_left.clear();

  vector<num_t> fv = this->getFeatureData(featureIdx,sampleIcs_right);
  vector<num_t> tv = this->getFeatureData(targetIdx,sampleIcs_right);

  vector<size_t> sortIcs = utils::range(sampleIcs_right.size());
  utils::sortDataAndMakeRef(true,fv,sortIcs);
  utils::sortFromRef(tv,sortIcs);
  utils::sortFromRef(sampleIcs_right,sortIcs);

  size_t n_tot = fv.size();
  //size_t n_right = n_tot;
  size_t n_left = 0;

  if(n_tot < 2 * minSamples) {
    DI_best = 0.0;
    return( DI_best );
  }

  size_t bestSplitIdx = datadefs::MAX_IDX;

  //If the target is numerical, we use the incremental squared error formula
  if ( this->feature(targetIdx)->isNumerical() ) {

    DI_best = utils::numericalFeatureSplitsNumericalTarget(tv,fv,minSamples,bestSplitIdx);

  } else { // Otherwise we use the iterative gini index formula to update impurity scores while we traverse "right"

    DI_best = utils::numericalFeatureSplitsCategoricalTarget(tv,fv,minSamples,bestSplitIdx);

  }


  if ( bestSplitIdx == datadefs::MAX_IDX ) {
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
  size_t n_right = sampleIcs_right.size();

  assert(n_left + n_right == n_tot);

  //cout << "N : " << n_left << " <-> " << n_right << " : fitness " << splitFitness << endl;

  return( DI_best );
  
}

// !! Inadequate Abstraction: Refactor me.
num_t Treedata::categoricalFeatureSplit(const size_t targetIdx,
					const size_t featureIdx,
					const vector<num_t>& catOrder,
					const size_t minSamples,
					vector<size_t>& sampleIcs_left,
					vector<size_t>& sampleIcs_right,
					unordered_set<num_t>& splitValues_left) {

  num_t DI_best = 0.0;

  sampleIcs_left.clear();

  vector<num_t> fv = this->getFeatureData(featureIdx,sampleIcs_right);
  vector<num_t> tv = this->getFeatureData(targetIdx,sampleIcs_right);

  size_t n_tot = fv.size();

  if(n_tot < 2 * minSamples) {
    DI_best = 0.0;
    return( DI_best );
  }
  
  unordered_map<num_t,vector<size_t> > fmap_right(catOrder.size());
  unordered_map<num_t,vector<size_t> > fmap_left(catOrder.size());

  if ( this->feature(targetIdx)->isNumerical() ) {

    DI_best = utils::categoricalFeatureSplitsNumericalTarget2(tv,fv,minSamples,catOrder,fmap_left,fmap_right);

  } else {

    DI_best = utils::categoricalFeatureSplitsCategoricalTarget2(tv,fv,minSamples,catOrder,fmap_left,fmap_right);

  }

  if ( fabs(DI_best) < datadefs::EPS ) {
    return(DI_best);
  }

  // Assign samples and categories on the left. First store the original sample indices
  vector<size_t> sampleIcs = sampleIcs_right;

  // Then populate the left side (sample indices and split values)
  sampleIcs_left.resize(n_tot);
  splitValues_left.clear();
  splitValues_left.reserve(fmap_right.size());
  size_t iter = 0;
  for ( unordered_map<num_t,vector<size_t> >::const_iterator it(fmap_left.begin()); it != fmap_left.end(); ++it ) {
    for ( size_t i = 0; i < it->second.size(); ++i ) {
      sampleIcs_left[iter] = sampleIcs[it->second[i]];
      ++iter;
    }
    splitValues_left.insert( it->first );
  }
  sampleIcs_left.resize(iter);
  //assert( iter == n_left);
  assert( splitValues_left.size() == fmap_left.size() );

  // Last populate the right side (sample indices and split values)
  sampleIcs_right.resize(n_tot);
  //unordered_set<num_t> splitValues_right;
  iter = 0;
  for ( unordered_map<num_t,vector<size_t> >::const_iterator it(fmap_right.begin()); it != fmap_right.end(); ++it ) {
    for ( size_t i = 0; i < it->second.size(); ++i ) {
      sampleIcs_right[iter] = sampleIcs[it->second[i]];
      ++iter;
    }
    //splitValues_right.insert( it->first );
  }
  sampleIcs_right.resize(iter);
  //assert( iter == n_right );
  //assert( splitValues_right.size() == fmap_right.size() );

  return( DI_best );

}

num_t Treedata::textualFeatureSplit(const size_t targetIdx,
				    const size_t featureIdx,
				    const uint32_t hashIdx,
				    const size_t minSamples,
				    vector<size_t>& sampleIcs_left,
				    vector<size_t>& sampleIcs_right) {


  assert(features_[featureIdx].isTextual());

  size_t n_left = 0;
  size_t n_right = 0;
  size_t n_tot = sampleIcs_right.size();

  sampleIcs_left.resize(n_tot);

  num_t DI_best = 0.0;

  if ( this->feature(targetIdx)->isNumerical() ) {
  
    num_t mu_left = 0.0;
    num_t mu_right = 0.0;
    num_t mu_tot = 0.0;

    for ( size_t i = 0; i < n_tot; ++i ) {
      unordered_set<uint32_t>& hs = features_[featureIdx].hashSet[sampleIcs_right[i]];
      num_t x = features_[targetIdx].data[sampleIcs_right[i]];
      if ( hs.find(hashIdx) != hs.end() ) {
	sampleIcs_left[n_left++] = sampleIcs_right[i];
	mu_left += ( x - mu_left ) / n_left;
      } else {
	sampleIcs_right[n_right++] = sampleIcs_right[i];
	mu_right += ( x - mu_right ) / n_right;
      }
      mu_tot += x / n_tot;
    }
    
    DI_best = math::deltaImpurity_regr(mu_tot,n_tot,mu_left,n_left,mu_right,n_right);

  } else {

    unordered_map<num_t,size_t> freq_left,freq_right,freq_tot(n_tot);

    size_t sf_left = 0;
    size_t sf_right = 0;
    size_t sf_tot = 0;

    for ( size_t i = 0; i < sampleIcs_right.size(); ++i ) {
      unordered_set<uint32_t>& hs = features_[featureIdx].hashSet[sampleIcs_right[i]];
      num_t x = features_[targetIdx].data[sampleIcs_right[i]];
      if ( hs.find(hashIdx) != hs.end() ) {
        sampleIcs_left[n_left++] = sampleIcs_right[i];
	math::incrementSquaredFrequency(x,freq_left,sf_left);
      } else {
        sampleIcs_right[n_right++] = sampleIcs_right[i];
	math::incrementSquaredFrequency(x,freq_right,sf_right);
      }
      math::incrementSquaredFrequency(x,freq_tot,sf_tot);
    }

    DI_best = math::deltaImpurity_class(sf_tot,n_tot,sf_left,n_left,sf_right,n_right);

  }

  assert(n_tot == n_left + n_right);

  if ( n_left < minSamples || n_right < minSamples ) {
    return(0.0);
  }

  sampleIcs_left.resize(n_left);
  sampleIcs_right.resize(n_right);
  
  return(DI_best);
  
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
  if ( features_[featureIdx].isNumerical() ) {
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

  features_[featureIdx] = Feature( featureData, features_[featureIdx].name() );

}

void Treedata::replaceFeatureData(const size_t featureIdx, const vector<string>& rawFeatureData) {

  if(rawFeatureData.size() != features_[featureIdx].data.size() ) {
    cerr << "Treedata::replaceFeatureData(string) -- data dimension mismatch" << endl;
    exit(1);
  }

  features_[featureIdx] = Feature( rawFeatureData, features_[featureIdx].name() );

}



