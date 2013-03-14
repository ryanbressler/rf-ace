#include "feature.hpp"

#include "utils.hpp"

Feature::Feature():
  type_(Feature::Type::UNKNOWN) {
}

Feature::Feature(Feature::Type newType, const string& newName, const size_t nSamples):
  type_(newType),
  name_(newName) {

  if ( type_ == Feature::Type::NUM ) {
    numData.resize(nSamples);
    catData.clear();
    txtData.clear();
  } else if ( type_ == Feature::Type::CAT ) {
    numData.clear();
    catData.resize(nSamples);
    txtData.clear();
  } else {
    numData.clear();
    catData.clear();
    txtData.resize(nSamples);
  }

}

void Feature::setNumSampleValue(const size_t sampleIdx, const num_t val) {
  assert( type_ == Feature::Type::NUM );
  numData[sampleIdx] = val;
}

void Feature::setCatSampleValue(const size_t sampleIdx, const cat_t& val) {
  assert( type_ == Feature::Type::CAT );
  catData[sampleIdx] = val;
}

void Feature::setTxtSampleValue(const size_t sampleIdx, const string& str) {
  assert( type_ == Feature::Type::TXT );
  if ( datadefs::isNAN(str) ) {
    txtData[sampleIdx] = utils::hashText("");
  } else {
    txtData[sampleIdx] = utils::hashText(str);
  }
}

cat_t Feature::getCatData(const size_t sampleIdx) const {
  assert(type_ == Feature::Type::CAT);
  return(catData[sampleIdx]);
}

vector<cat_t> Feature::getCatData() const {
  assert(type_ == Feature::Type::CAT);
  return(catData);
}

vector<cat_t> Feature::getCatData(const vector<size_t>& sampleIcs) const {
  assert(type_ == Feature::Type::CAT);
  vector<cat_t> data(sampleIcs.size());
  for ( size_t i = 0; i < sampleIcs.size(); ++i ) {
    data[i] = catData[sampleIcs[i]];
  }
  return(data);
}

num_t Feature::getNumData(const size_t sampleIdx) const {
  assert(type_ == Feature::Type::NUM);
  return(numData[sampleIdx]);
}

vector<num_t> Feature::getNumData() const {
  assert(type_ == Feature::Type::NUM);
  return(numData);
}

vector<num_t> Feature::getNumData(const vector<size_t>& sampleIcs) const {
  assert(type_ == Feature::Type::NUM);
  vector<num_t> data(sampleIcs.size());
  for ( size_t i = 0; i < sampleIcs.size(); ++i ) {
    data[i] = numData[sampleIcs[i]];
  }
  return(data);
}

unordered_set<uint32_t> Feature::getTxtData(const size_t sampleIdx) const {
  assert(type_ == Feature::Type::TXT);
  return(txtData[sampleIdx]);
}

Feature::Feature(const vector<num_t>& newNumData, const string& newName):
  type_(Feature::Type::NUM),
  name_(newName) {
  numData = newNumData;
}

Feature::Feature(const vector<cat_t>& newCatData, const string& newName):
  type_(Feature::Type::CAT),
  name_(newName) {
  catData = newCatData;
  }

Feature::Feature(const vector<string>& newTxtData, const string& newName, const bool doHash):
  type_(Feature::Type::TXT),
  name_(newName) {
  
  assert(doHash);

  size_t nSamples = newTxtData.size();
  
  txtData.resize(nSamples);
  
  for ( size_t i = 0; i < nSamples; ++i ) {
    if ( datadefs::isNAN(newTxtData[i]) ) {
      txtData[i] = utils::hashText("");
    } else {
      txtData[i] = utils::hashText(newTxtData[i]);  
    }
  }
  
}

Feature::~Feature() { }

bool Feature::isNumerical() const {
  return( type_ == Feature::Type::NUM ? true : false );
}

bool Feature::isCategorical() const {
  return( type_ == Feature::Type::CAT ? true : false );
}

bool Feature::isTextual() const {
  return( type_ == Feature::Type::TXT ? true : false );
}

bool Feature::isMissing(const size_t sampleIdx) const {
  switch (type_) {
  case NUM:
    return( datadefs::isNAN(numData[sampleIdx]) );
  case CAT:
    return( datadefs::isNAN(catData[sampleIdx]) );
  case TXT:
    return( datadefs::isNAN(txtData[sampleIdx]) );
  case UNKNOWN:
    break;
  } 

  cerr << "Feature::isMissing() -- tried to use with unset feature object!" << endl;
  exit(1);
}

size_t Feature::nSamples() const {
  switch ( type_ ) {
  case NUM:
    return( numData.size() );
  case CAT:
    return( catData.size() );
  case TXT:
    return( txtData.size() );
  case UNKNOWN:
    break;
  }

  cerr << "Feature::nSamples() -- tried to use with unset feature object!" << endl;
  exit(1);
}

string Feature::name() const {
  return( name_ );
}

void Feature::setName(const string& newName) {
  assert( newName.length() > 0 );
  name_ = newName;
}


vector<cat_t> Feature::categories() const {
  
  vector<cat_t> categories;
  
  if( this->isNumerical() || this->isTextual() ) {
    return( categories );
  }

  unordered_set<cat_t> categoriesSet;
  
  for ( size_t i = 0; i < catData.size(); ++i ) {
    if ( !this->isMissing(i) ) {
      categoriesSet.insert(catData[i]);
    }
  }

  categories.resize(categoriesSet.size());

  copy(categoriesSet.begin(),categoriesSet.end(),categories.begin());
  
  return( categories );
  
}


uint32_t Feature::getHash(const size_t sampleIdx, const size_t integer) const {

  assert( type_ == Feature::Type::TXT );

  size_t pos = integer % this->txtData[sampleIdx].size();

  unordered_set<uint32_t>::const_iterator it(this->txtData[sampleIdx].begin());
  for ( size_t i = 0; i < pos; ++i ) {
    it++;
  }

  return(*it);

}

bool Feature::hasHash(const size_t sampleIdx, const uint32_t hashIdx) const {

  return( this->txtData[sampleIdx].find(hashIdx) != this->txtData[sampleIdx].end() );

}

unordered_map<uint32_t,size_t> Feature::getHashKeyFrequency() const {

  size_t nSamples = txtData.size();

  unordered_map<uint32_t,size_t> visitedKeys;
  
  for ( size_t i = 0; i < nSamples; ++i ) {
    for ( unordered_set<uint32_t>::const_iterator it(txtData[i].begin()); it != txtData[i].end(); ++it ) {
      visitedKeys[*it]++;
    }
  }
  
  return(visitedKeys);

}

num_t Feature::entropy() const {

  size_t nSamples = txtData.size();

  unordered_map<uint32_t,size_t> visitedKeys = getHashKeyFrequency();

  unordered_map<uint32_t,size_t>::const_iterator it(visitedKeys.begin());

  num_t entropy = 0.0;

  for ( ; it != visitedKeys.end(); ++it ) {
    num_t f = static_cast<num_t>(it->second) / static_cast<num_t>(nSamples);
    if ( fabs(f) > 1e-5 && fabs(1-f) > 1e-5 ) {
      entropy -= (f * log(f) + (1-f)*log(1-f))/log(2);
    }
  }

  return(entropy);

}

void Feature::removeFrequentHashKeys(num_t fThreshold) {

  size_t nSamples = txtData.size();

  const unordered_map<uint32_t,size_t> visitedKeys = this->getHashKeyFrequency();

  unordered_map<uint32_t,size_t>::const_iterator it(visitedKeys.begin());
  
  for ( ; it != visitedKeys.end(); ++it ) {
    num_t f = static_cast<num_t>(it->second) / static_cast<num_t>(nSamples);
    if ( f > fThreshold ) {
      for ( size_t sampleIdx = 0; sampleIdx < nSamples; ++sampleIdx ) {
	uint32_t hashKey = it->first;
	if ( txtData[sampleIdx].find(hashKey) != txtData[sampleIdx].end() ) {
	  txtData[sampleIdx].erase(it->first);
	}
      }
    }
  }
  
}
