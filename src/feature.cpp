#include "feature.hpp"

#include "utils.hpp"

Feature::Feature():
  type_(Feature::Type::UNKNOWN) {
}

Feature::Feature(Feature::Type newType, const string& newName, const size_t nSamples):
  type_(newType),
  name_(newName) {

  if ( type_ == Feature::Type::NUM || type_ == Feature::Type::CAT ) {
    data.resize(nSamples);
    hashSet.clear();
  } else if ( type_ == Feature::Type::TXT ) {
    data.clear();
    hashSet.resize(nSamples);
  }

  mapping.clear();
  backMapping.clear();

}

void Feature::setNumSampleValue(const size_t sampleIdx, const num_t val) {
  assert( type_ == Feature::Type::NUM );
  data[sampleIdx] = val;
}

void Feature::setCatSampleValue(const size_t sampleIdx, const string& str) {
  assert( type_ == Feature::Type::CAT );

  if ( datadefs::isNAN_STR(str) ) {
    data[sampleIdx] = datadefs::NUM_NAN;
    return;
  }

  map<string,num_t>::iterator it(mapping.find(str));

  if ( it != mapping.end() ) {
    data[sampleIdx] = it->second;
  } else {
    num_t val = static_cast<num_t>(mapping.size());
    mapping[str] = val;
    backMapping[val] = str;
    data[sampleIdx] = val;
  }
}

void Feature::setTxtSampleValue(const size_t sampleIdx, const string& str) {
  assert( type_ == Feature::Type::TXT );
  if ( datadefs::isNAN_STR(str) ) {
    hashSet[sampleIdx] = utils::hashText("");
  } else {
    hashSet[sampleIdx] = utils::hashText(str);
  }
}

Feature::Feature(const vector<num_t>& newData, const string& newName):
  type_(Feature::Type::NUM),
  name_(newName) {

  data = newData;
  mapping.clear();
  backMapping.clear();

}

Feature::Feature(const vector<string>& newStringData, const string& newName, bool doHash):
  type_( doHash ? Feature::Type::TXT : Feature::Type::CAT ),
  name_(newName) {

  if ( type_ == Feature::Type::CAT ) {

    utils::strv2catv(newStringData,data,mapping,backMapping);

  } else {

    size_t nSamples = newStringData.size();

    hashSet.resize(nSamples);

    for ( size_t i = 0; i < nSamples; ++i ) {

      hashSet[i] = utils::hashText(newStringData[i]);

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
  if ( type_ != TXT ) {
    return( datadefs::isNAN(data[sampleIdx]) );
  } else {
    return( hashSet[sampleIdx].size() == 0 );
  }
}

string Feature::name() const {
  return( name_ );
}

void Feature::setName(const string& newName) {
  assert( newName.length() > 0 );
  name_ = newName;
}

vector<string> Feature::categories() const {

  vector<string> categories;
  
  if( this->isNumerical() || this->isTextual() ) {
    return( categories );
  }
  
  for ( map<num_t,string>::const_iterator it( backMapping.begin() ) ; it != backMapping.end(); ++it ) {
    categories.push_back(it->second);
  }
  
  return( categories );
  
}

size_t Feature::nCategories() const {
  return( mapping.size() );
}

uint32_t Feature::getHash(const size_t sampleIdx, const size_t integer) const {

  assert( type_ == Feature::Type::TXT );

  size_t pos = integer % this->hashSet[sampleIdx].size();

  unordered_set<uint32_t>::const_iterator it(this->hashSet[sampleIdx].begin());
  for ( size_t i = 0; i < pos; ++i ) {
    it++;
  }

  return(*it);

}

bool Feature::hasHash(const size_t sampleIdx, const uint32_t hashIdx) const {

  return( this->hashSet[sampleIdx].find(hashIdx) != this->hashSet[sampleIdx].end() );

}

unordered_map<uint32_t,size_t> Feature::getHashKeyFrequency() const {

  size_t nSamples = hashSet.size();

  unordered_map<uint32_t,size_t> visitedKeys;
  
  for ( size_t i = 0; i < nSamples; ++i ) {
    for ( unordered_set<uint32_t>::const_iterator it(hashSet[i].begin()); it != hashSet[i].end(); ++it ) {
      visitedKeys[*it]++;
    }
  }
  
  return(visitedKeys);

}

num_t Feature::entropy() const {

  size_t nSamples = hashSet.size();

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

  size_t nSamples = hashSet.size();

  const unordered_map<uint32_t,size_t> visitedKeys = this->getHashKeyFrequency();

  unordered_map<uint32_t,size_t>::const_iterator it(visitedKeys.begin());
  
  for ( ; it != visitedKeys.end(); ++it ) {
    num_t f = static_cast<num_t>(it->second) / static_cast<num_t>(nSamples);
    if ( f > fThreshold ) {
      for ( size_t sampleIdx = 0; sampleIdx < nSamples; ++sampleIdx ) {
	uint32_t hashKey = it->first;
	if ( hashSet[sampleIdx].find(hashKey) != hashSet[sampleIdx].end() ) {
	  hashSet[sampleIdx].erase(it->first);
	}
      }
    }
  }
  
}
