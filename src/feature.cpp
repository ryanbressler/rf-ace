#include "feature.hpp"

#include "utils.hpp"

Feature::Feature():
  type_(Feature::Type::UNKNOWN) {
}

Feature::Feature(Feature::Type newType, const string& newName, const size_t nSamples):
  type_(newType) {

  name = newName;

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
  type_(Feature::Type::NUM) {

  data = newData;
  name = newName;
  mapping.clear();
  backMapping.clear();

}

Feature::Feature(const vector<string>& newStringData, const string& newName, bool doHash):
  type_( doHash ? Feature::Type::TXT : Feature::Type::CAT ) {

  name = newName;

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

num_t Feature::entropy() const {

  size_t nSamples = hashSet.size();

  unordered_map<uint32_t,size_t> visited_keys;

  for ( size_t i = 0; i < nSamples; ++i ) {
    for ( unordered_set<uint32_t>::const_iterator it(hashSet[i].begin()); it != hashSet[i].end(); ++it ) {
      visited_keys[*it]++;
    }
  }

  unordered_map<uint32_t,size_t>::const_iterator it(visited_keys.begin());

  num_t entropy = 0.0;

  for ( ; it != visited_keys.end(); ++it ) {
    num_t f = static_cast<num_t>(it->second) / static_cast<num_t>(nSamples);
    if ( fabs(f) > 1e-5 && fabs(1-f) > 1e-5 ) {
      entropy -= f * log(f) + (1-f)*log(1-f);
    }
  }

  return(entropy);

}
