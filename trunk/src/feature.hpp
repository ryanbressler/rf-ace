#ifndef FEATURE_HPP
#define FEATURE_HPP

#include <cstdlib>
#include <vector>
#include <unordered_set>
#include <map>
#include <string>

#include "datadefs.hpp"

using namespace std;
using datadefs::num_t;
using datadefs::cat_t;

class Feature {
public:

  enum Type { NUM, CAT, TXT, UNKNOWN };

  vector<num_t> numData;
  vector<cat_t> catData;
  vector<unordered_set<uint32_t> > txtData;

  Feature();
  Feature(Type newType, const string& newName, const size_t nSamples);
  Feature(const vector<num_t>& newNumData, const string& newName);
  Feature(const vector<cat_t>& newCatData, const string& newName);
  Feature(const vector<string>& newTxtData, const string& newName, const bool doHash);
  ~Feature();

  void setNumSampleValue(const size_t sampleIdx, const num_t   val);
  void setCatSampleValue(const size_t sampleIdx, const cat_t&  val);
  void setTxtSampleValue(const size_t sampleIdx, const string& str);

  num_t getNumData(const size_t sampleIdx) const;
  vector<num_t> getNumData() const;
  vector<num_t> getNumData(const vector<size_t>& sampleIcs) const;

  cat_t getCatData(const size_t sampleIdx) const;
  vector<cat_t> getCatData() const;
  vector<cat_t> getCatData(const vector<size_t>& sampleIcs) const;

  unordered_set<uint32_t> getTxtData(const size_t sampleIdx) const;

  bool isNumerical() const;
  bool isCategorical() const;
  bool isTextual() const;

  bool isMissing(const size_t sampleIdx) const;

  size_t nSamples() const;
  size_t nRealSamples() const;

  string name() const;
  void setName(const string& newName);
  
  vector<cat_t> categories() const;

  uint32_t getHash(const size_t sampleIdx, const size_t integer) const;
  bool hasHash(const size_t sampleIdx, const uint32_t hashIdx) const;

  num_t entropy() const;

  unordered_map<uint32_t,size_t> getHashKeyFrequency() const;

  void removeFrequentHashKeys(const num_t fThreshold);

#ifndef TEST__
private:
#endif

  Type type_;
  string name_;

};


#endif
