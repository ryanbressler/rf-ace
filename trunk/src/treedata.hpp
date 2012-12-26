//treedata.hpp
//
//

#ifndef TREEDATA_HPP
#define TREEDATA_HPP

#include <cstdlib>
#include <map>
#include <fstream>
#include <unordered_map>
#include <unordered_set>

#include "datadefs.hpp"
#include "distributions.hpp"
#include "options.hpp"

using namespace std;
using datadefs::num_t;


class Feature {
public:

  enum Type { NUM, CAT, TXT, UNKNOWN };
  
  vector<num_t> data;
  map<string,num_t> mapping;
  map<num_t,string> backMapping;
  vector<unordered_set<uint32_t> > hashSet;
  string name;
  
  Feature();
  Feature(const vector<num_t>& newData, const string& newName);
  Feature(const vector<string>& newStringData, const string& newName, bool doHash = false);
  ~Feature();
  
  bool isNumerical() const; 
  bool isCategorical() const;
  bool isTextual() const;

  uint32_t getHash(const size_t sampleIdx, const size_t integer) const;
  bool hasHash(const size_t sampleIdx, const uint32_t hashIdx) const;

  num_t entropy() const;
  
private:

  Type type_;
  
};


class Treedata  {
public:

  // Initializes the object 
  Treedata(const vector<Feature>& features, bool useContrasts = false, const vector<string>& sampleHeaders = vector<string>(0));

  // Initializes the object and reads in a data matrix
  // NOTE: will permute the contrasts, which is why it needs the RNG
  Treedata(string fileName, const char dataDelimiter, const char headerDelimiter, const bool useContrasts = false);

  ~Treedata();

  // Returns the number of features
  size_t nFeatures() const;

  // Calculates Pearson Correlation
  // TODO: WILL BECOME OBSOLETE
  num_t pearsonCorrelation(size_t featureidx1, size_t featureidx2);

  // Returns feature index, given the name
  size_t getFeatureIdx(const string& featureName) const;

  // A value denoting the "one-over-last" feature in matrix
  size_t end() const { return( datadefs::MAX_IDX ); }

  // Returns feature name, given the index
  string getFeatureName(const size_t featureIdx) const;

  num_t getFeatureEntropy(const size_t fatureIdx) const;

  // Returns sample name, given sample index
  string getSampleName(const size_t sampleIdx);

  // Returns the number of samples
  size_t nSamples() const;

  // Returns the number of real samples the feature has
  size_t nRealSamples(const size_t featureIdx);
  size_t nRealSamples(const size_t featureIdx1, const size_t featureIdx2);
  
  // Returns the number of categories a feature has
  size_t nCategories(const size_t featureIdx);

  vector<string> categories(const size_t featureIdx);

  // Returns the number of categories of the feature with the highest cardinality
  size_t nMaxCategories();

  // Prints the treedata matrix in its internal form
  void print();
  void print(const size_t featureIdx);

  vector<num_t> getFeatureData(const size_t featureIdx);
  num_t getFeatureData(const size_t featureIdx, const size_t sampleIdx);
  vector<num_t> getFeatureData(const size_t featureIdx, const vector<size_t>& sampleIcs);

  vector<num_t> getFilteredFeatureData(const size_t featureIdx,
				       vector<size_t>& sampleIcs);

  uint32_t getHash(const size_t featureIdx, const size_t sampleIdx, const size_t integer) const;
  bool hasHash(const size_t featureIdx, const size_t sampleIdx, const uint32_t hashIdx) const;

  vector<num_t> getFeatureWeights() const;

  num_t numericalFeatureSplit(const size_t targetIdx,
			      const size_t featureIdx,
			      const size_t minSamples,
			      vector<size_t>& sampleIcs_left,
			      vector<size_t>& sampleIcs_right,
			      num_t& splitValue);

  num_t categoricalFeatureSplit(const size_t targetIdx,
				const size_t featureIdx,
				const size_t minSamples,
				vector<size_t>& sampleIcs_left,
				vector<size_t>& sampleIcs_right,
				set<num_t>& splitValues_left,
				set<num_t>& splitValues_right);

  num_t textualFeatureSplit(const size_t targetIdx,
			    const size_t featureIdx,
			    const uint32_t hashIdx,
			    const size_t minSamples,
			    vector<size_t>& sampleIcs_left,
			    vector<size_t>& sampleIcs_right);
  
  void getFilteredFeatureDataPair(const size_t featureIdx1, 
				  const size_t featureIdx2, 
				  vector<size_t>& sampleIcs, 
				  vector<num_t>& featureData1, 
				  vector<num_t>& featureData2);
  
  /*
    void getFilteredAndSortedFeatureDataPair(const size_t targetIdx, 
    const size_t featureIdx, 
    vector<size_t>& sampleIcs, 
    vector<num_t>& targetData, 
    vector<num_t>& featureData);
        
    void getFilteredAndSortedFeatureDataPair2(const size_t targetIdx,
    const size_t featureIdx,
    vector<size_t>& sampleIcs,
    vector<num_t>& targetData,
    vector<num_t>& featureData);
  */    
    
  void getFilteredAndSortedFeatureDataPair3(const size_t targetIdx,
					    const size_t featureIdx,
					    vector<size_t>& sampleIcs,
					    vector<num_t>& targetData,
					    vector<num_t>& featureData);
  


  string getRawFeatureData(const size_t featureIdx, const size_t sampleIdx);
  string getRawFeatureData(const size_t featureIdx, const num_t data);
  vector<string> getRawFeatureData(const size_t featureIdx);
  
  // Generates a bootstrap sample from the real samples of featureIdx. Samples not in the bootstrap sample will be stored in oob_ics,
  // and the number of oob samples is stored in noob.
  void bootstrapFromRealSamples(distributions::Random* random,
				const bool withReplacement, 
                                const num_t sampleSize, 
                                const size_t featureIdx, 
                                vector<size_t>& ics, 
                                vector<size_t>& oobIcs);

  void createContrasts();
  void permuteContrasts(distributions::Random* random);

  bool isFeatureNumerical(const size_t featureIdx) const;
  bool isFeatureCategorical(const size_t featureIdx) const;
  bool isFeatureTextual(const size_t featureIdx) const;

  void replaceFeatureData(const size_t featureIdx, const vector<num_t>& featureData);
  void replaceFeatureData(const size_t featureIdx, const vector<string>& rawFeatureData);

  
#ifndef TEST__
private:
#endif
  
  enum FileType {UNKNOWN, AFM, ARFF};

  void readFileType(string& fileName, FileType& fileType);

  void readAFM(ifstream& featurestream, 
	       vector<vector<string> >& rawMatrix, 
	       vector<string>& featureHeaders, 
	       vector<string>& sampleHeaders, 
	       vector<Feature::Type>& featureTypes,
	       const char dataDelimiter,
	       const char headerDelimiter);
  
  void readARFF(ifstream& featurestream, 
		vector<vector<string> >& rawMatrix, 
		vector<string>& featureHeaders, 
		vector<Feature::Type>& featureTypes);

  void parseARFFattribute(const string& str, string& attributeName, bool& isFeatureNumerical);

  bool isValidNumericalHeader(const string& str, const char headerDelimiter);
  bool isValidCategoricalHeader(const string& str, const char headerDelimiter);
  bool isValidTextHeader(const string& str, const char headerDelimiter);
  bool isValidFeatureHeader(const string& str, const char headerDelimiter);


  template <typename T> void transpose(vector<vector<T> >& mat);

  bool useContrasts_;
  
  vector<Feature> features_;
  vector<string> sampleHeaders_;

  unordered_map<string,size_t> name2idx_;
  
};

#endif
