//densetreedata.hpp
//
//

#ifndef DENSETREEDATA_HPP
#define DENSETREEDATA_HPP

#include <cstdlib>
#include <map>
#include <fstream>
#include <unordered_map>
#include <unordered_set>

#include "datadefs.hpp"
#include "distributions.hpp"
#include "options.hpp"
#include "feature.hpp"
#include "reader.hpp"
#include "treedata.hpp"

using namespace std;
using datadefs::num_t;

class DenseTreeData : public TreeData {
public:

  // Initializes the object 
  DenseTreeData(const vector<Feature>& features, bool useContrasts = false, const vector<string>& sampleHeaders = vector<string>(0));

  // Initializes the object and reads in a data matrix
  DenseTreeData(string fileName, const char dataDelimiter, const char headerDelimiter, const bool useContrasts = false);

  ~DenseTreeData();

  // Reveals the Feature class interface to the user
  const Feature* feature(const size_t featureIdx) const {
    return( &features_[featureIdx] );
  }
  
  // Returns the number of features
  size_t nFeatures() const;
  
  // Returns feature index, given the name
  size_t getFeatureIdx(const string& featureName) const;
  
  // A value denoting the "one-over-last" feature in matrix
  size_t end() const { return( datadefs::MAX_IDX ); }
  
  // Returns sample name, given sample index
  string getSampleName(const size_t sampleIdx);
  
  // Returns the number of samples
  size_t nSamples() const;
  
  vector<num_t> getFeatureWeights() const;
  
  void separateMissingSamples(const size_t featureIdx,
			      vector<size_t>& sampleIcs,
			      vector<size_t>& missingIcs);
  
  num_t numericalFeatureSplit(const size_t targetIdx,
			      const size_t featureIdx,
			      const size_t minSamples,
			      vector<size_t>& sampleIcs_left,
			      vector<size_t>& sampleIcs_right,
			      num_t& splitValue);

  num_t categoricalFeatureSplit(const size_t targetIdx,
				const size_t featureIdx,
				const vector<cat_t>& catOrder,
				const size_t minSamples,
				vector<size_t>& sampleIcs_left,
				vector<size_t>& sampleIcs_right,
				unordered_set<cat_t>& splitValues_left);

  num_t textualFeatureSplit(const size_t targetIdx,
			    const size_t featureIdx,
			    const uint32_t hashIdx,
			    const size_t minSamples,
			    vector<size_t>& sampleIcs_left,
			    vector<size_t>& sampleIcs_right);
    
  //string getRawFeatureData(const size_t featureIdx, const size_t sampleIdx);
  //string getRawFeatureData(const size_t featureIdx, const num_t data);
  //vector<string> getRawFeatureData(const size_t featureIdx);
  
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

  void replaceFeatureData(const size_t featureIdx, const vector<num_t>& featureData);
  void replaceFeatureData(const size_t featureIdx, const vector<string>& rawFeatureData);

  
#ifndef TEST__
private:
#endif
  
  enum FileType {UNKNOWN, AFM, ARFF};

  FileType getFileType(const string& fileName);

  bool isRowsAsSamplesInAFM(Reader& reader, const char headerDelimiter);

  void readAFM(const string& fileName, const char dataDelimiter, const char headerDelimiter);
  //void readARFF(const string& fileName);

  //void parseARFFattribute(const string& str, string& attributeName, bool& isFeatureNumerical);

  bool isValidNumericalHeader(const string& str, const char headerDelimiter);
  bool isValidCategoricalHeader(const string& str, const char headerDelimiter);
  bool isValidTextHeader(const string& str, const char headerDelimiter);
  bool isValidFeatureHeader(const string& str, const char headerDelimiter);

  //template <typename T> void transpose(vector<vector<T> >& mat);

  bool useContrasts_;
  
  vector<Feature> features_;
  vector<string> sampleHeaders_;

  unordered_map<string,size_t> name2idx_;
  
};

#endif
