//treedata.hpp
//
//

#ifndef TREEDATA_HPP
#define TREEDATA_HPP

#include <cstdlib>
#include <map>
#include <fstream>

#include "datadefs.hpp"
#include "mtrand.h"

using namespace std;
using datadefs::num_t;

class Treedata  {
public:

  // Initializes the object and reads in a data matrix
  Treedata(string fileName, char dataDelimiter, char headerDelimiter, int seed = -1 );

  ~Treedata();

  // Takes a set of features that are to be retained in the Treedata object.
  // Others will be removed.
  void whiteList(const set<string>& featureNames);

  // Takes a set of features that are to be removed from the Treedata object.
  // Others will be retained.
  void blackList(const set<string>& featureNames);

  // Takes a boolean vector of flags that mark which feature is to be removed (0) and
  // which to be retained (1) in the Treedata object.
  void whiteList(const vector<bool>& featureIcs);

  // Returns the number of features
  size_t nFeatures();

  // Calculates Pearson Correlation
  // TODO: WILL BECOME OBSOLETE
  num_t pearsonCorrelation(size_t featureidx1, size_t featureidx2);

  // Returns feature index, given the name
  size_t getFeatureIdx(const string& featureName);

  // Returns feature name, given the index
  string getFeatureName(const size_t featureIdx);

  // Returns sample name, given sample index
  string getSampleName(const size_t sampleIdx);

  // Returns the number of samples
  size_t nSamples();

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

  // !! WILL BECOME OBSOLETE
  num_t getCategoricalSplitFitness(const num_t sf_tot,
                                   const num_t nsf,
                                   const size_t n);

  // !! WILL BECOME OBSOLETE
  num_t getNumericalSplitFitness(const num_t se_tot,
				 const num_t se_best);

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
  
  inline void getRandomData(const size_t featureIdx, num_t& data) {data = features_[featureIdx].data[randomInteger_() % sampleHeaders_.size() ]; }

  inline size_t getRandomIndex(const size_t n) { return( randomInteger_() % n ); }

  // 2^32 == 4294967296 == upper bound for randomInteger_
  inline num_t getRandomUnif() { return( static_cast<num_t>( 1.0 * randomInteger_() / 4294967296. ) ); }
  
  // Generates a bootstrap sample from the real samples of featureIdx. Samples not in the bootstrap sample will be stored in oob_ics,
  // and the number of oob samples is stored in noob.
  void bootstrapFromRealSamples(const bool withReplacement, 
                                const num_t sampleSize, 
                                const size_t featureIdx, 
                                vector<size_t>& ics, 
                                vector<size_t>& oobIcs);

  void permuteContrasts();

  bool isFeatureNumerical(const size_t featureIdx);

  void replaceFeatureData(const size_t featureIdx, const vector<num_t>& featureData);
  void replaceFeatureData(const size_t featureIdx, const vector<string>& rawFeatureData);

  template <typename T> 
  void permute(vector<T>& data) {

    // Data size
    size_t n = data.size();
    
    // Mapping indices
    vector<size_t> ics(n);
    
    // Permute indices
    for (size_t i = 0; i < n; ++i) {
      size_t j = randomInteger_() % (i + 1);
      ics[i] = ics[j];
      ics[j] = i;
    }
    
    // Re-map data based on permuted indices
    for(size_t i = 0; i < n; ++i) {
      T temp = data[i];
      data[i] = data[ics[i]];
      data[ics[i]] = temp;
    }

  }

#ifndef TEST__
private:
#endif

  struct Feature {
    vector<num_t> data;
    //vector<size_t> sortOrder;
    bool isNumerical;
    map<string,num_t> mapping;
    map<num_t,string> backMapping;
    string name;
  };
  
  enum FileType {UNKNOWN, AFM, ARFF};

  void readFileType(string& fileName, FileType& fileType);

  void readAFM(ifstream& featurestream, 
	       vector<vector<string> >& rawMatrix, 
	       vector<string>& featureHeaders, 
	       vector<string>& sampleHeaders, 
	       vector<bool>& isFeatureNumerical);
  
  void readARFF(ifstream& featurestream, 
		vector<vector<string> >& rawMatrix, 
		vector<string>& featureHeaders, 
		vector<bool>& isFeatureNumerical);

  void parseARFFattribute(const string& str, string& attributeName, bool& isFeatureNumerical);

  bool isValidNumericalHeader(const string& str);
  bool isValidCategoricalHeader(const string& str);
  bool isValidFeatureHeader(const string& str);

  // NOTE: original contents in ics will be replaced.
  void permute(vector<size_t>& ics);

  // Permutes data.
  void permute(vector<num_t>& data);

  // A helper function that creates sort indices for the feature for fast lookup
  //void updateSortOrder(const size_t featureIdx);
  
  template <typename T> void transpose(vector<vector<T> >& mat);
  
  char dataDelimiter_;
  char headerDelimiter_;
  vector<Feature> features_;
  vector<string> sampleHeaders_;

  map<string,size_t> name2idx_;

  MTRand_int32 randomInteger_;
  
};

#endif
