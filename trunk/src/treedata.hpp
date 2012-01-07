//treedata.hpp
//
//

#ifndef TREEDATA_HPP
#define TREEDATA_HPP

//#define __IBMCPP_TR1__
//#include <boost/unordered_map.hpp>

#include <cstdlib>
#include <map>
#include <fstream>
//#include <tr1/unordered_map>

#include "datadefs.hpp"
#include "mtrand.h"

using namespace std;
using datadefs::num_t;

class Treedata  {
public:

  // Initializes the object and reads in a data matrix
  Treedata(string fileName, char dataDelimiter, char headerDelimiter);
  Treedata(Treedata& treedata);
  Treedata(Treedata& treedata, const vector<size_t>& featureIcs);
  ~Treedata();
  
  void keepFeatures(const vector<size_t>& featureIcs);

  size_t nFeatures();

  num_t pearsonCorrelation(size_t featureidx1, size_t featureidx2);

  size_t getFeatureIdx(const string& featureName);
  string getFeatureName(const size_t featureIdx);
  string getSampleName(const size_t sampleIdx);

  // Returns the number of samples
  size_t nSamples();

  // Returns the number of real samples the target (resp. any feature) has
  size_t nRealSamples(const size_t featureIdx);
  size_t nRealSamples(const size_t featureIdx1, const size_t featureIdx2);
  
  // Returns the number of categories a feature has
  size_t nCategories(const size_t featureIdx);

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

  void getFilteredFeatureDataPair(const size_t featureIdx1, 
				  const size_t featureIdx2, 
				  vector<size_t>& sampleIcs, 
				  vector<num_t>& featureData1, 
				  vector<num_t>& featureData2);

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

#ifndef TEST__
private:
#endif

  struct Feature {
    vector<num_t> data;
    vector<size_t> sortOrder;
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

  //NOTE: original contents in ics will be replaced.
  void permute(vector<size_t>& ics);

  //Permutes data.
  void permute(vector<num_t>& data);

  // A helper function that creates sort indices for the feature for fast lookup
  void updateSortOrder(const size_t featureIdx);
  
  template <typename T> void transpose(vector<vector<T> >& mat);
  
  char dataDelimiter_;
  char headerDelimiter_;
  vector<Feature> features_;
  vector<string> sampleHeaders_;

  map<string,size_t> name2idx_;

  //map<string,size_t> name2idxHashTest_;

  MTRand_int32 randomInteger_;
};

#endif
