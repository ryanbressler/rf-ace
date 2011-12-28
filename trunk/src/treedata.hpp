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
  //Initializes the object and reads in a data matrix
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

  //Returns the number of samples
  size_t nSamples();

  //Returns the number of real samples the target (resp. any feature) has
  size_t nRealSamples(const size_t featureIdx);
  size_t nRealSamples(const size_t featureIdx1, const size_t featureIdx2);
  
  // Returns the number of categories a feature has
  size_t nCategories(const size_t featureIdx);
  //size_t nCategories(const string& featureName);

  // Returns the number of categories of the feature with the highest cardinality
  size_t nMaxCategories();

  //Prints the treedata matrix in its internal form
  void print();
  void print(const size_t featureIdx);

  vector<num_t> getFeatureData(const size_t featureIdx);
  num_t getFeatureData(const size_t featureIdx, const size_t sampleIdx);
  vector<num_t> getFeatureData(const size_t featureIdx, const vector<size_t>& sampleIcs);

  void getFilteredDataPair(const size_t featureIdx1, const size_t featureIdx2, vector<size_t>& sampleIcs, vector<num_t>& featureData1, vector<num_t>& featureData2);

  void getFilteredAndSortedDataPair(const size_t featureIdx1, const size_t featureIdx2, vector<size_t>& sampleIcs, vector<num_t>& featureData1, vector<num_t>& featureData2);

  //vector<num_t> operator[](size_t featureIdx);
  //vector<num_t> operator[](const string& featureName);

  string getRawFeatureData(const size_t featureIdx, const size_t sampleIdx);
  string getRawFeatureData(const size_t featureIdx, const num_t data);
  
  //string getRawFeatureData(const string& featureName, const size_t sampleIdx);
  //string getRawFeatureData(const string& featureName, const num_t data);

  string dataToRaw(const size_t featureIdx, const num_t data);

  inline void getRandomData(const size_t featureIdx, num_t& data) {data = features_[featureIdx].data[randomInteger_() % sampleHeaders_.size() ]; }
  inline size_t getRandomIndex(const size_t n) { return( randomInteger_() % n ); }

  //map<string,num_t> getDataMapping(const size_t featureIdx);
  //map<num_t,string> getDataBackMapping(const size_t featureIdx);
  
  //Generates a bootstrap sample from the real samples of featureIdx. Samples not in the bootstrap sample will be stored in oob_ics,
  //and the number of oob samples is stored in noob.
  void bootstrapFromRealSamples(const bool withReplacement, 
                                const num_t sampleSize, 
                                const size_t featureIdx, 
                                vector<size_t>& ics, 
                                vector<size_t>& oobIcs);

  void permuteContrasts();

  bool isFeatureNumerical(const size_t featureIdx);
  //bool isFeatureNumerical(const string& featureName);

  // DEPRECATED
  //void impurity(vector<num_t>& data, bool isFeatureNumerical, num_t& impurity, size_t& nreal);

  //WILL DECOME DEPRECATED
#ifndef TEST__  
protected: 
#endif

  //WILL BECOME DEPRECATED
  friend class StochasticForest;

#ifndef TEST__
private:
#endif

  
  struct Feature {
    vector<num_t> data;
    bool isNumerical;
    map<string,num_t> mapping;
    map<num_t,string> backMapping;
    size_t nCategories;
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

  //bool isPositiveInteger(const string& str, size_t& positiveInteger);

  //NOTE: original contents in ics will be replaced.
  void permute(vector<size_t>& ics);

  //Permutes data.
  void permute(vector<num_t>& data);
  
  template <typename T> void transpose(vector<vector<T> >& mat);
  
  char dataDelimiter_;
  char headerDelimiter_;
  vector<Feature> features_;
  vector<string> sampleHeaders_;

  map<string,size_t> name2idx_;

  MTRand_int32 randomInteger_;
};

#endif
