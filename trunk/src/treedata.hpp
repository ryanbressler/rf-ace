//treedata.hpp
//
//

#ifndef TREEDATA_HPP
#define TREEDATA_HPP

#include <cstdlib>
#include <map>
#include <fstream>
#include "datadefs.hpp"
#include "node.hpp"
#include "mtrand.h"

using namespace std;
//using datadefs::cat_t;
using datadefs::num_t;

class Treedata 
{
public:
  //Initializes the object and reads in a data matrix
  Treedata(string fileName);
  ~Treedata();

  //Returns the number of features
  size_t nfeatures();

  num_t corr(size_t featureidx1, size_t featureidx2);

  string get_featureheader(size_t featureidx);
  string get_targetheader();

  //Returns the number of samples
  size_t nsamples();

  //Returns the number of real samples the target (resp. any feature) has
  size_t nrealsamples();
  size_t nrealsamples(size_t featureidx);

  //Prints the treedata matrix in its internal form
  void print();
  void print(const size_t featureidx);

  void getFeatureData(const size_t featureIdx, vector<num_t>& data);
  void getFeatureData(const size_t featureIdx, const size_t sampleIdx, num_t& data);
  void getFeatureData(const size_t featureIdx, const vector<size_t>& sampleIcs, vector<num_t>& data);

  void getContrastData(const size_t featureIdx, vector<num_t>& data);
  void getContrastData(const size_t featureIdx, const size_t sampleIdx, num_t& data);
  void getContrastData(const size_t featureIdx, const vector<size_t>& sampleIcs, vector<num_t>& data);

  inline size_t sampleRandomIdx(const size_t n) { return(irand_() % n); }
  num_t sampleAtRandom(size_t featureIdx);

protected: 

  friend class Randomforest;
  friend class GBT;

  void killFeature(const size_t featureIdx);

  //Permutes integers in range 0,1,...,(ics.size()-1). 
  //NOTE: original contents in ics will be replaced.  
  void permute(vector<size_t>& ics);

  //Permutes data in x.
  void permute(vector<num_t>& x);

  //Generates a bootstrap sample. Samples not in the bootstrap sample will be stored in oob_ics, 
  //and the number of oob samples is stored in noob.
  //NOTE: ics.size() will depend on the number of non-NaN values the current target has
  void bootstrap(vector<size_t>& ics, vector<size_t>& oob_ics);

  //Selects target feature. 
  //NOTE: data will be sorted with respect to the target.
  void selectTarget(size_t targetidx);

  size_t getTarget();
  
  void permuteContrasts();

  bool isFeatureNumerical(size_t featureIdx);

  //size_t nrealvalues();
  //size_t nrealvalues(size_t featureidx);

  //void remove_nans(size_t featureidx, vector<size_t>& sampleics, size_t& nreal);
  

  //Given feature, finds and returns the optimal split point wrt. sampleics. 
  //Samples branching left and right will be stored in sampleics_left (resp. right)
  /*
    void split_target(size_t featureidx,
    const size_t min_split,
    vector<size_t>& sampleics,
    vector<size_t>& sampleics_left,
    vector<size_t>& sampleics_right,
    num_t& splitvalue,
    set<num_t>& values_left);
  */
    
  num_t splitFitness(vector<num_t> const& data,
		     bool const& isFeatureNumerical,
		     size_t const& minSplit,
		     vector<size_t> const& sampleIcs_left,
		     vector<size_t> const& sampleIcs_right);
  
  void impurity(vector<num_t>& data, bool isFeatureNumerical, num_t& impurity, size_t& nreal);
  



private:

  enum FileType {UNKNOWN, AFM, ARFF};

  void readFileType(string& fileName, FileType& fileType);

  void readAFM(ifstream& featurestream, vector<vector<string> >& rawMatrix, vector<string>& featureHeaders, vector<bool>& isFeatureNumerical);
  void readARFF(ifstream& featurestream, vector<vector<string> >& rawMatrix, vector<string>& featureHeaders, vector<bool>& isFeatureNumerical);

  void parseARFFattribute(const string& str, string& attributeName, bool& isFeatureNumerical);

  bool isValidNumericalHeader(const string& str);
  bool isValidCategoricalHeader(const string& str);
  bool isValidFeatureHeader(const string& str);

  /*
    void getFeatureData(size_t featureIdx, num_t& data);
    void getFeatureData(size_t featureIdx, size_t sampleIdx, num_t& data);
    void getFeatureData(size_t featureIdx, vector<size_t>& sampleIcs, vector<num_t>& data);
    
    void getContrastData(size_t featureIdx, num_t& data);
    void getContrastData(size_t featureIdx, size_t sampleIdx, num_t& data);
    void getContrastData(size_t featureIdx, vector<size_t>& sampleIcs, vector<num_t>& data);
    
    inline size_t sampleRandomIdx(const size_t n) { return(irand_() % n); }
    num_t sampleAtRandom(size_t featureIdx);
  */
  
  template <typename T> void transpose(vector<vector<T> >& mat);
  
  //Finds the best split for target with respect to selected feature splitter, which needs to be numerical.
  void numericalFeatureSplit(vector<num_t>& tv,
			     const bool isTargetNumerical,
			     vector<num_t>& fv,
			     const size_t min_split,
			     vector<size_t>& sampleIcs_left,
			     vector<size_t>& sampleIcs_right,
			     num_t& splitValue);
    
  void categoricalFeatureSplit(vector<num_t>& tv,
			       const bool isTargetNumerical,
			       vector<num_t>& fv,
			       vector<size_t>& sampleIcs_left,
			       vector<size_t>& sampleIcs_right,
			       set<num_t>& categories_left);
  
  
  //Splits a set of samples to "left" and "right", given a splitidx
  /*
    void split_samples(vector<size_t>& sampleics,
    size_t splitidx,
    vector<size_t>& sampleics_left,
    vector<size_t>& sampleics_right);  
  */    


  size_t targetidx_;

  struct Feature {
    vector<num_t> data;
    vector<num_t> contrast;
    bool isNumerical;
    string name;
  };

  vector<Feature> features_;
  //vector<bool> isfeaturenum_;
  size_t nrealsamples_; //WILL BE DEPRECATED

  size_t nsamples_;
  size_t nfeatures_;

  //vector<string> featureheaders_;
  //vector<string> sampleheaders_;

  MTRand_int32 irand_;
};

#endif
