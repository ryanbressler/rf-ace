//densetreedata.hpp
//
//

#ifndef TREEDATA_HPP
#define TREEDATA_HPP

#include <cstdlib>
#include <vector>

#include "datadefs.hpp"
#include "distributions.hpp"
#include "options.hpp"
#include "feature.hpp"
#include "reader.hpp"

using namespace std;
using datadefs::num_t;

class TreeData {
public:

  //TreeData();

  // Initializes the object 
  //TreeData(const vector<Feature>& features, bool useContrasts = false, const vector<string>& sampleHeaders = vector<string>(0));

  // Initializes the object and reads in a data matrix
  //TreeData(string fileName, const char dataDelimiter, const char headerDelimiter, const bool useContrasts = false);

  //~TreeData();

  // Reveals the Feature class interface to the user
  virtual const Feature* feature(const size_t featureIdx) const;
  
  // Returns the number of features
  virtual size_t nFeatures() const;
  
  // Returns feature index, given the name
  virtual size_t getFeatureIdx(const string& featureName) const;
  
  // A value denoting the "one-over-last" feature in matrix
  virtual size_t end() const;
  
  // Returns sample name, given sample index
  virtual string getSampleName(const size_t sampleIdx);
  
  // Returns the number of samples
  virtual size_t nSamples() const;
  
  virtual vector<num_t> getFeatureWeights() const;
  
  virtual void separateMissingSamples(const size_t featureIdx,
				      vector<size_t>& sampleIcs,
				      vector<size_t>& missingIcs);
  
  virtual num_t numericalFeatureSplit(const size_t targetIdx,
				      const size_t featureIdx,
				      const size_t minSamples,
				      vector<size_t>& sampleIcs_left,
				      vector<size_t>& sampleIcs_right,
				      num_t& splitValue);

  virtual num_t categoricalFeatureSplit(const size_t targetIdx,
					const size_t featureIdx,
					const vector<cat_t>& catOrder,
					const size_t minSamples,
					vector<size_t>& sampleIcs_left,
					vector<size_t>& sampleIcs_right,
					unordered_set<cat_t>& splitValues_left);
  
  virtual num_t textualFeatureSplit(const size_t targetIdx,
				    const size_t featureIdx,
				    const uint32_t hashIdx,
				    const size_t minSamples,
				    vector<size_t>& sampleIcs_left,
				    vector<size_t>& sampleIcs_right);
    
  // Generates a bootstrap sample from the real samples of featureIdx. Samples not in the bootstrap sample will be stored in oob_ics,
  // and the number of oob samples is stored in noob.
  virtual void bootstrapFromRealSamples(distributions::Random* random,
					const bool withReplacement, 
					const num_t sampleSize, 
					const size_t featureIdx, 
					vector<size_t>& ics, 
					vector<size_t>& oobIcs);

  virtual void createContrasts();
  virtual void permuteContrasts(distributions::Random* random);
  
};

#endif
