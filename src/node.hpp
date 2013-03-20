//node.hpp
//
//A node class for CARTs

#ifndef NODE_HPP
#define NODE_HPP

#include <cstdlib>
#include <vector>
#include <map>
#include <set>
#include <unordered_set>
#include <string>
#include "datadefs.hpp"
#include "treedata.hpp"
#include "options.hpp"
#include "utils.hpp"
#include "distributions.hpp"

using namespace std;
using datadefs::num_t;

class Node {
public:

  struct Prediction {

    Feature::Type type;
    num_t numTrainPrediction;
    cat_t catTrainPrediction;
    vector<num_t> numTrainData;
    vector<cat_t> catTrainData;
  };
  
  struct Splitter {
    
    num_t fitness;
    string name;
    Feature::Type type;
    uint32_t hashValue;
    num_t leftLeqValue;
    unordered_set<cat_t> leftValues;
    
    Splitter(): fitness(0.0), name(""), type(Feature::Type::UNKNOWN) {}
    
  };

  //Initializes node.
  Node();
  ~Node();
  
  //Gets the splitter for the node
  const string& splitterName() const { return( splitter_.name ); }
  
  //Sets a splitter feature for the node.
  //NOTE: splitter can be assigned only once! Subsequent setter calls will raise an assertion failure.
  void setSplitter(const num_t splitFitness,
		   const string& splitterName,
                   const num_t splitLeftLeqValue,
		   Node& leftChild,
		   Node& rightChild);
  
  void setSplitter(const num_t splitFitness,
		   const string& splitterName,
                   const unordered_set<cat_t>& leftSplitValues,
		   Node& leftChild,
		   Node& rightChild);
  
  void setSplitter(const num_t splitFitness,
		   const string& splitterName,
		   const uint32_t hashIdx,
		   Node& leftChild,
		   Node& rightChild);
  
  void setMissingChild(Node& missingChild);
  
  //Given a value, descends to either one of the child nodes, if existing, otherwise returns a pointer to the current node
  Node* percolate(Treedata* testData, const size_t sampleIdx, const size_t scrambleFeatureIdx = datadefs::MAX_IDX);
  
  void setNumTrainPrediction(const num_t& numTrainPrediction);
  void setCatTrainPrediction(const cat_t& catTrainPrediction);
  
  //Logic test whether the node has children or not
  inline bool hasChildren() const { return( this->leftChild() || this->rightChild() ); }

  Node* leftChild() const;
  Node* rightChild() const;
  Node* missingChild() const;

  vector<Node*> getSubTreeLeaves();

  void setNumTrainData(const vector<num_t>& numTrainData);
  void setCatTrainData(const vector<cat_t>& catTrainData);

  const Prediction& getPrediction();

  const Splitter& getSplitter();

  void recursiveWriteTree(string& traversal, ofstream& toFile);

  enum PredictionFunctionType { MEAN, MODE, GAMMA };

#ifndef TEST__
protected:
#endif

  struct SplitCache {

    size_t nSamples;
    vector<size_t> featureSampleIcs;

    vector<size_t> sampleIcs_left;
    vector<size_t> sampleIcs_right;
    vector<size_t> sampleIcs_missing;
    uint32_t hashIdx;
    size_t splitFeatureIdx;
    num_t splitValue;
    unordered_set<cat_t> splitValues_left;
    num_t splitFitness;

    vector<size_t> newSampleIcs_left;
    vector<size_t> newSampleIcs_right;
    vector<size_t> newSampleIcs_missing;
    uint32_t newHashIdx;
    size_t newSplitFeatureIdx;
    num_t newSplitValue;
    unordered_set<cat_t> newSplitValues_left;
    num_t newSplitFitness;

  };

  void recursiveNodeSplit(Treedata* treeData,
                          const size_t targetIdx,
			  const ForestOptions* forestOptions,
			  distributions::Random* random,
			  const PredictionFunctionType& predictionFunctionType,
			  const distributions::PMF* pmf,
			  const vector<size_t>& sampleIcs,
                          size_t* nLeaves,
			  size_t& childIdx,
			  vector<Node>& children,
			  SplitCache& splitCache);

  bool regularSplitterSeek(Treedata* treeData,
			   const size_t targetIdx,
			   const ForestOptions* forestOptions,
			   distributions::Random* random,
			   const vector<size_t>& sampleIcs,
			   size_t& childIdx,
			   vector<Node>& children,
			   SplitCache& splitCache);


  void recursiveGetSubTreeLeaves(vector<Node*>& leaves);

#ifndef TEST__
private:
#endif

  Splitter splitter_;

  Prediction prediction_;

  Node* leftChild_;
  Node* rightChild_;
  Node* missingChild_;

};

#endif
