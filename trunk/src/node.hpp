//node.hpp
//
//A node class for CARTs

#ifndef NODE_HPP
#define NODE_HPP

#include <cstdlib>
#include <vector>
#include <set>
#include "datadefs.hpp"
#include "treedata.hpp"

using namespace std;
using datadefs::num_t;

class Node {
public:
  //Initializes node.
  Node();
  ~Node();

  //Sets a splitter feature for the node.
  //NOTE: splitter can be assigned only once! Subsequent setter calls will raise an assertion failure.
  void setSplitter(size_t splitterIdx, set<num_t> classSet);
  void setSplitter(size_t splitterIdx, num_t threshold);

  //Gets the splitter for the node
  inline size_t getSplitter() { return(splitter_); }

  //Given a value, descends to either one of the child nodes, if existing, otherwise returns a pointer to the current node
  Node* percolateData(num_t value);

  void setLeafTrainPrediction(const num_t trainPrediction);
  num_t getLeafTrainPrediction();

  void accumulateLeafTestPredictionError(const num_t newTestData);
  void eraseLeafTestPredictionError();

  //Logic test whether the node has children or not
  inline bool hasChildren() { return(hasChildren_); }

  Node* leftChild();
  Node* rightChild();

  size_t nNodes();

  //Helper functions
  //void print();
  //void print_compact();

protected:

  void recursiveNodeSplit(Treedata* treeData,
			  const size_t targetIdx,
			  const vector<size_t>& sampleIcs,
			  const bool sampleWithReplacement,
			  const num_t sampleSizeFraction,
			  const size_t maxNodesToStop,
			  const size_t minNodeSizeToStop,
			  const bool isRandomSplit,
			  const size_t nFeaturesForSplit,
			  const bool useContrasts,
			  set<size_t>& featuresInTree,
			  size_t& nNodes);
    

private:

  void numericalFeatureSplit(vector<num_t> tv,
			     const bool isTargetNumerical,
			     vector<num_t> fv,
			     const size_t min_split,
			     vector<size_t>& sampleIcs_left,
			     vector<size_t>& sampleIcs_right,
			     num_t& splitValue);

  void categoricalFeatureSplit(vector<num_t> tv,
			       const bool isTargetNumerical,
			       vector<num_t> fv,
			       vector<size_t>& sampleIcs_left,
			       vector<size_t>& sampleIcs_right,
			       set<num_t>& categories_left);

  num_t splitFitness(vector<num_t> const& data,
                     bool const& isFeatureNumerical,
                     size_t const& minSplit,
                     vector<size_t> const& sampleIcs_left,
                     vector<size_t> const& sampleIcs_right);

  
  void setLeafTrainPrediction(const bool isTargetNumerical, const vector<num_t>& trainData);
  //void accumulateLeafTestPredictionError(const num_t newTestData);
  //void eraseLeafTestPredictionError();

  //A recursive function that will accumulate the number of descendant nodes the current nodes has
  void recursiveNDescendantNodes(size_t& n);

  bool isSplitterNumerical_;

  size_t splitter_;
  num_t threshold_;
  set<num_t> classSet_;

  bool isTrainPredictionSet_;
  num_t trainPrediction_;
  size_t nTestSamples_;
  num_t testPredictionError_;

  // WILL BECOME DEPRECATED
  //num_t prediction_; // saves the prediction of the node
    
  bool hasChildren_;
  Node* leftChild_;
  Node* rightChild_;

};

#endif
