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

class Node;

typedef void (Node::*LeafPredictionFunction)(const vector<num_t>&, const size_t);

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

  //void setLeafTrainPrediction(const num_t trainPrediction);
  num_t getLeafTrainPrediction();

  //void accumulateLeafTestPredictionError(const num_t newTestData);
  //void eraseLeafTestPredictionError();

  //Logic test whether the node has children or not
  inline bool hasChildren() { return(hasChildren_); }

  Node* leftChild();
  Node* rightChild();

  size_t nNodes();

  //Leaf prediction functions
  //void leafMean(const vector<num_t>& data, const size_t numClasses);
  //void leafMode(const vector<num_t>& data, const size_t numClasses);
  //void leafGamma(const vector<num_t>& data, const size_t numClasses);

  enum LeafPredictionFunctionType { LEAF_MEAN, LEAF_MODE, LEAF_GAMMA };

  //Helper functions
  //void print();
  //void print_compact();

protected:

  struct GrowInstructions {
    bool sampleWithReplacement;
    num_t sampleSizeFraction;
    size_t maxNodesToStop;
    size_t minNodeSizeToStop;
    bool isRandomSplit;
    size_t nFeaturesForSplit;
    bool useContrasts;
    bool isOptimizedNodeSplit;
    //LeafPredictionFunctionType leafPredictionFunctionType; //void (*leafPredictionFunction)(const vector<num_t>& data, const size_t numClasses);
    LeafPredictionFunction leafPredictionFunction;//void (*leafPredictionFunction)(const vector<num_t>&, const size_t);
    size_t numClasses;
  };

  //Leaf prediction functions
  void leafMean(const vector<num_t>& data, const size_t numClasses);
  void leafMode(const vector<num_t>& data, const size_t numClasses);
  void leafGamma(const vector<num_t>& data, const size_t numClasses);

  void recursiveNodeSplit(Treedata* treeData,
                          const size_t targetIdx,
                          const vector<size_t>& sampleIcs,
                          const GrowInstructions& GI,
                          set<size_t>& featuresInTree,
                          size_t& nNodes);
    

private:

  inline void cleanPairVectorFromNANs(const vector<num_t>& v1_copy,
                                      const vector<num_t>& v2_copy,
                                      vector<num_t>& v1,
                                      vector<num_t>& v2,
                                      vector<size_t>& mapIcs);

  void numericalFeatureSplit(const vector<num_t>& tv_copy,
                             const bool isTargetNumerical,
                             const vector<num_t>& fv_copy,
                             const size_t min_split,
                             vector<size_t>& sampleIcs_left,
                             vector<size_t>& sampleIcs_right,
                             num_t& splitValue,
                             num_t& splitFitness);

  void categoricalFeatureSplit(const vector<num_t>& tv_copy,
                               const bool isTargetNumerical,
                               const vector<num_t>& fv_copy, /* MIN SPLIT MISSING (CHECK numericalFeatureSplit()) */
                               vector<size_t>& sampleIcs_left,
                               vector<size_t>& sampleIcs_right,
                               set<num_t>& categories_left,
                               num_t& splitFitness);

  num_t splitFitness(vector<num_t> const& data,
                     bool const& isFeatureNumerical,
                     size_t const& minSplit,
                     vector<size_t> const& sampleIcs_left,
                     vector<size_t> const& sampleIcs_right);

  
  //void setLeafTrainPrediction(const vector<num_t>& trainData, const GrowInstructions& GI);
  

  //A recursive function that will accumulate the number of descendant nodes the current nodes has
  void recursiveNDescendantNodes(size_t& n);

  //Leaf prediction functions
  //void leafMean(const vector<num_t>& data);
  //void leafMode(const vector<num_t>& data);
  //void leafGamma(const vector<num_t>& data, const size_t numClasses);

  bool isSplitterNumerical_;

  size_t splitter_;
  num_t threshold_;
  set<num_t> classSet_;

  bool isTrainPredictionSet_;
  num_t trainPrediction_;
  size_t nTestSamples_;
  num_t testPredictionError_;

  // WILL BECOME DEPRECATED !! Yes, it did. Can we remove it from the code-base?
  //num_t prediction_; // saves the prediction of the node
    
  bool hasChildren_;
  Node* leftChild_;
  Node* rightChild_;

};

#endif
