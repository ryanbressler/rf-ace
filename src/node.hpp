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
#include "splitter.hpp"
#include "partitionsequence.hpp"

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
  void setSplitter(size_t splitterIdx, const string& splitterName, num_t splitLeftLeqValue);
  void setSplitter(size_t splitterIdx, const string& splitterName, const set<num_t>& leftSplitValues, const set<num_t>& rightSplitValues);

  //Gets the splitter for the node
  inline size_t splitterIdx() { return(splitterIdx_); }

  //Given a value, descends to either one of the child nodes, if existing, otherwise returns a pointer to the current node
  Node* percolateData(num_t value);

  num_t getLeafTrainPrediction();

  //Logic test whether the node has children or not
  inline bool hasChildren() { return( splitter_ ); }

  Node* leftChild();
  Node* rightChild();

  size_t nNodes();

  void print(ofstream& toFile);

  enum LeafPredictionFunctionType { LEAF_MEAN, LEAF_MODE, LEAF_GAMMA };

#ifndef TEST__
protected:
#endif

  struct GrowInstructions {
    bool sampleWithReplacement;
    num_t sampleSizeFraction;
    size_t maxNodesToStop;
    size_t minNodeSizeToStop;
    bool isRandomSplit;
    size_t nFeaturesForSplit;
    bool useContrasts;
    LeafPredictionFunction leafPredictionFunction;
    size_t numClasses;
    PartitionSequence* partitionSequence;
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

  bool regularSplitterSeek(Treedata* treeData,
			   const size_t targetIdx,
			   const vector<size_t>& sampleIcs,
			   const vector<size_t>& featureSampleIcs,
			   const GrowInstructions& GI,
			   size_t& splitFeatureIdx,
			   vector<size_t>& sampleIcs_left,
			   vector<size_t>& sampleIcs_right,
			   num_t& splitValue,
			   set<num_t>& splitValues_left,
			   set<num_t>& splitValues_right,
			   num_t& splitFitness);

#ifndef TEST__
private:
#endif

  void numericalFeatureSplit(Treedata* treedata,
                             const size_t targetIdx,
                             const size_t featureIdx,
                             const GrowInstructions& GI,
                             vector<size_t>& sampleIcs_left,
                             vector<size_t>& sampleIcs_right,
			     num_t& splitValue,
                             num_t& splitFitness);

  void categoricalFeatureSplit(Treedata* treedata,
                               const size_t targetIdx,
                               const size_t featureIdx,
			       const GrowInstructions& GI,
			       vector<size_t>& sampleIcs_left,
                               vector<size_t>& sampleIcs_right,
                               set<num_t>& splitValues_left,
                               set<num_t>& splitValues_right,
                               num_t& splitFitness);

  num_t splitFitness(vector<num_t> const& data,
                     bool const& isFeatureNumerical,
                     size_t const& minSplit,
                     vector<size_t> const& sampleIcs_left,
                     vector<size_t> const& sampleIcs_right);

  //A recursive function that will accumulate the number of descendant nodes the current nodes has
  void recursiveNDescendantNodes(size_t& n);

  size_t splitterIdx_;
  Splitter* splitter_;

  num_t trainPrediction_;
  size_t nTestSamples_;
  num_t testPredictionError_;
    
  Node* leftChild_;
  Node* rightChild_;

};

#endif
