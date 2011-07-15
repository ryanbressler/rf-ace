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

  void setPrediction(num_t value);
  num_t getPrediction();

  //void set_trainidx(size_t trainidx);
  //vector<size_t>* get_trainics();
  //void clearTrainIcs();

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
			  const size_t maxNodesToStop,
			  const size_t minNodeSizeToStop,
			  const bool isRandomSplit,
			  const size_t nFeaturesInSample,
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



  //A recursive function that will accumulate the number of descendant nodes the current nodes has
  void recursiveNDescendantNodes(size_t& n);

  bool isSplitterNumerical_;

  size_t splitter_;
  num_t threshold_;
  set<num_t> classSet_;

  //num_t impurity_;
  num_t prediction_; // saves the prediction of the node
    
  bool hasChildren_;
  Node* leftChild_;
  Node* rightChild_;

};

#endif
