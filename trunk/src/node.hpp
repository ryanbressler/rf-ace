//node.hpp
//
//A node class for CARTs

#ifndef NODE_HPP
#define NODE_HPP

#include <cstdlib>
#include <vector>
#include <set>
#include <string>
#include "datadefs.hpp"
#include "treedata.hpp"
//#include "splitter.hpp"
#include "partitionsequence.hpp"

using namespace std;
using datadefs::num_t;

class Node;

//typedef void (Node::*LeafPredictionFunction)(const vector<num_t>&, const size_t);

class Node {
public:
  //Initializes node.
  Node();
  ~Node();

  //Gets the splitter for the node
  inline size_t splitterIdx() { return( splitter_.idx ); }

  //Sets a splitter feature for the node.
  //NOTE: splitter can be assigned only once! Subsequent setter calls will raise an assertion failure.
  void setSplitter(const size_t splitterIdx, const string& splitterName, num_t splitLeftLeqValue);
  void setSplitter(const size_t splitterIdx, const string& splitterName, const set<string>& leftSplitValues, const set<string>& rightSplitValues);

  //Given a value, descends to either one of the child nodes, if existing, otherwise returns a pointer to the current node
  // TODO: get rid of these as there's too much redundancy hidden within
  Node* percolateData(const num_t data);
  Node* percolateData(const string& data);

  
  //Node* percolate(Treedata* treeData, const size_t sampleIdx);
  //Node* percolateRandomly(Treedata* treeData, const size_t sampleIdx);
  //void  percolateData(Treedata* treeData, const size_t sampleIdx, Node** nodep);

  void setTrainPrediction(const num_t trainPrediction);
  num_t getTrainPrediction();

  //Logic test whether the node has children or not
  inline bool hasChildren() { return( leftChild_ || rightChild_ ); }

  Node* leftChild();
  Node* rightChild();

  void deleteTree();

  void print(ofstream& toFile);
  void print(string& traversal, ofstream& toFile);

  enum PredictionFunctionType { MEAN, MODE, GAMMA };

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
    PredictionFunctionType predictionFunctionType;
    size_t numClasses;
    //PartitionSequence* partitionSequence;
  };

  void recursiveNodeSplit(Treedata* treeData,
                          const size_t targetIdx,
                          const vector<size_t>& sampleIcs,
                          const GrowInstructions& GI,
                          set<size_t>& featuresInTree,
                          size_t* nNodes);

  bool regularSplitterSeek(Treedata* treeData,
			   const size_t targetIdx,
			   const vector<size_t>& sampleIcs,
			   const vector<size_t>& featureSampleIcs,
			   const GrowInstructions& GI,
			   size_t& splitFeatureIdx,
			   vector<size_t>& sampleIcs_left,
			   vector<size_t>& sampleIcs_right,
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

  
  inline num_t getSplitFitness(const size_t n_left, 
			       const size_t sf_left,
			       const size_t n_right,
			       const size_t sf_right,
			       const size_t n_tot,
			       const size_t sf_tot) {
    
    return( ( -1.0*n_left*n_right*sf_tot + 
	       1.0*n_tot*n_right*sf_left + 
	       1.0*n_tot*n_left*sf_right ) / 
	    ( 1.0*n_left*n_right * (1.0*n_tot*n_tot - 1.0*sf_tot) ) ); }



  //size_t splitterIdx_;
  //Splitter* splitter_;

  struct Splitter {

    size_t idx;
    string name;
    bool isNumerical;
    num_t leftLeqValue;
    set<string> leftValues;
    set<string> rightValues;

  } splitter_;

  num_t trainPrediction_;

  Node* leftChild_;
  Node* rightChild_;

};

#endif
