//node.hpp
//
//A node class for CARTs

#ifndef NODE_HPP
#define NODE_HPP

#include <cstdlib>
#include <vector>
//#include <queue>
#include <map>
#include <set>
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
  //Initializes node.
  Node();
  ~Node();

  //Gets the splitter for the node
  //inline size_t splitterIdx() { return( splitter_.idx ); }
  inline string splitterName() { return( splitter_.name ); }

  //Sets a splitter feature for the node.
  //NOTE: splitter can be assigned only once! Subsequent setter calls will raise an assertion failure.
  void setSplitter(const string& splitterName,
                   const num_t splitLeftLeqValue);

  void setSplitter(const string& splitterName,
                   const set<string>& leftSplitValues,
                   const set<string>& rightSplitValues);

  void setSplitter(const string& splitterName,
		   const uint32_t hashIdx);

  //Given a value, descends to either one of the child nodes, if existing, otherwise returns a pointer to the current node
  Node* percolate(Treedata* testData, const size_t sampleIdx, const size_t scrambleFeatureIdx = datadefs::MAX_IDX);

  void setTrainPrediction(const num_t trainPrediction, const string& rawTrainPrediction );
  
  num_t getTrainPrediction();
  string getRawTrainPrediction();

  //Logic test whether the node has children or not
  inline bool hasChildren() { return( this->leftChild() || this->rightChild() ); }

  Node* leftChild();
  Node* rightChild();

  void deleteTree();

  void print(ofstream& toFile);
  void print(string& traversal, ofstream& toFile);

  enum PredictionFunctionType { MEAN, MODE, GAMMA };

#ifndef TEST__
protected:
#endif

  void recursiveNodeSplit(Treedata* treeData,
                          const size_t targetIdx,
			  const ForestOptions* forestOptions,
			  distributions::Random* random,
			  const PredictionFunctionType& predictionFunctionType,
			  const distributions::PMF* pmf,
			  const vector<size_t>& sampleIcs,
			  const size_t treeDepth,
			  set<size_t>& featuresInTree,
			  vector<size_t>& minDistToRoot,
                          size_t* nLeaves);

  bool regularSplitterSeek(Treedata* treeData,
			   const size_t targetIdx,
			   const ForestOptions* forestOptions,
			   distributions::Random* random,
			   const vector<size_t>& sampleIcs,
			   const vector<size_t>& featureSampleIcs,
			   size_t& splitFeatureIdx,
			   vector<size_t>& sampleIcs_left,
			   vector<size_t>& sampleIcs_right,
			   num_t& splitFitness);

#ifndef TEST__
private:
#endif

  
  struct Splitter {

    string name;
    //bool isNumerical;
    Feature::Type type;
    uint32_t hashValue;
    num_t leftLeqValue;
    set<string> leftValues;
    set<string> rightValues;
    
  } splitter_;

  num_t  trainPrediction_;
  string rawTrainPrediction_;

  Node* leftChild_;
  Node* rightChild_;

};

#endif
