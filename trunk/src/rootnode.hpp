#ifndef ROOTNODE_HPP
#define ROOTNODE_HPP

#include <cstdlib>
#include <map>
#include <vector>
#include <set>
#include "node.hpp"
#include "treedata.hpp"
#include "options.hpp"

class RootNode : public Node {
public:

  RootNode(Treedata* treeData,
	   options::General_options* parameters,
	   const size_t threadIdx);

  ~RootNode();

  void growTree();
  
  size_t nNodes();

  num_t getTestPrediction(Treedata* treeData, const size_t sampleIdx);
  string getRawTestPrediction(Treedata* treeData, const size_t sampleIdx);

  num_t getTrainPrediction(const size_t sampleIdx);
  num_t getPermutedTrainPrediction(const size_t featureIdx,
				   const size_t sampleIdx);

  vector<num_t> getTrainPrediction();

  vector<size_t> getOobIcs();

  //num_t getOobError();
  size_t nOobSamples();

  set<size_t> getFeaturesInTree() { return(featuresInTree_); }

  //Treedata* trainData() { return(&treeData_); }

#ifndef TEST__
private:
#endif

  Node* percolateSampleIdx(Treedata* testData, const size_t sampleIdx);
  Node* percolateSampleIdx(const size_t sampleIdx);

  //map<Node*,vector<size_t> > percolateSampleIcsAtRandom(const size_t featureIdx, const vector<size_t>& sampleIcs);
  Node* percolateSampleIdxAtRandom(const size_t featureIdx, const size_t sampleIdx);

  //num_t getPredictionError(const map<Node*,vector<size_t> >& percolatedSampleIcs);

  // Required parameters
  Treedata* treeData_;
  options::General_options* parameters_;
  size_t threadIdx_;

  // Parameters that are generated only when a tree is grown
  size_t nNodes_;
  vector<size_t> bootstrapIcs_;
  vector<size_t> oobIcs_;

  set<size_t> featuresInTree_;

  vector<num_t> trainPredictionCache_;

};

#endif
