#ifndef ROOTNODE_HPP
#define ROOTNODE_HPP

#include <cstdlib>
#include <map>
#include <vector>
#include <set>
#include <utility>
#include <fstream>
#include "node.hpp"
#include "treedata.hpp"
#include "options.hpp"
#include "distributions.hpp"
#include "datadefs.hpp"

using datadefs::num_t;
using datadefs::cat_t;

class RootNode : public Node {
public:

  // Empty tree
  RootNode();
  
  // Learn Tree from data
  RootNode(Treedata* trainData, const size_t targetIdx, const distributions::PMF* pmf, const ForestOptions* forestOptions, distributions::Random* random);

  // Load tree from file
  RootNode(ifstream& treeStream);

  ~RootNode();

  void reset(const size_t nNodes);

  void loadTree(ifstream& treeStream);

  void writeTree(ofstream& toFile);

  void growTree(Treedata* trainData, const size_t targetIdx, const distributions::PMF* pmf, const ForestOptions* forestOptions, distributions::Random* random);

  Node& childRef(const size_t childIdx);
  
  size_t nNodes() const;
  
  size_t nLeaves() const;

  const Prediction& getPrediction(Treedata* treeData, const size_t sampleIdx);

  vector<num_t> getChildLeafNumTrainData(Treedata* treeData, const size_t sampleIdx);
  vector<cat_t> getChildLeafCatTrainData(Treedata* treeData, const size_t sampleIdx);

  vector<size_t> getOobIcs();

  size_t nOobSamples();

  set<size_t> getFeaturesInTree() { return( featuresInTree_ ); }

  string getTargetName() const { return( targetName_ ); }
  bool isTargetNumerical() const { return( isTargetNumerical_ ); }

  vector<pair<size_t,size_t> > getMinDistFeatures();

  void verifyIntegrity() const;

#ifndef TEST__
private:
#endif

  size_t getTreeSizeEstimate(const size_t nSamples, const size_t nMaxLeaves, const size_t nodeSize) const;

  forest_t forestType_;
  string targetName_;
  bool isTargetNumerical_;

  // Parameters that are generated only when a tree is grown
  vector<Node> children_;

  size_t nLeaves_;

  vector<size_t> bootstrapIcs_;
  vector<size_t> oobIcs_;

  set<size_t> featuresInTree_;

  vector<size_t> minDistToRoot_;

  SplitCache splitCache_;

};

#endif
