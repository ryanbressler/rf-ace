#ifndef ROOTNODE_HPP
#define ROOTNODE_HPP

#include <cstdlib>
#include "node.hpp"
#include "treedata.hpp"

class RootNode : public Node {
public:

  
  RootNode(bool sampleWithReplacement,
           num_t sampleSizeFraction,
           size_t maxNodesToStop,
           size_t minNodeSizeToStop,
           bool isRandomSplit,
           size_t nFeaturesForSplit,
           bool useContrasts,
           bool isOptimizedNodeSplit,
           size_t numClasses);

  void growTree(Treedata* treeData,
                const size_t targetIdx,
                const LeafPredictionFunctionType leafPredictionFunctionType,
                vector<size_t>& oobIcs,
                set<size_t>& featuresInTree,
                size_t& nNodes);
  
private:

  GrowInstructions GI_;

};

#endif
