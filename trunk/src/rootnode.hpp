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
           size_t numClasses,
	   PartitionSequence* partitionSequence);

  ~RootNode();

  void growTree(Treedata* treeData,
                const size_t targetIdx,
                const LeafPredictionFunctionType leafPredictionFunctionType,
                vector<size_t>& oobIcs,
                set<size_t>& featuresInTree,
                size_t& nNodes);
  
  size_t nNodes();

private:

  size_t nNodes_;

  GrowInstructions GI_;

  

};

#endif
