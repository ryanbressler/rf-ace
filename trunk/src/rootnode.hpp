#ifndef ROOTNODE_HPP
#define ROOTNODE_HPP

#include <cstdlib>
#include "node.hpp"
#include "treedata.hpp"

class RootNode : public Node 
{
public:

  /*
    RootNode();
    ~RootNode();
  */

  void growTree(Treedata* treeData,
		size_t targetIdx,
		bool sampleWithReplacement,
		num_t sampleSize,
		size_t maxNodesToStop,
		size_t minNodeSizeToStop,
		bool isRandomSplit,
		size_t nFeaturesInSample,
		bool useContrasts,
		vector<size_t>& oobIcs,
		set<size_t>& featuresInTree,
		size_t& nNodes);
  
private:


};

#endif
