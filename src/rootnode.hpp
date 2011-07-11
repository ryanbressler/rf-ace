#ifndef ROOTNODE_HPP
#define ROOTNODE_HPP

#include <cstdlib>
#include "node.hpp"
#include "treedata.hpp"

class RootNode : public Node 
{
public:

  RootNode();
  ~RootNode();

  void growTree(Treedata* treeData,
		const size_t targetIdx,
		const bool sampleWithReplacement,
		const num_t sampleSize,
		const size_t maxNodesToStop,
		const size_t minNodeSizeToStop,
		const bool isRandomSplit,
		const size_t nFeaturesInSample,
		vector<size_t>& oobIcs,
		set<size_t>& featuresInTree,
		size_t& nNodes);
  
private:


};

#endif
