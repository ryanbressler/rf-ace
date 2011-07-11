#include <vector>
#include <string>
#include "rootnode.hpp"
#include "datadefs.hpp"


RootNode::RootNode(): Node()
{
}
  
RootNode::~RootNode()
{ 
  //Node::~Node();
}


void RootNode::growTree(Treedata* treeData, 
			const size_t targetIdx, 
			const bool sampleWithReplacement, 
			const num_t sampleSize, 
			const size_t nMaxNodes, 
			const size_t minNodeSize, 
			const bool isRandomSplit, 
			const size_t nFeaturesInSample, 
			vector<size_t>& oobIcs,
			set<size_t>& featuresInTree,
			size_t& nNodes)
{
  //Generate the vector for bootstrap indices
  vector<size_t> bootstrapIcs;
  //bool withReplacement = true;
  //num_t sampleSize = 1.0;

  //Generate bootstrap indices and oob-indices
  treeData->bootstrapFromRealSamples(sampleWithReplacement, sampleSize, targetIdx, bootstrapIcs, oobIcs);

  //This is to check that the bootstrap sample doesn't contain any missing values (it shouldn't!)
  if(false)
    {
      vector<num_t> targetData;
      treeData->getFeatureData(targetIdx,bootstrapIcs,targetData);
      for(size_t i = 0; i < targetData.size(); ++i)
        {
          assert(!datadefs::isNAN(targetData[i]));
        }

      treeData->getFeatureData(targetIdx,oobIcs,targetData);
      for(size_t i = 0; i < targetData.size(); ++i)
        {
          assert(!datadefs::isNAN(targetData[i]));
        }
      cout << "the generated bootstrap sample for tree looks ok" << endl;
    }



  if(false)
    {
      cout << "tree bootstrap indices [";
      for(size_t i = 0; i < bootstrapIcs.size(); ++i)
        {
          cout << " " << bootstrapIcs[i];
        }
      cout << " ]  oob [";
      for(size_t i = 0; i < oobIcs.size(); ++i)
        {
          cout << " " << oobIcs[i];
        }
      cout << " ]" << endl << endl;
    }

  nNodes = 1;
  featuresInTree.clear();

  //Start the recursive node splitting from the root node. This will generate the tree.
  Node::recursiveNodeSplit(treeData,targetIdx,bootstrapIcs,nMaxNodes,minNodeSize,isRandomSplit,nFeaturesInSample,nNodes,featuresInTree);

}

