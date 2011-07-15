#include <vector>
#include <string>
#include "rootnode.hpp"
#include "datadefs.hpp"


RootNode::RootNode(): Node()
{
}
  
RootNode::~RootNode()
{ 
}


void RootNode::growTree(Treedata* treeData, 
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
			size_t& nNodes)
{

  if(false)
    {
      cout << "Growing a tree: samplewithReplacement=" << sampleWithReplacement << " sampleSize=" << sampleSize << " maxNodesToStop=" << maxNodesToStop
	   << " minNodeSizeToStop=" << minNodeSizeToStop << " isRandomSplit=" << isRandomSplit << " nFeaturesInSample=" << nFeaturesInSample << endl; 
    }

  //Generate the vector for bootstrap indices
  vector<size_t> bootstrapIcs;
  
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
      cout << "bootstrap samples look ok, no missing values detected" << endl;
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
  Node::recursiveNodeSplit(treeData,
			   targetIdx,
			   bootstrapIcs,
			   maxNodesToStop,
			   minNodeSizeToStop,
			   isRandomSplit,
			   nFeaturesInSample,
			   useContrasts,
			   featuresInTree,
			   nNodes);

}

