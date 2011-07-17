#include <vector>
#include <string>
#include "rootnode.hpp"
#include "datadefs.hpp"



RootNode::RootNode(bool sampleWithReplacement,
		   num_t sampleSizeFraction,
		   size_t maxNodesToStop,
		   size_t minNodeSizeToStop,
		   bool isRandomSplit,
		   size_t nFeaturesForSplit,
		   bool useContrasts): 
  Node(),
  sampleWithReplacement_(sampleWithReplacement),
  sampleSizeFraction_(sampleSizeFraction),
  maxNodesToStop_(maxNodesToStop),
  minNodeSizeToStop_(minNodeSizeToStop),
  isRandomSplit_(isRandomSplit),
  nFeaturesForSplit_(nFeaturesForSplit),
  useContrasts_(useContrasts)
{
  // EMPTY
}

  
/*
  RootNode::~RootNode()
  { 
  
  ~Node();
  
  }
*/


void RootNode::growTree(Treedata* treeData,
			const size_t targetIdx,
			vector<size_t>& oobIcs,
			set<size_t>& featuresInTree,
			size_t& nNodes)
{

  if(false)
    {
      cout << "Growing a tree: samplewithReplacement=" << sampleWithReplacement_ 
	   << " sampleSizeFraction=" << sampleSizeFraction_ << " maxNodesToStop=" << maxNodesToStop_
	   << " minNodeSizeToStop=" << minNodeSizeToStop_ << " isRandomSplit=" << isRandomSplit_ << " nFeaturesForSplit=" << nFeaturesForSplit_ << endl; 
    }

  //Generate the vector for bootstrap indices
  vector<size_t> bootstrapIcs;
  
  //Generate bootstrap indices and oob-indices
  treeData->bootstrapFromRealSamples(sampleWithReplacement_, sampleSizeFraction_, targetIdx, bootstrapIcs, oobIcs);

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
			   sampleWithReplacement_,
			   sampleSizeFraction_,
			   maxNodesToStop_,
			   minNodeSizeToStop_,
			   isRandomSplit_,
			   nFeaturesForSplit_,
			   useContrasts_,
			   featuresInTree,
			   nNodes);

}

/*
  num_t RootNode::targetIdx()
  {
  return(targetIdx_);
  }
  
  bool RootNode::isTargetNumerical()
  {
  return(isTargetNumerical_);
  }
  
  bool RootNode::sampleWithReplacement()
  {
  return(sampleWithReplacement_);
  }
  
  size_t RootNode::maxNodesToStop()
  {
  return(maxNodesToStop_);
  }
  
  size_t RootNode::minNodeSizeToStop()
  {
  return(monNodeSizeToStop_);
  }
  
  bool RootNode::isRandomSplit()
  {
  return(isRandomSplit_);
  }
  
  size_t RootNode::featureSampleSize()
  {
  return(featureSampleSize_);
  }
  
  bool RootNode::useContrasts()
  {
  return(useContrasts_);
  }
*/


