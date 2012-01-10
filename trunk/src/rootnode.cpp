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
                   bool useContrasts,
                   size_t numClasses,
		   PartitionSequence* partitionSequence):
  Node() {
  GI_.sampleWithReplacement = sampleWithReplacement;
  GI_.sampleSizeFraction = sampleSizeFraction;
  GI_.maxNodesToStop = maxNodesToStop;
  GI_.minNodeSizeToStop = minNodeSizeToStop;
  GI_.isRandomSplit = isRandomSplit;
  GI_.nFeaturesForSplit = nFeaturesForSplit;
  GI_.useContrasts = useContrasts;
  GI_.numClasses = numClasses;
  GI_.partitionSequence = partitionSequence; 
}


RootNode::~RootNode() { 
  /* Empty dtor */
}



void RootNode::growTree(Treedata* treeData,
                        const size_t targetIdx,
                        //void (*leafPredictionFunction)(const vector<num_t>&, const size_t),
                        const LeafPredictionFunctionType leafPredictionFunctionType,
                        vector<size_t>& oobIcs,
                        set<size_t>& featuresInTree,
                        size_t& nNodes) {

  if ( this->hasChildren() ) {
    this->deleteTree();
  }

  if ( leafPredictionFunctionType == LEAF_MEAN ) {
    GI_.leafPredictionFunction = &RootNode::leafMean;
  } else if ( leafPredictionFunctionType == LEAF_MODE ) {
    GI_.leafPredictionFunction = &RootNode::leafMode;
  } else {
    GI_.leafPredictionFunction = &RootNode::leafGamma;
  }

  if ( false ) {
    cout << "Growing a tree: samplewithReplacement=" << GI_.sampleWithReplacement 
         << " sampleSizeFraction=" << GI_.sampleSizeFraction << " maxNodesToStop=" << GI_.maxNodesToStop
         << " minNodeSizeToStop=" << GI_.minNodeSizeToStop << " isRandomSplit=" << GI_.isRandomSplit << " nFeaturesForSplit=" << GI_.nFeaturesForSplit << endl; 
  }

  //Generate the vector for bootstrap indices
  vector<size_t> bootstrapIcs;
  
  //Generate bootstrap indices and oob-indices
  treeData->bootstrapFromRealSamples(GI_.sampleWithReplacement, GI_.sampleSizeFraction, targetIdx, bootstrapIcs, oobIcs);

  //This is to check that the bootstrap sample doesn't contain any missing values (it shouldn't!)
  if(false) {
    vector<num_t> targetData = treeData->getFeatureData(targetIdx,bootstrapIcs);
    for(size_t i = 0; i < targetData.size(); ++i) {
      assert(!datadefs::isNAN(targetData[i]));
    }

    targetData = treeData->getFeatureData(targetIdx,oobIcs);
    for(size_t i = 0; i < targetData.size(); ++i) {
      assert(!datadefs::isNAN(targetData[i]));
    }
    cout << "bootstrap samples look ok, no missing values detected" << endl;
  }



  if(false) {
    cout << "tree bootstrap indices [";
    for(size_t i = 0; i < bootstrapIcs.size(); ++i) {
      cout << " " << bootstrapIcs[i];
    }
    cout << " ]  oob [";
    for(size_t i = 0; i < oobIcs.size(); ++i) {
      cout << " " << oobIcs[i];
    }
    cout << " ]" << endl << endl;
  }

  nNodes = 1;
  featuresInTree.clear();

  //Start the recursive node splitting from the root node. This will generate the tree.
  this->recursiveNodeSplit(treeData,targetIdx,bootstrapIcs,GI_,featuresInTree,nNodes);
  
}
