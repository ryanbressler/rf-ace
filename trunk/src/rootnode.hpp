#ifndef ROOTNODE_HPP
#define ROOTNODE_HPP

#include <cstdlib>
#include <map>
#include <vector>
#include <set>
#include <utility>
#include "node.hpp"
#include "treedata.hpp"
#include "options.hpp"

class RootNode : public Node {
public:

  RootNode(Treedata* trainData,
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

  size_t nOobSamples();

  set<size_t> getFeaturesInTree() { return( featuresInTree_ ); }

    //    set<size_t> featuresInTree;
    //for ( minDistToRoot_t::const_iterator it( minDistToRoot_.begin() ); it != minDistToRoot_.end(); ++it ) {
    // featuresInTree.insert( it->first );
    // }
    //return( featuresInTree );
    //} // { assert(false); } //return(featuresInTree_); }

  vector<pair<size_t,size_t> > getMinDistFeatures();

#ifndef TEST__
private:
#endif

  // Parameters that are generated only when a tree is grown
  Treedata* trainData_;
  size_t nNodes_;
  vector<size_t> bootstrapIcs_;
  vector<size_t> oobIcs_;

  set<size_t> featuresInTree_;

  minDistToRoot_t minDistToRoot_;

  vector<num_t> trainPredictionCache_;

};

#endif
