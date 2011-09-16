#include<iostream>
#include<cassert>
#include "node.hpp"

Node::Node():
  isSplitterNumerical_(true),
  isTrainPredictionSet_(false),
  trainPrediction_(0.0),
  nTestSamples_(0),
  testPredictionError_(0.0),
  hasChildren_(false),
  leftChild_(NULL),
  rightChild_(NULL) {
}

// !! Documentation: establishes a recursive relationship that deletes all
// !! children, unless their pointers are orphaned.

// !! Correctness: Consider replacing hasChildren_ with a NULL check for
// !! leftChild_ and rightChild_, or at the very least, including an assert for
// !! such here.
Node::~Node() {
  if(hasChildren_) {
    delete leftChild_;
    delete rightChild_;
  }
}

// !! Documentation: consider combining with the documentation in the header
// !! file, fleshing it out a bit. Ideally, implementation notes should fall
// !! here; notes on the abstraction should fall in the header file.
void Node::setSplitter(size_t splitter, set<num_t> classSet) {
  assert(!hasChildren_);
  
  isSplitterNumerical_ = false;

  splitter_ = splitter;
  classSet_ = classSet;

  leftChild_ = new Node;
  rightChild_ = new Node;
  hasChildren_ = true;
}

// !! Documentation: consider combining with the documentation in the header
// !! file, fleshing it out a bit. Ideally, implementation notes should fall
// !! here; notes on the abstraction should fall in the header file.
void Node::setSplitter(size_t splitter, num_t threshold) {
  assert(!hasChildren_);
  isSplitterNumerical_ = true;

  splitter_ = splitter;
  threshold_ = threshold;

  leftChild_ = new Node;
  rightChild_ = new Node;
  hasChildren_ = true;
}


// !! Cruft: remove this after verifying the inline representation in our
// !! header file is sufficient.
/*
  int Node::getSplitter() {
  assert(hasChildren_);
  return(splitter_);
  }
*/

Node* Node::percolateData(num_t value) {

  if(!hasChildren_) { return(this); }

  if(isSplitterNumerical_) {
    if(value <= threshold_) { return(leftChild_); } else { return(rightChild_); }
  } else {
    if(classSet_.find(value) != classSet_.end()) { return(leftChild_); } else { return(rightChild_); }
  }
}

// !! Inefficient: just use NULL as a sentinel here. Proxying through another
// !! variable (hasChildren_) is just asking for bugs. If you really do want to
// !! use hasChildren as a resource lock, at the very least drop an assert here.
Node* Node::leftChild() {
  if(hasChildren_) {
    return(leftChild_);
  } else {
    return(NULL);
  }
}

// !! Inefficient: just use NULL as a sentinel here. Proxying through another
// !! variable (hasChildren_) is just asking for bugs. If you really do want to
// !! use hasChildren as a resource lock, at the very least drop an assert here.
Node* Node::rightChild() {
  if(hasChildren_) {
    return(rightChild_);
  } else {
    return(NULL);
  }
}

// !! Documentation: walks the number of nodes downwards, in linear time and
// !! log2 space. Not much to see here.
size_t Node::nNodes() {
  size_t n = 1;
  this->recursiveNDescendantNodes(n);
  return(n);
}

// !! Documentation: recursive function for the above.

// !! Inefficient: just rewrite these as a single function? 
void Node::recursiveNDescendantNodes(size_t& n) {
  if(!hasChildren_) {
    return;
  } else {
    n += 2;
    leftChild_->recursiveNDescendantNodes(n);
    rightChild_->recursiveNDescendantNodes(n);
  }
}

// !! Documentation: just your usual accessor, returning a copy of
// !! trainPrediction_.
num_t Node::getLeafTrainPrediction() {
  assert(!hasChildren_ && isTrainPredictionSet_);
  return(trainPrediction_);
}

void Node::recursiveNodeSplit(Treedata* treeData,
                              const size_t targetIdx,
                              const vector<size_t>& sampleIcs,
                              const GrowInstructions& GI,
                              set<size_t>& featuresInTree,
                              size_t& nNodes) {

  //const bool isTargetNumerical = treeData->isFeatureNumerical(targetIdx);
  size_t nSamples = sampleIcs.size();

  if(nSamples < 2 * GI.minNodeSizeToStop || nNodes >= GI.maxNodesToStop) {
    //cout << "Too few samples to start with, quitting" << endl;
    vector<num_t> leafTrainData;
    treeData->getFeatureData(targetIdx,sampleIcs,leafTrainData);
    //Node::setLeafTrainPrediction(leafTrainData,GI);

    // !! Potential Crash: This is unsafe. Add asserts or runtime checks.
    // !! Correctness: Violates the Principle of Least Knowledge. Refactor.
    (this->*GI.leafPredictionFunction)(leafTrainData,GI.numClasses);
    return;
  }
  
  vector<size_t> featureSampleIcs(GI.nFeaturesForSplit);

  if(GI.isRandomSplit) {
    if(GI.useContrasts) {
      for(size_t i = 0; i < GI.nFeaturesForSplit; ++i) {
        treeData->getRandomIndex(2*treeData->nFeatures(),featureSampleIcs[i]);
      }
    } else {
      for(size_t i = 0; i < GI.nFeaturesForSplit; ++i) {
        treeData->getRandomIndex(treeData->nFeatures(),featureSampleIcs[i]);
      }
    }
  } else {
    for(size_t i = 0; i < GI.nFeaturesForSplit; ++i) {
      featureSampleIcs[i] = i;
    }
  }

  vector<num_t> targetData;
  vector<num_t> featureData;
  
  //const bool isTargetNumerical = rootNode->isTargetNumerical();
  //treeData->getFeatureData(targetIdx,sampleIcs,targetData);

  vector<size_t> sampleIcs_left,sampleIcs_right;
  num_t splitValue;
  set<num_t> splitValues_left;

  num_t splitFitness;
  size_t splitFeatureIdx;

  bool isSplitSuccessful;
  //bool isSplitFeatureNumerical;
  

  if(GI.isOptimizedNodeSplit) {
    
    isSplitSuccessful = Node::optimizedSplitterSeek(treeData,
						    targetIdx,
						    sampleIcs,
						    featureSampleIcs,
						    GI.minNodeSizeToStop,
						    splitFeatureIdx,
						    sampleIcs_left,
						    sampleIcs_right,
						    splitValue,
						    splitValues_left,
						    splitFitness);

  } else {

    isSplitSuccessful = Node::regularSplitterSeek(treeData,
						  targetIdx,
						  sampleIcs,
						  featureSampleIcs,
						  GI.minNodeSizeToStop,
						  splitFeatureIdx,
						  sampleIcs_left,
						  sampleIcs_right,
						  splitValue,
						  splitValues_left,
						  splitFitness);
  }

  if(!isSplitSuccessful) {

    vector<num_t> leafTrainData;
    treeData->getFeatureData(targetIdx,sampleIcs,leafTrainData);

    // !! Potential Crash: This is unsafe. Add asserts or runtime checks.
    // !! Correctness: Violates the Principle of Least Knowledge. Refactor.
    (this->*GI.leafPredictionFunction)(leafTrainData,GI.numClasses);
    return;
  }

     
  if ( treeData->isFeatureNumerical(splitFeatureIdx) ) {
    //cout << "num splitter" << endl;
    Node::setSplitter(splitFeatureIdx,splitValue);
  } else {
    //cout << "cat splitter" << endl;
    Node::setSplitter(splitFeatureIdx,splitValues_left);
  }

  featuresInTree.insert(splitFeatureIdx);
  nNodes += 2;

  leftChild_->recursiveNodeSplit(treeData,targetIdx,sampleIcs_left,GI,featuresInTree,nNodes);
  rightChild_->recursiveNodeSplit(treeData,targetIdx,sampleIcs_right,GI,featuresInTree,nNodes);
  
}

bool Node::regularSplitterSeek(Treedata* treeData,
			       const size_t targetIdx,
			       const vector<size_t>& sampleIcs,
			       const vector<size_t>& featureSampleIcs,
			       const size_t minNodeSizeToStop,
			       size_t& splitFeatureIdx,
			       vector<size_t>& sampleIcs_left,
			       vector<size_t>& sampleIcs_right,
			       num_t& splitValue,
			       set<num_t>& splitValues_left,
			       num_t& splitFitness) {

  bool isTargetNumerical = treeData->isFeatureNumerical(targetIdx);
  vector<num_t> targetData;
  treeData->getFeatureData(targetIdx,sampleIcs,targetData);

  size_t nFeaturesForSplit = featureSampleIcs.size();
  splitFeatureIdx = nFeaturesForSplit;
  splitFitness = 0.0;
  
  for ( size_t i = 0; i < nFeaturesForSplit; ++i ) {
    
    //vector<num_t> newSplitFeatureData;
    size_t newSplitFeatureIdx = featureSampleIcs[i];
    //bool isFeatureNumerical = treeData->isFeatureNumerical(newSplitFeatureIdx);

    num_t newSplitFitness;
    vector<size_t> newSampleIcs_left, newSampleIcs_right;
    num_t newSplitValue;
    set<num_t> newSplitValues_left;

    vector<num_t> featureData;
    treeData->getFeatureData(newSplitFeatureIdx,sampleIcs,featureData);

    if ( isTargetNumerical ) {
      Node::numericalFeatureSplit(targetData,
				  isTargetNumerical,
				  featureData,
				  minNodeSizeToStop,
				  newSampleIcs_left,
				  newSampleIcs_right,
				  newSplitValue,
				  newSplitFitness);
    } else {
      Node::categoricalFeatureSplit(targetData,
				    isTargetNumerical,
				    featureData,
				    newSampleIcs_left,
				    newSampleIcs_right,
				    newSplitValues_left,
				    newSplitFitness);
    }

    if( newSplitFitness > splitFitness &&
	newSampleIcs_left.size() >= minNodeSizeToStop &&
	newSampleIcs_right.size() >= minNodeSizeToStop ) {
      
      splitFitness = newSplitFitness;
      splitFeatureIdx = newSplitFeatureIdx;
      splitValue = newSplitValue;
      splitValues_left = newSplitValues_left;
      sampleIcs_left = newSampleIcs_left;
      sampleIcs_right = newSampleIcs_right;

    }    

  }
  
  if ( splitFeatureIdx == nFeaturesForSplit ) {
    return(false);
  }

  for(size_t i = 0; i < sampleIcs_left.size(); ++i) {
    sampleIcs_left[i] = sampleIcs[sampleIcs_left[i]];
  }

  for(size_t i = 0; i < sampleIcs_right.size(); ++i) {
    sampleIcs_right[i] = sampleIcs[sampleIcs_right[i]];
  }

  return(true);

}

/** Optimized node split utilizes the bijection idea, such that if 
 *
 * y <- f(x)
 *
 * yields a good split, so will the inverse (assuming it exists) 
 *
 * x <- f^(-1)(y) 
 *
 * yield a good split also. This is not always true, but on average it's a good approximation. 
 * Assuming that our goal is to perform a binary split on the target, y, we can first specify 
 * an optimal split of y, i.e., split y with itself, and then project that split onto the 
 * candidates x_1 , x_2 , ... , x_i , ... , x_n. Whichever x_i splits the best is chosen to
 * make the split. The true benefit with this approximation is that the data needn't be sorted 
 * prior to testing each of the candidate splitter.
 *
 * If a splitter isn't found or there's something bad with it, the function will return false,
 * otherwise true.
 */
bool Node::optimizedSplitterSeek(Treedata* treeData, 
				 const size_t targetIdx, 
				 const vector<size_t>& sampleIcs, 
				 const vector<size_t>& featureSampleIcs, 
				 const size_t minNodeSizeToStop, 
				 size_t& splitFeatureIdx,
				 vector<size_t>& sampleIcs_left,
				 vector<size_t>& sampleIcs_right,
				 num_t& splitValue,
				 set<num_t>& splitValues_left,
				 num_t& splitFitness) {

  bool isTargetNumerical = treeData->isFeatureNumerical(targetIdx);
  vector<num_t> targetData;
  treeData->getFeatureData(targetIdx,sampleIcs,targetData);
  
  if ( isTargetNumerical ) {
    Node::numericalFeatureSplit(targetData,
				isTargetNumerical,
				targetData,
				minNodeSizeToStop,
				sampleIcs_left,
				sampleIcs_right,
				splitValue,
				splitFitness);
  } else {
    Node::categoricalFeatureSplit(targetData,
				  isTargetNumerical,
				  targetData,
				  sampleIcs_left,
				  sampleIcs_right,
				  splitValues_left,
				  splitFitness);
  }        
  
  size_t nFeaturesForSplit = featureSampleIcs.size();
  splitFeatureIdx = nFeaturesForSplit;
  splitFitness = 0.0;
  
  for ( size_t i = 0; i < nFeaturesForSplit; ++i ) {
    
    vector<num_t> newSplitFeatureData;
    size_t newSplitFeatureIdx = featureSampleIcs[i];
    bool isFeatureNumerical = treeData->isFeatureNumerical(newSplitFeatureIdx);
    
    //Neither the real nor the contrast feature can appear in the tree as splitter
    if ( newSplitFeatureIdx == targetIdx ) {
      continue;
    }
    
    treeData->getFeatureData(newSplitFeatureIdx,sampleIcs,newSplitFeatureData);
    
    num_t newSplitFitness = Node::splitFitness(newSplitFeatureData,
					       isFeatureNumerical,
					       minNodeSizeToStop,
					       sampleIcs_left,
					       sampleIcs_right);
    
    if( newSplitFitness > splitFitness && 
       sampleIcs_left.size() >= minNodeSizeToStop && 
       sampleIcs_right.size() >= minNodeSizeToStop ) {
      
      splitFitness = newSplitFitness;
      splitFeatureIdx = newSplitFeatureIdx;
    }
    
  }

  if ( splitFeatureIdx == nFeaturesForSplit ) {
    return(false);
  } 

  vector<num_t> featureData;
  treeData->getFeatureData(splitFeatureIdx,sampleIcs,featureData);

  if ( treeData->isFeatureNumerical(splitFeatureIdx) ) {
    
    Node::numericalFeatureSplit(targetData,
				isTargetNumerical,
				featureData,
				minNodeSizeToStop,
				sampleIcs_left,
				sampleIcs_right,
				splitValue,
				splitFitness);
  
  } else {
  
    Node::categoricalFeatureSplit(targetData,
				  isTargetNumerical,
				  featureData,
				  sampleIcs_left,
				  sampleIcs_right,
				  splitValues_left,
				  splitFitness);

  }

  if ( sampleIcs_left.size() < minNodeSizeToStop || sampleIcs_right.size() < minNodeSizeToStop ) {
    return(false);
  }

  for ( size_t i = 0; i < sampleIcs_left.size(); ++i ) {
    sampleIcs_left[i] = sampleIcs[sampleIcs_left[i]];
  }

  for ( size_t i = 0; i < sampleIcs_right.size(); ++i ) {
    sampleIcs_right[i] = sampleIcs[sampleIcs_right[i]];
  }

  return(true);
  
}


// !! Documentation: hey, it's a vector NAN-fixer!

// !! Correctness: Um. This should be removed once the checks for NAN go in
// !! place and a tolerable policy is created.
inline void Node::cleanPairVectorFromNANs(const vector<num_t>& v1_copy, 
                                          const vector<num_t>& v2_copy, 
                                          vector<num_t>& v1, 
                                          vector<num_t>& v2, 
                                          vector<size_t>& mapIcs) {
  v1.clear();
  v2.clear();
  mapIcs.clear();
  for(size_t i = 0; i < v1_copy.size(); ++i) {
    if(!datadefs::isNAN(v1_copy[i]) && !datadefs::isNAN(v2_copy[i])) {
      mapIcs.push_back(i);
      v1.push_back(v1_copy[i]);
      v2.push_back(v2_copy[i]);
    }
  }

}

// !! Correctness, Inadequate Abstraction: kill this method with fire. Refactor, REFACTOR, _*REFACTOR*_.
void Node::numericalFeatureSplit(const vector<num_t>& tv_copy,
                                 const bool isTargetNumerical,
                                 const vector<num_t>& fv_copy,
                                 const size_t min_split,
                                 vector<size_t>& sampleIcs_left,
                                 vector<size_t>& sampleIcs_right,
                                 num_t& splitValue,
                                 num_t& splitFitness) {

  assert(tv_copy.size() == fv_copy.size());

  vector<num_t> tv,fv;
  vector<size_t> mapIcs;
 
  Node::cleanPairVectorFromNANs(tv_copy,fv_copy,tv,fv,mapIcs);

  size_t n_tot = tv.size();
  size_t n_right = n_tot;
  size_t n_left = 0;

  if(n_tot < 2 * min_split) {
    splitFitness = datadefs::NUM_NAN;
    return;
  }

  sampleIcs_left.clear();
  sampleIcs_right.clear();

  //Check that there are enough samples to make the split in the first place
  //assert(n_tot >= 2*min_split);

  //Make reference indices that define the sorting wrt. feature
  bool isIncreasingOrder = true;//vector<size_t> refIcs;

  //Sort feature vector and collect reference indices
  datadefs::sortDataAndMakeRef(isIncreasingOrder,fv,sampleIcs_right);

  //Use the reference indices to sort sample indices
  //datadefs::sortFromRef<size_t>(sampleics,refIcs);
  datadefs::sortFromRef<num_t>(tv,sampleIcs_right);

  //Count how many real values the feature and target has
  //size_t nreal_f,nreal_t;
  //datadefs::countRealValues(fv,nreal_f);
  //datadefs::countRealValues(tv,nreal_t);
  //cout << nreal_f << " " << nreal_t << " " << n_tot << endl;
  //assert(nreal_t == n_tot && nreal_f == n_tot);

  int bestSplitIdx = -1;
  

  //If the target is numerical, we use the iterative squared error formula to update impurity scores while we traverse "right"
  if(isTargetNumerical) {
    num_t mu_right = 0.0;
    num_t se_right = 0.0;
    num_t mu_left = 0.0;
    num_t se_left = 0.0;
    num_t se_best = 0.0;
    num_t se_tot = 0.0;
    //size_t nreal_right = 0;

    datadefs::sqerr(tv,mu_right,se_right,n_right);
    assert(n_tot == n_right);
    se_best = se_right;
    se_tot = se_right;

    size_t idx = 0;
    while(n_left < n_tot - min_split) {
      datadefs::forward_backward_sqerr(tv[idx],n_left,mu_left,se_left,n_right,mu_right,se_right);
      if( se_left + se_right < se_best && n_left >= min_split) {
        bestSplitIdx = idx;
        se_best = se_left + se_right;
      }
      ++idx;
    }
    splitFitness = (se_tot - se_best) / se_tot;
  } else { //Otherwise we use the iterative gini index formula to update impurity scores while we traverse "right"
    map<num_t,size_t> freq_left;
    map<num_t,size_t> freq_right;
    size_t sf_left = 0;
    size_t sf_right = 0;
    //size_t nreal_right = 0;

    datadefs::sqfreq(tv,freq_right,sf_right,n_right);
    num_t sf_tot = sf_right;
    num_t nsf_best = 1.0 * sf_right / n_right;
    assert(n_tot == n_right);

    size_t idx = 0;
    while(n_left < n_tot - min_split) {
      datadefs::forward_backward_sqfreq(tv[idx],n_left,freq_left,sf_left,n_right,freq_right,sf_right);
      if(1.0 * n_right * sf_left + 1.0 * n_left * sf_right > n_left * n_right * nsf_best && n_left >= min_split) {
        bestSplitIdx = idx;
        nsf_best = 1.0 * sf_left / n_left + 1.0 * sf_right / n_right;
      }
      ++idx;
    }
    splitFitness = ( -1.0 * n_left*n_right*sf_tot + n_tot*n_right*sf_left + n_tot*n_left*sf_right ) / ( n_left*n_right * (1.0*n_tot*n_tot - sf_tot) );
  }

  splitValue = fv[bestSplitIdx];
  n_left = bestSplitIdx + 1;
  sampleIcs_left.resize(n_left);

  for(size_t i = 0; i < n_left; ++i) {
    sampleIcs_left[i] = mapIcs[sampleIcs_right[i]];
  }
  sampleIcs_right.erase(sampleIcs_right.begin(),sampleIcs_right.begin() + n_left);
  n_right = sampleIcs_right.size();
  for(size_t i = 0; i < n_right; ++i) {
    sampleIcs_right[i] = mapIcs[sampleIcs_right[i]];
  }

  //cout << sampleIcs_left.size() << " " << sampleIcs_right.size() << " == " << n_tot << endl;
  assert(n_left + n_right == n_tot);

  if(false) {
    cout << "Numerical feature splits target [";
    for(size_t i = 0; i < sampleIcs_left.size(); ++i) {
      cout << " " << tv_copy[sampleIcs_left[i]];
    }
    cout << " ] <==> [";
    for(size_t i = 0; i < sampleIcs_right.size(); ++i) {
      cout << " " << tv_copy[sampleIcs_right[i]];
    }
    cout << " ]" << endl;
  }
}

// !! Inadequate Abstraction: Refactor me.
void Node::categoricalFeatureSplit(const vector<num_t>& tv_copy,
                                   const bool isTargetNumerical,
                                   const vector<num_t>& fv_copy,
                                   vector<size_t>& sampleIcs_left,
                                   vector<size_t>& sampleIcs_right,
                                   set<num_t>& categories_left,
                                   num_t& splitFitness) {

  assert(tv_copy.size() == fv_copy.size());

  vector<num_t> tv,fv;
  vector<size_t> mapIcs;

  Node::cleanPairVectorFromNANs(tv_copy,fv_copy,tv,fv,mapIcs);

  size_t n_tot = tv.size();
  size_t n_right = n_tot;
  size_t n_left = 0;

  if(n_tot < 2) {
    splitFitness = datadefs::NUM_NAN;
    return;
  }

  sampleIcs_left.clear();
  sampleIcs_right.clear();
  categories_left.clear();

  //Check that sample size is positive
  //assert(n_tot > 0);

  //Count how many real values the feature and target has (this is just to make sure there are no missing values thrown in)
  //size_t nreal_f,nreal_t;
  //datadefs::countRealValues(fv,nreal_f);
  //datadefs::countRealValues(tv,nreal_t);
  //assert(nreal_t == n_tot && nreal_f == n_tot);

  //Map all feature categories to the corresponding samples and represent it as map. The map is used to assign samples to left and right branches
  map<num_t,vector<size_t> > fmap_right;
  size_t n_f;
  datadefs::map_data(fv,fmap_right,n_f);
  assert(n_tot == n_f);

  if(isTargetNumerical) {

    num_t mu_right;
    num_t mu_left = 0.0;
    num_t se_right;
    num_t se_left = 0.0;
    datadefs::sqerr(tv,mu_right,se_right,n_right);
    assert(n_tot == n_right);
    num_t se_best = se_right;
    num_t se_tot = se_right;

    while(fmap_right.size() > 0) {

      map<num_t,vector<size_t> >::iterator it_best(fmap_right.end());

      for(map<num_t,vector<size_t> >::iterator it(fmap_right.begin()); it != fmap_right.end(); ++it) {

        //Take samples from right and put them left
        for(size_t i = 0; i < it->second.size(); ++i) {
          datadefs::forward_backward_sqerr(tv[it->second[i]],n_left,mu_left,se_left,n_right,mu_right,se_right);
          //cout << n_left << "\t" << n_right << "\t" << se_left << "\t" << se_right << endl;
        }

        if(se_left+se_right < se_best) {
          it_best = it;
          se_best = se_left + se_right;
        }

        //Take samples from left back to right
        for(size_t i = 0; i < it->second.size(); ++i) {
          //cout << tv[it->second[i]] << ": ";
          datadefs::forward_backward_sqerr(tv[it->second[i]],n_right,mu_right,se_right,n_left,mu_left,se_left);
          //cout << n_left << "\t" << n_right << "\t" << se_left << "\t" << se_right << endl;
        }
      }

      //If we found a categorical split that reduces impurity...
      if(it_best == fmap_right.end()) {
        break;
      }


      //Take samples from right and put them left
      for(size_t i = 0; i < it_best->second.size(); ++i) {
        categories_left.insert(it_best->first);
        //cout << "removing index " << it_best->second[i] << " (value " << tv[it_best->second[i]] << ") from right: ";
        datadefs::forward_backward_sqerr(tv[it_best->second[i]],n_left,mu_left,se_left,n_right,mu_right,se_right);
        sampleIcs_left.push_back(it_best->second[i]);
        //cout << n_left << "\t" << n_right << "\t" << se_left << "\t" << se_right << endl;
      }
      fmap_right.erase(it_best->first);

    }
    splitFitness = (se_tot - se_best) / se_tot;
  } else {
    map<num_t,size_t> freq_left,freq_right;
    size_t sf_left = 0;
    size_t sf_right;
    datadefs::sqfreq(tv,freq_right,sf_right,n_right);
    assert(n_tot == n_right);

    num_t sf_tot = sf_right;
    num_t nsf_best = 1.0 * sf_right / n_right;

    while(fmap_right.size() > 1) {

      map<num_t,vector<size_t> >::iterator it_best(fmap_right.end());

      for(map<num_t,vector<size_t> >::iterator it(fmap_right.begin()); it != fmap_right.end(); ++it) {

        size_t n_left_c = n_left;
        size_t n_right_c = n_right;
        map<num_t,size_t> freq_left_c = freq_left;
        map<num_t,size_t> freq_right_c = freq_right;
        //Take samples from right and put them left
        //cout << "Moving " << it->second.size() << " samples corresponding to feature category " << it->first << " from right to left" << endl;
        for(size_t i = 0; i < it->second.size(); ++i) {
          //cout << " " << tv[it->second[i]];
          datadefs::forward_backward_sqfreq(tv[it->second[i]],n_left,freq_left,sf_left,n_right,freq_right,sf_right);
          //cout << n_left << "\t" << n_right << "\t" << sf_left << "\t" << sf_right << endl;
        }
        //cout << endl;

        if(1.0*n_right*sf_left + n_left*sf_right > n_left*n_right*nsf_best) {
          it_best = it;
          nsf_best = 1.0*sf_left/n_left + 1.0*sf_right/n_right;
        }

        //cout << "Moving " << it->second.size() << " samples corresponding to feature category " << it->first << " from left to right" << endl;
        //Take samples from left back to right
        for(size_t i = 0; i < it->second.size(); ++i) {
          //cout << " " << tv[it->second[i]];
          datadefs::forward_backward_sqfreq(tv[it->second[i]],n_right,freq_right,sf_right,n_left,freq_left,sf_left);
          //cout << n_left << "\t" << n_right << "\t" << sf_left << "\t" << sf_right << endl;
        }
        //cout << endl;

        assert(n_left_c == n_left);
        assert(n_right_c == n_right);
        assert(freq_left_c.size() == freq_left.size());
        assert(freq_right_c.size() == freq_right.size());

      }

      //If we found a categorical split that reduces impurity...
      if(it_best == fmap_right.end()) {
        break;
      }


      //Take samples from right and put them left
      for(size_t i = 0; i < it_best->second.size(); ++i) {
        categories_left.insert(it_best->first);
        //cout << "removing index " << it_best->second[i] << " (value " << tv[it_best->second[i]] << ") from right: ";
        datadefs::forward_backward_sqfreq(tv[it_best->second[i]],n_left,freq_left,sf_left,n_right,freq_right,sf_right);
        sampleIcs_left.push_back(it_best->second[i]);
        //cout << n_left << "\t" << n_right << "\t" << sf_left << "\t" << sf_right << endl;
      }
      fmap_right.erase(it_best->first);

    }
    splitFitness = ( -1.0 * n_left*n_right*sf_tot + n_tot*n_right*sf_left + n_tot*n_left*sf_right ) / ( n_left*n_right * (1.0*n_tot*n_tot - sf_tot) );
  }

  //Next we assign samples remaining on right
  sampleIcs_right.clear();
  for(map<num_t,vector<size_t> >::const_iterator it(fmap_right.begin()); it != fmap_right.end(); ++it) {
    for(size_t i = 0; i < it->second.size(); ++i) {
      sampleIcs_right.push_back(mapIcs[it->second[i]]);
    }
  }

  for(size_t i = 0; i < sampleIcs_left.size(); ++i) {
    sampleIcs_left[i] = mapIcs[sampleIcs_left[i]];
  }

  if(false) {
    cout << "Categorical feature splits target [";
    for(size_t i = 0; i < sampleIcs_left.size(); ++i) {
      cout << " " << tv_copy[sampleIcs_left[i]];
    }
    cout << " ] <==> [";
    for(size_t i = 0; i < sampleIcs_right.size(); ++i) {
      cout << " " << tv_copy[sampleIcs_right[i]];
    }
    cout << " ]" << endl;
  }
}

// !! Legibility: Clean out all of the print statements.
num_t Node::splitFitness(vector<num_t> const& data,
                         bool const& isFeatureNumerical,
                         size_t const& minSplit,
                         vector<size_t> const& sampleIcs_left,
                         vector<size_t> const& sampleIcs_right) {

  //assert(data.size() == sampleIcs_left.size() + sampleIcs_right.size());

  size_t n_left = 0;
  size_t n_right = 0;
  if(isFeatureNumerical) {
    num_t mu_left = 0.0;
    num_t se_left = 0.0;
    num_t mu_right = 0.0;
    num_t se_right = 0.0;

    for(size_t i = 0; i < sampleIcs_left.size(); ++i) {
      datadefs::forward_sqerr(data[sampleIcs_left[i]],n_right,mu_right,se_right);
      //cout << "forward sqerr: " << featurematrix_[featureidx][sampleics_left[i]] << " " << n_right << " " << mu_right << " " << se_right << endl;
    }

    for(size_t i = 0; i < sampleIcs_right.size(); ++i) {
      datadefs::forward_sqerr(data[sampleIcs_right[i]],n_right,mu_right,se_right);
      //cout << "forward sqerr: " << featurematrix_[featureidx][sampleics_right[i]] << " " << n_right << " " << mu_right << " " << se_right << endl;
    }

    if(n_right < 2*minSplit) {
      return(0.0);
    }

    num_t se_tot = se_right;

    for(size_t i = 0; i < sampleIcs_left.size(); ++i) {
      datadefs::forward_backward_sqerr(data[sampleIcs_left[i]],n_left,mu_left,se_left,n_right,mu_right,se_right);
      //cout << "fw bw sqerr: " << featurematrix_[featureidx][sampleics_left[i]] << " " << n_left << " " << mu_left << " " << se_left << " " << n_right
      //  << " " << mu_right << " " << se_right << endl;
    }

    if(n_left < minSplit || n_right < minSplit) {
      return(0.0);
    }

    return(( se_tot - se_left - se_right ) / se_tot);

  } else {
    map<num_t,size_t> freq_left,freq_right;
    size_t sf_left = 0;
    size_t sf_right = 0;

    for(size_t i = 0; i < sampleIcs_left.size(); ++i) {
      datadefs::forward_sqfreq(data[sampleIcs_left[i]],n_right,freq_right,sf_right);
      //cout << "forward sqfreq: " << featurematrix_[featureidx][sampleics_left[i]] << " " << n_right << " " << sf_right << endl;
    }

    for(size_t i = 0; i < sampleIcs_right.size(); ++i) {
      datadefs::forward_sqfreq(data[sampleIcs_right[i]],n_right,freq_right,sf_right);
      //cout << "forward sqfreq: " << featurematrix_[featureidx][sampleics_right[i]] << " " << n_right << " " << sf_right << endl;
    }

    if(n_right < 2*minSplit) {
      return(0.0);
    }

    size_t n_tot = n_right;
    size_t sf_tot = sf_right;

    for(size_t i = 0; i < sampleIcs_left.size(); ++i) {
      datadefs::forward_backward_sqfreq(data[sampleIcs_left[i]],n_left,freq_left,sf_left,n_right,freq_right,sf_right);
      //cout << "fw bw sqfreq: " << featurematrix_[featureidx][sampleics_left[i]] << " " << n_left << " "<< sf_left << " " << n_right << " " << sf_right << endl;
    }

    if(n_left < minSplit || n_right < minSplit) {
      return(0.0);
    }

    //cout << n_left << " " << n_right << " " << sf_tot << " " << n_tot << " " << n_right << " " << sf_left << " " << n_tot*n_left*sf_right << " " << n_left*n_right << " " << pow(n_tot,2) - sf_tot << endl;

    //num_t fitness = (-1.0*(n_left*n_right*sf_tot) + n_tot*n_right*sf_left + n_tot*n_left*sf_right) / (n_left*n_right*(pow(n_tot,2) - sf_tot));
    //cout << "Fitness " << fitness << endl;

    return( ( -1.0 * n_left*n_right*sf_tot + n_tot*n_right*sf_left + n_tot*n_left*sf_right ) / ( n_left*n_right * (1.0*n_tot*n_tot - sf_tot) ) ) ;

  }

}

void Node::leafMean(const vector<datadefs::num_t>& data, const size_t numClasses) {
  
  assert(!hasChildren_);
  assert(!isTrainPredictionSet_);
  size_t n = data.size();
  assert(n > 0);
  trainPrediction_ = 0.0;

  for(size_t i = 0; i < data.size(); ++i) {
    trainPrediction_ += data[i];
  }

  trainPrediction_ /= n;
  isTrainPredictionSet_ = true;

}

void Node::leafMode(const vector<datadefs::num_t>& data, const size_t numClasses) {

  assert(!hasChildren_);
  assert(!isTrainPredictionSet_);
  size_t n = data.size();
  assert(n > 0);
  trainPrediction_ = 0.0;

  map<num_t,size_t> freq;
  
  datadefs::count_freq(data,freq,n);
  map<num_t,size_t>::iterator it(max_element(freq.begin(),freq.end(),datadefs::freqIncreasingOrder()));
  trainPrediction_ = it->first;
  isTrainPredictionSet_ = true;

}

// !! Document
void Node::leafGamma(const vector<datadefs::num_t>& data, const size_t numClasses) {

  assert(!hasChildren_);
  assert(!isTrainPredictionSet_);
  size_t n = data.size();
  assert(n > 0);
  trainPrediction_ = 0.0;

  num_t numerator = 0.0;
  num_t denominator = 0.0;

  for (size_t i = 0; i < n; ++i) {
    num_t abs_data_i = fabs( data[i] );
    denominator += abs_data_i * (1.0 - abs_data_i);
    numerator   += data[i];
  }
  if ( fabs(denominator) <= datadefs::EPS ) {
    trainPrediction_ = datadefs::LOG_OF_MAX_NUM * numerator;
  } else {
    trainPrediction_ = (numClasses - 1)*numerator / (numClasses * denominator);
  }
  isTrainPredictionSet_ = true;
}
  
// !! Cruft: remove me.
/*
  if(isTargetNumerical) {
  num_t trainPredictionSE = 0.0;
  size_t nTrainSamples = 0;
  datadefs::sqerr(trainData,trainPrediction_,trainPredictionSE,nTrainSamples);
  
  assert(nTrainSamples > 0);
  } else {
  map<num_t,size_t> freq;
  size_t nTrainSamples = 0;
  datadefs::count_freq(trainData,freq,nTrainSamples);
  map<num_t,size_t>::iterator it(max_element(freq.begin(),freq.end(),datadefs::freqIncreasingOrder()));
  trainPrediction_ = it->first;
  
  assert(nTrainSamples > 0);
  }
  
  isTrainPredictionSet_ = true;
  }
*/

// !! Cruft: remove me.

/*
  void Node::accumulateLeafTestPredictionError(const num_t newTestData) {
  assert(!hasChildren_);
  assert(false); 
  }
  
  void Node::eraseLeafTestPredictionError() {
  assert(!hasChildren_);
  
  testPredictionError_ = 0.0;
  nTestSamples_ = 0;
  }
  
*/
