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

Node* Node::percolateData(num_t value) {

  if ( !hasChildren_ ) { return( this ); }

  if ( isSplitterNumerical_ ) {
    if ( value <= threshold_ ) { return( leftChild_ ); } else { return( rightChild_ ); }
  } else {
    if ( classSet_.find(value) != classSet_.end() ) { return( leftChild_ ); } else { return( rightChild_ ); }
  }
}

// !! Inefficient: just use NULL as a sentinel here. Proxying through another
// !! variable (hasChildren_) is just asking for bugs. If you really do want to
// !! use hasChildren as a resource lock, at the very least drop an assert here.
Node* Node::leftChild() {
  //assert( hasChildren_ ) {
  return( leftChild_ );
  //} else {
  //return(NULL);
  //}
}

// !! Inefficient: just use NULL as a sentinel here. Proxying through another
// !! variable (hasChildren_) is just asking for bugs. If you really do want to
// !! use hasChildren as a resource lock, at the very least drop an assert here.
Node* Node::rightChild() {
  //assert( hasChildren_ );
  return( rightChild_ );
  //} else {
  //  return(NULL);
  //}
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
						    GI,
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
						  GI,
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
    //cout << "Stopping tree generation after creation of " << nNodes << " nodes" << endl;
    return;
  }

  if ( false ) {
    cout << "Out of ";
    for ( size_t i = 0; i < featureSampleIcs.size(); ++i ) {
      cout << " " << treeData->getFeatureName(featureSampleIcs[i]);
    }
    cout << endl;
    cout << " ---- Feature " << treeData->getFeatureName(splitFeatureIdx) << " splits the data: [";
    for ( size_t i = 0; i < sampleIcs_left.size(); ++i ) {
      cout << " " << treeData->getRawFeatureData(targetIdx,sampleIcs_left[i]);
    }
    cout << " ] <==> [ ";
    for ( size_t i = 0; i < sampleIcs_right.size(); ++i ) {
      cout << " " << treeData->getRawFeatureData(targetIdx,sampleIcs_right[i]);
    }
    cout << " ] FITNESS: " << splitFitness << endl;
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
			       const GrowInstructions& GI,
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
  
  num_t newSplitFitness;
  vector<size_t> newSampleIcs_left;
  vector<size_t> newSampleIcs_right;
  num_t newSplitValue;
  set<num_t> newSplitValues_left;
  
  for ( size_t i = 0; i < nFeaturesForSplit; ++i ) {
    
    //vector<num_t> newSplitFeatureData;
    size_t newSplitFeatureIdx = featureSampleIcs[i];
    bool isFeatureNumerical = treeData->isFeatureNumerical(newSplitFeatureIdx);

    //Neither the real nor the contrast feature can appear in the tree as splitter
    if ( newSplitFeatureIdx == targetIdx ) {
      continue;
    }

    vector<num_t> featureData;
    treeData->getFeatureData(newSplitFeatureIdx,sampleIcs,featureData);

    if ( isFeatureNumerical ) {
      Node::numericalFeatureSplit(targetData,
				  isTargetNumerical,
				  featureData,
				  GI,
				  newSampleIcs_left,
				  newSampleIcs_right,
				  newSplitValue,
				  newSplitFitness);
    } else {
      Node::categoricalFeatureSplit(targetData,
				    isTargetNumerical,
				    featureData,
				    GI,
				    newSampleIcs_left,
				    newSampleIcs_right,
				    newSplitValues_left,
				    newSplitFitness);
    }

    if( newSplitFitness > splitFitness &&
	newSampleIcs_left.size() >= GI.minNodeSizeToStop &&
	newSampleIcs_right.size() >= GI.minNodeSizeToStop ) {
      
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
				 const GrowInstructions& GI, 
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
				GI,
				sampleIcs_left,
				sampleIcs_right,
				splitValue,
				splitFitness);
  } else {
    Node::categoricalFeatureSplit(targetData,
				  isTargetNumerical,
				  targetData,
				  GI,
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
					       GI.minNodeSizeToStop,
					       sampleIcs_left,
					       sampleIcs_right);
   
    if( newSplitFitness > splitFitness && 
       sampleIcs_left.size() >= GI.minNodeSizeToStop && 
       sampleIcs_right.size() >= GI.minNodeSizeToStop ) {
      
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
				GI,
				sampleIcs_left,
				sampleIcs_right,
				splitValue,
				splitFitness);
  
  } else {
  
    Node::categoricalFeatureSplit(targetData,
				  isTargetNumerical,
				  featureData,
				  GI,
				  sampleIcs_left,
				  sampleIcs_right,
				  splitValues_left,
				  splitFitness);

  }

  if ( sampleIcs_left.size() < GI.minNodeSizeToStop || sampleIcs_right.size() < GI.minNodeSizeToStop ) {
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
inline void Node::cleanPairVectorFromNANs(//const vector<num_t>& v1_copy, 
                                          //const vector<num_t>& v2_copy, 
                                          vector<num_t>& v1, 
                                          vector<num_t>& v2, 
					  vector<size_t>& mapIcs) {
  size_t n = v1.size();
  //v1.resize(n);
  //v2.resize(n);
  mapIcs.resize(n);
  size_t nReal = 0;
  for(size_t i = 0; i < n; ++i) {
    if(!datadefs::isNAN(v1[i]) && !datadefs::isNAN(v2[i])) {
      mapIcs[nReal] = i;
      v1[nReal] = v1[i];
      v2[nReal] = v2[i];
      ++nReal;
    }
  }
  v1.resize(nReal);
  v2.resize(nReal);
  mapIcs.resize(nReal);
}

// !! Correctness, Inadequate Abstraction: kill this method with fire. Refactor, REFACTOR, _*REFACTOR*_.
void Node::numericalFeatureSplit(vector<num_t> tv,
                                 const bool isTargetNumerical,
                                 vector<num_t> fv,
                                 const GrowInstructions& GI,
                                 vector<size_t>& sampleIcs_left,
                                 vector<size_t>& sampleIcs_right,
                                 num_t& splitValue,
                                 num_t& splitFitness) {

  assert(tv.size() == fv.size());

  vector<size_t> mapIcs;
 
  Node::cleanPairVectorFromNANs(tv,fv,mapIcs);

  size_t n_tot = tv.size();
  size_t n_right = n_tot;
  size_t n_left = 0;

  if(n_tot < 2 * GI.minNodeSizeToStop) {
    splitFitness = datadefs::NUM_NAN;
    return;
  }

  sampleIcs_left.clear();
  sampleIcs_right.clear();

  //Make reference indices that define the sorting wrt. feature
  bool isIncreasingOrder = true;//vector<size_t> refIcs;

  //Sort feature vector and collect reference indices
  datadefs::sortDataAndMakeRef(isIncreasingOrder,fv,sampleIcs_right);

  //Use the reference indices to sort sample indices
  datadefs::sortFromRef<num_t>(tv,sampleIcs_right);

  int bestSplitIdx = -1;
  
  //If the target is numerical, we use the iterative squared error formula to update impurity scores while we traverse "right"
  if(isTargetNumerical) {
    num_t mu_right = 0.0;
    num_t se_right = 0.0;
    num_t mu_left = 0.0;
    num_t se_left = 0.0;
    num_t se_best = 0.0;
    num_t se_tot = 0.0;

    datadefs::sqerr(tv,mu_right,se_right,n_right);
    assert(n_tot == n_right);
    se_best = se_right;
    se_tot = se_right;

    size_t idx = 0;
    while(n_left < n_tot - GI.minNodeSizeToStop) {
      datadefs::forward_backward_sqerr(tv[idx],n_left,mu_left,se_left,n_right,mu_right,se_right);
      if( se_left + se_right < se_best && n_left >= GI.minNodeSizeToStop) {
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

    datadefs::sqfreq(tv,freq_right,sf_right,n_right);
    num_t sf_tot = sf_right;
    num_t nsf_best = 1.0 * sf_right / n_right;
    assert(n_tot == n_right);

    size_t idx = 0;
    while(n_left < n_tot - GI.minNodeSizeToStop) {
      datadefs::forward_backward_sqfreq(tv[idx],n_left,freq_left,sf_left,n_right,freq_right,sf_right);
      if(1.0 * n_right * sf_left + 1.0 * n_left * sf_right > n_left * n_right * nsf_best && n_left >= GI.minNodeSizeToStop) {
        bestSplitIdx = idx;
        nsf_best = 1.0 * sf_left / n_left + 1.0 * sf_right / n_right;
	splitFitness = ( -1.0 * n_left*n_right*sf_tot + n_tot*n_right*sf_left + n_tot*n_left*sf_right ) / ( n_left*n_right * (1.0*n_tot*n_tot - sf_tot) );
      }
      ++idx;
    }
    //splitFitness = ( -1.0 * n_left*n_right*sf_tot + n_tot*n_right*sf_left + n_tot*n_left*sf_right ) / ( n_left*n_right * (1.0*n_tot*n_tot - sf_tot) );
  }

  if(bestSplitIdx == -1) {
    splitFitness = 0.0;
    return;
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
      cout << " " << sampleIcs_left[i] << ":" << tv[sampleIcs_left[i]];
    }
    cout << " ] <==> [";
    for(size_t i = 0; i < sampleIcs_right.size(); ++i) {
      cout << " " << sampleIcs_right[i] << ":" << tv[sampleIcs_right[i]];
    }
    cout << " ]" << endl;
  }
  
}

// !! Inadequate Abstraction: Refactor me.
void Node::categoricalFeatureSplit(vector<num_t> tv,
                                   const bool isTargetNumerical,
                                   vector<num_t> fv,
				   const GrowInstructions& GI,
                                   vector<size_t>& sampleIcs_left,
                                   vector<size_t>& sampleIcs_right,
                                   set<num_t>& categories_left,
                                   num_t& splitFitness) {

  //cout << "Node::categoricalFeatureSplit..." << endl;

  assert(tv.size() == fv.size());

  vector<size_t> mapIcs;
  Node::cleanPairVectorFromNANs(tv,fv,mapIcs);

  // Map all feature categories to the corresponding samples and represent it as map. The map is used to assign samples to left and right branches
  map<num_t,vector<size_t> > fmap_right;
  map<num_t,vector<size_t> > fmap_left;
  map<num_t,vector<size_t> > fmap_right_best;
  map<num_t,vector<size_t> > fmap_left_best;

  size_t n_tot = 0;
  datadefs::map_data(fv,fmap_right,n_tot);
  size_t n_right = n_tot;
  size_t n_left = 0;

  if(n_tot < GI.minNodeSizeToStop) {
    splitFitness = datadefs::NUM_NAN;
    return;
  }
  
  map<size_t,num_t> int2num;
  size_t iter = 0;
  for ( map<num_t,vector<size_t> >::const_iterator it( fmap_right.begin() ); it != fmap_right.end(); ++it ) {
    int2num.insert( pair<size_t,num_t>( iter, it->first ) );
    //cout << iter << "->" << it->first << endl;
    ++iter;
  }

  // A variable to determine the index for the last sample in the partition sequence: lastSample = GI.partitionSequence->at(psMax)
  size_t psMax = 0;
  if ( fmap_right.size() > 2 ) {
    psMax = ( 1 << (fmap_right.size() - 2) ); // 2^( fmap_right.size() - 2 )
  }

  //cout << "Splitter has " << fmap_right.size() << " categories => psMax is " << psMax << endl;

  bool foundSplit = false;

  if ( isTargetNumerical ) {

    //cout << "Target is numerical" << endl;

    num_t mu_right;
    num_t mu_left = 0.0;
    num_t se_right;
    num_t se_left = 0.0;
    datadefs::sqerr(tv,mu_right,se_right,n_right);
    assert(n_tot == n_right);
    num_t se_best = se_right;
    num_t se_tot = se_right;

    for ( size_t psIdx = 0; psIdx <= psMax; ++psIdx ) {

      //cout << "psIdx = " << psIdx << " <= " << psMax << ", PS(psIdx) = " << GI.partitionSequence->at(psIdx) << " is thrown " << flush;
            
      // If the category is added from right to left
      if ( GI.partitionSequence->isAdded(psIdx) ) {

	//cout << "from right to left: ics [";
      
	//Take samples from right and put them left
	map<num_t,vector<size_t> >::iterator it( fmap_right.find( int2num[ GI.partitionSequence->at(psIdx) ] ) );
	for(size_t i = 0; i < it->second.size(); ++i) {
	  //cout << " " << it->second[i];
	  datadefs::forward_backward_sqerr(tv[ it->second[i] ],n_left,mu_left,se_left,n_right,mu_right,se_right);
	  //cout << n_left << "\t" << n_right << "\t" << se_left << "\t" << se_right << endl;
	}
	//cout << " ]" << endl;

	fmap_left.insert( *it );
	fmap_right.erase( it->first );

      } else {
	
	//cout << "from left to right: ics [";

        //Take samples from left back to right
	map<num_t,vector<size_t> >::iterator it( fmap_left.find( int2num[ GI.partitionSequence->at(psIdx) ] ) );
        for(size_t i = 0; i < it->second.size(); ++i) {
	  //cout << " " << it->second[i];
          //cout << tv[it->second[i]] << ": ";
          datadefs::forward_backward_sqerr(tv[ it->second[i] ],n_right,mu_right,se_right,n_left,mu_left,se_left);
          //cout << n_left << "\t" << n_right << "\t" << se_left << "\t" << se_right << endl;
        }
	//cout << " ]" << endl;

	fmap_right.insert( *it );
	fmap_left.erase( it->first );

      }

      if ( se_left+se_right < se_best && n_left >= GI.minNodeSizeToStop && n_right >= GI.minNodeSizeToStop) {
	foundSplit = true;
	fmap_left_best = fmap_left;
	fmap_right_best = fmap_right;
	se_best = se_left + se_right;
      }
    }

    splitFitness = ( se_tot - se_best ) / se_tot;

  } else {

    //cout << "Target is categorical" << endl;
    
    map<num_t,size_t> freq_left,freq_right;
    size_t sf_left = 0;
    size_t sf_right;
    datadefs::sqfreq(tv,freq_right,sf_right,n_right);
    assert(n_tot == n_right);
    
    num_t sf_tot = sf_right;
    num_t nsf_best = 1.0 * sf_right / n_right;
    
    for ( size_t psIdx = 0; psIdx <= psMax; ++psIdx ) {

      //cout << "psIdx = " << psIdx << ", PS(psIdx) = " << GI.partitionSequence->at(psIdx) << " is thrown " << flush;
      
      // If the samples corresponding to the next shifted category is from right to left 
      if ( GI.partitionSequence->isAdded(psIdx) ) {
	
	//cout << "from right to left: ics [";

	// Take samples from right and put them left
	map<num_t,vector<size_t> >::iterator it( fmap_right.find( int2num[ GI.partitionSequence->at(psIdx) ] ) );
	for(size_t i = 0; i < it->second.size(); ++i) {
	  //cout << " " << tv[it->second[i]];
	  //cout << " " << it->second[i];
	  datadefs::forward_backward_sqfreq(tv[ it->second[i] ],n_left,freq_left,sf_left,n_right,freq_right,sf_right);
	  //cout << "<-" << tv[it->second[i]] << "   :" << n_left << "," << n_right << "," << sf_left << "," << sf_right << endl;
	}
	//cout << " ]" << endl;

	fmap_left.insert( *it );
        fmap_right.erase( it->first );
	
      } else {

	//cout << "from left to right: ics [";

        //Take samples from left back to right
	map<num_t,vector<size_t> >::iterator it( fmap_left.find( int2num[ GI.partitionSequence->at(psIdx) ] ) );
        for(size_t i = 0; i < it->second.size(); ++i) {
          //cout << " " << tv[it->second[i]];
	  //cout << " " << it->second[i];
          datadefs::forward_backward_sqfreq(tv[ it->second[i] ],n_right,freq_right,sf_right,n_left,freq_left,sf_left);
          //cout << "  " << tv[it->second[i]] << "-> :" << n_left << "," << n_right << "," << sf_left << "," << sf_right << endl;
        }
        //cout << " ]" << endl;

	fmap_right.insert( *it );
        fmap_left.erase( it->first );

      }
      
      if ( 1.0*n_right*sf_left + n_left*sf_right > n_left*n_right*nsf_best && n_left >= GI.minNodeSizeToStop && n_right >= GI.minNodeSizeToStop) {
	foundSplit = true;
	fmap_left_best = fmap_left;
	fmap_right_best = fmap_right;
	nsf_best = 1.0*sf_left/n_left + 1.0*sf_right/n_right;
	splitFitness = ( -1.0 * n_left*n_right*sf_tot + n_tot*n_right*sf_left + n_tot*n_left*sf_right ) / ( n_left*n_right * (1.0*n_tot*n_tot - sf_tot) );
      }
            
    }
    //cout << n_left << "," << sf_left << " <-> " << n_right << "," << "," << sf_right << endl;
    //splitFitness = ( -1.0 * n_left*n_right*sf_tot + n_tot*n_right*sf_left + n_tot*n_left*sf_right ) / ( n_left*n_right * (1.0*n_tot*n_tot - sf_tot) );
  }
  
  if(!foundSplit) {
    splitFitness = 0.0;
    return;
  }

  // Assign samples and categories on the left
  sampleIcs_left.resize(n_tot);
  categories_left.clear();
  iter = 0;
  for ( map<num_t,vector<size_t> >::const_iterator it(fmap_left_best.begin()); it != fmap_left_best.end(); ++it ) {
    for ( size_t i = 0; i < it->second.size(); ++i ) {
      sampleIcs_left[iter] = mapIcs[ it->second[i] ];
      ++iter;
    }
    categories_left.insert( it->first );
  }
  sampleIcs_left.resize(iter);

  // Assign samples on the right
  sampleIcs_right.resize(n_tot);
  iter = 0;
  for ( map<num_t,vector<size_t> >::const_iterator it(fmap_right_best.begin()); it != fmap_right_best.end(); ++it ) {
    for ( size_t i = 0; i < it->second.size(); ++i ) {
      sampleIcs_right[iter] = mapIcs[ it->second[i] ];
      ++iter;
    }
  }
  sampleIcs_right.resize(iter);
  
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
  
