#include<iostream>
#include<cassert>
#include<iomanip>
#include "node.hpp"

Node::Node():
  splitterIdx_(0),
  splitter_(NULL),
  trainPrediction_(datadefs::NUM_NAN),
  leftChild_(NULL),
  rightChild_(NULL) {
  
}

/** 
 * Destructor for the node: initiates a recursive destruction for all child nodes
 */
Node::~Node() {
  
  //If the node has children, some moery cleanup needs to be performed
  if ( this->hasChildren() ) {

    // Deallocates dynamically allocated memory
    this->deleteTree();
  }

}

/**
 * Deletes child nodes, which will cascade all the way to the leaf nodes 
 */
void Node::deleteTree() {
  
  delete leftChild_;
  delete rightChild_;
  delete splitter_;
  
}

// !! Documentation: consider combining with the documentation in the header
// !! file, fleshing it out a bit. Ideally, implementation notes should fall
// !! here; notes on the abstraction should fall in the header file.
void Node::setSplitter(const size_t splitterIdx, const string& splitterName, num_t splitLeftLeqValue) {
  
  if ( this->hasChildren() ) {
    cerr << "Cannot set a splitter to a node twice!" << endl;
    exit(1);
  }
  
  splitterIdx_ = splitterIdx;
  splitter_ = new Splitter(splitterName,splitLeftLeqValue);

  leftChild_ = new Node();
  rightChild_ = new Node();

}

// !! Documentation: consider combining with the documentation in the header
// !! file, fleshing it out a bit. Ideally, implementation notes should fall
// !! here; notes on the abstraction should fall in the header file.
void Node::setSplitter(const size_t splitterIdx, const string& splitterName, const set<string>& leftSplitValues, const set<string>& rightSplitValues) {

  if ( this->hasChildren() ) {
    cerr << "Cannot set a splitter to a node twice!" << endl;
    exit(1);
  }

  splitterIdx_ = splitterIdx;
  splitter_ = new Splitter(splitterName,leftSplitValues,rightSplitValues);

  leftChild_ = new Node();
  rightChild_ = new Node();

}


Node* Node::percolateData(const num_t data) {

  // Return this if the node doesn't have children ( == is a leaf node )
  if ( !this->hasChildren() ) {
    return( this );
  }

  // Return left child if splits left
  if ( splitter_->splitsLeft( data ) ) {
    return( leftChild_ );
  }

  // Return right child if splits right
  if ( splitter_->splitsRight( data ) ) {
    return( rightChild_ );
  }

  // Return this if splits neither left nor right, which can happen if 
  // the splitter is categorical or data is NaN
  return( this );

}

Node* Node::percolateData(const string& data) {

  // Return this if the node doesn't have children ( == is a leaf node )
  if ( !this->hasChildren() ) {
    return( this );
  }

  // Return left child if splits left
  if ( splitter_->splitsLeft( data ) ) {
    return( leftChild_ );
  }

  // Return right child if splits right
  if ( splitter_->splitsRight( data ) ) {
    return( rightChild_ );
  }

  // Return this if splits neither left nor right, which can happen if
  // the splitter is categorical
  return( this );

}

Node* Node::percolateData(Treedata* treeData, const size_t sampleIdx) {

  // Start percolating from this node
  Node* nodep(this);

  // Percolate until stopped, at which point "nodep" stores the end-node
  this->percolateData(treeData,sampleIdx,&nodep);

  // Return the end-node
  return( nodep );

}

void Node::percolateData(Treedata* treeData, const size_t sampleIdx, Node** nodep) {

  
  if( !this->hasChildren() ) {
    // If the node does not have children, set "nodep" to this and stop recursion
    *nodep = this;
  } else if ( splitter_->splitsLeft(treeData,sampleIdx) ) {
    // Continue percolating left
    *nodep = leftChild_;
    leftChild_->percolateData(treeData,sampleIdx,nodep);
  } else if ( splitter_->splitsRight(treeData,sampleIdx) ) {
    // Continue percolating right
    *nodep = rightChild_;
    rightChild_->percolateData(treeData,sampleIdx,nodep);
  } else {
    // If no split can be made, set "nodep" to this and stop recursion
    *nodep = this;
  }

}

Node* Node::leftChild() {
  return( leftChild_ );
}


Node* Node::rightChild() {
  return( rightChild_ );
}

/**
 * Recursively prints a tree to a stream (file)
 */
void Node::print(string& traversal, ofstream& toFile) {

  toFile << "NODE=" << traversal << ",PRED=" << setprecision(3) << trainPrediction_;
  
  if ( this->hasChildren() ) {
    
    toFile << "," << splitter_->print() << endl;

    string traversalLeft = traversal;
    traversalLeft.append("L");
    string traversalRight = traversal;
    traversalRight.append("R");
    
    leftChild_->print(traversalLeft,toFile);
    rightChild_->print(traversalRight,toFile);

  } else {
   toFile << endl;
  }
}

void Node::print(ofstream& toFile) {
  
  string traversal("*");
  this->print(traversal,toFile);

}

void Node::setTrainPrediction(const num_t trainPrediction) {
  trainPrediction_ = trainPrediction;
}

// !! Documentation: just your usual accessor, returning a copy of
// !! trainPrediction_.
num_t Node::getTrainPrediction() {
  return( trainPrediction_ );
}

void Node::recursiveNodeSplit(Treedata* treeData,
                              const size_t targetIdx,
                              const vector<size_t>& sampleIcs,
                              const GrowInstructions& GI,
                              set<size_t>& featuresInTree,
                              size_t* nNodes) {

  size_t nSamples = sampleIcs.size();

  if(nSamples < 2 * GI.minNodeSizeToStop || *nNodes >= GI.maxNodesToStop) {
    //cout << "Too few samples to start with, quitting" << endl;
    vector<num_t> leafTrainData = treeData->getFeatureData(targetIdx,sampleIcs);

    switch ( GI.predictionFunctionType ) {
    case MEAN: 
      this->leafMean(leafTrainData);
      break;
    case MODE:
      this->leafMode(leafTrainData);
      break;
    case GAMMA:
      this->leafGamma(leafTrainData,GI.numClasses);
      break;
    }

    return;
  }
  
  vector<size_t> featureSampleIcs(GI.nFeaturesForSplit);

  if(GI.isRandomSplit) {
    if(GI.useContrasts) {
      for(size_t i = 0; i < GI.nFeaturesForSplit; ++i) {
	
        featureSampleIcs[i] = treeData->getRandomIndex(treeData->nFeatures());
	// If the sampled feature is a contrast... 
	if( treeData->getRandomIndex( 35535 ) < 355.35 ) { // %1 sampling rate
	  //cout << "Generated CONTRAST \n" << endl;
	  featureSampleIcs[i] += treeData->nFeatures();
	}
      }
    } else {
      for(size_t i = 0; i < GI.nFeaturesForSplit; ++i) {
        featureSampleIcs[i] = treeData->getRandomIndex(treeData->nFeatures());
      }
    }
  } else {
    for(size_t i = 0; i < GI.nFeaturesForSplit; ++i) {
      //cout << " " << i;
      featureSampleIcs[i] = i;
    }
  }
  //cout << endl;

  vector<size_t> sampleIcs_left,sampleIcs_right;

  size_t splitFeatureIdx;
  num_t splitFitness;
  
  bool foundSplit = this->regularSplitterSeek(treeData,
					      targetIdx,
					      sampleIcs,
					      featureSampleIcs,
					      GI,
					      splitFeatureIdx,
					      sampleIcs_left,
					      sampleIcs_right,
					      splitFitness);
  
  if ( !foundSplit ) {

    vector<num_t> leafTrainData = treeData->getFeatureData(targetIdx,sampleIcs);

    switch ( GI.predictionFunctionType ) {
    case MEAN:
      this->leafMean(leafTrainData);
      break;
    case MODE:
      this->leafMode(leafTrainData);
      break;
    case GAMMA:
      this->leafGamma(leafTrainData,GI.numClasses);
      break;
    }

    return;
  }
  
  vector<num_t> trainData = treeData->getFeatureData(targetIdx,sampleIcs);

  switch ( GI.predictionFunctionType ) {
  case MEAN:
    this->leafMean(trainData);
    break;
  case MODE:
    this->leafMode(trainData);
    break;
  case GAMMA:
    this->leafGamma(trainData,GI.numClasses);
    break;
  }


  featuresInTree.insert(splitFeatureIdx);
  *nNodes += 2;

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
			       num_t& splitFitness) {
  
  num_t splitValue = datadefs::NUM_NAN;
  set<num_t> splitValues_left,splitValues_right;
  
  size_t nFeaturesForSplit = featureSampleIcs.size();
  splitFeatureIdx = nFeaturesForSplit;
  //splitFitness = 0.0;
  
  num_t newSplitFitness;
  num_t newSplitValue;
  set<num_t> newSplitValues_left, newSplitValues_right;

  splitFitness = 0.0;
  
  for ( size_t i = 0; i < nFeaturesForSplit; ++i ) {
    
    //vector<num_t> newSplitFeatureData;
    size_t newSplitFeatureIdx = featureSampleIcs[i];
    bool isFeatureNumerical = treeData->isFeatureNumerical(newSplitFeatureIdx);

    //Neither the real nor the contrast feature can appear in the tree as splitter
    if ( newSplitFeatureIdx == targetIdx ) {
      continue;
    }

    // vector<num_t> featureData = treeData->getFeatureData(newSplitFeatureIdx,sampleIcs);

    vector<size_t> newSampleIcs_left(0);
    vector<size_t> newSampleIcs_right = sampleIcs;

    //cout << "Splitting with feature " << newSplitFeatureIdx << endl;

    if ( isFeatureNumerical ) {
      Node::numericalFeatureSplit(treeData,
				  targetIdx,
				  newSplitFeatureIdx,
				  GI,
				  newSampleIcs_left,
				  newSampleIcs_right,
				  newSplitValue,
				  newSplitFitness);
    } else {
      Node::categoricalFeatureSplit(treeData,
				    targetIdx,
				    newSplitFeatureIdx,
				    GI,
				    newSampleIcs_left,
				    newSampleIcs_right,
				    newSplitValues_left,
				    newSplitValues_right,
				    newSplitFitness);
    }

    if( newSplitFitness > splitFitness &&
	newSampleIcs_left.size() >= GI.minNodeSizeToStop &&
	newSampleIcs_right.size() >= GI.minNodeSizeToStop ) {
      
      splitFitness = newSplitFitness;
      splitFeatureIdx = newSplitFeatureIdx;
      splitValue = newSplitValue;
      splitValues_left = newSplitValues_left;
      splitValues_right = newSplitValues_right;
      sampleIcs_left = newSampleIcs_left;
      sampleIcs_right = newSampleIcs_right;
    }    

  }
  
  if ( splitFeatureIdx == nFeaturesForSplit ) {
    return(false);
  }

  if ( treeData->isFeatureNumerical(splitFeatureIdx) ) {

    this->setSplitter(splitFeatureIdx,treeData->getFeatureName(splitFeatureIdx),splitValue);

  } else {
    
    set<string> rawSplitValues_left,rawSplitValues_right;

    for ( set<num_t>::const_iterator it(splitValues_left.begin()); it != splitValues_left.end(); ++it ) {
      rawSplitValues_left.insert( treeData->getRawFeatureData(splitFeatureIdx,*it) );
    }

    for ( set<num_t>::const_iterator it(splitValues_right.begin()); it != splitValues_right.end(); ++it ) {
      rawSplitValues_right.insert( treeData->getRawFeatureData(splitFeatureIdx,*it) );
    }

    this->setSplitter(splitFeatureIdx,treeData->getFeatureName(splitFeatureIdx),rawSplitValues_left,rawSplitValues_right);

  }

  return(true);

}

// !! Correctness, Inadequate Abstraction: kill this method with fire. Refactor, REFACTOR, _*REFACTOR*_.
void Node::numericalFeatureSplit(Treedata* treedata,
				 const size_t targetIdx,
                                 const size_t featureIdx,
                                 const GrowInstructions& GI,
                                 vector<size_t>& sampleIcs_left,
                                 vector<size_t>& sampleIcs_right,
                                 num_t& splitValue,
                                 num_t& splitFitness) {

  splitFitness = datadefs::NUM_NAN;

  vector<num_t> tv,fv;

  sampleIcs_left.clear();

  treedata->getFilteredAndSortedFeatureDataPair3(targetIdx,featureIdx,sampleIcs_right,tv,fv);

  size_t n_tot = fv.size();
  size_t n_right = n_tot;
  size_t n_left = 0;

  if(n_tot < 2 * GI.minNodeSizeToStop) {
    splitFitness = datadefs::NUM_NAN;
    return;
  }

  int bestSplitIdx = -1;
  
  //If the target is numerical, we use the iterative squared error formula to update impurity scores while we traverse "right"
  if ( treedata->isFeatureNumerical(targetIdx) ) {
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
      if(1.0*n_right*sf_left + 1.0*n_left*sf_right > n_left*n_right*nsf_best && n_left >= GI.minNodeSizeToStop) {
        bestSplitIdx = idx;
        nsf_best = 1.0*sf_left / n_left + 1.0*sf_right / n_right;
	splitFitness = this->getSplitFitness(n_left,sf_left,n_right,sf_right,n_tot,sf_tot);
	//splitFitness = ( -1.0 * n_left*n_right*sf_tot + 1.0*n_tot*n_right*sf_left + 1.0*n_tot*n_left*sf_right ) / ( 1.0*n_left*n_right * (1.0*n_tot*n_tot - 1.0*sf_tot) );
      }
      ++idx;
    }
    //splitFitness = ( -1.0 * n_left*n_right*sf_tot + n_tot*n_right*sf_left + n_tot*n_left*sf_right ) / ( n_left*n_right * (1.0*n_tot*n_tot - sf_tot) );
  }

  if(bestSplitIdx == -1) {
    //cout << "N : " << n_left << " <-> " << n_right << " : fitness " << splitFitness << endl;
    return;
  }

  splitValue = fv[bestSplitIdx];
  n_left = bestSplitIdx + 1;
  sampleIcs_left.resize(n_left);

  for(size_t i = 0; i < n_left; ++i) {
    sampleIcs_left[i] = sampleIcs_right[i];
  }
  sampleIcs_right.erase(sampleIcs_right.begin(),sampleIcs_right.begin() + n_left);
  n_right = sampleIcs_right.size();

  assert(n_left + n_right == n_tot);

  //cout << "N : " << n_left << " <-> " << n_right << " : fitness " << splitFitness << endl;

}

// !! Inadequate Abstraction: Refactor me.
void Node::categoricalFeatureSplit(Treedata* treedata,
                                   const size_t targetIdx,
                                   const size_t featureIdx,
				   const GrowInstructions& GI,
                                   vector<size_t>& sampleIcs_left,
                                   vector<size_t>& sampleIcs_right,
				   set<num_t>& splitValues_left,
				   set<num_t>& splitValues_right,
                                   num_t& splitFitness) {

  splitFitness = datadefs::NUM_NAN;

  vector<num_t> tv,fv;

  //cout << " -- sampleIcs_right.size() = " << sampleIcs_right.size();

  sampleIcs_left.clear();
  treedata->getFilteredFeatureDataPair(targetIdx,featureIdx,sampleIcs_right,tv,fv);

  //cout << " => " << sampleIcs_right.size() << endl;

  // Map all feature categories to the corresponding samples and represent it as map. The map is used to assign samples to left and right branches
  map<num_t,vector<size_t> > fmap_right;
  map<num_t,vector<size_t> > fmap_left;

  size_t n_tot = 0;
  datadefs::map_data(fv,fmap_right,n_tot);
  size_t n_right = n_tot;
  size_t n_left = 0;

  if(n_tot < GI.minNodeSizeToStop) {
    return;
  }

  size_t nCategories = fmap_right.size();
  size_t nCategoriesMoved = 0; 

  bool keepSplitting = true;

  if ( treedata->isFeatureNumerical(targetIdx) ) {

    num_t mu_right;
    num_t mu_left = 0.0;
    num_t se_right;
    num_t se_left = 0.0;
    datadefs::sqerr(tv,mu_right,se_right,n_right);
    assert(n_tot == n_right);
    num_t se_best = se_right;
    num_t se_tot = se_right;

    while ( nCategoriesMoved < nCategories - 1 && keepSplitting ) {

      map<num_t,vector<size_t> >::iterator it( fmap_right.begin() );

      // We test each category one by one and see if the fitness becomes improved
      while ( it != fmap_right.end() ) {

        // Take samples from right and put them left
        //cout << "from right to left: [";
        for(size_t i = 0; i < it->second.size(); ++i) {
          //cout << " " << it->second[i];
	  datadefs::forward_backward_sqerr(tv[ it->second[i] ],n_left,mu_left,se_left,n_right,mu_right,se_right);
        }
        //cout << " ]" << endl;

        //If the fitness becomes improved, make the proposed change, otherwise move samples back
        if ( se_left + se_right < se_best ) { //&& n_left >= GI.minNodeSizeToStop && n_right >= GI.minNodeSizeToStop )

          fmap_left.insert( *it );
          fmap_right.erase( it->first );

	  ++nCategoriesMoved;

	  se_best = se_left + se_right;

	  splitFitness = ( se_tot - se_best ) / se_tot;

          break;

        } else {

          //Take samples from left and put them right
          //cout << "From left to right: [";
          for(size_t i = 0; i < it->second.size(); ++i) {
            //cout << " " << it->second[i];
	    datadefs::forward_backward_sqerr(tv[ it->second[i] ],n_right,mu_right,se_right,n_left,mu_left,se_left);
          }
          //cout << " ]" << endl;

        }

        ++it;

        if ( it == fmap_right.end() || nCategories == 2 ) {
          keepSplitting = false;
	  break;
        }

      }

    }

  } else {

    map<num_t,size_t> freq_left,freq_right;
    size_t sf_left = 0;
    size_t sf_right;
    datadefs::sqfreq(tv,freq_right,sf_right,n_right);
    assert(n_tot == n_right);
    
    num_t sf_tot = sf_right;
    num_t nsf_best = 1.0 * sf_right / n_right;
    
    while ( nCategoriesMoved < nCategories - 1 && keepSplitting ) {
      
      map<num_t,vector<size_t> >::iterator it( fmap_right.begin() );

      // We test each category one by one and see if the fitness becomes improved
      while ( it != fmap_right.end() ) {
	
	// Take samples from right and put them left
	//cout << "from right to left: [";
	for(size_t i = 0; i < it->second.size(); ++i) {
	  //cout << " " << it->second[i];
	  datadefs::forward_backward_sqfreq(tv[ it->second[i] ],n_left,freq_left,sf_left,n_right,freq_right,sf_right);
	}
	//cout << " ]" << endl;
	
	//If the fitness becomes improved, make the proposed change, otherwise move samples back
	if ( 1.0*n_right*sf_left + 1.0*n_left*sf_right > 1.0*n_left*n_right*nsf_best ) { //&& n_left >= GI.minNodeSizeToStop && n_right >= GI.minNodeSizeToStop )
	  
	  nsf_best = 1.0*sf_left/n_left + 1.0*sf_right/n_right;
	  // splitFitness = ( -1.0 * n_left*n_right*sf_tot + n_tot*n_right*sf_left + n_tot*n_left*sf_right ) / ( n_left*n_right * (1.0*n_tot*n_tot - sf_tot) );

	  fmap_left.insert( *it );
	  fmap_right.erase( it->first );
	  
	  ++nCategoriesMoved;

	  splitFitness = this->getSplitFitness(n_left,sf_left,n_right,sf_right,n_tot,sf_tot); 
	  //( -1.0 * n_left*n_right*sf_tot + n_tot*n_right*sf_left + n_tot*n_left*sf_right ) / ( n_left*n_right * (1.0*n_tot*n_tot - sf_tot) );

	  break;
	  
	} else {
	  
	  //Take samples from left and put them right
	  //cout << "From left to right: [";
	  for(size_t i = 0; i < it->second.size(); ++i) {
	    //cout << " " << it->second[i];
	    datadefs::forward_backward_sqfreq(tv[ it->second[i] ],n_right,freq_right,sf_right,n_left,freq_left,sf_left);
	  }
	  //cout << " ]" << endl;
	  
	}
	
	++it;

	if ( it == fmap_right.end() || nCategories == 2 ) {
          keepSplitting = false;
          break;
        }

	
      }

    }
    
  } 

  //cout << "C " << nCategoriesMoved << ": " << n_left << " <-> " << n_right << " : fitness " << splitFitness << endl;

  if( nCategoriesMoved == 0 || n_left < GI.minNodeSizeToStop || n_right < GI.minNodeSizeToStop ) {
    return;
  }

  // Assign samples and categories on the left. First store the original sample indices 
  vector<size_t> sampleIcs = sampleIcs_right;

  assert( n_left + n_right == n_tot );

  // Then populate the left side (sample indices and split values)
  sampleIcs_left.resize(n_left);
  splitValues_left.clear();
  size_t iter = 0;
  for ( map<num_t,vector<size_t> >::const_iterator it(fmap_left.begin()); it != fmap_left.end(); ++it ) {
    for ( size_t i = 0; i < it->second.size(); ++i ) {
      sampleIcs_left[iter] = sampleIcs[it->second[i]];
      ++iter;
    }
    splitValues_left.insert( it->first );
  }
  assert( iter == n_left);
  assert( splitValues_left.size() == fmap_left.size() );

  // Last populate the right side (sample indices and split values) 
  sampleIcs_right.resize(n_right);
  splitValues_right.clear();
  iter = 0;
  for ( map<num_t,vector<size_t> >::const_iterator it(fmap_right.begin()); it != fmap_right.end(); ++it ) {
    for ( size_t i = 0; i < it->second.size(); ++i ) {
      sampleIcs_right[iter] = sampleIcs[it->second[i]];
      ++iter;
    }
    splitValues_right.insert( it->first );
  }
  assert( iter == n_right );
  assert( splitValues_right.size() == fmap_right.size() );

}


void Node::leafMean(const vector<datadefs::num_t>& data) {
  
  if ( !datadefs::isNAN(trainPrediction_) ) {
    cerr << "Tried to set node prediction twice!" << endl;
    exit(1);
  }

  size_t n = data.size();
  assert(n > 0);
  trainPrediction_ = 0.0;

  for(size_t i = 0; i < data.size(); ++i) {
    trainPrediction_ += data[i];
  }

  trainPrediction_ /= n;

}

void Node::leafMode(const vector<datadefs::num_t>& data) {

  if ( !datadefs::isNAN(trainPrediction_) ) {
    cerr << "Tried to set node prediction twice!" << endl;
    exit(1);
  }

  size_t n = data.size();
  assert(n > 0);
  trainPrediction_ = 0.0;

  map<num_t,size_t> freq;
  
  datadefs::count_freq(data,freq,n);
  map<num_t,size_t>::iterator it(max_element(freq.begin(),freq.end(),datadefs::freqIncreasingOrder()));
  trainPrediction_ = it->first;

}

// !! Document
void Node::leafGamma(const vector<datadefs::num_t>& data, const size_t numClasses) {

  if ( !datadefs::isNAN(trainPrediction_) ) {
    cerr << "Tried to set node prediction twice!" << endl;
    exit(1);
  }

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
}
  
