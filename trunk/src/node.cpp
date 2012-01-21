#include<iostream>
#include<cassert>
#include<iomanip>
#include "node.hpp"

Node::Node():
  splitterIdx_(0),
  splitter_(NULL),
  //isTrainPredictionSet_(false),
  trainPrediction_(datadefs::NUM_NAN),
  nTestSamples_(0),
  testPredictionError_(0.0),
  //hasChildren_(false),
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
  
  if ( leftChild_ || rightChild_ ) {
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

  if ( leftChild_ || rightChild_ ) {
    cerr << "Cannot set a splitter to a node twice!" << endl;
    exit(1);
  }

  splitterIdx_ = splitterIdx;
  splitter_ = new Splitter(splitterName,leftSplitValues,rightSplitValues);

  leftChild_ = new Node();
  rightChild_ = new Node();

}


Node* Node::percolateData(const num_t data) {

  //if ( datadefs::isNAN( data ) ) {
  //  cerr << "Node::percolateData(num_t) does not accept NaNs (" << data << ") to be percolated!" << endl;
  //  exit(1);
  //}

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

  //if ( datadefs::isNAN( data ) ) {
  //  cerr << "Node::percolateData(string) does not accept NaNs (" << data << ") to be percolated!" << endl;
  //  exit(1);
  //}

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
  
  Node* nodep(this);

  this->percolateData(treeData,sampleIdx,&nodep);

  return( nodep );

}

void Node::percolateData(Treedata* treeData, const size_t sampleIdx, Node** nodep) {

  if( !this->hasChildren() ) {
    *nodep = this;
  } else if ( splitter_->splitsLeft(treeData,sampleIdx) ) {
    *nodep = leftChild_;
    leftChild_->percolateData(treeData,sampleIdx,nodep);
  } else if ( splitter_->splitsRight(treeData,sampleIdx) ) {
    *nodep = rightChild_;
    rightChild_->percolateData(treeData,sampleIdx,nodep);
  } else {
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

  if ( !this->hasChildren() ) {
    //toFile << "\tLEAF" << endl;
    return;
  }

  string traversalLeft = traversal;
  traversalLeft.append("L");
  string traversalRight = traversal;
  traversalRight.append("R");

  bool isLeftNodeLeaf = false;
  bool isRightNodeLeaf = false;

  if ( !leftChild_->hasChildren() ) {
    isLeftNodeLeaf = true;
    traversalLeft.append("*");
  }

  if ( !rightChild_->hasChildren() ) {
    isRightNodeLeaf = true;
    traversalRight.append("*");
  }

  // TODO: Node::print() needs to be completed
  toFile << ",L_PRED=" << setprecision(3) << leftChild_->getTrainPrediction() << ",R_PRED=" 
	 << setprecision(3) << rightChild_->getTrainPrediction() << ",SPLITTER=" << splitter_->name() << endl;

  if ( !isLeftNodeLeaf ) {
    leftChild_->print(traversalLeft,toFile);
  }

  if ( !isRightNodeLeaf ) {
    rightChild_->print(traversalRight,toFile);
  }
}

void Node::print(ofstream& toFile) {
  
  string traversal("");
  this->print(traversal,toFile);

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
      featureSampleIcs[i] = i;
    }
  }

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


  vector<num_t> tv,fv;

  sampleIcs_left.clear();

  treedata->getFilteredAndSortedFeatureDataPair3(targetIdx,featureIdx,sampleIcs_right,tv,fv);

  //datadefs::print(tv);
  //datadefs::print(fv);

  //treedata->getFilteredFeatureDataPair(targetIdx,featureIdx,sampleIcs_right,tv,fv);
  //vector<num_t> fv = treedata->getFilteredFeatureData(featureIdx,sampleIcs_right);

  size_t n_tot = fv.size();
  size_t n_right = n_tot;
  size_t n_left = 0;

  if(n_tot < 2 * GI.minNodeSizeToStop) {
    splitFitness = datadefs::NUM_NAN;
    return;
  }

  //Make reference indices that define the sorting wrt. feature
  //bool isIncreasingOrder = true;

  //Sort feature vector and collect reference indices
  //vector<size_t> refIcs;
  //datadefs::sortDataAndMakeRef(isIncreasingOrder,fv,refIcs);

  //datadefs::sortFromRef<size_t>(sampleIcs_right,refIcs);

  //vector<num_t> tv = treedata->getFeatureData(targetIdx,sampleIcs_right);

  //Use the reference indices to sort sample indices
  //datadefs::sortFromRef<num_t>(tv,refIcs);
  //datadefs::sortFromRef<size_t>(sampleIcs_right,refIcs);

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
    sampleIcs_left[i] = sampleIcs_right[i];
  }
  sampleIcs_right.erase(sampleIcs_right.begin(),sampleIcs_right.begin() + n_left);
  n_right = sampleIcs_right.size();

  assert(n_left + n_right == n_tot);

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

  vector<num_t> tv,fv;

  sampleIcs_left.clear();
  treedata->getFilteredFeatureDataPair(targetIdx,featureIdx,sampleIcs_right,tv,fv);

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

  if ( treedata->isFeatureNumerical(targetIdx) ) {

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

      if ( se_left+se_right < se_best && n_left >= GI.minNodeSizeToStop && n_right >= GI.minNodeSizeToStop ) {
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
      
      if ( 1.0*n_right*sf_left + n_left*sf_right > n_left*n_right*nsf_best && n_left >= GI.minNodeSizeToStop && n_right >= GI.minNodeSizeToStop ) {
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

  // Assign samples and categories on the left. First store the original sample indices 
  vector<size_t> sampleIcs = sampleIcs_right;

  sampleIcs_left.resize(n_tot);
  splitValues_left.clear();

  // Then populate the left side (sample indices and split values)
  iter = 0;
  for ( map<num_t,vector<size_t> >::const_iterator it(fmap_left_best.begin()); it != fmap_left_best.end(); ++it ) {
    for ( size_t i = 0; i < it->second.size(); ++i ) {
      sampleIcs_left[iter] = sampleIcs[it->second[i]];
      ++iter;
    }
    splitValues_left.insert( it->first );
  }
  sampleIcs_left.resize(iter);

  // Last populate the right side (sample indices and split values) 
  iter = 0;
  for ( map<num_t,vector<size_t> >::const_iterator it(fmap_right_best.begin()); it != fmap_right_best.end(); ++it ) {
    for ( size_t i = 0; i < it->second.size(); ++i ) {
      sampleIcs_right[iter] = sampleIcs[it->second[i]];
      ++iter;
    }
    splitValues_right.insert( it->first );
  }
  sampleIcs_right.resize(iter);

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
  
