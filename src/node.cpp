#include<iostream>
#include<cassert>
#include<iomanip>
#include "node.hpp"
#include "math.hpp"

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

  vector<num_t> trainData = treeData->getFeatureData(targetIdx,sampleIcs);
  if ( GI.predictionFunctionType == MEAN ) {
    trainPrediction_ = math::mean(trainData);
  } else if ( GI.predictionFunctionType == MODE ) {
    trainPrediction_ = math::mode<num_t>(trainData);
  } else if ( GI.predictionFunctionType == GAMMA ) {
    trainPrediction_ = math::gamma(trainData,GI.numClasses);
  } else {
    cerr << "Node::recursiveNodeSplit() -- unknown prediction function!" << endl;
    exit(1);
  }
  
  assert( !datadefs::isNAN(trainPrediction_) );
  
  if(nSamples < 2 * GI.minNodeSizeToStop || *nNodes >= GI.maxNodesToStop) {
    //cout << "Too few samples to start with, quitting" << endl;
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
    return;
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
  
  // This many features will be tested for splitting the data
  size_t nFeaturesForSplit = featureSampleIcs.size();
  
  num_t splitValue = datadefs::NUM_NAN;
  set<num_t> splitValues_left; 
  set<num_t> splitValues_right;

  // Initialize split fitness to lowest possible value
  splitFitness = 0.0;

  // Loop through candidate splitters
  for ( size_t i = 0; i < nFeaturesForSplit; ++i ) {
    
    // Get split feature index
    size_t newSplitFeatureIdx = featureSampleIcs[i];

    // Get type of the splitter
    bool isFeatureNumerical = treeData->isFeatureNumerical(newSplitFeatureIdx);

    //If the splitter equals to the target feature, skip splitting
    if ( newSplitFeatureIdx == targetIdx ) {
      continue;
    }

    vector<size_t> newSampleIcs_left(0);
    vector<size_t> newSampleIcs_right = sampleIcs;
    num_t newSplitValue;
    set<num_t> newSplitValues_left;
    set<num_t> newSplitValues_right;
    num_t newSplitFitness = 0.0;

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
  
  // If none of the splitter candidates worked as a splitter
  if ( fabs(splitFitness) < datadefs::EPS ) {
    return(false);
  }
  
  if ( treeData->isFeatureNumerical(splitFeatureIdx) ) {

    this->setSplitter(splitFeatureIdx,treeData->getFeatureName(splitFeatureIdx),splitValue);

  } else {
    
    //splitValues_left = newSplitValues_left[splitFeatureIdx];
    //splitValues_right = newSplitValues_right[splitFeatureIdx];

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

  //splitFitness = datadefs::NUM_NAN;

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
  
  //If the target is numerical, we use the incremental squared error formula
  if ( treedata->isFeatureNumerical(targetIdx) ) {

    // We start with one sample on the left branch
    size_t n = 1;
    num_t mu_left = tv[0];
    vector<num_t> se_left(n_tot, 0.0);

    // Add one sample at a time, from right to left, while updating the 
    // mean and squared error
    // NOTE1: leftmost sample has been added already, so we start i from 1
    // NOTE2: rightmost sample needs to be included since we want to compute
    //        squared error for the whole data vector ( i.e. se_tot )
    for( size_t i = 1; i < n_tot; ++i ) {

      // We need to manually transfer the previous squared error to the 
      // present slot
      se_left[i] = se_left[i-1];

      // This function takes n'th data point tv[i] from right to left
      // and updates the mean and squared error
      math::incrementSquaredError(tv[i],++n,mu_left,se_left[i]);
    }

    // Make sure we accumulated the right number of samples
    assert( n == n_tot );

    // Make sure the squared error didn't become corrupted by NANs
    assert( !datadefs::isNAN(se_left.back()) );

    // Total squared error now equals to the squared error of left branch
    // wince all samples were taken to the left one by one
    num_t se_tot = se_left.back();

    // The best squared error is set to the worst case
    num_t se_best = se_left.back();

    // Now it's time to start adding samples from left to right
    // NOTE: it's intentional to set n to 0
    n = 0;
    num_t mu_right = 0.0;
    num_t se_right = 0.0;
    
    // Add samples one by one from left to right until we hit the 
    // minimum allowed size of the branch
    for( int i = n_tot-1; i >= static_cast<int>(GI.minNodeSizeToStop); --i ) {
      
      // Add n'th sample tv[i] from left to right and update 
      // mean and squared error
      math::incrementSquaredError(tv[i],++n,mu_right,se_right);
 
      // If the sample is repeated and we can continue, continue
      if ( i-1 >= static_cast<int>(GI.minNodeSizeToStop) && tv[i-1] == tv[i] ) {
	continue;
      }

      // If the split point "i-1" yields a better split than the previous one,
      // update se_best and bestSplitIdx
      if ( se_left[i-1] + se_right < se_best ) {

	bestSplitIdx = i-1;
	se_best = se_left[i-1] + se_right;
	
      }
      
    }
    
    // Make sure there were no NANs to corrupt the results
    assert( !datadefs::isNAN(se_right) );
  
    // Calculate split fitness
    splitFitness = (se_tot - se_best) / se_tot;
  
  } else { // Otherwise we use the iterative gini index formula to update impurity scores while we traverse "right"

    // We start with one sample on the left branch
    size_t n_left = 1;
    map<num_t,size_t> freq_left;
    freq_left[ tv[0] ] = 1;
    vector<size_t> sf_left(n_tot, 0);
    sf_left[0] = 1;

    // Add one sample at a time, from right to left, while updating the
    // squared frequency
    // NOTE1: leftmost sample has been added already, so we start i from 1
    // NOTE2: rightmost sample needs to be included since we want to compute
    //        squared error for the whole data vector ( i.e. sf_tot )
    for( size_t i = 1; i < n_tot; ++i ) {

      // We need to manually transfer the previous squared frequency to the
      // present slot
      sf_left[i] = sf_left[i-1];

      // This function takes n'th data point tv[i] from right to left
      // and updates the squared frequency
      math::incrementSquaredFrequency(tv[i],freq_left,sf_left[i]);
      ++n_left;
    }

    // Make sure we accumulated the right number of samples
    assert( n_left == n_tot );

    // Total squared frequency now equals to the squared frequency of left branch
    // since all samples were taken to the left one by one
    size_t sf_tot = sf_left.back();

    // The best normalized squared frequency is set to the worst case
    num_t nsf_best = sf_left.back() / n_left;

    // Now it's time to start adding samples from right to left
    // NOTE: it's intentional to set n_right to 0
    size_t n_right = 0;
    map<num_t,size_t> freq_right;
    size_t sf_right = 0;

    // Add samples one by one from left to right until we hit the
    // minimum allowed size of the branch
    for( int i = n_tot-1; i >= static_cast<int>(GI.minNodeSizeToStop); --i ) {

      // Add n'th sample tv[i] from right to left and update
      // mean and squared frequency
      math::incrementSquaredFrequency(tv[i],freq_right,sf_right);
      ++n_right;
      --n_left;

      // If we have repeated samples and can continue, continue
      if ( i-1 > static_cast<int>(GI.minNodeSizeToStop) && tv[i-1] == tv[i] ) {
	continue;
      }

      // If the split point "i-1" yields a better split than the previous one,
      // update se_best and bestSplitIdx
      if(1.0*n_right*sf_left[i-1] + 1.0*n_left*sf_right > n_left*n_right*nsf_best && n_left >= GI.minNodeSizeToStop) {
        bestSplitIdx = i-1;
        nsf_best = 1.0*sf_left[i-1] / n_left + 1.0*sf_right / n_right;
	splitFitness = this->getSplitFitness(n_left,sf_left[i-1],n_right,sf_right,n_tot,sf_tot);
      }


    }

  }
  
  if(bestSplitIdx == -1) {
    //cout << "N : " << n_left << " <-> " << n_right << " : fitness " << splitFitness << endl;
    splitFitness = datadefs::NUM_NAN;
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

  if(n_tot < 2 * GI.minNodeSizeToStop) {
    splitFitness = datadefs::NUM_NAN;
    return;
  }

  if ( treedata->isFeatureNumerical(targetIdx) ) {

    num_t mu_right = 0.0;
    num_t mu_left = 0.0;
    num_t se_right = 0.0;
    num_t se_left = 0.0;

    math::squaredError(tv,mu_right,se_right);

    //    for ( size_t i = 0; i < n_tot; ++i ) {
    // math::incrementSquaredError(tv[i],i+1,mu_right,se_right);
    //}

    //    datadefs::sqerr(tv,mu_right,se_right,n_right);
    assert(n_tot == n_right);
    num_t se_best = se_right;
    num_t se_tot = se_right;
    
    while ( fmap_right.size() > 1 ) {
      
      map<num_t,vector<size_t> >::iterator it_best( fmap_right.end() );
      
      // We test each category one by one and see if the fitness becomes improved
      for ( map<num_t,vector<size_t> >::iterator it( fmap_right.begin() ); it != fmap_right.end() ; ++it ) {

	//cout << "Testing to split with feature '" << treedata->getRawFeatureData(featureIdx,it->first) << "'" << endl;
	
        // Take samples from right and put them left
        //cout << "from right to left: [";
        for(size_t i = 0; i < it->second.size(); ++i) {
          //cout << " " << it->second[i];
	  
	  // Add sample to left
	  ++n_left;
	  math::incrementSquaredError(tv[ it->second[i] ], n_left, mu_left, se_left);

	  // Remove sample from right
	  --n_right;
	  math::decrementSquaredError(tv[ it->second[i] ], n_right, mu_right, se_right);

        }
        //cout << " ]" << endl;
	
        //If the split reduces impurity even further, save the point
        if ( se_left + se_right < se_best ) { //&& n_left >= GI.minNodeSizeToStop && n_right >= GI.minNodeSizeToStop )
	  
	  //cout << " => BETTER" << endl;
	  it_best = it;
	  se_best = se_left + se_right;
	}
	
	//Take samples from left and put them right
	//cout << "From left to right: [";
	for(size_t i = 0; i < it->second.size(); ++i) {
	  //cout << " " << it->second[i];
	  
	  // Add sample to right
	  ++n_right;
	  math::incrementSquaredError(tv[ it->second[i] ], n_right, mu_right, se_right);

	  // Remove sample from left
	  --n_left;
	  math::decrementSquaredError(tv[ it->second[i] ], n_left, mu_left, se_left);

	}
	//cout << " ]" << endl;
	
      }
      
      // After testing all categories,
      // if we couldn't find any split that would reduce impurity,
      // we'll exit the loop
      if ( it_best == fmap_right.end() ) {
	//cout << " -- STOP --" << endl;
	break;
      }

      // Otherwise move samples from right to left
      for(size_t i = 0; i < it_best->second.size(); ++i) {
	//cout << " " << it->second[i];

	// Add sample to left
	++n_left;
	math::incrementSquaredError(tv[ it_best->second[i] ], n_left, mu_left, se_left);

	// Remove sample from right
	--n_right;
	math::decrementSquaredError(tv[ it_best->second[i] ], n_right, mu_right, se_right);

      }

      // Update the maps
      fmap_left.insert( *it_best );
      fmap_right.erase( it_best->first );

    }

    // Calculate the final split fitness
    splitFitness = ( se_tot - se_best ) / se_tot;    
    
  } else {
    
    map<num_t,size_t> freq_left,freq_right;
    size_t sf_left = 0;
    size_t sf_right = 0;

    for( size_t i = 0; i < n_tot; ++i ) {
      math::incrementSquaredFrequency(tv[i], freq_right, sf_right);
    }

    //datadefs::sqfreq(tv,freq_right,sf_right,n_right);
    assert(n_tot == n_right);
    
    size_t sf_tot = sf_right;
    num_t nsf_best = 1.0 * sf_right / n_right;
    
    while ( fmap_right.size() > 1 ) {
      
      map<num_t,vector<size_t> >::iterator it_best( fmap_right.end() );
      //cout << "There are " << fmap_right.size() << " categories on right" << endl;

      // We test each category one by one and see if the fitness becomes improved
      for ( map<num_t,vector<size_t> >::iterator it( fmap_right.begin() ); it != fmap_right.end() ; ++it ) {

	//cout << "Testing to split with feature '" << treedata->getRawFeatureData(featureIdx,it->first) << "'" << endl;
	
	// Take samples from right and put them left
	//cout << "from right to left: [";
	for(size_t i = 0; i < it->second.size(); ++i) {
	  
	  // Add sample to left
	  ++n_left;
	  math::incrementSquaredFrequency(tv[ it->second[i] ], freq_left, sf_left);

	  // Remove sample from right
	  --n_right;
	  math::decrementSquaredFrequency(tv[ it->second[i] ], freq_right, sf_right);

	}
	//cout << " ]" << endl;
	
	//If the impurity becomes reduced even further, save the point
	if ( 1.0*n_right*sf_left + 1.0*n_left*sf_right > 1.0*n_left*n_right*nsf_best ) { //&& n_left >= GI.minNodeSizeToStop && n_right >= GI.minNodeSizeToStop )
	
	  nsf_best = 1.0*sf_left/n_left + 1.0*sf_right/n_right;
	  it_best = it;
	  //cout << "nsf_best is now " << nsf_best << endl;
	}
	
	// Take samples from left and put them right
	//cout << "From left to right: [";
	for(size_t i = 0; i < it->second.size(); ++i) {
	  
	  // Add sample to right
	  ++n_right;
	  math::incrementSquaredFrequency(tv[ it->second[i] ], freq_right, sf_right);

	  // Remove sample from left
	  --n_left;
	  math::decrementSquaredFrequency(tv[ it->second[i] ], freq_left, sf_left);

	}
	//cout << " ]" << endl;
	
      }
      
      // After testing all categories,
      // if we couldn't find any split that would reduce impurity,
      // we'll exit the loop      
      if ( it_best == fmap_right.end() ) {
	//cout << " -- STOP --" << endl;
	break;
      }

      // Take samples from right and put them left
      for(size_t i = 0; i < it_best->second.size(); ++i) {
        
	// Add sample to left
	++n_left;
	math::incrementSquaredFrequency(tv[ it_best->second[i] ], freq_left, sf_left);
	
	// Remove sample from right
	--n_right;
	math::decrementSquaredFrequency(tv[ it_best->second[i] ], freq_right, sf_right);
	
      }
      
      // Update the maps
      fmap_left.insert( *it_best );
      fmap_right.erase( it_best->first );
      
    }
    
    // Calculate the final split fitness
    splitFitness = this->getSplitFitness(n_left,sf_left,n_right,sf_right,n_tot,sf_tot);
    
  } 

  if( n_left < GI.minNodeSizeToStop || n_right < GI.minNodeSizeToStop ) {
    splitFitness = datadefs::NUM_NAN;
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

/*
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
*/

/*
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
*/

// !! Document
/*
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
*/
  
