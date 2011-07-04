#include "GBT.hpp"
#include <cmath>
#include <iostream>

GBT::GBT(Treedata* treeData, size_t nTrees, size_t nMaxLeaves, float shrinkage, float subSampleSize):
treeData_(treeData),
nTrees_(nTrees),
nMaxLeaves_(nMaxLeaves),
shrinkage_(shrinkage),
subSampleSize_(subSampleSize)
{
	numClasses_ = 0; // Denotes regression forest by default
	nodeSize_ = 1;   // Allow split down to individual example
	// Currently, we allocate and grow the forest outside of the constructor
}




GBT::~GBT()
{
}


// reserve memory
void GBT::allocateForest()
{
	// Determine the number of nodes we need for the GBT-forest.
	//
	// nMaxLeaves determines the max number of leaf nodes
	// The depth of the last complete level of the tree:
	int depthComplete = int(floor(log10(float(nMaxLeaves_))/log10(2.0)));
	// Number of nodes of the last complete level of the tree:
	int nTerminalNodesComplete = int(pow(2.0, depthComplete));
	//
	// The number of nodes in a complete binary tree of depth depthComplete
	// plus the last incomplete level
	nMaxNodes_ = 2*nTerminalNodesComplete-1  + 2*( nMaxLeaves_-nTerminalNodesComplete );

	// for classification, we need the number of classes, 
	// because it will multiply the number of trees
	size_t iTarget = GBT::getTarget();
	if ( !treeData_->isFeatureNumerical( iTarget ) )
	{
		map<num_t,size_t> cat2freq;
		size_t nReal;
		datadefs::count_freq( treeData_->features_[iTarget].data, cat2freq, nReal);
		numClasses_ = cat2freq.size(); // set numClasses_ to denote a classification forest
		cout << "Target "<<iTarget<<" has "<<numClasses_<<" classes, "<<nReal<<" non-missing."<<endl;
		nTrees_ = nTrees_ * numClasses_;
		cout << "Target "<<iTarget<<" has "<<numClasses_<<" classes, "<<nReal<<" non-missing."<<endl;
	}

	// Reserve memory for vector of node counts in each tree
	nNodes_ = vector<size_t>(nTrees_, 0);
	// Reserve memory for the forest as a vector of vector of nodes
	forest_ = vector< vector<Node> >(nTrees_, vector<Node>( nMaxNodes_ ));

	cout << "Forest initialized. " << forest_.size() << " trees and "
			<< nMaxNodes_ << " max nodes per tree generated." << endl;
}




// Grow a GBT "forest"
void GBT::growForest()
{
	if ( treeData_->isFeatureNumerical( GBT::getTarget() ) )
		GBT::growForestNumerical();
	else
		GBT::growForestCategorical();
}



// Grow a GBT "forest" for a numerical target variable
void GBT::growForestNumerical()
{
	//---------------------------------------------------------------------------
	// The target attribute in GBT should have been set beforehand
	size_t iTarget  = getTarget();
	size_t nSamples = treeData_->nSamples();
	// save a copy of the target column because it will be overwritten
	vector<num_t> trueTarget = treeData_->features_[iTarget].data;

	// Target for GBT is different for each tree
	// reference to the target column, will overwrite it
	vector<num_t>& curTarget = treeData_->features_[iTarget].data;


	//---------------------------------------------------------------------------
	// Set the initial prediction to zero.
	vector<num_t> prediction(nSamples, 0.0);
	vector<num_t> curPrediction(nSamples);

	//---------------------------------------------------------------------------
	for(size_t t = 0; t < nTrees_; ++t)
	{
		cout << endl;

		// current target is the negative gradient of the loss function
		// for square loss, it is target minus current prediction
		for (size_t i=0; i<nSamples; i++)
		{
			curTarget[i] = trueTarget[i] - prediction[i];
		}

		// Grow a tree to predict the current target
		GBT::growTree(t);

		// What kind of a prediction does the new tree produce?
		GBT::predictDatasetByTree(t, curPrediction);

		// Calculate current total prediction adding the newly generated tree
		for (size_t i=0; i<nSamples; i++)
		{
			prediction[i] = prediction[i] + shrinkage_ * curPrediction[i];
		}
	}
	// GBT-forest is now done!

	// restore true target column
	treeData_->features_[iTarget].data = trueTarget;
}


// Grow a GBT "forest" for a categorical target variable
void GBT::growForestCategorical()
{
	//---------------------------------------------------------------------------
	// The target attribute in GBT should have been set beforehand
	size_t iTarget  = getTarget();
	size_t nSamples = treeData_->nSamples();
	
	// Save a copy of the target column because it will be overwritten later.
	// We also know that it must be categorical.
	vector<num_t> trueTarget = treeData_->features_[iTarget].data;

	// Target for GBT is different for each tree.
	// We use the original target column to save each temporary target.
	// Reference to the target column, will overwrite it:
	vector<num_t>& curTarget = treeData_->features_[iTarget].data;

	// For categorical target variables, we predict each category probability separately.
	// This target is numerical, thus we need to change the target variable type to 
	// numerical from categorical.
	treeData_->features_[iTarget].isNumerical = true;
	
	
	// Initialize class probability estimates and the predictions.
	// Note that dimensions in these two are reversed!
	vector< vector<num_t> > prediction(    nSamples, vector<num_t>( numClasses_, 0.0 ) );
	vector< vector<num_t> > curPrediction( numClasses_, vector<num_t>( nSamples, 0.0 ) );
	
	vector< vector<num_t> > curProbability( nSamples, vector<num_t>( numClasses_)  );
	// TODO don't really need all this space allocated
		
	
	// Each iteration consists of numClasses_ trees, 
	// each of those predicting the probability residual for each class.
	size_t numIterations = nTrees_ / numClasses_;
	
	for(size_t m=0; m < numIterations; ++m)
	{	
		// Multiclass logistic transform of class probabilities from current probability estimates.
		for (size_t i=0; i<nSamples; i++)
		{		
			GBT::transformLogistic( prediction[i],  curProbability[i]); 
			// each prediction[i] is a vector<num_t>(numClasses_)
		}		
		
		// construct a tree for each class
		for (size_t k = 0; k < numClasses_; ++k)
		{
			// target for class k is ...
			for (size_t i=0; i<nSamples; i++)
			{
				// ... the difference between true target and current prediction
				curTarget[i] = (k==trueTarget[i]) - curProbability[i][k];
			}

			// Grow a tree to predict the current target
			size_t t = m*numClasses_ + k; // tree index
			GBT::growTree(t);

			// What kind of a prediction does the new tree produce?
			GBT::predictDatasetByTree(t, curPrediction[k] );

			// Calculate the current total prediction adding the newly generated tree
			for (size_t i=0; i<nSamples; i++)
			{
				prediction[i][k] = prediction[i][k] + shrinkage_ * curPrediction[k][i];
			}
		}
	}
		
	// GBT-forest is now done!
	// restore the true target column and its true type
	treeData_->features_[iTarget].data = trueTarget;
	treeData_->features_[iTarget].isNumerical = false;
}



void GBT::transformLogistic(vector<num_t>& prediction, vector<num_t>& probability)
{
	// Multiclass logistic transform of class probabilities from current probability estimates.
	assert( numClasses_ == prediction.size() );	
	vector<num_t>& expPrediction = probability; // just using the space by a different name
	
	// find maximum prediction
	vector<num_t>::iterator maxPrediction = max_element( prediction.begin(), prediction.end() );
	// scale by maximum to prevent numerical errors
	


	num_t expSum = 0.0;
	size_t k;
	for(k=0; k < numClasses_; ++k)
	{
		expPrediction[k] = exp( prediction[k] - *maxPrediction ); // scale by maximum
		expSum += expPrediction[k];		
	}
	for(k = 0; k < numClasses_; ++k)
	{
		probability[k] = expPrediction[k] / expSum;	
	}			
}



// Predict using a GBT "forest" from an arbitrary data set. GBT::numClasses_ determines 
// whether the prediction is for a categorical or a numerical variable.
void GBT::predictForest(Treedata* treeData, vector<num_t>& prediction)
{
	if ( 0 == numClasses_ ) 
		predictForestNumerical(treeData, prediction);
	else 
		predictForestCategorical(treeData, prediction);
}


// Predict categorical target using a GBT "forest" from an arbitrary data set. 
void GBT::predictForestCategorical(Treedata* treeData, vector<num_t>& categoryPrediction)
{
	size_t nSamples = treeData->nSamples();
	size_t iTarget  = getTarget();
	cout << "Predicting "<<nSamples<<" samples. Target="<<iTarget<<". Classes="<<numClasses_<<endl;
	
	// for classification, each "tree" is actually numClasses_ trees in a row
	// each predicting the probability of its own class.
	size_t numIterations = nTrees_ / numClasses_;
	size_t errorCnt = 0;
	vector<num_t> probPrediction( numClasses_ );

	
	for (size_t i=0; i<nSamples; i++)
	{
		vector<num_t> prediction( numClasses_ , 0);
		for (size_t k=0; k<numClasses_; ++k) // 
		{
			for(size_t m = 0; m < numIterations; ++m)
			{
				size_t t =  m*numClasses_ + k; // tree index
				prediction[k] = prediction[k] + shrinkage_ * predictSampleByTree(i, t, treeData);
			}
		}
		// find index of maximum prediction, this is the predicted cateory
		vector<num_t>::iterator maxProb = max_element( prediction.begin(), prediction.end() );
		categoryPrediction[i] = maxProb - prediction.begin() ; // classes are 0,1,2,...
		
		
		
		// diagnostic print out the true and the prediction
		errorCnt += (treeData->features_[iTarget].data[i] != categoryPrediction[i]);
		cout << i << "\t" << treeData->features_[iTarget].data[i] << "\t" << categoryPrediction[i];
		GBT::transformLogistic(prediction, probPrediction); // predictions-to-probabilities
		for (size_t k=0; k<numClasses_; ++k)
		{
			cout <<"\t"<<probPrediction[k];
		}
		cout <<endl;
	}
	cout<<"Error "<<errorCnt<<"/"<<nSamples<<"="<<100*errorCnt/nSamples<<"%"<<endl;
}



// Predict numerical target using a GBT "forest" from an arbitrary data set
void GBT::predictForestNumerical(Treedata* treeData, vector<num_t>& prediction)
{
	size_t nSamples = treeData->nSamples();
	size_t iTarget  = getTarget();
	cout << "Predicting "<<nSamples<<" samples. Target="<<iTarget<<endl;
	for (size_t i=0; i<nSamples; i++)
	{
		prediction[i] = 0;
		for(size_t t = 0; t < nTrees_; ++t)
		{
			prediction[i] = prediction[i] + shrinkage_ * predictSampleByTree(i, t, treeData);
		}
		// diagnostic print out the true and the prediction
		cout << i << "\t" << treeData->features_[iTarget].data[i] << "\t" << prediction[i] <<endl;
	}
}



// Use a single GBT tree to produce a prediction from a single data sample
// of an arbitrary data set. This wouldn't have to be the train set.
num_t GBT::predictSampleByTree(size_t sampleIdx, size_t treeIdx, Treedata *treeData)
{
	// root of current tree
	Node *currentNode = &forest_[treeIdx][0];
	num_t value;
	// traverse to the leaf of the tree
	while(currentNode->hasChildren())
	{
		size_t featureIdx = currentNode->getSplitter();
		treeData->getFeatureData(featureIdx, sampleIdx, value);
		currentNode = currentNode->percolateData(value);
	}
	// now currentNode points to a leaf node: get the prediction
	return currentNode->getPrediction();
}



// Use a single GBT tree to produce predictions for a whole training data set
void GBT::predictDatasetByTree(size_t treeIdx, vector<num_t>& curPrediction)
{
	size_t nSamples = treeData_->nSamples();
	// predict for all samples
	for(size_t sampleIdx = 0; sampleIdx < nSamples; ++sampleIdx)
	{
		curPrediction[sampleIdx] = GBT::predictSampleByTree(sampleIdx, treeIdx, treeData_);

		cout << "Sample " << sampleIdx << ", prediction " << curPrediction[sampleIdx]  << endl;
	}
}



void GBT::growTree(size_t treeIdx)
{
	size_t nSamples = treeData_->nSamples();

	//Generate a new random  vector for sample indices
	vector<size_t> sampleIcs(nSamples);
	for(size_t i = 0; i < nSamples; ++i)
	{
		sampleIcs[i] = i;
	}
	
	// treeData_->permute(sampleIcs);
	// TODO implement training using only subSampleSize_ samples.
	
	cout << "Growing tree " << treeIdx << " with a sample:" << endl;

	size_t rootnode = 0;

	//Start the recursive node splitting from the root node. This will generate the tree.
	GBT::recursiveNodeSplit(treeIdx, rootnode, sampleIcs);
}



void GBT::SetNodePrediction( size_t treeIdx, size_t nodeIdx, vector<size_t>& sampleIcs)
{

	// calculate the prediction from leaf node training set
	vector<num_t> data(sampleIcs.size());
	size_t targetIdx = GBT::getTarget();
	//
	treeData_->getFeatureData(targetIdx, sampleIcs, data);
	
	num_t nodePred;
	if (0==numClasses_) 
	{
		// numerical target
		num_t se;
		size_t nreal;
		datadefs::sqerr( data, nodePred,se,nreal); // only need mu as the node Pred
		// TODO Simplify: we do not need everything that the above calculates
	}
	else
	{
		// categorical target, we are predicting one class probability, calculate gamma
		num_t numerator=0.0, denominator=0.0;
		for (size_t i = 0; i < data.size(); ++i)
		{
			num_t abs_data_i = abs( data[i] );
			denominator += abs_data_i * (1.0 - abs_data_i);
			numerator   += data[i];
		}
		if ( fabs(denominator) <= datadefs::EPS )
		{
			nodePred = LOG_OF_MAX_FLT * numerator;
		}
		else
		{
			nodePred = (numClasses_ - 1)*numerator / (numClasses_ * denominator);
		}
	}
		
	forest_[treeIdx][nodeIdx].setPrediction(nodePred);
	cout << "setting pred=" << nodePred << ", tree="<<treeIdx<<" node="<<nodeIdx<<endl;
}



void GBT::recursiveNodeSplit(size_t treeIdx, size_t nodeIdx, vector<size_t>& sampleIcs)
{
	size_t nSamples = sampleIcs.size();
	size_t nFeatures = treeData_->nFeatures();

	// exceeding max number of nodes if we split this node?
	if ( 2+nNodes_[treeIdx] >= nMaxNodes_)
	{
		cout << "tree=" << treeIdx << " nodeIdx=" << nodeIdx << " nMaxNodes=" << nMaxNodes_ << " nNodes=" << nNodes_[treeIdx] << endl;
		// Leaf node
		GBT::SetNodePrediction(treeIdx, nodeIdx, sampleIcs);
		return;
	}

	// not enough samples to split?
	if (sampleIcs.size() <= nodeSize_)
	{
		cout << "tree=" << treeIdx << "  nodeIdx=" << nodeIdx << ". Only " << sampleIcs.size() << " samples!" << endl;
		// Leaf node
		GBT::SetNodePrediction(treeIdx, nodeIdx, sampleIcs);
		return;
	}
	cout << "tree=" << treeIdx << "  nodeIdx=" << nodeIdx << endl;

	
	
	// Find the best target split
	vector<num_t> targetData;
	const bool isTargetNumerical = treeData_->isFeatureNumerical(targetIdx_);
	treeData_->getFeatureData(targetIdx_,sampleIcs,targetData);
	
	vector<size_t> sampleIcs_left, sampleIcs_right;
	num_t splitValue;
	set<num_t> values_left;	
	
	if(isTargetNumerical) 
	{
		treeData_->numericalFeatureSplit(targetData,isTargetNumerical,targetData,nodeSize_,sampleIcs_left,sampleIcs_right,splitValue);
	}
	else 
	{
		treeData_->categoricalFeatureSplit(targetData,isTargetNumerical,targetData,sampleIcs_left,sampleIcs_right,values_left);
	}
	cout << "Target split with itself." << endl;

	
	
	// Find the best feature split corresponding to the target split
	num_t bestFitness = 0.0;
	size_t bestFeatureIdx = nFeatures;
	vector<num_t> featureData;
	

	// TODO Instead of going through all features, sample features here
	// using sampling weights proportional to current importances
	for (size_t featureIdx = 0; featureIdx < nFeatures; ++featureIdx)
	{
		if (featureIdx == targetIdx_)
		{
			continue;
		}
		bool isFeatureNumerical = treeData_->isFeatureNumerical(featureIdx);
		treeData_->getFeatureData(featureIdx,sampleIcs,featureData);
		num_t fitness = treeData_->splitFitness(featureData,isFeatureNumerical,nodeSize_,sampleIcs_left,sampleIcs_right);

		if (fitness > bestFitness)
		{
			bestFitness = fitness;
			bestFeatureIdx = featureIdx;
			if (fabs(bestFitness - 1.0) < datadefs::EPS)
			{
				//cout << "Maximum fitness reached for splitter " << bestfeatureidx << ", stopped searching" << endl;
				break;
			}
		}
	}

	if (bestFeatureIdx == nFeatures)
	{
		cout << "No splitter found, quitting." << endl;
		// Leaf node
		GBT::SetNodePrediction(treeIdx, nodeIdx, sampleIcs);
		return;
	}
	cout << "Best splitter feature is " << bestFeatureIdx << " with fitness of " << bestFitness << endl; 

	
	treeData_->getFeatureData(bestFeatureIdx,sampleIcs,featureData);
	
	vector<size_t> NANIcs;

	size_t nRealSamples = 0;
	for(size_t i = 0; i < nSamples; ++i)
	{
		if(!datadefs::isNAN(featureData[i]))
		{
			featureData[nRealSamples] = featureData[i];
			targetData[nRealSamples] = targetData[i];
			++nRealSamples;
		}
	}
	featureData.resize(nRealSamples);
	targetData.resize(nRealSamples);
	
	if (nRealSamples < 2*nodeSize_)
    {
		cout << "Splitter has too few non-missing values, quitting" << endl;
		//TODO This needs to be fixed such that one of the surrogates will determine the split instead
		GBT::SetNodePrediction(treeIdx, nodeIdx, sampleIcs);
		return;
    }	
	
	
	size_t nodeIdx_left; 
	size_t nodeIdx_right;
	bool isBestFeatureNumerical = treeData_->isFeatureNumerical(bestFeatureIdx);

	if (isBestFeatureNumerical)
	{
		treeData_->numericalFeatureSplit(targetData,isTargetNumerical,featureData,nodeSize_,sampleIcs_left,sampleIcs_right,splitValue);
	}
	else
	{
		treeData_->categoricalFeatureSplit(targetData,isTargetNumerical,featureData,sampleIcs_left,sampleIcs_right,values_left);
	}

	assert(sampleIcs_left.size() + sampleIcs_right.size() == nRealSamples);

	if (sampleIcs_left.size() < nodeSize_ || sampleIcs_right.size() < nodeSize_)
	{
		// Leaf node
		GBT::SetNodePrediction(treeIdx, nodeIdx, sampleIcs);
		return;
	}

	nodeIdx_left  = ++nNodes_[treeIdx]; // we start numbering from zero
	nodeIdx_right = ++nNodes_[treeIdx];

	if (isBestFeatureNumerical)
	{
		forest_[treeIdx][nodeIdx].setSplitter(bestFeatureIdx,splitValue,forest_[treeIdx][nodeIdx_left],forest_[treeIdx][nodeIdx_right]);
	}
	else
	{
		forest_[treeIdx][nodeIdx].setSplitter(bestFeatureIdx,values_left,forest_[treeIdx][nodeIdx_left],forest_[treeIdx][nodeIdx_right]);
	}
	
	GBT::recursiveNodeSplit(treeIdx,nodeIdx_left, sampleIcs_left );
	GBT::recursiveNodeSplit(treeIdx,nodeIdx_right,sampleIcs_right);

/*
	if (sampleIcs_left.size() >= 2*nodeSize_)
	{
		//WILL BECOME DEPRECATED
		for (size_t i = 0; i < sampleIcs_left.size(); ++i)
		{
			sampleIcs_left[i] = sampleIcs[sampleIcs_left[i]];
		}

		GBT::recursiveNodeSplit(treeIdx,nodeIdx_left,sampleIcs_left);
	}

	if (sampleIcs_right.size() >= 2*nodeSize_)
	{
		//WILL BECOME DEPRECATED
		for (size_t i = 0; i < sampleIcs_right.size(); ++i)
		{
			sampleIcs_right[i] = sampleIcs[sampleIcs_right[i]];
		}

		GBT::recursiveNodeSplit(treeIdx,nodeIdx_right,sampleIcs_right);
	}
*/	
}




void GBT::setTarget(size_t targetIdx)
{
	targetIdx_ = targetIdx;
}


size_t GBT::getTarget()
{
	return(targetIdx_);
}
