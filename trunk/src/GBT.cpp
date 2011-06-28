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
	if ( !treeData_->isfeaturenum( iTarget ) )
	{
		map<num_t,size_t> cat2freq;
		size_t nReal;
		datadefs::count_freq( treeData_->featurematrix_[iTarget], cat2freq, nReal);
		numClasses_ = cat2freq.size(); // set numClasses_ to denote a classification forest
		cout << "Target "<<iTarget<<" has "<<numClasses_<<" classes, "<<nReal<<" non-missing."<<endl;
		nTrees_ = nTrees_ * numClasses_;
		cout << "Target "<<iTarget<<" has "<<numClasses_<<" classes, "<<nReal<<" non-missing."<<endl;
	}

	// Reserve memory for vector of node counts in each tree
	nNodes_ = vector<size_t>(nTrees_);
	// Reserve memory for the forest "pointer" vector
	forest_ = vector< vector<Node> >(nTrees_);
	
	// Memory for nodes of each individual tree
	for(size_t t = 0; t < nTrees_; ++t)
	{
		forest_[t] = vector<Node>(nMaxNodes_);
		nNodes_[t] = 0;
	}

	cout << "Forest initialized. " << forest_.size() << " trees and "
			<< nMaxNodes_ << " max nodes per tree generated." << endl;
}




// Grow a GBT "forest"
void GBT::growForest()
{
	if ( treeData_->isfeaturenum( GBT::getTarget() ) )
	{
		GBT::growForestNumerical();
	}
	else
	{
		GBT::growForestCategorical();
	}
}



// Grow a GBT "forest" for a numerical target variable
void GBT::growForestNumerical()
{
	//---------------------------------------------------------------------------
	// The target attribute in treedata should have been set beforehand
	size_t iTarget  = treeData_->get_target();
	size_t nSamples = treeData_->nsamples();
	// save a copy of the target column because it will be overwritten
	vector<num_t> trueTarget = treeData_->featurematrix_[iTarget];

	// Target for GBT is different for each tree
	// reference to the target column, will overwrite it
	vector<num_t>& curTarget = treeData_->featurematrix_[iTarget];


	//---------------------------------------------------------------------------
	// (Set the initial prediction as the mean of the target of the dataset)
	// ... OR set it to zero to avoid treating the 1st node differently
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


		// set target in order to re-sort
		//treeData_->select_target( iTarget );
		treeData_->targetidx_ = iTarget;

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
	treeData_->featurematrix_[iTarget] = trueTarget;
}


// Grow a GBT "forest" for a categorical target variable
void GBT::growForestCategorical()
{
	//---------------------------------------------------------------------------
	// The target attribute in treedata should have been set beforehand
	size_t iTarget  = treeData_->get_target();
	size_t nSamples = treeData_->nsamples();
	
	// Save a copy of the target column because it will be overwritten later.
	// We also know that it must be categorical.
	vector<num_t> trueTarget = treeData_->featurematrix_[iTarget];

	// Target for GBT is different for each tree.
	// We use the original target column to save each temporary target.
	// Reference to the target column, will overwrite it:
	vector<num_t>& curTarget = treeData_->featurematrix_[iTarget];

	// For categorical target variables, we predict each category probability separately.
	// This target is numerical, thus we need to change the target variable type to 
	// numerical from categorical.
	treeData_->isfeaturenum_[iTarget] = true;
	
	
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
	treeData_->featurematrix_[iTarget] = trueTarget;
	treeData_->isfeaturenum_[iTarget] = false;
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
	{
		predictForestNumerical(treeData, prediction);
	} 
	else 
	{
		predictForestCategorical(treeData, prediction);
	}
	
}


// Predict categorical target using a GBT "forest" from an arbitrary data set. 
void GBT::predictForestCategorical(Treedata* treeData, vector<num_t>& categoryPrediction)
{
	size_t nSamples = treeData->nsamples();
	size_t iTarget  = treeData->get_target();
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
		errorCnt += (treeData->featurematrix_[iTarget][i] != categoryPrediction[i]);
		cout << i << "\t" << treeData->featurematrix_[iTarget][i] << "\t" << categoryPrediction[i];
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
	size_t nSamples = treeData->nsamples();
	size_t iTarget  = treeData->get_target();
	cout << "Predicting "<<nSamples<<" samples. Target="<<iTarget<<endl;
	for (size_t i=0; i<nSamples; i++)
	{
		prediction[i] = 0;
		for(size_t t = 0; t < nTrees_; ++t)
		{
			prediction[i] = prediction[i] + shrinkage_ * predictSampleByTree(i, t, treeData);
		}
		// diagnostic print out the true and the prediction
		cout << i << "\t" << treeData->featurematrix_[iTarget][i] << "\t" << prediction[i] <<endl;
	}
}



// Use a single GBT tree to produce a prediction from a single data sample
// of an arbitrary data set. This wouldn't have to be the train set.
num_t GBT::predictSampleByTree(size_t sampleidx, size_t treeidx, Treedata *treeData)
{
	// root of current tree
	Node *currentNode = &forest_[treeidx][0];
	// traverse to the leaf of the tree
	while(currentNode->has_children())
	{
		size_t featureidx = currentNode->get_splitter();
		num_t value = treeData->at(featureidx,sampleidx);
		currentNode = currentNode->percolate(value);
	}
	// now currentNode points to a leaf node: get the prediction
	return currentNode->get_prediction();
}



// Use a single GBT tree to produce predictions for a whole training data set
void GBT::predictDatasetByTree(size_t treeidx, vector<num_t>& curPrediction)
{
	size_t nSamples = treeData_->nsamples();
	// predict for all samples
	for(size_t sampleidx = 0; sampleidx < nSamples; ++sampleidx)
	{
		curPrediction[sampleidx] = GBT::predictSampleByTree(sampleidx, treeidx, treeData_);

		cout << "Sample " << sampleidx << ", prediction " << curPrediction[sampleidx]  << endl;
	}
}



void GBT::growTree(size_t treeidx)
{
	size_t nSamples = treeData_->nsamples();

	//Generate a new random  vector for sample indices
	vector<size_t> sample_ics(nSamples);
	for(size_t i = 0; i < nSamples; ++i)
	{
		sample_ics[i] = i;
	}
	treeData_->permute(sample_ics);

	// TODO implement training using only subSampleSize_ samples.
	cout << "Growing tree " << treeidx << " with a sample:" << endl;

	size_t rootnode = 0;

	//Start the recursive node splitting from the root node. This will generate the tree.
	GBT::recursiveNodesplit(treeidx, rootnode, sample_ics);
}



void GBT::SetNodePrediction( size_t treeidx, size_t nodeidx, vector<size_t>& sampleics)
{

	// calculate the prediction from leaf node training set
	vector<num_t> data(sampleics.size());
	size_t targetidx = GBT::getTarget();
	//
	for(size_t i = 0; i < sampleics.size(); ++i)
	{
		data[i] = (treeData_->featurematrix_)[targetidx][ sampleics[i] ];
	}
	
	num_t nodePred;
	if (0==numClasses_) 
	{
		// numerical target
		num_t se;
		size_t nreal;
		datadefs::sqerr( data, nodePred,se,nreal); // only need mu as the node Pred
		// TODO we do no need everything that the above calculates
	}
	else
	{
		// categorical target, we are predicting one class probability, calculate gamma
		num_t numerator=0.0, denominator=0.0;
		for(size_t i = 0; i < data.size(); ++i)
		{
			num_t abs_data_i = abs( data[i] );
			denominator += abs_data_i * (1.0 - abs_data_i);
			numerator   += data[i];
		}
		if ( abs(denominator) <= datadefs::eps )
		{
			nodePred = LOG_OF_MAX_FLT * numerator;
		}
		else
		{
			nodePred = (numClasses_ - 1)*numerator / (numClasses_ * denominator);
		}
	}
		
	forest_[treeidx][nodeidx].set_prediction(nodePred);
	cout << "setting pred=" << nodePred << ", tree="<<treeidx<<" node="<<nodeidx<<endl;
}



void GBT::recursiveNodesplit(size_t treeidx, size_t nodeidx, vector<size_t>& sampleics)
{
	size_t n_tot = sampleics.size();
	size_t minSplit = 1; // smallest node size
	size_t nFeatures = treeData_->nfeatures();

	// exceeding max number of nodes if we split this node
	if ( 2+nNodes_[treeidx] >= nMaxNodes_)
	{
		cout << "tree=" << treeidx << " nodeidx=" << nodeidx << " nMaxNodes=" << nMaxNodes_ << " nNodes=" << nNodes_[treeidx] << endl;
		// Leaf node
		GBT::SetNodePrediction(treeidx, nodeidx, sampleics);
		return;
	}

	// not enough samples to split
	if (sampleics.size() <= minSplit)
	{
		cout << "tree=" << treeidx << "  nodeidx=" << nodeidx << ". Only " << sampleics.size() << " samples!" << endl;
		// Leaf node
		GBT::SetNodePrediction(treeidx, nodeidx, sampleics);
		return;
	}
	cout << "tree=" << treeidx << "  nodeidx=" << nodeidx << endl;


	// Find the best target split
	size_t targetidx = treeData_->get_target();
	vector<size_t> sampleics_left, sampleics_right;
	num_t splitvalue;
	set<num_t> values_left;	
	// treeData_->split_target(minSplit, sampleics, sampleics_left, sampleics_right);
	treeData_->split_target(targetidx,minSplit,sampleics,sampleics_left,sampleics_right,splitvalue,values_left);


	size_t nreal_tot;
	size_t n_left = sampleics_left.size();
	size_t n_right = sampleics_right.size();

	num_t bestrelativedecrease = 0;
	size_t bestsplitter_i = nFeatures;

	// TODO Instead of going through all features, sample features here
	// using sampling weights proportional to current importances
	for(size_t featureidx = 0; featureidx < nFeatures; ++featureidx)
	{

		if(featureidx == targetidx)
		{
			continue;
		}

		num_t impurity_tot,impurity_left,impurity_right;
		size_t nreal_left,nreal_right;
		treeData_->impurity(featureidx, sampleics, impurity_tot, nreal_tot);

		if(impurity_tot < datadefs::eps)
		{
			continue;
		}

		treeData_->impurity(featureidx,sampleics_left,impurity_left, nreal_left);
		treeData_->impurity(featureidx,sampleics_right,impurity_right, nreal_right);
		assert(sampleics_left.size() == n_left);
		assert(sampleics_right.size() == n_right);

		num_t relativedecrease = (impurity_tot-nreal_left*impurity_left/nreal_tot-nreal_right*impurity_right/nreal_tot)/impurity_tot;

		if(relativedecrease > bestrelativedecrease)
		{
			bestrelativedecrease = relativedecrease;
			bestsplitter_i = featureidx;
		}
	}

	if(bestsplitter_i == nFeatures)
	{
		cout << "No splitter found, quitting." << endl;
		// Leaf node
		GBT::SetNodePrediction(treeidx, nodeidx, sampleics);
		return;
	}

	size_t splitterfeatureidx = bestsplitter_i;

	cout << "Best splitter featureidx=" << splitterfeatureidx << " with relative decrease in impurity of " << bestrelativedecrease << endl;

	treeData_->remove_nans(splitterfeatureidx,sampleics,nreal_tot);
	cout << "Splitter feature has " << n_tot - nreal_tot << " missing values, which will be omitted in splitting" << endl;
	n_tot = nreal_tot;

	size_t nodeidx_left  = ++nNodes_[treeidx];
	size_t nodeidx_right = ++nNodes_[treeidx];


	treeData_->split_target(splitterfeatureidx,minSplit,sampleics,sampleics_left,sampleics_right,splitvalue,values_left);
	if(sampleics_left.size() == 0 || sampleics_right.size() == 0)
	{
		// Leaf node
		GBT::SetNodePrediction(treeidx, nodeidx, sampleics);
		return;
	}
	if(treeData_->isfeaturenum(splitterfeatureidx))
	{
		forest_[treeidx][nodeidx].set_splitter(splitterfeatureidx,splitvalue,forest_[treeidx][nodeidx_left],forest_[treeidx][nodeidx_right]);
	}
	else
	{
		forest_[treeidx][nodeidx].set_splitter(splitterfeatureidx,values_left,forest_[treeidx][nodeidx_left],forest_[treeidx][nodeidx_right]);
	}
	cout << forest_[treeidx][nodeidx].has_children() << endl;

	

	GBT::recursiveNodesplit(treeidx, nodeidx_left, sampleics_left);
	GBT::recursiveNodesplit(treeidx, nodeidx_right, sampleics_right);
}




void GBT::selectTarget(size_t targetidx)
{
	if(treeData_->get_target() != targetidx)
	{
		treeData_->select_target(targetidx);
	}
}


size_t GBT::getTarget()
{
	return(treeData_->get_target());
}
