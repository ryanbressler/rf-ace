#include<iostream>
#include<cassert>
#include "node.hpp"

Node::Node():
  isSplitterNumerical_(true),
  hasChildren_(false),
  leftChild_(NULL),
  rightChild_(NULL)
{
}

Node::~Node()
{
  if(hasChildren_)
    {
      delete leftChild_;
      delete rightChild_;
    }
}

void Node::setSplitter(size_t splitter, set<num_t> classSet)
{
  assert(!hasChildren_);
  
  isSplitterNumerical_ = false;

  splitter_ = splitter;
  classSet_ = classSet;

  leftChild_ = new Node;
  rightChild_ = new Node;
  hasChildren_ = true;
}

void Node::setSplitter(size_t splitter, num_t threshold)
{
  assert(!hasChildren_);
  isSplitterNumerical_ = true;

  splitter_ = splitter;
  threshold_ = threshold;

  leftChild_ = new Node;
  rightChild_ = new Node;
  hasChildren_ = true;
}

/*
  int Node::getSplitter()
  {
  assert(hasChildren_);
  return(splitter_);
  }
*/

Node* Node::percolateData(num_t value)
{

  if(!hasChildren_) { return(this); }

  if(isSplitterNumerical_)
    {
      if(value <= threshold_) { return(leftChild_); } else { return(rightChild_); }
    }
  else
    {
      if(classSet_.find(value) != classSet_.end()) { return(leftChild_); } else { return(rightChild_); }
    }
}

Node* Node::leftChild()
{
  if(hasChildren_)
    {
      return(leftChild_);
    }
  else
    {
      return(NULL);
    }
}

Node* Node::rightChild()
{
  if(hasChildren_)
    {
      return(rightChild_);
    }
  else
    {
      return(NULL);
    }
}

size_t Node::nNodes()
{
  size_t n = 1;
  this->recursiveNDescendantNodes(n);
  return(n);
}

void Node::recursiveNDescendantNodes(size_t& n)
{
  if(!hasChildren_)
    {
      return;
    }
  else
    {
      n += 2;
      leftChild_->recursiveNDescendantNodes(n);
      rightChild_->recursiveNDescendantNodes(n);
    }

}

void Node::setPrediction(num_t value)
{
  prediction_ = value;
}

num_t Node::getPrediction()
{
  return(prediction_);
}

void Node::recursiveNodeSplit(Treedata* treeData,
			      const size_t targetIdx,
			      const vector<size_t>& sampleIcs,
			      const size_t maxNodesToStop,
			      const size_t minNodeSizeToStop,
			      const bool isRandomSplit,
			      const size_t nFeaturesInSample,
			      set<size_t>& featuresInTree,
			      size_t& nNodes)
{



  size_t nSamples = sampleIcs.size();

  if(nSamples < 2 * minNodeSizeToStop || nNodes + 1 > maxNodesToStop)
    {
      //cout << "Too few samples to start with, quitting" << endl;
      return;
    }


  vector<size_t> featureSample(nFeaturesInSample);

  const size_t nFeatures = treeData->nFeatures();

  if(isRandomSplit)
    {
      for(size_t i = 0; i < nFeaturesInSample; ++i)
        {
          treeData->getRandomIndex(2*nFeatures,featureSample[i]);
        }
    }
  else
    {
      for(size_t i = 0; i < nFeaturesInSample; ++i)
	{
	  featureSample[i] = i;
	}
    }

  vector<num_t> targetData;
  const bool isTargetNumerical = treeData->isFeatureNumerical(targetIdx);
  treeData->getFeatureData(targetIdx,sampleIcs,targetData);

  vector<size_t> sampleIcs_left,sampleIcs_right;
  num_t splitValue;
  set<num_t> values_left;

  if(isTargetNumerical)
    {
      Node::numericalFeatureSplit(targetData,isTargetNumerical,targetData,minNodeSizeToStop,sampleIcs_left,sampleIcs_right,splitValue);
    }
  else
    {
      Node::categoricalFeatureSplit(targetData,isTargetNumerical,targetData,sampleIcs_left,sampleIcs_right,values_left);
    }

  vector<num_t> fitnessVector(nFeaturesInSample, 0.0);
  num_t bestFitness = 0.0;
  size_t bestFeatureIdx = nFeaturesInSample;

  for(size_t i = 0; i < nFeaturesInSample; ++i)
    {

      vector<num_t> featureData;
      size_t featureIdx = featureSample[i];
      bool isFeatureNumerical = treeData->isFeatureNumerical(featureIdx);

      //Neither the real nor the contrast feature can appear in the tree as splitter
      if(featureIdx == targetIdx)
        {
          continue;
        }

      treeData->getFeatureData(featureIdx,sampleIcs,featureData);

      fitnessVector[i] = Node::splitFitness(featureData,isFeatureNumerical,minNodeSizeToStop,sampleIcs_left,sampleIcs_right);

    }

  for(size_t i = 0; i < nFeaturesInSample; ++i)
    {
      if(fitnessVector[i] > bestFitness)
        {
          bestFitness = fitnessVector[i];
          bestFeatureIdx = featureSample[i];
        }
    }


  if(bestFeatureIdx == nFeaturesInSample)
    {
      //cout << "No splitter found, quitting" << endl << endl;
      return;
    }

  vector<num_t> featureData;
  treeData->getFeatureData(bestFeatureIdx,sampleIcs,featureData);

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

  if(nRealSamples < 2*minNodeSizeToStop)
    {
      //cout << "Splitter has too few non-missing values, quitting" << endl;
      //This needs to be fixed such that one of the surrogates will determine the split instead
      return;
    }

  bool isBestFeatureNumerical = treeData->isFeatureNumerical(bestFeatureIdx);

  if(isBestFeatureNumerical)
    {
      Node::numericalFeatureSplit(targetData,isTargetNumerical,featureData,minNodeSizeToStop,sampleIcs_left,sampleIcs_right,splitValue);
    }
  else
    {
      Node::categoricalFeatureSplit(targetData,isTargetNumerical,featureData,sampleIcs_left,sampleIcs_right,values_left);
    }

  assert(sampleIcs_left.size() + sampleIcs_right.size() == nRealSamples);

  if(sampleIcs_left.size() < minNodeSizeToStop || sampleIcs_right.size() < minNodeSizeToStop)
    {
      return;
    }

  if(isBestFeatureNumerical)
    {
      //cout << "num splitter" << endl;
      Node::setSplitter(bestFeatureIdx,splitValue);
    }
  else
    {
      //cout << "cat splitter" << endl;
      Node::setSplitter(bestFeatureIdx,values_left);
    }

  featuresInTree.insert(bestFeatureIdx);
  nNodes += 2;

  for(size_t i = 0; i < sampleIcs_left.size(); ++i)
    {
      sampleIcs_left[i] = sampleIcs[sampleIcs_left[i]];
    }

  //cout << "split left..." << endl;
  leftChild_->recursiveNodeSplit(treeData,targetIdx,sampleIcs_left,maxNodesToStop,minNodeSizeToStop,isRandomSplit,nFeaturesInSample,featuresInTree,nNodes);

  for(size_t i = 0; i < sampleIcs_right.size(); ++i)
    {
      sampleIcs_right[i] = sampleIcs[sampleIcs_right[i]];
    }

  //cout << "split right..." << endl;
  rightChild_->recursiveNodeSplit(treeData,targetIdx,sampleIcs_right,maxNodesToStop,minNodeSizeToStop,isRandomSplit,nFeaturesInSample,featuresInTree,nNodes);
  
}

void Node::numericalFeatureSplit(vector<num_t>& tv,
				 const bool isTargetNumerical,
				 vector<num_t>& fv,
				 const size_t min_split,
				 vector<size_t>& sampleIcs_left,
				 vector<size_t>& sampleIcs_right,
				 num_t& splitValue)
{

  assert(tv.size() == fv.size());

  size_t n_tot = tv.size();
  size_t n_right = n_tot;
  size_t n_left = 0;

  sampleIcs_left.clear();
  sampleIcs_right.clear();

  //Check that there are enough samples to make the split in the first place
  assert(n_tot >= 2*min_split);

  //Make reference indices that define the sorting wrt. feature
  vector<size_t> refIcs;

  //Sort feature vector and collect reference indices
  datadefs::sortDataAndMakeRef(fv,sampleIcs_right);

  //Use the reference indices to sort sample indices
  //datadefs::sortFromRef<size_t>(sampleics,refIcs);
  datadefs::sortFromRef<num_t>(tv,sampleIcs_right);

  //Count how many real values the feature and target has
  size_t nreal_f,nreal_t;
  datadefs::countRealValues(fv,nreal_f);
  datadefs::countRealValues(tv,nreal_t);
  //cout << nreal_f << " " << nreal_t << " " << n_tot << endl;
  assert(nreal_t == n_tot && nreal_f == n_tot);

  int bestSplitIdx = -1;

  //If the target is numerical, we use the iterative squared error formula to update impurity scores while we traverse "right"
  if(isTargetNumerical)
    {
      num_t mu_right = 0.0;
      num_t se_right = 0.0;
      num_t mu_left = 0.0;
      num_t se_left = 0.0;
      num_t se_best = 0.0;
      size_t nreal_right = 0;

      datadefs::sqerr(tv,mu_right,se_right,nreal_right);
      assert(n_tot == nreal_right);
      se_best = se_right;

      size_t idx = 0;
      while(n_left < n_tot - min_split)
        {
	  datadefs::forward_backward_sqerr(tv[idx],n_left,mu_left,se_left,n_right,mu_right,se_right);
          if( se_left + se_right < se_best && n_left >= min_split)
            {
              bestSplitIdx = idx;
              se_best = se_left + se_right;
            }
          ++idx;
        }
    }
  else //Otherwise we use the iterative gini index formula to update impurity scores while we traverse "right"
    {
      map<num_t,size_t> freq_left;
      map<num_t,size_t> freq_right;
      size_t sf_left = 0;
      size_t sf_right = 0;
      size_t nreal_right = 0;

      datadefs::sqfreq(tv,freq_right,sf_right,nreal_right);
      num_t nsf_best = 1.0 * sf_right / nreal_right;
      assert(n_tot == nreal_right);

      size_t idx = 0;
      while(n_left < n_tot - min_split)
        {
	  datadefs::forward_backward_sqfreq(tv[idx],n_left,freq_left,sf_left,n_right,freq_right,sf_right);
          if(1.0 * n_right * sf_left + 1.0 * n_left * sf_right > n_left * n_right * nsf_best && n_left >= min_split)
            {
              bestSplitIdx = idx;
              nsf_best = 1.0 * sf_left / n_left + 1.0 * sf_right / n_right;
            }
          ++idx;
        }
    }

  splitValue = fv[bestSplitIdx];
  n_left = bestSplitIdx + 1;
  sampleIcs_left.resize(n_left);

  for(size_t i = 0; i < n_left; ++i)
    {
      sampleIcs_left[i] = sampleIcs_right[i];
    }
  sampleIcs_right.erase(sampleIcs_right.begin(),sampleIcs_right.begin() + n_left);

  //cout << sampleIcs_left.size() << " " << sampleIcs_right.size() << " == " << n_tot << endl;
  assert(sampleIcs_left.size() + sampleIcs_right.size() == n_tot);

  if(false)
    {
      cout << "Feature splits target [";
      for(size_t i = 0; i < sampleIcs_left.size(); ++i)
        {
          cout << " " << tv[sampleIcs_left[i]];
        }
      cout << " ] <==> [";
      for(size_t i = 0; i < sampleIcs_right.size(); ++i)
        {
          cout << " " << tv[sampleIcs_right[i]];
        }
      cout << " ]" << endl;
    }
}

void Node::categoricalFeatureSplit(vector<num_t>& tv,
				   const bool isTargetNumerical,
				   vector<num_t>& fv,
				   vector<size_t>& sampleIcs_left,
				   vector<size_t>& sampleIcs_right,
				   set<num_t>& categories_left)
{

  assert(tv.size() == fv.size());

  size_t n_tot = tv.size();
  size_t n_right = n_tot;
  size_t n_left = 0;

  sampleIcs_left.clear();
  sampleIcs_right.clear();
  categories_left.clear();

  //Check that sample size is positive
  assert(n_tot > 0);

  //Count how many real values the feature and target has (this is just to make sure there are no missing values thrown in)
  size_t nreal_f,nreal_t;
  datadefs::countRealValues(fv,nreal_f);
  datadefs::countRealValues(tv,nreal_t);
  assert(nreal_t == n_tot && nreal_f == n_tot);

  //Map all feature categories to the corresponding samples and represent it as map. The map is used to assign samples to left and right branches
  map<num_t,vector<size_t> > fmap_right;
  datadefs::map_data(fv,fmap_right,nreal_f);
  assert(n_tot == nreal_f);

  if(isTargetNumerical)
    {

      num_t mu_right;
      num_t mu_left = 0.0;
      num_t se_right;
      num_t se_left = 0.0;
      datadefs::sqerr(tv,mu_right,se_right,n_right);
      assert(n_tot == n_right);
      num_t se_best = se_right;

      while(fmap_right.size() > 0)
        {

          map<num_t,vector<size_t> >::iterator it_best(fmap_right.end());

          for(map<num_t,vector<size_t> >::iterator it(fmap_right.begin()); it != fmap_right.end(); ++it)
            {

              //Take samples from right and put them left
              for(size_t i = 0; i < it->second.size(); ++i)
                {
		  datadefs::forward_backward_sqerr(tv[it->second[i]],n_left,mu_left,se_left,n_right,mu_right,se_right);
                  //cout << n_left << "\t" << n_right << "\t" << se_left << "\t" << se_right << endl;
                }

              if(se_left+se_right < se_best)
                {
                  it_best = it;
                  se_best = se_left + se_right;
                }

              //Take samples from left back to right
              for(size_t i = 0; i < it->second.size(); ++i)
                {
                  //cout << tv[it->second[i]] << ": ";
		  datadefs::forward_backward_sqerr(tv[it->second[i]],n_right,mu_right,se_right,n_left,mu_left,se_left);
                  //cout << n_left << "\t" << n_right << "\t" << se_left << "\t" << se_right << endl;
                }
            }

          //If we found a categorical split that reduces impurity...
          if(it_best == fmap_right.end())
            {
              break;
            }


          //Take samples from right and put them left
          for(size_t i = 0; i < it_best->second.size(); ++i)
            {
              categories_left.insert(it_best->first);
              //cout << "removing index " << it_best->second[i] << " (value " << tv[it_best->second[i]] << ") from right: ";
	      datadefs::forward_backward_sqerr(tv[it_best->second[i]],n_left,mu_left,se_left,n_right,mu_right,se_right);
              sampleIcs_left.push_back(it_best->second[i]);
              //cout << n_left << "\t" << n_right << "\t" << se_left << "\t" << se_right << endl;
            }
          fmap_right.erase(it_best->first);

        }
    }
  else
    {
      map<num_t,size_t> freq_left,freq_right;
      size_t sf_left = 0;
      size_t sf_right;
      datadefs::sqfreq(tv,freq_right,sf_right,n_right);
      assert(n_tot == n_right);

      num_t nsf_best = 1.0 * sf_right / n_right;

      while(fmap_right.size() > 1)
        {

          map<num_t,vector<size_t> >::iterator it_best(fmap_right.end());

          for(map<num_t,vector<size_t> >::iterator it(fmap_right.begin()); it != fmap_right.end(); ++it)
            {

              size_t n_left_c = n_left;
              size_t n_right_c = n_right;
              map<num_t,size_t> freq_left_c = freq_left;
              map<num_t,size_t> freq_right_c = freq_right;
              //Take samples from right and put them left
              //cout << "Moving " << it->second.size() << " samples corresponding to feature category " << it->first << " from right to left" << endl;
              for(size_t i = 0; i < it->second.size(); ++i)
                {
                  //cout << " " << tv[it->second[i]];
		  datadefs::forward_backward_sqfreq(tv[it->second[i]],n_left,freq_left,sf_left,n_right,freq_right,sf_right);
                  //cout << n_left << "\t" << n_right << "\t" << sf_left << "\t" << sf_right << endl;
                }
              //cout << endl;

              if(1.0*n_right*sf_left + n_left*sf_right > n_left*n_right*nsf_best)
                {
                  it_best = it;
                  nsf_best = 1.0*sf_left/n_left + 1.0*sf_right/n_right;
                }

              //cout << "Moving " << it->second.size() << " samples corresponding to feature category " << it->first << " from left to right" << endl;
              //Take samples from left back to right
              for(size_t i = 0; i < it->second.size(); ++i)
                {
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
          if(it_best == fmap_right.end())
            {
              break;
            }


          //Take samples from right and put them left
          for(size_t i = 0; i < it_best->second.size(); ++i)
            {
              categories_left.insert(it_best->first);
              //cout << "removing index " << it_best->second[i] << " (value " << tv[it_best->second[i]] << ") from right: ";
	      datadefs::forward_backward_sqfreq(tv[it_best->second[i]],n_left,freq_left,sf_left,n_right,freq_right,sf_right);
              sampleIcs_left.push_back(it_best->second[i]);
              //cout << n_left << "\t" << n_right << "\t" << sf_left << "\t" << sf_right << endl;
            }
          fmap_right.erase(it_best->first);

        }
    }

  //Next we assign samples remaining on right
  sampleIcs_right.clear();
  for(map<num_t,vector<size_t> >::const_iterator it(fmap_right.begin()); it != fmap_right.end(); ++it)
    {
      for(size_t i = 0; i < it->second.size(); ++i)
        {
          sampleIcs_right.push_back(it->second[i]);
        }
    }

  if(false)
    {
      cout << "Feature splits target [";
      for(size_t i = 0; i < sampleIcs_left.size(); ++i)
        {
          cout << " " << tv[sampleIcs_left[i]];
        }
      cout << " ] <==> [";
      for(size_t i = 0; i < sampleIcs_right.size(); ++i)
        {
          cout << " " << tv[sampleIcs_right[i]];
        }
      cout << " ]" << endl;
    }
}

num_t Node::splitFitness(vector<num_t> const& data,
			 bool const& isFeatureNumerical,
			 size_t const& minSplit,
			 vector<size_t> const& sampleIcs_left,
			 vector<size_t> const& sampleIcs_right)
{

  assert(data.size() == sampleIcs_left.size() + sampleIcs_right.size());

  size_t n_left = 0;
  size_t n_right = 0;
  if(isFeatureNumerical)
    {
      num_t mu_left = 0.0;
      num_t se_left = 0.0;
      num_t mu_right = 0.0;
      num_t se_right = 0.0;

      for(size_t i = 0; i < sampleIcs_left.size(); ++i)
        {
	  datadefs::forward_sqerr(data[sampleIcs_left[i]],n_right,mu_right,se_right);
          //cout << "forward sqerr: " << featurematrix_[featureidx][sampleics_left[i]] << " " << n_right << " " << mu_right << " " << se_right << endl;
        }

      for(size_t i = 0; i < sampleIcs_right.size(); ++i)
        {
	  datadefs::forward_sqerr(data[sampleIcs_right[i]],n_right,mu_right,se_right);
          //cout << "forward sqerr: " << featurematrix_[featureidx][sampleics_right[i]] << " " << n_right << " " << mu_right << " " << se_right << endl;
        }

      if(n_right < 2*minSplit)
        {
          return(0.0);
        }

      num_t se_tot = se_right;

      for(size_t i = 0; i < sampleIcs_left.size(); ++i)
        {
	  datadefs::forward_backward_sqerr(data[sampleIcs_left[i]],n_left,mu_left,se_left,n_right,mu_right,se_right);
          //cout << "fw bw sqerr: " << featurematrix_[featureidx][sampleics_left[i]] << " " << n_left << " " << mu_left << " " << se_left << " " << n_right
	  //  << " " << mu_right << " " << se_right << endl;
        }

      if(n_left < minSplit || n_right < minSplit)
        {
          return(0.0);
        }

      return(( se_tot - se_left - se_right ) / se_tot);

    }
  else
    {
      map<num_t,size_t> freq_left,freq_right;
      size_t sf_left = 0;
      size_t sf_right = 0;

      for(size_t i = 0; i < sampleIcs_left.size(); ++i)
        {
	  datadefs::forward_sqfreq(data[sampleIcs_left[i]],n_right,freq_right,sf_right);
          //cout << "forward sqfreq: " << featurematrix_[featureidx][sampleics_left[i]] << " " << n_right << " " << sf_right << endl;
        }

      for(size_t i = 0; i < sampleIcs_right.size(); ++i)
        {
	  datadefs::forward_sqfreq(data[sampleIcs_right[i]],n_right,freq_right,sf_right);
          //cout << "forward sqfreq: " << featurematrix_[featureidx][sampleics_right[i]] << " " << n_right << " " << sf_right << endl;
        }

      if(n_right < 2*minSplit)
	{
	  return(0.0);
	}

      size_t n_tot = n_right;
      size_t sf_tot = sf_right;

      for(size_t i = 0; i < sampleIcs_left.size(); ++i)
        {
	  datadefs::forward_backward_sqfreq(data[sampleIcs_left[i]],n_left,freq_left,sf_left,n_right,freq_right,sf_right);
	  //cout << "fw bw sqfreq: " << featurematrix_[featureidx][sampleics_left[i]] << " " << n_left << " "<< sf_left << " " << n_right << " " << sf_right << endl;
        }

      if(n_left < minSplit || n_right < minSplit)
        {
          return(0.0);
        }

      //cout << n_left << " " << n_right << " " << sf_tot << " " << n_tot << " " << n_right << " " << sf_left << " " << n_tot*n_left*sf_right << " " << n_left*n_right << " " << pow(n_tot,2) - sf_tot << endl;

	//num_t fitness = (-1.0*(n_left*n_right*sf_tot) + n_tot*n_right*sf_left + n_tot*n_left*sf_right) / (n_left*n_right*(pow(n_tot,2) - sf_tot));
	//cout << "Fitness " << fitness << endl;

	return( ( -1.0 * n_left*n_right*sf_tot + n_tot*n_right*sf_left + n_tot*n_left*sf_right ) / ( n_left*n_right * (1.0*n_tot*n_tot - sf_tot) ) ) ;

    }

}
