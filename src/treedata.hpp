//treedata.hpp
//
//

#ifndef TREEDATA_HPP
#define TREEDATA_HPP

#include<cstdlib>
#include<map>
#include "datadefs.hpp"

using namespace std;
//using datadefs::cat_t;
using datadefs::num_t;

class Treedata 
{
public:
  //Initializes the object and reads in a data matrix
  Treedata(string fname, bool is_featurerows);
  ~Treedata();

  //Returns the number of features
  size_t nfeatures();

  //Returns the number of samples
  size_t nsamples();

  //Permutes integers in range 0,1,...,(ics.size()-1). 
  //NOTE: original contents in ics will be replaced.  
  void permute(vector<size_t>& ics);

  //Permutes data in x.
  void permute(vector<num_t>& x);

  //Generates a bootstrap sample. Samples not in the bootstrap sample will be stored in oob_ics, 
  //and the number of oob samples is stored in noob. 
  void bootstrap(vector<size_t>& ics, vector<size_t>& oob_ics, size_t& noob);

  //Samples features
  void sample_features(size_t nfeatures, vector<size_t>& fsample);

  //Selects target feature. 
  //NOTE: data will be sorted with respect to the target.
  void select_target(size_t targetidx);

  //Sorts data with respect to a given feature.
  void sort_all_wrt_feature(size_t featureidx);
  
  //Sorts data with respect to target.
  void sort_all_wrt_target();

  //Given feature, finds and returns the optimal split point wrt. all samples.
  //**IMPLEMENTATION NOT READY**
  //void find_split(size_t featureidx,
  //               size_t& split_pos,
  //               num_t& impurity_left,
  //               num_t& impurity_right);
  
  //Given feature, finds and returns the optimal split point wrt. all sampleics.
  //**IMPLEMENTATION NOT READY**
  void find_split(size_t featureidx, 
		  vector<size_t>& sampleics,
		  vector<size_t>& sampleics_left,
		  vector<size_t>& sampleics_right,
		  size_t& n_left,
		  size_t& n_right,
		  num_t& impurity,
		  num_t& impurity_left, 
		  num_t& impurity_right);
  
  //Performs split at given split point.
  //**IMPLEMENTATION NOT READY**
  void split_at_pos(size_t featureidx,
		    vector<size_t>& sampleics,
		    num_t& impurity_left,
		    num_t& impurity_right);
  
  //Prints contents in Treedata
  void print();

  void range(vector<size_t>& ics);

private:
  
  template <typename T1,typename T2> void make_pairedv(vector<T1> const& v1, 
						       vector<T2> const& v2, 
						       vector<pair<T1,T2> >& p);

  template <typename T1,typename T2> void separate_pairedv(vector<pair<T1,T2> > const& p, 
							   vector<T1>& v1, 
							   vector<T2>& v2);

  template <typename T> void sort_and_make_ref(vector<T>& v, vector<size_t>& ref_ics);

  //Sorts a given input data vector of type T based on a given reference ordering of type vector<int>
  template <typename T> void sort_from_ref(vector<T>& in, vector<size_t> const& ref_ics);

  size_t targetidx_;

  vector<vector<num_t> > featurematrix_;
  vector<bool> isfeaturenum_;

  size_t nsamples_;
  size_t nfeatures_;

  vector<size_t> sampleics_; //Indices 0,1,...,(nsamples_-1)
  vector<size_t> featureics_; //Indices 0,1,...,(nfeatures_-1)

  vector<string> featureheaders_;
  vector<string> sampleheaders_;

};

#endif
