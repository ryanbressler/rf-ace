//treedata.hpp
//
//

#ifndef TREEDATA_HPP
#define TREEDATA_HPP

#include<cstdlib>
#include "datadefs.hpp"

using namespace std;
using datadefs::cat_t;
using datadefs::num_t;

class Treedata 
{
public:
  //Initializes the object and reads in a data matrix
  Treedata(string fname, bool is_featurerows);
  ~Treedata();

  size_t nfeatures();
  size_t nsamples();

  void select_target(size_t targetidx);

  void sort_all_wrt_feature(size_t featureidx);
  
  void sort_all_wrt_target();
  
  void find_split(int featureidx, 
		  vector<size_t> sampleics, 
		  size_t& split_pos, 
		  num_t& impurity_left, 
		  num_t& impurity_right) {/*LACKS IMPLEMENTATION*/};
  
  void split_at_pos(int featureidx,
		    vector<size_t> sampleics,
		    num_t& impurity_left,
		    num_t& impurity_right) {/*LACKS IMPLEMENTATION*/};
  
  void print();

private:

  void range(vector<size_t>& ics);
  template <typename T1,typename T2> void join_pairedv(vector<T1>& v1, vector<T2>& v2, vector<pair<T1,T2> >& p);
  template <typename T1,typename T2> void separate_pairedv(vector<pair<T1,T2> >& p, vector<T1>& v1, vector<T2>& v2);


  //Sorts a given input data vector of type T based on a given reference ordering of type vector<int>
  template <typename T> void sort_from_ref(vector<T>& in, vector<size_t> const& reference);

  //Sorts a given input data vector of type T1 based on a given reference ordering of type vector<pair<int,T2> >
  //template <typename T1, typename T2>
  //void sort_from_ref(vector<T1>& in,vector<pair<int,T2> > const& reference);

  void transpose();

  bool istarget_;
  size_t targetidx_;
  size_t internaltargetidx_;
  bool isnumtarget_;

  vector<vector<cat_t> > catmatrix_;
  vector<vector<num_t> > nummatrix_;

  size_t nsamples_;
  size_t nfeatures_;

  size_t ncatfeatures_;
  size_t nnumfeatures_;

  vector<string> catfeatureheaders_;
  vector<string> numfeatureheaders_;

  vector<size_t> catfeatureics_;
  vector<size_t> numfeatureics_;

  vector<string> sampleheaders_;

};

#endif
