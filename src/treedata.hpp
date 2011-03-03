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
  Treedata(string fname, bool is_featurerows); //PARTIAL IMPLEMENTATION
  ~Treedata();

  void sort_wrt_feature(size_t featureidx) {/*LACKS IMPLEMENTATION*/};
  
  void sort_wrt_target() {/*LACKS IMPLEMENTATION*/};
  
  void find_split(size_t featureidx, 
		  vector<size_t> sampleics, 
		  size_t& split_pos, 
		  num_t& impurity_left, 
		  num_t& impurity_right) {/*LACKS IMPLEMENTATION*/};
  
  void split_at_pos(size_t featureidx,
		    vector<size_t> sampleics,
		    num_t& impurity_left,
		    num_t& impurity_right) {/*LACKS IMPLEMENTATION*/};
  
private:

  void transpose();

  size_t targetidx_;

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
