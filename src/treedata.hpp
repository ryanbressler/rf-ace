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
  Treedata(string fname, bool is_featurerows);
  ~Treedata();
  
private:

  void transpose();

  size_t targetidx_;

  vector<vector<cat_t> > catmatrix_;
  vector<vector<num_t> > nummatrix_;

  size_t nsamples_;
  size_t nfeatures_;

  size_t ncatfeatures_;
  size_t nnumfeatures_;

  vector<size_t> catfeatureics_;
  vector<size_t> numfeatureics_;

  vector<string> featureheaders_;
  vector<string> sampleheaders_;

};

#endif
