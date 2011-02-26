#ifndef NODE_HPP
#define NODE_HPP

#include<cstdlib>
#include<vector>
#include<set>
#include "datadefs.hpp"

using namespace std;
using namespace datadefs;

class Node {
public:
  Node(int nsamples, bool isregr);
  ~Node();

  void set_splitter(int splitter, set<cat_t> classet, int leftchild, int rightchild);
  void set_splitter(int splitter, num_t threshold, int leftchild, int rightchild);
  int get_splitter();
  
  int percolate(cat_t value);
  int percolate(num_t value);

  //Helper functions
  void print();
  void print_compact();

private:
  bool isregr_;

  int splitter_;
  num_t threshold_;
  set<cat_t> classet_;

  vector<cat_t> cat_trainsamples_;
  vector<cat_t> cat_testsamples_;

  vector<num_t> num_trainsamples_;
  vector<num_t> num_testsamples_;

  size_t ntrainsamples_;
  size_t ntestsamples_;
  
  bool haschildren_;
  int leftchild_;
  int rightchild_;
};

#endif
