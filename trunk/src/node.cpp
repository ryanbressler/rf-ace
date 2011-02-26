#include<cassert>
#include "node.hpp"

Node::Node(int nsamples, bool isregr):
  isregr_(isregr),
  ntrainsamples_(0),
  ntestsamples_(0),
  haschildren_(false)
{
  if(isregr_)
    {
      vector<num_t> num_trainsamples_(nsamples);
      vector<num_t> num_testsamples_(nsamples);
    }
  else
    {
      vector<cat_t> cat_trainsamples_(nsamples);
      vector<cat_t> cat_testsamples_(nsamples);
    }
}

Node::~Node()
{

}

void Node::set_splitter(int splitter, set<cat_t> classet, int leftchild, int rightchild)
{
  assert(!isregr_);
  assert(!haschildren_);

  splitter_ = splitter;
  classet_ = classet;

  leftchild_ = leftchild;
  rightchild_ = rightchild;
  haschildren_ = true;
}

void Node::set_splitter(int splitter, num_t threshold, int leftchild, int rightchild)
{
  assert(isregr_);
  assert(!haschildren_);

  splitter_ = splitter;
  threshold_ = threshold;

  leftchild_ = leftchild;
  rightchild_ = rightchild;
  haschildren_ = true;
}

int Node::get_splitter()
{
  assert(haschildren_);

  return(splitter_);
}

int Node::percolate(cat_t value)
{
  if(classet_.find(value) != classet_.end()) {return(leftchild_);} else {return(rightchild_);}
}

int Node::percolate(num_t value)
{
  if(value <= threshold_) {return(leftchild_);} else {return(rightchild_);}
}

