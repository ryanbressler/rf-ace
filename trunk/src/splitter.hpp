#ifndef SPLITTER_HPP
#define SPLITTER_HPP

#include <set>

#include "datadefs.hpp"

using namespace std;
using datadefs::num_t;

class Splitter {

public:

  Splitter(num_t splitLeftLeqValue);
  Splitter(const set<num_t>& splitLeftValues, const set<num_t>& splitRightValues);

  ~Splitter();

  bool splitsLeft(num_t testValue);
  bool splitsRight(num_t testValue);


private:

  enum SplitterType {NUMERICAL_SPLITTER,CATEGORICAL_SPLITTER};

  SplitterType splitterType_; 

  num_t splitLeftLeqValue_;
  set<num_t> splitLeftValues_;
  set<num_t> splitRightValues_;

};


#endif
