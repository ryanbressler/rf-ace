#ifndef SPLITTER_HPP
#define SPLITTER_HPP

#include <set>
#include <string>

#include "datadefs.hpp"

using namespace std;
using datadefs::num_t;

class Splitter {

public:

  //Splitter();
  Splitter(const string& splitterName, const num_t splitLeftLeqValue);
  Splitter(const string& splitterName, const set<string>& splitLeftValues, const set<string>& splitRightValues);
  Splitter(const Splitter& splitter);

  ~Splitter();

  bool splitsLeft(const num_t testValue);
  bool splitsRight(const num_t testValue);

  bool splitsLeft(const string& testValue);
  bool splitsRight(const string& testValue);

  inline string name() { return(splitterName_); }

  void print();

#ifndef TEST__
private:
#endif

  enum SplitterType {NO_SPLITTER,NUMERICAL_SPLITTER,CATEGORICAL_SPLITTER};

  SplitterType splitterType_; 
  string splitterName_;

  num_t splitLeftLeqValue_;
  set<string> splitLeftValues_;
  set<string> splitRightValues_;

};


#endif
