#ifndef PARTITIONSEQUENCE_HPP
#define PARTITIONSEQUENCE_HPP

#include<cstdlib>
#include<vector>

using namespace std;

class PartitionSequence {

public:

  typedef size_t graycode_t;

  PartitionSequence(const size_t nMaxLength);
  ~PartitionSequence();

  bool isAdded(const size_t pos);
  size_t at(const size_t pos);

private:

  vector<graycode_t> graycode_;
  vector<size_t> bitSequence_;
  vector<bool> addBit_;

};

#endif
