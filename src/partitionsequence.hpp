#ifndef PARTITIONSEQUENCE_HPP
#define PARTITIONSEQUENCE_HPP

#include<cstdlib>
#include<vector>

using namespace std;

class PartitionSequence {

public:

  PartitionSequence();
  ~PartitionSequence();

  bool isAdded(const size_t pos);
  size_t at(const size_t pos);

private:

  size_t grayCode(const size_t pos);

  //typedef long int graycode_t;
    
  //unsigned bit_ : 1 << 63;
  //unsigned isAdded_ : 1 << 63;
    
  //vector<graycode_t> graycode_;
  // vector<size_t> bitSequence_;
  // vector<bool> addBit_;

};

#endif
