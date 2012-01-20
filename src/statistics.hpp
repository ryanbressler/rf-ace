#ifndef STATISTICS_HPP
#define STATISTICS_HPP

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <cassert>

#include "datadefs.hpp"

using namespace std;
using datadefs::num_t;
using datadefs::NUM_NAN;

namespace statistics {
  
  class RF_statistics {
    
  public:

    RF_statistics();
    RF_statistics(vector<vector<num_t> > importanceMat, vector<vector<num_t> > contrastImportanceMat, vector<vector<size_t> > nodeMat, num_t executionTime);

    void printContrastImportance(ofstream& toFile);
    
    void print(ofstream& toFile);

  private:

    vector<vector<num_t> > importanceMat_;
    vector<vector<num_t> > contrastImportanceMat_;

    vector<vector<size_t> > nodeMat_;

    num_t executionTime_;

  };
}


#endif
