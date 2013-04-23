#ifndef DISTRIBUTIONS_HPP
#define DISTRIBUTIONS_HPP

#include <cstdlib>
#include <vector>
#include <tr1/random>
#include <ctime>
#include "datadefs.hpp"

using namespace std;

namespace distributions {

  typedef tr1::mt19937 Engine;
  // typedef std::ranlux_base_01 Engine;
  
  inline unsigned int generateSeed() { return( clock() + time(0) ); }

  class Random {
  public:
    
    // Initialize the generator
    Random();
    Random(size_t seed);
    
    // Destructor
    ~Random();
    
    void seed(size_t seed);
    
    // Return random int
    size_t integer();

    // Generate and normalize random int
    datadefs::num_t uniform();
    
    size_t minIdx() { return( 0 ); }
    size_t maxIdx() { return( datadefs::MAX_IDX ); }
    
  private:
    
    Engine eng_;
    tr1::uniform_int<size_t> rand_;
    
  };
 
  class PMF {
  public:
    
    PMF(const vector<datadefs::num_t>& weights);
    ~PMF();
    
    size_t sample(Random* random) const;
    
#ifndef TEST__
  private:
#endif
    
    vector<datadefs::num_t> prob_;
    vector<size_t> alias_;
    
  };
  
}

#endif
