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

  class RandInt {
  public:
    
    // Initialize the generator
    RandInt();
    RandInt(size_t seed);
    
    // Destructor
    ~RandInt();
    
    void seed(size_t seed);
    
    // Return random int
    size_t operator()() {
      return( rand_(eng_) );
    }
    
    // Generate and normalize random int
    datadefs::num_t uniform();
    
    size_t minIdx() { return( 0 ); }
    size_t maxIdx() { return( datadefs::MAX_IDX ); }
    
  private:
    
    Engine eng_;
    tr1::uniform_int<size_t> rand_;
    
  };
  
  typedef map<datadefs::num_t,size_t> icdf_t;

  //icdf_t 

}

#endif
