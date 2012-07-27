#ifndef DISTRIBUTIONS_HPP
#define DISTRIBUTIONS_HPP

#include <cstdlib>
#include <vector>
#include <random>
#include <ctime>
#include "datadefs.hpp"

namespace distributions {

  //typedef std::mt19937 Engine;
  typedef std::ranlux_base_01 Engine;

  inline unsigned int generateSeed() { return( std::clock() + std::time(0) ); }

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

  private:

    Engine eng_;
    std::uniform_int<size_t> rand_;

  };

}

#endif
