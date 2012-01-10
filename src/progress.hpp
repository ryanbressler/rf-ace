#include <cstdlib>
#include <iostream>
#include <iomanip>

#include "datadefs.hpp"

using namespace std;
using datadefs::num_t;

class Progress {
public:
  Progress();
  ~Progress();

  void update(const num_t fraction);

private:

  void reset();

  size_t width_;

};


