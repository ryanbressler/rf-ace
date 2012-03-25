#include "progress.hpp"


Progress::Progress(): 
  width_(3) { 
  cout << setw(width_) << "0" << "%" << flush; 
}

Progress::~Progress() { 
  reset(); 
}

void Progress::update(const num_t fraction) { 
  
  reset(); 
  
  cout << setw(width_) << static_cast<size_t>(fraction*100) << "%" << flush; 

}

void Progress::reset() { 

  for(size_t i = 0; i <= width_; ++i) { 
    cout << "\b"; 
  } 

}




