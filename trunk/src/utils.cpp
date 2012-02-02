#include "utils.hpp"

vector<num_t> utils::removeNANs(const vector<num_t>& x) {

  vector<num_t> trimmed(x.size());

  size_t nRetained = 0;
  for(size_t i = 0; i < x.size(); ++i) {
    if( !datadefs::isNAN(x[i]) ) {
      trimmed[nRetained] = x[i];
      ++nRetained;
    }
  }
  trimmed.resize(nRetained);
  if ( nRetained == 0 ) {
    cout << "utils::removeNANs() -- data has 0 real values!" << endl;
  }

  return(trimmed);
}
