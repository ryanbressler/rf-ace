#include "statistics.hpp"

statistics::RF_statistics::RF_statistics() {

}

statistics::RF_statistics::RF_statistics(vector<vector<num_t> > importanceMat, 
					 vector<vector<num_t> > contrastImportanceMat,
					 vector<vector<size_t> > nodeMat):

    importanceMat_(importanceMat),
    contrastImportanceMat_(contrastImportanceMat),
    nodeMat_(nodeMat) {
  
}
