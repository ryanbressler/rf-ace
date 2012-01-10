#include <cstdlib>
#include <iostream>
#include <vector>
#include <cassert>

#include "datadefs.hpp"

using namespace std;
using datadefs::num_t;
using datadefs::NUM_NAN;

namespace statistics {
  
  struct RF_statistics {
    
    vector<vector<num_t> > importanceMat;
    vector<vector<num_t> > contrastImportanceMat;
    
    num_t meanNodesPerTree;
    num_t stdNodesPerTree;
    
    num_t meanNodesCreatedPerSecond;
    
    num_t meanContrastNodesPerTree;
    num_t stdContrastNodesPerTree;
    
    num_t meanTimePerForest;
    num_t stdTimePerForest;
    
    RF_statistics():
      importanceMat(0),
      contrastImportanceMat(0),
      
      meanNodesPerTree(NUM_NAN),
      stdNodesPerTree(NUM_NAN),
      
      meanNodesCreatedPerSecond(NUM_NAN),
      
      meanContrastNodesPerTree(NUM_NAN),
      stdContrastNodesPerTree(NUM_NAN),
      
      meanTimePerForest(NUM_NAN),
      stdTimePerForest(NUM_NAN) {}
    
    void print(ofstream& toFile) {
      
      size_t nPerms = importanceMat.size();
      
      assert( nPerms == contrastImportanceMat.size() );

      vector<num_t> importanceVec( nPerms );
      vector<num_t> contrastImportanceVec( nPerms );

      size_t nReal;
      for ( size_t permIdx = 0; permIdx < nPerms; ++permIdx ) {
	datadefs::mean(importanceMat[permIdx],importanceVec[permIdx],nReal);
	datadefs::mean(contrastImportanceMat[permIdx],contrastImportanceVec[permIdx],nReal);
      }

      num_t meanImportance;
      num_t meanContrastImportance;

      num_t stdImportance;
      num_t stdContrastImportance;

      datadefs::sqerr(importanceVec,meanImportance,stdImportance,nReal);
      stdImportance = sqrtf(stdImportance) / nReal;

      datadefs::sqerr(contrastImportanceVec,meanContrastImportance,stdContrastImportance,nReal);
      stdContrastImportance = sqrtf(stdContrastImportance) / nReal;


      toFile << "Random Forest statistics" << endl
	     << "------------------------" << endl
	     << "-- MEAN          IMPORTANCE = " << meanImportance << endl
	     << "-- STD           IMPORTANCE = " << stdImportance << endl
	     << "-- MEAN CONTRAST IMPORTANCE = " << meanContrastImportance << endl
	     << "-- STD  CONTRAST IMPORTANCE = " << stdContrastImportance << endl;
    }
    
    
  };
}
