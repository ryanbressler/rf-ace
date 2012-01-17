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
    RF_statistics(vector<vector<num_t> > importanceMat, vector<vector<num_t> > contrastImportanceMat, vector<vector<size_t> > nodeMat);

    void printContrastImportance(ofstream& toFile) {
      
      size_t nFeatures = contrastImportanceMat_[0].size();
      size_t nPerms = contrastImportanceMat_.size();

      for ( size_t featureIdx = 0; featureIdx < nFeatures; ++featureIdx ) {
	
	vector<num_t> fSample( nPerms );

	for( size_t permIdx = 0; permIdx < nPerms; ++permIdx ) {
	  fSample[permIdx] = contrastImportanceMat_[permIdx][featureIdx];
	}

	size_t nReal = 0;
	num_t mu = 0.0;

	datadefs::mean(fSample,mu,nReal);

	if ( nReal == 0 ) {
	  mu = datadefs::NUM_NAN;
	}

	toFile << mu << endl;
        
      }


    }
    
    void print(ofstream& toFile) {
      
      assert( nodeMat_.size() > 0 );

      size_t nPerms = importanceMat_.size();
      size_t nTrees = nodeMat_[0].size();

      assert( nPerms == contrastImportanceMat_.size() );
      assert( nPerms == nodeMat_.size() );

      vector<num_t> importanceVec( nPerms );
      vector<num_t> contrastImportanceVec( nPerms );
      
      size_t nNodes = 0;

      size_t nReal;
      for ( size_t permIdx = 0; permIdx < nPerms; ++permIdx ) {
	datadefs::mean(importanceMat_[permIdx],importanceVec[permIdx],nReal);
	datadefs::mean(contrastImportanceMat_[permIdx],contrastImportanceVec[permIdx],nReal);
	for ( size_t treeIdx = 0; treeIdx < nodeMat_[permIdx].size(); ++treeIdx ) {
	  nNodes += nodeMat_[permIdx][treeIdx];
	}
      }

      num_t meanNodesPerTree = 1.0 * nNodes / ( nPerms * nTrees );

      num_t meanNodesPerSecond = 1.0 * nNodes / executionTime;

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
	     << "-- NUMBER OF TREES PER FOREST = " << nTrees << endl
	     << "--          NUMBER OF FORESTS = " << nPerms << endl
	     << "--            MEAN IMPORTANCE = " << meanImportance << endl
	     << "--             STD IMPORTANCE = " << stdImportance << endl
	     << "--   MEAN CONTRAST IMPORTANCE = " << meanContrastImportance << endl
	     << "--    STD CONTRAST IMPORTANCE = " << stdContrastImportance << endl
	     << "--        MEAN NODES PER TREE = " << meanNodesPerTree << endl
	     << "--      MEAN NODES PER SECOND = " << meanNodesPerSecond << endl;
      

    }

  private:

    vector<vector<num_t> > importanceMat_;
    vector<vector<num_t> > contrastImportanceMat_;

    vector<vector<size_t> > nodeMat_;

    num_t executionTime;

    num_t meanContrastNodesPerTree;
    num_t stdContrastNodesPerTree;

    num_t meanTimePerForest;
    num_t stdTimePerForest;
    
    
  };
}


#endif
