#include "statistics.hpp"

statistics::RF_statistics::RF_statistics() {

}

statistics::RF_statistics::RF_statistics(vector<vector<num_t> > importanceMat, 
					 vector<vector<num_t> > contrastImportanceMat,
					 vector<vector<size_t> > nodeMat,
					 num_t executionTime):

    importanceMat_(importanceMat),
    contrastImportanceMat_(contrastImportanceMat),
    nodeMat_(nodeMat),
    executionTime_(executionTime) {
  
  

}

void statistics::RF_statistics::printContrastImportance(ofstream& toFile) {

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

void statistics::RF_statistics::print(ofstream& toFile) {

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

  num_t meanNodesPerSecond = 1.0 * nNodes / executionTime_;

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

