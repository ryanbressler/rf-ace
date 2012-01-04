#include "splitter.hpp"

Splitter::Splitter(const string& splitterName,
		   const num_t splitLeftLeqValue):
  splitterType_(NUMERICAL_SPLITTER),
  splitterName_(splitterName),
  splitLeftLeqValue_(splitLeftLeqValue) {

  /* Empty constructor */

}

Splitter::Splitter(const string& splitterName,
		   const set<string>& splitLeftValues,
		   const set<string>& splitRightValues):
  splitterType_(CATEGORICAL_SPLITTER),
  splitterName_(splitterName),
  splitLeftValues_(splitLeftValues),
  splitRightValues_(splitRightValues) {

  /* Empty constructor */

}


Splitter::Splitter(const Splitter& splitter) {
  
  splitterType_ = splitter.splitterType_;
  splitterName_ = splitter.splitterName_;
  splitLeftLeqValue_ = splitter.splitLeftLeqValue_;
  splitLeftValues_ = splitter.splitLeftValues_;
  splitRightValues_ = splitter.splitRightValues_;
  
}


Splitter::~Splitter() {
  /* Empty destructor */
}

bool Splitter::splitsLeft(const num_t testValue) {
  
  if ( splitterType_ == CATEGORICAL_SPLITTER ) {
    cerr << "Splitter::splitsLeft() -- tried to split with NUMERICAL data although splitter is CATEGORICAL" << endl;
    exit(1);
  }
  
  if ( testValue < splitLeftLeqValue_ ) {
    return( true );
  } else {
    return( false );
  }
  
}

bool Splitter::splitsLeft(const string& testValue) {
  
  if ( splitterType_ == NUMERICAL_SPLITTER ) {
    cerr << "Splitter::splitsLeft() -- tried to split with CATEGORICAL data although splitter is NUMERICAL" << endl;
    exit(1);
  }
  
  if ( splitLeftValues_.find( testValue ) != splitLeftValues_.end() ) {
    return( true );
  } else {
    return( false );
  }
  
}

bool Splitter::splitsRight(const num_t testValue) {

  if ( splitterType_ == CATEGORICAL_SPLITTER ) {
    cerr << "Splitter::splitsLeft() -- tried to split with NUMERICAL data although splitter is CATEGORICAL" << endl;
    exit(1);
  }

  if ( splitLeftLeqValue_ <= testValue) {
    return( true );
  } else {
    return( false );
  }

}

bool Splitter::splitsRight(const string& testValue) {

  if ( splitterType_ == NUMERICAL_SPLITTER ) {
    cerr << "Splitter::splitsRight() -- tried to split with CATEGORICAL data although splitter is NUMERICAL" << endl;
    exit(1);
  }

  if ( splitRightValues_.find( testValue ) != splitRightValues_.end() ) {
    return( true );
  } else {
    return( false );
  }

}

bool Splitter::splitsLeft(Treedata* treeData, const size_t sampleIdx) {
  
  size_t featureIdx = treeData->getFeatureIdx(splitterName_);

  if ( treeData->isFeatureNumerical(featureIdx) ) {
    return( this->splitsLeft(treeData->getFeatureData(featureIdx,sampleIdx)));
  } else {
    return( this->splitsLeft(treeData->getRawFeatureData(featureIdx,sampleIdx)));
  }

}

bool Splitter::splitsRight(Treedata* treeData, const size_t sampleIdx) {

  size_t featureIdx = treeData->getFeatureIdx(splitterName_);

  if ( treeData->isFeatureNumerical(featureIdx) ) {
    return( this->splitsRight(treeData->getFeatureData(featureIdx,sampleIdx)));
  } else {
    return( this->splitsRight(treeData->getRawFeatureData(featureIdx,sampleIdx)));
  }

}
