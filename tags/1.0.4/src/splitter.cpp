#include <sstream>
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

string Splitter::type() {

  switch ( splitterType_ ) {
  case NUMERICAL_SPLITTER:
    return("NUMERICAL");
  case CATEGORICAL_SPLITTER:
    return("CATEGORICAL");
  case NO_SPLITTER:
    return("");
  }
  
  return("");

}

string Splitter::leftSplitValues() {
  
  stringstream ss;

  //ss << "{";

  if ( splitterType_ == NUMERICAL_SPLITTER ) {
    ss << splitLeftLeqValue_;
  } else {

    //ss << "{";

    set<string>::const_iterator it(splitLeftValues_.begin());

    if ( it != splitLeftValues_.end() ) {
      ss << *it;
      ++it;
    }

    for ( ; it != splitLeftValues_.end(); ++it ) {
      ss << "," << *it;
    }

    //ss << "}";

  }
  
  //ss << "}";

  string str;

  ss >> str;

  return(str);
}

string Splitter::rightSplitValues() {

  stringstream ss;

  // ss << "{";

  if ( splitterType_ == NUMERICAL_SPLITTER ) {
    ss << splitLeftLeqValue_;
  } else {

    //ss << "{";

    set<string>::const_iterator it(splitRightValues_.begin());
    
    if ( it != splitRightValues_.end() ) {
      ss << *it;
      ++it;
    }
    
    for ( ; it != splitRightValues_.end(); ++it ) {
      ss << ":" << *it;
    }

    //ss << "}";
    
  }
  
  // ss << "}";
  
  string str;
  
  ss >> str;
  
  return(str);
}


string Splitter::print() {
  
  stringstream ss;

  ss << "SPLITTER=" << "\"" << this->name() << "\"" 
     << ",SPLITTERTYPE=" << this->type()
     << ",LVALUES=" << "\"" << this->leftSplitValues() << "\""
     << ",RVALUES=" << "\"" << this->rightSplitValues() << "\"";

  string str("");
  
  ss >> str;

  return(str);

}

