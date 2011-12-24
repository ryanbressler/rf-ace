#include "splitter.hpp"

/*
  Splitter::Splitter():
  splitterType_(NO_SPLITTER) {
  
  }
*/

Splitter::Splitter(const string& splitterName,
		   num_t splitLeftLeqValue):
  splitterType_(NUMERICAL_SPLITTER),
  splitterName_(splitterName),
  splitLeftLeqValue_(splitLeftLeqValue) {

  /* Empty constructor */

}

Splitter::Splitter(const string& splitterName,
		   const set<num_t>& splitLeftValues,
		   const set<num_t>& splitRightValues):
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

}

bool Splitter::splitsLeft(num_t testValue) {

  if ( splitterType_ == NO_SPLITTER ) {
    cerr << "Splitter::splitsLeft -- cannot split since the splitter is not set!" << endl;
    exit(1);
  }

  if ( splitterType_ == NUMERICAL_SPLITTER ) {
    
    if ( testValue < splitLeftLeqValue_ ) {
      return( true );
    } else {
      return( false );
    }

  } else if ( splitterType_ == CATEGORICAL_SPLITTER ) {
    
    if ( splitLeftValues_.find( testValue ) != splitLeftValues_.end() ) {
      return( true );
    } else {
      return( false );
    }

  } else {
    cerr << "Splitter::splitsLeft -- unknown splitter type!" << endl;
    exit(1);
  }

}

bool Splitter::splitsRight(num_t testValue) {

  if ( splitterType_ == NO_SPLITTER ) {
    cerr << "Splitter::splitsLeft -- cannot split since the splitter is not set!" << endl;
    exit(1);
  }

  if ( splitterType_ == NUMERICAL_SPLITTER ) {
    
    if ( splitLeftLeqValue_ <= testValue ) {
      return( true );
    } else {
      return( false );
    }
    
  } else if ( splitterType_ == CATEGORICAL_SPLITTER ) {
    
    if ( splitRightValues_.find( testValue ) != splitRightValues_.end() ) {
      return( true );
    } else {
      return( false );
    }
    
  } else {
    cerr << "Splitter::splitsLeft -- unknown splitter type!" << endl;
    exit(1);
  }
  
}

/*
  void Splitter::print(ofstream& toFile) {
  
  toFile << this->name();
  
  if ( splitterType_ == NO_SPLITTER ) {
  toFile << "Splitter: NOT SET" << endl;
  } else if ( splitterType_ == CATEGORICAL_SPLITTER ) {
  toFile << "Splitter: CATEGORICAL" << endl;
  toFile << "[";
  for ( set<num_t>::const_iterator it( splitLeftValues_.begin() ); it != splitLeftValues_.end(); ++it) {
  toFile << " " << *it;
  }
  toFile << " ] <==> [";
  for ( set<num_t>::const_iterator it( splitRightValues_.begin() ); it != splitRightValues_.end(); ++it) {
  toFile << " " << *it;
  }
  toFile << " ]" << endl;
  } else {
  toFile << "Splitter: NUMERICAL" << endl;
  toFile << " < " << splitLeftLeqValue_ << endl;  
  }
  
  }
*/
