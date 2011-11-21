#include "splitter.hpp"

Splitter::Splitter(num_t splitLeftLeqValue):
  splitterType_(NUMERICAL_SPLITTER),
  splitLeftLeqValue_(splitLeftLeqValue) {

  /* Empty constructor */

}

Splitter::Splitter(const set<num_t>& splitLeftValues,
		   const set<num_t>& splitRightValues):
  splitterType_(CATEGORICAL_SPLITTER),
  splitLeftValues_(splitLeftValues),
  splitRightValues_(splitRightValues) {

  /* Empty constructor */

}


Splitter::~Splitter() {

}

bool Splitter::splitsLeft(num_t testValue) {

  if ( splitterType_ == NUMERICAL_SPLITTER ) {
    
    if ( testValue < splitLeftLeqValue_ ) {
      return( true );
    } else {
      return( false );
    }

  } else {
    
    if ( splitLeftValues_.find( testValue ) != splitLeftValues_.end() ) {
      return( true );
    } else {
      return( false );
    }

  }

}

bool Splitter::splitsRight(num_t testValue) {

  if ( splitterType_ == NUMERICAL_SPLITTER ) {
    
    if ( splitLeftLeqValue_ <= testValue ) {
      return( true );
    } else {
      return( false );
    }
    
  } else {
    
    if ( splitRightValues_.find( testValue ) != splitRightValues_.end() ) {
      return( true );
    } else {
      return( false );
    }
    
  }
  
}

