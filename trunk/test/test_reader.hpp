#ifndef TEST_READER_HPP
#define TEST_READER_HPP

#include "newtest.hpp"
#include "reader.hpp"
#include "datadefs.hpp"
#include "treedata.hpp"

using namespace std;
using datadefs::num_t;

void test_readAFM();

void test_reader() {

  newtest( "Reading delimited data from AFM file", &test_readAFM );

}

void test_readAFM() {

  Reader reader("test/data/3by8_mixed_NA_matrix.afm",'\t',"NA");

  newassert( reader.nLines() == 4 );

  size_t nSamples = reader.nLines() - 1;

  vector<Feature> features;

  // Removing top-left corner from table having column and row headers
  reader.nextLine().skipField();

  size_t nVars = 0;

  // Check that all variable names are valid
  for ( ; ! reader.endOfLine(); ++nVars ) {
    string varName; reader >> varName;
    if ( varName.substr(0,2) == "N:" ) {
      features.push_back( Feature(Feature::Type::NUM,varName,nSamples) );
    } else if ( varName.substr(0,2) == "C:" ) {
      features.push_back( Feature(Feature::Type::CAT,varName,nSamples) );
    } else if ( varName.substr(0,2) == "T:" ) {
      features.push_back( Feature(Feature::Type::TXT,varName,nSamples) );
    } else {
      newassert( false );
    }
    
  }
  
  newassert( nVars == 8 );
  newassert( features.size() == 8 );

  // We should have reached end of the first line
  newassert( reader.endOfLine() );
  
  // Get the next line and start reading...
  reader.nextLine();
  string field;
  reader >> field; newassert( field == "s0" );  
  reader >> field; newassert( field == "NA" );  
  reader >> field; newassert( field == "foo" ); 
  reader >> field; newassert( field == "2.2" ); 
  reader >> field; newassert( field == "3.3" ); 
  reader >> field; newassert( field == "4.4" ); 
  reader >> field; newassert( field == "5.5" ); 
  reader >> field; newassert( field == "6.6" ); 
  reader >> field; newassert( field == "Ah, be so good. Yes, no?" ); 

  // Make sure that we reached end of line again
  newassert( reader.endOfLine() ); 

  // Go to the start of file and get first line
  reader.rewind().nextLine();

  vector<string> sampleNames(nSamples);

  // Go through lines 2,3,...
  for ( size_t i = 0; i < nSamples; ++i ) {
    reader.nextLine();
    // Sample name is the first field of the line
    reader >> sampleNames[i];
    for ( size_t j = 0; j < nVars; ++j ) {
      if ( features[j].isNumerical() ) {
	num_t val; reader >> val;
	features[j].setNumSampleValue(i,val);
      } else if ( features[j].isCategorical() ) {
	string str; reader >> str;
	features[j].setCatSampleValue(i,str);
      } else if ( features[j].isTextual() ) {
	string str; reader >> str;
	features[j].setTxtSampleValue(i,str);
      }
    }
    // By now, we should have reached end of line
    newassert( reader.endOfLine() );     
  }

  // Did we recover the correct sample names from the file
  newassert( sampleNames[0] == "s0" ); 
  newassert( sampleNames[1] == "s1" ); 
  newassert( sampleNames[2] == "s2" ); 

  // Rewind again to the start, and start reading from line 2
  reader.rewind().nextLine().nextLine();

  // Variables for storing all data on line 2
  string s0;
  num_t  v1,v3,v4,v5,v6,v7;
  string v2,v8;
  
  // Read the 2nd line in one pass
  reader >> s0 >> v1 >> v2 >> v3 >> v4 >> v5 >> v6 >> v7 >> v8; 

  // Again, end of line should have been reached
  newassert( reader.endOfLine() );

  // Make sure the content of the 2nd line is as expected
  newassert( s0 == "s0" ); 
  newassert( datadefs::isNAN(v1) ); 
  newassert( v2 == "foo" ); 
  newassert( fabs( v3 - 2.2 ) < datadefs::EPS ); 
  newassert( fabs( v4 - 3.3 ) < datadefs::EPS ); 
  newassert( fabs( v5 - 4.4 ) < datadefs::EPS ); 
  newassert( fabs( v6 - 5.5 ) < datadefs::EPS ); 
  newassert( fabs( v7 - 6.6 ) < datadefs::EPS ); 
  newassert( v8 == "Ah, be so good. Yes, no?" ); 

  // Go back to the beginning
  reader.rewind();

  // While reading the whole file line by line till the end, we should 
  // not reach end of line nor end of file, since we have the last line
  // stored in the linefeed...
  for ( size_t i = 0; i < reader.nLines(); ++i ) {
    reader.nextLine();
    newassert( ! reader.endOfLine() );
    newassert( ! reader.endOfFile() );
  }

  // ... that means that we can then read the last line, field by field, 
  // into string variables
  for ( size_t i = 0; i < nVars + 1; ++i ) {
    newassert( ! reader.endOfLine() );
    string field; reader >> field;
  }

  // After we are done reading the last line, we should have reached end of line 
  // and end of file
  newassert( reader.endOfLine() );
  newassert( reader.endOfFile() );

}

#endif
