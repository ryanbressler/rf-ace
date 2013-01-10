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

  newassert( reader.nRows() == 4 ); 
  newassert( reader.nCols() == 9 );

  size_t nVars = reader.nCols() - 1;
  size_t nSamples = reader.nRows() - 1;

  vector<Feature> features(nVars);

  // Removing top-left corner from table having column and row headers
  reader.nextLine().skipField();

  // Check that all variable names are valid
  for ( size_t i = 0; i < nVars; ++i ) {
    string varName; reader >> varName;
    if ( varName.substr(0,2) == "N:" ) {
      features[i] = Feature(Feature::Type::NUM,varName,nSamples);
    } else if ( varName.substr(0,2) == "C:" ) {
      features[i] = Feature(Feature::Type::CAT,varName,nSamples);
    } else if ( varName.substr(0,2) == "T:" ) {
      features[i] = Feature(Feature::Type::TXT,varName,nSamples);
    } else {
      newassert( false );
    }
    
  }
  
  newassert( reader.endOfLine() );
  
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

  newassert( reader.endOfLine() ); 

  reader.rewind().nextLine();

  vector<string> sampleNames(nSamples);

  for ( size_t i = 0; i < nSamples; ++i ) {
    reader.nextLine();
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
    newassert( reader.endOfLine() );     
  }

  newassert( sampleNames[0] == "s0" ); 
  newassert( sampleNames[1] == "s1" ); 
  newassert( sampleNames[2] == "s2" ); 

  reader.rewind().nextLine().nextLine();

  string s0;
  num_t  v1,v3,v4,v5,v6,v7;
  string v2,v8;
  
  reader >> s0 >> v1 >> v2 >> v3 >> v4 >> v5 >> v6 >> v7 >> v8; 

  newassert( s0 == "s0" ); 
  newassert( datadefs::isNAN(v1) ); 
  newassert( v2 == "foo" ); 
  newassert( fabs( v3 - 2.2 ) < datadefs::EPS ); 
  newassert( fabs( v4 - 3.3 ) < datadefs::EPS ); 
  newassert( fabs( v5 - 4.4 ) < datadefs::EPS ); 
  newassert( fabs( v6 - 5.5 ) < datadefs::EPS ); 
  newassert( fabs( v7 - 6.6 ) < datadefs::EPS ); 
  newassert( v8 == "Ah, be so good. Yes, no?" ); 

}

#endif
