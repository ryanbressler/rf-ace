#ifndef TEST_READER_HPP
#define TEST_READER_HPP

#include <regex>

#include "newtest.hpp"
#include "reader.hpp"
#include "datadefs.hpp"
#include "treedata.hpp"

using namespace std;
using datadefs::num_t;

extern size_t N_SUCCESS;
extern size_t N_FAIL;

void test_readAFM();

void test_reader() {

  newtest( "Reading delimited data from AFM file", &test_readAFM );

}

void test_readAFM() {

  Reader reader("test/data/3by8_mixed_NA_matrix.afm",'\t',"NA");

  assert( reader.nRows() == 4 ); N_SUCCESS++;
  assert( reader.nCols() == 9 ); N_SUCCESS++;

  /*
            N:var0  C:var1  N:var2  N:var3  N:var4  N:var5  N:var6  T:var7
    s0      NA      foo     2.2     3.3     4.4     5.5     6.6     Ah, be so good. Yes, no?
    s1      0.00    NA      2.22    3.33    4.44    5.55    NA      NA
    s2      0.000   bar     2.222   3.333   4.444   5.555   6.666   Some more text, but not much.
  */

  size_t nVars = reader.nCols() - 1;
  size_t nSamples = reader.nRows() - 1;

  vector<Feature> features(nVars);

  // Removing top-left corner from table having column and row headers
  reader.nextLine().skipField();

  // Regular expressions for valid variable names
  // NOTE: these will be moved to datadefs
  regex numVarRE("(N:)(.*)");
  regex catVarRE("(C:)(.*)");
  regex txtVarRE("(T:)(.*)");
  
  // Check that all variable names are valid
  for ( size_t i = 0; i < nVars; ++i ) {
    string varName; reader >> varName;
    if ( regex_match(varName, numVarRE) ) {
      features[i] = Feature(Feature::Type::NUM,varName,nSamples);
    } else if ( regex_match(varName, catVarRE) ) {
      features[i] = Feature(Feature::Type::CAT,varName,nSamples);
    } else if ( regex_match(varName, txtVarRE) ) {
      features[i] = Feature(Feature::Type::TXT,varName,nSamples);
    } else {
      N_FAIL++;
    }
    
  }
  
  assert( reader.endOfLine() ); N_SUCCESS++;
  
  reader.nextLine();

  string field;
  
  reader >> field; assert( field == "s0" );  N_SUCCESS++;
  reader >> field; assert( field == "NA" );  N_SUCCESS++;
  reader >> field; assert( field == "foo" ); N_SUCCESS++;
  reader >> field; assert( field == "2.2" ); N_SUCCESS++;
  reader >> field; assert( field == "3.3" ); N_SUCCESS++;
  reader >> field; assert( field == "4.4" ); N_SUCCESS++;
  reader >> field; assert( field == "5.5" ); N_SUCCESS++;
  reader >> field; assert( field == "6.6" ); N_SUCCESS++;
  reader >> field; assert( field == "Ah, be so good. Yes, no?" ); N_SUCCESS++;

  assert( reader.endOfLine() ); N_SUCCESS++;

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
    assert( reader.endOfLine() ); N_SUCCESS++;    
  }

  assert( sampleNames[0] == "s0" ); N_SUCCESS++;
  assert( sampleNames[1] == "s1" ); N_SUCCESS++;
  assert( sampleNames[2] == "s2" ); N_SUCCESS++;

  reader.rewind().nextLine().nextLine();

  string s0;
  num_t  v1,v3,v4,v5,v6,v7;
  string v2,v8;
  
  reader >> s0 >> v1 >> v2 >> v3 >> v4 >> v5 >> v6 >> v7 >> v8; 

  assert( s0 == "s0" ); N_SUCCESS++;
  assert( datadefs::isNAN(v1) ); N_SUCCESS++;
  assert( v2 == "foo" ); N_SUCCESS++;
  assert( fabs( v3 - 2.2 ) < datadefs::EPS ); N_SUCCESS++;
  assert( fabs( v4 - 3.3 ) < datadefs::EPS ); N_SUCCESS++;
  assert( fabs( v5 - 4.4 ) < datadefs::EPS ); N_SUCCESS++;
  assert( fabs( v6 - 5.5 ) < datadefs::EPS ); N_SUCCESS++;
  assert( fabs( v7 - 6.6 ) < datadefs::EPS ); N_SUCCESS++;
  assert( v8 == "Ah, be so good. Yes, no?" ); N_SUCCESS++;

}

#endif
