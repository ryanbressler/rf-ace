#ifndef UTILSTEST_HPP
#define UTILSTEST_HPP

#include <cppunit/extensions/HelperMacros.h>
#include <vector>
#include <string>
#include <algorithm>
#include "datadefs.hpp"
#include "utils.hpp"
#include "distributions.hpp"

class UtilsTest : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE( UtilsTest );
  CPPUNIT_TEST( test_removeNANs );
  CPPUNIT_TEST( test_parse );
  CPPUNIT_TEST( test_str2 );
  CPPUNIT_TEST( test_write );
  CPPUNIT_TEST( test_range );
  CPPUNIT_TEST( test_trim );
  CPPUNIT_TEST( test_chomp );
  CPPUNIT_TEST( test_split );
  CPPUNIT_TEST( test_permute );
  CPPUNIT_TEST( test_splitRange );

  // NEW
  CPPUNIT_TEST( test_strv2catv );
  CPPUNIT_TEST( test_strv2numv );
  CPPUNIT_TEST( test_sortDataAndMakeRef );
  CPPUNIT_TEST( test_sortFromRef );

  CPPUNIT_TEST_SUITE_END();
  
public:
  void setUp();
  void tearDown();
  
  void test_removeNANs();
  void test_parse();
  void test_str2();
  void test_write();
  void test_range();
  void test_trim();
  void test_chomp();
  void test_split();
  void test_permute();
  void test_splitRange();

  // NEW
  void test_strv2catv();
  void test_strv2numv();
  void test_sortDataAndMakeRef();
  void test_sortFromRef();
  
  
};

void UtilsTest::setUp() {}
void UtilsTest::tearDown() {}

void UtilsTest::test_removeNANs() {
  
}

void UtilsTest::test_parse() {

  string s1("KEY1=val1,KEY2='val2',key3='val3,which=continues\"here'");

  map<string,string> m1 = utils::parse(s1,',','=','\'');

  CPPUNIT_ASSERT( m1["KEY1"] == "val1" );
  CPPUNIT_ASSERT( m1["KEY2"] == "val2" );
  CPPUNIT_ASSERT( m1["key3"] == "val3,which=continues\"here");

  string s2("KEY1=\"\",KEY2=,KEY3=\"=\"");

  map<string,string> m2 = utils::parse(s2,',','=','"');

  CPPUNIT_ASSERT( m2["KEY1"] == "" );
  CPPUNIT_ASSERT( m2["KEY2"] == "" );
  CPPUNIT_ASSERT( m2["KEY3"] == "=" );
  
}

void UtilsTest::test_str2() {

  string a("0.0");
  string b("1.0");
  string c("-1.0");
  string d("-1.0e10");

  CPPUNIT_ASSERT(utils::str2<num_t>(a) == 0.0);
  CPPUNIT_ASSERT(utils::str2<num_t>(b) == 1.0);
  CPPUNIT_ASSERT(utils::str2<num_t>(c) == -1.0);
  CPPUNIT_ASSERT(utils::str2<num_t>(d) == -1.0e10);

}

void UtilsTest::test_write() {

  vector<string> foo;
  foo.push_back("a");
  foo.push_back("b");

  stringstream is;
  string result;

  utils::write(is,foo.begin(),foo.end(),',');
  is >> result;
  is.clear();

  CPPUNIT_ASSERT( result == "a,b" );

  utils::write(is,foo.begin(),foo.end(),'a');
  is >> result;
  is.clear();

  CPPUNIT_ASSERT( result == "aab" );

  foo.pop_back();

  utils::write(is,foo.begin(),foo.end(),',');
  is >> result;
  is.clear();

  CPPUNIT_ASSERT( result == "a" );

  foo.pop_back();

  // Nothing should be read into "is", which is why "result" stays unchanged
  utils::write(is,foo.begin(),foo.end(),',');
  is >> result;
  is.clear();

  CPPUNIT_ASSERT( result == "a" );

}

void UtilsTest::test_range() {

  size_t n = 50;

  vector<size_t> ics = utils::range(n);

  for (size_t i = 0; i < n; ++i) {
    CPPUNIT_ASSERT(ics[i] == i);
  }
  
}

void UtilsTest::test_trim() {
  
  string wh1 = " \t";
  string wh2 = " ";

  string str = " \t  a \tb\t ";

  CPPUNIT_ASSERT( utils::trim(str,wh1) == "a \tb" );
  CPPUNIT_ASSERT( utils::trim(str,wh2) == "\t  a \tb\t" );

  str = "";
  
  CPPUNIT_ASSERT( utils::trim(str,wh1) == "" );
  CPPUNIT_ASSERT( utils::trim(str,wh2) == "" );

  str = "\t";

  CPPUNIT_ASSERT( utils::trim(str,wh1) == "" );
  CPPUNIT_ASSERT( utils::trim(str,wh2) == "\t" );

  str = " ";

  CPPUNIT_ASSERT( utils::trim(str,wh1) == "" );
  CPPUNIT_ASSERT( utils::trim(str,wh2) == "" );

}

void UtilsTest::test_chomp() {
  
  string str = "\tb\t\r";

  CPPUNIT_ASSERT( utils::chomp(str) == "\tb\t" );

  str = "\r  \t \r";

  CPPUNIT_ASSERT( utils::chomp(str) == "" );

  str = "";

  CPPUNIT_ASSERT( utils::chomp(str) == "" );

  str = "\r\n";

  CPPUNIT_ASSERT( utils::chomp(str) == "" );

  str = "a\n\r\n";

  CPPUNIT_ASSERT( utils::chomp(str) == "a" );

  str = "a\n";

  CPPUNIT_ASSERT( utils::chomp(str) == "a" );

}

void UtilsTest::test_split() {

  string str = " ab, c  , def,gh,,i j, ";

  vector<string> vec = utils::split(str,','," ");

  CPPUNIT_ASSERT( vec.size() == 7 );

  CPPUNIT_ASSERT( vec[0] == "ab" );
  CPPUNIT_ASSERT( vec[1] == "c" );
  CPPUNIT_ASSERT( vec[2] == "def" );
  CPPUNIT_ASSERT( vec[3] == "gh" );
  CPPUNIT_ASSERT( vec[4] == "" );
  CPPUNIT_ASSERT( vec[5] == "i j" );
  CPPUNIT_ASSERT( vec[6] == "" );

}

void UtilsTest::test_permute() {

  distributions::Random rand(0);

  num_t initData[] = {1.0,3.1,2.2,4.2,4.1,6.5,7.5,3,2};

  // Repeat the test 5 times
  for ( size_t i = 0; i < 5; ++i ) {
    
    vector<datadefs::num_t> data(initData,initData+8);
    vector<datadefs::num_t> dataOrig = data;
    
    utils::permute(data,&rand);
    
    bool anyChange = false;
    
    for ( size_t i = 0; i < data.size(); ++i ) {
      if ( data[i] != dataOrig[i] ) anyChange = true;
    }
    
    CPPUNIT_ASSERT( anyChange );
    
    sort(data.begin(),data.end());
    sort(dataOrig.begin(),dataOrig.end());
    
    for ( size_t i = 0; i < data.size(); ++i ) {
      CPPUNIT_ASSERT( data[i] == dataOrig[i] );
    }
  }

}

void UtilsTest::test_splitRange() {

  //vector<size_t> ics = utils::range(10);

  vector<vector<size_t> > icsSets = utils::splitRange(10,1);

  CPPUNIT_ASSERT( icsSets.size() == 1 );

  for ( size_t i = 0; i < icsSets.size(); ++i ) {
    for ( size_t j = 0; j < icsSets[i].size(); ++j ) {
      CPPUNIT_ASSERT( icsSets[i][j] == j );
    }
  }


}

void UtilsTest::test_strv2catv() {
  vector<string> strvec(51,"");
  vector<datadefs::num_t> catvec(51,0.0);
  map<string,datadefs::num_t> mapping;
  map<datadefs::num_t,string> backMapping;
  strvec[0] = "a";
  strvec[1] = "b";
  strvec[2] = "c";
  strvec[3] = "A";
  strvec[50] = "NaN";
  utils::strv2catv(strvec, catvec, mapping, backMapping);

  CPPUNIT_ASSERT(catvec[0] == 0.0);
  CPPUNIT_ASSERT(catvec[1] == 1.0);
  CPPUNIT_ASSERT(catvec[2] == 2.0);
  CPPUNIT_ASSERT(catvec[3] == 3.0);
  for (int i = 4; i < 50; ++i) {
    CPPUNIT_ASSERT(4.0);
  }
  CPPUNIT_ASSERT(datadefs::isNAN(catvec[50]));

  CPPUNIT_ASSERT( mapping["a"] == 0.0 );
  CPPUNIT_ASSERT( mapping["b"] == 1.0 );
  CPPUNIT_ASSERT( mapping["c"] == 2.0 );
  CPPUNIT_ASSERT( mapping["A"] == 3.0 );
  CPPUNIT_ASSERT( mapping[""]  == 4.0 );

  CPPUNIT_ASSERT( backMapping[0.0] == "a" );
  CPPUNIT_ASSERT( backMapping[1.0] == "b" );
  CPPUNIT_ASSERT( backMapping[2.0] == "c" );
  CPPUNIT_ASSERT( backMapping[3.0] == "A" );
  CPPUNIT_ASSERT( backMapping[4.0]  == "" );

}

void UtilsTest::test_strv2numv() {
  vector<string> strvec(51,"3.0");
  vector<datadefs::num_t> catvec(51,0.0);
  strvec[0] = "0.0";
  strvec[1] = "1.0";
  strvec[2] = "2.0";
  strvec[3] = "0.00";
  strvec[50] = "NaN";
  utils::strv2numv(strvec, catvec);

  CPPUNIT_ASSERT(catvec[0] == 0.0);
  CPPUNIT_ASSERT(catvec[1] == 1.0);
  CPPUNIT_ASSERT(catvec[2] == 2.0);
  CPPUNIT_ASSERT(catvec[3] == 0.0);
  for (int i = 4; i < 50; ++i) {
    CPPUNIT_ASSERT(catvec[i] == 3.0);
  }
  CPPUNIT_ASSERT(datadefs::isNAN(catvec[50]));
}

void UtilsTest::test_sortDataAndMakeRef() {
  vector<datadefs::num_t> data;
  vector<size_t> refIcs;
  data.push_back(0.0);
  for (int i = 0; i < 50; ++i) {
    data.push_back(static_cast<datadefs::num_t>(i));
  }

  for (int i = 49; i > -1; --i) { // Deliberately of length data.size() - 1
    refIcs.push_back(static_cast<size_t>(i));
  } // This allocation should be irrelevant. We keep it here deliberately to
    //  ensure the results are flattened in a safe manner; we expect a bounds
    //  checker to complain violently should that not be the case.

  utils::sortDataAndMakeRef(true, data, refIcs);
  CPPUNIT_ASSERT(data[0] == 0.0);
  CPPUNIT_ASSERT(refIcs[0] == 0);
  for (int i = 1; i < 51; ++i) {
    CPPUNIT_ASSERT(data[i] == static_cast<datadefs::num_t>(i-1));
    CPPUNIT_ASSERT(refIcs[i] == i);
  }

  utils::sortDataAndMakeRef(false, data, refIcs);
  for (int i = 0; i < 50; ++i) {
    CPPUNIT_ASSERT(data[i] == static_cast<datadefs::num_t>(49-i));
    CPPUNIT_ASSERT(refIcs[i] == 50-i);
  }
  CPPUNIT_ASSERT(data[50] == 0.0);
  CPPUNIT_ASSERT(refIcs[50] == 0);

  // Check for correct behavior with an empty data list and arbitrary refIcs
  data.clear();
  utils::sortDataAndMakeRef(true, data, refIcs);
  CPPUNIT_ASSERT(data.size() == 0);
  CPPUNIT_ASSERT(refIcs.size() == 0);

  utils::sortDataAndMakeRef(false, data, refIcs);
  CPPUNIT_ASSERT(data.size() == 0);
  CPPUNIT_ASSERT(refIcs.size() == 0);

  // NaNs are not checked as sorting targets, as their behavior is currently undefined
}

void UtilsTest::test_sortFromRef() {
  vector<int> data(50,0);
  vector<size_t> refIcs(50,0);
  for (int i = 0; i < 50; ++i) {
    data[i] = i;
    refIcs[i] = 49-i;
  }

  utils::sortFromRef<int>(data,refIcs);
  for (int i = 0; i < 50; ++i) {
    CPPUNIT_ASSERT(data[i] == 49-i);
  }
}


// Registers the fixture into the test 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( UtilsTest );

#endif // UTILSTEST_HPP
