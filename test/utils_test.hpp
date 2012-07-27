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
  CPPUNIT_TEST( test_join );
  CPPUNIT_TEST( test_range );
  CPPUNIT_TEST( test_trim );
  CPPUNIT_TEST( test_chomp );
  CPPUNIT_TEST( test_split );
  CPPUNIT_TEST( test_permute );
  CPPUNIT_TEST_SUITE_END();
  
public:
  void setUp();
  void tearDown();
  
  void test_removeNANs();
  void test_parse();
  void test_str2();
  void test_join();
  void test_range();
  void test_trim();
  void test_chomp();
  void test_split();
  void test_permute();
  
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

void UtilsTest::test_join() {

  vector<string> foo;
  foo.push_back("a");
  foo.push_back("b");

  CPPUNIT_ASSERT( utils::join(foo.begin(),foo.end(),',') == "a,b" );
  CPPUNIT_ASSERT( utils::join(foo.begin(),foo.end(),'a') == "aab" );
  foo.pop_back();
  CPPUNIT_ASSERT( utils::join(foo.begin(),foo.end(),',') == "a" );
  CPPUNIT_ASSERT( utils::join(foo.begin(),foo.end(),'a') == "a" );
  foo.pop_back();
  CPPUNIT_ASSERT( utils::join(foo.begin(),foo.end(),',') == "" );
  CPPUNIT_ASSERT( utils::join(foo.begin(),foo.end(),'a') == "" );


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

  distributions::RandInt rand(0);

  num_t initData[] = {1.0,3.1,2.2,4.2,4.1,6.5,7.5,3,2};

  // Repeat the test 5 times
  for ( size_t i = 0; i < 5; ++i ) {
    
    vector<datadefs::num_t> data(initData,initData+8);
    vector<datadefs::num_t> dataOrig = data;
    
    utils::permute(data,rand);
    
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


// Registers the fixture into the test 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( UtilsTest );

#endif // UTILSTEST_HPP
