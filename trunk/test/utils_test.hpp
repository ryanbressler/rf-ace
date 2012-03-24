#ifndef UTILSTEST_HPP
#define UTILSTEST_HPP

#include <cppunit/extensions/HelperMacros.h>
#include "utils.hpp"

class UtilsTest : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE( UtilsTest );
  CPPUNIT_TEST( test_removeNANs );
  CPPUNIT_TEST( test_parse );
  CPPUNIT_TEST( test_str2 );
  CPPUNIT_TEST_SUITE_END();
  
public:
  void setUp();
  void tearDown();
  
  void test_removeNANs();
  void test_parse();
  void test_str2();
  
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

// Registers the fixture into the test 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( UtilsTest );

#endif // UTILSTEST_HPP
