#ifndef ARGPARSETEST_H
#define ARGPARSETEST_H

#include "argparse.hpp"
#include "errno.hpp"
#include <cppunit/extensions/HelperMacros.h>

class ArgParseTest : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE( ArgParseTest );
  CPPUNIT_TEST( test_getArgument );
  CPPUNIT_TEST( test_getFlag );
  CPPUNIT_TEST( test_spuriousArgc );
  CPPUNIT_TEST( test_spuriousArgv );
  CPPUNIT_TEST_SUITE_END();
  
public:
  void setUp();
  void tearDown();

  void test_getArgument();
  void test_getFlag();
  void test_spuriousArgc();
  void test_spuriousArgv();
};

void ArgParseTest::setUp() {}
void ArgParseTest::tearDown() {}


void ArgParseTest::test_getArgument() {
  const int argc = 6;
  char* const argv[] = {
    (char*)"test",
    (char*)"-f",
    (char*)"oo",
    (char*)"--bar",
    (char*)"boz",
    (char*)"--boz=1"};

  string val1("");
  string val2("");
  int val3(0);
  
  ArgParse ap(argc, argv);
  CPPUNIT_ASSERT(ap.getArgument<string>("f","foo",val1) == true);
  CPPUNIT_ASSERT(ap.getArgument<string>("b","bar",val2) == true);
  CPPUNIT_ASSERT(ap.getArgument<int>("z","boz",val3) == true);

  CPPUNIT_ASSERT( strcmp(val1.c_str(), "oo") == 0 );
  CPPUNIT_ASSERT( strcmp(val2.c_str(), "boz") == 0 );
  CPPUNIT_ASSERT( val3 == 1 );
}

void ArgParseTest::test_getFlag() {
  const int argc = 6;
  char* const argv[] = {
    (char*)"test",
    (char*)"-abcd",
    (char*)"-e",
    (char*)"-f",
    (char*)"g",
    (char*)"h"};

  bool a = false;
  bool b = false;
  bool c = false;
  bool d = false;
  bool e = false;
  bool f = false;
  bool g = false;
  bool h = false;
  bool i = false;
  bool z = true;

  ArgParse ap(argc, argv);
  CPPUNIT_ASSERT(ap.getFlag("a","aaa",a) == true);
  CPPUNIT_ASSERT(ap.getFlag("b","bbb",b) == true);
  CPPUNIT_ASSERT(ap.getFlag("c","ccc",c) == true);
  CPPUNIT_ASSERT(ap.getFlag("d","ddd",d) == true);
  CPPUNIT_ASSERT(ap.getFlag("e","eee",e) == true);
  CPPUNIT_ASSERT(ap.getFlag("f","fff",f) == true);
  CPPUNIT_ASSERT(ap.getFlag("g","ggg",g) == false);
  CPPUNIT_ASSERT(ap.getFlag("h","hhh",h) == false);
  CPPUNIT_ASSERT(ap.getFlag("z","zzz",z) == false);

  CPPUNIT_ASSERT(a == true);
  CPPUNIT_ASSERT(b == true);
  CPPUNIT_ASSERT(c == true);
  CPPUNIT_ASSERT(d == true);
  CPPUNIT_ASSERT(e == true);
  CPPUNIT_ASSERT(f == true);
  CPPUNIT_ASSERT(g == false);
  CPPUNIT_ASSERT(h == false);
  CPPUNIT_ASSERT(i == false);
  CPPUNIT_ASSERT(z == false);
}

void ArgParseTest::test_spuriousArgc() {
  const int argc = -200;
  char* const argv[] = {
    (char*)"test",
    (char*)"-f",
    (char*)"oo",
    (char*)"--bar",
    (char*)"boz",
    (char*)"--boz=1"};
  
  try {
    ArgParse ap(argc, argv);
    CPPUNIT_FAIL("ArgParse::ArgParse didn't throw the expected exception: ERRNO_INVALID_ARGUMENT");
  } catch (int e) {
    CPPUNIT_ASSERT(e == ERRNO_INVALID_ARGUMENT);
  }
}

void ArgParseTest::test_spuriousArgv() {
  const int argc = 6;
  char* const argv[] = {
    (char*)"test",
    (char*)"-abcd",
    (char*)"-e",
    (char*)"-f",
    (char*)"g",
    NULL};

  try {
    ArgParse ap(argc, argv);
    CPPUNIT_FAIL("ArgParse::ArgParse didn't throw the expected exception: ERRNO_ILLEGAL_MEMORY_ACCESS");
  } catch (int e) {
    CPPUNIT_ASSERT(e == ERRNO_ILLEGAL_MEMORY_ACCESS);
  }

  char* const argv2[] = {
    (char*)"test",
    (char*)"-abcd",
    (char*)"-e",
    (char*)"-f",
    (char*)"g"};
  
  try {
    ArgParse ap(argc, argv2); 
    CPPUNIT_FAIL("ArgParse::ArgParse didn't throw the expected exception: ERRNO_ILLEGAL_MEMORY_ACCESS");
  } catch (int e) {
    CPPUNIT_ASSERT(e == ERRNO_ILLEGAL_MEMORY_ACCESS);
  }
}

// Registers the fixture into the test 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( ArgParseTest );

#endif // ARGPARSETEST_H
