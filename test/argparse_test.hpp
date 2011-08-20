#ifndef ARGPARSETEST_HPP
#define ARGPARSETEST_HPP

#include "argparse.hpp"
#include "errno.hpp"
#include <cppunit/extensions/HelperMacros.h>

class ArgParseTest : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE( ArgParseTest );
  CPPUNIT_TEST( test_getArgument );
  CPPUNIT_TEST( test_getFlag );
  CPPUNIT_TEST( test_mixedArgumentsAndFlags );
  CPPUNIT_TEST( test_spuriousArgc );
  CPPUNIT_TEST( test_spuriousArgv );
  CPPUNIT_TEST_SUITE_END();
  
public:
  void setUp();
  void tearDown();

  void test_getArgument();
  void test_getFlag();
  void test_mixedArgumentsAndFlags();
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
    (char*)"efg",
    (char*)"-h",
    (char*)"i",
    (char*)"j"};

  bool a = false;
  bool b = false;
  bool c = false;
  bool d = false;
  bool e = false;
  bool f = false;
  bool g = false;
  bool h = false;
  bool i = false;
  bool j = false;
  bool z = false;

  ArgParse ap(argc, argv);
  CPPUNIT_ASSERT(ap.getFlag("a","aaa",a) == true);
  CPPUNIT_ASSERT(ap.getFlag("b","bbb",b) == true);
  CPPUNIT_ASSERT(ap.getFlag("c","ccc",c) == true);
  CPPUNIT_ASSERT(ap.getFlag("d","ddd",d) == true);
  CPPUNIT_ASSERT(ap.getFlag("e","eee",e) == false);
  CPPUNIT_ASSERT(ap.getFlag("f","fff",f) == false);
  CPPUNIT_ASSERT(ap.getFlag("g","ggg",g) == false);
  CPPUNIT_ASSERT(ap.getFlag("h","hhh",h) == true);
  CPPUNIT_ASSERT(ap.getFlag("i","iii",i) == false);
  CPPUNIT_ASSERT(ap.getFlag("j","jjj",j) == false);
  CPPUNIT_ASSERT(ap.getFlag("z","zzz",z) == false);

  CPPUNIT_ASSERT(a == true);
  CPPUNIT_ASSERT(b == true);
  CPPUNIT_ASSERT(c == true);
  CPPUNIT_ASSERT(d == true);
  CPPUNIT_ASSERT(e == false);
  CPPUNIT_ASSERT(f == false);
  CPPUNIT_ASSERT(g == false);
  CPPUNIT_ASSERT(h == true);
  CPPUNIT_ASSERT(i == false);
  CPPUNIT_ASSERT(j == false);
  CPPUNIT_ASSERT(z == false);
}

void ArgParseTest::test_mixedArgumentsAndFlags() {
  const int argc = 6;
  char* const argv[] = {
    (char*)"test",
    (char*)"-abcd",
    (char*)"e",
    (char*)"--fff=g",
    (char*)"-h",
    (char*)"i"};

  // These options are deliberately confusable, to expose cases of misparsing
  bool a = false;
  bool b = false;
  bool c = false;
  string d("");
  bool e = false;
  string f("");
  bool g = false;
  string h("");
  bool i = false;

  ArgParse ap(argc, argv);
  CPPUNIT_ASSERT(ap.getFlag("a","aaa",a) == true);
  CPPUNIT_ASSERT(ap.getFlag("b","bbb",b) == true);
  CPPUNIT_ASSERT(ap.getFlag("c","ccc",c) == true);
  CPPUNIT_ASSERT(ap.getArgument<string>("d","ddd",d) == true);
  CPPUNIT_ASSERT(ap.getFlag("e","eee",e) == false);
  CPPUNIT_ASSERT(ap.getArgument<string>("f","fff",f) == true);
  CPPUNIT_ASSERT(ap.getFlag("g","ggg",g) == false);
  CPPUNIT_ASSERT(ap.getArgument<string>("h","hhh",h) == true);
  CPPUNIT_ASSERT(ap.getFlag("i","iii",i) == false);

  CPPUNIT_ASSERT(a == true);
  CPPUNIT_ASSERT(b == true);
  CPPUNIT_ASSERT(c == true);
  CPPUNIT_ASSERT(strcmp(d.c_str(),"e") == 0);
  CPPUNIT_ASSERT(e == false);
  CPPUNIT_ASSERT(strcmp(f.c_str(),"g") == 0);
  CPPUNIT_ASSERT(g == false);
  CPPUNIT_ASSERT(strcmp(h.c_str(),"i") == 0);
  CPPUNIT_ASSERT(i == false);
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
    CPPUNIT_FAIL("ArgParse::ArgParse didn't throw any exception: expected 'ERRNO_INVALID_ARGUMENT'");
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
    CPPUNIT_FAIL("ArgParse::ArgParse didn't throw any exception; expected 'ERRNO_ILLEGAL_MEMORY_ACCESS'");
  } catch (int e) {
    CPPUNIT_ASSERT(e == ERRNO_ILLEGAL_MEMORY_ACCESS);
  }

  // Deprecated, since allowing this logic to crash in system-specific ways is
  //  preferrable to the heuristic alternatives. For discussion on one such
  //  strategy, see:
  //   http://blogs.msdn.com/b/oldnewthing/archive/2006/09/27/773741.aspx
  //
  //  For a more in-depth discussion, see:
  //   http://stackoverflow.com/questions/496034/most-efficient-replacement-for-isbadreadptr
  /**
  char* const argv2[] = {
    (char*)"test",
    (char*)"-abcd",
    (char*)"-e",
    (char*)"-f",
    (char*)"g"};
  
  try {
    ArgParse ap(argc, argv2); 
    CPPUNIT_FAIL("ArgParse::ArgParse didn't throw any exception; expected 'ERRNO_ILLEGAL_MEMORY_ACCESS'");
  } catch (int e) {
    CPPUNIT_ASSERT(e == ERRNO_ILLEGAL_MEMORY_ACCESS);
  } */
}

// Registers the fixture into the test 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( ArgParseTest );

#endif // ARGPARSETEST_HPP
