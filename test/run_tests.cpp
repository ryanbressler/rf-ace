#ifndef TEST__
#define TEST__
#endif

//#include "partitionsequence_test.hpp"
#include "argparse_test.hpp"
#include "datadefs_test.hpp"
#include "utils_test.hpp"
#include "math_test.hpp"
#include "node_test.hpp"
#include "rootnode_test.hpp"
#include "stochasticforest_test.hpp"
//#include "splitter_test.hpp"
#include "treedata_test.hpp"
#include <cppunit/CompilerOutputter.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>


// Adapted code from the CppUnit documentation, intermixed with updates 
//  by Larry Hosken, who seems to enjoy Comic Sans more than I do:
//  http://lahosken.san-francisco.ca.us/

int RunTests(void) {
  // Get the top level suite from the registry
  CppUnit::Test *suite = CppUnit::TestFactoryRegistry::getRegistry().makeTest();

  // Adds the test to the list of test to run
  // CppUnit::TextUi::TestRunner runner;
  CppUnit::TextTestRunner runner;
  runner.addTest( suite );

  // Change the default outputter to a compiler error format outputter
  runner.setOutputter( new CppUnit::CompilerOutputter( &runner.result(),
                                                       std::cerr ) );
  // Run the tests.
  bool wasSucessful = runner.run();

  // Return error code 1 if the one of test failed.
  return wasSucessful ? 0 : 1;
}

int main(void) {
  return RunTests();
}
