#ifndef NODETEST_HPP
#define NODETEST_HPP

#include <cppunit/extensions/HelperMacros.h>
#include "datadefs.hpp"
#include "treedata.hpp"
#include "node.hpp"
#include "errno.hpp"

class NodeTest : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE( NodeTest );
  CPPUNIT_TEST( test_setSplitter );
  CPPUNIT_TEST( test_percolateData );
  CPPUNIT_TEST( test_getLeafTrainPrediction );
  CPPUNIT_TEST( test_hasChildren );
  CPPUNIT_TEST( test_recursiveNodeSplit );
  CPPUNIT_TEST( test_cleanPairVectorFromNANs );
  CPPUNIT_TEST( test_recursiveNDescendantNodes );
  CPPUNIT_TEST_SUITE_END();
  
public:
  void setUp();
  void tearDown();

  void test_setSplitter();
  void test_percolateData();
  void test_getLeafTrainPrediction();
  void test_hasChildren();
  void test_recursiveNodeSplit();
  void test_cleanPairVectorFromNANs();
  void test_recursiveNDescendantNodes();

private:

  ForestOptions forestOptions_;
  size_t threadIdx_;

};

void NodeTest::setUp() {

  forestOptions_.setCARTDefaults();
  threadIdx_ = 0;

}

void NodeTest::tearDown() {}


void NodeTest::test_setSplitter() {
  
  size_t splitterIdx = 3;
  datadefs::num_t splitLeftLeqValue = 0.5; 
  datadefs::num_t leftFraction = 0.5;
  
  //Splitter::Splitter splitter(0.5);
  
  Node node,leftChild,rightChild;

  node.setSplitter("foo",splitLeftLeqValue,leftChild,rightChild);
  
  //CPPUNIT_ASSERT( node.splitterIdx() == splitterIdx );
  CPPUNIT_ASSERT( node.splitter_.type == Feature::Type::NUM );
  CPPUNIT_ASSERT( fabs(node.splitter_.leftLeqValue - splitLeftLeqValue) < datadefs::EPS );
  //CPPUNIT_ASSERT( fabs(node.splitter_.leftFraction - leftFraction) < datadefs::EPS );
  
}

void NodeTest::test_percolateData() {
  
  Node node0,node00,node01;
  //Splitter splitter("foo",0.1);
  node0.setSplitter("foo",0.1,node00,node01);
  //CPPUNIT_ASSERT( node0.leftChild() == node0.percolate(0.09) );
  //CPPUNIT_ASSERT( node0.rightChild() == node0.percolate(0.11) );
  //CPPUNIT_ASSERT( NULL == node0.percolateData(datadefs::NUM_NAN));

  Node node1,node10,node11;
  
  set<string> leftValues;
  set<string> rightValues;

  leftValues.insert("a");
  leftValues.insert("b");

  rightValues.insert("c");
  rightValues.insert("d");

  node1.setSplitter("foo",leftValues,rightValues,node10,node11);

  //CPPUNIT_ASSERT( node1.percolate("a") == node1.leftChild() );
  //CPPUNIT_ASSERT( node1.percolate("b") == node1.leftChild() );
  //CPPUNIT_ASSERT( node1.percolate("c") == node1.rightChild() );
  //CPPUNIT_ASSERT( node1.percolate("d") == node1.rightChild() );

  //CPPUNIT_ASSERT( node1.percolateData("foo") == NULL );

}

void NodeTest::test_getLeafTrainPrediction() {
}

void NodeTest::test_hasChildren() {
}

void NodeTest::test_recursiveNodeSplit() { 
}

void NodeTest::test_cleanPairVectorFromNANs() { 

}

void NodeTest::test_recursiveNDescendantNodes() {
  
}
// Registers the fixture into the test 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( NodeTest ); 

#endif // NODETEST_HPP
