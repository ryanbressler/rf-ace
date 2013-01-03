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
  void test_regularSplitterSeek();

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
  
  Treedata treeData("test_2by10_text_matrix.afm",'\t',':');

  uint32_t h;

  MurmurHash3_x86_32("c",1,0,&h);

  Node node,leftChild,rightChild;

  node.setSplitter("T:in",h,leftChild,rightChild);
  
  CPPUNIT_ASSERT( &leftChild == node.leftChild() );
  CPPUNIT_ASSERT( &rightChild == node.rightChild() );

  CPPUNIT_ASSERT( NULL == node.missingChild() );

  CPPUNIT_ASSERT( node.percolate(&treeData,0,1) == &rightChild );
  CPPUNIT_ASSERT( node.percolate(&treeData,1,1) == &rightChild );
  CPPUNIT_ASSERT( node.percolate(&treeData,2,1) == &rightChild );
  CPPUNIT_ASSERT( node.percolate(&treeData,3,1) == &rightChild );
  CPPUNIT_ASSERT( node.percolate(&treeData,4,1) == &rightChild );
  CPPUNIT_ASSERT( node.percolate(&treeData,5,1) == &leftChild );
  CPPUNIT_ASSERT( node.percolate(&treeData,6,1) == &leftChild );
  CPPUNIT_ASSERT( node.percolate(&treeData,7,1) == &leftChild );
  CPPUNIT_ASSERT( node.percolate(&treeData,8,1) == &leftChild );
  CPPUNIT_ASSERT( node.percolate(&treeData,9,1) == &leftChild );
  CPPUNIT_ASSERT( node.percolate(&treeData,10,1) == &leftChild );
  CPPUNIT_ASSERT( node.percolate(&treeData,11,1) == &leftChild );
  CPPUNIT_ASSERT( node.percolate(&treeData,12,1) == &leftChild );
  CPPUNIT_ASSERT( node.percolate(&treeData,13,1) == &leftChild );
  CPPUNIT_ASSERT( node.percolate(&treeData,14,1) == &leftChild );
  CPPUNIT_ASSERT( node.percolate(&treeData,15,1) == &rightChild );
  CPPUNIT_ASSERT( node.percolate(&treeData,16,1) == &rightChild );
  CPPUNIT_ASSERT( node.percolate(&treeData,17,1) == &rightChild );
  CPPUNIT_ASSERT( node.percolate(&treeData,18,1) == &rightChild );
  CPPUNIT_ASSERT( node.percolate(&treeData,19,1) == &rightChild );


}

void NodeTest::test_regularSplitterSeek() {

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
