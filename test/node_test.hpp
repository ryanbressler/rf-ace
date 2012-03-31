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
  CPPUNIT_TEST( test_numericalFeatureSplit );
  CPPUNIT_TEST( test_categoricalFeatureSplit );
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
  void test_numericalFeatureSplit();
  void test_categoricalFeatureSplit();
  void test_recursiveNDescendantNodes();

};

void NodeTest::setUp() {}
void NodeTest::tearDown() {}


void NodeTest::test_setSplitter() {
  
  size_t splitterIdx = 3;
  datadefs::num_t splitLeftLeqValue = 0.5; 
  
  //Splitter::Splitter splitter(0.5);
  
  Node::Node node;
    
  node.setSplitter(splitterIdx,"foo",splitLeftLeqValue);
  
  CPPUNIT_ASSERT( node.splitter_.idx == splitterIdx );
  CPPUNIT_ASSERT( node.splitter_.isNumerical );
  CPPUNIT_ASSERT( fabs(node.splitter_.leftLeqValue - splitLeftLeqValue) < datadefs::EPS );
  
  
}

void NodeTest::test_percolateData() {
  
  Node::Node node0;
  //Splitter splitter("foo",0.1);
  node0.setSplitter(1,"foo",0.1);
  CPPUNIT_ASSERT( node0.leftChild() == node0.percolateData(0.09) );
  CPPUNIT_ASSERT( node0.rightChild() == node0.percolateData(0.11) );
  CPPUNIT_ASSERT( NULL == node0.percolateData(datadefs::NUM_NAN));

  Node::Node node1;
  
  set<string> leftValues;
  set<string> rightValues;

  leftValues.insert("a");
  leftValues.insert("b");

  rightValues.insert("c");
  rightValues.insert("d");

  node1.setSplitter(1,"foo",leftValues,rightValues);

  CPPUNIT_ASSERT( node1.percolateData("a") == node1.leftChild() );
  CPPUNIT_ASSERT( node1.percolateData("b") == node1.leftChild() );
  CPPUNIT_ASSERT( node1.percolateData("c") == node1.rightChild() );
  CPPUNIT_ASSERT( node1.percolateData("d") == node1.rightChild() );

  CPPUNIT_ASSERT( node1.percolateData("foo") == NULL );

}

void NodeTest::test_getLeafTrainPrediction() {
}

void NodeTest::test_hasChildren() {
}

void NodeTest::test_recursiveNodeSplit() { 
}

void NodeTest::test_cleanPairVectorFromNANs() { 

}

void NodeTest::test_numericalFeatureSplit() { 

  Treedata treedata("test_2by8_numerical_matrix.tsv",'\t',':');

  CPPUNIT_ASSERT( treedata.getFeatureName(0) == "N:F1" );
  CPPUNIT_ASSERT( treedata.getFeatureName(1) == "N:F2" );

  vector<size_t> sampleIcs_left(0);
  vector<size_t> sampleIcs_right(8);
  datadefs::range(sampleIcs_right);
  
  datadefs::num_t splitValue;
  datadefs::num_t splitFitness;

  Node node;

  Node::GrowInstructions GI;
  GI.minNodeSizeToStop = 2;

  size_t targetIdx = 1;
  size_t featureIdx = 0;

  node.numericalFeatureSplit(&treedata,
			     targetIdx,
			     featureIdx,
			     GI,
			     sampleIcs_left,
			     sampleIcs_right,
			     splitValue,
			     splitFitness);

  //cout << splitValue << endl;

  //Splitter::Splitter splitter("foo",splitValue);

  //CPPUNIT_ASSERT( splitter.splitsLeft(3.1 - 0.00001) );
  //CPPUNIT_ASSERT( splitter.splitsRight(3.1) );

  CPPUNIT_ASSERT( fabs(splitFitness - 0.530492285084497) < 1e-10 );

  //datadefs::print(sampleIcs_left);
  //datadefs::print(sampleIcs_right);

  CPPUNIT_ASSERT(sampleIcs_left.size() == 4);
  for(size_t i = 0; i < sampleIcs_left.size(); ++i) {
    size_t idx = sampleIcs_left[i];
    CPPUNIT_ASSERT(idx == 1 || idx == 6 || idx == 0 || idx == 7);
  }

  CPPUNIT_ASSERT(sampleIcs_right.size() == 2);
  for(size_t i = 0; i < sampleIcs_right.size(); ++i) {
    size_t idx = sampleIcs_right[i];
    CPPUNIT_ASSERT(idx == 5 || idx == 4);
  }

}

void NodeTest::test_categoricalFeatureSplit() { 

  size_t n = 10;

  Treedata treedata("test_3by10_categorical_matrix.tsv",'\t',':');

  size_t featureIdx = 0;

  size_t targetIdx1 = 1; 
  size_t targetIdx2 = 2;

  Node node;
  //PartitionSequence PS;

  Node::GrowInstructions GI;
  GI.minNodeSizeToStop = 2;
  //GI.partitionSequence = &PS;

  vector<size_t> sampleIcs_left(0);
  vector<size_t> sampleIcs_right(n);
  datadefs::range(sampleIcs_right);
  
  set<num_t> splitValues_left,splitValues_right;
  datadefs::num_t splitFitness;
  
  node.categoricalFeatureSplit(&treedata,
			       targetIdx1,
			       featureIdx,
			       GI,
			       sampleIcs_left,
			       sampleIcs_right,
			       splitValues_left,
			       splitValues_right,
			       splitFitness);
   
  //datadefs::print<num_t>(splitValues_left);
  //datadefs::print<num_t>(splitValues_right);
  //datadefs::print<size_t>(sampleIcs_left);
  //datadefs::print<size_t>(sampleIcs_right);
  //datadefs::print<num_t>(treedata.getFeatureData(targetIdx1));
  //datadefs::print<num_t>(treedata.getFeatureData(featureIdx));
  //datadefs::print<num_t>(treedata.getFeatureData(targetIdx1,sampleIcs_left));
  //datadefs::print<num_t>(treedata.getFeatureData(targetIdx1,sampleIcs_right));
  
  set<string> rawSplitValues_left,rawSplitValues_right;

  for(set<num_t>::const_iterator it(splitValues_left.begin()); it != splitValues_left.end(); ++it ) {
    rawSplitValues_left.insert(treedata.getRawFeatureData(featureIdx,*it));
  }
  
  for(set<num_t>::const_iterator it(splitValues_right.begin()); it != splitValues_right.end(); ++it ) {
    rawSplitValues_right.insert(treedata.getRawFeatureData(featureIdx,*it));
  }

  node.setSplitter(0,"foo",rawSplitValues_left,rawSplitValues_right);
  
  CPPUNIT_ASSERT( sampleIcs_left.size() == sampleIcs_right.size() );
  
  CPPUNIT_ASSERT( node.percolateData("0") == node.leftChild() );
  CPPUNIT_ASSERT( node.percolateData("1") == node.leftChild() );
  CPPUNIT_ASSERT( node.percolateData("4") == node.leftChild() );

  //cout << iter << endl;
  
  for(size_t i = 0; i < sampleIcs_left.size(); ++i ) {
    //cout << " " << treedata.getRawFeatureData(targetIdx1,sampleIcs_left[i]); 
    CPPUNIT_ASSERT( treedata.getRawFeatureData(targetIdx1,sampleIcs_left[i]) == "1" );
    CPPUNIT_ASSERT( treedata.getRawFeatureData(targetIdx1,sampleIcs_right[i]) == "0" );
  }

  //cout << endl;
  
  //cout << splitFitness << endl;
  
  CPPUNIT_ASSERT( fabs(splitFitness - 1) < datadefs::EPS );
   
  sampleIcs_left.clear();
  sampleIcs_right.resize(n);
  datadefs::range(sampleIcs_right);
  
  splitValues_left.clear();
  splitValues_right.clear();
  //datadefs::num_t splitFitness;

  Node node2;

  node2.categoricalFeatureSplit(&treedata,
				targetIdx2,
				featureIdx,
				GI,
				sampleIcs_left,
				sampleIcs_right,
				splitValues_left,
				splitValues_right,
				splitFitness);
  
  //datadefs::print<num_t>(splitValues_left);
  //datadefs::print<num_t>(splitValues_right);
  //datadefs::print<size_t>(sampleIcs_left);
  //datadefs::print<size_t>(sampleIcs_right);
  //datadefs::print<num_t>(treedata.getFeatureData(targetIdx2));
  //datadefs::print<num_t>(treedata.getFeatureData(featureIdx));
  //datadefs::print<num_t>(treedata.getFeatureData(targetIdx2,sampleIcs_left));
  //datadefs::print<num_t>(treedata.getFeatureData(targetIdx2,sampleIcs_right));

  rawSplitValues_left.clear();
  rawSplitValues_right.clear();

  for(set<num_t>::const_iterator it(splitValues_left.begin()); it != splitValues_left.end(); ++it ) {
    rawSplitValues_left.insert(treedata.getRawFeatureData(featureIdx,*it));
  }

  for(set<num_t>::const_iterator it(splitValues_right.begin()); it != splitValues_right.end(); ++it ) {
    rawSplitValues_right.insert(treedata.getRawFeatureData(featureIdx,*it));
  }

  
  node2.setSplitter(0,"foo",rawSplitValues_left,rawSplitValues_right);
  
  CPPUNIT_ASSERT( !node2.splitter_.isNumerical );
  
  CPPUNIT_ASSERT( sampleIcs_left.size() == 4 );
  CPPUNIT_ASSERT( sampleIcs_right.size() == 6 );
  CPPUNIT_ASSERT( node2.percolateData("0") == node2.leftChild() );
  CPPUNIT_ASSERT( node2.percolateData("1") == node2.leftChild() );

  
  //if ( iter == 0 ) {
  CPPUNIT_ASSERT( fabs( splitFitness - 0.642857142857143 ) < 1e-10 );
  //} else {
  //CPPUNIT_ASSERT( fabs(splitFitness - 0.843750000000000 ) < 1e-10 );
  //}
  
}

//void NodeTest::test_splitFitness() { 
//}

void NodeTest::test_recursiveNDescendantNodes() {
  
}
// Registers the fixture into the test 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( NodeTest ); 

#endif // NODETEST_HPP
