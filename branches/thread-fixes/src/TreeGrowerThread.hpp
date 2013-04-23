#ifndef TREEGROWERTHREAD_HPP
#define TREEGROWERTHREAD_HPP

#include <vector>

#include "rootnode.hpp"

using namespace std;

class TreeGrowerThread {
public:
  vector<RootNode*>& rootNodes_;
  TreeData* trainData_;
  const size_t targetIdx_;
  const ForestOptions* forestOptions_;
  const distributions::PMF* PMF_;
  distributions::Random* randomDist_;

  TreeGrowerThread(vector<RootNode*>& rootNodes,
		   TreeData* trainData,
		   const size_t targetIdx,
		   const ForestOptions* forestOptions,
		   const distributions::PMF* pmf,
		   distributions::Random* random);
  
  ~TreeGrowerThread();

  void operator()();
};

#endif
