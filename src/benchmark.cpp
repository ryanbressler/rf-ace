#include<cstdlib>
#include<cstdio>
#include<iostream>
#include<vector>
#include<string>
#include<algorithm>
//#include "node.hpp"
#include "datadefs.hpp"
//#include "treedata.hpp"

using namespace std;
using datadefs::num_t;

int main()
{
  vector<int> foo(10);
  for(size_t i = 0; i < foo.size(); ++i)
    {
      foo[i] = i;
    }

  vector<int> sub(3);
  sub[0] = 4;
  sub[1] = 6;
  sub[2] = 1;

  vector<float> result(10);

  vector<float>::iterator it;

  it = set_difference(foo.begin(),foo.end(),sub.begin(),sub.end(),result.begin());
    
  cout << distance(it,result.begin()) << endl;
  for(vector<int> )

  return(EXIT_SUCCESS);
}
