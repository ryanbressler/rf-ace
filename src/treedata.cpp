#include "treedata.hpp"
#include<cstdlib>
#include<fstream>
#include<cassert>
#include<iostream>
#include<sstream>

using namespace std;

Treedata::Treedata(string fname, bool is_featurerows):
  targetidx_(-1)
{

  cout << "Treedata: reading matrix from file '" << fname << "'" << endl;

  //Initialize stream to read from file
  ifstream featurestream;
  featurestream.open(fname.c_str());
  assert(featurestream.good());

  string field;
  string row;

  //Remove upper left element from the matrix as useless
  getline(featurestream,field,'\t');

  //Count the number of columns
  getline(featurestream,row);
  stringstream ss(row);
  size_t ncols = 0;
  while(getline(ss,field,'\t'))
    {
      ++ncols;
    }

  //Count the number of rows
  size_t nrows = 0;
  while(getline(featurestream,row))
    {
      ++nrows;
    }

  //Reset streams and remove upper left element from the matrix as useless
  featurestream.clear();
  featurestream.seekg(0, ios::beg);
  getline(featurestream,field,'\t');
  ss.clear();
  ss.str("");

  cout << ncols << " columns and " << nrows << " rows read." << endl;

  if(is_featurerows)
    {
      this->transpose();
    }
}

Treedata::~Treedata()
{
}

void Treedata::transpose()
{
  
}



