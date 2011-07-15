#ifndef ARGPARSE_HPP
#define ARGPARSE_HPP

#include <cstdlib>
#include <string>
#include <vector>
#include <set>

#include <sstream>
#include <cassert>

using namespace std;

class ArgParse 
{
public:

  ArgParse(const int argc, char* argv[]);
  ~ArgParse();

  template <typename T> void getArgument(const string& shortName, const string& longName, T& argument)
  {
    bool foundShort = false;
    bool foundLong = false;
    size_t foundAtIdx = tokens_.size();
    for(size_t i = 0; i < tokens_.size(); ++i)
      {

	if(tokens_[i].name.compare(shortName) == 0)
	  {
	    foundShort = true;
	    foundAtIdx = i;
	  }

	if(tokens_[i].name.compare(longName) == 0)
	  {
	    foundLong = true;
	    foundAtIdx = i;
	  }
      }

    if(!foundShort && !foundLong)
      {
	//stringstream ss("0");
	//ss >> argument;
	//argument = false;
	return;
      }

    assert(foundShort^foundLong);

    if(tokens_[foundAtIdx].value == "")
      {
	if(sizeof(T) != sizeof(bool))
	  {
	    cerr << "error: arguments without value must be handles with logical (bool) type. -" << shortName << " / --" << longName << " is not. Quitting..." << endl;
	    exit(1);
	  }
	stringstream ss("1");
	ss >> argument;
	//argument = true;
	return;
      }

    stringstream ss(tokens_[foundAtIdx].value);
    ss >> argument;
    if(!ss.eof())
      {
	cerr << "error: incorrect value (" << tokens_[foundAtIdx].value <<") for the argument -" << shortName << " / --" << longName << ". Quitting..." << endl;
	exit(1);
      }
  }

private:

  bool getKeys(string& str, vector<string>& keys);
  bool foundRepeats();

  //JUST FOR TESTING
  void printTokens();

  struct Token
  {
    string name;
    string value;
  };

  vector<Token> tokens_;

};

#endif
