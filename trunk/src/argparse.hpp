#ifndef ARGPARSE_HPP
#define ARGPARSE_HPP

#include <getopt.h>
#include <cassert>
#include <iostream>
#include <sstream>
#include <string.h>

using namespace std;

/**
 * Generic argument parser that allows options to be checked
 *  speculatively. Currently wraps GNU getopt_long, which is completely
 *  thread-unsafe. Let the maintainer beware.
 */
class ArgParse {

public:
  ArgParse(const int argc, char* const argv[]) {
    argc_ = argc;
    argv_ = argv;
    opterr = 0; // Suppress invalid print messages from getopt, which trigger
                //  spuriously from this class when long options are issued.
  }
  ~ArgParse() {}

  /**
   * Shims the local version of ArgParse to GNU getopt_long. Changes the input
   *  assumptions for long arguments slightly, but makes the whole system more
   *  robust as a result.
   *
   *  Note that due to interface imparity, each pass of getArgument takes linear
   *  time with respect to the original command line. It also incurs the
   *  overhead of a stringstream should an option be found.
   *
   *  Contractual guarantees:
   *
   *  + The behavior of the extraction operator (>>) will be used with input
   *     type T. This may cause unexpected results if your type explicitly
   *     specifies an append instead of an overwrite for this
   *     operation. Declare the contents of returnVal carefully or redefine
   *     your type for these cases.
   *
   *  + Certain types specifiable for T may cause memory access violations that
   *     are difficult to debug. For example, specifying char* may throw a
   *     memory access violation at 'SOptarg >> returnVal'. It is expected and
   *     indeed required that your types be well-defined before passing them
   *     to this method.
   *
   *  + Attempting to pass a NULL pointer for any input value will not
   *     work. Don't do it.
   *
   *  Sets returnVal and returns true if an argument was found; false
   *  otherwise. !! TODO return a unified status code instead
   */
  template <typename T> bool getArgument(const char* shortName, const char* longName, T& returnVal) {

    assert(shortName != NULL);
    assert(longName != NULL);
    assert(strlen(shortName) == 1);
    assert(*longName != 0);

    const char cShortName = *shortName;
    const struct option long_options[] = {
      {longName, 1, NULL, cShortName}
    };
    
    int option_index = 0;
    int c = 0;
    bool rst = false;
    while (c != -1) {
      c = getopt_long (argc_, argv_, ALLOWED_ARGS, long_options, &option_index);
      if (c == cShortName) {
        stringstream sOptarg(optarg); // Streamify the arg received from getopt
        sOptarg >> returnVal;
        rst = true;
        break;
      }
    }
    optind = 0; // Reset the getopt parsing index to 0
    return rst;
  }

  template <typename T> bool getArgument(string& shortName, string& longName, T& returnVal) {
    return getArgument<T>(shortName.c_str(), longName.c_str(), returnVal);
  }

  /**
   * Shims the local version of ArgParse to GNU getopt_long. Changes the input
   *  assumptions for long arguments slightly, but makes the whole system more
   *  robust as a result.
   *
   *  Note that due to interface imparity, each pass of getArgument takes linear
   *  time with respect to the original command line.
   *
   *  Tests for the presence of a flag. Unlike the above, returnVal must be of
   *  type bool. !! TODO return a unified status code instead
   *
   *  Contractual guarantees:
   *
   *  + Attempting to pass a NULL pointer for any input value will not
   *     work. Don't do it.
   *
   */
  bool getFlag(const char* shortName, const char* longName, bool& returnVal) {
    
    assert(strlen(shortName) == 1);
    assert(*longName != 0);

    const char cShortName = *shortName;
    const struct option long_options[] = {
      {longName, 0, NULL, cShortName}
    };
    int option_index = 0;
    int c = 0;
    bool rst = false;
    returnVal = false;
    while (c != -1) {
      c = getopt_long (argc_, argv_, ALLOWED_FLAGS, long_options, &option_index);
      if (c == cShortName) {
        returnVal = true;
        rst = true;
        break;
      }
    }

    // Set returnVal to false to signal the non-presence of this flag
    if (!rst) {
      returnVal = false;
    }

    optind = 0; // Reset the getopt parsing index to 0
    return rst;
  }

  bool getFlag(string& shortName, string& longName, bool& returnVal) {
    return getFlag(shortName.c_str(), longName.c_str(), returnVal);
  }

  
private:
  int argc_;
  char* const* argv_;
  static const char* ALLOWED_FLAGS;
  static const char* ALLOWED_ARGS;
};

const char* ArgParse::ALLOWED_FLAGS =
"ABCDEFGHIJKLMNOPQRSTUVWXYZ\
abcdefghijklmnopqrstuvwxyz\
0123456789";

const char* ArgParse::ALLOWED_ARGS =
"A:B:C:D:E:F:G:H:I:J:K:L:M:N:O:P:Q:R:S:T:U:V:W:X:Y:Z:\
a:b:c:d:e:f:g:h:i:j:k:l:m:n:o:p:q:r:s:t:u:v:w:x:y:z:\
0:1:2:3:4:5:6:7:8:9:";

#endif
