#ifndef EXCEPTIONS_HPP
#define EXCEPTIONS_HPP

#include <exception>
#include <string>
#include <sstream>

#include "errno.hpp"

class RFACE_EXCEPTION : public exception {

public:

  RFACE_EXCEPTION(const ERRNO& errno, const std::string& note = "") throw():errno_(errno) {
    stringstream ss;
    ss << "ERRNO (" << errno_ << "): ";
    switch ( errno_ ) {
    case ERRNO::INVALID_ARGUMENT:
      ss << "invalid command-line argument.";
      break;
    case ERRNO::INVALID_VALUE:
      ss << "invalid command-line value.";
      break;
    case ERRNO::ILLEGAL_MEMORY_ACCESS:
      ss << "illegal memory access.";
      break;
    default:
      ss << "unknown exception!";
      break;
    }
    ss << " " << note;
    msg_ = ss.str();
  }

  ~RFACE_EXCEPTION() throw() {}

  virtual const char* what() const throw() {
    return msg_.c_str();
  }
  
private:
  std::string msg_;
  ERRNO errno_;
  
};

#endif
