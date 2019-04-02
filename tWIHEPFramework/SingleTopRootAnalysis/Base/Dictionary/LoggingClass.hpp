#ifndef _STRALOGGING_
#define _STRALOGGING_

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
const std::string errorTypes[4]={"ERROR:","WARNING:","INFO:","DEBUG:"};
class LoggingClass {
public:
  enum {
    ERROR=0,
    WARNING=1,
    INFO=2,
    DEBUG=3
  };
  LoggingClass();
  LoggingClass(std::string DefaultTAG);
  ~LoggingClass();
  void SetTAGWidth(int width);
  static bool SetGlobalLogLevel(int logLevel);
  bool SetLocalLogLevel(int logLevel);
  bool e(std::string TAG,std::string message);
  bool e(std::string message);
  bool w(std::string TAG,std::string message);
  bool w(std::string message);
  bool i(std::string TAG,std::string message);
  bool i(std::string message);
  bool d(std::string TAG,std::string message);
  bool d(std::string message);
  static std::string tostr(long longarg);
  static std::string tostr(double doublearg);
			  
private:
  std::string _defaultTAG;
  static int _TAGWidth;
  static int _globalLogLevel;
  int _localLogLevel;
  bool generic (int type, std::string TAG, std::string message);
  
};
#endif
