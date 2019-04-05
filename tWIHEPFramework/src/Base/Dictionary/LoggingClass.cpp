#include "SingleTopRootAnalysis/Base/Dictionary/LoggingClass.hpp"
int LoggingClass::_globalLogLevel=1;
int LoggingClass::_TAGWidth=20;
LoggingClass::LoggingClass():_defaultTAG(""),_localLogLevel(-1){

}
LoggingClass::LoggingClass(std::string defaultTAG):_defaultTAG(defaultTAG),_localLogLevel(-1){
}
LoggingClass::~LoggingClass(){
}
bool LoggingClass::e(std::string message){
  return e(_defaultTAG,message);
}
bool LoggingClass::e(std::string TAG,std::string message){
  return generic(0,TAG,message);
}
bool LoggingClass::w(std::string message){
  return w(_defaultTAG,message);
}
bool LoggingClass::w(std::string TAG,std::string message){
  return generic(1,TAG,message);
}
bool LoggingClass::i(std::string message){
  return i(_defaultTAG,message);
}
bool LoggingClass::i(std::string TAG,std::string message){
  return generic(2,TAG,message);
}
bool LoggingClass::d(std::string message){
  return d(_defaultTAG,message);
}
bool LoggingClass::d(std::string TAG,std::string message){
  return generic(3,TAG,message);
}
bool LoggingClass::generic(int logLevel,std::string TAG,std::string message){
  int thisLogLevel=(_localLogLevel==-1?_globalLogLevel:_localLogLevel);
  if (logLevel<=thisLogLevel){
    try{
      std::cout << std::left<<std::setw(_TAGWidth) <<TAG;
      std::cout << std::left<<std::setw(8)<<errorTypes[logLevel];
      std::cout << message <<std::endl;
    }
    catch(char * str){
      return false;
    }
  }
  return true;
}
void LoggingClass::SetTAGWidth(int width){
  _TAGWidth=width;
}
bool LoggingClass::SetGlobalLogLevel(int logLevel){
  if (logLevel>=0 && logLevel<=3){
      _globalLogLevel=logLevel;
  }
  else{
    return false;
  }
  return true;
}
bool LoggingClass::SetLocalLogLevel(int logLevel){
  if (logLevel>=0 && logLevel<=3){
      _localLogLevel=logLevel;
  }
  else{
    return false;
  }
  return true;
}
std::string LoggingClass::tostr(long longarg){
  std::stringstream placeholder;
  placeholder<<longarg;
  return placeholder.str();
}

std::string LoggingClass::tostr(double doublearg){
  std::stringstream placeholder;
  placeholder<<doublearg;
  return placeholder.str();
}

