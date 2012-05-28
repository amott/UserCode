#ifndef ReadConfig_h
#define ReadConfig_h

#include<fstream>
#include<string>
#include<vector>
#include<map>
#include<cstring>

#include <iostream>
using namespace std;

class ReadConfig{
 public:
  ReadConfig();
  ReadConfig(string);
  bool is_init(){return isInit;}
  int read(string);
  string getParameter(string s);
  vector<string> getTokens(string s,const char *tok);
  void printAll();
 private:
  int parseLine(string);  // returns 0 if no error, >0 for specific errors (defined in function)
  string stripComments(string);  //removes comments from the string
  string stripSpaces(string); //strips leading and trailing spaces
  bool isInit;
  map<string,string> parameters_;
};

ReadConfig::ReadConfig():
  isInit(false)
{

}

ReadConfig::ReadConfig(string s):
  isInit(false)
{
  read(s);
}

int ReadConfig::read(string s){
  ifstream file(s.c_str());
  string line;
  if(!file.is_open()) return -1;
  
  int lineNo = 0;
  int ret = 0;
  while(file.good()){
    getline (file,line);
    ret = parseLine(line);
    if(ret) break;  // if parseLine returns non-zero, error parsing line; raise error flag and don't set isInit
    lineNo++;
  }
  isInit = !ret; // isInit only true if return value is 0
  if(!isInit) cout << "ERROR reading cfg" << endl;
  return (isInit==true?0:10*lineNo+ret); // returns 0 or 10*linenumber+errorcode
}

int ReadConfig::parseLine(string s){
  /*
    return codes: [errorcode]
    0 -- no errors
    1 -- no variable 
    2 -- no value
  */

  //cout << s <<endl;
  string sc = stripComments(s);
  //cout << "sc: " << sc <<endl << "   " << sc.find_first_not_of(' ') << endl;
  if(sc.find_first_not_of(' ') == string::npos) return 0; // if there was nothing other than spaces on the line

  int eqPos = sc.find_first_of('='); // location of the assignment
  //cout << eqPos << endl;
  if(eqPos == sc.find_first_not_of(' ')) return 1; // can't start the line with =
  if(eqPos == string::npos) return 2; // and it must exist

  string var = stripSpaces(sc.substr(0,eqPos)); //variable
  string val = stripSpaces(sc.substr(eqPos+1));  //value
  //cout << var <<endl;
  //cout << val <<endl;

  if(var == "") return 1;
  if(val == "") return 2;

  parameters_[var] = val;
  return 0;
}

string ReadConfig::stripComments(string s){  //walk down the string and remove everything after the first # _unless its in a string 
  //return s.substr(0,s.find_first_of('#'));
  bool inQuote = false;
  string out;
  for(int i=0;i<s.size();i++){
    if(s[i] == '"') inQuote = !inQuote;
    if(s[i] == '#' && !inQuote) break;
    out+=s[i];
  }
  return out;
}

string ReadConfig::stripSpaces(string s){
  return s.substr(s.find_first_not_of(' '),1+s.find_last_not_of(' ')-s.find_first_not_of(' '));
}

string ReadConfig::getParameter(string s){
  if(parameters_.find(s) == parameters_.end()) return "";
  else return parameters_.find(s)->second;
  return "";
}

vector<string> ReadConfig::getTokens(string s,const char *tok){
  vector<string> out;
  char line[4000];
  strcpy(line,this->getParameter(s).c_str());
  char *chunk = strtok(line,tok);
  while(chunk){
    out.push_back(string(chunk));
    chunk = strtok(NULL,","); 
  }
  return out;
}

void ReadConfig::printAll(){
  cout << "Config Parameters are:" <<endl;
  map<string,string>::const_iterator it;
  for(it = parameters_.begin();it !=parameters_.end(); it++){
    cout << ">> k: " << it->first << "   v: " << it->second << endl;
  }
}
#endif
