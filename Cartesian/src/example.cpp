#include "Parameter.h"
#include <iostream>
#include <vector>
#include <string>

using std::cout;
using std::vector;
using std::string;

Parameter<string> jobname("jobname","Run1");
Parameter<double> lambda_nm("lambda_nm",32, "Wavelength in nm");
//                                          ^^ optional comment
Parameter<vector<int> > point("point", vector<int>(3,0) );

main(){
// read configuration file
ConfigFileParser cfg( "example.cfg" );

ParameterMap pMap;
// initialize all Parameter<T> type objects from the config file
pMap.init( cfg );


cout<<"Wavenumber: "<< 1/( lambda_nm * 1e-9 ) <<"\n";

// string/map access:
double x; cfg.getValueInto(x, "lambda_nm"); // ConfigFileParser style
double y = pMap["lambda_nm"];               // ParameterMap style
string s = pMap["jobname"];                 // kind of "chameleon"

try {
  cout << "Trying to pass a Parameter<double> to an int:\n";
  int i = pMap["lambda_nm"];
} catch ( std::runtime_error e ) {
  cerr << "   ->    " << e.what() <<"\n";
}

cout<<pMap;
}
