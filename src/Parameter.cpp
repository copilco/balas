/**
 * @file   Parameter.cpp
 * @author Arvid Requate <arvid.requate@gmail.com>
 * @date   Mo Jul  3 2006  
 *
 * Copyright (c) 2006 Arvid Requate (http://www.pks.mpg.de/cpp/Parameter/)
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 * 
 * Arvid Requate                                    \n
 * <arvid.requate@gmail.com>                        \n
 */
#include <string>
#include <map>
#include <iostream>
#include "Parameter.h"
#include <iomanip>

using namespace std;

/* ======================= static data ============================ */

map<string,ParameterMap::converter> ParameterMap::_contentMap; /**< the static core of the ParameterMap Class*/ 
string ParameterMap::_delimiter="="; // defaults
string ParameterMap::_comment="#";   // defaults

/* ----------------------- Friends -------------------------------- */

inline void printParameterMap( std::ostream& os, const ParameterMap & f, string _prefix) {
  const int offset = 20;
  os << _prefix << f._comment << " --------------------- <Configuration> ---------------------" << endl;
  // print date of run and directory
  for (ParameterMap::mapci it = f.begin(); it != f.end(); it++) {
    Parameter_ABC &p(it->second);   // short
    os << _prefix << it->first << " " << f._delimiter << " " << left << setw(offset) << p.getValueString();
    string comment = p.getComment();
    if (comment.size() != 0)
      os << " " << f._comment << " " << comment;
    os << endl;

  }
  os << _prefix << f._comment << " --------------------- </Configuration> --------------------" << endl;
}

/* ================== output operators ============================ */

/** 
 * This operator is used to save (redirect) a ParameterMap to (out)file stream
 * \todo : add comment management
 */
ofstream& operator<<( std::ofstream& os, const ParameterMap & f) {
  printParameterMap(os, f, "# ");
  return os;
}

/** 
 * This operator is used to save (redirect) a ParameterMap to output stream 
 * \todo : add comment management
 */
ostream& operator<<( std::ostream& os, const ParameterMap & f) {
  printParameterMap(os, f);
  return os;
}

/**
 * output stream operator for Parameter_ABC (Abstract Base Class)
 */
ostream& operator<<( ostream& os, const Parameter_ABC& rh ) {
  //return os << rh.getLabel() << " = " << rh.getValueString();
  return os << rh.getLabel() << " = " << rh.getValueString() << "\t\t# " << rh.getComment();
};

/* ======================= Implementations ======================== */

void ParameterMap::saveToFile(const string & filename)
{
    ofstream cfgFile(filename.c_str(), ios::out);
    if ( !cfgFile) throw FileNotCreated(filename);	
    cfgFile << (*this);
    cfgFile.close();	
}
void ParameterMap::init(ConfigFileParser& _cfg)
{
  for (mapci it = begin(); it != end(); it++) {
    ((Parameter_ABC&) it->second).init(_cfg);
  }
}

/** 
 * Faciliate access to a value associated to a key. Throw KeyNotFound when the key 
 * doesn't exist.
 * \note
 * Difference to map\<string,string\>::operator[] is the behaviour on KeyNotFound
 * @param key    the name of the key
 * @return       the value associated (may be I shall return an empty string if 
 *               the key does not exist). 
 */
/*
const Parameter_ABC* ParameterMap::operator[](const string & key) const {
    mapci it=find(key);
    if ( end()==it ) throw KeyNotFound(key); 
    return it->second; 

}
*/
ParameterMap::converter ParameterMap::operator[](const string & key) {
    mapi it=find(key);
    if ( end()==it ) throw KeyNotFound(key); 
    return it->second; 

}

