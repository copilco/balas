// $Id: Parameter.h 232 2006-06-16 16:56:28Z arvid $
/**
 * @file   Parameter.h
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
#ifndef PARAMETER_H
#define PARAMETER_H
#include <string>
#include <sstream>
#include <map>
#include <iostream>
#include <vector>
#include <complex>
#include <bits/cpp_type_traits.h>
#include <stdexcept>
#include "ConfigFileParser.h"
#include "DownConverter.h"

using namespace std;

namespace ParameterHelper{

/**
 * Exceptions for the FromStringTo<T> and FromStringToVector<T> conversion
 * functions
 */
  class Error : public runtime_error {
    public:
    Error(const string & str = "") : runtime_error(str) {}
  };
  class ConversionError : public Error {
    public:
    ConversionError(const string & str = "") : Error(str) {}
  };
  class SizeChangeWarning : public Error {
    public:
    SizeChangeWarning(const string & str = "") : Error(str) {}
  };

/* ----------------------- Utility functions ---------------------- */

/** 
 * \class FromStringTo
 * try to convert a string to a value of type T 
 * \remarks Example : int i = FromStringTo<int>(s)
 * @param str string to convert 
 * @return  value of type T corresponding to the conversion of s to this type
 */
template<typename T> 
inline T FromStringTo(const std::string & str) {
    istringstream is(str);
    T t;
    is >> t;

    bool failFlag=0;
    if(is.fail()) failFlag=1;

    // check if every character of the string has been read
    string buf;
    is >> buf;
    if(failFlag || buf.size() !=0) {
      string type=typeid(t).name();
      throw ConversionError(
        "\"" + str + "\" has a bad format for conversion to type " + type \
        +" in FromStringTo<" + type + ">");
    }

    return t;
}
template<> 
inline string FromStringTo<string>(const std::string & str) {
    return str;
}
/** specialization for bool
 * \li "false", "F", "no", "n", "0" are interpreted as false
 * \li "true", "T", "yes", "y", "1", "-1", or anything else are interpreted as true
*/
template<> 
inline bool FromStringTo<bool>(const std::string & str) {
  /*
  static const char whitespace[] = " \n\t\v\r\f";
  string::size_type i=str.find_first_not_of(whitespace);
  string::size_type f=str.find_last_not_of(whitespace) + 1U;
  string tmp=str.substr(i,f-i);  // trim
  */

  // capitalisation to faciliate comparisons
  string sub=str;
  for(string::iterator it = sub.begin(); it != sub.end(); it++) *it = toupper(*it);
  if( sub==string("FALSE") || sub==string("F") || sub==string("NO") || 
    sub==string("N") || sub==string("0") || sub==string("NONE") )
    return false;
  if( sub==string("TRUE") || sub==string("T") || sub==string("YES") || 
    sub==string("Y") || sub==string("1") )
    return true;
  throw ConversionError(
      "\"" + str + "\" in not a boolean expression");
}

/** 
 * \class FromStringToVector
 * try to convert a string to a vector containing type T 
 * \remarks Example : vector<int> i = FromStringToVector<int>(s)
 * @param s string to convert 
 * @return  vector container of type T corresponding to the conversion of s to this type
 */
template<typename T> 
class FromStringToVector : public vector<T> {
public:
  FromStringToVector(const std::string & str) {
  string::size_type loc = str.find( "(", 0 );
  if( loc == string::npos ) {
    throw ConversionError("opening bracket of vector missing");
  }
  loc+=1;
  //cout<<"loc="<<loc<<endl;
  string::size_type _end = str.find( ")", loc );
  if( _end == string::npos ) {
    throw ConversionError("closing bracket of vector missing");
  }
  //cout<<"end="<<_end<<endl;
  string::size_type last=loc;
  string::size_type index = str.find( ",", loc );
  try {
  while(index != string::npos){
    //cout<<"last="<<last<<endl;
    //cout<<"index="<<index<<endl;
      push_back(FromStringTo<T>(str.substr(last, index-last)));
    //cout<<back()<<endl;
    last=index+1;
    index = str.find( ",", last );
  }
  push_back(FromStringTo<T>(str.substr(last, _end-last)));
  } catch (ConversionError & e) {
    ostringstream s; s << this->size();
    throw ConversionError("[" + s.str() + "]: " + e.what());
  }		
  }
};

/** 
 * \class ToStringFrom
 * convert any value v of type T to string 
 * \remarks Example: string s = ToStringFrom(14)
 * @param v value to convert
 * @return  a string corresponding to v
 */
template<typename T> 
inline string ToStringFrom(const T & v) {
    ostringstream s;
    s << v;
    return s.str();
}
}	/* </namespace ParameterHelper> */

/* ======================= Template tools ========================= */

/* -------------------- ternery _if template ---------------------- */
///
/// generic ternery _if template:
/// used by select_type
///
template<bool,typename T,typename R>
struct _if {
   typedef T type;
};
template <typename T, typename R>
struct _if<false,T,R> {
   typedef R type;
};

/* -------------------- dereference operator ---------------------- */
///
/// dereference operator
/// used by select_type
///
/// the default case: do nothing, define type = R
template<typename R, template<typename> class T>
struct dereference {
   typedef R type;
};
/// the special case: if first parameter is T<R>, define type = R
///                   (instead of T<R>)
template<typename R, template<typename> class T>
struct dereference<T<R>, T > {
   typedef R type;
};

/* ======================= ParameterMap =========================== */

class Parameter_ABC;  // forward declaration
template<typename T> class Parameter;  // forward declaration
/**
 * ParameterMap is the registry for all Parameter Objects
 * this is coded as a Wrapper object containing a static map.
 * Since the internal map is static, all ParameterMap instances refer to the
 * same one. i.e. it is guaranteed that only one map exists, but several
 * ParameterMap objects can be instanciated locally in the code to access it
 */
class ParameterMap {
public:
  typedef DownConverter<Parameter_ABC, Parameter> converter; /**<The type of the objects in the static _contentMap */
  // Write configuration
  friend void printParameterMap( std::ostream& os, const ParameterMap & f, string _prefix="");
private:
  static map<string,converter> _contentMap; /**< the mapping (key, value) */ 
  static string _delimiter;	  /**< separator between key and value          */
  static string _comment;	  /**< character(s) used to introduce a comment */

public:
  /*  wrapper functions to access methods of contained _contentMap */
  // Return a string value associated to key (eventually throw KeyNotFound)
  //const Parameter_ABC * operator[](const string & key) const;
  converter operator[](const string & key);
  /// Redefinition of iterators types to faciliate navigation in _contentMap
  typedef map<string,converter>::iterator       mapi;
  typedef map<string,converter>::const_iterator mapci;
  mapci find(const string& str) const { return _contentMap.find(str); }
  mapi find(const string& str) { return _contentMap.find(str); }
  /// Return iterators to the beginning and the end of _contentMap
  mapci begin() const { return _contentMap.begin(); }
  mapci end()   const { return _contentMap.end();   }
  void add(string &label, converter ptr) {
    if ( _contentMap.count(label) ) {
      cerr<<"Error: parameter with label \""<<label<<"\" doubly defined"<<endl;
      throw ParameterHelper::Error();
    } else {
      /// register in static _contentMap
      _contentMap.insert(pair<string, converter>(label,ptr));
    }
  }
  void remove(string &label) {
    _contentMap.erase(label);
  }

  /// Update the configuration file with the content of the map
  void saveToFile(const string & filename);
  // late initialization of all map entries
  void init(ConfigFileParser& _cfg);
  /// Exception types when the file can't be created
  struct FileNotCreated {
    string filename;
    FileNotCreated(const string & f = "") : filename(f) {}
  };
  /// Exception types when a key is not found
  struct KeyNotFound {
    KeyNotFound(const string & k = "") {
      cerr<<"No Parameter with label \""<<k<<"\" found! Please make sure to declare it."<<endl;
    }
  };
};

/* ------------------------- Parameter_ABC ------------------------ */
/* -------------- Abstract Base Class for Parameter<T> ------------ */

/**
 * Abstract Base Class(ABC) for the Parameter\<T\> class,
 * necessary to be able to store the Parameters of different type
 * in the ParameterMap.
 *
 * @see Parameter\<T\>
 */
class Parameter_ABC { 
  ParameterMap pm;		/**<Access Object for static ParameterMap */

  Parameter_ABC(const Parameter_ABC&);          /**<hidden, never used */
  //Parameter_ABC operator=(const Parameter_ABC&);/**<hidden, never used */
protected: // the derived classes can access these things
  string label, comment;
public:
  //Parameter_DownConverter value;        /**< the DownConverter */
  template<typename T>
  Parameter_ABC(T& _vtref, string _label, string _comment) :
  label(_label), comment(_comment) {
    pm.add(label, ParameterMap::converter(this, _vtref));
  };
  virtual ~Parameter_ABC(){		// virtual destructor : This way delete(Parameter_ABC* x) calls the destructor of the derived class that x belongs to. Yet this routine is called implicitely after the call to the destructor of the derived class.
    //cerr<<"~Parameter_ABC("<<label<<")"<<endl;
    pm.remove(label);
  };

  string getLabel() const { return label; };
  string getComment() const { return comment; };
  virtual void setValue(const string & ) = 0;
  virtual void init(ConfigFileParser& _cfg) = 0;
  virtual string getValueString() const = 0;

  virtual void convertMsg(const char* f, const char* t) const {
    cerr<<"Cannot convert Parameter_ABC[\""+label+"\"] from type "<<f<<" to "<<t<<"\n"; 
    cerr<<"In fact the porgramm should never have printed this.. \n";
  }
};

/* ======================= Parameter<T> =========================== */

/* -------- additional cpp_type_traits for complex<double> -------- */
///
/// additional cpp_type_traits for complex\<double\>
/// used by select_type
///
//#undef complex
#if __GNUC__ == 4
namespace std
{
  template<>
    struct __is_integer<complex<double> >
    {
      enum { __value = 0 };
      typedef __false_type __type;
    };
  template<typename _Tp>
  struct is_complex
  {
    enum { __value = 0 };
    typedef __false_type __type;
  };
  template<>
    struct is_complex<complex<double> >
    {
      enum { __value = 1 };
      typedef __true_type __type;
    };
}
#elif __GNUC__ == 3
namespace std
{
  template<>
    struct __is_integer<complex<double> >
    {
      enum { _M_type = 0 };
    };
  template<typename _Tp>
  struct is_complex
  {
    enum { _M_type = 0 };
  };
  template<>
    struct is_complex<complex<double> >
    {
      enum { _M_type = 1 };
    };
}
#endif
//#define complex complex<double>

/* ---------------- select_type template tool --------------------- */
///
/// select_type template tool
/// determines the return types of the templated binary
/// arithemtic operators of Parameter\<T\>
///
#if __GNUC__ == 4
template<typename T, typename R>
struct select_type {
  typedef typename dereference<R, Parameter>::type real_R;
  typedef typename _if<!(
    (__is_integer<T>::__value&&!__is_integer<real_R>::__value) ||
    (__is_floating<T>::__value&&is_complex<real_R>::__value))
    , T, real_R >::type type;
};
#elif __GNUC__ == 3
template<typename T, typename R>
struct select_type {
  typedef typename dereference<R, Parameter>::type real_R;
  typedef typename _if<!(
    (__is_integer<T>::_M_type&&!__is_integer<real_R>::_M_type) ||
    (__is_floating<T>::_M_type&&is_complex<real_R>::_M_type))
    , T, real_R >::type type;
};
#endif

/* ------------------ best_type template tool --------------------- */
/**
 * general type selection template for binary arithmetic operators:
 * generally its good to return the type of the left (lhs) argument:
template<typename T, typename R>
struct best_type {
    typedef T type;
};
 **/
/**
 * special type selection templates for binary arithmetic operators:
 * if the right (rhs) argument is more general than the lhs one, return lhs type
template<>
struct best_type<int, double> {
    typedef double type;
};
template<>
struct best_type<int, Parameter<double> > {
    typedef double type;
};
#undef complex
template<>
struct best_type<double, complex<double> > {
    typedef complex<double> type;
};
#define complex complex<double>
 **/

/* ------------------- finally: Parameter<T> ---------------------- */

/**
 * \class Parameter
 */
template<typename T>
class Parameter : public Parameter_ABC{ 
  T _data;
public:
  Parameter(string _label, T _value, string _comment = "");
  Parameter(ConfigFileParser& _cfg, string _label, T _value, string _comment = "");
  virtual void init(ConfigFileParser& _cfg);
  /** 
   * Set the _data from a string
   * @param str  the string representation of the value
   * @return       nothing
   */
  virtual void setValue(const string & str );
  /** 
   * Get a string from the _data
   * @return       the string representation of the value
   */
  virtual string getValueString() const;
  /// \name Access operators
  //@{ 
  //const T operator()() const { return _data; };
  //T& operator()() { return _data; };
  /// Conversion operator, non-virtual, See C++FAQ-lite [23.8]
  operator T() const { return _data; };
  /// reference conversion operator (nonconst = mutator), useful?
  operator T&() { return _data; };
  //@} 
  /// \name Assignment operators
  //@{ 
  void operator()(T _value){ _data=_value; };
  void set(T _value){ _data=_value; };
  //@} 
  /// \name Binary Arithmetic operators
  //@{ 
  template <typename R>
    typename select_type<T, R>::type
      operator+(const R&rhs) const { return _data + rhs; }
  template <typename R>
    typename select_type<T, R>::type
      operator-(const R&rhs) const { return _data - rhs; }
  template <typename R>
    typename select_type<T, R>::type
      operator*(const R&rhs) const { return _data*rhs; }
  template <typename R>
    typename select_type<T, R>::type
      operator/(const R&rhs) const { return _data/rhs; }
  template <typename R> T operator% (const R& rhs) { return _data % rhs; };
  //@} 

  /// \name Relational operators
  //@{ 
  template <typename R> bool operator== (const R& rhs) const { return _data == rhs; };
  template <typename R> bool operator!= (const R& rhs) const { return _data != rhs; };
  template <typename R> bool operator> (const R& rhs) const { return _data > rhs; };
  template <typename R> bool operator>= (const R& rhs) const { return _data >= rhs; };
  template <typename R> bool operator< (const R& rhs) const { return _data < rhs; };
  template <typename R> bool operator<= (const R& rhs) const { return _data <= rhs; };
  //@} 
  /// \name Unary Arithmetic operators
  //@{ 
  T operator-() const { return -_data; }
  //@} 

  virtual void convertMsg(const char* f, const char* t) const {
    cerr<<"Cannot convert Parameter<"<<f<<">[\""+label+"\"] to "<<t<<"\n"; 
  }
};

/**
 *  \defgroup BinArithTempl Binary Arithmetic operators for templated types (viz complex<T>)
 *  for templated classes as complex<T> the following template templates help:
 *   \relates Parameter
 */
/*@{*/
template< template<class> class C, typename T >
  C<T> operator+( C<T> const& lhs, Parameter<T> const& rhs )
    { return lhs+rhs(); }
template< template<class> class C, typename T >
  C<T> operator-( C<T> const& lhs, Parameter<T> const& rhs )
    { return lhs-rhs(); }
template< template<class> class C, typename T >
  C<T> operator*( C<T> const& lhs, Parameter<T> const& rhs )
    { return lhs*rhs(); }
template< template<class> class C, typename T >
  C<T> operator/( C<T> const& lhs, Parameter<T> const& rhs )
    { return lhs/rhs(); }
/*@}*/

/* =================== Parameter<vector<T> > ====================== */

/**
 * specialized Parameter\<vector\<T\> \> class
 */
template<typename T>
class Parameter<vector<T> > : public Parameter_ABC{ 
  vector<T> _data;
  static bool sizeChangeException;
  static bool sizeChangeMessage;
public:
  Parameter(string _label, vector<T> _value, string _comment = "");
  Parameter(ConfigFileParser& _cfg, string _label, vector<T> _value, string _comment = "");
  virtual void init(ConfigFileParser& _cfg);
  /** 
   * Set the _data from a string
   * @param str  the string representation of the vector of values
   * @return       nothing
   */
  virtual void setValue(const string & str );
  /** 
   * Get a string from the _data
   * @return       the string representation of the vector of values
   */
  virtual string getValueString() const;
  /// \name Access operators
  //@{ 
  const vector<T>& operator()() const { return _data; };
  /// Conversion operators
  operator vector<T>() { return _data; };
  /// reference conversion operator (nonconst = mutator), useful?
  operator vector<T>&() const { return _data; };
  /// Subscript operators
  const T& operator[](const int i) const { return _data[i]; };
  // (nonconst = mutator !)
  T& operator[](const int i) { return _data[i]; };
  //@} 
  /// \name Assignment operators
  //@{ 
  void operator()(vector<T> _value){ _data=_value; };
  void set(vector<T> _value){ _data=_value; };
  //@} 
  size_t size() const { return _data.size(); };

  virtual void convertMsg(const char* f, const char* t) const {
    cerr<<"Cannot convert Parameter<vector<"<<f<<"> >[\""+label+"\"] to "<<t<<"\n"; 
  }
};
template<typename T> bool Parameter<vector<T> >::sizeChangeException=false;
template<typename T> bool Parameter<vector<T> >::sizeChangeMessage=true;

/* ================== output operators ============================ */

ofstream& operator<<( std::ofstream& , const ParameterMap &);
ostream& operator<<( std::ostream& , const ParameterMap &);

ostream& operator<<( ostream& , const Parameter_ABC& );
//template<typename T> ostream& operator<<( ostream& os, const Parameter<T>& rh ) { return os << rh(); };
    
/* ------------------ general conveniences ------------------------ */
/**
 * output stream operator for vector\<T\> datatypes
 */
template<typename T>
inline ostream& operator<<( ostream& os, const vector<T>& rh ) {
  os << "(";
  typename vector<T>::const_iterator it;
  for (it=rh.begin(); it != rh.end() - 1 ; ++it){
    os << *it << ", ";
  }
  os << *it << ")";
  return os;
};

/**
 * for "const char*"+string this helps
 */
inline string operator+(const char *lhs, const string rhs) { string t(lhs); return t+=rhs; };


/* ======================= Implementations ======================== */

/**
 * Parameter\<T\> class constructors
 */
template<typename T>
Parameter<T>::Parameter(string _label, T _value, string _comment) :
  Parameter_ABC(_value, _label, _comment), _data(_value) {
  };
template<typename T>
Parameter<T>::Parameter(ConfigFileParser& _cfg, string _label, T _value, string _comment) :
  Parameter_ABC(_value, _label, _comment), _data(_value) {
  ConfigFileParser::mapci it=_cfg.find(_label);
  if ( _cfg.end() != it ) {
    setValue(it->second);
  } else {
    cerr << "Warning: Parameter \"" << _label << "\" not found in config file " << _cfg.getFilename() <<endl;
  }
};

/**
 * Parameter\<T\> init(ConfigFileParser) for late initialization
 */
template<typename T>
void Parameter<T>::init(ConfigFileParser& _cfg) {
  ConfigFileParser::mapci it=_cfg.find(label);
  if ( _cfg.end() != it ) {
    setValue(it->second);
  } else {
    // no message for late initialization
    //cerr << "Warning: Parameter \"" << _label << "\" not found in config file " << _cfg.getFilename() <<endl;
  }
};

/** 
 * Set the _data from a string
 * @param str  the string representation of the _data
 * @return       nothing
 */
template<typename T>
inline void Parameter<T>::setValue(const string & str ) {
  try {
    _data = ParameterHelper::FromStringTo<T>( str );
  } catch (ParameterHelper::ConversionError & e) {
    cerr<<"Error: "<<e.what()<<endl;
    throw ParameterHelper::ConversionError();
  }		
};
/** 
 * Get a string from the _data
 * @return       the string representation of the _data
 */
template<typename T>
inline string Parameter<T>::getValueString() const {
  return ToStringFrom<T>(_data);
};

/**
 * Parameter\<vector\<T\> \> class constructors
 */
template<typename T>
Parameter<vector<T> >::Parameter(string _label, vector<T> _value, string _comment) :
  Parameter_ABC(_value, _label, _comment), _data(_value) {};
template<typename T>
Parameter<vector<T> >::Parameter(ConfigFileParser& _cfg, string _label, vector<T> _value, string _comment) :
  Parameter_ABC(_value, _label, _comment), _data(_value) {
  ConfigFileParser::mapci it=_cfg.find(_label);
  if ( _cfg.end() != it ) {
    sizeChangeException=true;
    try {
      setValue(it->second);
    } catch (ParameterHelper::SizeChangeWarning & e) {
      cerr<<"Warning: Reading config file "<<_cfg.getFilename()<<": "<<endl;
      cerr<<"\t"<<e.what()<<endl;
    }
    sizeChangeException=false;
  } else {
    cerr << "Warning: Parameter \"" << _label << "\" not found in config file " << _cfg.getFilename() <<endl;
  }
};

/**
 * Parameter\<vector\<T\>\> init(ConfigFileParser) for late initialization
 */
template<typename T>
void Parameter<vector<T> >::init(ConfigFileParser& _cfg) {
  ConfigFileParser::mapci it=_cfg.find(label);
  if ( _cfg.end() != it ) {
    sizeChangeException=true;
    try {
      setValue(it->second);
    } catch (ParameterHelper::SizeChangeWarning & e) {
      cerr<<"Warning: Reading config file "<<_cfg.getFilename()<<": "<<endl;
      cerr<<"\t"<<e.what()<<endl;
    }
    sizeChangeException=false;
  } else {
    // no message for late initialization
    //cerr << "Warning: Parameter \"" << _label << "\" not found in config file " << _cfg.getFilename() <<endl;
  }
};

/** 
 * Set the _data from a string
 * @param str  the string representation of the vector of values
 * @return       nothing
 */
template<typename T>
inline void Parameter<vector<T> >::setValue(const string & str) {
  try {
    size_t vs=_data.size();
    _data=ParameterHelper::FromStringToVector<T>( str );
    if (_data.size() != vs) {
      ostringstream size1; size1 << vs;
      ostringstream size2; size2 << _data.size();
      string warn("size of vector \""+getLabel()+"\" changed from "+size1.str()+" to "+size2.str());
      if(sizeChangeException) {
	throw ParameterHelper::SizeChangeWarning(warn);
      } else if(sizeChangeMessage)
        cerr<<"Warning: "<<warn<<endl;
    }
  } catch (ParameterHelper::ConversionError & e) {
    cerr<<"Error in ParameterHelper::FromStringToVector: "<<getLabel()<<e.what()<<endl;
    throw /*ParameterHelper::ConversionError()*/;
  }		
};
/** 
 * Get a string from the _data
 * @return       the string representation of the vector of values
 */
template<typename T>
inline string Parameter<vector<T> >::getValueString() const {
  return ParameterHelper::ToStringFrom<vector<T> >(_data);
};

#endif /* PARAMETER_H */
