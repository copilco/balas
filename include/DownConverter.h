/**
 * @file   DownConverter.h
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
#ifndef DOWNCONVERTER_H
#define DOWNCONVERTER_H
#include <iostream>
#include <typeinfo>
#include <stdexcept>
using std::cerr;

/**
 * \class DownConverter
 * \brief template class: converts a pointer to B into a pointer to R if T<R> is derived from B
 *
 * \remarks
 * \li Situation: you have a base class B and a bunch of derived
 *            classes T\<R\> (where R is not void ;-)
 * \li Access: by means of conversion operator R()
 * \li Template Parameters: on construction it is given a _pointer_ to the base class and a reference to a variable of type R (so it knows what type of T<R> it points to).
 * \exception Incompatible_Type_Exception if ( typeid(R) != *myType )
 * \return (* static_cast\< *T\<R\> \> _pointer_)
 *
 * \sa Parameter_ABC
 */
template<typename Base, template<typename> class Derived>
class DownConverter {
  Base* base_ptr;
  const type_info* myType;
public:
  template<typename R>
  DownConverter(Base* _base_ptr, R& _vtref)
    : base_ptr(_base_ptr), myType(&typeid(_vtref)) {};

  struct Incompatible_Type_Exception : public std::runtime_error {
    Incompatible_Type_Exception()
      : std::runtime_error("Incompatible_Type_Exception") {};
  };

  /// the templated conversion operator, determined by l-value
  template<typename R>
  operator R() {
    //return(*dynamic_cast<Derived<R>*>(base_ptr)); //segfaults
    if ((* myType) == typeid(R)) {
      Derived<R>* ptr=static_cast<Derived<R>*>(base_ptr);
      if (ptr != base_ptr){ // assert that conversion worked.
        // this is kind of unnecessary since we checked the datatype R
        // contained in T<R>, but to be sure..
        base_ptr->convertMsg(myType->name(), typeid(R).name());
        throw Incompatible_Type_Exception();
      }
      return (*ptr);    // return the beast itself
    } else {
      base_ptr->convertMsg(myType->name(), typeid(R).name());
      throw Incompatible_Type_Exception();
    }
  }
  operator Base&() const { return *base_ptr; }
};

#endif /* DOWNCONVERTER_H */
