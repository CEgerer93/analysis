
#ifndef _operators_h_
#define _operators_h_

#include "io/adat_xmlio.h"
#include "io/adat_xml_group_reader.h"

namespace Pseudo
{
  template<typename T> std::ostream& operator<<(std::ostream& os, const std::vector<T>& d);
  std::ostream& operator<<(std::ostream& os, const XMLArray::Array<int>& x);
  /* template<typename T> std::ostream& operator<<(std::ostream& os, const XMLArray::Array<T>& x); */
  int operator*(const XMLArray::Array<int>& p, const std::vector<int>& d);
  std::string operator+(const std::string& s, const int n);
}
#endif
