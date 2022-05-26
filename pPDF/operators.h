
#ifndef _operators_h_
#define _operators_h_

#include<complex>
#include "io/adat_xmlio.h"
#include "io/adat_xml_group_reader.h"

namespace Pseudo
{
  template<typename T> std::ostream& operator<<(std::ostream& os, const std::vector<T>& d);
  /* std::ostream& operator<<(std::ostream& os, const XMLArray::Array<int>& x); */
  template<typename T> std::ostream& operator<<(std::ostream& os, const XMLArray::Array<T>& x);

  template<typename T1, typename T2>
    double operator*(const XMLArray::Array<T1>& p, const std::vector<T2>& d);

  int operator*(const XMLArray::Array<int>& l, const XMLArray::Array<int>& r);


  template<typename T> std::ostream& operator<<(std::ostream& os, const XMLArray::Array<T>& x)
    {
      if (x.size() > 0)
	{
	  os << x[0];

	  for ( int i = 1; i < x.size(); ++i )
	    {
	      os << " " << x[i];
	    }
	}
      return os;
    }


  template<typename T> std::vector<T> operator*=(std::vector<T>& v, T i)
    {
      for ( auto it = v.begin(); it != v.end(); ++it )
	(*it)*=i;
      return v;
    }


  template<typename T>
    std::vector<std::complex<T> > operator*=(std::vector<std::complex<T> >& v, std::complex<T> c);
  std::string operator+(const std::string& s, const int n);
}
#endif
