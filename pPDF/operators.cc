
#include "operators.h"

namespace Pseudo
{
  template<typename T>
  std::ostream& operator<<(std::ostream& os, const std::vector<T>& d)
  {
    if (d.size() > 0)
      {
  	os << d[0];
      
  	for ( int i = 1; i < d.size(); ++i )
  	  {
  	    os << " " << d[i];
  	  }
      }
    return os;
  }

  // template<typename T>
  // std::ostream& operator<<(std::ostream& os, const XMLArray::Array<T>& x)
  // {
  //   if (x.size() > 0)
  //     {
  // 	os << x[0];
	
  // 	for ( int i = 1; i < x.size(); ++i )
  // 	  {
  // 	    os << " " << x[i];
  // 	  }
  //     }
  //   return os;
  // }

#warning "FIX ME TEMPLATE!"
  template<typename T1, typename T2>
  double operator*(const XMLArray::Array<T1>& p, const std::vector<T2>& d)
  {
#warning "No elegant check for ioffe time! Only capturing the sign of ioffe time!"
    double ioffe = p[2]*d[0];
    return ioffe;
  }

  // Scale an XMLArray::Array<T1>
  // template<typename T>
  // XMLArray::Array<T> operator*(T scale, const XMLArray::Array<int>& a)
  XMLArray::Array<double> operator*(double scale, const XMLArray::Array<int>& a)
  {
    XMLArray::Array<double> rescaled; rescaled.resize(a.size());
    for ( int i = 0; i < a.size(); ++i )
      rescaled[i] = scale*a[i];
    return rescaled;
  }

  // Handle inner product of two XMLArrays
  int operator*(const XMLArray::Array<int>& l, const XMLArray::Array<int>& r)
  {
    int inner = 0;
    if ( l.size() != r.size() )
      {
	std::cerr << "XML Arrays of different dim!" << std::endl;
	exit(3);
      }
    else
      {
	for ( int i = 0; i < l.size(); ++i )
	  inner += l[i]*r[i];
	  // {
	  //   for ( int j = 0; j < r.size(); ++j )
	  //     inner += l[i]*r[j];
	  // }
      }
    return inner;
  }

  // template<typename T>
  // std::vector<T> operator*=(std::vector<T>& v, int i)
  // {
  //   for ( std::vector<T>::iterator it = v.begin(); it != v.end(); ++it )
  //     (*it)*i;
  //   return v;
  // }

  template<typename T>
  std::vector<std::complex<T> > operator*=(std::vector<std::complex<T> >& v, std::complex<T> c)
  {
    for ( auto vi = v.begin(); vi != v.end(); ++vi )
      (*vi) *= c;
    return v;
  }
      

  std::string operator+(const std::string& s, const int n)
  {
    std::string s1 = s + std::to_string(n);
    return s1;
  }
}
