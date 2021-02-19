
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
  std::ostream& operator<<(std::ostream& os, const XMLArray::Array<int>& x)
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

  // template<typename p, typename d>
  int operator*(const XMLArray::Array<int>& p, const std::vector<int>& d)
  {
#warning "No elegant check for ioffe time! Only capturing the sign of ioffe time!"
    int ioffe = p[2]*d[0];
    return ioffe;
  }

  std::string operator+(const std::string& s, const int n)
  {
    std::string s1 = s + std::to_string(n);
    return s1;
  }
}
