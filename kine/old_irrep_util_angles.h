/*
  For cases where momenta used irrep_util version that did not split Rref from Rlatt

  The below is pulled from hadron/irrep_util.h at https://github.com/JeffersonLab/adat on commit 182c0295b86e9c7b40399e47c5c6edc08fcafa95
*/
#ifndef _old_irrep_util_angles_h_
#define _old_irrep_util_angles_h_

#include "io/adat_xmlio.h"
#include "io/adat_xml_group_reader.h"
#include "hadron/irrep_util.h"

#include<map>

namespace SingleRotations
{
  struct angles_t
  {
    double alpha, beta, gamma;
  };
  // Momenta that were supported
  const std::map<XMLArray::Array<int>, int> supportedSingleRotationMoms = {
    { Hadron::Mom3d(1, 1, 0), 0 },
    { Hadron::Mom3d(1, 1, 0), 0 },
    { Hadron::Mom3d(0, 1, 1), 0 },
    { Hadron::Mom3d(1, 0, 1), 0 },
    { Hadron::Mom3d(1, -1, 0), 0 },
    { Hadron::Mom3d(0, 1, -1), 0 },
    { Hadron::Mom3d(-1, 0, 1), 0 },
    { Hadron::Mom3d(-1, 1, 0), 0 },
    { Hadron::Mom3d(0, -1, 1), 0 },
    { Hadron::Mom3d(1, 0, -1), 0 },
    { Hadron::Mom3d(-1, -1, 0), 0 },
    { Hadron::Mom3d(0, -1, -1), 0 },
    { Hadron::Mom3d(-1, 0, -1), 0 },

    { Hadron::Mom3d(1, 1, 1), 0},
    { Hadron::Mom3d(-1, 1, 1), 0},
    { Hadron::Mom3d(1, -1, 1), 0},
    { Hadron::Mom3d(1, 1, -1), 0},
    { Hadron::Mom3d(-1, -1, 1), 0},
    { Hadron::Mom3d(1, -1, -1), 0},
    { Hadron::Mom3d(-1, 1, -1), 0},
    { Hadron::Mom3d(-1, -1, -1), 0},

    { Hadron::Mom3d(0, 1, 2), 0},
    { Hadron::Mom3d(1, 2, 0), 0},
    { Hadron::Mom3d(2, 0, 1), 0},
    { Hadron::Mom3d(0, 2, 1), 0},
    { Hadron::Mom3d(2, 1, 0), 0},
    { Hadron::Mom3d(1, 0, 2), 0},
    { Hadron::Mom3d(0, 1, -2), 0},
    { Hadron::Mom3d(1, -2, 0), 0},
    { Hadron::Mom3d(-2, 0, 1), 0},
    { Hadron::Mom3d(0, -2, 1), 0},
    { Hadron::Mom3d(-2, 1, 0), 0},
    { Hadron::Mom3d(1, 0, -2), 0},
    { Hadron::Mom3d(0, -1, 2), 0},
    { Hadron::Mom3d(-1, 2, 0), 0},
    { Hadron::Mom3d(2, 0, -1), 0},
    { Hadron::Mom3d(0, 2, -1), 0},
    { Hadron::Mom3d(2, -1, 0), 0},
    { Hadron::Mom3d(-1, 0, 2), 0},
    { Hadron::Mom3d(0, -1, -2), 0},
    { Hadron::Mom3d(-1, -2, 0), 0},
    { Hadron::Mom3d(-2, 0, -1), 0},
    { Hadron::Mom3d(0, -2, -1), 0},
    { Hadron::Mom3d(-2, -1, 0), 0},
    { Hadron::Mom3d(-1, 0, -2), 0},

    { Hadron::Mom3d(1, 1, 2), 0},
    { Hadron::Mom3d(1, 2, 1), 0},
    { Hadron::Mom3d(2, 1, 1), 0},
    { Hadron::Mom3d(-1, 1, 2), 0},
    { Hadron::Mom3d(-1, 2, 1), 0},
    { Hadron::Mom3d(1, -1, 2), 0},
    { Hadron::Mom3d(1, 2, -1), 0},
    { Hadron::Mom3d(2, -1, 1), 0},
    { Hadron::Mom3d(2, 1, -1), 0},
    { Hadron::Mom3d(1, 1, -2), 0},
    { Hadron::Mom3d(1, -2, 1), 0},
    { Hadron::Mom3d(-2, 1, 1), 0},
    { Hadron::Mom3d(-1, -1, 2), 0},
    { Hadron::Mom3d(-1, 2, -1), 0},
    { Hadron::Mom3d(2, -1, -1), 0},
    { Hadron::Mom3d(-1, 1, -2), 0},
    { Hadron::Mom3d(-1, -2, 1), 0},
    { Hadron::Mom3d(1, -1, -2), 0},
    { Hadron::Mom3d(1, -2, -1), 0},
    { Hadron::Mom3d(-2, -1, 1), 0},
    { Hadron::Mom3d(-2, 1, -1), 0},
    { Hadron::Mom3d(-1, -1, -2), 0},
    { Hadron::Mom3d(-1, -2, -1), 0},
    { Hadron::Mom3d(-2, -1, -1), 0},

    { Hadron::Mom3d(1, 1, 3), 0},
    { Hadron::Mom3d(1, 3, 1), 0},
    { Hadron::Mom3d(3, 1, 1), 0},
    { Hadron::Mom3d(-1, 1, 3), 0},
    { Hadron::Mom3d(-1, 3, 1), 0},
    { Hadron::Mom3d(1, -1, 3), 0},
    { Hadron::Mom3d(1, 3, -1), 0},
    { Hadron::Mom3d(3, -1, 1), 0},
    { Hadron::Mom3d(3, 1, -1), 0},
    { Hadron::Mom3d(1, 1, -3), 0},
    { Hadron::Mom3d(1, -3, 1), 0},
    { Hadron::Mom3d(-3, 1, 1), 0},
    { Hadron::Mom3d(-1, -1, 3), 0},
    { Hadron::Mom3d(-1, 3, -1), 0},
    { Hadron::Mom3d(3, -1, -1), 0},
    { Hadron::Mom3d(-1, 1, -3), 0},
    { Hadron::Mom3d(-1, -3, 1), 0},
    { Hadron::Mom3d(1, -1, -3), 0},
    { Hadron::Mom3d(1, -3, -1), 0},
    { Hadron::Mom3d(-3, -1, 1), 0},
    { Hadron::Mom3d(-3, 1, -1), 0},
    { Hadron::Mom3d(-1, -1, -3), 0},
    { Hadron::Mom3d(-1, -3, -1), 0},
    { Hadron::Mom3d(-3, -1, -1), 0},

    { Hadron::Mom3d(2,2,1), 0},
    { Hadron::Mom3d(2,1,2), 0},
    { Hadron::Mom3d(1,2,2), 0},
    { Hadron::Mom3d(-2,2,1), 0},
    { Hadron::Mom3d(-2,1,2), 0},
    { Hadron::Mom3d(2,-2,1), 0},
    { Hadron::Mom3d(2,1,-2), 0},
    { Hadron::Mom3d(1,-2,2), 0},
    { Hadron::Mom3d(1,2,-2), 0},
    { Hadron::Mom3d(2,2,-1), 0},
    { Hadron::Mom3d(2,-1,2), 0},
    { Hadron::Mom3d(-1,2,2), 0},
    { Hadron::Mom3d(-2,-2,1), 0},
    { Hadron::Mom3d(-2,1,-2), 0},
    { Hadron::Mom3d(1,-2,-2), 0},
    { Hadron::Mom3d(-2,2,-1), 0},
    { Hadron::Mom3d(-2,-1,2), 0},
    { Hadron::Mom3d(2,-2,-1), 0},
    { Hadron::Mom3d(2,-1,-2), 0},
    { Hadron::Mom3d(-1,-2,2), 0},
    { Hadron::Mom3d(-1,2,-2), 0},
    { Hadron::Mom3d(-2,-2,-1), 0},
    { Hadron::Mom3d(-2,-1,-2), 0},
    { Hadron::Mom3d(-1,-2,-2), 0}
  };

  //----------------------------------------------------------------------------------
  //! Return rotation angles for a vector to a canonical direction
  angles_t momSingleRotationAngles(const XMLArray::Array<int>& mom);
}
#endif
