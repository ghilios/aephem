//! \file math_util.c
//! This file contains mathematical utility functions.
//==============================================================================
// AEPHEM - an astronomical ephemeris and reduction library.
// Copyright 2008 Adam Hincks.
//
// This file is part of AEPHEM.
//
// AEPHEM is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// AEPHEM is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with AEPHEM.  If not, see <http://www.gnu.org/licenses/>.
//==============================================================================

#include <math.h>
#include <stdlib.h>

#include "aephem.h"


//------------------------------------------------------------------------------
//! Reduce x modulo 360 degrees.
//!
//! \param x The number to be modulo'd.
//!
//! \return A number in the range [0: 360.0].
//------------------------------------------------------------------------------

double ae_mod_360(double x) {
  long k;
  double y;

  k = (long)(x / 360.0);
  y = x - k * 360.0;
  while (y < 0.0)
    y += 360.0;
  while (y >= 360.0)
    y -= 360.0;

  return(y);
}


//------------------------------------------------------------------------------
//! Reduce x modulo 360 degrees, but centre at 0 degrees.
//!
//! \param x The number to be modulo'd.
//!
//! \return A number in the range (-180: 180].
//------------------------------------------------------------------------------

double ae_mod_180(double x) {
  double y;

  if ((y = ae_mod_360(x)) > 180)
    y -= 360;

  return(y);
}


//------------------------------------------------------------------------------
//! Reduce x modulo 2 pi.
//!
//! \param x The number to be modulo'd.
//!
//! \return The number in the range [0: 2 pi].
//------------------------------------------------------------------------------

double ae_mod_2pi(double x) {
  double y;

  y = floor(x / 2 / M_PI);
  y = x - y * 2 * M_PI;
  while (y < 0.0)
    y += 2 * M_PI;
  while (y >= 2 * M_PI)
    y -= 2 * M_PI;

  return(y);
}


//------------------------------------------------------------------------------
//! Quadrant correct inverse circular tangent.
//! Returns radian angle between 0 and +2PI whose tangent is y/x.
//! 
//! Cephes Math Library Release 2.0:  April, 1987
//! Copyright 1984, 1987 by Stephen L. Moshier.
//! Direct inquiries to 30 Frost Street, Cambridge, MA 02140
//! Certain routines from the Library, including this one, may be used and 
//! distributed freely provided this notice is retained and source code is 
//! included with all distributions.
//!
//! \param x The denominator.
//! \param y The numerator.
//!
//! \return The quadrant correct inverse circular tangent.
//------------------------------------------------------------------------------

double ae_zatan2(double x, double y) {
  double z, w;
  short code;

  code = 0;

  if (x < 0.0)
    code = 2;
  if (y < 0.0)
    code |= 1;

  if (x == 0.0) {
    if (code & 1)
      return 1.5 * M_PI;
    if (y == 0.0)
      return 0.0;
    return 0.5 * M_PI;
  }

  if (y == 0.0) {
    if (code & 2)
      return M_PI;
    return 0.0;
  }

  switch (code) {
    default: case 0:
      w = 0.0;
      break;
    case 1:
      w = 2.0 * M_PI;
      break;
    case 2: case 3:
      w = M_PI;
      break;
    }

  z = atan(y / x);

  return w + z;
}
