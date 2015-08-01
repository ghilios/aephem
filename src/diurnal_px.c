//! \file diurnal_px.c
//! Correct for diurnal parallax.
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
//! Correct for diurnal parallax.
//! See AA page D3.
//!
//! This function does not bother to calculate anything unless the equatorial
//! horizontal parallax is at least 0.005".
//!
//! \param last The local apparent sidereal time, in degrees.
//! \param tlat The geocentric latitude of the observer, in degrees.
//! \param trho The distance from the centre of the earth to the observer, in
//!             earth radii.
//! \param dist The earth-object distance, in AU.
//! \param ra The right ascension of the object, in degrees; this routine
//!           modifies this parameter to correct for diurnal parallax.
//! \param dec The declination of the object, in degrees; this routine modifies
//!            this parameter to correct for diurnal parallax.
//------------------------------------------------------------------------------

void ae_diurnal_parallax(double last, double tlat, double trho, double dist, 
                         double *ra, double *dec) {
  double cosdec, sindec, coslat, sinlat;
  double p[3], dp[3], x, y, z, D;

  // Don't bother with this unless the equatorial horizontal parallax
  // is at least 0.005".
  if (dist > 1758.8)
    return;

  cosdec = cos(AE_DTR * (*dec));
  sindec = sin(AE_DTR * (*dec));

  // Observer's astronomical latitude.
  x = tlat * AE_DTR;
  coslat = cos(x);
  sinlat = sin(x);

  // Convert to equatorial rectangular coordinates in which unit distance == 
  // earth radius.
  D = dist * AE_AU / (0.001 * AE_R_EARTH);
  p[0] = D * cosdec * cos(AE_DTR * (*ra));
  p[1] = D * cosdec * sin(AE_DTR * (*ra));
  p[2] = D * sindec;

  dp[0] = -trho * coslat * cos(AE_DTR * last);
  dp[1] = -trho * coslat * sin(AE_DTR * last);
  dp[2] = -trho * sinlat;

  x = p[0] + dp[0];
  y = p[1] + dp[1];
  z = p[2] + dp[2];
  D = sqrt(x * x + y * y + z * z);  // Topocentric distance.

  /* Recompute ra and dec */
  *ra = AE_RTD * ae_zatan2(x, y);
  *dec = AE_RTD * asin(z / D);

  return;
}
