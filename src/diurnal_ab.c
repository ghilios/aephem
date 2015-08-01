//! \file diurnal_ab.c
//! Correct for diurnal aberration.
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
//! Correct for diurnal aberration.
//! This formula is less rigorous than the method used for annual aberration
//! (see ae_annual_aberration()).  However, the correction is small.
//!
//! \param last The local apparent sidereal time, in degrees.
//! \param tlat The geocentric latitude of the observer, in degrees.
//! \param trho The distance from the centre of the earth to the observer, in
//!             earth radii.
//! \param direction To add aberration to \p p, pass #AE_ADD_ABERRATION; to
//!                  remove it, pass #AE_REMOVE_ABERRATION; if neither is
//!                  passed, then #AE_REMOVE_ABERRATION is assumed.
//!                  #AE_ADD_ABERRATION is appropriate, for example, if one
//!                  wishes to compute the apparent position of a star from a 
//!                  catalogue position, while #AE_REMOVE_ABERRATION is 
//!                  appropriate, for example, if one wants to compute the true 
//!                  position of an observed star.
//! \param ra The right ascension of the object, in degrees; this routine 
//!           modifies this parameter to correct for diurnal aberration.
//! \param dec The declination of the object, in degrees; this routine modifies 
//!            this parameter to correct for diurnal aberration.
//------------------------------------------------------------------------------

void ae_diurnal_aberration(double last, double tlat, double trho, int direction,
                           double *ra, double *dec) {
  double lha, coslha, sinlha, cosdec, sindec;
  double coslat, N, D;

  lha = last - *ra;
  coslha = cos(AE_DTR * lha);
  sinlha = sin(AE_DTR * lha);
  cosdec = cos(AE_DTR * (*dec));
  sindec = sin(AE_DTR * (*dec));
  coslat = cos(AE_DTR * tlat);

  if (cosdec != 0.0)
    N = 1.5472e-6 * trho * coslat * coslha / cosdec;
  else
    N = 0.0;
  D = 1.5472e-6 * trho * coslat * sinlha * sindec;
  
  if (direction == AE_ADD_ABERRATION) {
    *ra += AE_RTD * N;
    *dec += AE_RTD * D;
  }
  else {
    *ra -= AE_RTD * N;
    *dec -= AE_RTD * D;
  }

  return;
}
