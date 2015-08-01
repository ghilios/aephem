//! \file annual_ab.c
//! Correct for annual aberration.
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

#include <stdlib.h>
#include <math.h>

#include "aephem.h"

//------------------------------------------------------------------------------
//! Correct for the annual aberration.
//! See AA pages B17, B37, C24.
//!
//! \param v_earth The heliocentric rectangular velocity of Earth in AU per day.
//! \param p A unit vector pointing from the earth to the object; the corrected
//!          position is returned in this parameter.
//! \param direction To add aberration to \p p, pass #AE_ADD_ABERRATION; to
//!                  remove it, pass #AE_REMOVE_ABERRATION; if neither is
//!                  passed, then #AE_REMOVE_ABERRATION is assumed.
//!                  #AE_ADD_ABERRATION is appropriate, for example, if one
//!                  wishes to compute the apparent position of a star from a 
//!                  catalogue position, while #AE_REMOVE_ABERRATION is 
//!                  appropriate, for example, if one wants to compute the true 
//!                  position of an observed star.
//------------------------------------------------------------------------------

void ae_annual_aberration(double v_earth[], double p[], int direction) {
  double A, B, C;
  double betai, pV;
  double x[3], V[3];
  int i;

  // Calculate the velocity of the earth.
  betai = 0.0;
  pV = 0.0;
  for (i = 0; i < 3; i++) {
	A = v_earth[i] / AE_CLIGHT_AUD;
	V[i] = A;
	betai += A * A;
	pV += p[i] * A;
  }
  
  // Make the adjustment for aberration.
  betai = sqrt(1.0 - betai);
  if (direction == AE_ADD_ABERRATION)
    C = 1.0 + pV;
  else
    C = 1.0 - pV;
  A = betai / C;
  
  if (direction == AE_ADD_ABERRATION)
    B = (1.0 + pV / (1.0 + betai)) / C;
  else
    B = -(1.0 - pV / (1.0 + betai)) / C;

  for (i = 0; i < 3; i++) {
    C = A * p[i] + B * V[i];
    x[i] = C;
  }

  for(i = 0; i < 3; i++)
    p[i] = x[i];
}
