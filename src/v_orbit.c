//! \file v_orbit.c
//! Calculate the velocity of an orbit.
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
//! Calculate the velocity of an object as it moves in its orbit.
//! The difference in velocity between this and assuming a circular orbit is 
//! only about 1 part in 10**4.  Note that this gives heliocentric, not 
//! barycentric, velocity.
//!
//! \param jd_tt The Julian date in TT.
//! \param orb Orbital elements for the object.
//! \param v For returning the heliocentric, rectangular velocity, in
//!          AU per day.
//------------------------------------------------------------------------------

void ae_v_orbit(double jd_tt, const struct ae_orbit_t *orb, double v[]) {
  double e[3], r[3], t;
  int i;

  // Calculate the current position of the earth.
  ae_kepler(jd_tt, orb, r);

  // Calculate heliocentric position of the earth, as of a short time ago.
  t = 0.005;
  ae_kepler(jd_tt - t, orb, e);

  // Calculate the velocity.
  for (i = 0; i < 3; i++)
    v[i] = (r[i] - e[i]) / t;
}
