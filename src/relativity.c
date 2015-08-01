//! \file relativity.c
//! Correct for light deflection due to solar gravitation.
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
//! Correct for light deflection due to solar gravitation.
//! Note that this blows up if the object is very near (or is) the sun.
//! Therefore, if the sun-object or sun-object distance is zero, no correction
//! is made.
//!
//! See AA page B37.
//!
//! \param p The unit vector from observer to an object; this routine returns
//!          the corrected vector in this parameter.
//! \param q The heliocentric ecliptic rectangular coordinates of the object.
//! \param o The heliocentric ecliptic rectangular coordinates of the observer.
//------------------------------------------------------------------------------

void ae_relativity(double p[], double q[], double o[]) {
  int i;
  double C, d_sun_obs, d_sun_obj, d_obj_obs, p_dot_q, o_dot_p, q_dot_o;
  double dp[3];

  d_sun_obs = 0;
  d_sun_obj = 0;
  d_obj_obs = 0;
  p_dot_q = 0;
  o_dot_p = 0;
  q_dot_o = 0;
  for (i = 0; i < 3; i++) {
    d_sun_obs += o[i] * o[i];
    d_sun_obj += q[i] * q[i];
    d_obj_obs += p[i] * p[i];
    p_dot_q += p[i] * q[i];
    o_dot_p += p[i] * o[i];
    q_dot_o += q[i] * o[i];
  }
  d_sun_obs = sqrt(d_sun_obs);
  d_sun_obj = sqrt(d_sun_obj);
  d_obj_obs = sqrt(d_obj_obs);
  if (d_sun_obj && d_obj_obs)
    p_dot_q /= d_sun_obj * d_obj_obs;
  if (d_sun_obs && d_obj_obs)
    o_dot_p /= d_sun_obs * d_obj_obs;
  if (d_sun_obj && d_sun_obs)
    q_dot_o /= d_sun_obj * d_sun_obs;

  // Only proceed if the distance between sun and observer and object is
  // non-zero.
  if (d_sun_obs && d_sun_obj) {
    C = 1.974e-8 / (d_sun_obs * (1.0 + q_dot_o));
    for (i = 0; i < 3; i++) {
      dp[i] = C * (p_dot_q * o[i] / d_sun_obs - o_dot_p * q[i] / d_sun_obj);
      p[i] += dp[i];
    }
  }
}

