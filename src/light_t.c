//! \file light_t.c
//! Correct for light-travel time.
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
#include <stdio.h>
#include <stdlib.h>

#include "aephem.h"


//------------------------------------------------------------------------------
//! Get the apparent position of an object, corrected for light-travel time.
//! The routine does three iterations of light-time correction.  Gravitational
//! retardation from the sun is neglected.  (However, it can be added if
//! necessary by using ae_light_t() directly.)
//!
//! This routine uses orbital elements.  For JPL ephemerides, use
//! ae_light_t_jpl().
//!
//! \param jd_tt The Julian date in TT.
//! \param o The heliocentric rectangular position of the observer, in AU.
//! \param orb The orbital elements of the object.
//! \param q For returning the heliocentric rectangular position of the object 
//!          in AU. 
//------------------------------------------------------------------------------

void ae_light_t_orbit(double jd_tt, double o[], const struct ae_orbit_t *orb, 
                      double q[]) {
  int i;
  double t;

  ae_kepler(jd_tt, orb, q);
  for (i = 0; i < 2; i++) {
    t = ae_light_t(o, q, 0);
    ae_kepler(jd_tt - t, orb, q);
  }

  return;
}


//------------------------------------------------------------------------------
//! Get the apparent position of an object, corrected for light-travel time.
//! The routine does three iterations of light-time correction.  Gravitational
//! retardation from the sun is neglected.  (However, it can be added if
//! necessary by using ae_light_t() directly.)
//!
//! This function uses JPL ephemerides.  For orbital elements, use
//! ae_light_t_orbit().
//!
//! \param jh A handle initialised by ae_jpl_init().
//! \param jd_tt The Julian date in TT.
//! \param o The heliocentric rectangular position of the earth, in AU.
//! \param obj_num The object for which to get coordinates.  If a planetary 
//!                ephemeris file is being used, one can use the predefined
//!                #ae_ss_bodies_t constants for clarity.
//! \param q For returning the heliocentric rectangular position of the object 
//!          in AU. 
//! \param v_q For returning the velocity of the object in AU.  Set to NULL if
//!            you do not require the velocity.
//! \param is_planetary If 0, assume these are not planetary ephemerides and
//!                     barycentric coordinates are returned.  Otherwise, it 
//!                     will be assumed that a Solar System ephemeris file is 
//!                     being used and that the object numbers are as in 
//!                     #ae_ss_bodies_t; heliocentric coordinates will be
//!                     returned. 
//!
//! \return 0 on success; on failure, a negative error code from #ae_retcode_t.
//------------------------------------------------------------------------------

int ae_light_t_jpl(struct ae_jpl_handle_t *jh, double jd_tt, double o[], 
                    int obj_num, double q[], double v_q[], char is_planetary) {
  int i, j;
  double t;

  if ((i = ae_jpl_get_coords(jh, jd_tt, obj_num, q, v_q, is_planetary)))
    return i;
  for (i = 0; i < 2; i++) {
    t = ae_light_t(o, q, 0);
    if ((j = ae_jpl_get_coords(jh, jd_tt - t, obj_num, q, v_q, is_planetary)))
      return j;
  }

  return 0;
}

//------------------------------------------------------------------------------
//! Get light-travel time between two objects.
//!
//! Given the position of two objects in heliocentric equatorial rectangular
//! coordinates, returns the light-travel time between them.  To correct for
//! light-travel time, this should be done iteratively a couple of times.
//! The routines ae_light_t_jpl() and ae_light_t_orbit() do this automatically.
//!
//! This routine can approximate the effect due to the sun's gravitational
//! retardation (AA, p. 36).  However, it is usually a small effect.  It will
//! blow up if one of the objects is the sun.
//!
//! \param p The first object, in heliocentric rectangular coordinates, in AU.
//! \param q The second object, in heliocentric rectangular coordinates, in AU.
//! \param do_retardation If 0, neglect gravitational retardation from the sun;
//!                       otherwise, compute gravitational retardation from the 
//!                       sun.
//!
//! \return The light-travel time, in days.
//------------------------------------------------------------------------------

double ae_light_t(double p[], double q[], int do_retardation) {
  int i;
  double p_r, q_r, pq_r, x, t;

  pq_r = 0;
  for (i = 0; i < 3; i++) {
    x = p[i] - q[i];
    pq_r += x * x;
  }
  t = sqrt(pq_r) / AE_CLIGHT_AUD;
 
  // Do gravitational retardation, if requested.  It blows up if p or q is the
  // sun (i.e., at the origin).
  if (do_retardation) {
    p_r = 0;
    q_r = 0;
    for (i = 0; i < 3; i++) {
      p_r += p[i] * p[i];
      q_r += q[i] * q[i];
    }
    p_r = sqrt(p_r);
    q_r = sqrt(q_r);
    pq_r = sqrt(pq_r);
    t += 1.97e-8 * log((p_r + pq_r + q_r) / (p_r - pq_r + q_r)) / AE_CLIGHT_AUD;
  }

  return t;
}
