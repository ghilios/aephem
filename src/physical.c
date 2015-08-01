//! \file physical.c
//! This file contains functionns for computing physical ephemerides.
//==============================================================================
// AEPHEM - an astronomical ephemeris and reduction library.
// Copyright 2012 Adam Hincks.
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
#include <time.h>

#include "aephem.h"


//------------------------------------------------------------------------------
//! Calculate the solid angle of an observed disc.
//! It is assumed that the disc is elliptical.  The semi-minor axis of the
//! projected ellipse can be obtained from ae_disc_semiminor().  Note that this
//! calculation assumes that the size of the disc is small compared to the 
//! distance of the observer.
//!
//! \param a The semi-major axis of the disc, in kilometres.
//! \param b The semi-minor axis of the disc, in kilometres.
//! \param dist The distance from the observer to the planet, in AU.
//!
//! \return The solid angle of the disc, in steradians.
//------------------------------------------------------------------------------

double ae_disc_solid_angle(double a, double b, double dist) {
  return M_PI * a * b / dist / dist / AE_AU / AE_AU;
}

//------------------------------------------------------------------------------
//! Calculate the solid angle of an observed disc.
//! This simplified routine is potentially less efficient, since the distance to
//! the planet is calculated as well as the sub-earth latitude to find the
//! semi-minor axis of the projected disc---ae_disc_solid_angle() and
//! ae_disc_semiminor() can make use of precomputed values.
//!
//! \param jd_ut1 The Julian date in UT1.
//! \param o_orb The orbital elements of the observer.
//! \param q_orb The orbital elements of the object being observed.
//! \param phys The physical elements of the object being observed.
//!
//! \return The solid angle of the disc, in steradians.
//------------------------------------------------------------------------------

double aes_disc_solid_angle(double jd_ut1, const struct ae_orbit_t *o_orb,
                            const struct ae_orbit_t *q_orb,
                            const struct ae_physical_t *phys) {
  double lat, lon, dist, a, b;
  
  // Get the sub-observer point.
  aes_subobs_point(jd_ut1, o_orb, q_orb, phys, &lat, &lon, &dist);

  // Get the semi-minor axis of the projected disc.
  a = phys->r_eq;
  b = ae_disc_semiminor(a, phys->r_pole, lat);

  // Finally, compute the solid angle.
  return ae_disc_solid_angle(a, b, dist);
}

//------------------------------------------------------------------------------
//! Calculate the semi-minor axis of an observed disc.
//! When an oblate spheroid is viewed, it appears as an ellipse (to a
//! sufficiently distant observer).  This routine calculates the semi-minor axis
//! (i.e., polar radius) of the projected ellipse.
//!
//! \param a The major axis (i.e., equatorial radius) of the planet.
//! \param c The minor axis (i.e., polar radius) of the planet.
//! \param lat The sub-observer latitude, in degrees.
//!
//! \return The semi-minor axis of the projected ellipse.  Note that the
//!         semi-major axis of the ellipse is the same as the input value.
//------------------------------------------------------------------------------

double ae_disc_semiminor(double a, double c, double lat) {
  double cos_lat, sin_lat;

  cos_lat = cos(lat * AE_DTR);
  sin_lat = sin(lat * AE_DTR);

  return sqrt(a * a * sin_lat * sin_lat + c * c * cos_lat * cos_lat);
}


//------------------------------------------------------------------------------
//! Determine whether a body's rotation is retrograde.
//! This function merely looks at the sign of the variable w_d in the struct
//! #ae_physical_t:  negative values are retrograde.
//!
//! \param phys The body to examine.
//!
//! \return 1 if the rotation is retrograde; 0 if the rotation is prograde.
//------------------------------------------------------------------------------

int ae_is_retrograde(const struct ae_physical_t *phys) {
  if (phys->w_d < 0)
    return 1;
  else
    return 0;
}

//------------------------------------------------------------------------------
//! Calculate the sub-observer point.
//! This is a simplified version of ae_subobs_point().  It automatically
//! calculates the requisite unit vectors, given the orbital and physical
//! elements of the bodies.
//!
//! The sub-observer point is the latitude and longitude on an observed body at
//! the observer is directly over-head.
//!
//! This routine includes corrections for light-travel time, annual 
//! aberration, precession and nutation.
//!
//! WARNING:  the longitude computed by this function is, for some reason, only
//! accurate to tens of arc minutes when compared to HORIZONS (and I don't know
//! why. Latitude is accurate to a fraction of an arc minute.
//!
//! \param jd_ut1 The Julian date in UT1.
//! \param o_orb The orbital elements of the observer.
//! \param q_orb The orbital elements of the object being observed.
//! \param phys The physical elements of the object being observed.
//! \param lat For returning the sub-observer latitude, in degrees.
//! \param lon For returning the sub-observer longitude, in degress.
//! \param dist For returning the distance to the object, in AU.  Set to NULL if
//!             you do not require this value.  Note that this returns the
//!             distance between the centres of the two bodies.
//------------------------------------------------------------------------------

void aes_subobs_point(double jd_ut1, const struct ae_orbit_t *o_orb, 
                      const struct ae_orbit_t *q_orb, 
                      const struct ae_physical_t *phys,
                      double *lat, double *lon, double *dist) {
  int i;
  double jd_tt, t_light, j[3], o[3], v_o[3], n[3], q[3], ra, dec, w, f;

  // Get the date in TT.
  jd_tt = jd_ut1 + ae_delta_t(jd_ut1) * AE_D_PER_S;

  // Get the unit vector pointing from the centre of the observed body to the
  // observer.  We could use ae_geocentric_from_orbit() here, but we also need
  // to know the light-travel time.  First, the observer coordinates.
  ae_kepler(jd_tt, o_orb, o);
  ae_v_orbit(jd_tt, o_orb, v_o);

  // Get the light-time corrected object coordinates.
  ae_kepler(jd_tt, q_orb, q);
  t_light = 0;
  for (i = 0; i < 2; i++) {
    t_light = ae_light_t(o, q, 0);
    ae_kepler(jd_tt - t_light, q_orb, q);
  }

  // The rest of the reduction is generic.
  ae_geocentric(jd_tt, q, o, v_o, &ra, &dec, dist);
  ae_polar_to_rect(ra, dec, -1.0, j);

  // Get the coordinates of the observed body's north pole, and prime meridean.
  // Account for light-travel time.
  ae_phys_pole(jd_ut1 - t_light, jd_tt, phys, n, &w);

  // Get the flattening.
  f = ae_flattening(phys);

  // Get the subobserver latitude and longitude.
  ae_subobs_point(j, n, w, f, ae_is_retrograde(phys), lat, lon);

  return;
}

//------------------------------------------------------------------------------
//! Calculate the sub-observer point.
//! The sub-observer point is the latitude and longitude on an observed body at
//! the observer is directly over-head.
//!
//! Note that all input values should be precessed and nutated to the epoch of
//! observation.
//!
//! WARNING:  the longitude computed by this function is, for some reason, only
//! accurate to tens of arc minutes when compared to HORIZONS (and I don't know
//! why. Latitude is accurate to a fraction of an arc minute.
//!
//! \param j The unit vector pointing from the centre of the observed body to 
//!          the observer.
//! \param n The unit vector of the north pole of the observed body, in
//!          equatorial coordinates of date.
//! \param w The prime meridean of the observed body, in degreees.
//! \param f The flattening of the observed body.
//! \param retrograde Pass non-zero value if the rotation is retrograde, i.e.,
//!                   the body rotates in the opposite sense of the earth's
//!                   rotation.  Pass 0 for prograde rotation.
//! \param lat For returning the sub-observer latitude, in degrees.
//! \param lon For returning the sub-observer longitude, in degrees.
//------------------------------------------------------------------------------

void ae_subobs_point(double j[], double n[], double w, double f, int retrograde,
                     double *lat, double *lon) {
  double u[3], y[3], p[3], pole_ra, pole_dec;
  double sin_g, cos_g, sin_s, cos_s, u_ra, u_dec;

  ae_rect_to_polar(n, &pole_ra, &pole_dec, NULL);
  
  sin_g = sin(w * AE_DTR) * cos(pole_dec * AE_DTR);
  cos_g = cos(asin(sin_g));
  cos_s = cos(w * AE_DTR) / cos_g;
  sin_s = sin(w * AE_DTR) * sin(pole_dec * AE_DTR) / cos_g;

  u_ra = pole_ra + 90.0 + ae_zatan2(cos_s, sin_s) * AE_RTD;
  u_dec = ae_zatan2(cos_g, sin_g) * AE_RTD;
  ae_polar_to_rect(u_ra, u_dec, 1.0, u);

  // Compute y = u x n.
  y[0] = u[1] * n[2] - u[2] * n[1];
  y[1] = u[2] * n[0] - u[0] * n[2];
  y[2] = u[0] * n[1] - u[1] * n[0];

  // Calculate the subearth point.
  p[0] = j[0] * u[0] + j[1] * u[1] + j[2] * u[2];
  p[1] = j[0] * y[0] + j[1] * y[1] + j[2] * y[2];
  p[2] = j[0] * n[0] + j[1] * n[1] + j[2] * n[2];

  // Find the longitude and latitude.  For some reason the longitude needs to be
  // flipped around here (i.e., the 360 - x) to agree with JPL---probably a 
  // minus sign missing somewhere above.
  *lat = ae_mod_180(ae_zatan2((1.0 - f) * (1.0 - f), tan(asin(p[2]))) * 
                    AE_RTD);
  if (retrograde)
    *lon = ae_mod_360(-ae_zatan2(p[0], p[1]) * AE_RTD);
  else
    *lon = ae_mod_360(ae_zatan2(p[0], p[1]) * AE_RTD);

  return;
}


//------------------------------------------------------------------------------
//! Calculate the axis of a body's north pole and its prime meridean.
//! Get the precessed, nutated direction of the body's north pole; also return
//! the prime meridean.
//!
//! \param jd_ut1 The Julian date in UT1.  For reasonable accuracy, pass a date
//!               that has been corrected for light travel time, to get the
//!               position of the north pole and meridean as seen by the
//!               observer.
//! \param jd_tt The Julian date in TT.  This should not be corrected for
//!              light-travel time.
//! \param phys The physical parameters of the body.
//! \param n For returning the unit vector of the north pole, in equatorial
//!          units of date, in degrees.
//! \param w For returning the prime meridean, in degrees.
//------------------------------------------------------------------------------

void ae_phys_pole(double jd_ut1, double jd_tt, const struct ae_physical_t *phys,
                  double n[], double *w) {
  int i;
  double d, t, ra_0, dec_0, w_0, x;

  // Calculate number of days and years since J2000.0.
  d = jd_ut1 - AE_J2000;
  t = d / 36525.0;

  // Calculate the ra, dec and w in J2000 coordinates.
  ra_0 = phys->pole_ra + phys->pole_ra_t * t;
  dec_0 = phys->pole_dec + phys->pole_dec_t * t;
  w_0 = phys->w + phys->w_d * d + phys->w_d_sq * d * d;
  for (i = 0; phys->ra_sin_term[i].time_var != AE_PHYSICAL_END; i++) {
    x = (phys->ra_sin_term[i].time_var == AE_PHYSICAL_D) ? d : t;
    ra_0 += phys->ra_sin_term[i].a * 
            sin(AE_DTR * phys->ra_sin_term[i].b +
                AE_DTR * phys->ra_sin_term[i].c * x +
                AE_DTR * phys->ra_sin_term[i].d * x * x);
  }
  for (i = 0; phys->dec_cos_term[i].time_var != AE_PHYSICAL_END; i++) {
    x = (phys->dec_cos_term[i].time_var == AE_PHYSICAL_D) ? d : t;
    dec_0 += phys->dec_cos_term[i].a * 
             cos(AE_DTR * phys->dec_cos_term[i].b +
                 AE_DTR * phys->dec_cos_term[i].c * x +
                 AE_DTR * phys->dec_cos_term[i].d * x * x);
  }
  for (i = 0; phys->w_sin_term[i].time_var != AE_PHYSICAL_END; i++) {
    x = (phys->w_sin_term[i].time_var == AE_PHYSICAL_D) ? d : t;
    w_0 += phys->w_sin_term[i].a * 
           sin(AE_DTR * phys->w_sin_term[i].b +
               AE_DTR * phys->w_sin_term[i].c * x +
               AE_DTR * phys->w_sin_term[i].d * x * x);
  }

  // Precess & nutate.
  ae_polar_to_rect(ra_0, dec_0, 1.0, n);
  ae_precess(jd_tt, n, AE_FROM_J2000);
  aes_nutate(jd_tt, n, AE_FROM_J2000);

  // Copy w.
  *w = ae_mod_360(w_0);

  return;
}


//------------------------------------------------------------------------------
//! Calculate a body's flattening.
//! Takes the equatorial radius, \f$a\f$, and polar radius, \f$c\f$, and
//! calculates the flattening of the oblate sphere:
//! \f[f \equiv \frac{a - b}{a}\f]
//!
//! \param phys The physical parameters of the body.
//!
//! \return The flattening.
//------------------------------------------------------------------------------

double ae_flattening(const struct ae_physical_t *phys) {
  return (phys->r_eq - phys->r_pole) / phys->r_eq;
}
