//! \file geocentric.c
//! Reduce heliocentric position to geocentric right ascension and declination.
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
#include <stdio.h>
#include <string.h>

#include "aephem.h"

/*
//! This function optionally calculates the apparent visual magnitude (see
//! parameter \p app_mag.  The magnitude is
//!   \f[V(1,0) + 2.5 \log_{10} \left( \frac{SE^2 SO^2}{\Phi} \right) ,\f]
//! where \f$SE\f$ is the sun-earth distance, \f$SO\f$ is the sun-object
//! distance, \f$\Phi\f$ is the phase and V(1,0) is the magnitude at 1 AU from
//! both earth and sun and 100% illumination.  The phase is the geometric
//! fraction of the disc illuminated and is given by
//!   \f[\Phi = \frac{1 + pq}{2},\f]
//! where \f$pq = \cos(sun-object-earth\;\;angle)\f$.  However, this routine
//! `fudges' the phase for light leakage in magnitude estimation and uses
//!   \f[\Phi_{fudge} = \frac{1.01 + 0.99 pq}{2}.\f]
//! Note this phase term estimate does not reflect reality well (no pun 
//! intended).  In particular, calculated magnitudes of Mercury and Venus are 
//! inaccurate.
//!
  // Calculate visual magnitude and/or phase.
  if (app_mag != NULL) {
    b = 0.5 * (1.01 + 0.99 * pq);
    *app_mag = el->mag + 2.1715 * log(EO * SO) - 1.085 * log(b);
  }
  if (phase != NULL)
    *phase = 0.5 * (1.0 + pq);
*/


//------------------------------------------------------------------------------
//! Given heliocentric coordinates, reduce to geocentric coordinates.
//! This routine reduces heliocentric coordinates to equitorial coordinates.
//! Unlike ae_geocentric_from_orbit() or ae_geocentric_from_jpl(), it does not
//! automatically compute the heliocentric coordinates.  This means that
//! light travel-time is <em>not</em> automatically accounted for!  The input
//! position of the planet is what is used.
//!
//! Corrections done during the reduction are:  light deflection, annual
//! aberration, precession and nutation.
//!
//! This routine should be thought of as a macro.  One can call the functions 
//! that perform the corrections listed above individually.  This routine simply
//! wraps them into one call.
//!
//! \param jd_tt The Julian date in TT.
//! \param q The heliocentric rectangular coordinates of the object being
//!          observed, in AU.
//! \param e The heliocentric rectangular coordinates of the earth, in AU.
//! \param v_e The velocity of the earth, in AU per day.
//! \param ra For returning the right ascension, in degrees.
//! \param dec For returning the declination, in degrees.
//! \param dist For returning the geocentric distance between earth and the
//!             object, in AU.  Set to NULL if you do not require this value.
//------------------------------------------------------------------------------

void ae_geocentric(double jd_tt, double q[], double e[], double v_e[], 
                   double *ra, double *dec, double *dist) {
  int i;
  double p[3], q_e_norm;
  
  // Find vector between object and earth; also calculate the norm of this
  // vector.
  for (i = 0, q_e_norm = 0; i < 3; i++) {
    p[i] = q[i] - e[i];
    q_e_norm += p[i] * p[i];
  }
  q_e_norm = sqrt(q_e_norm);

  // Reduce to a unit vector.
  for (i = 0; i < 3; i++)
    p[i] /= q_e_norm;

  // Correct position for light deflection.
  ae_relativity(p, q, e);

  // Correct for aberration.
  ae_annual_aberration(v_e, p, AE_ADD_ABERRATION);

  // Precession of the equinox and ecliptic from J2000.0 to ephemeris date.
  ae_precess(jd_tt, p, AE_FROM_J2000);

  // Nutation.
  aes_nutate(jd_tt, p, AE_FROM_J2000);

  // Convert to ra/dec.
  ae_rect_to_polar(p, ra, dec, NULL);
  if (dist != NULL)
    *dist = q_e_norm;
}


//------------------------------------------------------------------------------
//! Determine geocentric coordinates using orbital elements.
//! This routine reduces the heliocentric equitorial coordinates,
//! obtained by calling ae_kepler(), to the geocentric right ascension and
//! declination.  The position is corrected for light travel time, light
//! deflection, annual aberration, precession and nutation.
//!
//! This routine should be thought of as a macro.  One can call ae_kepler() and
//! the functions that perform the corrections listed above individually.  This
//! routine simply wraps them into one call.
//!
//! To obtain geocentric coordinates from a JPL ephemeris file, use
//! ae_geocentric_from_jpl().
//!
//! \param jd_tt The Julian date in TT.
//! \param o_orb Orbital elements of the observer.
//! \param q_orb The orbital elements of the object.
//! \param ra For returning the geocentric right ascension, in degrees.
//! \param dec For returning the geocentric declination, in degrees.
//! \param dist For returning the geocentric distance between earth and the
//!             object, in AU.  Set to NULL if you do not require this value.
//------------------------------------------------------------------------------

void ae_geocentric_from_orbit(double jd_tt, const struct ae_orbit_t *o_orb,
                              const struct ae_orbit_t *q_orb, double *ra,
                              double *dec, double *dist) {
  double q[3], o[3], v_o[3];

  // Get the position of the observer, and the light-corrected position of the
  // object.
  ae_kepler(jd_tt, o_orb, o);
  ae_light_t_orbit(jd_tt, o, q_orb, q);

  // Get the velocity of the earth.
  ae_v_orbit(jd_tt, o_orb, v_o);

  // The rest of the reduction is very generic.
  ae_geocentric(jd_tt, q, o, v_o, ra, dec, dist);

  return;
}


//------------------------------------------------------------------------------
//! Determine geocentric coordinates using JPL ephemerides.
//! This routine reduces the heliocentric equitorial coordinates,
//! obtained by calling ae_jpl_get_coords(), to the geocentric right ascension 
//! and declination.  The position is corrected for light travel time, light
//! deflection, annual aberration, precession and nutation.
//!
//! It is assumed that the JPL planetary ephemeris file is being used,
//! since it automatically polls the file for Earth's position and velocity.  If
//! you are using a different ephemeris file, you cannot use this routine.
//!
//! This routine should be thought of as a macro.  One can call
//! ae_jpl_get_coords() and the functions that perform the corrections listed 
//! above individually.  This routine simply wraps them into one call.
//!
//! To obtain geocentric coordinates from an orbital elements struct, use
//! ae_geocentric_from_orbit().
//!
//! \param jh A handle initialised by ae_jpl_init().
//! \param jd_tt The Julian date in TT.
//! \param obj_num The object for which to get coordinates.  If a planetary 
//!                ephemeris file is being used, one can use the predefined
//!                #ae_ss_bodies_t constants for clarity.
//! \param ra For returning the geocentric right ascension, in degrees.
//! \param dec For returning the geocentric declination, in degrees.
//! \param dist For returning the geocentric distance between earth and the
//!             object, in AU.  Set to NULL if you do not require this value.
//!
//! \return 0 on success; on failure, a negative error code from #ae_retcode_t.
//------------------------------------------------------------------------------

int ae_geocentric_from_jpl(struct ae_jpl_handle_t *jh, double jd_tt, 
                           int obj_num, double *ra, double *dec, double *dist) {
  int i;
  double q[3], e[3], v_e[3];

  // Get the position and velocity of the earth, and the light-corrected 
  // position of the object.
  if ((i = ae_jpl_get_coords(jh, jd_tt, AE_SS_EARTH, e, v_e, 1)))
    return i;
  if ((i = ae_light_t_jpl(jh, jd_tt, e, obj_num, q, NULL, 1)))
    return i;

  // The rest of the reduction is very generic.
  ae_geocentric(jd_tt, q, e, v_e, ra, dec, dist);

  return 0;
}


//------------------------------------------------------------------------------
//! Reduce a star's catalogue coordinates to its apparent place.
//! This routine puts the coordinates in J2000.0 (if not already there), and
//! performs precession, proper motion, parallax, light deflection, aberration 
//! and nutation.
//!
//! This routine should be thought of as a macro.  One can call the functions 
//! that perform the corrections listed above individually.  This routine 
//! simply wraps them into one call.
//!
//! \param jd_tt The Julian date in TT.
//! \param e The heliocentric rectangular coordinates of the earth, in AU.  
//!          These can be obtained from ae_kepler() or ae_jpl_get_coords().
//! \param v_e The velocity of the earth, in AU per day.  This can be obtained
//!            from ae_v_orbit() or ae_jpl_get_coords().
//! \param star The catalogue object.  If it is in the FK4 B1950.0 coordinates,
//!             this routine will convert it to FK5 J2000.0 coordinates.
//! \param ra For returning the apparent (geocentric) right ascension, in
//!           degrees.
//! \param dec For returning the apparent (geocentric) declination, in degrees.
//------------------------------------------------------------------------------

void ae_geocentric_from_cat(double jd_tt, double e[], double v_e[],
                            const struct ae_star_t *star, double *ra,
                            double *dec) {
  int i;
  double p[3], q[3], m[3], temp[3], T, vpi, r_earth[3],  cosdec, sindec, cosra;
  double sinra, norm;
  struct ae_star_t ls;

  // Make a local copy, since we might need to modify the original.
  memcpy(&ls, star, sizeof(struct ae_star_t));

  // Convert to FK5, if necessary.
  if (ls.epoch == AE_B1950)
    ae_fk4_to_fk5(&ls);

  // Convert to rectangular coordinates.
  cosdec = cos(ls.dec * AE_DTR);
  sindec = sin(ls.dec * AE_DTR);
  cosra = cos(ls.ra * AE_DTR);
  sinra = sin(ls.ra * AE_DTR);
  q[0] = cosra * cosdec;
  q[1] = sinra * cosdec;
  q[2] = sindec;
  vpi = 21.094952663 * ls.v * ls.px;// * AE_STR;
  m[0] = -ls.mura * cosdec * sinra - ls.mudec * sindec * cosra + vpi * q[0];
  m[1] = ls.mura * cosdec * cosra - ls.mudec * sindec * sinra + vpi * q[1];
  m[2] = ls.mudec * cosdec + vpi * q[2];

  // Calculate the current position of the earth.
  for (i = 0; i < 3; i++)
    r_earth[i] = e[i];

  // Precess the earth to the star epoch.
  ae_precess(ls.epoch, e, AE_FROM_J2000);

  // Correct for proper motion and parallax.
  T = (jd_tt - ls.epoch) / 36525.0;
  for (i = 0; i < 3; i++)
    p[i] = q[i] + (T * m[i] - ls.px * e[i]) * AE_STR;

  // Precess the star to J2000.
  ae_precess(ls.epoch, p, AE_TO_J2000);

  // Reset the earth to J2000.
  for (i = 0; i < 3; i++)
    e[i] = r_earth[i];

  // Find unit vector from earth in direction of object.
  norm = sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
  for (i = 0; i < 3; i++) {
    p[i] /= norm;
    temp[i] = p[i];
  }

  // Correct position for light deflection.
  ae_relativity(p, p, e);

  // Correct for annual aberration.
  ae_annual_aberration(v_e, p, AE_ADD_ABERRATION);

  // Precession of the equinox and ecliptic from J2000.0 to ephemeris date.
  ae_precess(jd_tt, p, AE_FROM_J2000);

  // Adjust for nutation at current ecliptic.
  aes_nutate(jd_tt, p, AE_FROM_J2000);

  // Update the star object.
  ae_rect_to_polar(p, ra, dec, NULL);

  return;
}

//------------------------------------------------------------------------------
//! Calculate the geocentric coordinates of the sun from orbital elements.
//! Corrects for precession, nutation, annual aberration.
//!
//! Note that if using JPL ephemeris, one can use ae_geocentric_from_jpl() to do
//! this calculation since the sun is an object in those ephemerides.
//!
//! \param jd_tt The Julian date in TT.
//! \param o_orb Orbital elements of the observer.
//! \param ra For returning the geocentric right ascension, in degrees.
//! \param dec For returning the geocentric declination, in degrees.
//! \param dist For returning the geocentric distance between earth and the
//!             object, in AU.  Set to NULL if you do not require this value.
//------------------------------------------------------------------------------

void ae_geocentric_sun_from_orbit(double jd_tt, const struct ae_orbit_t *o_orb,
                                  double *ra, double *dec, double *dist) {

  double q[3], o[3], v_o[3];
  
  // Get the position of the observer and its velocity.
  ae_kepler(jd_tt, o_orb, o);
  ae_v_orbit(jd_tt, o_orb, v_o);

  // The position is at the origin in heliocentric coordinates.
  q[0] = 0;
  q[1] = 0;
  q[2] = 0;

  // The rest of the reduction is very generic.
  ae_geocentric(jd_tt, q, o, v_o, ra, dec, dist);

  return;
}

//------------------------------------------------------------------------------
//! Calculate the geocentric coordinates of the moon from orbital elements.
//! Corrects for precession, nutation, annual aberration.  Neglects light travel
//! time.
//!
//! Note that if using JPL ephemeris, one can use ae_geocentric_from_jpl() to do
//! this calculation since the moon is an object in those ephemerides.
//!
//! The Astronomical Almanac (Section D, Daily Polynomial Coefficients)
//! seems to omit the annual aberration at one point, even though the reference
//! ephemeris is inertial . . .
//!
//! \param jd_tt The Julian date in TT.
//! \param o_orb Orbital elements of the observer.
//! \param ra For returning the geocentric right ascension, in degrees.
//! \param dec For returning the geocentric declination, in degrees.
//! \param dist For returning the geocentric distance between earth and the
//!             object, in AU.  Set to NULL if you do not require this value.
//------------------------------------------------------------------------------
// This routine used to do more . . . maybe I'll add it back in at some point.
// -AH
//
// The moon's elongation from the sun and illumination factor can be optionally
// returned.
//
// The phase of the moon can be optionally returned.  It is calculated by 
// comparing the moon's longitude with earth's longitude.  The number of days 
// before or past indicated phase is estimated by assuming the true longitudes 
// change linearly with time.  These rates are estimated for the date, but do 
// not stay constant.  The error can exceed 0.15 day in 4 days.
//
// \param elong For returning the moon's elongation from the sun.  Set to NULL
//              if you do not need this value.
// \param illum For returning the moon's illumination fraction.  Set to NULL
//              if you do not need this value.
// \param phase For returning the moon's phase.  The buffer should be able to
//              hold at least 28 characters.  By a strange coincidence, this is
//              also the orbital period of the moon, in days.  Set to NULL if 
//              you do not need this value.
//------------------------------------------------------------------------------

void ae_geocentric_moon_from_orbit(double jd_tt, const struct ae_orbit_t *o_orb,
                                   double *ra, double *dec, double *dist) {
  double cosB, sinB, cosL, sinL, e[3];
  double rect[3], pol[3], eps, sineps, coseps, nutl, nuto;
  //int i, x, y, z;
  //double pp[3], qq[3], pe[3], re[3], rearth[3], ep, pq, v_e[3];
  //double yy, qq_norm, re_norm, rect_norm, pp_norm;

  // Geometric equatorial coordinates of the earth.
  //for (i = 0; i < 3; i++) {
  //  rearth[i] = e[i];
  //  re[i] = rearth[i];
  //}

  // Get the geocentric coordinates of the observer.
  ae_kepler(jd_tt, o_orb, e);
  //ae_v_orbit(jd_tt, ae_orb_earth, e_vel);

  // Get geometric coordinates of the Moon.
  ae_gmoon(jd_tt, rect, pol);

  // Light time correction to longitude, about 0.7".
  pol[0] -= 0.0118 * AE_DTR * AE_R_EARTH_AU / pol[2];

  // Compute obliquity of the ecliptic.
  eps = ae_epsilon(jd_tt);
  coseps = cos(AE_STR * eps);
  sineps = sin(AE_STR * eps);

  // Convert to equatorial system of date.
  cosB = cos(pol[1]);
  sinB = sin(pol[1]);
  cosL = cos(pol[0]);
  sinL = sin(pol[0]);
  rect[0] = cosB * cosL;
  rect[1] = coseps * cosB * sinL - sineps * sinB;
  rect[2] = sineps * cosB * sinL + coseps * sinB;

  /*
  // Rotate to J2000.
  ae_precess(jd_tt, rect, AE_TO_J2000);

  // Find Euclidean vectors and angles between earth, object, and the sun.
  pp_norm = 0;
  for (i = 0; i < 3; i++) {
    pp[i] = rect[i] * pol[2];
    qq[i] = rearth[i] + pp[i];
    pp_norm += pp[i] * pp[i];
  }
  pp_norm = sqrt(pp_norm);

  // Aberration of light.  The Astronomical Almanac (Section D, Daily Polynomial
  // Coefficients) seems to omit this, even though the reference ephemeris is 
  // inertial.

  // Precess to date.
  ae_precess(jd_tt, rect, AE_FROM_J2000); */

  // Correct for nutation at date jd_tt.
  ae_nutation_lon_ob(jd_tt, &nutl, &nuto);
  ae_nutate(nutl, nuto, eps, rect, AE_FROM_J2000);

  // Apparent geocentric right ascension and declination.
  ae_rect_to_polar(rect, ra, dec, dist);
  if (dist != NULL)
    *dist = pol[2];

  return;

  /* PERHAPS LATER THESE PARTS CAN BE RESTORED SOMEHOW.  FOR NOW I WILL REMOVE
   * THEM FROM THE FUNCTIONALITY OF THIS FUNCTION.  THE FUNCTION CURRENTLY ENDS
   * HERE. -ADH
   *
  // For apparent ecliptic coordinates, rotate from the true equator into the 
  // ecliptic of date.
  cosL = cos(AE_STR * (eps + nuto));
  sinL = sin(AE_STR * (eps + nuto));
  y = cosL * rect[1] + sinL * rect[2];
  z = -sinL * rect[1] + cosL * rect[2];
  pol[0] = ae_zatan2( rect[0], y);
  pol[1] = asin(z);

  // Restore earth-moon distance.
  for (i = 0; i < 3; i++)
    rect[i] *= pp_norm;

  // Get apparent coordinates for the earth.
  z = re[0] * re[0] + re[1] * re[1] + re[2] * re[2];
  z = sqrt(z);
  for (i = 0; i < 3; i++)
    re[i] /= z;

  // Aberration of light.
  ae_annual_aberration(v_e, re, AE_ADD_ABERRATION);

  // Precession.
  ae_precess(jd_tt, re, AE_FROM_J2000);
  ae_nutate(jd_tt, nutl, nuto, re, AE_FROM_J2000);

  for (i = 0; i < 3; i++)
    re[i] *= z;
  ae_rect_to_polar_with_jd(re, jd_tt, 0, pe);

  // Find sun-moon-earth angles.
  qq_norm = 0;
  rect_norm = 0;
  re_norm = 0;
  pq = 0;
  ep = 0;
  for (i = 0; i < 3; i++) {
    qq[i] = re[i] + rect[i];
    pq = rect[i] * qq[i];
    ep = rect[i] * re[i];
    qq_norm += qq[i] * qq[i];
    rect_norm += rect[i] * rect[i];
    re_norm += re[i] * re[i];
  }
  qq_norm = sqrt(qq_norm);
  rect_norm = sqrt(rect_norm);
  re_norm = sqrt(re_norm);
  pq /= rect_norm * qq_norm;
  ep /= rect_norm * re_norm;

  if (elong != NULL)
    *elong = AE_RTD * acos(-ep);
  if (illum != NULL)
    *illum = 0.5 * (1.0 + pq);

  // Find phase of the Moon by comparing Moon's longitude with Earth's 
  // longitude.  The number of days before or past indicated phase is estimated
  // by assuming the true longitudes change linearly with time.  These rates are
  // estimated for the date, but do not stay constant.  The error can exceed 
  // 0.15 day in 4 days.
  if (phase != NULL) {
    x = pol[0] - pe[0];
    x = ae_mod_2pi(x ) * AE_RTD; // Difference in longitude.
    i = (int)(x / 90);          // Number of quarters.
    x = (x - i * 90.0);         // Phase angle mod 90 degrees.

    // Days per degree of phase angle.
    z = pol[2] / (12.3685 * 0.00257357);

    if( x > 45.0 ) {
      yy = -(x - 90.0) * z;
      if (yy > 1.0)
        sprintf(phase, "%.1f days before ", yy);
      else
        sprintf(phase, "%.2f days before ", yy);
      i = (i + 1) & 3;
    }
    else {
      yy = x * z;
      if (yy > 1.0)
        sprintf(phase, "%.1f days past ", yy);
      else
        sprintf(phase, "%.2f days past ", yy);
    }

    switch(i) {
      case 0:
        strcat(phase, "full moon");
        break;
      case 1:
        strcat(phase,  "third quarter");
        break;
      case 2:
        strcat(phase,  "new moon");
        break;
      case 3:
        strcat(phase,  "first quarter");
        break;
    }
  }*/
}
