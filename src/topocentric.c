//! \file topocentric.c
//! Reduce from geocentric to topocentric coordinates.
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
//! Apply diurnal aberrations and calculate topocentric altitude and azimuth.
//! Apply diurnal aberrations and calculate topocentric altitude and azimuth, 
//! given the geocentric apparent right ascension and declination.  From AA 
//! page B60 and D3.
//!
//! The geocentric coordinates can be obtained from ae_geocentric_from_orbit()
//! or a similar routine.
//!
//! \param jd_tt The Julian date in TT.
//! \param jd_ut1 The Julian date in UT1 measure.
//! \param tlat The geocentric latitude of the observer, in degrees.
//! \param glat The geodetic latitude of the observer, in degrees.
//! \param tlong The longitude of the observer, in degrees.
//! \param trho The distance from the centre of the earth to the observer, in
//!             earth radii.
//! \param refrac A routine that calculates refraction (pass NULL to perform
//!               no refraction correction).  Its first argument is the
//!               altitude, in degrees.  Its second argument is a void pointer
//!               for passing other parameters, which will be particular to the
//!               refraction algorithm.  These parameters are passed via \p
//!               refrac_param.
//! \param refrac_param For passing auxiliary parameters to the \p refrac
//!                     function.  Pass NULL if \p refrac is not being used.
//! \param dist The geocentric distance to the object being observed, in AU.
//!             Pass 0 for extra-solar-system objects.
//! \param ra The right ascension of the object, in degrees; this routine
//!           returns the topocentric right ascension in this parameter.
//! \param dec The declination of the object, in degrees; this routine returns
//!            the topocentric declination in this parameter.
//------------------------------------------------------------------------------

void ae_topocentric(double jd_tt, double jd_ut1, double tlat, double glat, 
                    double tlong, double trho, double (*refrac)(double, void *),
                    void *refrac_param, double dist, double *ra, double *dec) {
  double last, alt, az, D, nutl, nuto, eps;

  // Get nutation and epsilon.
  eps = ae_epsilon(jd_tt);
  ae_nutation_lon_ob(jd_tt, &nutl, &nuto);

  // Local apparent sidereal time.
  last = ae_last(jd_ut1, jd_tt, tlong, nutl, eps);

  // Correct for diurnal aberration.
  ae_diurnal_aberration(last, tlat, trho, AE_ADD_ABERRATION, ra, dec);

  // Diurnal parallax.
  if (dist)
    ae_diurnal_parallax(last, tlat, trho, dist, ra, dec);

  // Convert to alt/az.
  ae_radec_to_altaz(last, glat, *ra, *dec, &alt, &az);

  // Correction for atmospheric refraction (unit = degrees).
  if (refrac != NULL)
    D = refrac(alt, refrac_param);
  else
    D = 0;
  alt += D;

  // Convert back to ra/dec.
  ae_altaz_to_radec(last, glat, alt, az, ra, dec);
}


//------------------------------------------------------------------------------
//! Apply diurnal aberrations and calculate topocentric altitude and azimuth.
//! This is a simplified version of ae_topocentric(), which automatically
//! calculates \p jd_tt, \p tlat and \p trho.  It assumes that the altitude 
//! of the observer is at sea level.  It also uses no atmospheric
//! model.  It might be less efficient to use this routine if speed is required.
//!
//! \param jd_ut1 The Julian date in UT1 measure.
//! \param glat The geodetic latitude of the observer, in degrees.
//! \param tlong The longitude of the observer, in degrees.
//! \param dist The geocentric distance to the object being observed, in AU.
//!             Pass 0 for extra-solar-system objects.
//! \param ra The right ascension of the object, in degrees; this routine
//!           returns the topocentric right ascension in this parameter.
//! \param dec The declination of the object, in degrees; this routine returns
//!            the topocentric declination in this parameter.
//------------------------------------------------------------------------------

void aes_topocentric(double jd_ut1, double glat, double tlong, double dist, 
                    double *ra, double *dec) {
  double jd_tt, tlat, trho;

  jd_tt = jd_ut1 + ae_delta_t(jd_ut1) / AE_S_PER_D;
  ae_geocentric_lat(glat, 0, &tlat, &trho);

  ae_topocentric(jd_tt, jd_ut1, tlat, glat, tlong, trho, NULL, NULL, dist, ra,
                 dec);
}
