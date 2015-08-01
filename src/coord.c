//! \file coord.c
//! This file contains functions for performing coordinate transformations.
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

#include <ctype.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "aephem.h"

//------------------------------------------------------------------------------
//! Convert an equatorial rectangular unit vector to ra/dec/radius.
//!
//! \param rect An equatorial rectangular vector.
//! \param ra For returning the right ascension, in degrees.
//! \param dec For returning the declination, in degrees.
//! \param radius For returning the radius.  Set to NULL if you do not need this
//!               value.
//------------------------------------------------------------------------------

void ae_rect_to_polar(const double rect[], double *ra, double *dec,
                      double *radius) {
  int i;
  double r;

  r = 0.0;
  for (i = 0; i < 3; i++)
    r += rect[i] * rect[i];
  r = sqrt(r);

  *ra = AE_RTD * ae_zatan2(rect[0], rect[1]);
  *dec = AE_RTD * asin(rect[2] / r);
  if (radius != NULL)
    *radius = r;
}



//------------------------------------------------------------------------------
//! Convert equatorial rectangular coordinates to polar coordinates.
//! This routine is similar to ae_rect_to_polar(), except that it also performs //! any necessary precession.
//!
//! \param pp The equatorial rectangular coordinates.
//! \param jd The Julian date.
//! \param ofdate If 1, precess from J2000.0 to \p jd.
//! \param polar For returning the polar coordinates.
//------------------------------------------------------------------------------

void ae_rect_to_polar_with_jd(double pp[], double jd, char ofdate,
                              double polar[]) {
  double s[3], x, y, z, yy, zz, r, eps, coseps, sineps;
  int i;

  // Make local copy of position vector and calculate radius.
  r = 0.0;
  for (i = 0; i < 3; i++) {
    x = pp[i];
    s[i] = x;
    r += x * x;
  }
  r = sqrt(r);

  // Precess to equinox of date jd, if requested.
  if (ofdate)
    ae_precess(jd, s, AE_FROM_J2000);

  // Convert from equatorial to ecliptic coordinates.
  eps = ae_epsilon(jd);
  coseps = cos(AE_STR * eps);
  sineps = sin(AE_STR * eps);
  yy = s[1];
  zz = s[2];
  x = s[0];
  y = coseps * yy + sineps * zz;
  z = -sineps * yy + coseps * zz;

  yy = ae_zatan2(x, y);
  zz = asin(z / r);

  polar[0] = yy;
  polar[1] = zz;
  polar[2] = r;

  return;
}


//------------------------------------------------------------------------------
//! Convert ra/dec/radius to a rectangular vector.
//!
//! \param ra  The right ascension, in degrees.
//! \param dec The declination, in degrees.
//! \param radius The radius; the units of \p radius will determine the units of
//!               the output.
//! \param rect For returning the rectangular vector.
//------------------------------------------------------------------------------

void ae_polar_to_rect(double ra, double dec, double radius, double rect[]) {
  double cosra, sinra, cosdec, sindec;

  cosra = cos(ra * AE_DTR);
  sinra = sin(ra * AE_DTR);
  cosdec = cos(dec * AE_DTR);
  sindec = sin(dec * AE_DTR);

  rect[0] = radius * cosdec * cosra;
  rect[1] = radius * cosdec * sinra;
  rect[2] = radius * sindec;

  return;
}


//------------------------------------------------------------------------------
//! Convert ra/dec to alt/az.
//!
//! \param last The local aparent sidereal time, in degrees.
//! \param glat The geodetic latitude, in degrees.
//! \param ra The right ascension, in degrees.
//! \param dec The declination, in degrees.
//! \param alt For returning the altitude, in degrees.
//! \param az For returning the azimuth, in degrees.
//------------------------------------------------------------------------------

void ae_radec_to_altaz(double last, double glat, double ra, double dec,
                       double *alt, double *az) {
  double lha, coslha, sinlha, cosdec, sindec, coslat, sinlat, N, D;

  cosdec = cos(AE_DTR * dec);
  sindec = sin(AE_DTR * dec);
  lha = AE_DTR * (last - ra);
  coslha = cos(lha);
  sinlha = sin(lha);

  // Use the geodetic latitude for altitude and azimuth.
  coslat = cos(AE_DTR * glat);
  sinlat = sin(AE_DTR * glat);

  N = -cosdec * sinlha;
  D =  sindec * coslat  -  cosdec * coslha * sinlat;
  *az = AE_RTD * ae_zatan2(D, N);
  *alt = sindec * sinlat  +  cosdec * coslha * coslat;
  *alt = AE_RTD * asin(*alt);

  return;
}

  
//------------------------------------------------------------------------------
//! Convert alt/az to ra/dec.
//!
//! \param last The local aparent sidereal time, in degrees.
//! \param glat The geodetic latitude, in degrees.
//! \param alt The altitude, in degrees.
//! \param az The azimuth, in degrees.
//! \param ra For returning the right ascension, in degrees.
//! \param dec For returning the declination, in degrees.
//------------------------------------------------------------------------------

void ae_altaz_to_radec(double last, double glat, double alt, double az, 
                       double *ra, double *dec) {
  double x, y, z, sinlha, coslha, sindec, lha;
  static double coslat = 2, sinlat = 2, glat_save = 1000;

  // Use the geodetic latitude for altitude and azimuth.  Don't calculate what
  // we don't need to. 
  if (glat_save != glat) {
    coslat = cos(AE_DTR * glat);
    sinlat = sin(AE_DTR * glat);
  }
  
  x = cos(AE_DTR * alt);
  y = sin(AE_DTR * alt);
  z = cos(AE_DTR * az);
  sinlha = -x * sin(AE_DTR * az);
  coslha = y * coslat - x * z * sinlat;
  sindec = y * sinlat + x * z * coslat;

  lha = ae_zatan2(coslha, sinlha) * AE_RTD;
  (*ra) = ae_mod_360(last - lha);
  (*dec) = asin(sindec) * AE_RTD;

  return;
}

//------------------------------------------------------------------------------
//! Convert ra/dec to l/b (galactic coordinates).
//!
//! \param ra The right ascension, in degrees.
//! \param dec The declination, in degrees.
//! \param l For returning the galactic longitude, in degrees.
//! \param b For returning the galactic latitude, in degrees.
//! \param fk5 If non-zero, assume FK5 system (i.e., J2000 compatible); if
//!            zero, use the old FK4 (B1950 compatible) galactic pole.
//------------------------------------------------------------------------------

void ae_radec_to_gal(double ra, double dec, double *l, double *b, int fk5) {
  double cos_b, x1, y1, z1, x2, y2, z2, r;

  ra *= AE_DTR;
  dec *= AE_DTR;

  cos_b = cos(dec);
  x1 = cos(ra) * cos_b;
  y1 = sin(ra) * cos_b;
  z1 = sin(dec);

  if (fk5) {
    x2 = -0.054875539726 * x1 - 0.873437108010 * y1 - 0.483834985808 * z1;
    y2 =  0.494109453312 * x1 - 0.444829589425 * y1 + 0.746982251810 * z1;
    z2 = -0.867666135858 * x1 - 0.198076386122 * y1 + 0.455983795705 * z1;
  }
  else {
    x2 = -0.066988739415 * x1 - 0.872755765852 * y1 - 0.483538914632 * z1;
    y2 =  0.492728466075 * x1 - 0.450346958020 * y1 + 0.744584633283 * z1;
    z2 = -0.867600811151 * x1 - 0.188374601723 * y1 + 0.460199784784 * z1;
  }

  if ((r = sqrt(x2 * x2 + y2 * y2)))
    *l = ae_mod_360(atan2(y2, x2) * AE_RTD);
  else
    *l = 0;
  if (z2)
    *b = ae_mod_360(atan2(z2, r) * AE_RTD);
  else
    *b = 0;

  return;
}


//------------------------------------------------------------------------------
//! Convert l/b (galactic coordinates) to ra/dec.
//!
//! \param l The galactic longitude, in degrees.
//! \param b The galactic latitude, in degrees.
//! \param ra For returning the right ascension, in degrees.
//! \param dec For returning the declination, in degrees.
//! \param fk5 If non-zero, assume FK5 system (i.e., J2000 compatible); if
//!            zero, use the old FK4 (B1950 compatible) galactic pole.
//------------------------------------------------------------------------------

void ae_gal_to_radec(double l, double b, double *ra, double *dec, int fk5) {
  double cos_b, x1, y1, z1, x2, y2, z2, r;

  l *= AE_DTR;
  b *= AE_DTR;

  cos_b = cos(b);
  x1 = cos(l) * cos_b;
  y1 = sin(l) * cos_b;
  z1 = sin(b);

  if (fk5) {
    x2 = -0.054875539692 * x1 + 0.494109453288 * y1 - 0.867666135842 * z1;
    y2 = -0.873437107998 * x1 - 0.444829589425 * y1 - 0.198076386097 * z1;
    z2 = -0.483834985832 * x1 + 0.746982251827 * y1 + 0.455983795747 * z1;
  }
  else {
    x2 = -0.066988739415 * x1 + 0.492728466076 * y1 - 0.867600811152 * z1;
    y2 = -0.872755765852 * x1 - 0.450346958020 * y1 - 0.188374601723 * z1;
    z2 = -0.483538914632 * x1 + 0.744584633283 * y1 + 0.460199784784 * z1;
  }

  if ((r = sqrt(x2 * x2 + y2 * y2)))
    *ra = ae_mod_360(atan2(y2, x2) * AE_RTD);
  else
    *ra = 0;
  if (z2)
    *dec = ae_mod_360(atan2(z2, r) * AE_RTD);
  else
    *dec = 0;

  return;
}


//------------------------------------------------------------------------------
//! Calculate the geocentric latitude and the distance to the earth's centre.
//! This function uses the reciprocal of flattening and the earth's radius (or,
//! more precisely, semi-major axis) from from the WGS84 definition.  
//! (See #AE_FLAT and #AE_R_EARTH.)
//!
//! \param glat The geodetic latitude of the observer, in degrees.
//! \param height The height of the observer above sea level, in metres.
//! \param tlat For returning the geocentric latitude, in degrees.
//! \param trho For returning the distance to the earth's centre, in earth
//!             radii.
//------------------------------------------------------------------------------
  
void ae_geocentric_lat(double glat, double height, double *tlat, double *trho) {
  double a, b, fl, co, si, u;
  
  u = glat * AE_DTR;
  co = cos(u);
  si = sin(u);
  fl = 1.0 - 1.0 / AE_FLAT;
  fl = fl * fl;
  si = si * si;
  u = 1.0 / sqrt(co * co + fl * si);
  a = AE_R_EARTH * u + height;
  b = AE_R_EARTH * fl * u + height;
  *trho = sqrt(a * a * co * co + b * b * si);
  *tlat = AE_RTD * acos(a * co / (*trho));
  
  if (glat < 0.0) 
    *tlat = -(*tlat);
  *trho /= AE_R_EARTH;

  return;
}


//------------------------------------------------------------------------------
//! Convert change in rectangular coordinatates to change in ra/dec.
//! For changes greater than about 0.1 degree, the coordinates are converted 
//! directly to ra and dec and the results subtracted.  For small changes, the 
//! change is calculated to first order by differentiating \f$tan(ra) = y/x\f$
//! to obtain \f[\frac{dra}{\cos^2(ra)} = \frac{dy}{x} - \frac{y dx}{x^2},\f]
//! where \f[\cos^2(ra) = \frac{1}{1 + (y/x)^2}.\f]
//!
//! The change in declination \f$arcsin(z/R)\f$ is 
//! \f[d\arcsin(u) = \frac{du}{\sqrt{1-u^2}},\f] where \f$u = z/R.\f$
//!
//! \param q0 The initial object - earth vector.
//! \param q1 The vector after motion or aberration.
//! \param dra The change in right ascension.
//! \param ddec The change in declination.
//------------------------------------------------------------------------------

void ae_delta_q(double q0[], double q1[], double *dra, double *ddec) {
  double dp[3], A, B, P, Q, x, y, z;
  int i;

  P = 0.0;
  Q = 0.0;
  z = 0.0;
  for (i = 0; i < 3; i++) {
    x = q0[i];
    y = q1[i];
    P += x * x;
    Q += y * y;
    y = y - x;
    dp[i] = y;
    z += y*y;
  }

  A = sqrt(P);
  B = sqrt(Q);

  if ((A < 1.e-7) || (B < 1.e-7) || (z/(P + Q)) > 5.e-7) {
    P = ae_zatan2(q0[0], q0[1]);
    Q = ae_zatan2(q1[0], q1[1]);
    Q = Q - P;
    while (Q < -M_PI)
      Q += 2.0 * M_PI;
    while (Q > M_PI)
      Q -= 2.0 * M_PI;
    *dra = Q;
    P = asin(q0[2] / A);
    Q = asin(q1[2] / B);
    *ddec = Q - P;

    return;
  }

  x = q0[0];
  y = q0[1];
  if (x == 0.0)
    *dra = 1.0e38;
  else {
    Q = y / x;
    Q = (dp[1]  -  dp[0] * y / x)/(x * (1.0 + Q * Q));
    *dra = Q;
  }

  x = q0[2] / A;
  P = sqrt(1.0 - x * x);
  *ddec = (q1[2] / B - x) / P;

  return;
}
