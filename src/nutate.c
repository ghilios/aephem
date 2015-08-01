//! \file nutate.c
//! Correct for nutation.
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

void ae_nut_sscc(int k, double arg, int n);

//! Nutation parameters.
//! Each term in the expansion has a trigonometric argument given by 
//! \f[W = i MM + j MS + k FF + l DD + m OM\f]
//! where the variables are defined below.  The nutation in longitude is a sum 
//! of terms of the form \f$(a + bT) \sin(W)\f$. The terms for nutation in 
//! obliquity are of the form \f$(c + dT) \cos(W)\f$.  The coefficients are 
//! arranged in the tabulation as follows:
//!
//! <table>
//!   <tr><td>i </td><td> j </td><td> k </td><td> l </td><td> m </td><td> 
//!           a </td><td> b </td><td> c </td><td> d </td></tr>
//!   <tr><td> 0 </td><td> 0 </td><td> 0 </td><td> 0 </td><td> 1 </td><td>
//!            -171996 </td><td> -1742 </td><td> 92025 </td><td> 89 </td></tr>
//! </table>
//!
//! The first line of the table, above, is done separately since two of the 
//! values do not fit into 16 bit integers.  The values a and c are arc seconds
//! times 10000.  B and d are arc seconds per Julian century times 100000.  I 
//! through m are integers.  See the program for interpretation of MM, MS, etc.,
//! which are mean orbital elements of the Sun and Moon.
//!
//! If terms with coefficient less than X are omitted, the peak errors will be:
//!
//! <table>
//!   <tr><td> omit a &lt; </td><td> error, longitude </td><td>
//!            omit c &lt; </td><td> error, obliquity </td></tr>
//!   <tr><td>.0005"</td><td>.0100"</td><td>.0008" </td><td> .0094" </td></tr>
//!   <tr><td>.0046 </td><td> .0492	 </td><td> .0095 </td><td> .0481 </td></tr>
//!   <tr><td>.0123 </td><td> .0880	 </td><td> .0224 </td><td> .0905 </td></tr>
//!   <tr><td>.0386 </td><td> .1808	 </td><td> .0895 </td><td> .1129 </td></tr>
//! </table>
short ae_nt_param[105*9] = {
 0, 0, 0, 0, 2, 2062, 2,-895, 5,
-2, 0, 2, 0, 1, 46, 0,-24, 0,
 2, 0,-2, 0, 0, 11, 0, 0, 0,
-2, 0, 2, 0, 2,-3, 0, 1, 0,
 1,-1, 0,-1, 0,-3, 0, 0, 0,
 0,-2, 2,-2, 1,-2, 0, 1, 0,
 2, 0,-2, 0, 1, 1, 0, 0, 0,
 0, 0, 2,-2, 2,-13187,-16, 5736,-31,
 0, 1, 0, 0, 0, 1426,-34, 54,-1,
 0, 1, 2,-2, 2,-517, 12, 224,-6,
 0,-1, 2,-2, 2, 217,-5,-95, 3,
 0, 0, 2,-2, 1, 129, 1,-70, 0,
 2, 0, 0,-2, 0, 48, 0, 1, 0,
 0, 0, 2,-2, 0,-22, 0, 0, 0,
 0, 2, 0, 0, 0, 17,-1, 0, 0,
 0, 1, 0, 0, 1,-15, 0, 9, 0,
 0, 2, 2,-2, 2,-16, 1, 7, 0,
 0,-1, 0, 0, 1,-12, 0, 6, 0,
-2, 0, 0, 2, 1,-6, 0, 3, 0,
 0,-1, 2,-2, 1,-5, 0, 3, 0,
 2, 0, 0,-2, 1, 4, 0,-2, 0,
 0, 1, 2,-2, 1, 4, 0,-2, 0,
 1, 0, 0,-1, 0,-4, 0, 0, 0,
 2, 1, 0,-2, 0, 1, 0, 0, 0,
 0, 0,-2, 2, 1, 1, 0, 0, 0,
 0, 1,-2, 2, 0,-1, 0, 0, 0,
 0, 1, 0, 0, 2, 1, 0, 0, 0,
-1, 0, 0, 1, 1, 1, 0, 0, 0,
 0, 1, 2,-2, 0,-1, 0, 0, 0,
 0, 0, 2, 0, 2,-2274,-2, 977,-5,
 1, 0, 0, 0, 0, 712, 1,-7, 0,
 0, 0, 2, 0, 1,-386,-4, 200, 0,
 1, 0, 2, 0, 2,-301, 0, 129,-1,
 1, 0, 0,-2, 0,-158, 0,-1, 0,
-1, 0, 2, 0, 2, 123, 0,-53, 0,
 0, 0, 0, 2, 0, 63, 0,-2, 0,
 1, 0, 0, 0, 1, 63, 1,-33, 0,
-1, 0, 0, 0, 1,-58,-1, 32, 0,
-1, 0, 2, 2, 2,-59, 0, 26, 0,
 1, 0, 2, 0, 1,-51, 0, 27, 0,
 0, 0, 2, 2, 2,-38, 0, 16, 0,
 2, 0, 0, 0, 0, 29, 0,-1, 0,
 1, 0, 2,-2, 2, 29, 0,-12, 0,
 2, 0, 2, 0, 2,-31, 0, 13, 0,
 0, 0, 2, 0, 0, 26, 0,-1, 0,
-1, 0, 2, 0, 1, 21, 0,-10, 0,
-1, 0, 0, 2, 1, 16, 0,-8, 0,
 1, 0, 0,-2, 1,-13, 0, 7, 0,
-1, 0, 2, 2, 1,-10, 0, 5, 0,
 1, 1, 0,-2, 0,-7, 0, 0, 0,
 0, 1, 2, 0, 2, 7, 0,-3, 0,
 0,-1, 2, 0, 2,-7, 0, 3, 0,
 1, 0, 2, 2, 2,-8, 0, 3, 0,
 1, 0, 0, 2, 0, 6, 0, 0, 0,
 2, 0, 2,-2, 2, 6, 0,-3, 0,
 0, 0, 0, 2, 1,-6, 0, 3, 0,
 0, 0, 2, 2, 1,-7, 0, 3, 0,
 1, 0, 2,-2, 1, 6, 0,-3, 0,
 0, 0, 0,-2, 1,-5, 0, 3, 0,
 1,-1, 0, 0, 0, 5, 0, 0, 0,
 2, 0, 2, 0, 1,-5, 0, 3, 0, 
 0, 1, 0,-2, 0,-4, 0, 0, 0,
 1, 0,-2, 0, 0, 4, 0, 0, 0,
 0, 0, 0, 1, 0,-4, 0, 0, 0,
 1, 1, 0, 0, 0,-3, 0, 0, 0,
 1, 0, 2, 0, 0, 3, 0, 0, 0,
 1,-1, 2, 0, 2,-3, 0, 1, 0,
-1,-1, 2, 2, 2,-3, 0, 1, 0,
-2, 0, 0, 0, 1,-2, 0, 1, 0,
 3, 0, 2, 0, 2,-3, 0, 1, 0,
 0,-1, 2, 2, 2,-3, 0, 1, 0,
 1, 1, 2, 0, 2, 2, 0,-1, 0,
-1, 0, 2,-2, 1,-2, 0, 1, 0,
 2, 0, 0, 0, 1, 2, 0,-1, 0,
 1, 0, 0, 0, 2,-2, 0, 1, 0,
 3, 0, 0, 0, 0, 2, 0, 0, 0,
 0, 0, 2, 1, 2, 2, 0,-1, 0,
-1, 0, 0, 0, 2, 1, 0,-1, 0,
 1, 0, 0,-4, 0,-1, 0, 0, 0,
-2, 0, 2, 2, 2, 1, 0,-1, 0,
-1, 0, 2, 4, 2,-2, 0, 1, 0,
 2, 0, 0,-4, 0,-1, 0, 0, 0,
 1, 1, 2,-2, 2, 1, 0,-1, 0,
 1, 0, 2, 2, 1,-1, 0, 1, 0,
-2, 0, 2, 4, 2,-1, 0, 1, 0,
-1, 0, 4, 0, 2, 1, 0, 0, 0,
 1,-1, 0,-2, 0, 1, 0, 0, 0,
 2, 0, 2,-2, 1, 1, 0,-1, 0,
 2, 0, 2, 2, 2,-1, 0, 0, 0,
 1, 0, 0, 2, 1,-1, 0, 0, 0,
 0, 0, 4,-2, 2, 1, 0, 0, 0,
 3, 0, 2,-2, 2, 1, 0, 0, 0,
 1, 0, 2,-2, 0,-1, 0, 0, 0,
 0, 1, 2, 0, 1, 1, 0, 0, 0,
-1,-1, 0, 2, 1, 1, 0, 0, 0,
 0, 0,-2, 0, 1,-1, 0, 0, 0,
 0, 0, 2,-1, 2,-1, 0, 0, 0,
 0, 1, 0, 2, 0,-1, 0, 0, 0,
 1, 0,-2,-2, 0,-1, 0, 0, 0,
 0,-1, 2, 0, 1,-1, 0, 0, 0,
 1, 1, 0,-2, 1,-1, 0, 0, 0,
 1, 0,-2, 2, 0,-1, 0, 0, 0,
 2, 0, 0, 2, 0, 1, 0, 0, 0,
 0, 0, 2, 4, 2,-1, 0, 0, 0,
 0, 1, 0, 1, 0, 1, 0, 0, 0,
};

double ae_nut_ss[5][8];    //!< An array to hold the sine of multiple angles. 
double ae_nut_cc[5][8];    //!< An array to hold the cosine of multiple angles.
//! A handy macro.
#define ae_mod_3600(x) ((x) - 1296000. * floor ((x)/1296000.))

//------------------------------------------------------------------------------
//! Calculate the nutation in longitude and oblation.
//! References:
//! - "Summary of 1980 IAU Theory of Nutation (Final Report of the IAU Working
//!   Group on Nutation)", P. K. Seidelmann et al., in Transactions of the IAU
//!   Vol. XVIII A, Reports on Astronomy, P. A. Wayman, ed.; D. Reidel Pub. 
//!   Co., 1982.
//! - "Nutation and the Earth's Rotation", I.A.U. Symposium No. 78, May, 1977, 
//!   page 256. I.A.U., 1980.
//! - Woolard, E.W., "A redevelopment of the theory of nutation", The 
//!   Astronomical Journal, 58, 1-3 (1953).
//!
//! This program implements all of the 1980 IAU nutation series.  Results 
//! checked at 100 points against the 1986 AA; all agreed.
//!
//! \param jd_tt The Julian date in TT.
//! \param nutl For returning the nutation in longitude, in seconds of arc.
//! \param nuto For returning the nutation in oblation, in seconds of arc.
//------------------------------------------------------------------------------

void ae_nutation_lon_ob(double jd_tt, double *nutl, double *nuto) {
  double f, g, T, T2, T10;
  double MM, MS, FF, DD, OM;
  double cu, su, cv, sv, sw;
  double C, D;
  int i, j, k, k1, m;
  short *p;

  // Julian centuries from 2000 January 1.5, barycentric dynamical time.
  T = (jd_tt - 2451545.0) / 36525.0;
  T2 = T * T;
  T10 = T / 10.0;

  // Fundamental arguments in the FK5 reference system.
  // Longitude of the mean ascending node of the lunar orbit on the ecliptic, 
  // measured from the mean equinox of date.
  OM = (ae_mod_3600(-6962890.539 * T + 450160.280) + (0.008 * T + 7.455) * 
        T2) * AE_STR;

  // Mean longitude of the Sun minus the mean longitude of the Sun's perigee.
  MS = (ae_mod_3600(129596581.224 * T + 1287099.804) - (0.012 * T + 0.577) * 
        T2) * AE_STR;

  // Mean longitude of the Moon minus the mean longitude of the Moon's perigee.
  MM = (ae_mod_3600(1717915922.633 * T + 485866.733) + (0.064 * T + 31.310) * 
        T2) * AE_STR;

  // Mean longitude of the Moon minus the mean longitude of the Moon's node.
  FF = (ae_mod_3600(1739527263.137 * T + 335778.877) + (0.011 * T - 13.257) * 
        T2) * AE_STR;

  // Mean elongation of the Moon from the Sun.
  DD = (ae_mod_3600(1602961601.328 * T + 1072261.307) + (0.019 * T - 6.891) * 
        T2) * AE_STR;

  // Calculate sin(i * MM), etc. for needed multiple angles.
  ae_nut_sscc(0, MM, 3);
  ae_nut_sscc(1, MS, 2);
  ae_nut_sscc(2, FF, 4);
  ae_nut_sscc(3, DD, 4);
  ae_nut_sscc(4, OM, 2);

  C = 0.0;
  D = 0.0;
  p = &ae_nt_param[0]; // Point to start of table.

  for (i = 0; i < 105; i++)	{
    // Argument of sine and cosine.
	k1 = 0;
	cv = 0.0;
	sv = 0.0;
	for (m=0; m < 5; m++) {
      j = *p++;
      if (j) {
        k = j;
        if (j < 0)
          k = -k;
        su = ae_nut_ss[m][k - 1]; // sin(k*angle).
        if (j < 0)
          su = -su;
        cu = ae_nut_cc[m][k - 1];
        if( k1 == 0 ) { // Set first angle.
          sv = su;
          cv = cu;
          k1 = 1;
        }
        else { // Combine angle.
          sw = su*cv + cu*sv;
          cv = cu*cv - su*sv;
          sv = sw;
        }
      }
    }
    // Longitude coefficient.
    f = *p++;
    if ((k = *p++) != 0)
      f += T10 * k;

    // Obliquity coefficient.
    g = *p++;
    if((k = *p++) != 0)
      g += T10 * k;

    // Accumulate the terms.
    C += f * sv;
    D += g * cv;
  }

  // First terms, not in table.
  C += (-1742.0 * T10 - 171996.0) * ae_nut_ss[4][0];
  D += (89.0 * T10 + 92025.0) * ae_nut_cc[4][0];

  // Return values in seconds of arc.
  *nutl = 0.0001 * C;
  *nuto = 0.0001 * D;
}


//------------------------------------------------------------------------------
//! Correct for nutation.
//! From AA page B20.
//!
//! \param nutl The nutation in longitude in seconds of arc 
//!             (see ae_nutation_lon_ob()).
//! \param nuto The nutation in oblation in seconds of arc
//!             (see ae_nutation_lon_ob()).
//! \param eps The obliquity of the ecliptic in seconds of arc 
//!            (see ae_epsilon()).
//! \param p The equitorial rectangular position vector for mean ecliptic and
//!          equinox of date.  This is modified by the routine for nutation.
//! \param direction #AE_TO_J2000 to nutate from \p jd_tt to J2000; 
//!                  #AE_FROM_J2000 to nutate from J2000 to \p jd_tt.  (If 
//!                  neither is passed, #AE_TO_J2000 is assumed.)  To nutate 
//!                  from jd_1 to jd_2, first go from jd_1 to J2000, and then 
//!                  from J2000 to jd_2.
//------------------------------------------------------------------------------

void ae_nutate(double nutl, double nuto, double eps, double p[], 
               int direction) {
  double ce, se, cl, sl, sino, f, coseps, sineps;
  double p1[3];
  int i;

  // Initial trignometry.
  coseps = cos(AE_STR * eps);
  sineps = sin(AE_STR * eps);
  f = AE_STR * (eps + nuto);
  ce = cos(f);
  se = sin(f);
  sino = sin(AE_STR * nuto);
  cl = cos(AE_STR * nutl);
  sl = sin(AE_STR * nutl);

  // Apply adjustment to equatorial rectangular coordinates of object.  This is 
  // a composite of three rotations: rotate about x axis to ecliptic of date; 
  // rotate about new z axis by the nutation in longitude; rotate about new x 
  // axis back to equator of date plus nutation in obliquity.  For nutating back
  // to J2000, the inverse operation is performed.
  if (direction == AE_FROM_J2000) {
    p1[0] = cl * p[0] -
            sl * coseps * p[1] -
            sl * sineps * p[2];
    p1[1] = sl * ce * p[0] +
            (cl * coseps * ce + sineps * se) * p[1] -
            (sino + (1.0 - cl) * sineps * ce) * p[2];
    p1[2] = sl * se * p[0] + 
            (sino + (cl - 1.0) * se * coseps) * p[1] +
            (cl * sineps * se + coseps * ce) * p[2];
  }
  else {
    p1[0] = cl * p[0] +
            sl * ce * p[1] +
            sl * se * p[2];
    p1[1] = -sl * coseps * p[0] +
            (cl * coseps * ce + sineps * se) * p[1] +
            (coseps * cl * se - ce * sineps) * p[2];
    p1[2] = -sl * sineps * p[0] + 
            (cl * sineps * ce - se * coseps) * p[1] +
            (sineps * cl * se + ce * coseps) * p[2];
  }

  for(i = 0; i < 3; i++)
    p[i] = p1[i];
}


//------------------------------------------------------------------------------
//! Correct for nutation.
//! This is a simplified version of ae_nutate(), which automatically calculates
//! \p nutl, \p nuto and \p eps.  It might be less efficient to use this routine
//! if speed is required.
//!
//! \param jd_tt The Julian date in TT.
//! \param p The equitorial rectangular position vector for mean ecliptic and
//!          equinox of date.  This is modified by the routine for nutation.
//!
//! \param direction #AE_TO_J2000 to nutate from \p jd_tt to J2000; 
//!                  #AE_FROM_J2000 to nutate from J2000 to \p jd_tt.  (If 
//!                  neither is passed, #AE_TO_J2000 is assumed.)  To nutate 
//!                  from jd_1 to jd_2, first go from jd_1 to J2000, and then 
//!                  from J2000 to jd_2.
//------------------------------------------------------------------------------

void aes_nutate(double jd_tt, double p[], int direction) {
  double nutl, nuto, eps;

  ae_nutation_lon_ob(jd_tt, &nutl, &nuto);
  eps = ae_epsilon(jd_tt);

  ae_nutate(nutl, nuto, eps, p, direction);
}


//------------------------------------------------------------------------------
//! Do trigonometric calculations for ae_nutation_lon_ob().
//------------------------------------------------------------------------------

void ae_nut_sscc(int k, double arg, int n) {
  double cu, su, cv, sv, s;
  int i;

  su = sin(arg);
  cu = cos(arg);
  ae_nut_ss[k][0] = su;			/* sin(L) */
  ae_nut_cc[k][0] = cu;			/* cos(L) */
  sv = 2.0*su*cu;
  cv = cu*cu - su*su;
  ae_nut_ss[k][1] = sv;			/* sin(2L) */
  ae_nut_cc[k][1] = cv;

  for( i=2; i<n; i++ ) {
    s =  su*cv + cu*sv;
    cv = cu*cv - su*sv;
    sv = s;
    ae_nut_ss[k][i] = sv;		/* sin( i+1 L ) */
    ae_nut_cc[k][i] = cv;
  }
}
