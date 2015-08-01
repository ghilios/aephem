//! \file precess.c
//! Correct for precession.
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

//! Use Williams's coefficients for calculating precession.
//! Laskar's terms of order higher than \f$t^4\f$ have been retained, because 
//! Simon et al. mention that the solution is the same except for the lower 
//! order terms.
//!
//! Reference:
//! - James G. Williams, "Contributions to the Earth's obliquity rate, 
//!   precession, and nutation,"  Astron. J. 108, 711-724 (1994).
#define AE_PRECESS_USE_WILLIAMS 1

//! Use Simon's et al coefficients for calculating precession.
//! Laskar's terms of order higher than \f$t^4\f$ have been retained, because 
//! Simon et al. mention that the solution is the same except for the lower 
//! order terms.
//!
//! Reference:
//! - J. L. Simon, P. Bretagnon, J. Chapront, M. Chapront-Touze', G. Francou,
//!   and J. Laskar, "Numerical Expressions for precession formulae and mean 
//!   elements for the Moon and the planets," Astronomy and Astrophysics 282, 
//!   663-683 (1994).
#define AE_PRECESS_USE_SIMON 0

//! Use Lieske's et al IAU coefficients for calculating precession.
//! Reference:
//! - J. H. Lieske, T. Lederle, W. Fricke, and B. Morando, "Expressions for the 
//!   Precession Quantities Based upon the IAU (1976) System of Astronomical 
//!   Constants,"  Astronomy and Astrophysics 58, 1-16 (1977).
#define AE_PRECESS_USE_IAU 0

//! Use Laskar's IAU coefficients for calculating precession.
//! Newer formulas that cover a much longer time span are from:
//!   - J. Laskar, "Secular terms of classical planetary theories using the 
//!     results of general theory," Astronomy and Astrophysics 157, 59070 
//!     (1986).
//!
//! See also:
//!   - P. Bretagnon and G. Francou, "Planetary theories in rectangular and 
//!     spherical variables. VSOP87 solutions," Astronomy and Astrophysics 202,
//!     309-315 (1988).
//!
//! Laskar's expansions are said by Bretagnon and Francou to have "a precision 
//! of about 1" over 10000 years before and after J2000.0 in so far as the 
//! precession constants \f$p^0_A\f$ and \f$\epsilon^0_A\f$ are perfectly known.
//!
//! Bretagnon and Francou's expansions for the node and inclination of the 
//! ecliptic were derived from Laskar's data but were truncated after the term 
//! in \f$T^6\f$. I have recomputed these expansions from Laskar's data, 
//! retaining powers up to \f$T^10\f$ in the result.
//!
//! The following table indicates the differences between the result of the IAU 
//! formula and Laskar's formula using four different test vectors, checking at 
//! J2000 plus and minus the indicated number of years.
//!
//! <table>
//!   <tr><td>Years from J2000</td><td>Arcseconds</td></tr>
//!   <tr><td>0     </td><td>0      </td></tr>
//!   <tr><td>100   </td><td>0.006  </td></tr>
//!   <tr><td>200   </td><td>0.006  </td></tr>
//!   <tr><td>500   </td><td>0.015  </td></tr>
//!   <tr><td>1000  </td><td>0.28   </td></tr>
//!   <tr><td>2000  </td><td>6.4    </td></tr>
//!   <tr><td>3000  </td><td>38.0   </td></tr>
//!   <tr><td>10000 </td><td>9400.  </td></tr>
#define AE_PRECESS_USE_LASKAR 0

//! \var double ae_precess_coef[]
//! Precession coefficients.
//! The exact definition will depend on which of #AE_PRECESS_USE_WILLIAMS, 
//! #AE_PRECESS_USE_SIMON, #AE_PRECESS_USE_IAU or #AE_PRECESS_USE_LASKAR is set.

//! \var double ae_precess_node_coef[]
//! Precession node coefficients.
//! The exact definition will depend on which of #AE_PRECESS_USE_WILLIAMS, 
//! #AE_PRECESS_USE_SIMON, #AE_PRECESS_USE_IAU or #AE_PRECESS_USE_LASKAR is set.

//! \var double ae_precess_incl_coef[]
//! Precession inclination coefficients.
//! The exact definition will depend on which of #AE_PRECESS_USE_WILLIAMS, 
//! #AE_PRECESS_USE_SIMON, #AE_PRECESS_USE_IAU or #AE_PRECESS_USE_LASKAR is set.

#if AE_PRECESS_USE_WILLIAMS
  double ae_precess_coef[] = {
    #if 1
      // Corrections to Williams (1994) introduced in DE403.
      -8.66e-10, -4.759e-8, 2.424e-7, 1.3095e-5, 1.7451e-4, -1.8055e-3,
      -0.235316, 0.076, 110.5414, 50287.91959
    #else
      -8.66e-10, -4.759e-8, 2.424e-7, 1.3095e-5, 1.7451e-4, -1.8055e-3,
      -0.235316, 0.076, 110.5407, 50287.70000
    #endif
  };
  
  /* Pi from Williams' 1994 paper, in radians.  No change in DE403.  */
  double ae_precess_node_coef[] = {
    6.6402e-16, -2.69151e-15, -1.547021e-12, 7.521313e-12, 1.9e-10, -3.54e-9, 
    -1.8103e-7,  1.26e-7,  7.436169e-5, -0.04207794833,  3.052115282424};
  
  /* pi from Williams' 1994 paper, in radians.  No change in DE403.  */
  double ae_precess_incl_coef[] = {
    1.2147e-16, 7.3759e-17, -8.26287e-14, 2.503410e-13, 2.4650839e-11, 
    -5.4000441e-11, 1.32115526e-9, -6.012e-7, -1.62442e-5, 0.00227850649, 0.0};
#endif

#if AE_PRECESS_USE_SIMON
  double ae_precess_coef[] = {
    -8.66e-10, -4.759e-8, 2.424e-7, 1.3095e-5, 1.7451e-4, -1.8055e-3, -0.235316,
    0.07732, 111.2022, 50288.200};

  double ae_precess_node_coef[] = {
    6.6402e-16, -2.69151e-15, -1.547021e-12, 7.521313e-12, 1.9e-10, -3.54e-9, 
    -1.8103e-7, 2.579e-8, 7.4379679e-5, -0.0420782900, 3.0521126906};

  double ae_precess_incl_coef[] = {
    1.2147e-16, 7.3759e-17, -8.26287e-14, 2.503410e-13, 2.4650839e-11, 
    -5.4000441e-11, 1.32115526e-9, -5.99908e-7, -1.624383e-5, 0.002278492868, 
    0.0};
#endif

#if AE_PRECESS_USE_LASKAR
  double ae_precess_coef[] = {
    -8.66e-10, -4.759e-8, 2.424e-7, 1.3095e-5, 1.7451e-4, -1.8055e-3, -0.235316,
    0.07732, 111.1971, 50290.966 };

  double ae_precess_node_coef[] = {
    6.6402e-16, -2.69151e-15, -1.547021e-12, 7.521313e-12, 6.3190131e-10, 
    -3.48388152e-9, -1.813065896e-7, 2.75036225e-8, 7.4394531426e-5,
    -0.042078604317, 3.052112654975 };

  double ae_precess_incl_coef[] = {
    1.2147e-16, 7.3759e-17, -8.26287e-14, 2.503410e-13, 2.4650839e-11, 
    -5.4000441e-11, 1.32115526e-9, -5.998737027e-7, -1.6242797091e-5,
    0.002278495537, 0.0 };
#endif


//------------------------------------------------------------------------------
//! Precess a coordinate.
//! The precession is calcuated using #ae_precess_coef, #ae_precess_node_coef
//! and #ae_precess_incl_coef.
//!
//! \param jd_tt The Julian date in TT.
//! \param r A rectangular coordinate vector to be precessed.  This routine
//!          returns the precessed coordinate in this parameter.
//! \param direction #AE_TO_J2000 to precess from \p jd_tt to J2000; 
//!                  #AE_FROM_J2000 to precess from J2000 to \p jd_tt.  (If 
//!                  neither is passed, #AE_TO_J2000 is assumed.)  To precess 
//!                  from jd_1 to jd_2, first go from jd_1 to J2000, and then 
//!                  from J2000 to jd_2.
//------------------------------------------------------------------------------

void ae_precess(double jd_tt, double r[], int direction) {
  int i;
  double A, B, T, pA, W, z, x[3], *p, eps, coseps, sineps;
  #if AE_PRECESS_USE_IAU
    double sinth, costh, sinZ, cosZ, sinz, cosz, Z, TH;
  #endif

  if (direction != AE_TO_J2000 && direction != AE_FROM_J2000)
    direction = AE_TO_J2000;

  // Dumb case.
  if (jd_tt == AE_J2000)
    return;

  // Each precession angle is specified by a polynomial in T = Julian centuries 
  // from J2000.0.  See AA page B18.
  T = (jd_tt - AE_J2000) / 36525.0;

  #if AE_PRECESS_USE_IAU
    // Use IAU formula only for a few centuries, if at all.
    if (fabs(T) > Two)
      goto laskar;

    Z =  (( 0.017998 * T + 0.30188) * T + 2306.2181) * T * AE_STR;
    z =  (( 0.018203 * T + 1.09468) * T + 2306.2181) * T * AE_STR;
    TH = ((-0.041833 * T - 0.42665) * T + 2004.3109) * T * AE_STR;

    sinth = sin(TH);
    costh = cos(TH);
    sinZ = sin(Z);
    cosZ = cos(Z);
    sinz = sin(z);
    cosz = cos(z);
    A = cosZ * costh;
    B = sinZ * costh;

    if (direction == AE_FROM_J2000) { 
      // From J2000.0 to jd.
      x[0] = (A * cosz - sinZ * sinz) * r[0] -
	         (B * cosz + cosZ * sinz) * r[1] -
	                     sinth * cosz * r[2];

      x[1] = (A * sinz + sinZ * cosz) * r[0] -
             (B * sinz - cosZ * cosz) * r[1] -
                         sinth * sinz * r[2];

      x[2] = cosZ * sinth * r[0] -
	         sinZ * sinth * r[1] +
	                costh * r[2];
    }
    else { 
      // From jd to J2000.0.
      x[0] = (A * cosz - sinZ * sinz) * r[0] +
             (A * sinz + sinZ * cosz) * r[1] +
                         cosZ * sinth * r[2];

      x[1] = -(B * cosz + cosZ * sinz) * r[0] -
              (B * sinz - cosZ * cosz) * r[1] -
	                      sinZ * sinth * r[2];

      x[2] = -sinth * cosz * r[0] -
             sinth * sinz * r[1] +
	         costh * r[2];
    }	
    goto done;

    laskar:
  #endif

  // Implementation by elementary rotations using Laskar's expansions.  First 
  // rotate about the x axis from the initial equator to the ecliptic. 
  // (The input is equatorial.)
  if (direction == AE_TO_J2000)
    eps = ae_epsilon(jd_tt);  // To jd.
  else
	eps = ae_epsilon(AE_J2000); // From J2000.0.
  coseps = cos(AE_STR * eps);
  sineps = sin(AE_STR * eps);

  x[0] = r[0];
  z = coseps * r[1] + sineps * r[2];
  x[2] = -sineps * r[1] + coseps * r[2];
  x[1] = z;

  // Precession in longitude.
  T /= 10.0; // Thousands of years.
  p = ae_precess_coef;
  pA = *p++;
  for (i = 0; i < 9; i++)
    pA = pA * T + *p++;
  pA *= AE_STR * T;

  // Node of the moving ecliptic on the J2000 ecliptic.
  p = ae_precess_node_coef;
  W = *p++;
  for (i = 0; i < 10; i++)
    W = W * T + *p++;

  // Rotate about z axis to the node.
  if (direction == AE_TO_J2000)
    z = W + pA;
  else
    z = W;
  B = cos(z);
  A = sin(z);
  z = B * x[0] + A * x[1];
  x[1] = -A * x[0] + B * x[1];
  x[0] = z;

  // Rotate about new x axis by the inclination of the moving ecliptic on the 
  // J2000 ecliptic.
  p = ae_precess_incl_coef;
  z = *p++;
  for (i = 0; i < 10; i++)
    z = z * T + *p++;
  if (direction == AE_TO_J2000)
    z = -z;
  B = cos(z);
  A = sin(z);
  z = B * x[1] + A * x[2];
  x[2] = -A * x[1] + B * x[2];
  x[1] = z;

  // Rotate about new z axis back from the node.
  if (direction == AE_TO_J2000)
    z = -W;
  else
    z = -W - pA;
  B = cos(z);
  A = sin(z);
  z = B * x[0] + A * x[1];
  x[1] = -A * x[0] + B * x[1];
  x[0] = z;

  // Rotate about x axis to final equator.
  if (direction == AE_TO_J2000)
    eps = ae_epsilon(AE_J2000);
  else
    eps = ae_epsilon(jd_tt);
  coseps = cos(AE_STR * eps);
  sineps = sin(AE_STR * eps);
  z = coseps * x[1] - sineps * x[2];
  x[2] = sineps * x[1] + coseps * x[2];
  x[1] = z;

  #if AE_PRECESS_USE_IAU
    done:
  #endif

  for (i = 0; i < 3; i++)
    r[i] = x[i];

  return;
}
