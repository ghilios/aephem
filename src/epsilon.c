//! \file epsilon.c
//! Calculate the obliquity of the ecliptic.
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

#include <stdio.h>

#include "aephem.h"

//! Use William's obliquity calculations (c.f. #AE_PRECESS_USE_WILLIAMS).
#define AE_EPSILON_USE_WILLIAMS 1
//! Use Simon's et al obliquity calculations (c.f. #AE_PRECESS_USE_SIMON).
#define AE_EPSILON_USE_SIMON 0

//------------------------------------------------------------------------------
//! Calculate the obliquity of the ecliptic.
//! J. L. Simon, P. Bretagnon, J. Chapront, M. Chapront-Touze', G. Francou,
//! and J. Laskar, "Numerical Expressions for precession formulae and
//! mean elements for the Moon and the planets," Astronomy and Astrophysics
//! 282, 663-683 (1994).
//!
//! IAU Coefficients are from:
//!   J. H. Lieske, T. Lederle, W. Fricke, and B. Morando,  "Expressions for 
//!   the Precession Quantities Based upon the IAU (1976) System of Astronomical
//!   Constants,"  Astronomy and Astrophysics 58, 1-16 (1977).
//!
//! Before or after 200 years from J2000, the formula used is from:
//!   J. Laskar, "Secular terms of classical planetary theories using the 
//!   results of general theory," Astronomy and Astrophysics 157, 59070 (1986).
//!
//!  See precess.c and page B18 of the Astronomical Almanac.
//!
//! \param jd_tt The Julian date in TT.
//!
//! \return The obliquity, in seconds of arc.
//------------------------------------------------------------------------------

double ae_epsilon(double jd_tt) {
  double eps;
  double T;

  T = (jd_tt- 2451545.0) / 36525.0;

  #if AE_EPSILON_USE_WILLIAMS
    // DE403 values.
    T /= 10.0;
    eps = (((((((((2.45e-10 * T + 5.79e-9) * T + 2.787e-7) * T + 7.12e-7) * T - 
                   3.905e-5) * T - 2.4967e-3) * T - 5.138e-3) * T + 1.9989) *
                   T - 0.0175) * T - 468.33960) * T + 84381.406173;
  #else
    // This expansion is from the AA.
    // Note the official 1976 IAU number is 23d 26' 21.448", but the JPL 
    // numerical integration found 21.4119".
    #if AE_EPSILON_USE_SIMON
	  T /= 10.0;
      eps = (((((((((2.45e-10 * T + 5.79e-9) * T + 2.787e-7) * T + 7.12e-7) * 
                     T - 3.905e-5) * T - 2.4967e-3) * T	- 5.138e-3) * T + 
                     1.9989) * T - 0.0152) * T - 468.0927) * T + 84381.412;
    #else
      if (fabs(T) < 2.0)
        eps = ((1.813e-3 * T - 5.9e-4) * T - 46.8150) * T + 84381.448;
      else {
        // This expansion is from Laskar, cited above.  Bretagnon and Simon say,
        // in Planetary Programs and Tables, that it is accurate to 0.1" over a 
        // span of 6000 years. Laskar estimates the precision to be 0.01" after 
        // 1000 years and a few seconds of arc after 10000 years.
        eps = (((((((((2.45e-10 * T + 5.79e-9) * T + 2.787e-7) * T + 7.12e-7) * 
                       T - 3.905e-5)*T - 2.4967e-3) * T - 5.138e-3) * T + 
                       1.99925) * T - 0.0155) * T - 468.093) * T + 84381.448;
      }
    #endif
  #endif

  return eps;
}
