//! \file refrac.c
//! Correction for atmospheric refraction in visible bands.
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
//! Atmospheric refraction correction (in visible bands).
//! This function is suitable for passing to ae_topocentric().
//!
//! For high altitude angle, AA page B61 is used.  The accuracy is `usually 
//! about 0.1 arcsecond'.
//!
//! The formula for low altitude is from the Almanac for Computers. It gives 
//! the correction for observed altitude, so has to be inverted numerically to 
//! get the observed from the true.  The accuracy is about 0.2' for 
//! -20C < T < +40C and 970mb < P < 1050mb.
//!
//! \param alt The altitude of the observation, in degrees.
//! \param param For passing the atmospheric pressure, in millibar, and the
//!              temperature, in degrees centigrade.  These should be stored as
//!              doubles, such that the pressure is in ((double *) \p param)[0]
//!              and the temperature in ((double *) \p param)[1].
//!
//! \return The correction in degrees to be added to true altitude to obtain 
//!         apparent altitude.
//------------------------------------------------------------------------------

double ae_refrac_visible(double alt, void *param) {
  int i;
  double y, y0, D0, N, D, P, Q, atpress, attemp;

  atpress = ((double *)param)[0];
  attemp = ((double *)param)[1];

  if ((alt < -2.0) || (alt >= 90.0))
    return(0.0);

  // Low angle correction.
  if(alt > 15.0) {
    D = 0.00452 * atpress / ((273.0 + attemp) * tan(AE_DTR * alt));
    return(D);
  }
  // Start iteration assuming correction = 0.
  y = alt;
  D = 0.0;
  // Invert Almanac for Computers formula numerically.
  P = (atpress - 80.0) / 930.0;
  Q = 4.8e-3 * (attemp - 10.0);
  y0 = y;
  D0 = D;

  for (i=0; i < 4; i++) {
    N = y + (7.31 / (y + 4.4));
    N = 1.0 / tan(AE_DTR * N);
    D = N * P / (60.0 + Q * (N + 39.0));
    N = y - y0;
    y0 = D - D0 - N; // Denominator of derivative.

    // Newton iteration with numerically estimated derivative.
    if ((N != 0.0) && (y0 != 0.0))
		N = y - N*(alt + D - y)/y0;
    else // Can't do it on first pass.
      N = alt + D;

    y0 = y;
    D0 = D;
    y = N;
  }

  return D;
}


//------------------------------------------------------------------------------
//! Atmospheric refraction correction (in millimetre bands).
//! This function is suitable for passing to ae_topocentric().
//!
//! This function uses the equations from Ulich 1981.  It does not cite
//! which millimetre wavelength ranges it is expected to be good for.  The ALMA
//! memo 366 does analysis on another model and finds it good to 2% up to
//! frequencies of 1000 GHz.  (The worst is at 500 GHz.)
//!
//! Ulich claims accuracy to 2" above altitudes of 3 degrees.
//!
//! References:
//! -B. L. Ulich.  "Millimeter Wave Radio Telescopes: Gain and Pointing 
//! Charactersitics".  <em>International Journal of Infrared and Millimeter 
//! Waves</em>, <strong>2</strong>, 2 (1981).
//! - J. G. Mangum. "ALMA Memo 366:  A Telescope Pointing Algorithm for ALMA".
//!   <em>Available at 
//!   http://www.alma.nrao.edu/memos/html-memos/abstracts/abs366.html.</em>
//!
//! \param alt The altitude of the observation, in degrees.
//! \param param For passing the atmospheric pressure, in hPa, the
//!              temperature, in degrees centigrade and the surface humidity, as
//!              a percentage.  These should be stored as
//!              doubles, such that:
//!              -pressure is in ((double *) \p param)[0],
//!              -temperature is in ((double *) \p param)[1].
//!              -and humidity is in ((double *) \p param)[2].
//!
//! \return The correction in degrees to be added to true altitude to obtain 
//!         apparent altitude.
//------------------------------------------------------------------------------

double ae_refrac_ulich(double alt, void *param) {
  double p_t, t_k, p_w, r_h, r, f, t, p, h;

  p = ((double *)param)[0];
  t = ((double *)param)[1];
  h = ((double *)param)[2];
  
  t_k = t + 273.15;
  p_t = p * 100.0 / 101325.0 * 760.0;
  r_h = h;
  
  p_w = r_h * pow(10.0, 23.636 - 2948.0 / t_k - 5.0 * log10(t_k)) / 100;
  
  r = 21.36 * p_t / t_k - 1.66 * p_w / t_k + 103030 * p_w / t_k / t_k;

  f = cos(alt * AE_DTR) / (sin(alt * AE_DTR) + 0.00175 * 
                           tan((87.5 - alt) * AE_DTR));
  
  return r * f * AE_STD;
}
