//! \file calendar.c
//! This file contains functions for calendar and time conversions.
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
#include <time.h>

#include "aephem.h"

//------------------------------------------------------------------------------
//! Get the Greenwich mean sidereal time.
//! Get the mean sidereal time at Greenwich (i.e., longitude 0).  The mean
//! sidereal time does not include nutation corrections.  Coefficients are from 
//! the IAU and can be found in:
//! - George H. Kaplan, "The IAU Resolutions on Astronomical Reference Systems,
//!   Time Scales, and Earth Rotation Models," United States Naval Observatory
//!   Circular No. 179, 2005.
//!
//! \param jd_ut1 The Julian date, in UT1.
//! \param jd_tt The Julian date, in TT.  (See ae_delta_t().)
//!
//! \return The Greenwich mean sidereal time, in degrees.
//------------------------------------------------------------------------------
  
double ae_gmst(double jd_ut1, double jd_tt) {
  double jd0, gmst, T0, T, theta, d_u;
    
  // Julian day at given UT1.
  jd0 = floor(jd_ut1);
    
  // Julian centuries from standard epoch J2000.0.
  T0 = (jd0 - AE_J2000) / 36525.0;
  
  // Julian days since J2000.0.
  d_u = jd_ut1 - AE_J2000;
  
  // Get the earth rotation angle.
  theta = 0.7790572732640 + 0.00273781191135448 * d_u + (jd_ut1 - jd0);
    
    
  // Centuries since J2000, in TT.
  T = (jd_tt - AE_J2000) / 36525.0;

  gmst = theta * 360.0 + (0.014506 + T * (4612.156534 + T * (1.3915817 -
                 T * (4.4e-7 - T * (2.9956e-5 - T * 3.68e-8))))) * AE_STD;
    
  return ae_mod_360(gmst);
}   
    

//------------------------------------------------------------------------------
//! Get Greenwich apparent sidereal time.
//! The Greenwich apparent sidereal time is the Greenwich mean sidereal time
//! corrected for nutation.  See ae_gmst() for reference.
//!
//! \param jd_ut1 The Julian date, in UT1.
//! \param jd_tt The Julian date, in TT.  (See ae_delta_t().)
//! \param nutl The nutation in longitude in seconds of arc (see 
//!             ae_nutation_lon_ob()).
//! \param eps The obliquity of the ecliptic in seconds of arc 
//!            (see ae_epsilon()).
//!
//! \return Greenwich apparent sidereal time in degrees.
//------------------------------------------------------------------------------

double ae_gast(double jd_ut1, double jd_tt, double nutl, double eps) {
  double gast;

  gast = ae_gmst(jd_ut1, jd_tt) + nutl * AE_STD * cos(AE_STR * eps);

  return ae_mod_360(gast);
}


//------------------------------------------------------------------------------
//! Get the local mean sidereal time.
//! The local mean sidereal time is simply the Greenwich mean sidereal time
//! plus the local longitude.  See ae_gmst().
//!
//! \param jd_ut1 The Julian date in UT1.
//! \param jd_tt The Julian date in TT.  (See ae_delta_t().)
//! \param tlong The longitude of the observer in degrees.
//!
//! \return Local sidereal time in degrees.
//------------------------------------------------------------------------------

double ae_lmst(double jd_ut1, double jd_tt, double tlong) {
  return ae_mod_360(ae_gmst(jd_ut1, jd_tt) + tlong);
}


//------------------------------------------------------------------------------
//! Get the local apparent sidereal time.
//! The local apparent sidereal time is simply the Greenwich apparent sidereal 
//! time plus the local longitude.  See ae_gast().
//!
//! \param jd_ut1 The Julian date in UT1.
//! \param jd_tt The Julian date in TT.  (See ae_gast().)
//! \param tlong The longitude of the observer in degrees.
//! \param nutl The nutation in longitude in seconds of arc (see 
//!             ae_nutation_lon_ob()).
//! \param eps The obliquity of the ecliptic in seconds of arc 
//!            (see ae_epsilon()).
//!
//! \return Local sidereal time in degrees.
//------------------------------------------------------------------------------

double ae_last(double jd_ut1, double jd_tt, double tlong, double nutl,
               double eps) {
  return ae_mod_360(ae_gast(jd_ut1, jd_tt, nutl, eps) + tlong);
}


//------------------------------------------------------------------------------
//! Get the local apparent sidereal time.
//! This is a simplified version of ae_last(), which automatically calculates
//! \p jd_tt, \p nutl and \p eps.  It might be less efficient to use this 
//! routine if speed is required.
//!
//! \param jd_ut1 The Julian date in UT1.
//! \param tlong The longitude of the observer in degrees.
//!
//! \return Local sidereal time in degrees.
//------------------------------------------------------------------------------

double aes_last(double jd_ut1, double tlong) {
  double jd_tt, nutl, nuto, eps;
  
  jd_tt = jd_ut1 + ae_delta_t(jd_ut1) * AE_D_PER_S;
  eps = ae_epsilon(jd_tt);
  ae_nutation_lon_ob(jd_tt, &nutl, &nuto);
  
  return ae_last(jd_ut1, jd_tt, tlong, nutl, eps);
}


//------------------------------------------------------------------------------
//! Find Barycentric Dynamical Time from Terrestrial Dynamical Time.
//! See AA page B5.
//!
//! \param jd A Julian date, in TDT.
//!
//! \return The corresponding time in TDB.
//------------------------------------------------------------------------------

double ae_tdb(double jd) {
  double M, T;

  // Find time T in Julian centuries from J2000.
  T = (jd - AE_J2000) / 36525.0;

  // Mean anomaly of sun = l' (J. Laskar).
  M = 129596581.038354 * T +  1287104.76154;

  // Reduce arc seconds mod 360 degrees.
  M = M - 1296000.0 * floor(M / 1296000.0);

  M += ((((((((1.62e-20 * T - 1.0390e-17 ) * T - 3.83508e-15 ) * T +
               4.237343e-13 ) * T + 8.8555011e-11 ) * T - 4.77258489e-8 ) * T -
               1.1297037031e-5 ) * T + 1.4732069041e-4 ) * T -
               0.552891801772 ) * T * T;
  M *= AE_STR;

  // TDB - TDT, in seconds.
  T = 0.001658 * sin(M) + 0.000014 * sin(M + M);

  T = jd + T / 86400.0;

  return T;
}


//-----------------------------------------------------------------------------
//! Calculate Julian date from a Gregorian calendar date.
//! The Julian date is double precision floating point with the origin used by
//! astronomers. 
//!
//! There is no year 0.  Enter B.C. years as negative; i.e., 2 B.C. = -2.
//! 
//! The approximate range of dates handled is 4713 B.C. to 54,078 A.D.  This 
//! should be adequate for most applications.
//!
//! B.C. dates are calculated by extending the Gregorian sequence of leap years 
//! and century years into the past.  This seems the only sensible definition, 
//! but it may not be the official one.
//!
//! Note that the astronomical Julian day starts at noon on the previous 
//! calendar day.  Thus at midnight in the morning of the present calendar day 
//! the Julian date ends in 0.5; it rolls over to tomorrow at noon today.
//!
//! The month-finding algorithm is attributed to Meeus.
//!
//! \param year The Gregorian year.
//! \param month The month of the year, with January = 1, February = 2, etc.
//! \param day The fractional day of the month, starting at 1.
//!
//! \return The Julian day number.
//------------------------------------------------------------------------------

double ae_cal_to_jd(long year, int month, double day) {
  long y, a, b, c, e, m;

  // The origin should be chosen to be a century year that is also a leap year. 
  // We pick 4801 B.C.
  y = year + 4800;
  if (year < 0)
    y += 1;

  // The following magic arithmetic calculates a sequence whose successive terms
  // differ by the correct number of days per calendar month.  It starts at 
  // 122 = March; January and February come after December.
  m = month;
  if (m <= 2) {
    m += 12;
    y -= 1;
  }
  e = (306 * (m + 1)) / 10;
  
  a = y / 100; // Number of centuries.
  if (year <= 1582L) {
    if (year == 1582L) {
      if (month < 10)
        goto julius;
      if (month > 10)
        goto gregor;
      if (day >= 15)
        goto gregor;
    }

    julius:
      // printf("Aephem: ae_caltojd: Julian Calendar assumed!\n");
      b = -38;
  }
  else { 
    // -Number of century years that are not leap years.
    gregor:
      b = (a / 4) - a;
  }

  c = (36525L * y) / 100; // Julian calendar years and leap years.

  // Add up these terms, plus offset from J 0 to 1 Jan 4801 B.C.  Also fudge for
  // the 122 days from the month algorithm.
  return b + c + e + day - 32167.5;
}


//------------------------------------------------------------------------------
//! Calculate month, day, and year from a Julian date.
//!
//! \param jd A Julian day number.
//! \param year For returning the year of \p jd.
//! \param month For returning the month of \p jd, in the range [1:12].
//! \param day For returning the fractional day of \p jd, with the first day of
//!            the month being 1.
//------------------------------------------------------------------------------

void ae_jd_to_cal(double jd, int *year, int *month, double *day) {
  long a, c, d, x, y, jd_adj;
  int BC;

  if (jd < 1721425.5) // January 1.0, 1 A.D.
    BC = 1;
  else
    BC = 0;

  jd_adj = (long) (jd + 0.5); // Round Julian date up to integer.

  // Find the number of Gregorian centuries since March 1, 4801 B.C.
  a = (100 * jd_adj + 3204500L) / 3652425L;

  // Transform to Julian calendar by adding in Gregorian century years that 
  // are not leap years. Subtract 97 days to shift origin of JD to March 1.
  // Add 122 days for magic arithmetic algorithm.  Add four years to ensure the 
  // first leap year is detected.
  c = jd_adj + 1486;
  if( jd_adj >= 2299160.5 )
    c += a - a / 4;
  else
    c += 38;

  // Offset 122 days, which is where the magic arithmetic month formula sequence
  // starts (March 1 = 4 * 30.6 = 122.4).
  d = (100 * c - 12210L) / 36525L;
  // Days in that many whole Julian years.
  x = (36525L * d) / 100L;

  // Find month and day.
  y = ((c - x) * 100L) / 3061L;
  *day = (double)((int)(c - x - ((306L * y)/ 10L)));
  *month = (int) (y - 1);
  if (y > 13)
    *month -= 12;

  // Get the year right.
  *year = d - 4715;
  if (*month > 2)
    *year -= 1;

  // Fractional part of day.
  *day += jd - jd_adj + 0.5;

  if (BC)
    *year = -(*year) + 1;
}




//------------------------------------------------------------------------------
//! Convert a C time to a Julian date.
//! The C time is the time since the Epoch (00:00:00 UTC, January 1, 1970),
//! measured in seconds.
//!
//! \param t The C time at which to calculate the Julian date.
//!
//! \return The Julian date.
//------------------------------------------------------------------------------

double ae_ctime_to_jd(double t) {
  time_t t2;
  double d;
  struct tm *t3;

  t2 = (time_t)t;
  t3 = gmtime(&t2);
  d = ((double)(t3->tm_hour * 3600 + t3->tm_min * 60 + t3->tm_sec) +
       t - (double)t2) * AE_STDAY + (double)t3->tm_mday;

  return ae_cal_to_jd(t3->tm_year + 1900, t3->tm_mon + 1, d);
}


//------------------------------------------------------------------------------
//! Calculate the local apparent sidereal time from a C time.
//! The C time is the time since the Epoch (00:00:00 UTC, January 1, 1970),
//! measured in seconds.
//!
//! \param t The UNIX time in UT1 measure at which to calculate the LST.
//! \param delta_t TT - UT1 (see ae_delta_t()), in seconds.
//! \param tlong The longitude of the observer, in degrees.
//! \param nutl The nutation in longitude in seconds of arc 
//!             (see ae_nutation_lon_ob()).
//! \param eps The obliquity of the ecliptic in seconds of arc 
//!            (see ae_epsilon()).
//!
//! \return The LAST in degrees.
//------------------------------------------------------------------------------
    
double ae_ctime_to_last(double t, double delta_t, double tlong, double nutl,
                        double eps) {
  double jd_ut1, jd_tt;
  
  jd_ut1 = ae_ctime_to_jd(t);
  jd_ut1 += ae_dut1(jd_ut1) * AE_D_PER_S;
  jd_tt = jd_ut1 + delta_t * AE_D_PER_S;

  return ae_last(jd_ut1, jd_tt, tlong, nutl, eps);
} 
  
  
//------------------------------------------------------------------------------
//! Calculate the local sidereal time from a C time.
//! This is a simplification of ae_ctime_to_last().  It automatically does
//! the calculations of \p delta_t, \p nutl and \p eps.  Use ae_ctime_to_last()
//! for more efficiency, especially with repeated calls.
//!
//! \param t The UNIX time at which to calculate the LST.
//! \param tlong The longitude of the observer, in degrees.
//!
//! \return The LAST in degrees.
//------------------------------------------------------------------------------
  
double aes_ctime_to_last(double t, double tlong) {
  double jd_ut1;

  jd_ut1 = ae_ctime_to_jd(t);
  jd_ut1 += ae_dut1(jd_ut1) * AE_D_PER_S;

  return aes_last(jd_ut1, tlong);
}


//------------------------------------------------------------------------------
//! Given a Julian Date, get the Modified Julian Date.
//! The difference between the two is simply a constant, #AE_MJD_START.
//!
//! \return The Modified Julian Date.
//------------------------------------------------------------------------------

double ae_mjd(double jd) {
  return jd - AE_MJD_START;
}


extern short ae_delta_t_tab[];
extern const int ae_delta_t_start_year;
extern const int ae_delta_t_end_year;

//------------------------------------------------------------------------------
//! Get the value of delta T == Ephemeris Time - Universal Time.
//! This routine uses the table #ae_delta_t_tab for historic values and
//! near-future predictions.
//!
//! The program adjusts for a value of secular tidal acceleration ndot. It is 
//! -25.8 arcsec per century squared for JPL's DE403 ephemeris.  ELP2000 and 
//! DE200 use the value -23.8946.
//!
//! For dates earlier than the tabulated range, the program calculates
//! approximate formulae of Stephenson and Morrison or K. M. Borkowski.  These 
//! approximations have an estimated error of 15 minutes at 1500 B.C.  They are
//! not adjusted for small improvements in the current estimate of ndot because 
//! the formulas were derived from studies of ancient eclipses and other 
//! historical information, whose interpretation depends only partly on ndot.
//! 
//! A quadratic extrapolation formula, that agrees in value and slope with
//! current data, predicts future values of deltaT.
//!
//! References:
//!
//!   - Stephenson, F. R., and L. V. Morrison, "Long-term changes in the
//!     rotation of the Earth: 700 B.C. to A.D. 1980," Philosophical 
//!     Transactions of the Royal Society of London Series A 313, 47-70 (1984)
//!   - Borkowski, K. M., "ELP2000-85 and the Dynamical Time - Universal Time 
//!     relation," Astronomy and Astrophysics 205, L8-L10 (1988)
//!     - Borkowski's formula is derived from eclipses going back to 2137 BC
//!       and uses lunar position based on tidal coefficient of -23.9 
//!       arcsec/cy^2.
//!   - Chapront-Touze, Michelle, and Jean Chapront, _Lunar Tables and Programs 
//!     from 4000 B.C. to A.D. 8000_, Willmann-Bell 1991
//!     - Their table agrees with the one here, but the entries are rounded to 
//!       the nearest whole second.
//!   - Stephenson, F. R., and M. A. Houlden, _Atlas of Historical Eclipse
//!     Maps_, Cambridge U. Press (1986)
//!
//! \param jd_ut1 The Julian date.
//!
//! \return The ET - UT (delta T) in seconds.
//------------------------------------------------------------------------------

double ae_delta_t(double jd_ut1) {
  double ans, p, B, y;
  int d[6], delta_t_range, i, iy, k;

  delta_t_range = (ae_delta_t_end_year - ae_delta_t_start_year + 1);

  // Convert to centuries.
  y = 2000.0 + (jd_ut1 - AE_J2000) / 365.25;

  if (y > ae_delta_t_end_year) {
    #if 0
      // Morrison, L. V. and F. R. Stephenson, "Sun and Planetary System"
      // vol 96,73 eds. W. Fricke, G. Teleki, Reidel, Dordrecht (1982)
      B = 0.01*(y-1800.0) - 0.1;
      ans = -15.0 + 32.5*B*B;
      return(ans);
    #else
      // Extrapolate forward by a second-degree curve that agrees with the most
      // recent data in value and slope, and vaguely fits over the past century.
      // This idea communicated by Paul Muller, who says NASA used to do 
      // something like it.
      B = y - ae_delta_t_end_year;
      // Slope.
      p = ae_delta_t_tab[delta_t_range - 1] -
          ae_delta_t_tab[delta_t_range - 2];
      // Square term.
      ans = (ae_delta_t_tab[delta_t_range - 101] -
            (ae_delta_t_tab[delta_t_range - 1] - 100.0 * p)) * 1e-4;
      ans = 0.01 * (ae_delta_t_tab[delta_t_range-1] + p * B + ans * B * B);

      return ans;
    #endif
  }

  if (y < ae_delta_t_start_year) {
    if (y >= 948.0) {
      // Stephenson and Morrison, stated domain is 948 to 1600: 
      //   25.5 * (centuries from 1800)^2 - 1.9159 * (centuries from 1955)^2
      B = 0.01 * (y - 2000.0);
      ans = (23.58 * B + 100.3) * B + 101.6;
    }
    else {
      // Borkowski.
      B = 0.01 * (y - 2000.0) + 3.75;
      ans = 35.0 * B * B + 40.0;
    }

    return ans;
  }

  // Besselian interpolation from tabulated values.  See AA page K11.

  // Index into the table.
  p = floor(y);
  iy = (int)(p - ae_delta_t_start_year);

  // Zeroth order estimate is value at start of year.
  ans = ae_delta_t_tab[iy];
  k = iy + 1;
  if (k >= delta_t_range)
    goto done;

  // The fraction of tabulation interval.
  p = y - p;

  // First order interpolated value.
  ans += p * (ae_delta_t_tab[k] - ae_delta_t_tab[iy]);
  if ((iy - 1 < 0) || (iy + 2 >= delta_t_range))
    goto done;

  // Make table of first differences.
  k = iy - 2;
  for (i = 0; i < 5; i++) {
    if ((k < 0) || (k + 1 >= delta_t_range))
      d[i] = 0;
    else
      d[i] = ae_delta_t_tab[k + 1] - ae_delta_t_tab[k];
    k += 1;
  }

  // Compute second differences.
  for (i = 0; i < 4; i++)
    d[i] = d[i + 1] - d[i];
  B = 0.25 * p * (p - 1.0);
  ans += B * (d[1] + d[2]);
  if (iy + 2 >= delta_t_range)
    goto done;

  // Compute third differences.
  for (i = 0; i < 3; i++)
    d[i] = d[i + 1] - d[i];
  B = 2.0 * B / 3.0;
  ans += (p - 0.5) * B * d[1];
  if ((iy - 2 < 0) || (iy + 3 > delta_t_range))
    goto done;

  // Compute fourth differences.
  for (i = 0; i < 2; i++)
    d[i] = d[i + 1] - d[i];
  B = 0.125 * B * (p + 1.0) * (p - 2.0);
  ans += B * (d[0] + d[1]);
  
  done:
  
  // Astronomical Almanac table is corrected by adding the expression
  //   -0.000091 (ndot + 26)(year-1955)^2  seconds
  // to entries prior to 1955 (AA page K8), where ndot is the secular tidal term
  // in the mean motion of the Moon.  Entries after 1955 are referred to atomic 
  // time standards and are not affected by errors in Lunar or planetary theory.
  ans *= 0.01;
  if (y < 1955.0) {
    B = (y - 1955.0);
    #if 1 
      ans += -0.000091 * (-25.8 + 26.0) * B * B;
    #else
      ans += -0.000091 * (-23.8946 + 26.0) * B * B;
    #endif 
  }

  return ans;
}
  
extern const double ae_dut1_first_mjd;
extern const double ae_dut1_last_mjd;
extern short ae_dut1_tab[];
extern const double ae_dut1_predict_offset;
extern const double ae_dut1_predict_scale;
extern const double ae_dut1_predict_mjd;

//------------------------------------------------------------------------------
//! Given a UTC, return DUT1.
//! DUT1 is the offset between UTC and UT1, i.e., \f$ DUT1 \equiv UT1 - UTC \f$.
//! Functions that find sidereal time expect UT1, so this function will be
//! required to convert UTC if reasonable precision (better than roughly 10
//! seconds of arc).
//!
//! Values are taken from a look-up table.  Past values start on 19 May 1976
//! are are tabulated daily until the release date of this version of aephem.
//! After this, predicted values are used through about one year.
//!
//! If DUT1 is requested for a date before the tabulation starts, 0 will be 
//! returned.  After the last predicted value, a value will be estimated using 
//! an interpolation formula.
//!
//! Past values, predicted values and the interpolation formula were taken from 
//! the U.S. Naval Observatory IERS Bulletin A.  This is available online at
//! http://maia.usno.navy.mil/.
//!
//! \param jd_utc The Julian date in UTC.
//!
//! \return DUT1, in seconds. If jd_utc is before the first tabulated value,
//!         0 will be returned.
//------------------------------------------------------------------------------

double ae_dut1(double jd_utc) {
  int i;
  double mjd, dt, t, dut2, B;

  mjd = ae_mjd(jd_utc);

  // Case 1:  before first tabulated values.
  if (mjd < ae_dut1_first_mjd)
    return 0;

  // Case 2:  a tabulated value.
  if (mjd < ae_dut1_last_mjd) {
    i = (int)(mjd - ae_dut1_first_mjd);
    dt = mjd - floor(mjd);
    t = (double)ae_dut1_tab[i] * (1.0 - dt) + 
        (double)ae_dut1_tab[i + 1] * dt;

    return t / 1e4;
  }

  // Case 3:  beyond the tabulated predicted values.
  // First get UT2 - UT1.
  B = 1900.0 + (jd_utc - 2415020.31352) / 365.242198781; // Besselian year.

  dut2 = 0.022 * sin(2 * M_PI * B) - 0.012 * cos(2 * M_PI * B) -
         0.006 * sin(4 * M_PI * B) + 0.007 * cos(4 * M_PI * B);

  t = ae_dut1_predict_offset - 
      ae_dut1_predict_scale * (mjd - ae_dut1_predict_mjd) - dut2;
  
  return t;
}
