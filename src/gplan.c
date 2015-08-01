//! \file gplan.c
//! Routines to chew through tables of perturbations.
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

//! A modulo macro.
#define ae_gplan_mods3600(x) ((x) - 1.296e6 * floor ((x)/1.296e6))

//! The length of #ae_gplan_freq and #ae_gplan_phase.
#define AE_GPLAN_N_HARMONIC 9

//! Orbital frequency harmonics.
//! Units are seconds of arc per 10000 Julian years.  From Simon et al. (1994).
double ae_gplan_freq[AE_GPLAN_N_HARMONIC] = {
  53810162868.8982,
  21066413643.3548,
  12959774228.3429,
  6890507749.3988,
  1092566037.7991,
  439960985.5372,
  154248119.3933,
  78655032.0744,
  52272245.1795
};

//! Orbital phase harmonics.
//! Units are seconds of arc.  From Simon et al. (1994).
double ae_gplan_phase[AE_GPLAN_N_HARMONIC] = {
  252.25090552 * 3600.,
  181.97980085 * 3600.,
  100.46645683 * 3600.,
  355.43299958 * 3600.,
  34.35151874 * 3600.,
  50.07744430 * 3600.,
  314.05500511 * 3600.,
  304.34866548 * 3600.,
  860492.1546
};

//! Sine look-up table for the gplan() routines.
double ae_gplan_ss[AE_PLANTBL_N_HARMONIC][31];
//! Cosine look-up table for the gplan() routines.
double ae_gplan_cc[AE_PLANTBL_N_HARMONIC][31];


//------------------------------------------------------------------------------
//! Prepare lookup table of sin and cos (i * Lj) for required multiple angles.
//!
//! \param k The first index.
//! \param arg The sin/cos argument.
//! \param n The second index.
//------------------------------------------------------------------------------

void ae_gplan_sscc(int k, double arg, int n) {
  double cu, su, cv, sv, s;
  int i;

  su = sin(arg);
  cu = cos(arg);
  ae_gplan_ss[k][0] = su;  // Sin(L).
  ae_gplan_cc[k][0] = cu;  // Cos(L).
  sv = 2.0 * su * cu;
  cv = cu * cu - su * su;
  ae_gplan_ss[k][1] = sv;  // Sin(2L).
  ae_gplan_cc[k][1] = cv;

  for (i = 2; i < n; i++) {
    s = su * cv + cu * sv;
    cv = cu * cv - su * sv;
    sv = s;
    ae_gplan_ss[k][i] = sv;  // Sin((i + 1)L).
    ae_gplan_cc[k][i] = cv;
  }
}

//------------------------------------------------------------------------------
//! Use a table to find a planet's polar, heliocentric position.
//!
//! \param jd_tt The Julian date in TT.
//! \param plan The planet's orbital table.
//! \param pobj For returning the polar position of the planet.
//------------------------------------------------------------------------------

void ae_gplan(double jd_tt, struct ae_plantbl_t *plan, double pobj[]) {
  register double su, cu, sv, cv, T;
  double t, sl, sb, sr;
  int i, j, k, m, n, k1, ip, np, nt;
  char *p;
  double *pl, *pb, *pr;

  T = (jd_tt - AE_J2000) / plan->timescale;
  n = plan->maxargs;

  // Calculate sin(i * MM), etc. for needed multiple angles.
  for (i = 0; i < n; i++) {
    if ((j = plan->max_harmonic[i]) > 0) {
      if (i < AE_GPLAN_N_HARMONIC)
        sr = (ae_gplan_mods3600(ae_gplan_freq[i] * T) + 
                                ae_gplan_phase[i]) * AE_STR;
      else
        sr = 0;
      ae_gplan_sscc(i, sr, j);
    }
  }

  // Point to start of table of arguments.
  p = plan->arg_tbl;

  // Point to tabulated cosine and sine amplitudes.
  pl = (double *)plan->lon_tbl;
  pb = (double *)plan->lat_tbl;
  pr = (double *)plan->rad_tbl;
  sl = 0.0;
  sb = 0.0;
  sr = 0.0;
  for (;;) {
    // Argument of sine and cosine.
    // Number of periodic arguments.
    np = *p++;
    if (np < 0)
      break;
    if (np == 0) {
      nt = *p++;

      // Longitude polynomial.
      cu = *pl++;

      for (ip = 0; ip < nt; ip++)
        cu = cu * T + *pl++;
      sl +=  ae_gplan_mods3600(cu);

      // Latitude polynomial.
      cu = *pb++;
      for (ip = 0; ip < nt; ip++)
        cu = cu * T + *pb++;
      sb += cu;

	  // Radius polynomial.
	  cu = *pr++;
	  for (ip = 0; ip < nt; ip++)
        cu = cu * T + *pr++;
      sr += cu;

      continue;
    }

    k1 = 0;
    cv = 0.0;
    sv = 0.0;

    for (ip = 0; ip < np; ip++)	{
      // Which harmonic.
	  j = *p++;

	  // Which planet.
      m = *p++ - 1;

      if (j) {
        k = j;
        if (j < 0)
          k = -k;
        k -= 1;
        su = ae_gplan_ss[m][k];	  // Sin(k*angle).
        if (j < 0)
          su = -su;
        cu = ae_gplan_cc[m][k];
        if (k1 == 0) {
          // Set first angle.
          sv = su;
          cv = cu;
          k1 = 1;
        }
        else {
          // Combine angles.
          t = su * cv + cu * sv;
          cv = cu * cv - su * sv;
          sv = t;
        }
      }
	}

    // Highest power of T.
    nt = *p++;

    // Longitude.
    cu = *pl++;
    su = *pl++;
    for (ip = 0; ip < nt; ip++) {
      cu = cu * T + *pl++;
      su = su * T + *pl++;
    } 
    sl += cu * cv + su * sv;
     
    // Latitiude.
    cu = *pb++;
    su = *pb++;
    for (ip = 0; ip < nt; ip++)	{
      cu = cu * T + *pb++;
      su = su * T + *pb++;
    }
    sb += cu * cv + su * sv;

    // Radius.
    cu = *pr++;
    su = *pr++;
    for (ip = 0; ip < nt; ip++)	{
      cu = cu * T + *pr++;
      su = su * T + *pr++;
    }
    sr += cu * cv + su * sv;
  }

  pobj[0] = AE_STR * sl;
  pobj[1] = AE_STR * sb;
  pobj[2] = AE_STR * plan->distance * sr + plan->distance;

  return;
}

//------------------------------------------------------------------------------
//! Compute mean elements.
//! The first nine arguments are for mean longitudes of planets (Simon et al, 
//! 1994), with 0.047" subtracted from constant term for offset to DE403 origin.
//!
//! \param jd_tt The Julian date in TT.
//! \param arg The mean element arguments.
//------------------------------------------------------------------------------

void ae_gplan_mean_element(double jd_tt, double *arg) {
  double x, T, T2;

  // Time variables.  T is in Julian centuries.
  T = (jd_tt - 2451545.0) / 36525.0;
  T2 = T * T;

  // Mean longitudes of planets (Simon et al, 1994).
  // 0.047" subtracted from constant term for offset to DE403 origin.

  // Mercury.
  x = ae_gplan_mods3600(538101628.6889819 * T + 908103.213);
  x += (6.39e-6 * T - 0.0192789) * T2;
  arg[0] = AE_STR * x;

  // Venus.
  x = ae_gplan_mods3600(210664136.4335482 * T + 655127.236);
  x += (-6.27e-6 * T + 0.0059381) * T2;
  arg[1] = AE_STR * x;

  // Earth.
  x = ae_gplan_mods3600(129597742.283429 * T + 361679.198);
  x += (-5.23e-6 * T - 2.04411e-2) * T2;
  arg[2] = AE_STR * x;

  // Mars.
  x = ae_gplan_mods3600(68905077.493988 * T + 1279558.751);
  x += (-1.043e-5 * T + 0.0094264) * T2;
  arg[3] = AE_STR * x;

  // Jupiter.
  x = ae_gplan_mods3600(10925660.377991 * T + 123665.420);
  x += ((((-3.4e-10 * T + 5.91e-8) * T + 4.667e-6) * T + 5.706e-5) * T - 
           3.060378e-1) * T2;
  arg[4] = AE_STR * x;

  // Saturn.
  x = ae_gplan_mods3600(4399609.855372 * T + 180278.752);
  x += ((((8.3e-10 * T - 1.452e-7) * T - 1.1484e-5) * T - 1.6618e-4) * T + 
           7.561614E-1) * T2;
  arg[5] = AE_STR * x;

  // Uranus.
  x = ae_gplan_mods3600(1542481.193933 * T + 1130597.971) + 
      (0.00002156 * T - 0.0175083)*T2;
  arg[6] = AE_STR * x;

  // Neptune.
  x = ae_gplan_mods3600(786550.320744 * T + 1095655.149) + 
      (-0.00000895 * T + 0.0021103) * T2;
  arg[7] = AE_STR * x;

  // Copied from cmoon.c, DE404 version. (?)
  // Mean elongation of moon = D.
  x = ae_gplan_mods3600(1.6029616009939659e+09 * T + 1.0722612202445078e+06);
  x += (((((-3.207663637426e-013 * T + 2.555243317839e-011) * T + 
            2.560078201452e-009) * T - 3.702060118571e-005) * T + 
            6.9492746836058421e-03) * T - 6.7352202374457519e+00) * T2;
  arg[9] = AE_STR * x;

  // Mean distance of moon from its ascending node = F.
  x = ae_gplan_mods3600(1.7395272628437717e+09 * T + 3.3577951412884740e+05);
  x += (((((4.474984866301e-013 * T + 4.189032191814e-011) * T - 
            2.790392351314e-009) * T - 2.165750777942e-006) * T - 
            7.5311878482337989e-04) * T - 1.3117809789650071e+01) * T2;
  arg[10] = AE_STR * x;

  // Mean anomaly of sun = l' (J. Laskar)
  x = ae_gplan_mods3600(1.2959658102304320e+08 * T + 1.2871027407441526e+06);
  x += ((((((((1.62e-20 * T - 1.0390e-17 ) * T - 3.83508e-15 ) * T + 
               4.237343e-13 ) * T + 8.8555011e-11 ) * T - 4.77258489e-8 ) * T -
               1.1297037031e-5 ) * T + 8.7473717367324703e-05) * T -
               5.5281306421783094e-01) * T2;
  arg[11] = AE_STR * x;

  // Mean anomaly of moon = l.
  x = ae_gplan_mods3600(1.7179159228846793e+09 * T + 4.8586817465825332e+05);
  x += (((((-1.755312760154e-012 * T + 3.452144225877e-011) * T - 
            2.506365935364e-008) * T - 2.536291235258e-004) * T + 
            5.2099641302735818e-02) * T + 3.1501359071894147e+01) * T2;
  arg[12] = AE_STR * x;

  // Mean longitude of moon, re mean ecliptic and equinox of date = L.
  x = ae_gplan_mods3600(1.7325643720442266e+09 * T + 7.8593980921052420e+05);
  x += (((((7.200592540556e-014 * T + 2.235210987108e-010) * T - 
            1.024222633731e-008) * T - 6.073960534117e-005) * T + 
            6.9017248528380490e-03) * T - 5.6550460027471399e+00) * T2;
  arg[13] = AE_STR * x;

  // Precession of the equinox.
  x = (((((((((-8.66e-20 * T - 4.759e-17) * T + 2.424e-15) * T + 
               1.3095e-12) * T + 1.7451e-10) * T - 1.8055e-8) * T - 
               0.0000235316) * T + 0.000076) * T + 1.105414) * T + 
               5028.791959) * T;
  // Moon's longitude re fixed AE_J2000 equinox.
  // arg[13] -= x;

  // Free librations.  Longitudinal libration 2.891725 years.
  x = ae_gplan_mods3600(4.48175409e7 * T + 8.060457e5);
  arg[14] = AE_STR * x;

  // Libration P, 24.2 years.
  x = ae_gplan_mods3600(5.36486787e6 * T - 391702.8);
  arg[15] = AE_STR * x;

  // REMOVE???
  #if 0
    arg[16] = 0.0;
  #endif

  // Libration W, 74.7 years.
  x = ae_gplan_mods3600(1.73573e6 * T);
  arg[17] = AE_STR * x;
}


//------------------------------------------------------------------------------
//! Accumulate the sum of trigonometric series in three variables.
//! This is a generic program to accumulate sum of trigonometric series in 
//! three variables (e.g., longitude, latitude, radius) of the same list of 
//! arguments.
//!
//! \param jd_tt The Julian date in TT.
//! \param plan A table of planetary orbit data.
//! \param pobj For returning the polar position of the planet.
//! \param objnum The planet number.
//------------------------------------------------------------------------------

void ae_g3plan(double jd_tt, struct ae_plantbl_t *plan, double pobj[], 
               int objnum) {
  char *p;
  int i, j, k, m, n, k1, ip, np, nt;
  long *pl, *pb, *pr;
  double su, cu, sv, cv, T, t, sl, sb, sr, arg[AE_PLANTBL_N_HARMONIC];  

  ae_gplan_mean_element(jd_tt, arg);
  
  T = (jd_tt - AE_J2000) / plan->timescale;
  n = plan->maxargs;

  // Calculate sin( i*MM ), etc. for needed multiple angles.
  for (i = 0; i < n; i++) {
    if ((j = plan->max_harmonic[i]) > 0)
	  ae_gplan_sscc(i, arg[i], j);
  }

  // Point to start of table of arguments.
  p = plan->arg_tbl;
  // Point to tabulated cosine and sine amplitudes.
  pl = (long *)plan->lon_tbl;
  pb = (long *)plan->lat_tbl;
  pr = (long *)plan->rad_tbl;
  sl = 0.0;
  sb = 0.0;
  sr = 0.0;

  for (;;) {
    // Argument of sine and cosine.
    np = *p++;
    if (np < 0)
      break;
    if (np == 0) {
      nt = *p++;

      // "Longitude" polynomial (phi).
      cu = *pl++;
      for (ip = 0; ip < nt; ip++)
        cu = cu * T + *pl++;
      // sl +=  ae_gplan_mods3600(cu);
      sl += cu;

	  // "Latitude" polynomial (theta).
      cu = *pb++;
	  for (ip = 0; ip < nt; ip++)
        cu = cu * T + *pb++;
      sb += cu;

	  // "Radius" polynomial (psi).
      cu = *pr++;
      for (ip = 0; ip < nt; ip++)
        cu = cu * T + *pr++;
	  sr += cu;

      continue;
    }

    k1 = 0;
    cv = 0.0;
    sv = 0.0;
    for (ip = 0; ip < np; ip++) {
	  // Which harmonic.
	  j = *p++;
	  // Which planet.
	  m = *p++ - 1;
	  if (j) {
        // k = abs (j);
        if (j < 0)
          k = -j;
        else
          k = j;
        k -= 1;
        su = ae_gplan_ss[m][k];	// Sin(k*angle).
        if (j < 0)
          su = -su;
        cu = ae_gplan_cc[m][k];
        if (k1 == 0) { // Set first angle.
          sv = su;
          cv = cu;
          k1 = 1;
        }
        else { // Combine angles.
          t = su * cv + cu * sv;
          cv = cu * cv - su * sv;
          sv = t;
        }
      }
    }
   
    // Highest power of T.
    nt = *p++;

    // "Longitude".
    cu = *pl++;
    su = *pl++;
    for (ip = 0; ip < nt; ip++)	{
      cu = cu * T + *pl++;
      su = su * T + *pl++;
    }
    sl += cu * cv + su * sv;
    
    // "Latitiude".
    cu = *pb++;
    su = *pb++;
    for (ip = 0; ip < nt; ip++)	{
      cu = cu * T + *pb++;
      su = su * T + *pb++;
    }
    sb += cu * cv + su * sv;

    // "Radius".
    cu = *pr++;
    su = *pr++;
    for (ip = 0; ip < nt; ip++) {
      cu = cu * T + *pr++;
      su = su * T + *pr++;
    }
    sr += cu * cv + su * sv;
  }

  t = plan->trunclvl;
  pobj[0] = arg[objnum - 1] + AE_STR * t * sl;
  pobj[1] = AE_STR * t * sb;
  pobj[2] = plan->distance * (1.0 + AE_STR * t * sr);
  
  return;
}


//------------------------------------------------------------------------------
//! Accumulate the sum of trigonometric series in two variables.
//! This is a generic routine to accumulate sum of trigonometric series in two 
//! variables (e.g., longitude, radius) of the same list of arguments.
//!
//! \param jd_tt The Julian date in TT.
//! \param plan A table of planetary orbit data.
//! \param pobj For returning the polar position of the planet.
//! \param lp_equinox For returning the lp_equinox (?); set to NULL if it is not
//!                   not needed
//------------------------------------------------------------------------------

void ae_g2plan(double jd_tt, struct ae_plantbl_t *plan, double pobj[], 
               double *lp_equinox) {
  char *p;
  int i, j, k, m, n, k1, ip, np, nt;
  long *pl, *pr;
  double su, cu, sv, cv, T, t, sl, sr, arg[AE_PLANTBL_N_HARMONIC];

  ae_gplan_mean_element(jd_tt, arg);
  if (lp_equinox != NULL)
    *lp_equinox = arg[13] / AE_STR;

  T = (jd_tt - AE_J2000) / plan->timescale;
  n = plan->maxargs;
  
  // Calculate sin(i * MM), etc. for needed multiple angles.
  for (i = 0; i < n; i++) {
    if ((j = plan->max_harmonic[i]) > 0)
      ae_gplan_sscc(i, arg[i], j);
  }

  // Point to start of table of arguments.
  p = plan->arg_tbl;
  // Point to tabulated cosine and sine amplitudes.
  pl = (long *)plan->lon_tbl;
  pr = (long *)plan->rad_tbl;
  sl = 0.0;
  sr = 0.0;

  for (;;) {
    // Argument of sine and cosine.
    np = *p++;
    
    if (np < 0)
      break;
    if (np == 0) {
      nt = *p++;

      // "Longitude" polynomial.
      cu = *pl++;
      for (ip = 0; ip < nt; ip++)
	      cu = cu * T + *pl++;
	  // sl +=  ae_gplan_mods3600 (cu);
      sl += cu;

	  // "Radius" polynomial.
      cu = *pr++;
      for (ip = 0; ip < nt; ip++)
        cu = cu * T + *pr++;
      sr += cu;

      continue;
    }

    k1 = 0;
    cv = 0.0;
    sv = 0.0;
    for (ip = 0; ip < np; ip++)	{
      // Which harmonic.
      j = *p++;
      // Which planet.
      m = *p++ - 1;
      if (j) {
	    // k = abs (j);
        if (j < 0)
          k = -j;
        else
          k = j;
        k -= 1;
        su = ae_gplan_ss[m][k];	// Sin(k * angle).
        if (j < 0)
          su = -su;
        cu = ae_gplan_cc[m][k];
        if (k1 == 0) { // Set first angle.
          sv = su;
          cv = cu;
          k1 = 1;
        }
        else { // Combine angles.
          t = su * cv + cu * sv;
          cv = cu * cv - su * sv;
          sv = t;
        }
      }
	}

    // Highest power of T.
    nt = *p++;
   
    // Longitude.
    cu = *pl++;
    su = *pl++;
    for (ip = 0; ip < nt; ip++)	{
      cu = cu * T + *pl++;
	  su = su * T + *pl++;
    }
    sl += cu * cv + su * sv;
   
    // Radius.
    cu = *pr++;
    su = *pr++;
    for (ip = 0; ip < nt; ip++)	{
      cu = cu * T + *pr++;
      su = su * T + *pr++;
    }
    sr += cu * cv + su * sv;
  }

  t = plan->trunclvl;
  pobj[0] = t * sl;
  pobj[2] = t * sr;
  
  return;
}


//------------------------------------------------------------------------------
//! Accumulate the sum of trigonometric series in one variable.
//! This is a generic routine to accumulate sum of trigonometric series in one 
//! variable of the same list of arguments.
//!
//! \param jd_tt The Julian date in TT.
//! \param plan A table of planetary orbit data.
//!
//! \return The accumlated sum.
//------------------------------------------------------------------------------

double ae_g1plan(double jd_tt, struct ae_plantbl_t *plan) {
  char *p;
  int i, j, k, m, k1, ip, np, nt;
  long *pl;
  double su, cu, sv, cv, T, t, sl, arg[AE_PLANTBL_N_HARMONIC];  

  T = (jd_tt - AE_J2000) / plan->timescale;
  ae_gplan_mean_element(jd_tt, arg);

  // Calculate sin(i * MM), etc. for needed multiple angles.
  for (i = 0; i < AE_PLANTBL_N_HARMONIC; i++) {
    if ((j = plan->max_harmonic[i]) > 0)
      ae_gplan_sscc(i, arg[i], j);
  }

  // Point to start of table of arguments.
  p = plan->arg_tbl;
  // Point to tabulated cosine and sine amplitudes.
  pl = (long *)plan->lon_tbl;
  sl = 0.0;

  for (;;) {
    // Argument of sine and cosine.
    np = *p++;
    
    if (np < 0)
      break;
    if (np == 0) {
      nt = *p++;
      cu = *pl++;
      for (ip = 0; ip < nt; ip++)
        cu = cu * T + *pl++;
      // sl += ae_gplan_mods3600(cu);
      sl += cu;

      continue;
	}

    k1 = 0;
    cv = 0.0;
    sv = 0.0;
    for (ip = 0; ip < np; ip++)	{
      // Which harmonic.
      j = *p++;
      // Which planet.
      m = *p++ - 1;
      if (j) {
        // k = abs (j);
        if (j < 0)
          k = -j;
        else
          k = j;
        k -= 1;
        su = ae_gplan_ss[m][k];	// Sin(k * angle).
        if (j < 0)
          su = -su;
        cu = ae_gplan_cc[m][k];
        if (k1 == 0) { // Set first angle.
          sv = su;
          cv = cu;
          k1 = 1;
        }
        else { // Combine angles.
          t = su * cv + cu * sv;
          cv = cu * cv - su * sv;
          sv = t;
        }
      }
    }

    // Highest power of T.
    nt = *p++;

    // Cosine and sine coefficients.
    cu = *pl++;
    su = *pl++;
    for (ip = 0; ip < nt; ip++) {
      cu = cu * T + *pl++;
      su = su * T + *pl++;
    }
    sl += cu * cv + su * sv;
  }
 
  return plan->trunclvl * sl;
}


//------------------------------------------------------------------------------
//! Compute geocentric moon position.
//!
//! \param jd_tt The Julian date in TT.
//! \param rect For returning the position in rectangular coordinates.
//! \param pol For returning the position in polar coordinates.
//------------------------------------------------------------------------------

void ae_gmoon(double jd_tt, double rect[], double pol[]) {
  double x, cosB, sinB, cosL, sinL, lp_equinox, eps, coseps, sineps;

  ae_g2plan(jd_tt, &ae_mlr404, pol, &lp_equinox);
  x = pol[0];
  x += lp_equinox;
  if (x < -6.48e5)
    x += 1.296e6;
  if (x > 6.48e5)
    x -= 1.296e6;
  pol[0] = AE_STR * x;

  x = ae_g1plan(jd_tt, &ae_mlat404);
  pol[1] = AE_STR * x;
  x = (1.0 + AE_STR * pol[2]) * ae_mlr404.distance;
  pol[2] = x;

  // Convert ecliptic polar to equatorial rectangular coordinates.
  eps = ae_epsilon(jd_tt);
  coseps = cos(AE_STR * eps);
  sineps = sin(AE_STR * eps);
  cosB = cos(pol[1]);
  sinB = sin(pol[1]);
  cosL = cos(pol[0]);
  sinL = sin(pol[0]);
  rect[0] = cosB * cosL * x;
  rect[1] = (coseps * cosB * sinL - sineps * sinB) * x;
  rect[2] = (sineps * cosB * sinL + coseps * sinB) * x;
}
