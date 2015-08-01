//! \file init.c
//! Define default orbit structures and do other initialisation.
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

#include <stdlib.h>
#include <string.h>

#include "aephem.h"

//! Default orbital parameters for Mercury.
struct ae_orbit_t ae_orb_mercury = {
  "Mercury",
  2446800.5,    // January 5.0, 1987.
  7.0048,
  48.177,
  29.074,
  0.387098,
  4.09236,
  0.205628,
  198.7199,
  2446800.5,
  &ae_mer404,
  0.0,
  0.0,
  0.0
};

//! Default orbital parameters for Venus.
struct ae_orbit_t ae_orb_venus = {
  "Venus",
  2446800.5,
  3.3946,
  76.561,
  54.889,
  0.723329,
  1.60214,
  0.006757,
  9.0369,
  2446800.5,
  &ae_ven404,
  0.0,
  0.0,
  0.0
};

//! Default orbital parameters for Earth.
struct ae_orbit_t ae_orb_earth = {
  "Earth",
  2446800.5,
  0.0,
  0.0,
  102.884,
  0.999999,
  0.985611,
  0.016713,
  1.1791,
  2446800.5,
  &ae_ear404,
  0.0,
  0.0,
  0.0
};

//! Default orbital parameters for Mars.
struct ae_orbit_t ae_orb_mars = {
  "Mars",
  2446800.5,
  1.8498,
  49.457,
  286.343,
  1.523710,
  0.524023,
  0.093472,
  53.1893,
  2446800.5,
  &ae_mar404,
  0.0,
  0.0,
  0.0
};

//! Default orbital parameters for Jupiter.
struct ae_orbit_t ae_orb_jupiter = {
  "Jupiter",
  2446800.5,
  1.3051,
  100.358,
  275.129,
  5.20265,
  0.0830948,
  0.048100,
  344.5086,
  2446800.5,
  &ae_jup404,
  0.0,
  0.0,
  0.0
};

//! Default orbital parameters for Saturn.
struct ae_orbit_t ae_orb_saturn = {
  "Saturn",
  2446800.5,
  2.4858,
  113.555,
  337.969,
  9.54050,
  0.0334510,
  0.052786,
  159.6327,
  2446800.5,
  &ae_sat404,
  0.0,
  0.0,
  0.0
};

//! Default orbital parameters for Uranus.
struct ae_orbit_t ae_orb_uranus = {
  "Uranus",
  2446800.5,
  0.7738,
  73.994,
  98.746,
  19.2233,
  0.0116943,
  0.045682,
  84.8516,
  2446800.5,
  &ae_ura404,
  0.0,
  0.0,
  0.0
};

//! Default orbital parameters for Neptune.
struct ae_orbit_t ae_orb_neptune = {
  "Neptune",
  2446800.5,
  1.7697,
  131.677,
  250.623,
  30.1631,
  0.00594978,
  0.009019,
  254.2568,
  2446800.5,
  &ae_nep404,
  0.0,
  0.0,
  0.0
};

//! Default orbital parameters for Pluto.
struct ae_orbit_t ae_orb_pluto = {
  "Pluto",
  2446640.5,
  17.1346,
  110.204,
  114.21,
  39.4633,
  0.00397570,
  0.248662,
  355.0554,
  2446640.5,
  &ae_plu404,
  0.0,
  0.0,
  0.0
};

//! Pointers to the planetary orbit objects.
//! The indices correspond to #ae_ss_bodies_t.
struct ae_orbit_t *ae_orb_planet[] = {
  &ae_orb_mercury,
  &ae_orb_venus,
  NULL,
  &ae_orb_mars,
  &ae_orb_jupiter,
  &ae_orb_saturn,
  &ae_orb_uranus,
  &ae_orb_neptune,
  &ae_orb_pluto,
  NULL,
  NULL,
  NULL,
  NULL,
  &ae_orb_earth,
  NULL,
  NULL
};

//! Pointers to the planetary physical information objects.
//! The indices correspond to #ae_ss_bodies_t.
struct ae_physical_t *ae_phys_planet[] = {
  &ae_phys_mercury,
  &ae_phys_venus,
  NULL,
  &ae_phys_mars,
  &ae_phys_jupiter,
  &ae_phys_saturn,
  &ae_phys_uranus,
  &ae_phys_neptune,
  &ae_phys_pluto,
  NULL,
  &ae_phys_sun,
  NULL,
  NULL,
  &ae_phys_earth,
  NULL,
  NULL
};

//! The names of the solar system bodies.
//! The indices correspond to #ae_ss_bodies_t.
const char *ae_ss_name[] = {
  "Mercury",
  "Venus",
  "Earth-moon barycentre",
  "Mars",
  "Jupiter",
  "Saturn",
  "Uranus",
  "Neptune",
  "Pluto",
  "moon, relative to earth-moon barycentre",
  "sun",
  "nutation",
  "libration",
  "Earth",
  "moon",
  "Solar System barycentre",
};

//! Human readable strings describing return codes.
const char *ae_retcode_name[] = {
  "function excecuted normally",
  "could not open file",
  "unexpected end-of-file encountered",
  "unknown read error occured",
  "bad JPL file header",
  "bad JPL file header (KSIZE != 2 * N_COEFF)",
  "bad JPL file header (could not find a required group)",
  "bad JPL file header (differing numbers of constants and values)",
  "bad JPL file header (bad date range)",
  "unrecognised JPL file format (expecting ASCII)",
  "unrecognised JPL file format (expecting binary)",
  "bad JPL data block",
  "JPL data file seems to be corrupt",
  "Julian date outside of JPL file date range",
  "unable to parse expected float in JPL ASCII file",
  "bad Chebyshev parameter when interpolating data",
  "unknown object number",
  "object not found in catalogue file",
  "function requires aephem to be built with cfitsio support",
  "a required FITS keyword does not exist",
  "a required FITS keyword is corrupt",
  "unknown WCS projection type",
  "an error occurred in a call to a libcfitsio function",
  "could not compute the native pole of the projection",
  "a matrix is not invertible",
  "functionality under construction",
  NULL
};
//------------------------------------------------------------------------------
//! For returning human readable return codes.
//!
//! \param retcode A return code.  The absolute value will be used, so either
//!                a negative or positive number may be passed.
//!
//! \return A string describing the return code.
//------------------------------------------------------------------------------

const char *ae_retcode_string(int retcode) {
  static char *ret;
  int i;
  static int n_ret = -1;

  if (n_ret < 0) {
    i = strlen("unknown return code (9223372036854775808) pad");
    for (n_ret = 0; ae_retcode_name[n_ret] != NULL; n_ret++) {
      if (strlen(ae_retcode_name[n_ret]) > i)
        i = strlen(ae_retcode_name[n_ret]);
    }
    ret = (char *)malloc((i + 1) * sizeof(char));
  }

  retcode = abs(retcode);
  if (retcode < 0 || retcode >= n_ret)
    sprintf(ret, "unknown return code (%d)", retcode);
  else
    strcpy(ret, ae_retcode_name[retcode]);

  return ret;
}
