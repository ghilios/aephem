//==============================================================================
// AEPHEM - an astronomical ephemeris and reductpann library.
// Copyright 2012 Adam Hincks.
//
// This file is part of AEPHEM.
//
// AEPHEM is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundatpann, either verspann 3 of the License, or
// (at your optpann) any later verspann.
//
// AEPHEM is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with AEPHEM.  If not, see <http://www.gnu.org/licenses/>.
//==============================================================================
// In this program, we first determine the sub-observer points, positions and 
// solid angles of all the planets at the current time. Then, to show more 
// functions, we define the physical parameters for Pan (the moon of Saturn) 
// and then find out some stuff about it at the current time.

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <aephem.h>

int main(int argc, char *argv[]) {
  int i;
  double jd_ut, jd_ut1, jd_tt, ra, dec, dist, sa, w, n[3], p[3], e[3], v_e[3];
  double j[3], light_t, lat, lon, f, a, b;
  struct ae_physical_t pan;

  // Get the current Julian Date in the required time definitions.
  jd_ut = ae_ctime_to_jd(time(NULL));
  jd_ut1 = jd_ut + ae_dut1(jd_ut) * AE_D_PER_S;
  jd_tt = jd_ut1 + ae_delta_t(jd_ut1) * AE_D_PER_S;
  printf("The current Julian Date in UT1 is %0.5f.\n\n", jd_ut1);

  // Get the planets' solid angles.
  printf("Here are the planets' data (all geocentric):\n");
  printf("  -------------------------------------------------------------------"
         "-------\n");
  printf("  Planet      R. Ascension      Declination   Sub-lat.(°)  "
         "S.-Ang.(arcsec^2)\n");
  printf("  -------------------------------------------------------------------"
         "-------\n");
  for (i = 0; i <= AE_SS_NEPTUNE; i++) {
    if (i == AE_SS_EMBARY)
      continue; // Don't do the E-M barycentre!

    aes_subobs_point(jd_ut1, &ae_orb_earth, ae_orb_planet[i],
                     ae_phys_planet[i], &lat, &lon, NULL);
    sa = aes_disc_solid_angle(jd_ut1, &ae_orb_earth, ae_orb_planet[i],
                              ae_phys_planet[i]);
    ae_geocentric_from_orbit(jd_tt, &ae_orb_earth, ae_orb_planet[i], &ra,
                             &dec, NULL);
    printf("  %-8s " ae_hms_sfmt " " ae_dms_sfmt "        %6.2f            "
           "%7.2f\n", ae_ss_name[i], ae_hms_arg(ra), ae_dms_arg(dec), lat,
           sa * AE_RTS * AE_RTS);
  }

  // Now do some heavy-lifting on Pan.  Define its physical parameters.
  pan.r_mean       =  17.2;
  pan.r_eq         =  10.4;
  pan.r_pole       =   9.1;
  pan.pole_ra      =  40.6;
  pan.pole_ra_t    =  -0.036;
  pan.ra_sin_term  = ae_phys_no_term;
  pan.pole_dec     =  83.5;
  pan.pole_dec_t   =  -0.004;
  pan.dec_cos_term = ae_phys_no_term;
  pan.w            =   48.8;
  pan.w_d          =  626.044;
  pan.w_d_sq       = 0;
  pan.w_sin_term   = ae_phys_no_term;

  // Get the unit vector pointing from us to Pan. We'll use Saturn as a proxy
  // since we don't have ephemerides for Pan itself. We could use
  // ae_geocentric_from_orbit() here, but we also need to know the light-travel
  // times. First, our coordinates.
  ae_kepler(jd_tt, &ae_orb_earth, e);
  ae_v_orbit(jd_tt, &ae_orb_earth, v_e);

  // Get the light-corrected coordinates of Saturn.
  ae_kepler(jd_tt, &ae_orb_saturn, p);
  light_t = 0;
  for (i = 0; i < 2; i++) {
    light_t = ae_light_t(e, p, 0);
    ae_kepler(jd_tt - light_t, &ae_orb_saturn, p);
  }

  // Now get the geocentric position of Saturn. (We won't be precise and reduce
  // to topocentric coordinates. After all, we are pretending that Pan is at
  // Saturn's exact position.)
  ae_geocentric(jd_tt, p, e, v_e, &ra, &dec, &dist);

  // Get the unit vector pointing from us to Pan. We pass a distance of -1 so
  // that the unit vector is pointing _towards_ us, not away from us.
  ae_polar_to_rect(ra, dec, -1.0, j);

  // Get the coordinates of Pan's pole and prime meridean, accounting for
  // light-travel time.
  ae_phys_pole(jd_ut1 - light_t, jd_tt, &pan, n, &w);

  // Get the flattening of the spheroid.
  f = ae_flattening(&pan);

  // Get the subobserver point.
  ae_subobs_point(j, n, w, f, ae_is_retrograde(&pan), &lat, &lon);

  // Get the semi-minor axis of the projected disc.
  a = pan.r_eq;
  b = ae_disc_semiminor(a, pan.r_pole, lat);
  
  // Get its solid angle.
  sa = ae_disc_solid_angle(a, b, dist);

  // Now print everything out.
  printf("\n");
  printf("Here are some data on Saturn's moon Pan:\n");
  printf("  Ra/dec:            " ae_hms_fmt ", " ae_dms_fmt "\n",
         ae_hms_arg(ra), ae_dms_arg(dec));
  printf("  Distance:          %.3f AU\n", dist);
  printf("  Light-travel time: %.2f min\n", light_t * 24.0 * 60.0);
  printf("  North-pole coords: [%.2f %.2f %.2f]\n", n[0], n[1], n[2]);
  printf("  Sub-obs. lon/lat:  %.1f, %.1f (deg)\n", lon, lat);
  printf("  Projected axes:    %.1f, %.1f (km)\n", a, b);
  printf("  Solid angle:       %.2g arcsec^2\n", sa * AE_RTS * AE_RTS);

  return 0;
}
