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
// This interactive program gets input from the user about which object he 
// would like to view, using which ephemerides and which time-spans. It then
// prints everything in a big ol' table.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <aephem.h>

int main(int argc, char *argv[]) {
  int i, j, n_tab, obj_num, dummy;
  char ob_type[15], cat_path[128], star_name[32], full_star_name[128];
  char name[128];
  double jd_tt, jd_ut, jd_ut1, tlong, glat, tlat, trho, height, y, m, d;
  double e_pos[3], e_vel[3], geo_ra, geo_dec, topo_ra, topo_dec, alt, az;
  double dist, last, tab_int, dut1, delta_t;
  struct ae_orbit_t *orb;
  struct ae_star_t star;
  struct ae_jpl_handle_t jh;

  printf("Enter the observer east longitude (degrees): ");
  dummy = scanf("%lf", &tlong);
  printf("Enter the observer north geodetic latitude (degrees): ");
  dummy = scanf("%lf", &glat);
  printf("Enter the observer altitude (m above sea level): ");
  dummy = scanf("%lf", &height);
  printf("Enter the starting UT date (JD 0 0 or year month fractional_day): ");
  dummy = scanf("%lf %lf %lf", &y, &m, &d);
  if (m == 0 && d == 0)
    jd_ut = y;
  else
    jd_ut = ae_cal_to_jd((long)y, (long)m, d);
  printf("Enter tabulation intervals (days) and number of tabulations: ");
  dummy = scanf("%lf %d", &tab_int, &n_tab);
  printf("\nHere is a list of options:\n");
  printf("  P . . . . . use built-in planetary ephemeris\n");
  printf("  J . . . . . use a JPL ephemeris file\n");
  printf("  O . . . . . use an aephem orbital elements file\n");
  printf("  C . . . . . use a star catalogue\n");
  printf("Please enter one of these options: ");
  dummy = scanf("%15s", ob_type);

  switch(ob_type[0]) {
    case 'P': case 'p':
      printf("Here is a list of available objects:\n");
      for (i = 0; i < AE_N_SS_BODIES; i++) {
        if (ae_orb_planet[i] != NULL)
          printf("  %2d . . . . %s\n", i, ae_ss_name[i]);
      }
      printf("Enter the object number: ");
      dummy = scanf("%d", &i);
      if (i < 0 || i >= AE_N_SS_BODIES || ae_orb_planet[i] == NULL) {
        printf("That is not a valid option.\n");
        exit(0);
      }
      orb = ae_orb_planet[i];
      strcpy(name, orb->name);

      break;

    case 'J': case 'j':
      printf("Enter the name of the JPL file: ");
      dummy = scanf("%127s", cat_path);

      if ((i = ae_jpl_init(cat_path, &jh))) {
        printf("Could not load JPL ephemeris file:  %s.\n", 
               ae_retcode_string(i));
        exit(0);
      }
      
      printf("Here is a list of JPL objects:\n");
      for (i = 0; i < AE_N_SS_BODIES; i++)
        printf("  %2d . . . . %s\n", i, ae_ss_name[i]);
      printf("Enter the object number: ");
      dummy = scanf("%d", &obj_num);
      if (obj_num < 0 || obj_num >= AE_N_SS_BODIES) {
        printf("That is not a valid option.\n");
        exit(0);
      }
      strcpy(name, ae_ss_name[obj_num]);
      break;

    case 'O': case 'o':
      printf("Enter the name of the orbital file: ");
      dummy = scanf("%127s", cat_path);
      printf("Enter the name of the orbit, as it appears in %s: ", cat_path);
      dummy = scanf("%31s", star_name);

      if ((i = ae_read_orbit_from_cat(cat_path, star_name, orb))) {
        printf("Could not load orbit:  %s.\n", ae_retcode_string(i));
        exit(0);
      }

      ob_type[0] = 'p';  // From now on, treat this like any other orbit.
      strcpy(name, orb->name);

      break;

    case 'C': case 'c':
      printf("Enter the name of the catalogue file: ");
      dummy = scanf("%127s", cat_path);
      printf("Enter the name of the object, as it appears in %s: ", cat_path);
      dummy = scanf("%31s", star_name);

      if ((i = ae_read_star_from_cat(cat_path, star_name, &star))) {
        printf("Could not load catalogue item: %s.\n", ae_retcode_string(i));
        exit(0);
      }

      break;
    
    default:
      printf("You must enter 'P', 'O' or 'C'.\n");
      exit(0);
  }

  // Get the geocentric latitude.
  ae_geocentric_lat(glat, height, &tlat, &trho);

  printf("\n");
  printf("Observer coordinates (geocentric): " ae_hms_fmt ", " ae_dms_fmt "\n", 
         ae_hms_arg(tlong), ae_dms_arg(tlat));
  printf("Observer R_{earth}: %g\n", trho);

  for (i = 0; i < n_tab; i++, jd_ut += tab_int) {
    // Get delta_t and DUT1.
    delta_t = ae_delta_t(jd_ut);
    dut1 = ae_dut1(jd_ut);
    
    jd_ut1 = jd_ut + dut1 / AE_S_PER_D;
    jd_tt = jd_ut1 + delta_t / AE_S_PER_D;

    switch (ob_type[0]) {
      case 'P': case 'p':
        // Do geocentric reduction.
        ae_geocentric_from_orbit(jd_tt, &ae_orb_earth, orb, &geo_ra, &geo_dec,
                                 &dist);
        break;
      case 'J': case 'j':
        // Do geocentric reduction.
        printf("%.10f\n", jd_tt);
        if ((j = ae_geocentric_from_jpl(&jh, jd_tt, obj_num, &geo_ra, &geo_dec,
                                        &dist))) {
          printf("Could not do geocentric reduction:  %s.", 
                 ae_retcode_string(j));
          exit(0);
        }

        break;
      case 'C': case 'c':
        // Do geocentric reduction.
        ae_kepler(jd_tt, &ae_orb_earth, e_pos);
        ae_v_orbit(jd_tt, &ae_orb_earth, e_vel);
        ae_geocentric_from_cat(jd_tt, e_pos, e_vel, &star, &geo_ra, &geo_dec);
        dist = 0;
        break;
    }

    // Do topocentric reduction.
    topo_ra = geo_ra;
    topo_dec = geo_dec;
    ae_topocentric(jd_tt, jd_ut1, tlat, glat, tlong, trho, NULL, NULL, dist, 
                   &topo_ra, &topo_dec);

    if (!i) {
      // Print header.  First get the constellation.
      switch (ob_type[0]) {
        case 'J': case 'j': case 'P': case 'p': 
          printf("Object is %s, located in %s.\n", name,
                 ae_constel_name[ae_coord_to_constel_index
                                 (topo_ra, topo_dec, jd_tt)] + 4);
          printf("+--------------+--------+----------+-----------------+-------"
                 "-----------+------------+-----------------+------------------"
                 "+------------------+------------------+\n");
          printf("| JD (UT1)     | dT (s) | DUT1 (s) | Geocentric RA   | "
                 "Geocentric Dec.  | Dist. (AU) | Topocentric RA  | "
                 "Topocentric Dec. | Azimuth          | Altitude         |\n");
          printf("+--------------+--------+----------+-----------------+-------"
                 "-----------+------------+-----------------+------------------"
                 "+------------------+------------------+\n");
          break;

        case 'C': case 'c':
          ae_cat_to_constel_index(star_name, full_star_name, 128);
          printf("Star is %s with vis. mag. of %g.\n", full_star_name, 
                 star.mag);
          printf("+--------------+--------+----------+-----------------+-------"
                 "-----------+-----------------+------------------+------------"
                 "------+------------------+\n");
          printf("| JD (UT1)     | dT (s) | DUT1 (s) | Geocentric RA   | "
                 "Geocentric Dec.  | Topocentric RA  | Topocentric Dec. | "
                 "Azimuth          | Altitude         |\n");
          printf("+--------------+--------+----------+-----------------+-------"
                 "-----------+-----------------+------------------+------------"
                 "------+------------------+\n");
          break;
      }
    }

    // Find local alt/az.
    last = aes_last(jd_ut1, tlong);
    ae_radec_to_altaz(last, glat, topo_ra, topo_dec, &alt, &az);

    // Print the results.
    switch (ob_type[0]) {
      case 'p': case 'P': case 'j': case 'J':
        printf("| %12.4f | %5.2f  | %5.2f    | " ae_hms_sfmt " | " ae_dms_sfmt 
               " | %6.4f     | " ae_hms_sfmt " | " ae_dms_sfmt  " | " 
               ae_dms_sfmt " | " ae_dms_sfmt " |\n", jd_ut1, delta_t, dut1, 
               ae_hms_arg(geo_ra), ae_dms_arg(geo_dec), dist, 
               ae_hms_arg(topo_ra), ae_dms_arg(topo_dec), ae_dms_arg(az), 
               ae_dms_arg(alt));
        break;
      case 'c': case 'C':
        printf("| %12.4f | %5.2f  | %5.2f    | " ae_hms_sfmt " | " ae_dms_sfmt 
               " | "  ae_hms_sfmt " | " ae_dms_sfmt  " | " ae_dms_sfmt " | " 
               ae_dms_sfmt " |\n", jd_ut1, delta_t, dut1, ae_hms_arg(geo_ra), 
               ae_dms_arg(geo_dec), ae_hms_arg(topo_ra), ae_dms_arg(topo_dec), 
               ae_dms_arg(az), ae_dms_arg(alt));
        break;
    }
  }
  switch (ob_type[0]) {
    case 'p': case 'P': case 'j': case 'J':
       printf("+--------------+--------+----------+-----------------+-------"
              "-----------+------------+-----------------+------------------"
              "+------------------+------------------+\n");
      break;
    case 'c': case 'C':
        printf("+--------------+--------+----------+-----------------+-------"
               "-----------+-----------------+------------------+------------"
               "------+------------------+\n");
      break;
  }
  
  return 0;
}
