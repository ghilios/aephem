//==============================================================================
// AEPHEM - an astronomical ephemeris and reduction library.
// Copyright 2012 Adam Hincks.
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
#include <stdlib.h>
#include <string.h>

#include <aephem.h>

void usage() {
  fprintf(stderr, "Usage: test <FILE>\n");
  fprintf(stderr, "See README file.\n");
  exit(0);
}

int main(int argc, char *argv[]) {
  char target_name[128], tmpstr[128], *ptr, minus[128];
  int i, n, target_id, start_ok, jd_pos, radec_pos, azel_pos, subpnt_pos;
  long rec_start;
  double *jd, *ra_in, *dec_in, *az_in, *alt_in, *subpnt_lon_in, *subpnt_lat_in;
  double *ra, *dec, *az, *alt, *subpnt_lon, *subpnt_lat, jd_tt, dist, last;
  double jd_ut1, longitude, latitude, altitude, d1, m1, s1, d2, m2, s2;
  FILE *fp;

  if (argc != 2)
    usage();

  // Read in HORIZONS file.
  if ((fp = fopen(argv[1], "r")) == NULL) {
    fprintf(stderr, "Could not open %s for reading.\n", argv[1]);
    exit(0);
  }

  target_id = -1;
  longitude = 1000;
  latitude = 1000;
  altitude = 0;
  start_ok = 0;
  jd_pos = -1;
  radec_pos = -1;
  azel_pos = -1;
  subpnt_pos = -1;
  while (fgets(tmpstr, 127, fp) != NULL) {
    if (target_id < 0) {
      // First get the target body name.
      if (!strncmp(tmpstr, "Target body name:", 17)) {
        if (sscanf(tmpstr, "%*s %*s %*s %127s", target_name) != 1) {
          fprintf(stderr, "Error parsing file:  could not read target body "
                          "name.\n");
          exit(0);
        }
       
        for (target_id = 0; target_id < AE_N_SS_BODIES; target_id++) {
          if (!strcasecmp(target_name, ae_ss_name[target_id]))
            break;
        }
        if (target_id == AE_N_SS_BODIES) {
          fprintf(stderr, "Error parsing file:  did not recognise target "
                          "name \"%s\".", target_name);
          exit(0);
        }
        printf("Target is %s.\n", ae_ss_name[target_id]);
      }
    }
    else if (latitude > 360) {
      // Now get the observer position.
      if (!strncmp(tmpstr, "Center geodetic :", 17)) {
        if (sscanf(tmpstr, "%*s %*s %*s %lf,%lf,%lf", 
                           &longitude, &latitude, &altitude) != 3) {
          fprintf(stderr, "Error parsing file:  could not read observer "
                          "coordinates.\n");
          exit(0);
        }
        printf("Observer position is:  lat " ae_dms_fmt ", lon " ae_dms_fmt ", "
               "alt %g km\n", ae_dms_arg(longitude), ae_dms_arg(latitude),
               altitude);
      }
    }
    else if (!start_ok) {
      if (!strncmp("Date_", tmpstr, 5)) {
        if ((ptr = strstr(tmpstr, "JDUT")) != NULL)
          jd_pos = 0;
        if ((ptr = strstr(tmpstr, "R.A.__(a-apparent)")) != NULL)
          radec_pos = (int)(ptr - tmpstr);
        if ((ptr = strstr(tmpstr, "Azi_(a-appr)")) != NULL)
          azel_pos = (int)(ptr - tmpstr);
        if ((ptr = strstr(tmpstr, "Ob-lon")) != NULL)
          subpnt_pos = (int)(ptr - tmpstr);
        if (strstr(tmpstr, "Ob-lat") == NULL)
          subpnt_pos = -1;

        if (jd_pos >= 0 && radec_pos >= 0 && azel_pos >= 0 && subpnt_pos >=0)
          start_ok = 1;
        else {
          fprintf(stderr, "Could not parse file.  The following columns are "
                          "required:\n"
                          "  JD date, apparent ra/dec, az/el, ob-lon, "
                          "ob-lat\n");
          exit(0);
        }
      }
    }
    else {
      if (!strncmp("$$SOE", tmpstr, 5)) {
        start_ok = 2;
        break;
      }
    }
  }

  if (start_ok == 0) {
    fprintf(stderr, "Could not parse file.  The following columns are "
                    "required:\n");
    fprintf(stderr, "  JD date, apparent ra/dec, az/el, ob-lon, ob-lat\n");
    exit(0);
  }
  else if (start_ok == 1) {
    fprintf(stderr, "Could not find record start (\"$$SOE\").");
    exit(0);
  }

  // Count number of records.
  rec_start = ftell(fp);
  for (n = 0, i = 0; fgets(tmpstr, 127, fp); n++) {
    if (!strncmp("$$EOE", tmpstr, 5)) {
      i = 1;
      break;
    }
  }
  if (!i) {
    fprintf(stderr, "Could not find record end (\"$$SOE\").");
    exit(0);
  }
  fseek(fp, rec_start, SEEK_SET);
  
  printf("%d records found.\n", n);
  jd = (double *)malloc(n * sizeof(double));
  ra_in = (double *)malloc(n * sizeof(double));
  ra = (double *)malloc(n * sizeof(double));
  dec_in = (double *)malloc(n * sizeof(double));
  dec = (double *)malloc(n * sizeof(double));
  az_in = (double *)malloc(n * sizeof(double));
  az = (double *)malloc(n * sizeof(double));
  alt_in = (double *)malloc(n * sizeof(double));
  alt = (double *)malloc(n * sizeof(double));
  subpnt_lon_in = (double *)malloc(n * sizeof(double));
  subpnt_lon = (double *)malloc(n * sizeof(double));
  subpnt_lat_in = (double *)malloc(n * sizeof(double));
  subpnt_lat = (double *)malloc(n * sizeof(double));

  // Read in the HORIZONS values.
  for (i = 0; i < n; i++) {
    if (fgets(tmpstr, 127, fp) == NULL) {
      fprintf(stderr, "Wow, congratulations on triggering this error!\n");
      exit(0);
    }

    sscanf(tmpstr, "%lf", &jd[i]);
    sscanf(tmpstr + radec_pos, "%lf %lf %lf %lf %lf %lf", &d1, &m1, &s1, &d2,
                               &m2, &s2);
    ra_in[i] = (d1 + m1 / 60.0 + s1 / 3600.0) * 15.0;
    dec_in[i] = fabs(d2) + m2 / 60.0 + s2 / 3600.0;
    sscanf(tmpstr + radec_pos, "%*s %*s %*s %s", minus);
    if (minus[0] == '-')
      dec_in[i] *= -1.0;
    sscanf(tmpstr + azel_pos, "%lf %lf", &az_in[i], &alt_in[i]);
    sscanf(tmpstr + subpnt_pos, "%lf %lf", &subpnt_lon_in[i], 
                                &subpnt_lat_in[i]);
  }
  fclose(fp);

  // Now calculate the values with AEPHEM.
  for (i = 0; i < n; i++) {
    jd_ut1 = jd[i] + ae_dut1(jd[i]) * AE_D_PER_S;     // UT1 measure.
    jd_tt = jd_ut1 + ae_delta_t(jd_ut1) * AE_D_PER_S;   // TT measure.
    last = aes_last(jd_ut1, longitude);   // Local apparent sidereal time.
 
    // Get the geocentric ra/dec of the object.
    ae_geocentric_from_orbit(jd_tt, &ae_orb_earth, ae_orb_planet[target_id],
                             &ra[i], &dec[i], &dist);
 
    // Get the apparent ra/dec of the obect (no atmospheric model).
    aes_topocentric(jd_ut1, latitude, longitude, dist, &ra[i], &dec[i]);
 
    // Convert to alt/az.
    ae_radec_to_altaz(last, latitude, ra[i], dec[i], &alt[i], &az[i]);
  
    // Get the subobserver point.
    aes_subobs_point(jd_ut1, &ae_orb_earth, ae_orb_planet[target_id],
                     ae_phys_planet[target_id], &subpnt_lat[i],
                     &subpnt_lon[i], NULL);
  }

  printf("\nRESULTS:\n\n");
  printf("JD                  ??ra(\")  ??dec(\")  ??alt(\")   ??az(\")"
         " ??spnt_lon(') ??spnt_lat(')\n");
  printf("---------------------------------------------------------"
         "-----------------------\n");
  for (i = 0; i < n; i++) {
    printf("%17.9f %8.3f %8.3f %8.3f %8.3f     %8.3f     %8.3f\n", jd[i],
           (ra_in[i] - ra[i]) * 3600.0, (dec_in[i] - dec[i]) * 3600.0,
           (alt_in[i] - alt[i]) * 3600.0, (az_in[i] - az[i]) * 3600.0,
           (subpnt_lon_in[i] - subpnt_lon[i]) * 60.0,
           (subpnt_lat_in[i] - subpnt_lat[i]) * 60.0);
  }


  return 0;
}
