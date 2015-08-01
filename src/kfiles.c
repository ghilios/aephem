//! \file kfiles.c
//! Read orbital element files.
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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "aephem.h"

//------------------------------------------------------------------------------
//! Parse orbital elements from a catalogue file.
//! The catalogue file should be a columnated ASCII file.  Lines starting with
//! a hash ('#') are treated as comments.  The columns are:
//! - epoch of elements;
//! - inclination, in degrees;
//! - longitude of the ascending node, in degrees;
//! - argument of the perihelion, in degrees;
//! - mean distance, in AU;
//! - daily motion, in AU per day;
//! - eccentricity;
//! - mean anomaly, in degrees;
//! - epoch of equinox and ecliptic;
//! - name.
//!
//! \param path The location of the catalogue file.
//! \param name The name of the object, as it appears in the catalogue file.
//! \param orb For returning the orbital elements.  If the routine fails,
//!            the information in \p el will be undefined.
//!
//! \return 0 on successful load; on failure, a negative error code from
//!         #ae_retcode_t.
//------------------------------------------------------------------------------

int ae_read_orbit_from_cat(const char *path, const char *name,
                           struct ae_orbit_t *orb) {
  char line[512];
  FILE *fp;

  if ((fp = fopen(path, "r")) == NULL)
    return -AE_RET_BAD_PATH;

  while (fgets(line, 511, fp) != NULL) {
    if (line[0] == '#')
      continue;

    if (sscanf(line, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %31s",
               &orb->epoch, &orb->i, &orb->W, &orb->w, &orb->a, &orb->dm, 
               &orb->ecc, &orb->M, &orb->equinox, orb->name) != 10)
      continue;
    orb->name[31] = '\0';

    if (!strcmp(name, orb->name)) {
      // We found the object.  Clear out the rest of the orbit elements and
      // return.
      orb->ptable = NULL;
      orb->L = 0.0;
      orb->r = 0.0;
      orb->plat = 0.0;

      fclose(fp);

      return 0;
    }
  }

  // Unsuccessful search.
  fclose(fp);
  return -AE_RET_OBJ_NOT_FOUND;
}


//------------------------------------------------------------------------------
//! Parse star positional data from a catalogue file.
//! The catalogue file should be a columnated ASCII file.  Lines starting with
//! a hash ('#') are treated as comments.  The columns are:
//! - epoch
//! - ra hours
//! - ra minutes
//! - ra second
//! - dec degrees
//! - dec minutes
//! - dec seconds
//! - mu ra (in seconds/century)
//! - mu dec (in seconds of arc/century)
//! - v (in km/s) 
//! - px (in seconds of arc)
//! - mag 
//! - name
//!
//! \param path The location of the catalogue file.
//! \param name The name of the star, as it appears in the catalogue file.
//! \param star For returning the star position data.  If the routine fails,
//!             the information in \p el will be undefined.
//!
//! \return 0 on successful load; on failure, a negative error code from
//!         #ae_retcode_t.
//------------------------------------------------------------------------------

int ae_read_star_from_cat(const char *path, const char *name, 
                          struct ae_star_t *star) { 
  char line[512], *p;
  int i, sign;
  double rh, rm, rs, dd, dm, ds, x, z;
  FILE *fp;

  if ((fp = fopen(path, "r")) == NULL)
    return -AE_RET_BAD_PATH;

  while (fgets(line, 511, fp) != NULL) {
    if (line[0] == '#')
      continue;

    if (sscanf(line, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %31s",
          &x, &rh, &rm, &rs, &dd, &dm, &ds, &star->mura, &star->mudec, 
          &star->v, &star->px, &star->mag, star->name) != 13)
      continue;
    star->name[31] = '\0';

    if (!strcmp(name, star->name)) {
      // We found the object.  Format its data.

      // Parse the epoch.
      if (x == 2000.0)
        x = AE_J2000;
      else if (x == 1950.0)
        x = AE_B1950;
      else if (x == 1900.0)
        x = AE_J1900;
      else
        x = AE_J2000 + 365.25 * (x - 2000.0);
      star->epoch = x;

      // Parse the right ascension.
      star->ra = 15.0 * (rh + rm / 60.0 + rs / 3600.0);

      // Parse the declination.
      sign = 1;
      if ((dd < 0.0) || (dm < 0.0) || (ds < 0.0))
        sign = -1;
      z = fabs(dd) + fabs(dm) / 60.0 + fabs(ds) / 3600.0;

      if (dd == 0.0) {
        // Scan the text for possible minus sign in front of declination 0.
        p = line;
        // Skip over 4 fields.
        for (i = 0; i < 4; i++) {
          for (; *p == ' ' || *p == '\t'; p++);
          for (; *p != ' ' && *p != '\t'; p++);
        }
        for (; *p == ' ' || *p == '\t'; p++);
        p--;
        if (*p == '-')
          sign = -1;
      }
      if (sign < 0)
        z = -z;
      star->dec = z;

      // Put proper motion into "/century.
      star->mura *= 15.0; // s/century -> "/century

      z = star->px;
      if (z < 1.0) {
        if (z <= 0.0)
          star->px = 0.0;
        else
          star->px = z;  // Assume px in seconds of arc.
      }
      else
        star->px = 1.0 / z;    // Parsecs -> seconds of arc.
  
      fclose(fp);
      return 0;
    }
  }

  // Unsuccessful search.
  fclose(fp);
  return -AE_RET_OBJ_NOT_FOUND;
}
