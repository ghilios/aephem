//! \file brightness.c
//! This file contains functions for computing planetary brightnesses.
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

#include <fcntl.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include "aephem.h"

//------------------------------------------------------------------------------
//! Quickly fill a #ae_stokes_t struct.
//!
//! \param s The Stokes structure to be filled.
//! \param i The I parameter.
//! \param q The Q parameter.
//! \param u The U parameter.
//! \param v The V parameter.
//------------------------------------------------------------------------------

void ae_set_stokes(struct ae_stokes_t *s, double i, double q, double u,
                       double v) {
  s->i = i;
  s->q = q;
  s->u = u;
  s->v = v;

  return;
}


//------------------------------------------------------------------------------
//! Copy Stokes parameter.
//!
//! \param origin The struct to copy from.
//! \param dest The struct to copy to.
//------------------------------------------------------------------------------

void ae_copy_stokes(const struct ae_stokes_t *origin, 
                    struct ae_stokes_t *dest) {
  dest->i = origin->i;
  dest->q = origin->q;
  dest->u = origin->u;
  dest->v = origin->v;

  return;
}


//------------------------------------------------------------------------------
//! Allocate the elements of an ae_saturn_components_t struct; set to zero.
//!
//! \param n_ring The number of rings.
//! \param comp The struct to be allocated.
//------------------------------------------------------------------------------

void ae_saturn_components_alloc(int n_ring, 
                                struct ae_saturn_components_t *comp) {
  int i;

  comp->n_ring_disc = (double *)malloc(n_ring * sizeof(double));
  comp->n_ring_back = (double *)malloc(n_ring * sizeof(double));
  comp->exp_tau = (double *)malloc(n_ring * sizeof(double));

  comp->n_disc = 0;
  for (i = 0; i < n_ring; i++) {
    comp->n_ring_disc[i] = 0;
    comp->n_ring_back[i] = 0;
    comp->exp_tau[i] = 0;
  }

  return;
}

//------------------------------------------------------------------------------
//! Free the elements of an allocated ae_saturn_components_t struct.
//!
//! \param comp The struct to be freed.
//------------------------------------------------------------------------------

void ae_saturn_components_free(struct ae_saturn_components_t *comp) {
  free(comp->n_ring_disc);
  free(comp->n_ring_back);
  free(comp->exp_tau);

  return;
}


//------------------------------------------------------------------------------
//! Figures out what fraction of a pixel falls within a disc.
//!
//! \param x The x-coordinate of the pixel's centre.
//! \param y The y-coordinate of the pixel's centre.
//! \param res The resolution of the pixel.
//! \param a The semi-major axis of the disc.
//! \param b the semi-minor axis of the disc.
//!
//! \return The fraction of the pixel falling in the disc, from 0 to 1.
//------------------------------------------------------------------------------

double ae_brightness_in_disc(double x, double y, double res, double a,
                             double b) {
  double a2, b2, xs, ys, x1, x2, y1, y2, hres;

  // Symmetry allows us to consider only one quadrant.
  x = fabs(x);
  y = fabs(y);

  // Useful values.
  hres = res / 2.0;
  a2 = a * a;
  b2 = b * b;

  // Get the intersection of the disc with each of the four sides of the pixel,
  // in units of the pixel size, i.e., with 0 being the near corner and 1 being
  // the far corner:  hence, the disc intersects the side if it is in the range
  // [0, 1].
  // Lower side, parallel to x-axis.
  ys = y - hres;
  x1 = a2 * (1.0 - ys * ys / b2);
  if (x1 < 0)
    x1 = -1.0;
  else {
    x1 = (sqrt(x1) - (x - hres)) / res;
  }

  // Upper side, parallel to y-axis.
  ys = y + hres;
  x2 = a2 * (1.0 - ys * ys / b2);
  if (x2 < 0)
    x2 = -1.0;
  else
    x2 = (sqrt(x2) - (x - hres)) / res;

  // Left side, parallel to y-axis.
  xs = x - hres;
  y1 = b2 * (1.0 - xs * xs / a2);
  if (y1 < 0)
    y1 = -1.0;
  else
    y1 = (sqrt(y1) - (y - hres)) / res;

  // Right side, parallel to y-axis.
  xs = x + hres;
  y2 = b2 * (1.0 - xs * xs / a2);
  if (y2 < 0)
    y2 = -1.0;
  else
    y2 = (sqrt(y2) - (y - hres)) / res;

  if (x1 < 0 && y1 < 0)
    return 0;   // The whole pixel falls outside the disc.
  else if (x2 > 1.0 && y2 > 1.0)
    return 1.0; // The whole pixel falls inside the disc.
  else if ((0 <= x1 && x1 <= 1.0) && (0 <= x2 && x2 <= 1.0))
    return x1 + 0.5 * (x2 - x1);  // The pixel is cut in half vertically.
  else if ((0 <= y1 && y1 <= 1.0) && (0 <= y2 && y2 <= 1.0))
    return y1 + 0.5 * (y2 - y1);  // The pixel is cut in half horizontally.
  else if ((0 <= x1 && x1 <= 1.0) && (0 <= y1 && y1 <= 1.0))
    return 0.5 * x1 * y1; // Only the bottom left corner is filled.
  else if ((0 <= x2 && x2 <= 1.0) && (0 <= y2 && y2 <= 1.0))
    return 1.0 - 0.5 * (1.0 - x2) * 
                       (1.0 - y2); // All but the upper right corner is filled.
  else if ((0 <= x1 && x1 <= 1.0) || (0 <= y1 && y1 <= 1.0))
    return 0.0; // The rare case when the the pixel is slightly over-osculated.
  else if ((0 <= x2 && x2 <= 1.0) || (0 <= y2 && y2 <= 1.0))
    return 1.0; // The rare case when the the pixel is almost completely filled.
  else {
    fprintf(stderr, "Aephem bug detected!!.\n");
    return 2.0; // We should never reach this point!
  }
}

//------------------------------------------------------------------------------
//! Figures out what fraction of a pixel falls inside a ring.
//!
//! \param x The x-coordinate of the pixel centre.
//! \param y The y-coordinate of the pixel centre.
//! \param res The resolution of the pixel.
//! \param sin_ring_inc The sine of the inclination of the ring.
//! \param a_inner The inner radius of the ring.
//! \param a_outer The outer radius of the ring.
//! 
//! \return The fraction of the pixel falling in the ring, from 0 to 1.
//------------------------------------------------------------------------------

double ae_brightness_in_ring(double x, double y, double res, 
                             double sin_ring_inc, double a_inner,
                             double a_outer) {
  double r2_inner, r2_outer, a2, b2, out_frac, in_frac;

  // Get the fractions falling within the disc defined by the ring's outer
  // and inner radii.
  out_frac = ae_brightness_in_disc(x, y, res, a_outer, a_outer * sin_ring_inc);
  in_frac = ae_brightness_in_disc(x, y, res, a_inner, a_inner * sin_ring_inc);

  // The total fraction for the ring is the outer minus the inner.
  return out_frac - in_frac;

  a2 = a_inner * a_inner;
  b2 = a2 * sin_ring_inc * sin_ring_inc;
  r2_inner = x * x / a2 + y * y / b2;

  a2 = a_outer * a_outer;
  b2 = a2 * sin_ring_inc * sin_ring_inc;
  r2_outer = x * x / a2 + y * y / b2;

  if ((r2_inner > 1.0) && (r2_outer < 1.0))
    return 1;
  else
    return 0;
}


//------------------------------------------------------------------------------
//! Make a map of a planet with rings.
//! This function creates a bit map of a planet with rings, given the
//! temperatures of the disc as well as the rings.
//!
//! \param res The resolution for the map, in units of the planet's equatorial
//!            radius.
//! \param b The body to be mapped.
//! \param r The parameters of the rings of the body.
//! \param inclination The inclination of the rings relative to the viewer, in
//!                    degrees.
//! \param sm The model of Saturn to be used.
//! \param map For returning the map. The member \p data of this structure
//!            will be allocated by this subroutine.
//! \param comp Unless NULL, returns a \p ae_saturn_components_t
//!             struct, whhich records the number of pixels in the map
//!             containing just the disc, the disc seen through rings and the
//!             background seen through rings. This information can be used to
//!             rapidly compute the effective temperature for various disc and
//!             ring temperatures without needing to redraw the map, using
//!             ae_saturn_t_eff().
//------------------------------------------------------------------------------

void ae_map_saturn(double res, const struct ae_physical_t *b,
                   const struct ae_ring_geometry_t *r, double inclination,
                   const struct ae_saturn_temp_model_t *sm,
                   struct ae_planet_map_t *map,
                   struct ae_saturn_components_t *comp) {
  int i, j, k, n_pix, x_half;
  double x, y, *exp_tau, pathlength, sin_ring_inc, cos_ring_inc, r_b, span;
  double disc_frac, this_ring_frac, ring_frac, disc_q, disc_u;
  struct ae_stokes_t **m, t_on, t_off;

  // Calculate the opacities of the rings.
  sin_ring_inc = fabs(sin(inclination * AE_DTR));
  cos_ring_inc = fabs(cos(inclination * AE_DTR));
  exp_tau = (double *)malloc(r->n * sizeof(double));
  for (i = 0; i < r->n; i++) {
    if ((pathlength = 1.0 / sin_ring_inc) > 100.0)
      pathlength = 100.0;
    exp_tau[i] = exp(-r->tau[i] * pathlength);
  }

  // Calculate the semi-minor axis for our ring inclination, where the
  // semi-major axis is defined as unity.
  x = b->r_pole / b->r_eq;
  r_b = ae_disc_semiminor(1.0, x, inclination);

  // Figure out how big the map should be, using the size of the largest ring.
  span = r->a_outer[0] * 1.05;
  n_pix = (int)(span / res) * 2 + 1;
  x_half = (int)(span / res);
  span = (double)x_half * res;

  // Figure out the q and u factors for the disc.
  disc_q = sm->pol_disc * cos(sm->pol_disc_angle * AE_DTR);
  disc_u = sm->pol_disc * sin(sm->pol_disc_angle * AE_DTR);

  // Allocate.
  m = (struct ae_stokes_t **)malloc(n_pix * sizeof(struct ae_stokes_t *));
  for (i = 0; i < n_pix; i++)
    m[i] = (struct ae_stokes_t *)malloc(n_pix * sizeof(struct ae_stokes_t));
  if (comp) {
    ae_saturn_components_alloc(r->n, comp);
    for (i = 0; i < r->n; i++)
      comp->exp_tau[i] = exp_tau[i];
    comp->cos_ring_inc = cos_ring_inc;
  }

  // Now make the map.
  for (i = 0, x = -span; i < n_pix; i++, x += res) {
    for (j = 0, y = -span; j < n_pix; j++, y += res) {
      // Get fraction of pixel that falls in the main disc.
      disc_frac = ae_brightness_in_disc(x, y, res, 1.0, r_b);
      ring_frac = 0.0;

      // First:  process the portion on-disc.
      ae_set_stokes(&t_on, 0, 0, 0, 0);
 
      // We are looking towards the disc.  Is there a ring over it?
      if ((inclination < 0 && y > 0) || (inclination > 0 && y < 0)) {
        // Work from inner rings out.
        for (k = 0; k < r->n; k++) {
          this_ring_frac = ae_brightness_in_ring(x, y, res, sin_ring_inc, 
                                                 r->a_inner[k], r->a_outer[k]);
          t_on.i += (sm->t_disc * exp_tau[k] + sm->t_ring[k]) * this_ring_frac;
          t_on.q += (sm->t_disc * disc_q * exp_tau[k] +
                     sm->t_ring[k] * cos_ring_inc * sm->pol_ring[k]) *
                    this_ring_frac;
          t_on.u += sm->t_disc * disc_u * exp_tau[k] * this_ring_frac;
          ring_frac += this_ring_frac;
          if (comp)
            comp->n_ring_disc[k] += this_ring_frac * disc_frac;
        }
      }
      
      t_on.i += sm->t_disc * (1.0 - ring_frac);
      t_on.q += sm->t_disc * disc_q * (1.0 - ring_frac);
      t_on.u += sm->t_disc * disc_u * (1.0 - ring_frac);
      t_on.i *= disc_frac;
      t_on.q *= disc_frac;
      t_on.u *= disc_frac;
      if (comp)
        comp->n_disc += (1.0 - ring_frac) * disc_frac;

      // Now:  process the portion off-disc.
      ae_set_stokes(&t_off, 0, 0, 0, 0);
      ring_frac = 0;
        
      // We are looking off-disc.  Are there rings?
      for (k = 0; k < r->n; k++) {
        this_ring_frac = ae_brightness_in_ring(x, y, res, sin_ring_inc, 
                                               r->a_inner[k], r->a_outer[k]);
        t_off.i += (sm->t_background * exp_tau[k] + sm->t_ring[k]) * 
                   this_ring_frac;
        t_off.q += sm->t_ring[k] * cos_ring_inc * this_ring_frac * 
                   sm->pol_ring[k];
        ring_frac += this_ring_frac;
        if (comp)
          comp->n_ring_back[k] += this_ring_frac * (1.0 - disc_frac);
      }
      
      t_off.i += sm->t_background * (1.0 - ring_frac);
      t_off.i *= (1.0 - disc_frac);
      t_off.q *= (1.0 - disc_frac);

      m[i][j].i = t_on.i + t_off.i;
      m[i][j].q = t_on.q + t_off.q;
      m[i][j].u = t_on.u + t_off.u;
      m[i][j].v = 0.0;
    }
  }

  free(exp_tau);
  map->data = m;
  map->res = res;
  map->n_pix = n_pix;

  return;
}


//------------------------------------------------------------------------------
//! Get the effective temperature of Saturn.
//! This is an alternative to ae_map_t_eff() for Saturn's temperature. It makes
//! use of previously-calculated tabulations of how much of the map falls on
//! disc, on rings, etc., which can be returned in the parameter \p comp of
//! ae_map_saturn().
//!
//! This function is useful if effective temperatures are being found for many
//! different ring/disc temperatures but the same viewing geometery, as new maps
//! do not need to be calculated for each new combination. It is therefore much
//! faster.
//!
//! \param r The parameters of the rings of the body.
//! \param sm The model of Saturn to be used.
//! \param comp A \p ae_saturn_components_t struct calculated by 
//!             ae_map_saturn().
//! \param t_eff For returning the Stokes parameters of the effective
//!              temperature.
//------------------------------------------------------------------------------

void ae_saturn_t_eff(const struct ae_ring_geometry_t *r,
                     const struct ae_saturn_temp_model_t *sm,
                     const struct ae_saturn_components_t *comp,
                     struct ae_stokes_t *t_eff) {
  int i;
  double n_disc, n_non_disc, disc_q, disc_u;

  // Figure out the q and u factors for the disc.
  disc_q = sm->pol_disc * cos(sm->pol_disc_angle * AE_DTR);
  disc_u = sm->pol_disc * sin(sm->pol_disc_angle * AE_DTR);

  // First, the unobstructed disc.
  t_eff->i = sm->t_disc * comp->n_disc;
  t_eff->q = sm->t_disc * disc_q * comp->n_disc;
  t_eff->u = sm->t_disc * disc_u * comp->n_disc;
  t_eff->v = 0;
  n_disc = comp->n_disc;
  n_non_disc = 0;

  // Next, rings over the disc.
  for (i = 0; i < r->n; i++) {
    t_eff->i += (sm->t_disc * comp->exp_tau[i] + sm->t_ring[i]) *
                comp->n_ring_disc[i];
    t_eff->q += (sm->t_disc * disc_q * comp->exp_tau[i] +
                 sm->t_ring[i] * comp->cos_ring_inc * sm->pol_ring[i]) *
                comp->n_ring_disc[i];
    t_eff->u += sm->t_disc * disc_u * comp->exp_tau[i] * sm->pol_ring[i];
    n_disc += comp->n_ring_disc[i];
  }

  // Finally, rings off-disc.
  for (i = 0; i < r->n; i++) {
    t_eff->i += (sm->t_background * comp->exp_tau[i] + sm->t_ring[i]) * 
                comp->n_ring_back[i];
    t_eff->q += sm->t_ring[i] * comp->cos_ring_inc * sm->pol_ring[i] *
                comp->n_ring_back[i];
    n_non_disc += comp->n_ring_back[i];
  }

  // Subtract the background temperature from the off-ring discs.
  t_eff->i -= sm->t_background * n_non_disc;

  // Divide by the solid angle of the disc to get the effective temperature.
  t_eff->i /= n_disc;
  t_eff->q /= n_disc;
  t_eff->u /= n_disc;
  t_eff->v /= n_disc;

  return;
}


//------------------------------------------------------------------------------
//! Given a planet map, determine the effective temperature.
//! The effective temperature is the total temperature of the map divided by the
//! area of the disc, after the total background temperature has been 
//! subtracted.
//!
//! \param map The planet map.
//! \param b The physical parameters of the planet.
//! \param subobs_lat The sub-observer latitude of the planet.
//! \param background The value of the background, which will be removed from
//!                   the effective temperature.
//! \param t_eff For returning the Stokes parameters of the effective
//!              temperature.
//------------------------------------------------------------------------------

void ae_map_t_eff(const struct ae_planet_map_t *map,
                  const struct ae_physical_t *b, double subobs_lat,
                  double background, struct ae_stokes_t *t_eff) {
  int i, j;
  double disc_sa, map_sa;
  struct ae_stokes_t avg;

  // Sum all the pixels.
  ae_set_stokes(&avg, 0, 0, 0, 0);
  for (i = 0; i < map->n_pix; i++) {
    for (j = 0; j < map->n_pix; j++) {
      avg.i += map->data[i][j].i;
      avg.q += map->data[i][j].q;
      avg.u += map->data[i][j].u;
      avg.v += map->data[i][j].v;
    }
  }


  // Figure out the area of the disc, in units of the planet's equatorial
  // radius.
  disc_sa = M_PI * ae_disc_semiminor(1.0, b->r_pole / b->r_eq, subobs_lat);

  // The area of the map, in units of the planet's equatorial radius.
  map_sa = map->res * map->res * map->n_pix * map->n_pix;

  // Convert the sum to units where the planet's equatorial radius is unity,
  avg.i *= map->res * map->res;
  avg.q *= map->res * map->res;
  avg.u *= map->res * map->res;
  avg.v *= map->res * map->res;

  // Remove the background.  (Assume that it is not polarised, so only subtract
  // from the intensity.)
  avg.i -= background * (map_sa - disc_sa);

  // Divide by the area of the disc.
  avg.i /= disc_sa;
  avg.q /= disc_sa;
  avg.u /= disc_sa;
  avg.v /= disc_sa;

  ae_copy_stokes(&avg, t_eff);
}

//------------------------------------------------------------------------------
//! Write a planet map to disc in a format readable by gnuplot.
//! The maps can be plotted in gnuplot with the command:
//! \pre gnuplot> plot "filename.map" binary with image
//!
//! \param map The map to be plotted.
//! \param stokes The Stokes parameter to be plotted:  one of 'i', 'q', 'u' or
//!               'v'.
//! \param path The file to write to.
//!
//! \return 0 on success; -1 on error.
//------------------------------------------------------------------------------

int ae_map_to_gnuplot(const struct ae_planet_map_t *map, char stokes,
                      const char *path) {
  int i, j, fd, n_pix_half;
  ssize_t dummy;
  float *writer;

  if ((fd = creat(path, 0666)) < 0)
    return -1;

  writer = (float *)malloc((map->n_pix + 1) * sizeof(float));
  n_pix_half = (int)(map->n_pix / 2);

  writer[0] = (float)map->n_pix;
  for (i = 1; i <= map->n_pix; i++)
    writer[i] = (float)((i - 1 - n_pix_half) * map->res);
  
  dummy = write(fd, writer, (map->n_pix + 1) * sizeof(float));
  for (i = 0; i < map->n_pix; i++) {
    writer[0] = (float)((i - n_pix_half) * map->res);
    for (j = 0; j < map->n_pix; j++) {
      switch (stokes) {
        case 'i': case 'I':
          writer[j + 1] = (float)map->data[j][i].i;
          break;
        case 'q': case 'Q':
          writer[j + 1] = (float)map->data[j][i].q;
          break;
        case 'u': case 'U':
          writer[j + 1] = (float)map->data[j][i].u;
          break;
        case 'v': case 'V': default:
          writer[j + 1] = (float)map->data[j][i].v;
          break;
      }
    }
    dummy = write(fd, writer, (map->n_pix + 1) * sizeof(float));
  }
  
  close(fd);

  free(writer); 

  return 0;
}
