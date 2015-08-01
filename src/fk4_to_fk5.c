//! \file fk4_to_fk5.c
//! Convert FK4 B1950.0 catalogue coordinates to FK5 J2000.0 coordinates.
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

//! Factors to eliminate E terms of aberration in FK4 to FK5 conversion.
double ae_fk4_to_fk5_abb[3] = {-1.62557e-6, -3.1919e-7, - 1.3843e-7};
//! Factors to eliminate E terms of aberration in FK4 to FK5 conversion.
double ae_fk4_to_fk5_abb_d[3] = {1.244e-3, -1.579e-3, -6.60e-4};

//! Tranformation matrix for FK4 to FK5 conversion.
//! Transformation matrix for unit direction vector, and motion vector in
//! seconds of arc per century.
double ae_fk4_to_fk5_mat[36] = {
  0.9999256782,  -0.0111820611,  -4.8579477e-3,
  2.42395018e-6, -2.710663e-8,   -1.177656e-8,
  0.0111820610,  0.9999374784,   -2.71765e-5,
  2.710663e-8,   2.42397878e-6,  -6.587e-11,
  4.8579479e-3,  -2.71474e-5,    0.9999881997,
  1.177656e-8,   -6.582e-11,     2.42410173e-6,
  -5.51e-4,      -0.238565,      0.435739,
  0.99994704,    -0.01118251,    -4.85767e-3,
  0.238514,      -2.667e-3,      -8.541e-3,
  0.01118251,    0.99995883,     -2.718e-5,
  -0.435623,     0.012254,       2.117e-3,
  4.85767e-3,    -2.714e-5,      1.00000956
};

//------------------------------------------------------------------------------
//! Convert FK4 B1950.0 catalogue coordinates to FK5 J2000.0 coordinates.
//! See AA page B58.
//!
//! \param el The catalogue object to update to J2000.0.
//------------------------------------------------------------------------------

void ae_fk4_to_fk5(struct ae_star_t *el) {
  int i, j;
  double a, b, c, *u, *v, R[6], p[3], m[3], cosdec, sindec, cosra, sinra, vpi;

  // Convert to rectangular coordinates.
  cosdec = cos(el->dec * AE_DTR);
  sindec = sin(el->dec * AE_DTR);
  cosra = cos(el->ra * AE_DTR);
  sinra = sin(el->ra * AE_DTR);
  p[0] = cosra * cosdec;
  p[1] = sinra * cosdec;
  p[2] = sindec;
  vpi = 21.094952663 * el->v * el->px;
  m[0] = -el->mura * cosdec * sinra - el->mudec * sindec * cosra + vpi * p[0];
  m[1] = el->mura * cosdec * cosra - el->mudec * sindec * sinra + vpi * p[1];
  m[2] = el->mudec * cosdec + vpi * p[2];

  a = 0.0;
  b = 0.0;
  for (i = 0; i < 3; i++) {
    a += ae_fk4_to_fk5_abb[i] * p[i];
    b += ae_fk4_to_fk5_abb_d[i] * p[i];
  }

  // Remove E terms of aberration from FK4.
  for(i = 0; i < 3; i++) {
    R[i] = p[i] - ae_fk4_to_fk5_abb[i] + a * p[i];
    R[i + 3] = m[i] - ae_fk4_to_fk5_abb_d[i] + b * p[i];
  }

  // Perform matrix multiplication.
  v = &ae_fk4_to_fk5_mat[0];
  for (i = 0; i < 6; i++) {
    a = 0.0;
    u = &R[0];
    for (j = 0; j < 6; j++)
      a += *u++ * *v++;
    if (i < 3)
      p[i] = a;
    else
      m[i - 3] = a;
  }

  // Transform the answers into J2000 catalogue entries in radian measure.
  b = p[0] * p[0] + p[1] * p[1];
  a = b + p[2] * p[2];
  c = a;
  a = sqrt(a);

  el->ra = ae_zatan2(p[0], p[1]) * AE_RTD;
  el->dec = asin(p[2] / a) * AE_RTD;

  // Note motion converted back to radians per (Julian) century.
  el->mura = (p[0] * m[1] - p[1] * m[0]) / b;
  el->mudec = (m[2] * b - p[2] * (p[0] * m[0] + p[1] * m[1]) ) / 
              (c * sqrt(b));

  if (el->px > 0.0) {
    c = 0.0;
    for (i = 0; i < 3; i++)
      c += p[i] * m[i];

    // Divide by AE_RTS to deconvert m (and therefore c) from arc seconds back 
    // to radians.
    el->v = c / (21.094952663 * el->px * a);
  }
  el->px /= a;  // A is dimensionless.
  el->epoch = AE_J2000;

  return;
}
