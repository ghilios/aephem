//! \file phys_data.c
//! This file contains physical data for solar system data.
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

#include "aephem.h"

//! An empty series for ephemeride terms with no extra terms.
struct ae_physical_term_t ae_phys_no_term[] = {
  {AE_PHYSICAL_END, 0, 0, 0, 0}
};

//! Physical parameters for the Sun.
struct ae_physical_t ae_phys_sun = {
  696000,
  696000,
  696000,
  0,
  0,
  0,
  286.13,          0, ae_phys_no_term,
   63.87,          0, ae_phys_no_term,
   84.10, 14.1844000, 0, ae_phys_no_term
};

//! Physical parameters for Mercury:  sine terms for w (M1 through M5).
struct ae_physical_term_t ae_phys_mercury_w[] = {
  {AE_PHYSICAL_D,  0.00993822, 174.791086,   4.0923350, 0},
  {AE_PHYSICAL_D, -0.00104581, 349.582171,   8.1846700, 0},
  {AE_PHYSICAL_D, -0.00010280, 164.373257,   12.277005, 0},
  {AE_PHYSICAL_D, -0.00002364, 339.164343,   16.369340, 0},
  {AE_PHYSICAL_D, -0.00000532, 153.955429,  20.4616750, 0},
  {AE_PHYSICAL_END, 0, 0, 0, 0}
};

//! Physical parameters for Mercury.
struct ae_physical_t ae_phys_mercury = {
  2439.7,
  2439.7,
  2439.7,
  1.0,
  4.6,
  2.5,
  281.0100, -0.0328000, ae_phys_no_term,
   61.4500, -0.0049000, ae_phys_no_term,
  329.5469,  6.1385025, 0, ae_phys_mercury_w
};

//! Physical parameters for Venus.
struct ae_physical_t ae_phys_venus = {
  6051.8,
  6051.8,
  6051.8,
  1,
  11,
  2,
  272.76,  0.0000000, ae_phys_no_term,
   67.16,  0.0000000, ae_phys_no_term,
  160.20, -1.4813688, 0, ae_phys_no_term
};

//! Physical parameters for Earth.
struct ae_physical_t ae_phys_earth = {
  6371.00,
  6378.14,
  6356.75,
  3.57,
  8.85,
  11.52,
    0.000,  -0.6410000, ae_phys_no_term,
   90.000,  -0.5570000, ae_phys_no_term,
  190.147, 360.9856235, 0, ae_phys_no_term
};

//! Physical parameters for Mars.
struct ae_physical_t ae_phys_mars = {
  3389.50,
  3396.19,
  3376.20,
  3.0,
  22.64,
  7.55,
  317.68143,  -0.10610000, ae_phys_no_term,
   52.88650,  -0.06090000, ae_phys_no_term,
  176.63000, 350.89198226, 0, ae_phys_no_term
};

//! Physical parameters for Jupiter:  sine terms for ra (Ja through Je).
struct ae_physical_term_t ae_phys_jupiter_sincos[] = {
  {AE_PHYSICAL_T, 0.000117,  99.360714, 4850.4046, 0},
  {AE_PHYSICAL_T, 0.000938, 175.895369, 1191.9605, 0},
  {AE_PHYSICAL_T, 0.001432, 300.323162,  262.5475, 0},
  {AE_PHYSICAL_T, 0.000030, 114.012305, 6070.2476, 0},
  {AE_PHYSICAL_T, 0.002150,  49.511251,   64.3000, 0},
  {AE_PHYSICAL_END, 0, 0, 0, 0}
};

//! Physical parameters for Jupiter.
struct ae_physical_t ae_phys_jupiter = {
  69911,
  71492,
  66854,
  62.1,
  31,
  102,
  268.056595,  -0.0064990, ae_phys_jupiter_sincos,
   64.495303,   0.0024130, ae_phys_jupiter_sincos,
  284.950000, 870.5360000, 0, ae_phys_no_term
};

//! Physical parameters for Saturn.
struct ae_physical_t ae_phys_saturn = {
  58232,
  60268,
  54364,
  102.9,
  8,
  205,
  40.589,  -0.0360000, ae_phys_no_term,
  83.537,  -0.0040000, ae_phys_no_term,
  38.900, 810.7939024, 0, ae_phys_no_term
};

// DEVEL start
//! Physical parameters for Saturn's rings' optical depths.
const double ae_saturn_rings_tau[] = {
  0.7, 0.1, 2.0, 1.0, 0.1, 0.15, 0.08
};

//! Physical parameters for Saturn's rings' inner radii.
const double ae_saturn_rings_a_inner[] = {
  2.025, 1.95, 1.64, 1.525, 1.43, 1.29, 1.24
};

//! Physical parameters for Saturn's rings' outer radii.
const double ae_saturn_rings_a_outer[] = {
  2.27, 2.025, 1.95, 1.64, 1.525, 1.43, 1.29
};

//! Names of Saturn's rings.
const char *ae_saturn_rings_name[] = {
  "A", "Cassini Division", "Outer B", "Inner B", "Outer C", "Middle C",
  "Inner C"
};

//! Physical parameters for Saturn's rings.
struct ae_ring_geometry_t ae_saturn_rings = {
  7,
  ae_saturn_rings_name,
  ae_saturn_rings_tau,
  ae_saturn_rings_a_inner,
  ae_saturn_rings_a_outer,
};
// DEVEL stop.

//! Physical parameters for Uranus.
struct ae_physical_t ae_phys_uranus = {
  25362,
  25559,
  24973,
  16.8,
  28,
  0,
  257.311,    0.0000000, ae_phys_no_term,
  -15.175,    0.0000000, ae_phys_no_term,
  203.810, -501.1600928, 0, ae_phys_no_term
};

//! Physical parameters for Neptune:  sine term for ra.
struct ae_physical_term_t ae_phys_neptune_ra_sin[] = {
  {AE_PHYSICAL_T,  0.70, 357.85, 52.316, 0},
  {AE_PHYSICAL_END, 0, 0, 0, 0}
};

//! Physical parameters for Neptune:  cosine term for ra.
struct ae_physical_term_t ae_phys_neptune_dec_cos[] = {
  {AE_PHYSICAL_T, -0.51, 357.85, 52.316, 0},
  {AE_PHYSICAL_END, 0, 0, 0, 0}
};

//! Physical parameters for Neptune:  sine term for w.
struct ae_physical_term_t ae_phys_neptune_w_sin[] = {
  {AE_PHYSICAL_T, -0.48, 357.85, 52.316, 0},
  {AE_PHYSICAL_END, 0, 0, 0, 0}
};

//! Physical parameters for Neptune.
struct ae_physical_t ae_phys_neptune = {
  24622,
  24764,
  24341,
  8,
  14,
  0,
  299.36,   0.0000000, ae_phys_neptune_ra_sin,
   43.46,   0.0000000, ae_phys_neptune_dec_cos,
  253.18, 536.3128492, 0, ae_phys_neptune_w_sin
};

//! Physical parameters for Pluto.
struct ae_physical_t ae_phys_pluto = {
  1195,
  1195,
  1195,
  0,
  0,
  0,
  132.993,  0.0000000, ae_phys_no_term,
   -6.163,  0.0000000, ae_phys_no_term,
  237.305, 56.3625225, 0, ae_phys_no_term
};
