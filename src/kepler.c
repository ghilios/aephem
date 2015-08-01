//! \file kepler.c
//! Solve Keplerian orbits.
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
#include <string.h>

#include "aephem.h"

//------------------------------------------------------------------------------
//! Adjust position from earth-moon barycentre to earth.
//!
//! \param jd_tt The Julian date in TT.
//! \param rect Equatorial rectangular coordinates of EMB.
//!
//! \return Earth's distance to the sun, in AU.
//------------------------------------------------------------------------------

double ae_kepler_embofs(double jd_tt, double rect[]) {
  int i;
  double pm[3], polm[3], a, b;

  // Compute the vector Moon - Earth.
  ae_gmoon(jd_tt, pm, polm);

  // Precess the lunar position to ecliptic and equinox of AE_J2000.0.
  ae_precess(jd_tt, pm, AE_TO_J2000);

  // Adjust the coordinates of the Earth.
  a = 1.0 / (AE_E_M_RAT +  1.0);
  b = 0.0;
  for(i = 0; i < 3; i++) {
    rect[i] = rect[i] - a * pm[i];
    b = b + rect[i] * rect[i];
  }
  
  return sqrt(b);
}


//------------------------------------------------------------------------------
//! Find a Keplerian orbit.
//! This routine solves for a Keplerian orbit, given orbital parameters and the
//! time.
//!
//! The routine detects several cases of given orbital elements.  If a program 
//! for perturbations is pointed to, it is called to calculate all the
//! elements. If there is no program, then the mean longitude is calculated
//! from the mean anomaly and daily motion.  If the daily motion is not given, 
//! it is calculated by Kepler's law.  If the eccentricity is given to be 1.0, 
//! it means that meandistance is really the perihelion distance, as in a comet
//! specification, and the orbit is parabolic.
//! 
//! Reference:
//! - Taff, L.G., "Celestial Mechanics, A Computational Guide for the 
//!   Practitioner."  Wiley, 1985.
//!
//! \param jd_tt The Julian date in TT.
//! \param orb The orbital elements.
//! \param q For returning the heliocentric equatorial rectangular coordinates 
//!          of the object, in AU.
//------------------------------------------------------------------------------

void ae_kepler(double jd_tt, const struct ae_orbit_t *orb, double q[]) {
  double alat, E, M, W, temp, epoch, inclination, ascnode, argperih;
  double meandistance, dailymotion, eccent, meananomaly, r, coso, sino, cosa;
  double eps, coseps, sineps, polar[3], equinox, orb_l;

  // Call program to compute position, if one is supplied.
  if (orb->ptable) {
    if (!strcasecmp(orb->name, "earth"))
      ae_g3plan(jd_tt, orb->ptable, polar, 3);
    else
      ae_gplan(jd_tt, orb->ptable, polar);
    E = polar[0]; // Longitude.
    //orb.L = E;
    W = polar[1]; // Latitude.
    r = polar[2]; // Radius.
    //orb.r = r;
    //orb.epoch = jd_tt;
    equinox = AE_J2000;

    goto kepdon;
  }

  // Decant the parameters from the data structure.
  epoch = orb->epoch;
  equinox = orb->equinox;
  inclination = orb->i;
  ascnode = orb->W * AE_DTR;
  argperih = orb->w;
  meandistance = orb->a; // Semi-major axis.
  dailymotion = orb->dm;
  eccent = orb->ecc;
  meananomaly = orb->M;
  // Check for parabolic orbit.
  if (eccent == 1.0) {
    // Meandistance = perihelion distance, q.
    // Epoch = perihelion passage date.
    temp = meandistance * sqrt(meandistance);
    W = (jd_tt - epoch ) * 0.0364911624 / temp;
    // The constant above is 3 k / sqrt(2), where k = Gaussian gravitational 
    // constant = 0.01720209895.
    E = 0.0;
    M = 1.0;
    while (fabs(M) > 1.0e-11) {
      temp = E * E;
      temp = (2.0 * E * temp + W) / ( 3.0 * (1.0 + temp));
      M = temp - E;
      if (temp != 0.0)
        M /= temp;
      E = temp;
    }
    r = meandistance * (1.0 + E * E);
    M = atan(E);
    M = 2.0 * M;
    alat = M + AE_DTR * argperih;

    goto parabcon;
  }

  if (eccent > 1.0)	{
    // The equation of the hyperbola in polar coordinates r, theta is 
    // r = a(e^2 - 1)/(1 + e cos(theta)) so the perihelion distance q = a(e-1),
    // the "mean distance"  a = q/(e-1).
    meandistance = meandistance / (eccent - 1.0);
    temp = meandistance * sqrt(meandistance);
    W = (jd_tt - epoch ) * 0.01720209895 / temp;

    // Solve M = -E + e sinh E.
    E = W / (eccent - 1.0);
    M = 1.0;
    while (fabs(M) > 1.0e-11) {
      M = -E + eccent * sinh(E) - W;
      E += M / (1.0 - eccent * cosh(E));
    }

    r = meandistance * (-1.0 + eccent * cosh(E));
    temp = (eccent + 1.0)/(eccent - 1.0);
    M = sqrt(temp) * tanh( 0.5*E );
    M = 2.0 * atan(M);	
    alat = M + AE_DTR * argperih;

    goto parabcon;
  }

  // Calculate the daily motion, if it is not given.
  if (dailymotion == 0.0) {
    // The constant is 180 k / pi, k = Gaussian gravitational constant.
    // Assumes object in heliocentric orbit is massless.
    dailymotion = 0.9856076686/(orb->a * sqrt(orb->a));
  }

  dailymotion *= jd_tt - epoch;

  // M is proportional to the area swept out by the radius vector of a circular
  // orbit during the time between perihelion passage and Julian date jd_tt.  
  // It is the mean anomaly at time J.
  M = AE_DTR * (meananomaly + dailymotion);
  M = ae_mod_2pi(M);

  // If mean longitude was calculated, adjust it also for motion since epoch of
  // elements.
  if (orb->L) {
    orb_l = orb->L + dailymotion;
    orb_l = ae_mod_360(orb_l);
  }
  else
    orb_l = 0;

  // By Kepler's second law, M must be equal to the area swept out in the same 
  // time by an elliptical orbit of same total area.  Integrate the ellipse
  // expressed in polar coordinates r = a(1-e^2)/(1 + e cosW) with respect to 
  // the angle W to get an expression for the area swept out by the radius 
  // vector.  The area is given by the mean anomaly; the angle is solved 
  // numerically.
  //
  // The answer is obtained in two steps.  We first solve Kepler's equation
  // M = E - eccent * sin(E) for the eccentric anomaly E.  Then there is a
  // closed form solution for W in terms of E.

  E = M; // Initial guess is same as circular orbit.
  temp = 1.0;
  do {
    // The approximate area swept out in the ellipse . . .
    temp = E - eccent * sin(E)
           - M; // . . . minus the area swept out in the circle.
    // . . . should be zero.  Use the derivative of the error to converge to 
    // solution by Newton's method.
    E -= temp / (1.0 - eccent * cos(E));
  } while (fabs(temp) > 1.0e-11);

  // The exact formula for the area in the ellipse is
  //   2.0 * atan(c2 * tan(0.5 * W)) - c1 * eccent * sin(W) / (1 + e * cos(W))
  // where
  //   c1 = sqrt(1.0 - eccent * eccent)
  //   c2 = sqrt((1.0 - eccent) / (1.0 + eccent)).
  // Substituting the following value of W yields the exact solution.
  temp = sqrt((1.0 + eccent) / (1.0 - eccent));
  W = 2.0 * atan(temp * tan(0.5 * E));

  // The true anomaly.
  W = ae_mod_2pi(W);

  meananomaly *= AE_DTR;
  // Orbital longitude measured from node (argument of latitude).
  if (orb_l)
    alat = (orb_l) * AE_DTR + W - meananomaly - ascnode;
  else
    alat = W + AE_DTR * argperih; // Mean longitude not given.

  // From the equation of the ellipse, get the radius from central focus to the
  // object.
  r = meandistance * (1.0 - eccent * eccent) / (1.0 + eccent * cos(W));

  parabcon:

  // The heliocentric ecliptic longitude of the object is given by
  //   tan(longitude - ascnode) = cos(inclination) * tan(alat).
  coso = cos(alat);
  sino = sin(alat);
  inclination *= AE_DTR;
  W = sino * cos(inclination);
  E = ae_zatan2(coso, W) + ascnode;

  // The ecliptic latitude of the object.
  W = sino * sin(inclination);
  W = asin(W);

  kepdon:

  // Convert to rectangular coordinates, using the perturbed latitude.
  q[2] = r * sin(W);
  cosa = cos(W);
  q[1] = r * cosa * sin(E);
  q[0] = r * cosa * cos(E);

  // Convert from heliocentric ecliptic rectangular to heliocentric equatorial
  // rectangular coordinates by rotating eps radians about the x axis.
  eps = ae_epsilon(equinox);
  coseps = cos(AE_STR * eps);
  sineps = sin(AE_STR * eps);
  W = coseps * q[1] - sineps * q[2];
  M = sineps * q[1] + coseps * q[2];
  q[1] = W;
  q[2] = M;

  // Precess the position to ecliptic and equinox of AE_J2000.0 if not already 
  // there.
  ae_precess(equinox, q, AE_TO_J2000);

  // If earth, adjust from earth-moon barycenter to earth by AA page E2.
  if (!strcasecmp(orb->name, "earth"))
    r = ae_kepler_embofs(jd_tt, q);

  // Rotate back into the ecliptic.
  eps = ae_epsilon(AE_J2000);
  coseps = cos(AE_STR * eps);
  sineps = sin(AE_STR * eps);
  W = coseps * q[1] + sineps * q[2];
  M = -sineps * q[1] + coseps * q[2];
  
  return;
}
