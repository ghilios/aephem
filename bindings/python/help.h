// \file help.h
// Defines the docstrings for the python binding methods.
//==============================================================================
// AEPHEM - an astronomical ephemeris and reduction library.
// Copyright 2012 Adam Hincks, Canadian Institute for Theoretical Astrophysics.
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

#ifndef HELP_H
#define HELP_H

//==============================================================================
// GENERAL FORMATTING DEFINITIONS
//==============================================================================

#define AEPY_HRULE \
   "------------------------------------------------------------------------\n"
#define AEPYO_HRULE \
   "--------------------------------------------------------------------\n"
#define AEPY_DBL_HRULE \
   "========================================================================\n"
#define AEPYO_DBL_HRULE \
   "====================================================================\n"
#define AEPY_METHOD_NAME ""
#define AEPYO_METHOD_NAME ""
#define AEPY_SHORT_DESC AEPY_DBL_HRULE
#define AEPYO_SHORT_DESC AEPYO_DBL_HRULE
#define AEPY_LONG_DESC "\n"
#define AEPYO_LONG_DESC "\n"
#define AEPY_PARAM "\nParameters\n" AEPY_HRULE
#define AEPYO_PARAM "\nParameters\n" AEPYO_HRULE
#define AEPY_RETURN "\nReturns\n" AEPY_HRULE
#define AEPYO_RETURN "\nReturns\n" AEPYO_HRULE
#define AEPY_SEE "\nSee also\n" AEPY_HRULE
#define AEPYO_SEE "\nSee also\n" AEPYO_HRULE

//==============================================================================
// HELP FOR CALENDAR.C
//==============================================================================

#define AEPY_DS_ctime_to_jd \
  AEPY_METHOD_NAME \
   "ctime_to_jd(t)\n"\
  AEPY_SHORT_DESC \
   "Convert a C time to a Julian date.\n"\
  AEPY_LONG_DESC \
   "The C time is the time since the Epoch (00:00:00 UTC, January 1, 1970),\n"\
   "measured in seconds.\n"\
  AEPY_PARAM \
   "t : The C time at which to calculate the Julian date.\n"\
  AEPY_RETURN \
   "The Julian date.\n"

#define AEPY_DS_ctime_to_last \
  AEPY_METHOD_NAME \
   "ctime_to_last(t, delta_t, [tlon, nutl, eps])\n"\
  AEPY_SHORT_DESC \
   "Calculate the local apparent sidereal time from a C time.\n"\
  AEPY_LONG_DESC \
   "The C time is the time since the Epoch (00:00:00 UTC, January 1, 1970),\n"\
   "measured in seconds.\n"\
   "\n"\
   "If optional arguments delta_t, nutl, eps are not passed, then they are\n"\
   "calculated by this function. Note that this makes the call slower,\n"\
   "especially if the nutation terms need to be calculated. If calling this\n"\
   "function repeatedly for similar dates, consider precalculating these\n"\
   "terms once and including them as optional arguments to increase\n"\
   "efficiency.\n"\
  AEPY_PARAM \
   "      t : The C time at which to calculate the Julian date.\n"\
   "   tlon : The longitude of the observer, in degrees.\n"\
   "delta_t : TT - UT1, in seconds. Optional.\n"\
   "   nutl : The nutation in longitude in seconds of arc. Optional.\n"\
   "    eps : The obliquity of the ecliptic in seconds of arc. Optional.\n"\
  AEPY_RETURN \
   "The local apparent sidereal time, in degrees.\n"\
  AEPY_SEE \
   "delta_t(), nutation_lon_ob(), epsilon()\n"

#define AEPY_DS_cal_to_jd \
  AEPY_METHOD_NAME \
   "cal_to_jd(year, month, day)\n"\
  AEPY_SHORT_DESC \
   "Calculate Julian date from a Gregorian calendar date.\n"\
  AEPY_LONG_DESC \
   "The Julian date is double precision floating point with the origin used\n"\
   "by astronomers.\n"\
   "\n"\
   "There is no year 0.  Enter B.C. years as negative; i.e., 2 B.C. = -2.\n"\
   "\n"\
   "The approximate range of dates handled is 4713 B.C. to 54,078 A.D. This\n"\
   "should be adequate for most applications.\n"\
   "\n"\
   "B.C. dates are calculated by extending the Gregorian sequence of leap\n"\
   "years and century years into the past.  This seems the only sensible\n"\
   "definition, but it may not be the official one.\n"\
   "\n"\
   "Note that the astronomical Julian day starts at noon on the previous\n"\
   "calendar day.  Thus at midnight in the morning of the present calendar\n"\
   "day the Julian date ends in 0.5; it rolls over to tomorrow at noon\n"\
   "today.\n"\
   "\n"\
   "The month finding algorithm is attributed to Meeus.\n"\
  AEPY_PARAM \
   " year : The Gregorian year (integer).\n"\
   "month : The month of the year (integer), with January = 1, February = 2,\n"\
   "        etc.\n"\
   "  day : The fractional day of the month, starting at 1.\n"\
  AEPY_RETURN \
   "The Julian day number.\n"\
  AEPY_SEE \
   "jd_to_cal(), ctime_to_jd()\n"

#define AEPY_DS_jd_to_cal \
  AEPY_METHOD_NAME \
   "jd_to_cal(jd)\n"\
  AEPY_SHORT_DESC \
   "Calculate month, day, and year from a Julian date.\n"\
  AEPY_PARAM \
   "jd : A Julian day number.\n"\
  AEPY_RETURN \
   "A (year, month, day) tuplet.\n"\
  AEPY_SEE \
   "cal_to_jd()\n"

#define AEPY_DS_gmst \
  AEPY_METHOD_NAME \
   "gmst(jd_ut1, jd_tt)\n"\
  AEPY_SHORT_DESC \
   "Get the Greenwich mean sidereal time.\n"\
  AEPY_LONG_DESC \
   "Get the mean sidereal time at Greenwich (i.e., longitude 0).  The mean\n"\
   "sidereal time does not include nutation corrections.  Coefficients are\n"\
   "from the IAU and can be found in:\n"\
   "\n"\
   "     George H. Kaplan, \"The IAU Resolutions on Astronomical Reference\n"\
   "     Systems, Time Scales, and Earth Rotation Models,\" United States\n"\
   "     Naval Observatory Circular No. 179, 2005.\n"\
  AEPY_PARAM \
   "jd_ut1 : The Julian date, in UT1.\n"\
   " jd_tt : The Julian date, in TT.\n"\
  AEPY_RETURN \
   "The Greenwich mean sidereal time, in degrees.\n"\
  AEPY_SEE \
   "dut1(), delta_t(), gast(), lmst(), last(), ctime_to_last()\n"

#define AEPY_DS_gast \
  AEPY_METHOD_NAME \
   "gast(jd_ut1, jd_tt, nutl, eps)\n"\
  AEPY_SHORT_DESC \
   "Get Greenwich apparent sidereal time.\n"\
  AEPY_LONG_DESC \
   "The Greenwich apparent sidereal time is the Greenwich mean sidereal time\n"\
   "corrected for nutation.  See gmst() for reference.\n"\
  AEPY_PARAM \
   "jd_ut1 : The Julian date, in UT1.\n"\
   " jd_tt : The Julian date, in TT.\n"\
   "  nutl : The nutation in longitude in seconds of arc.\n"\
   "   eps : The obliquity of the ecliptic in seconds of arc.\n"\
  AEPY_RETURN \
   "Greenwich apparent sidereal time, in degrees.\n"\
  AEPY_SEE \
   "dut1(), delta_t(), nutation_lon_ob(), epsilon(), gmst(), lmst(), last(),\n"\
   "ctime_to_last()\n"

#define AEPY_DS_lmst \
  AEPY_METHOD_NAME \
   "lmst(jd_ut1, jd_tt, tlon)\n"\
  AEPY_SHORT_DESC \
   "Get the local mean sidereal time.\n"\
  AEPY_LONG_DESC \
   "The local mean sidereal time is simply the Greenwich mean sidereal time\n"\
   "plus the local longitude.  See gmst().\n"\
  AEPY_PARAM \
   "jd_ut1 : The Julian date in UT1.\n"\
   " jd_tt : The Julian date in TT.\n"\
   "  tlon : The longitude of the obsesver in degrees.\n"\
  AEPY_RETURN \
   "Local mean sidereal time, in degrees.\n"\
  AEPY_SEE \
   "dut1(), delta_t(), gmst() gast(), last(), ctime_to_last()\n"

#define AEPY_DS_last \
  AEPY_METHOD_NAME \
   "last(jd_ut1, tlon [, jd_tt, nutl, eps])\n"\
  AEPY_SHORT_DESC \
   "Get the local apparent sidereal time.\n"\
  AEPY_LONG_DESC \
   "The local apparent sidereal time is simply the Greenwich apparent\n"\
   "sidereal time plus the local longitude.  See gast().\n"\
   "\n"\
   "If optional arguments jd_tt, nutl, eps are not passed, then they are\n"\
   "calculated by this function. Note that this makes the call slower,\n"\
   "especially if the nutation terms need to be calculated. If calling this\n"\
   "function repeatedly for similar dates, consider precalculating these\n"\
   "terms once and including them as optional arguments to increase\n"\
   "efficiency.\n"\
  AEPY_PARAM \
   "jd_ut1 : The Julian date in UT1.\n"\
   "  tlon : The longitude of the observer in degrees.\n"\
   " jd_tt : The Julian date in TT. Optional.\n"\
   "  nutl : The nutation in longitude in seconds of arc. Optional.\n"\
   "   eps : The obliquity of the ecliptic in seconds of arc. Optional.\n"\
  AEPY_RETURN \
   "Local apparent sidereal time in degrees.\n"\
  AEPY_SEE \
   "dut1(), delta_t(), nutation_lon_ob(), epsilon(), gmst(), lmst(), lmst(),\n"\
   "ctime_to_last()\n"

#define AEPY_DS_tdb \
  AEPY_METHOD_NAME \
   "tdb(jd)\n"\
  AEPY_SHORT_DESC \
   "Find Barycentric Dynamical Time from Terrestrial Dynamical Time (TDT).\n"\
  AEPY_LONG_DESC \
   "See AA page B5.\n"\
  AEPY_PARAM \
   "jd : A Julian date, in TDT.\n"\
  AEPY_RETURN \
   "The corresponding time in TDB.\n"

#define AEPY_DS_mjd \
  AEPY_METHOD_NAME \
   "mjd(jd)\n"\
  AEPY_SHORT_DESC \
   "Given a Julian Date, get the Modified Julian Date.\n"\
  AEPY_LONG_DESC \
   "The difference between the two is simply a constant, MJD_START.\n"\
  AEPY_PARAM \
   "jd : A Julian date.\n"\
  AEPY_RETURN \
   "The Modified Julian Date.\n"\
  AEPY_SEE \
   "MJD_START\n"

#define AEPY_DS_delta_t \
  AEPY_METHOD_NAME \
   "delta_t(jd_ut1)\n"\
  AEPY_SHORT_DESC \
   "Get the value of Terrestrial Time (TT) - Universal Time (UT1).\n"\
  AEPY_LONG_DESC \
   "This routine uses a table for historic values and near-future\n"\
   "predictions.\n"\
   "\n"\
   "The program adjusts for a value of secular tidal acceleration ndot. It\n"\
   "is -25.8 arcsec per century squared for JPL's DE403 ephemeris.  ELP2000\n"\
   "and DE200 use the value -23.8946.\n"\
   "\n"\
   "For dates earlier than the tabulated range, the program calculates\n"\
   "approximate formulae of Stephenson and Morrison or K. M. Borkowski.\n"\
   "These approximations have an estimated error of 15 minutes at 1500 B.C.\n"\
   "They are not adjusted for small improvements in the current estimate of\n"\
   "ndot because the formulas were derived from studies of ancient eclipses\n"\
   "and other historical information, whose interpretation depends only\n"\
   "partly on ndot.\n"\
   "\n"\
   "A quadratic extrapolation formula, that agrees in value and slope with\n"\
   "current data, predicts future values of deltaT.\n"\
  AEPY_PARAM \
   "jd_ut1 : The Julian date.\n"\
  AEPY_RETURN \
   "The TT - UT (delta T) in seconds.\n"\
  AEPY_SEE \
   "dut1()\n"

#define AEPY_DS_dut1 \
  AEPY_METHOD_NAME \
   "dut1(jd_utc)\n"\
  AEPY_SHORT_DESC \
   "Given a UTC, return DUT1.\n"\
  AEPY_LONG_DESC \
   "DUT1 is the offset between UTC and UT1, i.e., DUT1 = UT1 - UTC.\n"\
   "Functions that find sidereal time expect UT1, so this function will be\n"\
   "required to convert UTC if reasonable precision (better than roughly 10\n"\
   "seconds of arc).\n"\
   "\n"\
   "Values are taken from a look-up table.  Past values start on 19 May 1976\n"\
   "are are tabulated daily until the release date of this version of\n"\
   "aephem.  After this, predicted values are used through about one year.\n"\
   "\n"\
   "If DUT1 is requested for a date before the tabulation starts, 0 will be\n"\
   "returned.  After the last predicted value, a value will be estimated\n"\
   "using an interpolation formula.\n"\
   "\n"\
   "Past values, predicted values and the interpolation formula were taken\n"\
   "from the U.S. Naval Observatory IERS Bulletin A.  This is available\n"\
   "online at http://maia.usno.navy.mil/.\n"\
  AEPY_PARAM \
   "jd_utc : The Julian date in UTC.\n"\
  AEPY_RETURN \
   "DUT1, in seconds. If jd_utc is before the first tabulated value, 0 will\n"\
   "be returned.\n"

//==============================================================================
// BINDINGS FOR COORD.C
//==============================================================================

#define AEPY_DS_rect_to_polar \
  AEPY_METHOD_NAME \
   "rect_to_polar(rect [, radius])\n"\
  AEPY_SHORT_DESC \
   "Convert an equatorial rectangular unit vector to ra/dec/radius.\n"\
  AEPY_PARAM \
   "  rect : A list, of length three, containing the rectantular vector.\n"\
   "radius : If true, then return the radius in addition to the ra/dec.\n"\
  AEPY_RETURN \
   "A tuplet containing (ra, dec), in degrees; or, if radius = True,\n"\
   "(ra, dec, radius).\n"\
  AEPY_SEE \
   "polar_to_rect()\n"

#define AEPY_DS_polar_to_rect \
  AEPY_METHOD_NAME \
   "polar_to_rect(ra, dec, radius)\n"\
  AEPY_SHORT_DESC \
   "Convert ra/dec/radius to a rectangular vector.\n"\
  AEPY_PARAM \
   "    ra : The right ascension, in degrees.\n"\
   "   dec : The declination, in degrees.\n"\
   "radius : The radius; the units of radius will determine the units of the\n"\
   "         output.\n"\
  AEPY_RETURN \
   "A list, of length three, containing the rectangular vector.\n"\
  AEPY_SEE \
   "rect_to_polar()\n"

#define AEPY_DS_radec_to_altaz \
  AEPY_METHOD_NAME \
   "radec_to_altaz(last, glat, ra, dec)\n"\
  AEPY_SHORT_DESC \
   "Convert ra/dec to alt/az.\n"\
  AEPY_PARAM \
   "last : The local aparent sidereal time, in degrees.\n"\
   "glat : The geodetic latitude, in degrees.\n"\
   "  ra : The right ascension, in degrees.\n"\
   " dec : The declination, in degrees.\n"\
  AEPY_RETURN \
   "A tuple of (alt, az), in degrees.\n"\
  AEPY_SEE \
   "altaz_to_radec()\n"

#define AEPY_DS_altaz_to_radec \
  AEPY_METHOD_NAME \
   "altaz_to_radec(last, glat, alt, az)\n"\
  AEPY_SHORT_DESC \
   "Convert alt/az to ra/dec.\n"\
  AEPY_PARAM \
   "last : The local aparent sidereal time, in degrees.\n"\
   "glat : The geodetic latitude, in degrees.\n"\
   " alt : The altitude, in degrees.\n"\
   "  az : The azimuth, in degrees.\n"\
  AEPY_RETURN \
   "A tuple of (ra, dec), in degrees.\n"\
  AEPY_SEE \
   "radec_to_altaz()\n"

#define AEPY_DS_radec_to_gal \
  AEPY_METHOD_NAME \
   "radec_to_gal(ra, dec [, is_fk4]))\n"\
  AEPY_SHORT_DESC \
   "Convert ra/dec to l/b galactic coordinates.\n"\
  AEPY_PARAM \
   "    ra : The right ascension, in degrees.\n"\
   "   dec : The declination, in degrees.\n"\
   "is_fk4 : If true, use FK4 coordinates (B1950) rather than the default, "\
   "         which is FK5 (J2000)."\
  AEPY_RETURN \
   "A tuple of (l, b), in degrees.\n"\
  AEPY_SEE \
   "gal_to_radec()\n"

#define AEPY_DS_gal_to_radec \
  AEPY_METHOD_NAME \
   "gal_to_radec(l, b [, is_fk4]))\n"\
  AEPY_SHORT_DESC \
   "Convert galactic l/b to ra/dec equatorial coordinates.\n"\
  AEPY_PARAM \
   "     l : The galactic longitude, in degrees.\n"\
   "   dec : The galactic latitude, in degrees.\n"\
   "is_fk4 : If true, use FK4 coordinates (B1950) rather than the default, "\
   "         which is FK5 (J2000)."\
  AEPY_RETURN \
   "A tuple of (ra, dec), in degrees.\n"\
  AEPY_SEE \
   "radec_to_gal()\n"

#define AEPY_DS_geocentric_lat \
  AEPY_METHOD_NAME \
   "geocentric_lat(glat, height)\n"\
  AEPY_SHORT_DESC \
   "Get the geocentric latitude and the distance to the earth's centre.\n"\
  AEPY_LONG_DESC \
   "This function uses the reciprocal of flattening and the earth's radius\n"\
   "(or, more precisely, semi-major axis) from from the WGS84 definition.\n"\
  AEPY_PARAM \
   "  glat : The geodetic latitude of the observer, in degrees.\n"\
   "height : The height of the observer above sea level, in metres.\n"\
  AEPY_RETURN \
   "A tuple of geocentric latitude, in degrees, and distance to the earth's\n"\
   "centre, in earth radii: (tlat, trho).\n"\
  AEPY_SEE \
   "FLAT, R_EARTH"

#define AEPY_DS_delta_q \
  AEPY_METHOD_NAME \
   "delta_q(q0, q1)\n"\
  AEPY_SHORT_DESC \
   "Convert change in rectangular coordinatates to change in ra/dec.\n"\
  AEPY_LONG_DESC \
   "For changes greater than about 0.1 degree, the coordinates are converted\n"\
   "directly to ra and dec and the results subtracted.  For small changes,\n"\
   "the change is calculated to first order by differentiating the\n"\
   "conversion formulae.\n"\
  AEPY_PARAM \
   "q0 : The first vector, as a list of length three.\n"\
   "q1 : The second vector, as a list of length three.\n"\
  AEPY_RETURN \
   "The change in right ascension and declination, in degrees: (dra, ddec).\n"\

#define AEPY_DS_decimal \
  AEPY_METHOD_NAME \
   "decimal(d, m, s [, hour])\n"\
  AEPY_SHORT_DESC \
   "Convert degrees/hours, (arc)minutes and (arc)seconds to a decimal angle.\n"\
  AEPY_LONG_DESC \
   "Note that if _any_ of parameters d, m, s are negative, the whole angle\n"\
   "will be considered negative, and the absolute values of d, m, s will be\n"\
   "used to determine the absolute value of the decimal value.\n"\
  AEPY_PARAM \
   "   d - The degrees or hours; any decimals will be truncated.\n"\
   "   m - The arcminutes or minutes; any decimals will be truncated.\n"\
   "   s - The arcseconds or seconds.\n"\
   "hour - If true, then interpret the arguments as hours, minutes and\n"\
   "       seconds, rather than the default degrees, arcminutes and\n"\
   "       arcseconds. Optional.\n"\
  AEPY_RETURN \
   "The decimal angle. Note that no range-checking is done: if an angle\n"\
   "outside the range 0-360 degrees is passed, an angle outside that range\n"\
   "will be returned.\n"\
  AEPY_SEE \
   "dms(), dms_str(), hms(), hms_str()\n"

#define AEPY_DS_dms \
  AEPY_METHOD_NAME \
   "dms(angle)\n"\
  AEPY_SHORT_DESC \
   "Calculate the degrees, minutes and seconds of an angle.\n"\
  AEPY_PARAM \
   "angle : The angle, in degrees.\n"\
  AEPY_RETURN \
   "A tuple of (d, m, s). All values are floats, though degrees and minutes\n"\
   "are always rounded to integer values.\n"\
  AEPY_SEE \
   "dms_str(), hms(), hms_str(), decimal()\n"

#define AEPY_DS_hms \
  AEPY_METHOD_NAME \
   "hms(angle)\n"\
  AEPY_SHORT_DESC \
   "Calculate the hours, minutes and seconds of an angle.\n"\
  AEPY_PARAM \
   "angle : The angle, in degrees.\n"\
  AEPY_RETURN \
   "A tuple of (h, m, s). All values are floats, though hours and minutes\n"\
   "are always rounded to integer values.\n"\
  AEPY_SEE \
   "dms(), dms_str(), hms_str(), decimal()\n"

#define AEPY_DS_dms_str \
  AEPY_METHOD_NAME \
   "dms_str(angle [, fixed, sign])\n"\
  AEPY_SHORT_DESC \
   "Return an angle in a string, expressed in degrees, minutes and seconds.\n"\
  AEPY_PARAM \
   "angle : The angle to express.\n"\
   "fixed : If True, then return a string of fixed width (useful for\n"\
   "        printing columns of data.\n"\
   " sign : If True, then always print the sign, even when positive.\n"\
  AEPY_RETURN \
   "A string of the format ±DDDd MM' SS.SS\".\n"\
  AEPY_SEE \
   "dms(), hms(), hms_str(), decimal()\n"

#define AEPY_DS_hms_str \
  AEPY_METHOD_NAME \
   "hms_str(angle [, fixed, sign])\n"\
  AEPY_SHORT_DESC \
   "Return an angle in a string, expressed in hours, minutes and seconds.\n"\
  AEPY_PARAM \
   "angle : The angle to express.\n"\
   "fixed : If True, then return a string of fixed width (useful for\n"\
   "        printing columns of data.\n"\
   " sign : If True, then always print the sign, even when positive.\n"\
  AEPY_RETURN \
   "A string of the format ±HHh MMm SS.SSs.\n"\
  AEPY_SEE \
   "dms(), dms_str(), hms(), decimal()\n"

//==============================================================================
// BINDINGS FOR PRECESS.C, NUTATE.C AND EPSILON.C.
//==============================================================================

#define AEPY_DS_precess \
  AEPY_METHOD_NAME \
   "precess(r, jd_tt [, to_j2000])\n"\
  AEPY_SHORT_DESC \
   "Precess a coordinate.\n"\
  AEPY_LONG_DESC \
   "By default, the coordinate is precessed from J2000 to jd_tt; to precess\n"\
   "in the opposite direction, set optional parameter to_j2000 = True. To\n"\
   "precess from jd_1 to jd_2, first go from jd_1 to J2000, and then from\n"\
   "J2000 to jd_2.\n"\
  AEPY_PARAM \
   "       r - A list, of length 3, giving the rectangular coordinates to be\n"\
   "   jd_tt - The Julian date in TT.\n"\
   "to_j2000 - If to_j2000 = True, then precess from jd_tt to J2000;\n"\
   "           otherwise, precess from J2000 to jd_tt. Optional.\n"\
  AEPY_RETURN \
   "A list, of length 3, containing the precessed, rectangular coordinates.\n"\
  AEPY_SEE \
   "J2000, rect_to_polar(), polar_to_rect()\n"

#define AEPY_DS_nutation_lon_ob \
  AEPY_METHOD_NAME \
   "nutation_lon_ob(jd_tt)\n"\
  AEPY_SHORT_DESC \
   "Calculate the nutation in longitude and oblation.\n"\
  AEPY_LONG_DESC \
   "This program implements all of the 1980 IAU nutation series.  Results\n"\
   "checked at 100 points against the 1986 AA; all agreed.\n"\
  AEPY_PARAM \
   "jd_tt - The Julian date in TT.\n"\
  AEPY_RETURN \
   "The tuple (nutl, nuto), i.e., the nutation in longitude and oblation,\n"\
   "respectively, each in seconds of arc.\n"\
  AEPY_SEE \
   "nutate(), epsilon()\n"

#define AEPY_DS_nutate \
  AEPY_METHOD_NAME \
   "nutate(r, {(nutl, nuto, eps) | jd_tt} [, to_j2000])\n"\
  AEPY_SHORT_DESC \
   "Correct for nutation.\n"\
  AEPY_LONG_DESC \
   "This function can receive either nutation parameters (nutl, nuto, eps)\n"\
   "or a Julian date (jd_tt). Note that the latter is much slower, since\n"\
   "the nutation parameters need to be calculated by the function. If making\n"\
   "repeated calls at roughly the same epoch, consider using\n"\
   "nutation_lon_ob() and epsilon() to calculate the nutuation parameters\n"\
   "once, and then reusing them multiple times with this function.\n"\
   "\n"\
   "By default, the coordinate is nutated from J2000 to jd_tt; to nutate\n"\
   "in the opposite direction, set optional parameter to_j2000 = True. To\n"\
   "nutate from jd_1 to jd_2, first go from jd_1 to J2000, and then from\n"\
   "J2000 to jd_2.\n"\
  AEPY_PARAM \
   "       r : The rectangular coordinate to be nutated.\n"\
   "    nutl : The nutation in longitude in seconds of arc.\n"\
   "    nuto : The nutation in oblation in seconds of arc.\n"\
   "     eps : The obliquity of the ecliptic in seconds of arc.\n"\
   "   jd_tt : The Julian date in TT. Only pass this if (nutl, nuto, eps)\n"\
   "           are not being passed.\n"\
   "to_j2000 : If to_j2000 = True, then nutate from jd_tt to J2000;\n"\
   "           otherwise, nutate from J2000 to jd_tt. Optional.\n"\
  AEPY_RETURN \
   "A list, of length 3, containing the nutated, rectangular coordinates.\n"\
  AEPY_SEE \
   "nutation_lon_ob(), epsilon(), J2000, rect_to_polar(), polar_to_rect()\n"

#define AEPY_DS_epsilon \
  AEPY_METHOD_NAME \
   "epsilon(jd_tt)\n"\
  AEPY_SHORT_DESC \
   "Calculate the obliquity of the ecliptic.\n"\
  AEPY_PARAM \
   "jd_tt : The Julian date in TT.\n"\
  AEPY_RETURN \
   "The obliquity, in seconds of arc.\n"\
  AEPY_SEE \
   "nutate(), nutation_lon_ob()\n"

//==============================================================================
// BINDINGS FOR CONSTEL.C
//==============================================================================

#define AEPY_DS_constellation \
  AEPY_METHOD_NAME \
   "constellation(ra, dec [, epoch, abbriev])\n"\
  AEPY_SHORT_DESC \
   "Given a sky coordinate, find the constellation it is in.\n"\
  AEPY_PARAM \
   "     ra : The right ascension, in degrees.\n"\
   "    dec : The declination, in degrees.\n"\
   "  epoch : The epoch (in Julian days) of the ra/dec; default is J2000.\n"\
   "          Optional.\n"\
   "abbriev : If true, then return the three-letter abbrieviation of the\n"\
   "          constellation. Optional.\n"\
  AEPY_RETURN \
   "A string containing the name of the constellation (or its three-letter\n"\
   "abbreviation if abbriev = true); if the constellation cannot be found\n"\
   "for some reason, \"unknown\" is returned.\n"\

//==============================================================================
// BINDINGS FOR ANNUAL_AB.C
//==============================================================================

#define AEPY_DS_annual_aberration \
  AEPY_METHOD_NAME \
   "annual_aberration(v, p [, add])\n"\
  AEPY_SHORT_DESC \
   "Correct for the annual aberration.\n"\
  AEPY_PARAM \
   "  v : The heliocentric rectangular velocity of Earth in AU per day.\n"\
   "  p : A unit vector pointing from the earth to the object; the corrected\n"\
   "      position is returned in this parameter.\n"\
   "add : If true, then add the aberration effect; otherwise remove it\n"\
   "      (default). Adding is appropriate, for example, if one wishes to\n"\
   "      compute the apparent position of a star from a catalogue position,\n"\
   "      while removing is appropriate, for example, if one wants to\n"\
   "      compute the true position of an observed star. Optional.\n"\
  AEPY_RETURN \
   "A list, of length 3, being the vector p, corrected for aberration.\n"\
  AEPY_SEE \
   "diurnal_aberration()\n"\

#define AEPY_DS_diurnal_aberration \
  AEPY_METHOD_NAME \
   "diurnal_aberration(last, tlat, trho, ra, dec [, add])\n"\
  AEPY_SHORT_DESC \
   "Correct for diurnal aberration.\n"\
  AEPY_LONG_DESC \
   "This formula is less rigorous than the method used for annual aberration\n"\
   "(see ae_annual_aberration()).  However, the correction is much smaller,\n"\
   "so it is adequate.\n"\
  AEPY_PARAM \
   "last : The local apparent sidereal time, in degrees.\n"\
   "tlat : The geocentric latitude of the observer, in degrees.\n"\
   "trho : The distance from the centre of the earth to the observer, in\n"\
   "       earth radii.\n"\
   "  ra : The right ascension to be corrected, in degrees.\n"\
   " dec : The declination to be corrected, in degrees.\n"\
   " add : If true, then add the aberration effect; otherwise remove it\n"\
   "       (default). Adding is appropriate, for example, if one wishes to\n"\
   "       compute the apparent position of a star from a catalogue\n"\
   "       position, while removing is appropriate, for example, if one\n"\
   "       wants to compute the true position of an observed star. Optional.\n"\
  AEPY_RETURN \
   "The corrected right ascension and declination as a tuple (ra, dec), in\n"\
   "degrees.\n"\
  AEPY_SEE \
   "geocentric_lat(), annual_aberration()\n"

#define AEPY_DS_diurnal_parallax \
  AEPY_METHOD_NAME \
   "diurnal_parallax(last, tlat, trho, dist, ra, dec)\n"\
  AEPY_SHORT_DESC \
   "Correct for diurnal parallax.\n"\
  AEPY_LONG_DESC \
   "This function does not bother to calculate anything unless the\n"\
   "equatorial horizontal parallax is at least 0.005\".\n"\
  AEPY_PARAM \
   "last : The local apparent sidereal time, in degrees.\n"\
   "tlat : The geocentric latitude of the observer, in degrees.\n"\
   "trho : The distance from the centre of the earth to the observer, in\n"\
   "       earth radii.\n"\
   "dist : The earth-object distance, in AU.\n"\
   "  ra : The right ascension to be corrected, in degrees.\n"\
   " dec : The declination to be corrected, in degrees.\n"\
  AEPY_RETURN \
   "The corrected right ascension and declination as a tuple (ra, dec), in\n"\
   "degrees.\n"\
  AEPY_SEE \
   "geocentric_lat()\n"

//==============================================================================
// BINDINGS FOR RELATIVITY.C
//==============================================================================

#define AEPY_DS_relativity \
  AEPY_METHOD_NAME \
   "relativity(p, q, o)\n"\
  AEPY_SHORT_DESC \
   "Correct for light deflection due to solar gravitation.\n"\
  AEPY_LONG_DESC \
   "Note that this blows up if the object is very near (or is) the sun.\n"\
   "Therefore, if the sun-object or sun-object distance is zero, no\n"\
   "correction is made.\n"\
  AEPY_PARAM \
   "p : The unit vector from observer to an object.\n"\
   "q : The heliocentric ecliptic rectangular coordinates of the object.\n"\
   "o : The heliocentric ecliptic rectangular coordinates of the observer.\n"\
  AEPY_RETURN \
   "A list, of length 3, of the corrected unit vector p.\n"\

//==============================================================================
// BINDINGS FOR TOPOCENTRIC.C AND REFRAC.C.
//==============================================================================

#define AEPY_DS_refrac_visible \
  AEPY_METHOD_NAME \
   "refrac_visible(alt, pres, temp)\n"\
  AEPY_SHORT_DESC \
   "Atmospheric refraction correction (in visible bands).\n"\
  AEPY_LONG_DESC \
   "This function is suitable for passing to ae_topocentric().\n"\
   "\n"\
   "For high altitude angle, AA page B61 is used.  The accuracy is `usually\n"\
   "about 0.1 arcsecond'.\n"\
   "\n"\
   "The formula for low altitude is from the Almanac for Computers. It gives\n"\
   "the correction for observed altitude, so has to be inverted numerically\n"\
   "to get the observed from the true.  The accuracy is about 0.2' for\n"\
   "-20C < T < +40C and 970mb < P < 1050mb.\n"\
  AEPY_PARAM \
   " alt : The altitude of the observation, in degrees.\n"\
   "pres : The atmospheric pressure, in millibar.\n"\
   "temp : The temperature, in degrees centigrade.\n"\
  AEPY_RETURN \
   "The correction in degrees to be added to true altitude to obtain\n"\
   "apparent altitude.\n"\
  AEPY_SEE \
   "topocentric(), refrac_ulich()\n"

#define AEPY_DS_refrac_ulich \
  AEPY_METHOD_NAME \
   "refrac_ulich(alt, pres, temp, humid)\n"\
  AEPY_SHORT_DESC \
   "Atmospheric refraction correction (in millimetre bands).\n"\
  AEPY_LONG_DESC \
   "This function is suitable for passing to ae_topocentric().\n"\
   "\n"\
   "This function uses the equations from:\n"\
   "    B. L. Ulich.  `Millimeter Wave Radio Telescopes: Gain and Pointing\n"\
   "    Charactersitics.' International Journal of Infrared and Millimeter\n"\
   "    Waves, 2, 2 (1981).\n"\
   "\n"\
   "Ulich does not cite which millimetre wavelength ranges the equations are\n"\
   "expected to hold good for. However, the ALMA Memo 366 does analysis on\n"\
   "another model and finds it good to 2\% up to frequencies of 1000 GHz.\n"\
   "(The worst is at 500 GHz.)\n"\
   "\n"\
   "Ulich claims accuracy to 2\" above altitudes of 3 degrees.\n"\
  AEPY_PARAM \
   "  alt : The altitude of the observation, in degrees.\n"\
   " pres : The atmospheric pressure, in millibar.\n"\
   " temp : The temperature, in degrees centigrade.\n"\
   "humid : The surface humidity, as a percentage.\n"\
  AEPY_RETURN \
   "The correction in degrees to be added to true altitude to obtain\n"\
   "apparent altitude.\n"\
  AEPY_SEE \
   "topocentric(), refrac_visible()\n"

#define AEPY_DS_topocentric \
  AEPY_METHOD_NAME \
   "topocentric(jd_ut1, glat, tlon, ra, dec [, dist, jd_tt, tlat, trho, \\\n"\
   "            refrac])\n"\
  AEPY_SHORT_DESC \
   "Apply corrections to obtain topocentric altitude and azimuth."\
  AEPY_LONG_DESC \
   "To go from geocentric to apparent, topocentric coordinates, the diurnal\n"\
   "aberration and parallax should, in many cases, be accounted for, as well\n"\
   "as atmospheric refraction. This function can do each of these.\n"\
   "\n"\
   "Note that, for repeated calls at the same, or similar, epochs, it may\n"\
   "be more efficient to make these corrections `by hand', as this function\n"\
   "needs to compute nutation factors; calling it many times can therefore\n"\
   "be inefficient.\n"\
   "\n"\
   "If refrac is not passed, then no atmospheric refraction correction is\n"\
   "performed.\n"\
   "\n"\
   "If dist is not passed, then no diurnal parallax correction is performed.\n"\
   "\n"\
   "The other optional arguments can be passed to speed up execution of this\n"\
   "function, which otherwise needs to calculate them itself.\n"\
  AEPY_PARAM \
   "jd_ut1 : The Julian date in UT1 measure.\n"\
   "  glat : The geodetic latitude of the observer, in degrees.\n"\
   "  tlon : The longitude of the observer, in degrees.\n"\
   "    ra : The right ascension of the object, in degrees.\n"\
   "   dec : The declination of the object, in degrees.\n"\
   "  dist : The distance to the object being observed, in AU. Required for\n"\
   "         diurnal parallax correction. Optional.\n"\
   "refrac : A tuple of the form (func, ...). The first argument is a\n"\
   "         function that can calculate atmospheric refraction. Its first\n"\
   "         argument should be the altitude, which this function will\n"\
   "         provide. Its other arguments (such as temperature or pressure)\n"\
   "         should be given by the remaining arguments of the tuple refrac.\n"\
   "         Example: if the refraction function is do_refrac(alt, a, b, c),\n"\
   "         then refrac should be (do_refrac, a, b, c). Optional.\n"\
   " jd_tt : The Julian date in TT measure. Optional.\n"\
   "  tlat : The geocentric latitude of the observer, in degrees. Optional.\n"\
   "  trho : The distance from the centre of the earth to the observer, in\n"\
   "         earth radii. Optional.\n"\
  AEPY_RETURN \
   "The tuple (ra, dec), containing the topographic right ascension and\n"\
   "declination, in degrees.\n"\
  AEPY_SEE \
   "diurnal_aberration(), diurnal_parallax(), geocentric_lat(), delta_t(),\n"\
   "refrac_visible(), refrac_ulich()\n"



//==============================================================================
// BINDINGS FOR PHYSICAL.C
//==============================================================================

#define AEPY_DS_disc_semiminor \
  AEPY_METHOD_NAME \
   "disc_semiminor(a, c, latitude)\n"\
  AEPY_SHORT_DESC \
   "Calculate the semi-minor axis of an observed disc.\n"\
  AEPY_LONG_DESC \
   "When an oblate spheroid is viewed, it appears as an ellipse (to a\n"\
   "sufficiently distant observer).  This routine calculates the semi-minor\n"\
   "axis (i.e., polar radius) of the projected ellipse.\n"\
  AEPY_PARAM \
   "       a : The major axis (i.e., equatorial radius) of the planet.\n"\
   "       c : The minor axis (i.e., polar radius) of the planet.\n"\
   "latitude : The sub-observer latitude, in degrees.\n"\
  AEPY_RETURN \
   "The semi-minor axis of the projected ellipse.  Note that the semi-major\n"\
   "axis of the ellipse is the same as the input value.\n"\

#define AEPY_DS_disc_solid_angle \
  AEPY_METHOD_NAME \
   "disc_solid_angle(jd_ut1, obs, obj, phys)\n"\
  AEPY_SHORT_DESC \
   "Calculate the solid angle of an observed disc.\n"\
  AEPY_LONG_DESC \
   "Given the semi-minor and -major axes along with the distance, the\n"\
   "solid angle is simply:\n"\
   "  pi x a x b / d^2\n"\
   "\n"\
   "This routine determines all of the requisite figures and then performs\n"\
   "this calculation."\
   "\n"\
   "Note that this function cannot take JPL ephemerides, only 'aephem.orb'\n"\
   "objects.\n"\
  AEPY_PARAM \
   "jd_ut1 : The Julian date in UT1.\n"\
   "   obs : The orbital elements of the observer.\n"\
   "   obj : The orbital elements of the object being observed.\n"\
   "  phys : The physical elements of the object being observed.\n"\
  AEPY_RETURN \
   "The solid angle of the disc, in steradians.\n"\
  AEPY_SEE \
   "disc_semiminor()\n"

#define AEPY_DS_subobs_point \
  AEPY_METHOD_NAME \
   "subobs_point(jd_ut1, o_orb, q_orb, phys)\n"\
  AEPY_SHORT_DESC \
   "Get the subobserver point of a body.\n"\
  AEPY_LONG_DESC \
   "The subobserver point is the latitude and longitude on an observed body\n"\
   "which appears in the centre of the disc seen by the observer.\n"\
   "\n"\
   "Note that this function cannot take JPL ephemerides, only 'aephem.orb'\n"\
   "objects.\n"\
   "\n"\
   "WARNING:  the longitude computed by this function is, for some reason,\n"\
   "only accurate to tens of arc minutes when compared to HORIZONS (and I\n"\
   "don't know why). Latitude is accurate to a fraction of an arc minute.\n"\
  AEPY_PARAM \
   "jd_ut1 : The Julian date in UT1.\n"\
   "   obs : The orbital elements of the observer.\n"\
   "   obj : The orbital elements of the object being observed.\n"\
   "  phys : The physical elements of the object being observed.\n"\
  AEPY_RETURN \
   "A tuple containing the subobserver latitude, longitude, both in degrees\n"\
   "and the distance between observer and observed body, in AU.\n"

/*
#define AEPY_DS_disc_semiminor \
  AEPY_METHOD_NAME \
  AEPY_SHORT_DESC \
  AEPY_LONG_DESC \
  AEPY_PARAM \
  AEPY_RETURN \
  AEPY_SEE \
*/


//==============================================================================
// OTHER MAIN METHODS
//==============================================================================

#define AEPY_DS_geocentric \
  AEPY_METHOD_NAME \
   "geocentric(jd_tt, obs, obj)\n"\
  AEPY_SHORT_DESC \
   "Get the geocentric coordinates of an object.\n"\
  AEPY_LONG_DESC \
   "This routine returns the coordinates of a body as seen from the centre\n"\
   "of the observing body. For highest accuracy, the results should be\n"\
   "passed to topocentric() to get the coordinates as seen from a specific\n"\
   "point on the observing body.\n"\
   "\n"\
   "Most of the time, the 'obs' parameter will be an ephemeris object for\n"\
   "the earth, but it is also possible to get the coordinates of an object\n"\
   "as seen from another solar system body.\n"\
   "\n"\
   "Corrections done during the reduction are:  light deflection, annual\n"\
   "aberration, precession and nutation.\n"\
  AEPY_PARAM \
   "jd_tt : The Julian date in TT measure.\n"\
   "  obs : An ephemeris object for the observer. It should be one of:\n"\
   "        - An 'aephem.orb' object, such as 'aephem.p_orb[\"earth\"]'.\n"\
   "        - A tuple of an 'aephem.jpl' object and the object number, such\n"\
   "          as '(jpl_obj, aephem.SS_EARTH)'.\n"\
   "  obj : An ephemeris object for the body being observed. It should be\n"\
   "        one of:\n"\
   "        - An 'aephem.orb' object, such as 'aephem.p_orb[\"saturn\"]'.\n"\
   "        - A tuple of an 'aephem.jpl' object and the body number, such\n"\
   "          as '(jpl_obj, aephem.SS_SATURN)'.\n"\
   "        - An 'aephem.star' object, for a body from an extra-Solar System\n"\
   "          catalogue.\n"\
   "        - The string 'sun' or 'moon'; these have special calculations\n"\
   "          when using orbital elements. If a JPL ephemeris is being used\n"\
   "          then you can should the sun and moon ephemeris in it instead\n"\
   "          of this.\n"\
  AEPY_RETURN \
   "The tuple (ra, dec, dist), containing the topographic right ascension\n"\
   "declination, in degrees, and the distance to the body, in AU. For\n"\
   "'aephem.star' objects, the distance returned is always -1.\n"\
  AEPY_SEE \
   "topocentric()\n"

#define AEPY_DS_geocentric_from_helio \
  AEPY_METHOD_NAME \
   "geocentric_from_helio(jd_tt, o, v_o, q)\n"\
  AEPY_SHORT_DESC \
   "Given heliocentric coordinates, reduce to geocentric coordinates.\n"\
  AEPY_LONG_DESC \
   "The function geocentric() is more suitable for most applications; only\n"\
   "use this function if you need for some reason to reduce directly from\n"\
   "heliocentric coordinates.\n"\
   "\n"\
   "Unlike geocentric(), this function cannot account for light travel-time.\n"\
   "The input position of the planet is assumed to be already corrected.\n"\
   "\n"\
   "Corrections done during the reduction are:  light deflection, annual\n"\
   "aberration, precession and nutation.\n"\
  AEPY_PARAM \
   " jd_tt : The Julian date in TT measure.\n"\
   "     o : The heliocentric rectangular coordinates of the observer, in\n"\
   "         AU.\n"\
   "   v_e : The heliocentric velocity of the observer, in AU per day.\n"\
   "     q : The heliocentric rectangulage coordinates of the observed\n"\
   "         body, in AU.\n"\
  AEPY_RETURN \
   "A tuple of (ra, dec, distance), with ra and dec in degrees and distance\n"\
   "in AU.\n"\
  AEPY_SEE \
   "geocentric(), orb.coord(), jpl.coord()\n"

#define AEPY_DS_light_t \
  AEPY_METHOD_NAME \
   "light_t(jd_tt, o, obj)\n"\
  AEPY_SHORT_DESC \
   "Get the apparent position of an object, corrected for light-travel time.\n"\
  AEPY_LONG_DESC \
   "The routine does three iterations of light-time correction.\n"\
   "Gravitational retardation from the sun is neglected.\n"\
  AEPY_PARAM \
   " jd_tt : The Julian date in TT measure.\n"\
   "     o : The heliocentric rectangular coordinates of the observer, in\n"\
   "         AU.\n"\
   "   obj : An ephemeris object for the body being observed. It should be\n"\
   "         one of:\n"\
   "         - An 'aephem.orb' object, such as 'aephem.p_orb[\"mars\"]'.\n"\
   "         - A tuple of an 'aephem.jpl' object and the body number, such\n"\
   "           as '(jpl_obj, aephem.SS_SATURN)'.\n"\
  AEPY_RETURN \
   "A list of the heliocentric, rectangular coordinates of the object, in\n"\
   "AU.\n"

/*#define AEPY_DS_geocentric_from_helio \
  AEPY_METHOD_NAME \
  AEPY_SHORT_DESC \
  AEPY_LONG_DESC \
  AEPY_PARAM \
  AEPY_RETURN \
  AEPY_SEE \*/

//==============================================================================
// HELP FOR STAR OBJECT
//==============================================================================

#define AEPY_DS_star \
   "star([cat, star])\n"\
   "\n"\
   "A star object contains information on an extra-solar body---see the data\n"\
   "descriptors. The user can either fill these in himself, or he can direct\n"\
   "for them to be loaded from a catalogue via read_from_cat(). The\n"\
   "documentation for this method describes the catalogue file format. The\n"\
   "data can also optionally be read in via the object prototype, which\n"\
   "takes the same arguments as read_from_cat().\n"\
   "\n"\
   "There is one additional method: to convert from the FK4 to FK5 system,\n"\
   "use fk4_to_fk5().\n"\

#define AEPY_DS_star_read_from_cat \
  AEPYO_METHOD_NAME \
       "read_from_cat(cat, star)\n"\
  AEPYO_SHORT_DESC \
       "Parse star positional data from a catalogue file.\n"\
  AEPYO_LONG_DESC \
       "The catalogue file should be a columnated ASCII file.  Lines\n"\
       "starting with a hash ('#') are treated as comments.  The columns\n"\
       "are:\n"\
       "o epoch\n"\
       "o ra hours\n"\
       "o ra minutes\n"\
       "o ra second\n"\
       "o dec degrees\n"\
       "o dec arcminutes\n"\
       "o dec arcseconds\n"\
       "o mu ra (in seconds/century)\n"\
       "o mu dec (in arcseconds/century)\n"\
       "o v (in km/s)\n"\
       "o px (in seconds of arc)\n"\
       "o visual magnitude\n"\
       "o name\n"\
  AEPYO_PARAM \
       " cat : The full path to the catalogue file.\n"\
       "star : The name of the star, as it appears in the catalogue.\n"\
  AEPYO_RETURN \
       "None.\n"\

#define AEPY_DS_star_fk4_to_fk5 \
  AEPYO_METHOD_NAME \
       "fk4_to_fk5()\n"\
  AEPYO_SHORT_DESC \
       "Convert FK4 B1950.0 coordinates to FK5 J2000.0 coordinates.\n"\
  AEPYO_LONG_DESC \
       "It is assumed that the beginning epoch is B1950: this method does\n"\
       "not check that this indeed the case. The epoch will be updated to\n"\
       "J2000.\n"\
  AEPYO_RETURN \
       "None.\n"\

//==============================================================================
// HELP FOR ORB OBJECT
//==============================================================================

#define AEPY_DS_orb \
   "orb([cat, orb])\n"\
   "\n"\
   "An orb object contains information on an orbit---see the data\n"\
   "descriptors. The user can either fill these in himself, or he can direct\n"\
   "for them to be loaded from a catalogue file via read_from_cat(). The\n"\
   "documentation for this method describes the catalogue file format. The\n"\
   "data can also optionally be read in via the object prototype, which\n"\
   "takes the same arguments as read_from_cat().\n"\
   "\n"\
   "The method kepler() returns heliocentric, rectangular coordinates of the\n"\
   "orbit at a given time.\n"\
   "\n"\
   "There are built-in orb objects for the planets, pluto and Earth's moon,\n"\
   "available in the 'p_orb' dictionary as p_orb['mercury'], p_orb['venus']\n"\
   "and so on. These objects have perturbation tables which make them much \n"\
   "more accurate than simple Keplerian elements; they give results which \n"\
   "are within less than a second of arc from the JPL ephemerides.\n"\
   "\n"\
   "Note that currently in the python version, other 'orb' objects cannot\n"\
   "be created or loaded from file with perturbation tables besides those\n"\
   "in the 'planet' dictionary."

#define AEPY_DS_orb_read_from_cat \
  AEPYO_METHOD_NAME \
       "read_from_cat(cat, orb)\n"\
  AEPYO_SHORT_DESC \
       "Parse orbital data from a catalogue file.\n"\
  AEPYO_LONG_DESC \
       "The catalogue file should be a columnated ASCII file.  Lines\n"\
       "starting with a hash ('#') are treated as comments.  The columns\n"\
       "are:\n"\
       "o epoch\n"\
       "o inclination (in degrees)\n"\
       "o longitude of the ascending node (in degrees)\n"\
       "o argument of the perihelion (in degrees)\n"\
       "o mean distance (in AU)\n"\
       "o daily motion (in AU per day)\n"\
       "o eccentricity\n"\
       "o mean anomaly (in degrees)\n"\
       "o epoch of equinox and ecliptic\n"\
       "o name\n"\
  AEPYO_PARAM \
       "cat : The full path to the catalogue file.\n"\
       "orb : The name of the orbit, as it appears in the catalogue.\n"\
  AEPYO_RETURN \
       "None.\n"

#define AEPY_DS_orb_kepler \
  AEPYO_METHOD_NAME \
       "kepler(jd_tt [, velocity])\n"\
  AEPYO_SHORT_DESC \
       "Find a Keplerian orbit, and, optionally, its velocity.\n"\
  AEPYO_LONG_DESC \
       "This routine solves for a Keplerian orbit, given orbital parameters\n"\
       "and the time. If the 'velocity' parameter is true, also get the\n"\
       "velocity of the orbit.\n"\
       "\n"\
       "The routine detects several cases of given orbital elements.  If a\n"\
       "program for perturbations is pointed to, it is called to calculate\n"\
       "all the elements. If there is no program, then the mean longitude is\n"\
       "calculated from the mean anomaly and daily motion.  If the daily\n"\
       "motion is not given, it is calculated by Kepler's law.  If the\n"\
       "eccentricity is given to be 1.0, it means that meandistance is\n"\
       "really the perihelion distance, as in a comet specification, and the\n"\
       "orbit is parabolic.\n"\
  AEPYO_PARAM \
       "   jd_tt : The Julian date in TT.\n"\
       "velocity : If true (default is false), then also return the\n"\
       "           velocity.\n"\
  AEPYO_RETURN \
       "A list of the heliocentric, rectangular coordinates of the\n"\
       "object, in AU; if 'velocity' is true, then a tuple is returned, the\n"\
       "first item being the position, and the second the velocity (also\n"\
       "in heliocentric, rectangular coordinates), in AU per day.\n"\
  AEPYO_SEE \
       "jpl.get_coords()\n"

//==============================================================================
// HELP FOR JPL OBJECT
//==============================================================================

#define AEPY_DS_jpl \
   "jpl(path [, ascii, header_path, trust_header])\n"\
   "\n"\
   "A JPL object is a pointer to a JPL ephemeris file, either binary or\n"\
   "ASCII. Once initialised, it can be passed to various methods to\n"\
   "calculate astronomical positions.\n"\
  AEPY_PARAM \
   "        path : The location of the ephemeris file.\n"\
   "       ascii : The default is to read a binary ephemeris file. Set\n"\
   "               'ascii' to true if you wish to read an ASCII ephemeris.\n"\
   "               Optional.\n"\
   " header_path : The location of the header file for an ASCII ephemeris,\n"\
   "               if it is separate from the main data file. This parameter\n"\
   "               is only used if 'ascii' is true. Optional.\n"\
   "trust_header : Since ASCII ephemerides are sometimes broken into\n"\
   "               multiple files, the initialiser normally searches the\n"\
 "               file for the date range. If, however, you wish to \"trust\"\n"\
   "               the header and assume the range listed there is correct,\n"\
   "               set 'trust_header' to true. This parameter is only used\n"\
   "               if 'ascii' is true. Optional.\n"

#define AEPY_DS_jpl_get_coords \
  AEPYO_METHOD_NAME \
       "get_coords(jd_tt, obj_num [, not_planetary)\n"\
  AEPYO_SHORT_DESC \
       "Get object positions from a JPL ephemeris file.\n"\
  AEPYO_PARAM \
       "        jd_tt : The Julian date in TT.\n"\
       "      obj_num : The object for which to get coordinates.  If a\n"\
       "                planetary ephemeris file is being used, one can use\n"\
       "                the predefined constants SS_SUN, SS_MERCURY,\n"\
       "                SS_VENUS, and so on.\n"\
       "not_planetary : If true (default false), assume these are not\n"\
       "                planetary ephemerides, in which case barycentric\n"\
       "                coordinates are returned. Otherwise, it will be\n"\
       "                assumed that a Solar System ephemeris file is being\n"\
       "                used and heliocentric coordinates will be returned.\n"\
  AEPYO_RETURN \
       "A tuple of the rectangular position and velocity coordinates. Units\n"\
       "are AU and AU per day, unless nutations or librations are being\n"\
       "returned, in which case units are seconds of arc and seconds of arc\n"\
       "per day.\n"\
  AEPYO_SEE \
       "orb.kepler()\n"

//==============================================================================
// HELP FOR PHYS OBJECT
//==============================================================================

#define AEPY_DS_phys \
   "phys()\n"\
   "\n"\
   "A 'phys' object contains information on the physical and rotational\n"\
   "parameters---see the data descriptors.\n"\
   "\n"\
   "There are built-in 'phys' objects for the planets and pluto, available\n"\
   "in the 'p_phys' dictionary as p_phys['mercury'], p_phys['venus'] and \n"\
   "so on.\n"\
   "\n"\
   "The method pole() returns the position of the body's pole and prime\n"\
   "meridean at a given epoch. The methods flattening() and is_retrograde()\n"\
   "return some basic derived parameters on the body.\n"

#define AEPY_DS_phys_pole \
  AEPYO_METHOD_NAME \
       "pole(jd_ut1, jd_tt)\n"\
  AEPYO_SHORT_DESC \
       "Calculate the axis of a body's north pole and its prime meridean.\n"\
  AEPYO_LONG_DESC \
       "Get the precessed, nutated direction of the body's north pole; also\n"\
       "return the prime meridean.\n"\
  AEPYO_PARAM \
       "jd_ut1 : The Julian date in UT1.  For reasonable accuracy, pass a\n"\
       "         date that has been corrected for light travel time, to get\n"\
       "         the position of the north pole and meridean as seen by the\n"\
       "         observer.\n"\
       " jd_tt : The Julian date in TT.  This should not be corrected for\n"\
       "         light-travel time.\n"\
  AEPYO_RETURN \
       "The unit vector 'n', in equatorial units of date, in degrees, and\n"\
       "the prime meridean 'w', in degrees.\n"

#define AEPY_DS_phys_is_retrograde \
  AEPYO_METHOD_NAME \
       "is_retrograde()\n"\
  AEPYO_SHORT_DESC \
       "Determine whether a body's rotation is retrograde.\n"\
  AEPYO_LONG_DESC \
       "This function merely looks at the sign of attribute 'w_d': negative\n"\
       "values are retrograde.\n"\
  AEPYO_RETURN \
       "True if the rotation is retrograde; false if the rotation is\n"\
       "prograde.\n"

#define AEPY_DS_phys_flattening \
  AEPYO_METHOD_NAME \
       "flattening()\n"\
  AEPYO_SHORT_DESC \
       "Calculate a body's flattening.\n"\
  AEPYO_LONG_DESC \
       "The flattening is f = (a - b) / a, where 'a' is the equatorial\n"\
       "radius ('r_eq') and 'b' is the polar radius ('r_pole').\n"\
  AEPYO_RETURN \
       "The flattening.\n"

/*
#define AEPY_DS_orb_read_from_cat \
  AEPYO_METHOD_NAME \
  AEPYO_SHORT_DESC \
  AEPYO_LONG_DESC \
  AEPYO_PARAM \
  AEPYO_RETURN \
  AEPYO_SEE
*/

#endif
