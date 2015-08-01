//! \file aephem.h
//! Header file for AEPHEM library.
//==============================================================================
// AEPHEM - an astronomical ephemeris and reduction library.
// Copyright 2012 Adam Hincks, Canadian Institute for Theoretical Astrophysics
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

#ifndef AEPHEM_H
#define AEPHEM_H

#ifdef __cplusplus
  extern "C" {
#endif

#include <math.h>
#include <stdio.h>

// Constants.
#define AE_DTR          1.7453292519943295769e-2 //!< Degrees to radians.
#define AE_RTD          5.7295779513082320877e1  //!< Radians to degrees.
#define AE_RTS          2.0626480624709635516e5  //!< Arcseconds per radian.
#define AE_STR          4.8481368110953599359e-6 //!< Radians per arcsecond.
#define AE_STD          (1.0 / 3600.0)           //!< Degrees per arcsecond.
#define AE_STDAY        0.00001157407407407407   //!< Days per second.
#define AE_COSSUN       -0.014543897651582657    //!< Cos(90d 50').
#define AE_COSZEN       -9.8900378587411476e-3   //!< Cos(90d 34').
#define AE_J2000        2451545.0                //!< 1.5 Jan 2000.
#define AE_B1950        2433282.423              //!< 0.923 Jan 1950, Besselian.
#define AE_J1900        2415020.0                //!< 0 Jan 1900.
//! The Julian Date when the Modified Julian Date begins.
#define AE_MJD_START    2400000.5
#define AE_AU           1.49597870691e8          //!< Astronomical unit in km.
#define AE_CLIGHT       2.99792458e5             //!< Speed of light in km/s.
#define AE_S_PER_D      86400.0                  //!< Number of seconds per day.
#define AE_D_PER_S      (1.0 / AE_S_PER_D)       //!< Number of days per second.
//! The speed of light in AU / day.
#define AE_CLIGHT_AUD   (86400.0 * AE_CLIGHT / AE_AU)
#define AE_R_EARTH      6378137.0                //!< Radius of earth in m.
//! Radius of the earth in AU.
#define AE_R_EARTH_AU   (0.001 * AE_R_EARTH / AE_AU) //!< Radius of earth in AU.
#define AE_E_M_RAT      81.300585                //!< Earth/moon mass ratio.
//! The reciprocal of flattening for the earth.
//! Taken from the World Geodetic System 1984 (WGS84) definitions.
#define AE_FLAT         298.257222101
//! The temperature of the cosmic background radiation, in degrees kelvin.
#define AE_T_CMB        2.725
//! Size of the JPL binary header.
#define AE_JPL_HEAD_SIZE (5 * sizeof(double) + 41 * sizeof(int32_t))
#define AE_JPL_CONST_LEN 6      //!< Constant name lengths in JPL ephemerides.

//! Error codes that functions might return.
//! Functions will generally return negative values, e.g., -AE_RET_BAD_PATH == 
//! -1, not AE_RET_BAD_PATH == 1.
enum ae_retcode_t {
  AE_RET_BAD_PATH = 1,          //!< A path could not be accessed for I/O.
  AE_RET_UNEXPECTED_EOF,        //!< An unexpected EOF was found in a file.
  AE_RET_READ_ERROR,            //!< An unknown read error occured in a file.
  AE_RET_BAD_JPL_HEADER,        //!< The JPL header was corrupted.
  AE_RET_BAD_JPL_HEADER_KSIZE,  //!< Bad JPL header; ksize != 2 * n_coeff.
  AE_RET_BAD_JPL_HEADER_GROUP,  //!< Bad JPL header; could not find a group.
  //! Bad JPL header; differing numbers of constants and values.
  AE_RET_BAD_JPL_HEADER_NUM, 
  AE_RET_BAD_JPL_HEADER_RANGE,  //!< Bad JPL header; could not read date range.
  AE_RET_BAD_JPL_NOT_ASCII,     //!< Was expecting a JPL ASCII file.
  AE_RET_BAD_JPL_NOT_BIN,       //!< Was expecting a JPL binary file.
  AE_RET_BAD_JPL_BLOCK,         //!< Bad JPL block.
  AE_RET_BAD_JPL_CORRUPT,       //!< Data file seems corrupt.
  AE_RET_BAD_JPL_DATE,          //!< Julian date outside JPL ephemeris range.
  AE_RET_BAD_JPL_FLOAT,         //!< Unable to parse float in JPL file.
  AE_RET_BAD_JPL_CHEBY_TIME,    //!< Bad chebyshev time when gettin JPL data.
  AE_RET_BAD_OBJ_INDEX,         //!< Unknown object number.
  AE_RET_OBJ_NOT_FOUND,         //!< Object not found in a catalogue file.
  //! A function requiring cfitsio was called, but cfitsio is not available.
  AE_RET_CFITSIO_NOT_ENABLED,
  AE_RET_CFITSIO_NO_KEYWORD,    //!< A required FITS keyword did not exist.
  AE_RET_CFITSIO_BAD_KEYWORD,   //!< A required FITS keyword was corrupt.
  AE_RET_WCS_UNKNOWN_PROJ,      //!< An unknown projection type was detected.
  AE_RET_CFITSIO_ERROR,         //!< An error occurred in a call to libcfitsio.
  AE_RET_BAD_NATIVE_POLE,       //!< Could not compute projection's native pole.
  AE_RET_NO_INVERSE,            //!< A matrix was not invertible.
  //! The function does not work as it is still being developed!
  AE_RET_UNDER_CONSTRUCTION
};

//! Number of entries in #ae_ss_bodies_t.
#define AE_N_SS_BODIES  16
//! Number of bodies given in the JPL ephemeris.
#define AE_N_SS_BODIES_JPL  13
//! Enumerate Solar System bodies.
//! The first #AE_N_SS_BODIES_JPL are in the the JPL ephemeris order.  The rest
//! are derived from these.
enum ae_ss_bodies_t {
  AE_SS_MERCURY = 0,    //!< Mercury, the Winged Messenger.
  AE_SS_VENUS,          //!< Venus, the Bringer of Peace.
  AE_SS_EMBARY,         //!< Earth-moon barycentre.
  AE_SS_MARS,           //!< Mars, the Bringer of War.
  AE_SS_JUPITER,        //!< Jove, the Bringer of Jollity.
  AE_SS_SATURN,         //!< Saturn, the Bringer of Old Age.
  AE_SS_URANUS,         //!< Uranus, the Magician.
  AE_SS_NEPTUNE,        //!< Neptune, the Mystic.
  AE_SS_PLUTO,          //!< Pluto, now a dwarf planet.
  AE_SS_MOON_EMB,       //!< Our moon, relative to the Earth-moon barycentre.
  AE_SS_SUN,            //!< The sun.
  AE_SS_NUTATION,       //!< Nutation.
  AE_SS_LIBRATION,      //!< Libration.
  AE_SS_EARTH,          //!< The earth.
  AE_SS_MOON,           //!< The moon.
  AE_SS_SSBARY          //!< The Solar System Barycenter.
};

//! Define directions for precession and nutation.
enum ae_precess_direction_t {
  AE_FROM_J2000 = -1,   //!< Transform from J2000 to JD.
  AE_TO_J2000 = 1       //!< Transform from JD to J2000.
};

//! Define directions for aberration.
enum ae_aberration_direction_t {
  AE_ADD_ABERRATION = -1,  //!< Add the effect of aberration.
  AE_REMOVE_ABERRATION = 1 //!< Remove the effect of aberration.
};

//! Define whether a term uses T or d for physical ephemeride terms.
enum ae_physical_use_d_or_t {
  AE_PHYSICAL_END = 0,   //!< For marking the end of the series.
  AE_PHYSICAL_D,         //!< Use d = # days since J2000.
  AE_PHYSICAL_T          //!< Use T = # Julian centuries since J2000.
};

#define AE_PLANTBL_N_HARMONIC     18  //!< Number of arguments for ae_plantbl_t.
//! For holding data on planetary ephemerides.
struct ae_plantbl_t {
  char maxargs;     //!< The length of arg_tbl, lon_tbl, lat_tbl and rad_tbl.
  char max_harmonic[AE_PLANTBL_N_HARMONIC]; //!< The maximum harmonic.
  char max_power_of_t;  //!< The maximum power in time.
  char *arg_tbl;        //!< Table of arguments (?).
  void *lon_tbl;        //!< Table of longitudinal thingies (?).
  void *lat_tbl;        //!< Table of latitudinal thingies (?).
  void *rad_tbl;        //!< Table of radial thingies (?).
  double distance;      //!< Distance.
  double timescale;     //!< Timescale.
  double trunclvl;      //!< ?
};

//! For holding information on an orbit.
struct ae_orbit_t {
  char name[32];  //!< Name of the object.
  double epoch;     //!< Epoch of orbital elements.
  double i;	        //!< Inclination, in degrees.
  double W;	        //!< Longitude of the ascending node, in degrees.
  double w;	        //!< Argument of the perihelion, in degrees.
  double a;	        //!< Mean distance (semimajor axis), in AU.
  double dm;	    //!< Daily motion, in AU per day.
  double ecc;	    //!< Eccentricity.
  double M;	        //!< Mean anomaly, in degrees.
  double equinox;	//!< Epoch of equinox and ecliptic.

  // The following used by perterbation formulas.
  struct ae_plantbl_t *ptable;  //!< The perturbation tables.
  double L;   	    //!< Computed mean longitude.
  double r;	        //!< Computed radius vector.
  double plat;	    //!< Perturbation in ecliptic latitude.
};

//! For holding information on a star/planet.
//! Note the items for a star are in different measurement units
//! in the ASCII file description.
struct ae_star_t {
  char name[32];  //!< Object name (31 chars). 
  double epoch;     //!< Epoch of coordinates.
  double ra;        //!< Right Ascension, degrees.
  double dec;       //!< Declination, degrees.
  double px;        //!< Parallax, seconds of arc.
  double mura;      //!< Proper motion in R.A., "/century.
  double mudec;     //!< Proper motion in Dec., "/century.
  double v;         //!< Radial velocity in units of km/s.
  double mag;       //!< Visual magnitude.
};

//! For holding information on a JPL ephemeris file.
//! When an ASCII or binary header is read using ae_jpl_init(), this structure 
//! is populated.
struct ae_jpl_handle_t {
  double start;             //!< JD at which the ephemeris begins.
  double end;               //!< JD at which the ephemeris ends.
  double step;              //!< Ephemeris step-size, in days.
  int n_const;              //!< Number of constants.
  double au;                //!< Definition of the AU.
  double em_ratio;          //!< Earth-moon ratio.
  int ipt[13][3];           //!< Index pointers to the Chebyshev coefficents.
  int ephemeris_version;    //!< Version number.
  char **const_name;     //!< Constant names.
  double *const_val;        //!< Constant values.

  // The following are computed within the code.
  FILE *fp;                 //!< The JPL file being read from.
  int is_bin;               //!< Is the file binary or ASCII?
  int n_coeff;              //!< Number of Chebyshev coefficients.
  int swop;                 //!< Do we need to swop bytes?
  int rec_size;             //!< Record size.
  int last_block_num;       //!< Last block number read.
  int max_cheby;            //!< Maximum Chebyshev coefficients for a body.
  long ascii_start;     //!< The beginning byte number of the data (ASCII only).
  long ascii_len_block;     //!< The length of a block, in bytes (ASCII only).
  double sun_pv[6];         //!< For storing barycentric sun position/velocity.
  double sun_pv_t;          //!< The time for which sun_pv was last calculated.
  double *block;            //!< A local cache of the most recent block read.
  double *pc;               //!< For caching position Chebyshev coefficients.
  double *vc;               //!< For caching velocity Chebyshev coefficients.
  double last_x;            //!< Value at which coefficients were evaluated.
};

//! The coefficients for sine/cosine terms in a body's rotation data.
struct ae_physical_term_t {
  //! One of #AE_PHYSICAL_D or #AE_PHYSICAL_T, or #AE_PHYSICAL_END to mark last
  //! term.
  int time_var;         
  double a;             //!< The a in \f$a\cos(b + ct + dt^2)\f$.
  double b;             //!< The b in \f$a\cos(b + ct + dt^2)\f$, in degrees.
  double c;             //!< The b in \f$a\cos(b + ct + dt^2)\f$, in degrees.
  double d;             //!< The c in \f$a\cos(b + ct + dt^2)\f$, in degrees.
};

//! For holding information on a body's physical and rotational properties.
//! The members \c pole_ra_sin, \c pole_dec_cos and \c w_sin point to arrays
//! which are of length \f$4 \times n\f$, where \f$n\f$ is \c n_pole_ra_t, etc.
//! The elements of the arrays are \f${a_1, b_1, c_1, t_1, ..., a_n, b_n, c_n,
//! t_n}\f$, so that the contribution to the rotation terms are as \f$ a 
//! \sin(b + cx)\f$, where \f$x\f$ is the number of Julian centuries since J2000
//! if \f$t\f$ is #AE_PHYSICAL_D or the number of days since J2000 if
//! \f$t\f$ is #AE_PHYSICAL_T, and the sine term is a cosine for the case of
//! \c pole_dec_cos.
struct ae_physical_t {
  double r_mean;                //!< Mean radius, km. Not required.
  double r_eq;                  //!< Equatorial radius, km.
  double r_pole;                //!< Polar radius, km.
  double rms_dev_spheroid;   //!< RMS deviation from spheroid, km. Not required.
  double max_elev;              //!< Maximum elevation, km. Not required.
  double max_depress;           //!< Maximum depression, km. Not required.
  double pole_ra;               //!< North pole ra, J2000.
  //! North pole ra coefficient for T = Julian centuries since J2000.
  double pole_ra_t;
  //! Sine series coefficents for the north pole ra.
  struct ae_physical_term_t *ra_sin_term;
  double pole_dec;              //!< North pole declination, J2000.
  //! North pole dec coefficient for T = Julian centuries since J2000.
  double pole_dec_t;
  //! Cosine series coefficents for the north pole dec.
  struct ae_physical_term_t *dec_cos_term;
  double w;                     //!< Prime meridean.
  //! Prime meridean coefficient for d = days since J2000.
  double w_d;
  //! Prime meridean coefficient for d^2.
  double w_d_sq;
  //! Sine series coefficents for the prime meridian.
  struct ae_physical_term_t *w_sin_term;
};

// Predefined orbits.
extern struct ae_orbit_t ae_orb_mercury;
extern struct ae_orbit_t ae_orb_venus;
extern struct ae_orbit_t ae_orb_earth;
extern struct ae_orbit_t ae_orb_mars;
extern struct ae_orbit_t ae_orb_jupiter;
extern struct ae_orbit_t ae_orb_saturn;
extern struct ae_orbit_t ae_orb_uranus;
extern struct ae_orbit_t ae_orb_neptune;
extern struct ae_orbit_t ae_orb_pluto;
extern struct ae_orbit_t *ae_orb_planet[];

// Predefined perturbation tables.
extern struct ae_plantbl_t ae_mer404;
extern struct ae_plantbl_t ae_ven404;
extern struct ae_plantbl_t ae_ear404;
extern struct ae_plantbl_t ae_mar404;
extern struct ae_plantbl_t ae_jup404;
extern struct ae_plantbl_t ae_sat404;
extern struct ae_plantbl_t ae_ura404;
extern struct ae_plantbl_t ae_nep404;
extern struct ae_plantbl_t ae_plu404;
extern struct ae_plantbl_t ae_mlat404;
extern struct ae_plantbl_t ae_mlr404;

// Predefined physical parameters.
extern struct ae_physical_t ae_phys_sun;
extern struct ae_physical_t ae_phys_mercury;
extern struct ae_physical_t ae_phys_venus;
extern struct ae_physical_t ae_phys_earth;
extern struct ae_physical_t ae_phys_mars;
extern struct ae_physical_t ae_phys_jupiter;
extern struct ae_physical_t ae_phys_saturn;
extern struct ae_ring_geometry_t ae_saturn_rings;
extern struct ae_physical_t ae_phys_uranus;
extern struct ae_physical_t ae_phys_neptune;
extern struct ae_physical_t ae_phys_pluto;
extern struct ae_physical_t *ae_phys_planet[];
extern struct ae_physical_term_t ae_phys_no_term[];

// Predefined names.
extern const char *ae_constel_name[];
extern const char *ae_ss_name[];

void ae_altaz_to_radec(double last, double glat, double alt, double az,
                       double *ra, double *dec);
void ae_annual_aberration(double v_earth[], double p[], int direction);
double ae_cal_to_jd(long year, int month, double day);
int ae_cat_to_constel_index(const char *cat_name, char *full_name, int len);
int ae_coord_to_constel_index(double ra, double dec, double epoch);
double ae_ctime_to_jd(double t);
double ae_ctime_to_last(double t, double delta_t, double tlong, double nutl, 
                        double eps);
double aes_ctime_to_last(double t, double tlong);
double ae_delta_t(double jd_ut1);
void ae_delta_q(double q0[], double q1[], double *dra, double *ddec);
double ae_disc_semiminor(double a, double c, double lat);
double ae_disc_solid_angle(double a, double b, double dist);
double aes_disc_solid_angle(double jd_ut1, const struct ae_orbit_t *o_orb,
                            const struct ae_orbit_t *q_orb,
                            const struct ae_physical_t *phys);
void ae_diurnal_aberration(double last, double tlat, double trho, int direction,
                           double *ra, double *dec);
void ae_diurnal_parallax(double last, double tlat, double trho, double dist,
                         double *ra, double *dec);
double ae_dut1(double jd_utc);
double ae_epsilon(double jd_tt);
void ae_fk4_to_fk5(struct ae_star_t *star);
void ae_geocentric(double jd_tt, double q[], double e[], double v_e[], 
                   double *ra, double *dec, double *dist);
void ae_geocentric_moon_from_orbit(double jd_tt, const struct ae_orbit_t *o_orb,
                                   double *ra, double *dec, double *dist);
void ae_geocentric_from_cat(double jd_tt, double e[], double v_e[],
                            const struct ae_star_t *star, double *ra,
                            double *dec);
int ae_geocentric_from_jpl(struct ae_jpl_handle_t *jh, double jd_tt, 
                           int obj_num, double *ra, double *dec, double *dist);
void ae_geocentric_from_orbit(double jd_tt, const struct ae_orbit_t *o_orb, 
                              const struct ae_orbit_t *q_orb, double *ra, 
                              double *dec, double *dist);
void ae_geocentric_sun_from_orbit(double jd_tt, const struct ae_orbit_t *o_orb,
                                   double *ra, double *dec, double *dist);
void ae_geocentric_lat(double glat, double height, double *tlat, double *trho);
int ae_read_orbit_from_cat(const char *path, const char *name, 
                           struct ae_orbit_t *orb);
int ae_read_star_from_cat(const char *path, const char *name, 
                          struct ae_star_t *star);
void ae_gmoon(double jd, double rect[], double pol[]);
double ae_gmst(double jd_ut1, double jd_tt);
double ae_gast(double jd_ut1, double jd_tt, double nutl, double eps);
double ae_lmst(double jd_ut1, double jd_tt, double tlong);
double ae_last(double jd_ut1, double jd_tt, double tlong, double nutl,
               double eps);
double aes_last(double jd_ut1, double tlong);
void ae_gplan(double jd_tt, struct ae_plantbl_t *plan, double pobj[]);
double ae_g1plan(double jd_tt, struct ae_plantbl_t *plan);
void ae_g2plan(double jd_tt, struct ae_plantbl_t *plan, double pobj[],
               double *lp_equinox);
void ae_g3plan(double jd_tt, struct ae_plantbl_t *plan, double pobj[], 
               int objnum);
void ae_jd_to_cal(double jd, int *year, int *month, double *day);
double ae_jpl_get_const_val(const struct ae_jpl_handle_t *jh,
                            const char *const_name);
int ae_jpl_init(const char *path, struct ae_jpl_handle_t *jh);
int ae_jpl_init_ascii(const char *path, const char *header_path,
                      int search_dates, struct ae_jpl_handle_t *jh);
int ae_jpl_init_bin(const char *path, struct ae_jpl_handle_t *jh);
int ae_jpl_get_coords(struct ae_jpl_handle_t *jh, double jd_tt, int obj_num,
                      double r[], double v[], char is_planetary);
void ae_jpl_close(struct ae_jpl_handle_t *jh);
void ae_kepler(double jd_tt, const struct ae_orbit_t *orb, double q[]);
double ae_light_t(double p[], double q[], int do_retardation);
int ae_light_t_jpl(struct ae_jpl_handle_t *jh, double jd_tt, double e[], 
                   int obj_num, double q[], double v_q[], char is_planetary);
void ae_light_t_orbit(double jd_tt, double e[], const struct ae_orbit_t *orb,
                      double q[]);
double ae_mjd(double jd);
double ae_mod_2pi(double x);
double ae_mod_360(double x);
double ae_mod_180(double x);
void ae_nutation_lon_ob(double jd_tt, double *nutl, double *nuto);
void ae_nutate(double nutl, double nuto, double eps, double p[], int direction);
void aes_nutate(double jd_tt, double p[], int direction);
void ae_phys_pole(double jd_ut1, double jd_tt, const struct ae_physical_t *phys,
                  double n[], double *w);
int ae_is_retrograde(const struct ae_physical_t *phys);
void ae_polar_to_rect(double ra, double dec, double radius, double rect[]);
void ae_precess(double jd_tt, double r[], int direction);
void ae_radec_to_altaz(double last, double glat, double ra, double dec,
                       double *alt, double *az);
void ae_radec_to_gal(double ra, double dec, double *l, double *b, int fk5);
void ae_gal_to_radec(double ra, double dec, double *l, double *b, int fk5);
void ae_rect_to_polar(const double rect[], double *ra, double *dec, 
                      double *radius);
void ae_rect_to_polar_with_jd(double pp[], double jd, char ofdate, 
                              double polar[]);
double ae_refrac_visible(double alt, void *param);
double ae_refrac_ulich(double alt, void *param);
void ae_relativity(double p[], double q[], double o[]);
const char *ae_retcode_string(int retcode);
double ae_tdb(double jd);
void ae_subobs_point(double j[], double n[], double w, double f, int retrograde,
                     double *lat, double *lon);
void aes_subobs_point(double jd_ut1, const struct ae_orbit_t *o_orb,
                      const struct ae_orbit_t *q_orb,
                      const struct ae_physical_t *phys,
                      double *lat, double *lon, double *dist);
double ae_flattening(const struct ae_physical_t *phys);
void ae_topocentric(double jd_tt, double jd_ut1, double tlat, double glat, 
                    double tlong, double trho, double (*refrac)(double, void *),
                    void *refrac_param, double dist, double *ra, double *dec);
void aes_topocentric(double jd_ut1, double glat, double tlong, double dist,
                     double *ra, double *dec);
void ae_v_orbit(double jd_tt, const struct ae_orbit_t *orb, double v[]);
double ae_zatan2(double x, double y);

// Useful macros.

//! Extract the degrees/hours from a decimal angle/time.
//! \param a The decimal angle.
//! \return The degrees/hours.
#define ae_get_d(a) ((double)((int)(a)))

//! Extract the minutes from a decimal angle/time.
//! \param a The decimal angle.
//! \return The minutes.
#define ae_get_m(a) ((double)fabs(((int)(((a) - ae_get_d(a)) * 60))))

//! Extract the seconds from a decimal angle/time.
//! \param a The decimal angle.
//! \return The seconds.
#define ae_get_s(a) (fabs(((a) - ae_get_d(a) - ((double)((int)(((a) - ae_get_d(a)) * 60))) / 60) * 3600))

//! Defines a format string for printing D:M:S.
//! Use in conjuction with #ae_dms_arg.
#define ae_dms_fmt "%dd %d' %.2f\""

//! Defines a format string for printing D:M:S with fixed spacing.
//! Use in conjuction with #ae_dms_arg.
#define ae_dms_sfmt "%4dd %02d' %05.2f\""

//! Defines the arguments for printing D:M:S.
//! Use in conjunction with #ae_dms_fmt.
#define ae_dms_arg(a) (int)ae_get_d(a), (int)ae_get_m(a), ae_get_s(a)

//! Defines a format string for printing H:M:S.
//! Use in conjuction with #ae_hms_arg.
#define ae_hms_fmt "%dh %dm %.2fs"

//! Defines a format string for printing H:M:S with fixed spacing.
//! Use in conjuction with #ae_hms_arg.
#define ae_hms_sfmt "%3dh %02dm %05.2fs"

//! Defines the arguments for printing D:M:S.
//! Use in conjunction with #ae_dms_fmt.
#define ae_hms_arg(a) (int)ae_get_d(a / 15.0), (int)ae_get_m(a / 15.0), ae_get_s(a / 15.0)

// DEVEL START
//! For storing information on a planet's (i.e., Saturn's) rings.
struct ae_ring_geometry_t {
  int n;              //!< Number of rings.
  const char **name;        //!< Names of rings.
  const double *tau;        //!< Opacities of rings.
  //! Inner radius of rings, in Saturn equatorial radii.
  const double *a_inner;
  //! Outer radius of rings, in Saturn equatorial radii.
  const double *a_outer;
};

//! For holding Stokes parameters.
struct ae_stokes_t {
  double i;     //!< The intensity Stokes parameter.
  double q;     //!< The Q stokes parameter.
  double u;     //!< The U stokes parameter.
  double v;     //!< The V stokes parameter.
};

//! For modelling Saturn's emission.
struct ae_saturn_temp_model_t {
  double t_disc;                //!< The temperature of the disc.
  double *t_ring;               //!< The temperature of the rings.
  double t_background;          //!< The background temperature.
  //! Fractional linear polarisation of rings, assumed normal to ring surface.
  double *pol_ring;
  //!< Fractional linear polarisation of disc.
  double pol_disc;              //!< Fractional linear polarisation of disc.
  //! Angle east of north of disc polarisation, in degrees.
  double pol_disc_angle;
};

//! For storing a planetary map.
//! Planet maps are always square.
struct ae_planet_map_t {
  int n_pix;        //!< The dimension of each side of the map.
  //! The resolution of the map, in units of the planet's semi-major axis.
  double res;
  struct ae_stokes_t **data;    //!< The bitmap.
};

void ae_map_saturn(double res, const struct ae_physical_t *b,
                   const struct ae_ring_geometry_t *r, double inclination,
                   const struct ae_saturn_temp_model_t *t_sat,
                   struct ae_planet_map_t *map);
int ae_map_to_gnuplot(const struct ae_planet_map_t *map, char stokes,
                      const char *path);
void ae_map_t_eff(const struct ae_planet_map_t *map,
                  const struct ae_physical_t *b, double subobs_lat,
                  double background, struct ae_stokes_t *t_eff);
void ae_set_stokes(struct ae_stokes_t *s, double i, double q, double u,
                   double v);
void ae_copy_stokes(const struct ae_stokes_t *origin, struct ae_stokes_t *dest);
// DEVEL END

#ifdef __cplusplus
  }
#endif

#endif
