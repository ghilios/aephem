dnl Need the math library!
AC_CHECK_LIB([m], [cos], [],
             AC_MSG_ERROR([Aephem requires the math library!]))

dnl Check for pthread library: cfitsio might need it.
AC_CHECK_LIB([pthread],[pthread_create], [pthread_lib="-lpthread"],
             [pthread_lib=""])

dnl Checks for cfitsio library.
AC_ARG_WITH([cfitsio], [AS_HELP_STRING([--with-cfitsio=PREFIX],
            [include support for FITS headers @<:@default=autodetect@:>@])],
            [
             case "${withval}" in
               no) cfitsio_prefix= ;;
               *) cfitsio_prefix="${withval}" ;;
             esac
            ], [with_cfitsio=])

saved_ldflags=$LDFLAGS
saved_libs=$LIBS
if test "x$cfitsio_prefix" != "x"; then
  LDFLAGS="$LDFLAGS -L$cfitsio_prefix/lib -lm ${pthread_lib}"
else
  LDFLAGS="$LDFLAGS -lm ${pthread_lib}"
fi
AC_CHECK_LIB([cfitsio],[fits_init_cfitsio], [have_cfitsio_lib=yes],
             [have_cfitsio_lib=no])
LDFLAGS=$saved_ldflags
LIBS=$saved_libs

saved_cppflags=$CPPFLAGS
if test "x$cfitsio_prefix" != "x"; then
  CPPFLAGS="$CPPFLAGS -I$cfitsio_prefix/include"
fi
AC_CHECK_HEADERS([fitsio.h], [have_cfitsio_head=yes], [have_cfitsio_head=no])
CPPFLAGS=$saved_cppflags

if test "x$have_cfitsio_lib" = "xyes" -a "x$have_cfitsio_head" = "xyes"; then
  if test "x$cfitsio_prefix" = "x"; then
    CFITSIO_LIB="-lcfitsio ${pthread_lib}"
  else
    CFITSIO_LIB="-L${cfitsio_prefix}/lib -lcfitsio ${pthread_lib}"
    CFITSIO_INC="-I${cfitsio_prefix}/include"
  fi
  AC_DEFINE([HAVE_CFITSIO], [], [The cfitsio library is present.])
else
  if test "x$cfitsio_prefix" != "x"; then
    AC_MSG_ERROR([cannot find cfitsio library at PREFIX=${cfitsio_prefix}.])
  fi
fi
AC_SUBST(CFITSIO_INC)
AC_SUBST(CFITSIO_LIB)

AC_CONFIG_FILES([Makefile src/Makefile demo/Makefile])
AC_OUTPUT


#define AE_WCS_N_PARAM  21     //!< Length of \p param in #ae_wcs_t.

//! Define map projection types.
enum ae_proj_type_t {
  AE_PROJ_AZP = 0,      //!< Zenithal perspective (NOT IMPLEMENTED).
  AE_PROJ_SZP,          //!< Slant zenithal perspective (NOT IMPLEMENTED).
  AE_PROJ_TAN,          //!< Gnomonic.
  AE_PROJ_STG,          //!< Stereographic (NOT IMPLEMENTED).
  AE_PROJ_SIN,          //!< Slant orthographic (NOT IMPLEMENTED).
  AE_PROJ_ARC,          //!< Zenithal equidistant (NOT IMPLEMENTED).
  AE_PROJ_ZPN,          //!< Zenithal polynomial (NOT IMPLEMENTED).
  AE_PROJ_ZEA,          //!< Zenithal equal-area (NOT IMPLEMENTED).
  AE_PROJ_AIR,          //!< Airy projection (NOT IMPLEMENTED).
  AE_PROJ_CYP,          //!< Cylindrical perspective (NOT IMPLEMENTED).
  AE_PROJ_CEA,          //!< Cylindrical equal-area (NOT IMPLEMENTED).
  AE_PROJ_CAR,          //!< Plate carrée (NOT IMPLEMENTED).
  AE_PROJ_MER,          //!< Mercator (NOT IMPLEMENTED).
  AE_PROJ_SFL,          //!< Sanson-Flamsteed (NOT IMPLEMENTED).
  AE_PROJ_PAR,          //!< Parabolic (NOT IMPLEMENTED).
  AE_PROJ_MOL,          //!< Mollweide's (NOT IMPLEMENTED).
  AE_PROJ_AIT,          //!< Hammer-Aitoff (NOT IMPLEMENTED).
  AE_PROJ_COP,          //!< Conic perspective (NOT IMPLEMENTED).
  AE_PROJ_COE,          //!< Conic equal area (NOT IMPLEMENTED).
  AE_PROJ_COD,          //!< Conic equidistant (NOT IMPLEMENTED).
  AE_PROJ_COO,          //!< Conic orthomorphic (NOT IMPLEMENTED).
  AE_PROJ_BON,          //!< Bonne's equal area (NOT IMPLEMENTED).
  AE_PROJ_PCO,          //!< Polyconic (NOT IMPLEMENTED).
  AE_PROJ_TSC,          //!< Tangential spherical cube (NOT IMPLEMENTED).
  AE_PROJ_CSC,      //!< COBE quadrilaterlised spherical cube (NOT IMPLEMENTED).
  AE_PROJ_QSC           //!< Quadrilateralised spherical cube (NOT IMPLEMENTED).
};


//! For storing information on a World Coordinate System (WCS).
//! Astronomical maps are stored as pixels which have often been projected from 
//! a spherical surface. This structure can hold the parameters describing the
//! projection.
//!
//! See Calabretta, R. and Greisen, E.W., `Representations of celestial 
//! coordinates in FITS (Paper II)', Astronomy & Astrophysics, 395, 1077-1122, 
//! 2002. 
struct ae_wcs_t {
  int proj_type;        //!< One of ae_proj_type_t.
  //! The celestial longitude of the fiducial point, in degrees.
  double cel_f_lon;
  //! The celestial latitude of the fiducial point, in degrees.
  double cel_f_lat;
  //! The celestial longitude of the native pole, in degrees.
  double cel_p_lon;
  //! The celestial latitude of the native pole, in degrees.
  double cel_p_lat;
  double nat_f_lon;     //!< Native longitude of the fiducial point, in degrees.
  double nat_f_lat;     //!< Native latitude of the fiducial point, in degrees.
  double nat_p_lon;     //!< Native longitude of the celestial pole, in degrees.
  double nat_p_lat;     //!< Native latitude of the celestial pole, in degrees.
  double pix_f_x;       //!< The x-pixel number of the fiducial point.
  double pix_f_y;       //!< The y-pixel number of the fiducial point.
  //! Matrix for transforming pixels (scaling, rotating, etc.) to native
  //! spherical coordinates.degrees.
  double trans[2][2];
  //! The inverse of \p trans_mat.
  double trans_inv[2][2];
  //! Additional parameters used for some projections.
  double param[AE_WCS_N_PARAM];
};

int ae_wcs_compute_native_pole(struct ae_wcs_t *wcs);
void ae_wcs_dump(const struct ae_wcs_t *wcs);
int ae_wcs_invert_matrix(struct ae_wcs_t *wcs);
void ae_wcs_deproject(const struct ae_wcs_t *wcs, double x, double y,
                      double *lon, double *lat);
int ae_wcs_parse_fits_header(void *fp, struct ae_wcs_t *wcs,
                             char suffix, char *bad_keyword);
int ae_wcs_write_fits_header(void *fits_fp, const struct ae_wcs_t *wcs,
                             char suffix, char *lon_keyword, char *lat_keyword);
void ae_wcs_project(const struct ae_wcs_t *wcs, double lon, double lat,
                    double *x, double *y);

//------------------------------------------------------------------------------
//! Initialise an #ae_wcs_t structure.
//! Not all of the elements of #ae_wcs_t are always relevant. This routine fills
//! them with default values, chiefly zero. Or, e.g., \p trans is set to 
//! the identity matrix.
//!
//! \param wcs An #ae_wcs_t struct.
//------------------------------------------------------------------------------

void ae_wcs_init(struct ae_wcs_t *wcs) {
  wcs->cel_f_lon = 0;
  wcs->cel_f_lat = 0;
  wcs->nat_f_lon = 0;
  wcs->nat_f_lat = 0;
  wcs->pix_f_x = 0;
  wcs->pix_f_y = 0;
  wcs->nat_p_lon = 180.0;    // This is normally the case for ra/dec.
  wcs->nat_p_lat = 0;        // Rarely used, I believe.
  wcs->trans[0][0] = 1.0;
  wcs->trans[0][1] = 0;
  wcs->trans[1][0] = 0;
  wcs->trans[1][1] = 1.0;
}


//------------------------------------------------------------------------------
//! Invert the transformation matrix of a WCS.
//! It is assumed that \p trans has been initialised in \p wcs.
//!
//! \param wcs The WCS containing the matrix to be inverted. If the inverse
//!            exists, then \p trans_inv will be set.
//!
//! \return 1 on success; 0 if the matrix is not invertible.
//------------------------------------------------------------------------------

int ae_wcs_invert_matrix(struct ae_wcs_t *wcs) {
  double a, b, c, d, denom;

  a = wcs->trans[0][0];
  b = wcs->trans[0][1];
  c = wcs->trans[1][0];
  d = wcs->trans[1][1];

  if ((denom = a * d - b * c) == 0)
    return 0;

  wcs->trans_inv[0][0] = d / denom;
  wcs->trans_inv[0][1] = -b / denom;
  wcs->trans_inv[1][0] = -c / denom;
  wcs->trans_inv[1][1] = a / denom;

  return 1;
}

//------------------------------------------------------------------------------
//! Calculate the celestial coordinates of the native pole for a WCS.
//! In general, this is a finickity calculation. See:  Calabretta, R. and
//! Greisen, E.W., `Representations of celestial coordinates in FITS (Paper 
//! II)', <i>Astronomy & Astrophysics</i>, <b>395</b>, 1077&ndash;1122, 2002.
//!
//! If the WCS is not properly defined, there may be no solution. In this case,
//! this function does not change wcs.
//!
//! \param wcs The WCS to work with. This function assumes that it has been
//!            populated with all the parameters <i>except</i> for
//!            \p cel_p_lon and \p cel_p_lat. If there is a solution, these 
//!            values are set by this function.
//!
//! \return 1 if there exists a solution, 0 if there is no solution.
//------------------------------------------------------------------------------

int ae_wcs_compute_native_pole(struct ae_wcs_t *wcs) {
  double lon, lat, a, l1, l2, lat1, lat2, sin_cel_f_lat, cos_cel_f_lat;
  double cos_nat_f_lat, sin_nat_f_lat, sin_nat_pf_lon, cos_nat_pf_lon;

  cos_cel_f_lat = cos(wcs->cel_f_lat * AE_DTR);
  sin_cel_f_lat = sin(wcs->cel_f_lat * AE_DTR);
  cos_nat_f_lat = cos(wcs->nat_f_lat * AE_DTR);
  sin_nat_f_lat = sin(wcs->nat_f_lat * AE_DTR);
  cos_nat_pf_lon = cos((wcs->nat_p_lon - wcs->nat_f_lon) * AE_DTR);
  sin_nat_pf_lon = sin((wcs->nat_p_lon - wcs->nat_f_lon) * AE_DTR);

  // If this is less than zero, then there is no solution.
  a = 1.0 - cos_nat_f_lat * cos_nat_f_lat * sin_nat_pf_lon * sin_nat_pf_lon;
  if (a < 0)
    return 0;

  // Calculate the two terms for the latitude and construct the two possible
  // solutions.
  l1 = atan2(sin_nat_f_lat, cos_nat_f_lat * cos_nat_pf_lon);
  l2 = acos(sin_cel_f_lat / sqrt(a));
  lat1 = l1 + l2;
  lat2 = l1 - l2;

  // Sort through the degeneracies in the latitude, which must lie in the
  // range [-90, +90].
  if (wcs->nat_f_lat == 0 && wcs->cel_f_lat == 0 &&
      fabs(wcs->nat_f_lon - wcs->nat_p_lon) == 90.0) {
    // A special case.
    lat = wcs->nat_p_lat;
  }
  else if (fabs(lat1) <= M_PI / 2 && fabs(lat2) <= M_PI / 2) {
    // Both solutions are correct. Take the one closest to the native latitude
    // of the celestial pole.
    if (fabs(lat1 - wcs->nat_p_lat) < fabs(lat2 - wcs->nat_p_lat))
      lat = lat1;
    else
      lat = lat2;
  }
  else if (fabs(lat1) <= M_PI / 2 && fabs(lat2) > M_PI / 2)
    lat = lat1;
  else if (fabs(lat1) > M_PI / 2 && fabs(lat2) <= M_PI / 2)
    lat = lat2;
  else
    return 0;  // No valid solution!

  // Figure out what the longitude is now.
  if (fabs(wcs->cel_f_lat) == 90.0)
    lon = wcs->cel_f_lon;
  else if (lat == M_PI / 2)
    lon = wcs->cel_f_lon + wcs->nat_p_lon - wcs->nat_f_lon - 180.0;
  else if (lat == -M_PI / 2)
    lon = wcs->cel_f_lon - wcs->nat_p_lon + wcs->nat_f_lon;
  else
    lon = wcs->cel_f_lon -
          atan2(sin_nat_pf_lon * cos_nat_f_lat / cos_cel_f_lat,
                (sin_nat_f_lat - sin(lat) * sin_cel_f_lat) /
                cos(lat) / cos_cel_f_lat) * AE_RTD;

  // Record the values.
  wcs->cel_p_lon = lon;
  wcs->cel_p_lat = lat * AE_RTD;

  return 1;
}


//------------------------------------------------------------------------------
//! Given a WCS, deproject a pixel back onto the sphere.
//! Reference: Calabretta, R. and Greisen, E.W., `Representations of celestial 
//! coordinates in FITS (Paper II)', <i>Astronomy & Astrophysics</i>,
//! <b>395</b>, 1077&ndash;1122, 2002.
//!
//! \param wcs The WCS to use.
//! \param x The x-coordinate of the pixel.
//! \param y The y-coordinate of the pixel.
//! \param lat For returning the latitude, in degrees.
//! \param lon For returning the longitude, in degrees.
//------------------------------------------------------------------------------

void ae_wcs_deproject(const struct ae_wcs_t *wcs, double x, double y,
                      double *lon, double *lat) {
  double xc, yc, xpp, ypp, nlon, nlat, cos_nlat, sin_nlat, cos_cplat;
  double sin_cplat, cos_nplon, sin_nplon, th_1, th_2, eta, th_a, r0, y0, c, gam;
  double sin_th_1, sin_th_2, sin_th_12;

  // Step 1: convert to projection plane coordinates.
  xc = (double)(x - wcs->pix_f_x);
  yc = (double)(y - wcs->pix_f_y);
  xpp = wcs->trans[0][0] * xc + wcs->trans[0][1] * yc;
  ypp = wcs->trans[1][0] * xc + wcs->trans[1][1] * yc;

  // Step 2: deproject to native spherical coordinates.
  switch (wcs->proj_type) {
    case AE_PROJ_TAN:   // Gnomonic.
      nlon = atan2(xpp, -ypp);
      nlat = atan(1.0 / sqrt(xpp * xpp + ypp * ypp) / AE_DTR);
      break;

    case AE_PROJ_COE:   // Conic equal-area.
      th_a = wcs->param[1] * AE_DTR;
      eta = wcs->param[2] * AE_DTR;
      th_1 = th_a - eta;
      th_2 = th_a + eta;
      sin_th_1 = sin(th_1);
      sin_th_2 = sin(th_2);
      sin_th_12 = sin((th_1 + th_2) / 2.0);
      gam = sin_th_1 + sin_th_2;
      c = gam / 2.0;
      y0 = 2.0 * AE_RTD * sqrt(1.0 + sin_th_1 * sin_th_2 - gam * sin_th_12) /
           gam;
      r0 = sqrt(xpp * xpp + (y0 - ypp) * (y0 - ypp));
      if (th_a < 0)
        r0 *= -1.0;
      else if (th_a == 0)
        r0 = 0;
      nlon = atan2(xpp / r0, (y0 - ypp) / r0) / c;
      nlat = asin(1.0 / gam + sin_th_1 * sin_th_2 / gam -
                  gam * r0 * r0 * AE_DTR * AE_DTR / 4.0);
      break;

    default:
      fprintf(stderr, "Projection type %d is not implemented!\n",
                      wcs->proj_type);
      exit(0);
  }

  // Step 3: rotate to celestial spherical coordinates.
  cos_nlat = cos(nlat);
  sin_nlat = sin(nlat);
  cos_cplat = cos(wcs->cel_p_lat * AE_DTR);
  sin_cplat = sin(wcs->cel_p_lat * AE_DTR);
  cos_nplon = cos(nlon - wcs->nat_p_lon * AE_DTR);
  sin_nplon = sin(nlon - wcs->nat_p_lon * AE_DTR);
  *lon = wcs->cel_p_lon + atan2(-cos_nlat * sin_nplon,
                                sin_nlat * cos_cplat -
                                cos_nlat * sin_cplat * cos_nplon) * AE_RTD;
  *lat = asin(sin_nlat * sin_cplat + cos_nlat * cos_cplat * cos_nplon) * AE_RTD;
}


//------------------------------------------------------------------------------
//! Given a WCS, project a point on the sphere into a pixel.
//! Reference: Calabretta, R. and Greisen, E.W., `Representations of celestial 
//! coordinates in FITS (Paper II)', <i>Astronomy & Astrophysics</i>,
//! <b>395</b>, 1077&ndash;1122, 2002.
//!
//! \param wcs The WCS to use.
//! \param lat The latitude, in degrees.
//! \param lon The longitude, in degrees.
//! \param x For returning the x-coordinate of the pixel.
//! \param y For returning the y-coordinate of the pixel.
//------------------------------------------------------------------------------

void ae_wcs_project(const struct ae_wcs_t *wcs, double lon, double lat,
                    double *x, double *y) {
  double xc, yc, xpp, ypp, a, r, nlon, nlat, cos_lat, sin_lat, cos_cplat;
  double sin_cplat, cos_cplon, sin_cplon, th_1, th_2, eta, th_a, r0, y0, c, gam;
  double sin_th, sin_th_1, sin_th_2, sin_th_12;

  // Step 1: rotate to native spherical coordinates.
  cos_lat = cos(lat * AE_DTR);
  sin_lat = sin(lat * AE_DTR);
  cos_cplat = cos(wcs->cel_p_lat * AE_DTR);
  sin_cplat = sin(wcs->cel_p_lat * AE_DTR);
  cos_cplon = cos((lon - wcs->cel_p_lon) * AE_DTR);
  sin_cplon = sin((lon - wcs->cel_p_lon) * AE_DTR);
  nlon = wcs->nat_p_lon * AE_DTR + atan2(-cos_lat * sin_cplon,
                                         sin_lat * cos_cplat -
                                         cos_lat * sin_cplat * cos_cplon);
  nlat = asin(sin_lat * sin_cplat + cos_lat * cos_cplat * cos_cplon);

  // Step 2: project to projection plane coordinates.
  switch (wcs->proj_type) {
    case AE_PROJ_TAN:   // Gnomonic.
      if ((a = tan(nlat)) == 0) {
        // Silly user! This point diverges.  Throw zeros at him.
        xpp = 0;
        ypp = 0;
      }
      else {
        r = AE_RTD / a;
        xpp = sin(nlon) * r;
        ypp = -cos(nlon) * r;
      }
      break;

    case AE_PROJ_COE:   // Conic equal-area.
      th_a = wcs->param[1] * AE_DTR;
      eta = wcs->param[2] * AE_DTR;
      th_1 = th_a - eta;
      th_2 = th_a + eta;
      sin_th_1 = sin(th_1);
      sin_th_2 = sin(th_2);
      sin_th_12 = sin((th_1 + th_2) / 2.0);
      sin_th = sin(nlat);
      gam = sin_th_1 + sin_th_2;
      c = gam / 2.0;
      y0 = 2.0 * AE_RTD * sqrt(1.0 + sin_th_1 * sin_th_2 - gam * sin_th_12) /
           gam;
      r0 = 2.0 * AE_RTD * sqrt(1.0 + sin_th_1 * sin_th_2 - gam * sin_th) / gam;
      xpp = r0 * sin(c * nlon);
      ypp = -r0 * cos(c * nlon) + y0;
      break;

    default:
      fprintf(stderr, "Projection type %d is not implemented!\n",
                      wcs->proj_type);
      exit(0);
  }

  // Step 3: convert to pixels.
  xc = wcs->trans_inv[0][0] * xpp + wcs->trans_inv[0][1] * ypp;
  yc = wcs->trans_inv[1][0] * xpp + wcs->trans_inv[1][1] * ypp;
  *x = xc + wcs->pix_f_x;
  *y = yc + wcs->pix_f_y;

  return;
}


#ifdef HAVE_CFITSIO
//------------------------------------------------------------------------------
//! A wrapper for reading FITS keywords and deciding what to do.
//!
//! \param fits_fp An open FITS file.
//! \param name The keyword to attempt to read.
//! \param suffix A suffix to append to \p name.
//! \param suffix A suffix to append to \p name.
//! \param value For returning the value.
//! \param required If non-zero, then return an error if the keyword is not 
//!                 found.
//! \param default_value The value to assign to \p value if \p required is
//!                      false and the keyword is not found.
//! \param bad_keyword If not NULL, then return the offending keyword here if it
//!                    was not found and should have been, or was corrupt.
//!
//! \return return 0 on success; otherwise, or one of #ae_retcode_t.
//------------------------------------------------------------------------------

int ae_read_fits_keyword(fitsfile *fits_fp, const char *name, char suffix,
                         double *value, int required, double default_value,
                         char *bad_keyword) {
  int status;
  char field_name[FLEN_KEYWORD + 1], comment[FLEN_COMMENT + 1];
  char value_str[FLEN_VALUE + 1];

  status = 0;
  if (suffix)
    snprintf(field_name, FLEN_KEYWORD, "%s%c", name, suffix);
  else
    strcpy(field_name, name);
  fits_read_keyword(fits_fp, field_name, value_str, comment, &status);

  switch (status) {
    case 0:
      // Successful call.
      *value = atof(value_str);
      return 0;

    case KEY_NO_EXIST:
      if (required) {
        if (bad_keyword != NULL)
          strcpy(bad_keyword, field_name);
        return AE_RET_CFITSIO_NO_KEYWORD;
      }
      else {
        *value = default_value;
        return 0;
      }

    case VALUE_UNDEFINED:
      if (bad_keyword != NULL)
        strcpy(bad_keyword, field_name);
      return AE_RET_CFITSIO_BAD_KEYWORD;

    default:
      return AE_RET_CFITSIO_ERROR;
  }

  // We should never reach here, but just in case . . .
  return AE_RET_CFITSIO_ERROR;
}


//------------------------------------------------------------------------------
//! Make a FITS field name.
//!
//! \param template The template for the name. In the template, a (lowercase)
//!                 'i' is replaced with the value passed in \p i, 'j' with \p j
//!                 and 'a' with \p suffix, unless \p suffix is zero, in which
//!                 case 'a' is ignored.
//! \param i The axis number.  Normally 1 == longitude and 2 == latitude.
//! \param j The dimension number.  Normally 1 == x and 2 == y.
//! \param suffix The suffix. Pass 0 for no suffix.
//! \param result For returning the field name. It should be long enough to
//!               accept a string of length FLEN_KEYWORD (as defined in
//!               fitsio.h).
//------------------------------------------------------------------------------


void ae_fits_keyname(const char *template, int i, int j, char suffix,
                          char *result) {
  int k, l, len;

  len = strlen(template);
  for (k = 0, l = 0; k < len && l < FLEN_KEYWORD; k++) {
    switch (template[k]) {
      case 'i':
        l += snprintf(result + l, FLEN_KEYWORD - l, "%d", i);
        break;
      case 'j':
        l += snprintf(result + l, FLEN_KEYWORD - l, "%d", j);
        break;
      case 'a':
        if (suffix)
          result[l++] = suffix;
        break;
      default:
        result[l++] = template[k];
        break;
    }
  }
  result[l] = '\0';

  return;
}


//------------------------------------------------------------------------------
//! Wrapper for writing a header item to a FITS file.
//!
//! \param fp An open FITS file.
//! \param template The template for the name. In the template, a (lowercase)
//!                 'i' is replaced with the value passed in \p i, 'j' with \p j
//!                 and 'a' with \p suffix, unless \p suffix is zero, in which
//!                 case 'a' is ignored.
//! \param i The axis number.  Normally 1 == longitude and 2 == latitude.
//! \param j The dimension number.  Normally 1 == x and 2 == y.
//! \param suffix The suffix. Pass 0 for no suffix.
//! \param value The value to write.
//------------------------------------------------------------------------------


int ae_write_fits_key(fitsfile *fp, const char *template, int i, int j,
                        char suffix, double value) {
  char field_name[FLEN_KEYWORD];
  int f_stat;

  f_stat = 0;
  ae_fits_keyname(template, i, j, suffix, field_name);

  return fits_write_key(fp, TDOUBLE, field_name, &value, NULL, &f_stat);
}
#endif


//------------------------------------------------------------------------------
//! Parse a FITS header for the WCS information.
//! Reads in a FITS header and parses it for the WCS information. If the
//! necessary WCS information is found, this function then proceeds to call
//! ae_wcs_compute_native_pole() and ae_wcs_invert_matrix().
//!
//! This function requires that the following keywords be present in the FITS 
//! file, and will fail if they are not:
//! - \p CTYPE1a: this is used purely to determine the projection type.
//! - \p CRPIXja: the reference pixel coordinates,
//! - \p CDELTia: coordinate scale, assumed to be in degrees,
//! - \p CRVALia: celestial longitude and latitude of fiducial point, assumed to
//!      be in degees,
//!
//! where <tt>j = 1,2</tt> for x- and y-axes, respectively, <tt>i = 1,2</tt>
//! for longitude and latitude, respectively, and \p a is any suffix passed by
//! \p keyword_suffix.
//!
//! The following keywords are searched for, and if not found, then default
//! values are used.
//! - \p CDi_ja or \p PCi_ja: the linear transformation matrix. Default 1 for 
//!      diagonal terms and 0 for off-diagonal terms.
//! - \p PV1_ia, \p PVi_2a: native longitude and latitude of fiducial point,
//!      where \p i corresponds to the longitude axis. Default:
//!      - for longitude, 0 degrees;
//!      - for latitude:
//!        - 90 degrees if \p wcs->proj_type is zenithal (i.e., one of 
//!          #AE_PROJ_AZP, #AE_PROJ_SZP, #AE_PROJ_TAN, #AE_PROJ_STG,
//!          #AE_PROJ_SIN, #AE_PROJ_ARC, #AE_PROJ_ZPN, #AE_PROJ_ZEA
//!          or #AE_PROJ_AIR)),
//!        - equal to \f$\theta_a\f$ (see below) if \p wcs->proj_type is conical
//!          (i.e., one of #AE_PROJ_COP, #AE_PROJ_COE, #AE_PROJ_COD or
//!          #AE_PROJ_COO),
//!        - 0 degrees in all other cases.
//! - \p LONGPOLEa (<tt>= PVi_3a</tt>), LATPOLEa (<tt>= PVi_4a</tt>): native
//!      longitude and latitude of the pole, where \p i corresponds to the
//!      longitude axis. Default:
//!      - for longitude:
//!        - 0 degrees if longitude of fiducial point in celestial coordinates
//!          is greater or equal to longitude of fiducial point in native 
//!          coordinates,
//!        - 180 degrees otherwise,
//!      - for latitude, 90.
//!
//! Some projections need additional parameters. These can be stored in
//! PVi_0a&ndash;PVi_20a, where \p i corresponds to the latitude axis, though
//! usually only PVi_1a and PVi_2a are required. This function does <i>not</i>
//! require that these be present in the header, and assigns default values of
//! 0. However, the projection will fail if required parameters are missing. 
//! Here are the projection types that require these variables:
//! - #AE_PROJ_AZP: \f$\mu\f$ = \p PVi_1a, \f$\gamma\f$ = \p PVi_2a;
//! - #AE_PROJ_SZP: \f$\mu\f$ = \p PVi_1a, \f$\phi_c\f$ = \p PVi_2a,
//!                 \f$\theta_c\f$ = \p PVi_3a;
//! - #AE_PROJ_SIN: \f$\xi\f$ = \p PVi_1a, \f$\eta\f$ = \p PVi_2a,
//! - #AE_PROJ_ZPN: \f$P_0\f$ = \p PVi_0a, . . ., \f$P_{20}\f$ = \p PVi_20a;
//! - #AE_PROJ_AIR: \f$\theta_b\f$ = \p PVi_1a;
//! - #AE_PROJ_CYP: \f$\mu\f$ = \p PVi_1a, \f$\lambda\f$ = \p PVi_2a;
//! - #AE_PROJ_CEA: \f$\lambda\f$ = \p PVi_1a;
//! - #AE_PROJ_COP, #AE_PROJ_COE, #AE_PROJ_COD, #AE_PROJ_COO: 
//!                 \f$\lambda\f$ = \p PVi_1a, \f$\eta\f$ = \p PVi_2a;
//! - #AE_PROJ_BON: \f$\theta_1\f$ = \p PVi_1a.
//!
//! \b Important: this function assumes that the first axis (CTYPE1) is 
//! longitude and the second is latitude (CTYPE2).
//!
//! This function requires that cfitsio was enabled when aephem was configured.
//!
//! \param fits_fp An open FITS file. It is cast as void in the prototype in
//!                case the user's build does not include the cfitsio library.
//! \param wcs For recording the WCS.
//! \param suffix A suffix (or 0 for none) to append to all the keywords in the 
//!               FITS header&mdash;the FITS specification allows more than one
//!               WCS definition per header, distinguished by a single trailing
//!               letter in all the keywords.
//! \param bad_keyword If not NULL, then this is for returning the name of a
//!                    keyword that was expected but could not be read; will be
//!                    set if either #AE_RET_CFITSIO_BAD_KEYWORD or
//!                    #AE_RET_CFITSIO_NO_KEYWORD is returned. It should be able
//!                    to hold FLEN_KEYWORD characters (as defined by fitsio.h).
//!
//! \return An error code from #ae_retcode_t.
//------------------------------------------------------------------------------

int ae_wcs_parse_fits_header(void *fits_fp, struct ae_wcs_t *wcs,
                             char suffix, char *bad_keyword) {
#ifdef HAVE_CFITSIO
  char field_name[FLEN_KEYWORD + 1], comment[FLEN_COMMENT + 1];
  char value[FLEN_VALUE + 1];
  int i, j, status;
  double res_x, res_y;
  fitsfile *fp;

  fp = (fitsfile *)fits_fp;
  status = 0;

  // Figure out what kind of projection this is.
  if (suffix)
    snprintf(field_name, FLEN_KEYWORD, "CTYPE1%c", suffix);
  else
    strcpy(field_name, "CTYPE1");
  fits_read_keyword(fp, field_name, value, comment, &status);
  if (status == KEY_NO_EXIST) {
    if (bad_keyword != NULL)
      strcpy(bad_keyword, field_name);
    return AE_RET_CFITSIO_NO_KEYWORD;
  }
  else if (status == VALUE_UNDEFINED) {
    if (bad_keyword != NULL)
      strcpy(bad_keyword, field_name);
    return AE_RET_CFITSIO_BAD_KEYWORD;
  }
  else if (status)
    return AE_RET_CFITSIO_ERROR;

  // The projection type is encoded in the last three characters of the
  // 8-character string, surrounded by single quotes.
  if (strlen(value) < 8)
    return AE_RET_WCS_UNKNOWN_PROJ;
  value[9] = '\0';  // Get rid of trailing single quote.
  if (!strcasecmp(value + 6, "AZP"))
    wcs->proj_type = AE_PROJ_AZP;
  else if (!strcasecmp(value + 6, "SZP"))
    wcs->proj_type = AE_PROJ_SZP;
  else if (!strcasecmp(value + 6, "TAN"))
    wcs->proj_type = AE_PROJ_TAN;
  else if (!strcasecmp(value + 6, "STG"))
    wcs->proj_type = AE_PROJ_STG;
  else if (!strcasecmp(value + 6, "SIN"))
    wcs->proj_type = AE_PROJ_SIN;
  else if (!strcasecmp(value + 6, "ARC"))
    wcs->proj_type = AE_PROJ_ARC;
  else if (!strcasecmp(value + 6, "ZPN"))
    wcs->proj_type = AE_PROJ_ZPN;
  else if (!strcasecmp(value + 6, "ZEA"))
    wcs->proj_type = AE_PROJ_ZEA;
  else if (!strcasecmp(value + 6, "AIR"))
    wcs->proj_type = AE_PROJ_AIR;
  else if (!strcasecmp(value + 6, "CYP"))
    wcs->proj_type = AE_PROJ_CYP;
  else if (!strcasecmp(value + 6, "CEA"))
    wcs->proj_type = AE_PROJ_CEA;
  else if (!strcasecmp(value + 6, "CAR"))
    wcs->proj_type = AE_PROJ_CAR;
  else if (!strcasecmp(value + 6, "MER"))
    wcs->proj_type = AE_PROJ_MER;
  else if (!strcasecmp(value + 6, "SFL"))
    wcs->proj_type = AE_PROJ_SFL;
  else if (!strcasecmp(value + 6, "PAR"))
    wcs->proj_type = AE_PROJ_PAR;
  else if (!strcasecmp(value + 6, "MOL"))
    wcs->proj_type = AE_PROJ_MOL;
  else if (!strcasecmp(value + 6, "AIT"))
    wcs->proj_type = AE_PROJ_AIT;
  else if (!strcasecmp(value + 6, "COP"))
    wcs->proj_type = AE_PROJ_COP;
  else if (!strcasecmp(value + 6, "COE"))
    wcs->proj_type = AE_PROJ_COE;
  else if (!strcasecmp(value + 6, "COD"))
    wcs->proj_type = AE_PROJ_COD;
  else if (!strcasecmp(value + 6, "COO"))
    wcs->proj_type = AE_PROJ_COO;
  else if (!strcasecmp(value + 6, "BON"))
    wcs->proj_type = AE_PROJ_BON;
  else if (!strcasecmp(value + 6, "PCO"))
    wcs->proj_type = AE_PROJ_PCO;
  else if (!strcasecmp(value + 6, "TSC"))
    wcs->proj_type = AE_PROJ_TSC;
  else if (!strcasecmp(value + 6, "CSC"))
    wcs->proj_type = AE_PROJ_CSC;
  else if (!strcasecmp(value + 6, "QSC"))
    wcs->proj_type = AE_PROJ_QSC;
  else
    return AE_RET_WCS_UNKNOWN_PROJ;

  // Now get the required values.
  if ((i = ae_read_fits_keyword(fp, "CRPIX1", suffix, &wcs->pix_f_x,
                                1, 0, bad_keyword)))
    return i;
  if ((i = ae_read_fits_keyword(fp, "CRPIX2", suffix, &wcs->pix_f_y,
                                1, 0, bad_keyword)))
    return i;
  if ((i = ae_read_fits_keyword(fp, "CDELT1", suffix, &res_x,
                                1, 0, bad_keyword)))
    return i;
  if ((i = ae_read_fits_keyword(fp, "CDELT2", suffix, &res_y,
                                1, 0, bad_keyword)))
    return i;
  if ((i = ae_read_fits_keyword(fp, "CRVAL1", suffix, &wcs->cel_f_lon,
                                1, 0, bad_keyword)))
    return i;
  if ((i = ae_read_fits_keyword(fp, "CRVAL2", suffix, &wcs->cel_f_lat,
                                1, 0, bad_keyword)))
    return i;

  // Determine the linear transformation matrix.
  if ((i = ae_read_fits_keyword(fp, "CD1_1", suffix, &wcs->trans[0][0],
                                0, 1.0, bad_keyword)))
    return i;
  if ((i = ae_read_fits_keyword(fp, "CD1_2", suffix, &wcs->trans[0][1],
                                0, 0.0, bad_keyword)))
    return i;
  if ((i = ae_read_fits_keyword(fp, "CD2_1", suffix, &wcs->trans[1][0],
                                0, 0.0, bad_keyword)))
    return i;
  if ((i = ae_read_fits_keyword(fp, "CD2_2", suffix, &wcs->trans[1][1],
                                0, 1.0, bad_keyword)))
    return i;

  // Multiply by the resolution.
  wcs->trans[0][0] *= res_x;
  wcs->trans[0][1] *= res_x;
  wcs->trans[1][0] *= res_y;
  wcs->trans[1][1] *= res_y;

  // Get the native longitudes and latitudes of the fiducial point and pole.
  if ((i = ae_read_fits_keyword(fp, "PV1_1", suffix, &wcs->nat_f_lon,
                                0, 0.0, bad_keyword)))
    return i;
  switch (wcs->proj_type) {
    case AE_PROJ_AZP: case AE_PROJ_SZP: case AE_PROJ_TAN: case AE_PROJ_STG:
    case AE_PROJ_SIN: case AE_PROJ_ARC: case AE_PROJ_ZPN: case AE_PROJ_ZEA:
    case AE_PROJ_AIR:
      if ((i = ae_read_fits_keyword(fp, "PV1_2", suffix, &wcs->nat_f_lat,
                                    0, 90.0, bad_keyword)))
        return i;
      break;

    default:
      if ((i = ae_read_fits_keyword(fp, "PV1_2", suffix, &wcs->nat_f_lat,
                                    0, 0.0, bad_keyword)))
        return i;
      break;
  }
  if (ae_mod_180(wcs->cel_f_lon) >= wcs->nat_f_lon) {
    if ((i = ae_read_fits_keyword(fp, "PV1_1", suffix, &wcs->nat_p_lon,
                                  0, 0.0, bad_keyword)))
      return i;
  }
  else {
    if ((i = ae_read_fits_keyword(fp, "PV1_1", suffix, &wcs->nat_p_lon,
                                  0, 180.0, bad_keyword)))
      return i;
  }
  if ((i = ae_read_fits_keyword(fp, "PV1_4", suffix, &wcs->nat_p_lat,
                                0, 90.0, bad_keyword)))
    return i;

  // Get up to 21 additional projection parameters.
  for (i = 0; i < AE_WCS_N_PARAM; i++) {
    if ((j = ae_read_fits_keyword(fp, "PV2_0", suffix, &wcs->param[i],
                                  0, 0, bad_keyword)))
    return j;
  }

  // Phew! We made it this far!
  if (!ae_wcs_compute_native_pole(wcs))
    return AE_RET_BAD_NATIVE_POLE;
  if (!ae_wcs_invert_matrix(wcs))
    return AE_RET_NO_INVERSE;

  return 0;
#else
  return AE_RET_CFITSIO_NOT_ENABLED;
#endif
}

//------------------------------------------------------------------------------
//! Write WCS parameters to a FITS header.
//! This function writes all the elements of a WCS structure appropriately to
//! an open FITS file.
//!
//! The linear transformation matrix is divided into the pixel scale
//! (\p CDELTia) matrix (\p CD/PCi_ja), unless it is the identity matrix, in
//! which case nothing is written.
//!
//! The elements \p param of \p wcs are only written when the projection type
//! requires them&mdash;see the list in the usage of ae_wcs_parse_fits_header().
//!
//! \param fits_fp An open FITS file. It is cast as void in the prototype in
//!                case the user's build does not include the cfitsio library.
//! \param wcs The WCS to use.
//! \param suffix A suffix (or 0 for none) to append to all the keywords in the 
//!               FITS header&mdash;the FITS specification allows more than one
//!               WCS definition per header, distinguished by a single trailing
//!               letter in all the keywords.
//! \param lon_keyword The name of the longitude axis (e.g., "RA"), maximum four
//!                    characters.
//! \param lat_keyword The name of the latitude axis (e.g., "DEC"), maximum four
//!                    characters.
//!
//! \return A status code from the cfitsio calls: 0 on success.
//------------------------------------------------------------------------------

int ae_wcs_write_fits_header(void *fits_fp, const struct ae_wcs_t *wcs,
                             char suffix, char *lon_keyword,
                             char *lat_keyword) {
#ifdef HAVE_CFITSIO
  char field_val[FLEN_VALUE + 1], field_name[FLEN_KEYWORD + 1], tmpstr[32];
  int i, j, k, ret;
  double x;
  fitsfile *fp;

  if (wcs->proj_type < 0 || wcs->proj_type >= AE_PROJ_QSC)
    return AE_RET_WCS_UNKNOWN_PROJ;

  ret = 0;
  fp = (fitsfile *)fits_fp;

  if (fits_write_comment(fp, "WCS parameters written by "
                             "ae_wcs_write_fits_header() routine of aephem.",
                             &ret))
    return ret;

  // Write the axis information. The format should be like RA---TAN.
  for (i = 0; i < strlen(lon_keyword) && i < 4; i++)
    field_val[i] = toupper(lon_keyword[i]);
  for (; i < 5; i++)
    field_val[i] = '-';
  strncpy(field_val + 5, ae_proj_type_name[wcs->proj_type], 3);
  field_val[8] = '\0';
  ae_fits_keyname("CTYPEia", 1, 0, suffix, field_name);
  if (fits_write_key(fp, TSTRING, field_name, field_val, NULL, &ret))
    return ret;

  for (i = 0; i < strlen(lat_keyword) && i < 4; i++)
    field_val[i] = toupper(lat_keyword[i]);
  for (; i < 5; i++)
    field_val[i] = '-';
  strncpy(field_val + 5, ae_proj_type_name[wcs->proj_type], 3);
  field_val[8] = '\0';
  ae_fits_keyname("CTYPEia", 2, 0, suffix, field_name);
  if (fits_write_key(fp, TSTRING, field_name, field_val, NULL, &ret))
    return ret;

  // Write all the other fields.
  if ((ret = ae_write_fits_key(fp, "CRPIXja", 0, 1, suffix, wcs->pix_f_x)))
    return ret;
  if ((ret = ae_write_fits_key(fp, "CRPIXja", 0, 2, suffix, wcs->pix_f_y)))
    return ret;
  if ((ret = ae_write_fits_key(fp, "CDELTia", 1, 0, suffix, wcs->trans[0][0])))
    return ret;
  if ((ret = ae_write_fits_key(fp, "CDELTia", 2, 0, suffix, wcs->trans[1][1])))
    return ret;
  if ((ret = ae_write_fits_key(fp, "PVi_1a", 1, 0, suffix, wcs->nat_f_lon)))
    return ret;
  if ((ret = ae_write_fits_key(fp, "PVi_1a", 1, 0, suffix, wcs->nat_f_lat)))
    return ret;
  if ((ret = ae_write_fits_key(fp, "CRVALia", 1, 0, suffix, wcs->cel_f_lon)))
    return ret;
  if ((ret = ae_write_fits_key(fp, "CRVALia", 2, 0, suffix, wcs->cel_f_lat)))
    return ret;
  if ((ret = ae_write_fits_key(fp, "LONPOLEa", 0, 0, suffix, wcs->nat_p_lon)))
    return ret;
  if ((ret = ae_write_fits_key(fp, "LATPOLEa", 0, 0, suffix, wcs->nat_p_lat)))
    return ret;

  // Only write the transformation matrix if there are off-axis elements.
  if (wcs->trans[0][1] || wcs->trans[1][0]) {
    if ((ret = ae_write_fits_key(fp, "CDi_ja", 1, 1, suffix, 1.0)))
      return ret;
    if ((ret = ae_write_fits_key(fp, "CDi_ja", 2, 2, suffix, 1.0)))
      return ret;
    x = wcs->trans[0][1] / wcs->trans[0][0];
    if ((ret = ae_write_fits_key(fp, "CDi_ja", 1, 2, suffix, x)))
      return ret;
    x = wcs->trans[1][0] / wcs->trans[1][1];
    if ((ret = ae_write_fits_key(fp, "CDi_ja", 2, 2, suffix, x)))
      return ret;
  }

  // Only write the extra parameters if they are required.
  i = 1;
  j = 0;
  switch (wcs->proj_type) {
    case AE_PROJ_AIR: case AE_PROJ_CEA: case AE_PROJ_BON:
      j = 1;    // One parameter.
      break;
    case AE_PROJ_AZP: case AE_PROJ_SIN: case AE_PROJ_CYP: case AE_PROJ_COP:
    case AE_PROJ_COE: case AE_PROJ_COD: case AE_PROJ_COO:
      j = 2;    // Two parameters.
      break;
    case AE_PROJ_SZP:
      j = 3;    // Three parameters.
      break;
    case AE_PROJ_ZPN:
      i = 0;
      j = 21;   // All parameters.
      break;
  }
  for (k = i; k < i + j; k++) {
    sprintf(tmpstr, "PV2_%d", k);
    if ((ret = ae_write_fits_key(fp, tmpstr, 0, 0, suffix, wcs->param[k])))
      return ret;
  }

  return 0;
#else
  return AE_RET_CFITSIO_NOT_ENABLED;
#endif
}

//------------------------------------------------------------------------------
//! Dump to \p stdout the contents of a WCS struct.
//!
//! \p wcs The WCS to be dumped.
//------------------------------------------------------------------------------

void ae_wcs_dump(const struct ae_wcs_t *wcs) {
  int i;

  if (wcs->proj_type < 0 || wcs->proj_type > AE_PROJ_QSC) {
    printf("Unknown projection id %d!\n", wcs->proj_type);
    return;
  }
  printf("Projection type:          %s\n", ae_proj_type_name[wcs->proj_type]);
  printf("\n");
  printf("Celestial fiducial point: %10.6f %10.6f\n", wcs->cel_f_lon,
                                                        wcs->cel_f_lat);
  printf("Native fiducial point:    %10.6f %10.6f\n", wcs->nat_f_lon,
                                                        wcs->nat_f_lat);
  printf("Pixel fiducial point:     %10.6f %10.6f\n", wcs->pix_f_x,
                                                        wcs->pix_f_y);
  printf("Celestial pole:           %10.6f %10.6f\n", wcs->cel_p_lon,
                                                        wcs->cel_p_lat);
  printf("Native pole:              %10.6f %10.6f\n", wcs->nat_p_lon,
                                                        wcs->nat_p_lat);
  printf("\n");
  printf("Transformation matrix:   (%10.6f %10.6f)\n", wcs->trans[0][0],
                                                         wcs->trans[0][1]);
  printf("                         (%10.6f %10.6f)\n", wcs->trans[1][0],
                                                         wcs->trans[1][1]);
  printf("Inverse matrix:          (%10.6f %10.6f)\n", wcs->trans_inv[0][0],
                                                         wcs->trans_inv[0][1]);
  printf("                         (%10.6f %10.6f)\n", wcs->trans_inv[1][0],
                                                         wcs->trans_inv[1][1]);
  printf("\n");
  printf("Parameters:               ");

  for (i = 0; i < AE_WCS_N_PARAM; i++) {
    if (i && !(i % 5))
      printf("\n                          ");
    printf("%10.6f ", wcs->param[i]);
  }
  printf("\n");

  return;
}
