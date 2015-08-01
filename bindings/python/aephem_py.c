//! \file aephem_py.c
//! This main file for making the python bindings.
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

#include "aephem_py.h"
#include "help.h"
#include "../src/config.h"

//==============================================================================
// CONVERSION AND OTHER MISCELLANEOUS FUNCTIONS
//==============================================================================

int aepy_list_to_3d(PyObject *obj, void *p) {
  int i;
  char err_string[128];
  struct aepy_3d_t *pp;

  pp = (struct aepy_3d_t *)p;

  if (!PyList_Check(obj)) {
    sprintf(err_string, "'%s' must be a list object", pp->keyword);
    PyErr_SetString(PyExc_TypeError, err_string);
    return 0;
  }
  if (PyList_Size(obj) != 3) {
    sprintf(err_string, "'%s' must have exactly 3 elements", pp->keyword);
    PyErr_SetString(PyExc_TypeError, err_string);
    return 0;
  }

  for (i = 0; i < 3; i++)
    pp->vec[i] = PyFloat_AsDouble(PyList_GetItem(obj, i));

  return 1;
}


PyObject *aepy_3d_to_list(double *v) {
  int i;
  PyObject *ret;

  ret = PyList_New(3);
  for (i = 0; i < 3; i++)
    PyList_SetItem(ret, i, PyFloat_FromDouble(v[i]));

  return ret;
}

PyObject *aepy_make_pytuple(const char *type, ...) {
  int i, n;
  va_list item;
  PyObject *ret;

  n = strlen(type);
  ret = PyTuple_New(n);
  va_start(item, type);
  
  for (i = 0; i < n; i++) {
    switch (type[i]) {
      case 'i':
        PyTuple_SetItem(ret, i, PyInt_FromLong((long)va_arg(item, int)));
        break;
      case 'l':
        PyTuple_SetItem(ret, i, PyInt_FromLong((long)va_arg(item, long)));
        break;
      case 'd': default:
        PyTuple_SetItem(ret, i, PyFloat_FromDouble(va_arg(item, double)));
        break;
    }
  }

  va_end(item);

  return ret;
}

struct aepy_orb_t *get_orb_obj(PyObject *obj) {
  if (!strcmp("aephem.orb", obj->ob_type->tp_name))
    return (struct aepy_orb_t *)obj;
  else
    return NULL;
}

struct aepy_star_t *get_star_obj(PyObject *obj) {
  if (!strcmp("aephem.star", obj->ob_type->tp_name))
    return (struct aepy_star_t *)obj;
  else
    return NULL;
}

struct aepy_jpl_t *get_jpl_obj(PyObject *obj, int *obj_num,
                               int *not_planetary) {
  if (!PyTuple_Check(obj))
    return NULL;
  if (PyTuple_Size(obj) != 2 && PyTuple_Size(obj) != 3)
    return NULL;

  if (strcmp("aephem.jpl", PyTuple_GetItem(obj, 0)->ob_type->tp_name))
    return NULL;
  if (!PyInt_Check(PyTuple_GetItem(obj, 1)))
    return NULL;
  if (PyTuple_Size(obj) == 3) {
    if (!PyInt_Check(PyTuple_GetItem(obj, 2)))
      return NULL;
  }
    
  if (PyTuple_Size(obj) == 3)
    *not_planetary = PyInt_AsLong(PyTuple_GetItem(obj, 2));
  else
    *not_planetary = 0;

  *obj_num = PyInt_AsLong(PyTuple_GetItem(obj, 1));

  return (struct aepy_jpl_t *)PyTuple_GetItem(obj, 0);
}

//==============================================================================
// BINDINGS FOR CALENDAR.C
//==============================================================================

static PyObject *aepy_ctime_to_jd(struct aepy_t *self, PyObject *args,
                                  PyObject *keys) {
  char *keywords[] = {"t", NULL};
  double t;

  if (!PyArg_ParseTupleAndKeywords(args, keys, "d", keywords, &t))
    return NULL;

  return PyFloat_FromDouble(ae_ctime_to_jd(t));
}

static PyObject *aepy_ctime_to_last(struct aepy_t *self, PyObject *args,
                                    PyObject *keys) {
  char *keywords[] = {"t", "tlon", "delta_t", "nutl", "eps", NULL};
  double jd_ut1, jd_tt, t, delta_t, tlon, nutl, nuto, eps;
  
  delta_t = 1e99;
  nutl = 1e99;
  eps = 1e99;
  if (!PyArg_ParseTupleAndKeywords(args, keys, "dd|ddd", keywords, &t, &tlon,
                                   &delta_t, &nutl, &eps))
    return NULL;

  // Get the Julian dates.
  jd_ut1 = ae_ctime_to_jd(t);
  jd_ut1 += ae_dut1(jd_ut1) * AE_D_PER_S;
  if (delta_t > 1e10)
    delta_t = ae_delta_t(jd_ut1);
  jd_tt = jd_ut1 + delta_t * AE_D_PER_S;

  // Get nutation parameters, if necessary.
  if (nutl > 1e10)
    ae_nutation_lon_ob(jd_tt, &nutl, &nuto);
  if (eps > 1e10)
    eps = ae_epsilon(jd_tt);

  return PyFloat_FromDouble(ae_last(jd_ut1, jd_tt, tlon, nutl, eps));
}

static PyObject *aepy_cal_to_jd(struct aepy_t *self, PyObject *args,
                                PyObject *keys) {
  char *keywords[] = {"month", "year", "day", NULL};
  int month;
  long year;
  double day;

  if (!PyArg_ParseTupleAndKeywords(args, keys, "lid", keywords, &year, &month,
                                   &day))
    return NULL;

  return PyFloat_FromDouble(ae_cal_to_jd(year, month, day));
}

static PyObject *aepy_jd_to_cal(struct aepy_t *self, PyObject *args,
                                PyObject *keys) {
  char *keywords[] = {"jd", NULL};
  int year, month;
  double jd, day;

  if (!PyArg_ParseTupleAndKeywords(args, keys, "d", keywords, &jd))
    return NULL;

  ae_jd_to_cal(jd, &year, &month, &day);

  return aepy_make_pytuple("iid", year, month, day);
}

static PyObject *aepy_gmst(struct aepy_t *self, PyObject *args,
                           PyObject *keys) {
  char *keywords[] = {"jd_ut1", "jd_tt", NULL};
  double jd_ut1, jd_tt;

  if (!PyArg_ParseTupleAndKeywords(args, keys, "dd", keywords, &jd_ut1, &jd_tt))
    return NULL;

  return PyFloat_FromDouble(ae_gmst(jd_ut1, jd_tt));
}

static PyObject *aepy_gast(struct aepy_t *self, PyObject *args,
                           PyObject *keys) {
  char *keywords[] = {"jd_ut1", "jd_tt", "nutl", "eps", NULL};
  double jd_ut1, jd_tt, nutl, eps;

  if (!PyArg_ParseTupleAndKeywords(args, keys, "dddd", keywords, &jd_ut1, 
                                   &jd_tt, &nutl, &eps))
    return NULL;

  return PyFloat_FromDouble(ae_gast(jd_ut1, jd_tt, nutl, eps));
}

static PyObject *aepy_lmst(struct aepy_t *self, PyObject *args,
                           PyObject *keys) {
  char *keywords[] = {"jd_ut1", "jd_tt", "tlon", NULL};
  double jd_ut1, jd_tt, tlon;

  if (!PyArg_ParseTupleAndKeywords(args, keys, "ddd", keywords, &jd_ut1,
                                   &jd_tt, &tlon))
    return NULL;

  return PyFloat_FromDouble(ae_lmst(jd_ut1, jd_tt, tlon));
}

static PyObject *aepy_last(struct aepy_t *self, PyObject *args,
                           PyObject *keys) {
  char *keywords[] = {"jd_ut1", "tlon", "jd_tt", "nutl", "eps", NULL};
  double jd_ut1, tlon, jd_tt, nutl, eps, nuto;

  jd_tt = 1e99;
  nutl = 1e99;
  eps = 1e99;
  if (!PyArg_ParseTupleAndKeywords(args, keys, "dd|ddd", keywords, &jd_ut1,
                                   &tlon, &jd_tt, &nutl, &eps))
    return NULL;
  
  // Get TT, if necessary.
  if (jd_tt > 1e10)
    jd_tt = jd_ut1 + ae_delta_t(jd_ut1) * AE_D_PER_S;

  // Get nutation parameters, if necessary.
  if (nutl > 1e10)
    ae_nutation_lon_ob(jd_tt, &nutl, &nuto);
  if (eps > 1e10)
    eps = ae_epsilon(jd_tt);

  return PyFloat_FromDouble(ae_last(jd_ut1, jd_tt, tlon, nutl, eps));
}

static PyObject *aepy_tdb(struct aepy_t *self, PyObject *args,
                          PyObject *keys) {
  char *keywords[] = {"jd", NULL};
  double jd;

  if (!PyArg_ParseTupleAndKeywords(args, keys, "d", keywords, &jd))
    return NULL;

  return PyFloat_FromDouble(ae_tdb(jd));
}

static PyObject *aepy_mjd(struct aepy_t *self, PyObject *args,
                          PyObject *keys) {
  char *keywords[] = {"jd", NULL};
  double jd;

  if (!PyArg_ParseTupleAndKeywords(args, keys, "d", keywords, &jd))
    return NULL;

  return PyFloat_FromDouble(ae_mjd(jd));
}


static PyObject *aepy_delta_t(struct aepy_t *self, PyObject *args,
                              PyObject *keys) {
  char *keywords[] = {"jd_ut1", NULL};
  double jd_ut1;

  if (!PyArg_ParseTupleAndKeywords(args, keys, "d", keywords, &jd_ut1))
    return NULL;

  return PyFloat_FromDouble(ae_delta_t(jd_ut1));
}

static PyObject *aepy_dut1(struct aepy_t *self, PyObject *args,
                           PyObject *keys) {
  char *keywords[] = {"jd_utc", NULL};
  double jd_utc;

  if (!PyArg_ParseTupleAndKeywords(args, keys, "d", keywords, &jd_utc))
    return NULL;

  return PyFloat_FromDouble(ae_dut1(jd_utc));
}




//==============================================================================
// BINDINGS FOR COORD.C
//==============================================================================

static PyObject *aepy_rect_to_polar(struct aepy_t *self, PyObject *args,
                                    PyObject *keys) {
  char *keywords[] = {"rect", "radius", NULL};
  int do_radius;
  double ra, dec, radius;
  struct aepy_3d_t rect;

  strcpy(rect.keyword, "rect");
  do_radius = 0;
  if (!PyArg_ParseTupleAndKeywords(args, keys, "O&|i", keywords,
                                   &aepy_list_to_3d, (void *)&rect, &do_radius))
    return NULL;

  if (do_radius) {
    ae_rect_to_polar(rect.vec, &ra, &dec, &radius);
    return aepy_make_pytuple("ddd", ra, dec, radius);
  }
  else {
    ae_rect_to_polar(rect.vec, &ra, &dec, NULL);
    return aepy_make_pytuple("dd", ra, dec);
  }
}

static PyObject *aepy_polar_to_rect(struct aepy_t *self, PyObject *args,
                                    PyObject *keys) {
  char *keywords[] = {"ra", "dec", "radius", NULL};
  double ra, dec, radius, rect[3];

  if (!PyArg_ParseTupleAndKeywords(args, keys, "ddd", keywords, &ra, &dec,
                                   &radius))
    return NULL;

  ae_polar_to_rect(ra, dec, radius, rect);

  return aepy_3d_to_list(rect);
}

static PyObject *aepy_radec_to_altaz(struct aepy_t *self, PyObject *args,
                                     PyObject *keys) {
  char *keywords[] = {"last", "glat", "ra", "dec", NULL};
  double last, glat, ra, dec, alt, az;

  if (!PyArg_ParseTupleAndKeywords(args, keys, "dddd", keywords, &last, &glat,
                                   &ra, &dec))
    return NULL;

  ae_radec_to_altaz(last, glat, ra, dec, &alt, &az);

  return aepy_make_pytuple("dd", alt, az);
}


static PyObject *aepy_altaz_to_radec(struct aepy_t *self, PyObject *args,
                                     PyObject *keys) {
  char *keywords[] = {"last", "glat", "alt", "az", NULL};
  double last, glat, ra, dec, alt, az;

  if (!PyArg_ParseTupleAndKeywords(args, keys, "dddd", keywords, &last, &glat,
                                   &alt, &az))
    return NULL;

  ae_altaz_to_radec(last, glat, alt, az, &ra, &dec);

  return aepy_make_pytuple("dd", ra, dec);
}

static PyObject *aepy_radec_to_gal(struct aepy_t *self, PyObject *args,
                                   PyObject *keys) {
  char *keywords[] = {"ra", "dec", "is_fk4", NULL};
  int is_fk4;
  double ra, dec, l, b;

  is_fk4 = 0;
  if (!PyArg_ParseTupleAndKeywords(args, keys, "dd|i", keywords, &ra, &dec,
                                   &is_fk4))
    return NULL;

  if (is_fk4)
    ae_radec_to_gal(ra, dec, &l, &b, 0);
  else
    ae_radec_to_gal(ra, dec, &l, &b, 1);

  printf("!!!!!!!!!!!!!!!!!\n");
  return aepy_make_pytuple("dd", l, b);
}

static PyObject *aepy_gal_to_radec(struct aepy_t *self, PyObject *args,
                                   PyObject *keys) {
  char *keywords[] = {"l", "b", "is_fk4", NULL};
  int is_fk4;
  double ra, dec, l, b;

  is_fk4 = 0;
  if (!PyArg_ParseTupleAndKeywords(args, keys, "dd|i", keywords, &l, &b,
                                   &is_fk4))
    return NULL;

  if (is_fk4)
    ae_gal_to_radec(l, b, &ra, &dec, 0);
  else
    ae_gal_to_radec(l, b, &ra, &dec, 1);

  printf("!!!!!!!!!!!!!!!!!\n");
  return aepy_make_pytuple("dd", ra, dec);
}

static PyObject *aepy_geocentric_lat(struct aepy_t *self, PyObject *args,
                                     PyObject *keys) {
  char *keywords[] = {"glat", "height", NULL};
  double glat, height, tlat, trho;

  if (!PyArg_ParseTupleAndKeywords(args, keys, "dd", keywords, &glat, &height))
    return NULL;

  ae_geocentric_lat(glat, height, &tlat, &trho);

  return aepy_make_pytuple("dd", tlat, trho);
}

static PyObject *aepy_delta_q(struct aepy_t *self, PyObject *args,
                              PyObject *keys) {
  char *keywords[] = {"q0", "q1", NULL};
  double dra, ddec;
  struct aepy_3d_t q0, q1;

  strcpy(q0.keyword, "q0");
  strcpy(q1.keyword, "q1");
  if (!PyArg_ParseTupleAndKeywords(args, keys, "O&O&", keywords,
                                   &aepy_list_to_3d, (void *)&q0,
                                   &aepy_list_to_3d, (void *)&q1))
    return NULL;

  ae_delta_q(q0.vec, q1.vec, &dra, &ddec);

  return aepy_make_pytuple("dd", dra, ddec);
}

void aepy_theta_to_dhms(double theta, double *dh, double *m, double *s,
                        int do_h) {
  double sign;

  if (theta < 0)
    sign = -1.0;
  else
    sign = 1.0;
  if (do_h)
    theta /= 15.0;

  *dh = (double)((int)theta);
  theta = fabs(theta - (*dh)) * 60.0;
  *m = (double)((int)theta);
  *s = fabs(theta - (*m)) * 60.0;

  if (*dh == 0 && sign < 0) {
    if (*m == 0)
      *s *= -1;
    else
      *m *= -1;
  }
  
  return;
}

static PyObject *aepy_decimal(struct aepy_t *self, PyObject *args,
                              PyObject *keys) {
  char *keywords[] = {"d", "m", "s", "hour", NULL};
  int hour;
  double d, m, s, mult, theta;

  hour = 0;
  if (!PyArg_ParseTupleAndKeywords(args, keys, "ddd|i", keywords, &d, &m, &s,
                                   &hour))
    return NULL;

  if (d < 0 || m < 0 || s < 0)
    mult = -1.0;
  else
    mult = 1.0;
  if (hour)
    mult *= 15.0;

  theta = (floor(fabs(d)) + floor(fabs(m)) / 60.0 + fabs(s) / 3600.0) * mult;

  return PyFloat_FromDouble(theta);
}

static PyObject *aepy_dms(struct aepy_t *self, PyObject *args,
                          PyObject *keys) {
  char *keywords[] = {"angle", NULL};
  double theta, d, m, s;

  if (!PyArg_ParseTupleAndKeywords(args, keys, "d", keywords, &theta))
    return NULL;

  aepy_theta_to_dhms(theta, &d, &m, &s, 0);

  return aepy_make_pytuple("ddd", d, m, s);
}

static PyObject *aepy_hms(struct aepy_t *self, PyObject *args,
                          PyObject *keys) {
  char *keywords[] = {"angle", NULL};
  double theta, h, m, s;

  if (!PyArg_ParseTupleAndKeywords(args, keys, "d", keywords, &theta))
    return NULL;
  
  aepy_theta_to_dhms(theta, &h, &m, &s, 1);

  return aepy_make_pytuple("ddd", h, m, s);
}

static PyObject *aepy_dms_str(struct aepy_t *self, PyObject *args,
                              PyObject *keys) {
  char *keywords[] = {"angle", "fixed", "sign", NULL};
  char str[24];
  int fixed, sign, i;
  double theta, d, m, s;

  fixed = 0;
  sign = 0;
  if (!PyArg_ParseTupleAndKeywords(args, keys, "d|ii", keywords, &theta, &fixed,
                                   &sign))
    return NULL;

  aepy_theta_to_dhms(theta, &d, &m, &s, 0);

  if (d < 0 || m < 0 || s < 0) {
    str[0] = '-';
    i = 1;
  }
  else if (sign) {
    str[0] = '+';
    i = 1;
  }
  else
    i = 0;

  if (fixed)
    sprintf(str + i, "%3dd %02d' %05.2f\"", abs((int)d), abs((int)m), fabs(s));
  else
    sprintf(str + i, "%dd %d' %.2f\"", abs((int)d), abs((int)m), fabs(s));

  return PyString_FromString(str);
}

static PyObject *aepy_hms_str(struct aepy_t *self, PyObject *args,
                              PyObject *keys) {
  char *keywords[] = {"angle", "fixed", "sign", NULL};
  char str[24];
  int fixed, sign, i;
  double theta, h, m, s;

  fixed = 0;
  sign = 0;
  if (!PyArg_ParseTupleAndKeywords(args, keys, "d|ii", keywords, &theta, &fixed,
                                   &sign))
    return NULL;

  aepy_theta_to_dhms(theta, &h, &m, &s, 1);

  if (h < 0 || m < 0 || s < 0) {
    str[0] = '-';
    i = 1;
  }
  else if (sign) {
    str[0] = '+';
    i = 1;
  }
  else
    i = 0;

  if (fixed)
    sprintf(str + i, "%3dh %02dm %05.2fs", abs((int)h), abs((int)m), fabs(s));
  else
    sprintf(str + i, "%dh %dm %.2fs", abs((int)h), abs((int)m), fabs(s));

  return PyString_FromString(str);
}




//==============================================================================
// BINDINGS FOR PRECESS.C, NUTATE.C AND EPSILON.C.
//==============================================================================

static PyObject *aepy_precess(struct aepy_t *self, PyObject *args,
                              PyObject *keys) {
  char *keywords[] = {"r", "jd_tt", "to_j2000", NULL};
  int direction;
  double jd_tt;
  struct aepy_3d_t r;

  direction = AE_FROM_J2000;
  strcpy(r.keyword, "r");
  if (!PyArg_ParseTupleAndKeywords(args, keys, "O&d|i", keywords,
                                   &aepy_list_to_3d, &r, &jd_tt, &direction))
    return NULL;

  if (direction != AE_FROM_J2000) {
    if (direction == 0)
      direction = AE_FROM_J2000;
    else
      direction = AE_TO_J2000;
  }

  ae_precess(jd_tt, r.vec, direction);

  return aepy_3d_to_list(r.vec);
}

static PyObject *aepy_nutation_lon_ob(struct aepy_t *self, PyObject *args,
                                      PyObject *keys) {
  char *keywords[] = {"jd_tt", NULL};
  double jd_tt, nutl, nuto;

  if (!PyArg_ParseTupleAndKeywords(args, keys, "d", keywords, &jd_tt))
    return NULL;

  ae_nutation_lon_ob(jd_tt, &nutl, &nuto);

  return aepy_make_pytuple("dd", nutl, nuto);
}

static PyObject *aepy_nutate(struct aepy_t *self, PyObject *args,
                             PyObject *keys) {
  char *keywords[] = {"r", "nutl", "nuto", "eps", "jd_tt", "to_j2000", NULL};
  int direction, has_tt, no_tt;
  double jd_tt, nutl, nuto, eps;
  struct aepy_3d_t r;

  jd_tt = 1e99;
  nutl = 1e99;
  nuto = 1e99;
  eps = 1e99;
  direction = AE_FROM_J2000;
  strcpy(r.keyword, "r");
  if (!PyArg_ParseTupleAndKeywords(args, keys, "O&|ddddi", keywords,
                                   &aepy_list_to_3d, (void *)&r, &nutl, &nuto,
                                   &eps, &jd_tt, &direction))
    return NULL;

  if (direction != AE_FROM_J2000) {
    if (direction == 0)
      direction = AE_FROM_J2000;
    else
      direction = AE_TO_J2000;
  }

  has_tt = 0;
  no_tt = 0;
  if (jd_tt < 1e50)
    has_tt = 1;
  if (nutl < 1e50 && nuto < 1e50 && eps < 1e50)
    no_tt = 1;


  if (has_tt == no_tt) {
    PyErr_SetString(PyExc_KeyError, "Either all of nutl, nuto and eps must be "
                                    "passed, or jd_tt must be passed.");
    return NULL;
  }

  if (has_tt)
    aes_nutate(jd_tt, r.vec, direction);
  else
    ae_nutate(nutl, nuto, eps, r.vec, direction);

  return aepy_3d_to_list(r.vec);
}

static PyObject *aepy_epsilon(struct aepy_t *self, PyObject *args,
                              PyObject *keys) {
  char *keywords[] = {"jd_tt", NULL};
  double jd_tt;

  if (!PyArg_ParseTupleAndKeywords(args, keys, "d", keywords, &jd_tt))
    return NULL;

  return PyFloat_FromDouble(ae_epsilon(jd_tt));
}




//==============================================================================
// BINDINGS FOR CONSTEL.C
//==============================================================================

static PyObject *aepy_constellation(struct aepy_t *self, PyObject *args,
                                    PyObject *keys) {
  char *keywords[] = {"ra", "dec", "epoch", "abbriev", NULL};
  char ret[32];
  int index, abbriev;
  double ra, dec, epoch;

  epoch = AE_J2000;
  abbriev = 0;
  if (!PyArg_ParseTupleAndKeywords(args, keys, "dd|di", keywords, &ra, &dec,
                                   &epoch, &abbriev))
    return NULL;

  if ((index = ae_coord_to_constel_index(ra, dec, epoch)) < 0)
    strcpy(ret, "unknown");
  else if (abbriev) {
    strncpy(ret, ae_constel_name[index], 3);
    ret[3] = '\0';
  }
  else
    strcpy(ret, ae_constel_name[index] + 4);

  return PyString_FromString(ret);
}




//==============================================================================
// BINDINGS FOR ANNUAL_AB.C, DIURNAL_AB.C AND DIURNAL_PX.C.
//==============================================================================

static PyObject *aepy_annual_aberration(struct aepy_t *self, PyObject *args,
                                        PyObject *keys) {
  char *keywords[] = {"v", "p", "add", NULL};
  int add, direction;
  struct aepy_3d_t v, p;

  add = 0;
  strcpy(v.keyword, "v");
  strcpy(p.keyword, "p");
  if (!PyArg_ParseTupleAndKeywords(args, keys, "O&O&|i", keywords,
                                   &aepy_list_to_3d, (void *)&v,
                                   &aepy_list_to_3d, (void *)&p, &add))
    return NULL;

  if (add)
    direction = AE_ADD_ABERRATION;
  else
    direction = AE_REMOVE_ABERRATION;

  ae_annual_aberration(v.vec, p.vec, direction);

  return aepy_3d_to_list(p.vec);
}

static PyObject *aepy_diurnal_aberration(struct aepy_t *self, PyObject *args,
                                         PyObject *keys) {
  char *keywords[] = {"last", "tlat", "trho", "ra", "dec", "add", NULL};
  int add, direction;
  double ra, dec, last, tlat, trho;

  add = 0;
  if (!PyArg_ParseTupleAndKeywords(args, keys, "ddddd|i", keywords, &last,
                                   &tlat, &trho, &ra, &dec, &add))
    return NULL;

  if (add)
    direction = AE_ADD_ABERRATION;
  else
    direction = AE_REMOVE_ABERRATION;

  ae_diurnal_aberration(last, tlat, trho, direction, &ra, &dec);

  return aepy_make_pytuple("dd", ra, dec);
}

static PyObject *aepy_diurnal_parallax(struct aepy_t *self, PyObject *args,
                                       PyObject *keys) {
  char *keywords[] = {"last", "tlat", "trho", "dist", "ra", "dec", NULL};
  double ra, dec, last, tlat, trho, dist;

  if (!PyArg_ParseTupleAndKeywords(args, keys, "dddddd", keywords, &last, &tlat,
                                   &trho, &dist, &ra, &dec))
    return NULL;

  ae_diurnal_parallax(last, tlat, trho, dist, &ra, &dec);

  return aepy_make_pytuple("dd", ra, dec);
}




//==============================================================================
// BINDINGS FOR RELATIVITY.C.
//==============================================================================

static PyObject *aepy_relativity(struct aepy_t *self, PyObject *args,
                                 PyObject *keys) {
  char *keywords[] = {"p", "q", "o", NULL};
  struct aepy_3d_t p, q, o;

  strcpy(p.keyword, "p");
  strcpy(q.keyword, "q");
  strcpy(o.keyword, "o");

  if (!PyArg_ParseTupleAndKeywords(args, keys, "O&O&O&", keywords,
                                   &aepy_list_to_3d, (void *)&p,
                                   &aepy_list_to_3d, (void *)&q,
                                   &aepy_list_to_3d, (void *)&o))
    return NULL;

  ae_relativity(p.vec, q.vec, o.vec);

  return aepy_3d_to_list(p.vec);
}





//==============================================================================
// BINDINGS FOR TOPOCENTRIC.C AND REFRAC.C
//==============================================================================

static PyObject *aepy_refrac_visible(struct aepy_t *self, PyObject *args,
                                     PyObject *keys) {
  char *keywords[] = {"alt", "pres", "temp", NULL};
  double alt, pt[2];

  if (!PyArg_ParseTupleAndKeywords(args, keys, "ddd", keywords, &alt, &pt[0],
                                   &pt[1]))
    return NULL;

  return PyFloat_FromDouble(ae_refrac_visible(alt, (void *)pt));
}

static PyObject *aepy_refrac_ulich(struct aepy_t *self, PyObject *args,
                                   PyObject *keys) {
  char *keywords[] = {"alt", "pres", "temp", "humid", NULL};
  double alt, pt[3];

  if (!PyArg_ParseTupleAndKeywords(args, keys, "dddd", keywords, &alt, &pt[0],
                                   &pt[1], &pt[2]))
    return NULL;

  return PyFloat_FromDouble(ae_refrac_ulich(alt, (void *)pt));
}

static PyObject *aepy_topocentric(struct aepy_t *self, PyObject *args,
                                  PyObject *keys) {
  char *keywords[] = {"jd_ut1", "glat", "tlon", "ra", "dec", "dist", "refrac",
                      "jd_tt", "tlat", "trho", NULL};
  int i, n;
  double jd_ut1, glat, tlon, dist, ra, dec, jd_tt, tlat, trho;
  double nutl, nuto, eps, last, alt, az;
  PyObject *refrac, *func, *param, *res;

  dist = 0;
  jd_tt = 1e99;
  tlat = 1e99;
  trho = 1e99;
  refrac = NULL;
  if (!PyArg_ParseTupleAndKeywords(args, keys, "ddddd|dO!ddd", keywords,
                                   &jd_ut1, &glat, &tlon, &ra, &dec, &dist,
                                   &PyTuple_Type, &refrac, &jd_tt, &tlat,
                                   &trho))
    return NULL;

  // Compute any necessary values.
  if (jd_tt > 1e50)
    jd_tt = jd_ut1 + ae_delta_t(jd_ut1) * AE_D_PER_S;
  if (tlat > 1e50 || trho > 1e50)
    ae_geocentric_lat(glat, 0, &tlat, &trho);

  // Get nutation and epsilon.
  eps = ae_epsilon(jd_tt);
  ae_nutation_lon_ob(jd_tt, &nutl, &nuto);

  // Local apparent sidereal time.
  last = ae_last(jd_ut1, jd_tt, tlon, nutl, eps);

  // Correct for diurnal aberration.
  ae_diurnal_aberration(last, tlat, trho, AE_ADD_ABERRATION, &ra, &dec);

  // Diurnal parallax.
  if (dist)
    ae_diurnal_parallax(last, tlat, trho, dist, &ra, &dec);

  // Apply atmospheric refraction correction, if requested.
  if (refrac != NULL) {
    // We need to convert to alt/az.
    ae_radec_to_altaz(last, glat, ra, dec, &alt, &az);

    // The first item of refrac should be the function.
    func = PyTuple_GetItem(refrac, 0);
    if (!PyFunction_Check(func) && !PyCFunction_Check(func)) {
      PyErr_SetString(PyExc_KeyError, "The first item of the tuple refrac must "
                                      "be a function, or a built-in function.");
      return NULL;
    }
    
    // The subsequent items should be the arguments for the function. The first
    // argument, though, is the altitude.
    n = (int)PyTuple_Size(refrac);
    param = PyTuple_New(n);
    PyTuple_SetItem(param, 0, PyFloat_FromDouble(alt));
    for (i = 1; i < n; i++) {
      Py_INCREF(PyTuple_GetItem(refrac, i));
      PyTuple_SetItem(param, i, PyTuple_GetItem(refrac, i));
    }

    // Call the function. The result should be a float.
    if (!PyFloat_Check((res = PyObject_CallObject(func, param)))) {
      PyErr_SetString(PyExc_KeyError, "The refraction function must return a "
                                      "float.");
      Py_DECREF(param);
      Py_DECREF(res);
      return NULL;
    }

    // Now apply the correction.
    alt += PyFloat_AsDouble(res);
    ae_altaz_to_radec(last, glat, alt, az, &ra, &dec);

    // Clean up references.  
    Py_DECREF(param);
    Py_DECREF(res);
  }

  return aepy_make_pytuple("dd", ra, dec);
}




//==============================================================================
// OTHER IMPORTANT BINDINGS
//==============================================================================

static PyObject *aepy_geocentric_from_helio(struct aepy_t *self, PyObject *args,
                                            PyObject *keys) {
  char *keywords[] = {"jd_tt", "o", "v_o", "q", NULL};
  double jd_tt, ra, dec, dist;
  struct aepy_3d_t o, v_o, q;

  strcpy(o.keyword, "o");
  strcpy(v_o.keyword, "v_o");
  strcpy(q.keyword, "q");
  if (!PyArg_ParseTupleAndKeywords(args, keys, "dO&|i", keywords, &jd_tt,
                                   &aepy_list_to_3d, (void *)&o,
                                   &aepy_list_to_3d, (void *)&v_o,
                                   &aepy_list_to_3d, (void *)&q))
    return NULL;

  ae_geocentric(jd_tt, q.vec, o.vec, v_o.vec, &ra, &dec, &dist);

  return aepy_make_pytuple("ddd", ra, dec, dist);
};

static PyObject *aepy_geocentric(struct aepy_t *self, PyObject *args,
                                 PyObject *keys) {
  char *keywords[] = {"jd_tt", "obs", "obj", NULL};
  int i, obj_num, not_planetary;
  double jd_tt, ra, dec, dist, o[3], v_o[3], q[3];
  PyObject *obs, *obj;
  struct aepy_orb_t *orb;
  struct aepy_star_t *star;
  struct aepy_jpl_t *jpl;

  if (!PyArg_ParseTupleAndKeywords(args, keys, "dOO", keywords, &jd_tt,
                                   &obs, &obj))
    return NULL;

  // Get the position and velocity of the observer.
  if ((orb = get_orb_obj(obs))) {
    ae_kepler(jd_tt, orb->o, o);
    ae_v_orbit(jd_tt, orb->o, v_o);
  }
  else if ((jpl = get_jpl_obj(obs, &obj_num, &not_planetary))) {
    if ((i = ae_jpl_get_coords(jpl->j, jd_tt, obj_num, o, v_o,
                               !not_planetary))) {
      PyErr_SetString(PyExc_KeyError, ae_retcode_string(i));
      return NULL;
    }
  }
  else {
    PyErr_SetString(PyExc_KeyError, "The argument 'obs' must be either an "
                                    "'aephem.orb' object, or a tuple of an "
                                    "'aephem.jpl' object and the planet index "
                                    "number.");
    return NULL;
  }

  // Get the light-corrected position of the object and then the geocentric
  // thingy.
  if ((orb = get_orb_obj(obj))) {
    ae_light_t_orbit(jd_tt, o, orb->o, q);
    ae_geocentric(jd_tt, q, o, v_o, &ra, &dec, &dist);
  }
  else if ((jpl = get_jpl_obj(obj, &obj_num, &not_planetary))) {
    if ((i = ae_light_t_jpl(jpl->j, jd_tt, o, obj_num, q, NULL, 1))) {
      PyErr_SetString(PyExc_KeyError, ae_retcode_string(i));
      return NULL;
    }
    ae_geocentric(jd_tt, q, o, v_o, &ra, &dec, &dist);
  }
  else if ((star = get_star_obj(obj))) {
    dist = -1;
    ae_geocentric_from_cat(jd_tt, o, v_o, star->s, &ra, &dec);
  }
  else if (PyString_Check(obj)) {
    if (!strcasecmp(PyString_AsString(obj), "sun")) {
      q[0] = 0;
      q[1] = 0;
      q[2] = 0;
      ae_geocentric(jd_tt, q, o, v_o, &ra, &dec, &dist);
    }
    else if (!strcasecmp(PyString_AsString(obj), "moon")) {
      if ((orb = get_orb_obj(obs)))
        ae_geocentric_moon_from_orbit(jd_tt, orb->o, &ra, &dec, &dist);
      else {
        PyErr_SetString(PyExc_KeyError, "If key 'obj' is set to 'moon', then "
                                        "the 'obs' key must be an 'aephem.orb' "
                                        "object.");
        return NULL;
      }
    }
    else 
      PyErr_SetString(PyExc_KeyError, "The key 'obj' must be one of the "
                                      "following: an 'aephem.orb' object, an "
                                      "'aephem.jpl' object, an 'aephem.star' "
                                      "object, or the string 'sun' or 'moon'.");
  }
  else
    PyErr_SetString(PyExc_KeyError, "The key 'obj' must be one of the "
                                    "following: an 'aephem.orb' object, an "
                                    "'aephem.jpl' object, an 'aephem.star' "
                                    "object, or the string 'sun' or 'moon'.");

  return aepy_make_pytuple("ddd", ra, dec, dist);
}

static PyObject *aepy_light_t(struct aepy_t *self, PyObject *args, 
                              PyObject *keys) {
  char *keywords[] = {"jd_tt", "o", "obj", NULL};
  int obj_num, not_planetary;
  double jd_tt, q[3];
  struct aepy_3d_t o;
  struct aepy_orb_t *orb;
  struct aepy_jpl_t *jpl;
  PyObject *obj;

  strcpy(o.keyword, "o");
  if (!PyArg_ParseTupleAndKeywords(args, keys, "dO&!|i", keywords, &jd_tt,
                                   &aepy_list_to_3d, (void *)&o, &obj))
    return NULL;

  if ((orb = get_orb_obj(obj)))
    ae_light_t_orbit(jd_tt, o.vec, orb->o, q);
  else if ((jpl = get_jpl_obj(obj, &obj_num, &not_planetary)))
    ae_light_t_jpl(jpl->j, jd_tt, o.vec, obj_num, q, NULL, !not_planetary);
  else {
    PyErr_SetString(PyExc_KeyError, "The key 'obj' must be either an "
                                    "'aephem.orb' object or an ('aephem.jpl', "
                                    "obj_num) tuple.");
    return NULL;
  }

  return aepy_3d_to_list(q);
}




//==============================================================================
// BINDINGS FOR PHYSICAL.C
//==============================================================================

static PyObject *aepy_disc_semiminor(struct aepy_t *self, PyObject *args,
                                     PyObject *keys) {
  char *keywords[] = {"a", "c", "latitude", NULL};
  double a, c, lat;

  if (!PyArg_ParseTupleAndKeywords(args, keys, "ddd", keywords, &a, &c, &lat))
    return NULL;

  return PyFloat_FromDouble(ae_disc_semiminor(a, c, lat));
}

static PyObject *aepy_disc_solid_angle(struct aepy_t *self, PyObject *args,
                                       PyObject *keys) {
  char *keywords[] = {"jd_ut1", "obs", "obj", "phys", NULL};
  double jd_ut1;
  struct aepy_orb_t *o_orb, *q_orb;
  struct aepy_phys_t *phys;

  if (!PyArg_ParseTupleAndKeywords(args, keys, "dO!O!O!", keywords, &jd_ut1,
                                   &aepy_orb, &o_orb, &aepy_orb, &q_orb,
                                   &aepy_phys, &phys))
    return NULL;

  return PyFloat_FromDouble(aes_disc_solid_angle(jd_ut1, o_orb->o, q_orb->o,
                                                 phys->p));
}

static PyObject *aepy_subobs_point(struct aepy_t *self, PyObject *args,
                                   PyObject *keys) {
  char *keywords[] = {"jd_ut1", "o_orb", "q_orb", "phys", NULL};
  double jd_ut1, lat, lon, dist;
  struct aepy_orb_t *o_orb, *q_orb;
  struct aepy_phys_t *phys;

  if (!PyArg_ParseTupleAndKeywords(args, keys, "dO!O!O!", keywords, &jd_ut1,
                                   &aepy_orb, &o_orb, &aepy_orb, &q_orb,
                                   &aepy_phys, &phys))
    return NULL;

  aes_subobs_point(jd_ut1, o_orb->o, q_orb->o, phys->p, &lat, &lon, &dist);

  return aepy_make_pytuple("ddd", lat, lon, dist);
}




//==============================================================================
// MODULE SETUP
//==============================================================================

static PyMethodDef aepy_methods[] = {
  {"ctime_to_jd", (PyCFunction)aepy_ctime_to_jd,
   METH_VARARGS | METH_KEYWORDS, AEPY_DS_ctime_to_jd},
  {"ctime_to_last", (PyCFunction)aepy_ctime_to_last,
   METH_VARARGS | METH_KEYWORDS, AEPY_DS_ctime_to_last},
  {"cal_to_jd", (PyCFunction)aepy_cal_to_jd,
   METH_VARARGS | METH_KEYWORDS, AEPY_DS_cal_to_jd},
  {"jd_to_cal", (PyCFunction)aepy_jd_to_cal,
   METH_VARARGS | METH_KEYWORDS, AEPY_DS_jd_to_cal},
  {"gmst", (PyCFunction)aepy_gmst,
   METH_VARARGS | METH_KEYWORDS, AEPY_DS_gmst},
  {"gast", (PyCFunction)aepy_gast,
   METH_VARARGS | METH_KEYWORDS, AEPY_DS_gast},
  {"lmst", (PyCFunction)aepy_lmst,
   METH_VARARGS | METH_KEYWORDS, AEPY_DS_lmst},
  {"last", (PyCFunction)aepy_last,
   METH_VARARGS | METH_KEYWORDS, AEPY_DS_last},
  {"tdb", (PyCFunction)aepy_tdb,
   METH_VARARGS | METH_KEYWORDS, AEPY_DS_tdb},
  {"mjd", (PyCFunction)aepy_mjd,
   METH_VARARGS | METH_KEYWORDS, AEPY_DS_mjd},
  {"delta_t", (PyCFunction)aepy_delta_t,
   METH_VARARGS | METH_KEYWORDS, AEPY_DS_delta_t},
  {"dut1", (PyCFunction)aepy_dut1,
   METH_VARARGS | METH_KEYWORDS, AEPY_DS_dut1},
  {"rect_to_polar", (PyCFunction)aepy_rect_to_polar,
    METH_VARARGS | METH_KEYWORDS, AEPY_DS_rect_to_polar},
  {"polar_to_rect", (PyCFunction)aepy_polar_to_rect,
    METH_VARARGS | METH_KEYWORDS, AEPY_DS_polar_to_rect},
  {"radec_to_altaz", (PyCFunction)aepy_radec_to_altaz,
    METH_VARARGS | METH_KEYWORDS, AEPY_DS_radec_to_altaz},
  {"altaz_to_radec", (PyCFunction)aepy_altaz_to_radec,
    METH_VARARGS | METH_KEYWORDS, AEPY_DS_altaz_to_radec},
  {"radec_to_gal", (PyCFunction)aepy_radec_to_gal,
    METH_VARARGS | METH_KEYWORDS, AEPY_DS_radec_to_gal},
  {"gal_to_radec", (PyCFunction)aepy_gal_to_radec,
    METH_VARARGS | METH_KEYWORDS, AEPY_DS_gal_to_radec},
  {"geocentric_lat", (PyCFunction)aepy_geocentric_lat,
    METH_VARARGS | METH_KEYWORDS, AEPY_DS_geocentric_lat},
  {"delta_q", (PyCFunction)aepy_delta_q,
    METH_VARARGS | METH_KEYWORDS, AEPY_DS_delta_q},
  {"dms", (PyCFunction)aepy_dms,
    METH_VARARGS | METH_KEYWORDS, AEPY_DS_dms},
  {"dms_str", (PyCFunction)aepy_dms_str,
    METH_VARARGS | METH_KEYWORDS, AEPY_DS_dms_str},
  {"hms", (PyCFunction)aepy_hms,
    METH_VARARGS | METH_KEYWORDS, AEPY_DS_hms},
  {"hms_str", (PyCFunction)aepy_hms_str,
    METH_VARARGS | METH_KEYWORDS, AEPY_DS_hms_str},
  {"precess", (PyCFunction)aepy_precess,
    METH_VARARGS | METH_KEYWORDS, AEPY_DS_precess},
  {"nutation_lon_ob", (PyCFunction)aepy_nutation_lon_ob,
    METH_VARARGS | METH_KEYWORDS, AEPY_DS_nutation_lon_ob},
  {"nutate", (PyCFunction)aepy_nutate,
    METH_VARARGS | METH_KEYWORDS, AEPY_DS_nutate},
  {"epsilon", (PyCFunction)aepy_epsilon,
    METH_VARARGS | METH_KEYWORDS, AEPY_DS_epsilon},
  {"constellation", (PyCFunction)aepy_constellation,
    METH_VARARGS | METH_KEYWORDS, AEPY_DS_constellation},
  {"decimal", (PyCFunction)aepy_decimal,
    METH_VARARGS | METH_KEYWORDS, AEPY_DS_decimal},
  {"annual_aberration", (PyCFunction)aepy_annual_aberration,
    METH_VARARGS | METH_KEYWORDS, AEPY_DS_annual_aberration},
  {"diurnal_aberration", (PyCFunction)aepy_diurnal_aberration,
    METH_VARARGS | METH_KEYWORDS, AEPY_DS_diurnal_aberration},
  {"diurnal_parallax", (PyCFunction)aepy_diurnal_parallax,
    METH_VARARGS | METH_KEYWORDS, AEPY_DS_diurnal_parallax},
  {"relativity", (PyCFunction)aepy_relativity,
    METH_VARARGS | METH_KEYWORDS, AEPY_DS_relativity},
  {"refrac_visible", (PyCFunction)aepy_refrac_visible,
    METH_VARARGS | METH_KEYWORDS, AEPY_DS_refrac_visible},
  {"refrac_ulich", (PyCFunction)aepy_refrac_ulich,
    METH_VARARGS | METH_KEYWORDS, AEPY_DS_refrac_ulich},
  {"topocentric", (PyCFunction)aepy_topocentric,
    METH_VARARGS | METH_KEYWORDS, AEPY_DS_topocentric},
  {"disc_semiminor", (PyCFunction)aepy_disc_semiminor,
    METH_VARARGS | METH_KEYWORDS, AEPY_DS_disc_semiminor},
  {"disc_solid_angle", (PyCFunction)aepy_disc_solid_angle,
    METH_VARARGS | METH_KEYWORDS, AEPY_DS_disc_solid_angle},
  {"geocentric_from_helio", (PyCFunction)aepy_geocentric_from_helio,
    METH_VARARGS | METH_KEYWORDS, AEPY_DS_geocentric_from_helio},
  {"geocentric", (PyCFunction)aepy_geocentric,
    METH_VARARGS | METH_KEYWORDS, AEPY_DS_geocentric},
  {"light_t", (PyCFunction)aepy_light_t,
    METH_VARARGS | METH_KEYWORDS, AEPY_DS_light_t},
  {"subobs_point", (PyCFunction)aepy_subobs_point,
    METH_VARARGS | METH_KEYWORDS, AEPY_DS_subobs_point},
  {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC initaephem(void) {
  int i, j;
  char name[128];
  PyObject *mod, *orbit_dict, *phys_dict, *curr_orb, *curr_phys;

  if (PyType_Ready(&aepy_jpl) < 0)
    return;
  if (PyType_Ready(&aepy_orb) < 0)
    return;
  if (PyType_Ready(&aepy_phys) < 0)
    return;
  if (PyType_Ready(&aepy_star) < 0)
    return;

  // Methods.
  mod = Py_InitModule3("aephem", aepy_methods,
                       "Python bindings to the aephem C library.");
  if (mod == NULL)
    return;

  Py_INCREF(&aepy_jpl);
  PyModule_AddObject(mod, "jpl", (PyObject *)&aepy_jpl);
  Py_INCREF(&aepy_orb);
  PyModule_AddObject(mod, "orb", (PyObject *)&aepy_orb);
  Py_INCREF(&aepy_phys);
  PyModule_AddObject(mod, "phys", (PyObject *)&aepy_phys);
  Py_INCREF(&aepy_star);
  PyModule_AddObject(mod, "star", (PyObject *)&aepy_star);

  PyModule_AddObject(mod, "__version__", Py_BuildValue("(iiis)",
                     AE_MAJOR, AE_MINOR, AE_REVISION, AE_VERSION_SUFFIX));
  PyModule_AddStringConstant(mod, "__author__",
            "Adam D. Hincks <adam.hincks@jesuits.ca>");

  // Add a dictionary of built-in orbits and physical data.
  orbit_dict = PyDict_New();
  phys_dict = PyDict_New();
  for (i = 0; i < AE_N_SS_BODIES; i++) {
    if (ae_orb_planet[i] == NULL)
      continue;

    strcpy(name, ae_ss_name[i]);
    for (j = 0; ae_ss_name[i][j] != '\0' && j < 128; j++)
      name[j] = tolower(ae_ss_name[i][j]);

    curr_orb = PyObject_CallObject((PyObject *)&aepy_orb, NULL);
    free(((struct aepy_orb_t *)curr_orb)->o);
    ((struct aepy_orb_t *)curr_orb)->o = ae_orb_planet[i];
    PyDict_SetItem(orbit_dict, PyString_FromString(name), curr_orb);
    Py_DECREF(curr_orb);

    curr_phys = PyObject_CallObject((PyObject *)&aepy_phys, NULL);
    free(((struct aepy_phys_t *)curr_phys)->p->ra_sin_term);
    free(((struct aepy_phys_t *)curr_phys)->p->dec_cos_term);
    free(((struct aepy_phys_t *)curr_phys)->p->w_sin_term);
    free(((struct aepy_phys_t *)curr_phys)->p);
    ((struct aepy_phys_t *)curr_phys)->p = ae_phys_planet[i];
    PyDict_SetItem(phys_dict, PyString_FromString(name), curr_phys);
    Py_DECREF(curr_phys);
  }

  Py_INCREF(orbit_dict);
  PyModule_AddObject(mod, "p_orb", orbit_dict);
  
  Py_INCREF(phys_dict);
  PyModule_AddObject(mod, "p_phys", phys_dict);

  // Add constants.
  for (i = 0; aepy_const_list[i].name != NULL; i++)
    PyModule_AddObject(mod, aepy_const_list[i].name,
                       PyFloat_FromDouble(aepy_const_list[i].value));
  for (i = 0; aepy_enum_list[i].name != NULL; i++)
    PyModule_AddIntConstant(mod, aepy_enum_list[i].name,
                            aepy_enum_list[i].value);

  return;
}
