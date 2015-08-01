//! \file phys.
//! Implements bindings for routines that use ae_physical_t.
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

//==============================================================================
// MODULE METHODS
//==============================================================================

static PyObject *aepy_phys_create(PyTypeObject *type, PyObject *args,
                                  PyObject *keys) {
  struct aepy_phys_t *self;

  self = (struct aepy_phys_t *)type->tp_alloc(type, 0);

  return (PyObject *)self;
}

static int aepy_phys_init(struct aepy_phys_t *self, PyObject *args,
                          PyObject *keys) {
  char *keywords[] = {NULL};
  
  if (!PyArg_ParseTupleAndKeywords(args, keys, "", keywords))
    return -1;

  self->p = (struct ae_physical_t *)malloc(sizeof(struct ae_physical_t));

  self->p->r_mean = 0;
  self->p->r_eq = 0;
  self->p->r_pole = 0;
  self->p->rms_dev_spheroid = 0;
  self->p->max_elev = 0;
  self->p->max_depress = 0;
  self->p->pole_ra = 0;
  self->p->pole_ra_t = 0;
  self->p->ra_sin_term =
    (struct ae_physical_term_t *)malloc(sizeof(struct ae_physical_term_t));
  self->p->pole_dec = 0;
  self->p->pole_dec_t = 0;
  self->p->dec_cos_term =
    (struct ae_physical_term_t *)malloc(sizeof(struct ae_physical_term_t));
  self->p->w = 0;
  self->p->w_d = 0;
  self->p->w_d_sq = 0;
  self->p->w_sin_term =
    (struct ae_physical_term_t *)malloc(sizeof(struct ae_physical_term_t));

  self->p->ra_sin_term->time_var = AE_PHYSICAL_END;
  self->p->dec_cos_term->time_var = AE_PHYSICAL_END;
  self->p->w_sin_term->time_var = AE_PHYSICAL_END;

  return 0;
}

static void aepy_phys_delete(struct aepy_phys_t *self) {
  free(self->p->ra_sin_term);
  free(self->p->dec_cos_term);
  free(self->p->w_sin_term);
  free(self->p);

  return;
}

static PyObject *aepy_phys_pole(struct aepy_phys_t *self, PyObject *args,
                                PyObject *keys) {
  char *keywords[] = {"jd_ut1", "jd_tt", NULL};
  double jd_ut1, jd_tt, n[3], w;
  PyObject *ret;

  if (!PyArg_ParseTupleAndKeywords(args, keys, "dd", keywords, &jd_ut1, &jd_tt))
    return NULL;

  ae_phys_pole(jd_ut1, jd_tt, self->p, n, &w);
  ret = PyTuple_New(2);
  PyTuple_SetItem(ret, 0, aepy_3d_to_list(n));
  PyTuple_SetItem(ret, 1, PyFloat_FromDouble(w));

  return ret;
}

static PyObject *aepy_phys_is_retrograde(struct aepy_phys_t *self,
                                         PyObject *args, PyObject *keys) {
  char *keywords[] = {NULL};

  if (!PyArg_ParseTupleAndKeywords(args, keys, "", keywords))
    return NULL;

  if (ae_is_retrograde(self->p)) {
    Py_INCREF(Py_True);
    return Py_True;
  }
  else {
    Py_INCREF(Py_False);
    return Py_False;
  }
}

static PyObject *aepy_phys_flattening(struct aepy_phys_t *self, PyObject *args,
                                      PyObject *keys) {
  char *keywords[] = {NULL};

  if (!PyArg_ParseTupleAndKeywords(args, keys, "", keywords))
    return NULL;

  return PyFloat_FromDouble(ae_flattening(self->p));
}

int aepy_phys_fill_term(PyObject *list, struct ae_physical_term_t **t,
                        const char *var_name) {
  char tmpstr[128], type;
  int i, j, n;
  double coef[4];
  PyObject *tup;

  if (!PyList_Check(list)) {
    sprintf(tmpstr, "'aephem.phys' attribute '%s' must of of type 'list'",
            var_name);
    PyErr_SetString(PyExc_TypeError, tmpstr);
    return -1;
  }

  n = PyList_Size(list);
  *t = (struct ae_physical_term_t *)realloc(*t, n *
                                            sizeof(struct ae_physical_term_t));
  for (i = 0; i < n; i++) {
    tup = PyList_GetItem(list, i);
    if (!PyTuple_Check(tup)) {
      sprintf(tmpstr, "item %d of 'aephem.phys' attribute '%s' must be of type "
                      "'tuple'", i + 1, var_name);
      PyErr_SetString(PyExc_TypeError, tmpstr);
      (*t)[0].time_var = AE_PHYSICAL_END;
      return -1;
    }

    if (!PyString_Check(PyTuple_GetItem(tup, 0))) {
      sprintf(tmpstr, "first item of tuple %d of 'aephem.phys' attribute "
                      "'%s' must be either 't' or 'd'", i + 1, var_name);
      PyErr_SetString(PyExc_TypeError, tmpstr);
      (*t)[0].time_var = AE_PHYSICAL_END;
      return -1;
    }
    sscanf(PyString_AsString(PyTuple_GetItem(tup, 0)), "%c", &type);
    switch (type) {
      case 'd':
        (*t)[i].time_var = AE_PHYSICAL_D;
        break;
      case 't':
        (*t)[i].time_var = AE_PHYSICAL_T;
        break;
      default:
        sprintf(tmpstr, "first item of tuple %d of 'aephem.phys' attribute "
                        "'%s' must be either 't' or 'd'", i + 1, var_name);
        PyErr_SetString(PyExc_TypeError, tmpstr);
        (*t)[0].time_var = AE_PHYSICAL_END;
        return -1;
    }

    for (j = 0; j < 4; j++) {
      if (!PyFloat_Check(PyTuple_GetItem(tup, j + 1)) &&
          !PyInt_Check(PyTuple_GetItem(tup, j + 1))) {
        sprintf(tmpstr, "item %d of tuple %d of 'aephem.phys' attribute "
                        "'%s' must be of type 'float'", j + 2, i + 1, 
                        var_name);
        PyErr_SetString(PyExc_TypeError, tmpstr);
        (*t)[0].time_var = AE_PHYSICAL_END;
        return -1;
      }
      coef[j] = PyFloat_AsDouble(PyTuple_GetItem(tup, j + 1));
    }

    (*t)[i].a = coef[0];
    (*t)[i].b = coef[1];
    (*t)[i].c = coef[2];
    (*t)[i].d = coef[3];
  }

  return 0;
}
static PyObject *aepy_phys_get_term(struct ae_physical_term_t *t) {
  int n;
  char type[2];
  PyObject *ret, *tup;

  // Count number of terms.
  for (n = 0; t[n].time_var != AE_PHYSICAL_END; n++);

  ret = PyList_New(n);
  for (n = 0; t[n].time_var != AE_PHYSICAL_END; n++) {
    tup = PyTuple_New(5);
    strcpy(type, t[n].time_var == AE_PHYSICAL_T ? "t" : "d");
    PyTuple_SetItem(tup, 0, PyString_FromString(type));
    PyTuple_SetItem(tup, 1, PyFloat_FromDouble(t[n].a));
    PyTuple_SetItem(tup, 2, PyFloat_FromDouble(t[n].b));
    PyTuple_SetItem(tup, 3, PyFloat_FromDouble(t[n].c));
    PyTuple_SetItem(tup, 4, PyFloat_FromDouble(t[n].d));
    PyList_SetItem(ret, n, tup);
  }

  return ret;
}

static int aepy_phys_set_ra_sin_term(struct aepy_phys_t *self, PyObject *val,
                                     void *closure) {
  return aepy_phys_fill_term(val, &self->p->ra_sin_term, "ra_sin_term");
}

static int aepy_phys_set_dec_cos_term(struct aepy_phys_t *self, PyObject *val,
                                     void *closure) {
  return aepy_phys_fill_term(val, &self->p->dec_cos_term, "dec_cos_term");
}

static int aepy_phys_set_w_sin_term(struct aepy_phys_t *self, PyObject *val,
                                    void *closure) {
  return aepy_phys_fill_term(val, &self->p->w_sin_term, "w_sin_term");
}

static PyObject *aepy_phys_get_ra_sin_term(struct aepy_phys_t *self,
                                           PyObject *val, void *closure) {
  return aepy_phys_get_term(self->p->ra_sin_term);
}

static PyObject *aepy_phys_get_dec_cos_term(struct aepy_phys_t *self,
                                            PyObject *val, void *closure) {
  return aepy_phys_get_term(self->p->dec_cos_term);
}

static PyObject *aepy_phys_get_w_sin_term(struct aepy_phys_t *self,
                                          PyObject *val, void *closure) {
  return aepy_phys_get_term(self->p->w_sin_term);
}

GETSET_FUNCTIONS(phys, p,           r_mean, Float, Double)
GETSET_FUNCTIONS(phys, p,             r_eq, Float, Double)
GETSET_FUNCTIONS(phys, p,           r_pole, Float, Double)
GETSET_FUNCTIONS(phys, p, rms_dev_spheroid, Float, Double)
GETSET_FUNCTIONS(phys, p,         max_elev, Float, Double)
GETSET_FUNCTIONS(phys, p,      max_depress, Float, Double)
GETSET_FUNCTIONS(phys, p,          pole_ra, Float, Double)
GETSET_FUNCTIONS(phys, p,        pole_ra_t, Float, Double)
GETSET_FUNCTIONS(phys, p,         pole_dec, Float, Double)
GETSET_FUNCTIONS(phys, p,       pole_dec_t, Float, Double)
GETSET_FUNCTIONS(phys, p,                w, Float, Double)
GETSET_FUNCTIONS(phys, p,              w_d, Float, Double)
GETSET_FUNCTIONS(phys, p,           w_d_sq, Float, Double)

//==============================================================================
// MODULE SETUP
//==============================================================================

static PyMethodDef aepy_phys_methods[] = {
  {"pole", (PyCFunction)aepy_phys_pole,
   METH_VARARGS | METH_KEYWORDS, AEPY_DS_phys_pole},
  {"is_retrograde", (PyCFunction)aepy_phys_is_retrograde,
   METH_VARARGS | METH_KEYWORDS, AEPY_DS_phys_is_retrograde},
  {"flattening", (PyCFunction)aepy_phys_flattening,
   METH_VARARGS | METH_KEYWORDS, AEPY_DS_phys_flattening},
  {NULL, NULL, 0, NULL}
};

static PyGetSetDef aepy_phys_getset[] = {
  GETSET_STRUCT(phys,      r_mean, "Mean radius, km. Not required.")
  GETSET_STRUCT(phys,        r_eq, "Equatorial radius, km.")
  GETSET_STRUCT(phys,      r_pole, "Polar radius, km.")
  GETSET_STRUCT(phys, rms_dev_spheroid, "RMS deviation from spheroid, km. Not "
                                        "required.")
  GETSET_STRUCT(phys,    max_elev, "Maximum elevation, km. Not required.")
  GETSET_STRUCT(phys, max_depress, "Maximum depression, km. Not required.")
  GETSET_STRUCT(phys,     pole_ra, "North pole ra, J2000.")
  GETSET_STRUCT(phys,   pole_ra_t, "North pole ra coefficient for T = Julian "
                                   "centuries since J2000.")
  GETSET_STRUCT(phys,    pole_dec, "North pole declination, J2000.")
  GETSET_STRUCT(phys,  pole_dec_t, "North pole dec coefficient for T = Julian "
                                   "centuries since J2000.")
  GETSET_STRUCT(phys,           w, "Prime meridean.")
  GETSET_STRUCT(phys,         w_d, "Prime meridean coefficient for d = days "
                                   "since J2000.")
  GETSET_STRUCT(phys,      w_d_sq, "Prime meridean coefficient for d^2.")
  GETSET_STRUCT(phys, ra_sin_term, "Sine terms for ra: a list of tuples "
                                   "('d'/'t', a, b, c, d), where the\nfirst "
                                   "item indicates whether the time variable x "
                                   "is days ('d') or\nyears ('t') since J2000, "
                                   "and the terms are: a sin(b + c x + d x^2).")
  GETSET_STRUCT(phys, dec_cos_term, "Cosine terms for dec: a list of tuples "
                                   "as in 'ra_sin_term', though\nwith cosine "
                                   "instead of sine.")
  GETSET_STRUCT(phys,  w_sin_term, "Sine terms for w: a list of tuples "
                                   "as in 'ra_sin_term'.")
  {NULL, NULL, NULL, NULL}
};

PyTypeObject aepy_phys = {
  PyObject_HEAD_INIT(NULL)
  0,                               // ob_size
  "aephem.phys",                   // tp_name
  sizeof(struct aepy_phys_t),      // tp_basicsize
  0,                               // tp_itemsize
  (destructor)aepy_phys_delete,    // tp_dealloc
  0,                               // tp_print
  0,                               // tp_getattr
  0,                               // tp_setattr
  0,                               // tp_compare
  0,                               // tp_repr
  0,                               // tp_as_number
  0,                               // tp_as_sequence
  0,                               // tp_as_mapping
  0,                               // tp_hash
  0,                               // tp_call
  0,                               // tp_str
  0,                               // tp_getattro
  0,                               // tp_setattro
  0,                               // tp_as_buffer
  Py_TPFLAGS_DEFAULT,              // tp_flags
  AEPY_DS_phys,                    // tp_doc
  0,                               // tp_traverse
  0,                               // tp_clear
  0,                               // tp_richcompare
  0,                               // tp_weaklistoffset
  0,                               // tp_iter
  0,                               // tp_iternext
  aepy_phys_methods,                // tp_methods
  0,                               // tp_members
  aepy_phys_getset,                 // tp_getset
  0,                               // tp_base
  0,                               // tp_dict
  0,                               // tp_descr_get
  0,                               // tp_descr_set
  0,                               // tp_dictoffset
  (initproc)aepy_phys_init,         // tp_init
  0,                               // tp_alloc
  aepy_phys_create                  // tp_new
};
