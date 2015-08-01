//! \file orb.
//! Implements bindings for routines that use ae_orbit_t.
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

static PyObject *aepy_orb_create(PyTypeObject *type, PyObject *args,
                                 PyObject *keys) {
  struct aepy_orb_t *self;

  if ((self = (struct aepy_orb_t *)type->tp_alloc(type, 0)) != NULL)
    self->o = NULL;

  return (PyObject *)self;
}

static int aepy_orb_init(struct aepy_orb_t *self, PyObject *args,
                         PyObject *keys) {
  char *keywords[] = {"cat", "orb", NULL};
  const char *cat, *orb;
  int i;
  
  cat = NULL;
  orb = NULL;
  if (!PyArg_ParseTupleAndKeywords(args, keys, "|ss", keywords, &cat, &orb))
    return -1;

  self->o = (struct ae_orbit_t *)malloc(sizeof(struct ae_orbit_t));

  if (cat != NULL && orb != NULL) {
    if ((i = ae_read_orbit_from_cat(cat, orb, self->o))) {
      PyErr_SetString(PyExc_KeyError, ae_retcode_string(i));
      free(self->o);
      self->o = NULL;
      return -1;
    }
  }
  else if (cat != NULL || orb != NULL) {
    PyErr_SetString(PyExc_KeyError, "Either both or neither of optional "
                                    "arguments 'cat' and 'orb' must be "
                                    "passed.");
    free(self->o);
    self->o = NULL;
    return -1;
  }
  else {
    strcpy(self->o->name, "none");
    self->o->epoch = AE_J2000;
    self->o->i = 0;
    self->o->W = 0;
    self->o->w = 0;
    self->o->a = 0;
    self->o->dm = 0;
    self->o->ecc = 0;
    self->o->M = 0;
    self->o->equinox = 0;
    self->o->ptable = NULL;
    self->o->L = 0;
    self->o->r = 0;
    self->o->plat = 0;
  }

  return 0;
}

static void aepy_orb_delete(struct aepy_orb_t *self) {
  if (self->o != NULL)
    free(self->o->ptable);
  free(self->o);

  return;
}

static PyObject *aepy_orb_read_from_cat(struct aepy_orb_t *self,
                                        PyObject *args, PyObject *keys) {
  char *keywords[] = {"cat", "orb", NULL};
  const char *cat, *orb;
  int i;

  if (!PyArg_ParseTupleAndKeywords(args, keys, "ss", keywords, &cat, &orb))
    return NULL;

  if ((i = ae_read_orbit_from_cat(cat, orb, self->o))) {
    PyErr_SetString(PyExc_KeyError, ae_retcode_string(i));
    return NULL;
  }

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *aepy_orb_kepler(struct aepy_orb_t *self, PyObject *args,
                                 PyObject *keys) {
  char *keywords[] = {"jd_tt", "velocity", NULL};
  int i, velocity;
  double jd_tt, dt, q[3], q2[3], v[3];
  PyObject *ret;

  velocity = 0;
  if (!PyArg_ParseTupleAndKeywords(args, keys, "d|i", keywords, &jd_tt,
                                   &velocity))
    return NULL;

  ae_kepler(jd_tt, self->o, q);

  if (velocity) {
    dt = 0.005;
    ae_kepler(jd_tt - dt, self->o, q2);
    for (i = 0; i < 3; i++)
      v[i] = (q[i] - q2[i]) / dt;
  }

  if (velocity) {
    ret = PyTuple_New(2);
    PyTuple_SetItem(ret, 0, aepy_3d_to_list(q));
    PyTuple_SetItem(ret, 1, aepy_3d_to_list(v));
    return ret;
  }
  else
    return aepy_3d_to_list(q);
}

GETSET_STR_FUNCTIONS(orb, o, name)
GETSET_FUNCTIONS(orb, o,   epoch, Float, Double)
GETSET_FUNCTIONS(orb, o,       i, Float, Double)
GETSET_FUNCTIONS(orb, o,       W, Float, Double)
GETSET_FUNCTIONS(orb, o,       w, Float, Double)
GETSET_FUNCTIONS(orb, o,       a, Float, Double)
GETSET_FUNCTIONS(orb, o,      dm, Float, Double)
GETSET_FUNCTIONS(orb, o,     ecc, Float, Double)
GETSET_FUNCTIONS(orb, o,       M, Float, Double)
GETSET_FUNCTIONS(orb, o, equinox, Float, Double)

//==============================================================================
// MODULE SETUP
//==============================================================================

static PyMethodDef aepy_orb_methods[] = {
  {"read_from_cat", (PyCFunction)aepy_orb_read_from_cat,
   METH_VARARGS | METH_KEYWORDS, AEPY_DS_orb_read_from_cat},
  {"kepler", (PyCFunction)aepy_orb_kepler,
   METH_VARARGS | METH_KEYWORDS, AEPY_DS_orb_kepler},
  {NULL, NULL, 0, NULL}
};

static PyGetSetDef aepy_orb_getset[] = {
  GETSET_STRUCT(orb, name, "Object name (up to 31 chars).")
  GETSET_STRUCT(orb,  epoch, "Epoch of coordinates (default J2000).")
  GETSET_STRUCT(orb,       i, "Inclination.")
  GETSET_STRUCT(orb,       W, "Longitude of the ascending node.")
  GETSET_STRUCT(orb,       w, "Argument of the perihelion.")
  GETSET_STRUCT(orb,       a, "Mean distance (semimajor axis).")
  GETSET_STRUCT(orb,      dm, "Daily motion.")
  GETSET_STRUCT(orb,     ecc, "Eccentricity.")
  GETSET_STRUCT(orb,       M, "Mean anomaly.")
  GETSET_STRUCT(orb, equinox, "Epoch of equinox and ecliptic.")
  {NULL, NULL, NULL, NULL}
};

PyTypeObject aepy_orb = {
  PyObject_HEAD_INIT(NULL)
  0,                               // ob_size
  "aephem.orb",                    // tp_name
  sizeof(struct aepy_orb_t),       // tp_basicsize
  0,                               // tp_itemsize
  (destructor)aepy_orb_delete,     // tp_dealloc
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
  AEPY_DS_orb,                     // tp_doc
  0,                               // tp_traverse
  0,                               // tp_clear
  0,                               // tp_richcompare
  0,                               // tp_weaklistoffset
  0,                               // tp_iter
  0,                               // tp_iternext
  aepy_orb_methods,                // tp_methods
  0,                               // tp_members
  aepy_orb_getset,                 // tp_getset
  0,                               // tp_base
  0,                               // tp_dict
  0,                               // tp_descr_get
  0,                               // tp_descr_set
  0,                               // tp_dictoffset
  (initproc)aepy_orb_init,         // tp_init
  0,                               // tp_alloc
  aepy_orb_create                  // tp_new
};
