//! \file star.c
//! Implements bindings for routines that use ae_star_t.
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

static PyObject *aepy_star_create(PyTypeObject *type, PyObject *args,
                                  PyObject *keys) {
  struct aepy_star_t *self;

  if ((self = (struct aepy_star_t *)type->tp_alloc(type, 0)) != NULL)
    self->s = NULL;

  return (PyObject *)self;
}

static int aepy_star_init(struct aepy_star_t *self, PyObject *args,
                          PyObject *keys) {
  char *keywords[] = {"cat", "star", NULL};
  const char *cat, *star;
  int i;
  
  cat = NULL;
  star = NULL;
  if (!PyArg_ParseTupleAndKeywords(args, keys, "|ss", keywords, &cat, &star))
    return -1;

  self->s = (struct ae_star_t *)malloc(sizeof(struct ae_star_t));;
  if (cat != NULL && star != NULL) {
    if ((i = ae_read_star_from_cat(cat, star, self->s))) {
      PyErr_SetString(PyExc_KeyError, ae_retcode_string(i));
      free(self->s);
      self->s = NULL;
      return -1;
    }
  }
  else if (cat != NULL || star != NULL) {
    PyErr_SetString(PyExc_KeyError, "Either both or neither of optional "
                                    "arguments 'cat' and 'star' must be "
                                    "passed.");
    free(self->s);
    self->s = NULL;
    return -1;
  }
  else {
    strcpy(self->s->name, "none");
    self->s->epoch = AE_J2000;
    self->s->ra = 0;
    self->s->dec = 0;
    self->s->px = 0;
    self->s->mura = 0;
    self->s->mudec = 0;
    self->s->v = 0;
    self->s->mag = 0;
  }

  return 0;
}

static void aepy_star_delete(struct aepy_star_t *self) {
  free(self->s);

  return;
}

static PyObject *aepy_star_read_from_cat(struct aepy_star_t *self,
                                         PyObject *args, PyObject *keys) {
  char *keywords[] = {"cat", "star", NULL};
  const char *cat, *star;
  int i;

  cat = NULL;
  star = NULL;
  if (!PyArg_ParseTupleAndKeywords(args, keys, "ss", keywords, &cat, &star))
    return NULL;

  if ((i = ae_read_star_from_cat(cat, star, self->s))) {
    PyErr_SetString(PyExc_KeyError, ae_retcode_string(i));
    return NULL;
  }

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *aepy_star_fk4_to_fk5(struct aepy_star_t *self, PyObject *args,
                                      PyObject *keys) {
  char *keywords[] = {NULL};

  if (!PyArg_ParseTupleAndKeywords(args, keys, "", keywords))
    return NULL;

  ae_fk4_to_fk5(self->s);

  Py_INCREF(Py_None);
  return Py_None;
}

GETSET_STR_FUNCTIONS(star, s, name)
GETSET_FUNCTIONS(star, s, epoch, Float, Double)
GETSET_FUNCTIONS(star, s,    ra, Float, Double)
GETSET_FUNCTIONS(star, s,   dec, Float, Double)
GETSET_FUNCTIONS(star, s,    px, Float, Double)
GETSET_FUNCTIONS(star, s,  mura, Float, Double)
GETSET_FUNCTIONS(star, s, mudec, Float, Double)
GETSET_FUNCTIONS(star, s,     v, Float, Double)
GETSET_FUNCTIONS(star, s,   mag, Float, Double)


//==============================================================================
// MODULE SETUP
//==============================================================================

static PyMethodDef aepy_star_methods[] = {
  {"read_from_cat", (PyCFunction)aepy_star_read_from_cat,
   METH_VARARGS | METH_KEYWORDS, AEPY_DS_star_read_from_cat},
  {"fk4_to_fk5", (PyCFunction)aepy_star_fk4_to_fk5,
   METH_VARARGS | METH_KEYWORDS, AEPY_DS_star_fk4_to_fk5},
  {NULL, NULL, 0, NULL}
};

static PyGetSetDef aepy_star_getset[] = {
  GETSET_STRUCT(star, name, "Object name (up to 31 chars).")
  GETSET_STRUCT(star,  epoch, "Epoch of coordinates (default J2000).")
  GETSET_STRUCT(star,     ra, "Right ascension, degrees.")
  GETSET_STRUCT(star,    dec, "Declination, degrees.")
  GETSET_STRUCT(star,     px, "Parallax, seconds of arc.")
  GETSET_STRUCT(star,   mura, "Proper motion in ra, \"/century.")
  GETSET_STRUCT(star,  mudec, "Proper motion in dec, \"/century.")
  GETSET_STRUCT(star,      v, "Radial velocity, km/s.")
  GETSET_STRUCT(star,    mag, "Visual magnitude.")
  {NULL, NULL, NULL, NULL}
};

PyTypeObject aepy_star = {
  PyObject_HEAD_INIT(NULL)
  0,                               // ob_size
  "aephem.star",                   // tp_name
  sizeof(struct aepy_star_t),      // tp_basicsize
  0,                               // tp_itemsize
  (destructor)aepy_star_delete,    // tp_dealloc
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
  AEPY_DS_star,                    // tp_doc
  0,                               // tp_traverse
  0,                               // tp_clear
  0,                               // tp_richcompare
  0,                               // tp_weaklistoffset
  0,                               // tp_iter
  0,                               // tp_iternext
  aepy_star_methods,               // tp_methods
  0,                               // tp_members
  aepy_star_getset,                // tp_getset
  0,                               // tp_base
  0,                               // tp_dict
  0,                               // tp_descr_get
  0,                               // tp_descr_set
  0,                               // tp_dictoffset
  (initproc)aepy_star_init,        // tp_init
  0,                               // tp_alloc
  aepy_star_create                 // tp_new
};
