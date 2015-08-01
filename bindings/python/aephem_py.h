//! \file aephem_py.h
//! Header file for python bindings of AEPHEM library.
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

#ifndef AEPHEM_PY_H
#define AEPHEM_PY_H

#include <Python.h>

#include "../../src/aephem.h"

struct aepy_const_t {
  char *name;
  double value;
};

struct aepy_t {
  PyObject_HEAD
};

struct aepy_star_t {
  PyObject_HEAD
  struct ae_star_t *s;
};

struct aepy_orb_t {
  PyObject_HEAD
  struct ae_orbit_t *o;
};

struct aepy_jpl_t {
  PyObject_HEAD
  struct ae_jpl_handle_t *j;
};

struct aepy_phys_t {
  PyObject_HEAD
  struct ae_physical_t *p;
};

struct aepy_3d_t {
  char keyword[64];
  double vec[3];
};

extern const struct aepy_const_t aepy_const_list[];
extern const struct aepy_const_t aepy_enum_list[];
extern PyTypeObject aepy_jpl;
extern PyTypeObject aepy_orb;
extern PyTypeObject aepy_phys;
extern PyTypeObject aepy_star;

int aepy_list_to_3d(PyObject *obj, void *p);
PyObject *aepy_3d_to_list(double *v);
PyObject *aepy_make_pytuple(const char *type, ...);

#define GET_FUNCTION(type, memb, var, pyvar_type, cvar_type) \
  static PyObject *aepy_##type##_get_##var(struct aepy_##type##_t *self,\
                                           void *closure) { \
    return Py##pyvar_type##_From##cvar_type(self->memb->var); \
  }

#define SET_FUNCTION(type, memb, var, pyvar_type, cvar_type) \
  static int aepy_##type##_set_##var(struct aepy_##type##_t *self,\
                                     PyObject *val, void *closure) { \
    if (!Py##pyvar_type##_Check(val)) { \
      PyErr_SetString(PyExc_TypeError, "'aephem.type' attribute 'val' must " \
                                       "be of type " #pyvar_type); \
      return -1; \
    } \
    self->memb->var = Py##pyvar_type##_As##cvar_type(val); \
    return 0; \
  }

#define GETSET_FUNCTIONS(type, memb, var, pyvar_type, cvar_type) \
  GET_FUNCTION(type, memb, var, pyvar_type, cvar_type) \
  SET_FUNCTION(type, memb, var, pyvar_type, cvar_type)

#define GETSET_STR_FUNCTIONS(type, memb, var) \
  static PyObject *aepy_##type##_get_##var(struct aepy_##type##_t *self,\
                                           void *closure) { \
    return PyString_FromString(self->memb->var); \
  } \
  static int aepy_##type##_set_##var(struct aepy_##type##_t *self,\
                                     PyObject *val, void *closure) { \
    if (!PyString_Check(val)) { \
      PyErr_SetString(PyExc_TypeError, "'aephem.type' attribute 'val' must " \
                                       "be of type String"); \
      return -1; \
    } \
    strcpy(self->memb->var, PyString_AsString(val)); \
    return 0; \
  }

#define GETSET_STRUCT(type, var, docstring) \
  {#var, (getter)aepy_##type##_get_##var, (setter)aepy_##type##_set_##var, \
   docstring},

#define GET_STRUCT(type, var, docstring) \
  {#var, (getter)aepy_##type##_get_##var, NULL, \
   docstring},

#endif
