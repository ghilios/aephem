//! \file jpl.c
//! Implements bindings for routines that use ae_jpl_handle_t.
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

static PyObject *aepy_jpl_create(PyTypeObject *type, PyObject *args,
                                 PyObject *keys) {
  struct aepy_jpl_t *self;

  if ((self = (struct aepy_jpl_t *)type->tp_alloc(type, 0)) != NULL)
    self->j = NULL;

  return (PyObject *)self;
}

static int aepy_jpl_init(struct aepy_jpl_t *self, PyObject *args,
                         PyObject *keys) {
  char *keywords[] = {"path", "ascii", "header_path", "trust_header", NULL};
  char *path, *header_path;
  int i, is_ascii, trust_header;
  
  header_path = NULL;
  is_ascii = 0;
  trust_header = 0;
  if (!PyArg_ParseTupleAndKeywords(args, keys, "s|isi", keywords, &path,
                                   &is_ascii, &header_path, &trust_header))
    return -1;

  self->j = (struct ae_jpl_handle_t *)malloc(sizeof(struct ae_jpl_handle_t));

  if (is_ascii)
    i = ae_jpl_init_ascii(path, header_path, !trust_header, self->j);
  else
    i = ae_jpl_init_bin(path, self->j);

  if (i) {
    PyErr_SetString(PyExc_KeyError, ae_retcode_string(i));
    free(self->j);
    self->j = NULL;
    return -1;
  }

  return 0;
}

static void aepy_jpl_delete(struct aepy_jpl_t *self) {
  if (self->j != NULL)
    ae_jpl_close(self->j);
  free(self->j);

  return;
}

static PyObject *aepy_jpl_get_coords(struct aepy_jpl_t *self, PyObject *args,
                                     PyObject *keys) {
  char *keywords[] = {"jd_tt", "obj_num", "not_planetary", NULL};
  int obj_num, not_planetary;
  double jd_tt, r[3], v[3];
  PyObject *ret;

  not_planetary = 0;
  if (!PyArg_ParseTupleAndKeywords(args, keys, "di|i", keywords, &jd_tt,
                                   &obj_num, &not_planetary))
    return NULL;

  ae_jpl_get_coords(self->j, jd_tt, obj_num, r, v, !not_planetary);

  ret = PyTuple_New(2);
  PyTuple_SetItem(ret, 0, aepy_3d_to_list(r));
  PyTuple_SetItem(ret, 1, aepy_3d_to_list(v));
  
  return ret;
}

static PyObject *aepy_jpl_get_ipt(struct aepy_jpl_t *self, void *closure) {
  int i, j;
  PyObject *ret, *item;

  ret = PyList_New(13);
  for (i = 0; i < 13; i++) {
    item = PyList_New(3);
    for (j = 0; j < 3; j++)
      PyList_SetItem(item, j, PyInt_FromLong((long)self->j->ipt[i][j]));
    PyList_SetItem(ret, i, item);
  }

  return ret;
}

static PyObject *aepy_jpl_get_const(struct aepy_jpl_t *self, void *closure) {
  int i, j;
  char trunc_name[AE_JPL_CONST_LEN + 1];
  PyObject *ret;

  ret = PyDict_New();
  for (i = 0; i < self->j->n_const; i++) {
    for (j = 0; j < AE_JPL_CONST_LEN && self->j->const_name[i][j] != ' '; j++)
      trunc_name[j] = self->j->const_name[i][j];
    trunc_name[j] = '\0';
    PyDict_SetItem(ret, PyString_FromString(trunc_name),
                        PyFloat_FromDouble(self->j->const_val[i]));
  }

  return ret;
}

GET_FUNCTION(jpl, j,             start, Float, Double)
GET_FUNCTION(jpl, j,               end, Float, Double)
GET_FUNCTION(jpl, j,              step, Float, Double)
GET_FUNCTION(jpl, j,                au, Float, Double)
GET_FUNCTION(jpl, j,          em_ratio, Float, Double)
GET_FUNCTION(jpl, j,           n_const,   Int,   Long)
GET_FUNCTION(jpl, j, ephemeris_version,   Int,   Long)


//==============================================================================
// MODULE SETUP
//==============================================================================

static PyMethodDef aepy_jpl_methods[] = {
  {"get_coords", (PyCFunction)aepy_jpl_get_coords,
   METH_VARARGS | METH_KEYWORDS, AEPY_DS_jpl_get_coords},
  {NULL, NULL, 0, NULL}
};

static PyGetSetDef aepy_jpl_getset[] = {
  GET_STRUCT(jpl,             start, "The starting Julian Date.")
  GET_STRUCT(jpl,               end, "The final Julian Date.")
  GET_STRUCT(jpl,              step, "The ephemeris step size (in days).")
  GET_STRUCT(jpl,                au, "The definition of the AU in the "
                                     "ephemeris (in km).")
  GET_STRUCT(jpl,          em_ratio, "The definition of the earth/moon ratio "
                                     "in the ephemeris.")
  GET_STRUCT(jpl,           n_const, "The number of defined constants.")
  GET_STRUCT(jpl, ephemeris_version, "The ephemeris version.")
  GET_STRUCT(jpl,               ipt, "Index pointers to the Chebyshev "
                                     "coefficients.")
  GET_STRUCT(jpl,             const, "The constants defined by the ephemeris.")
  {NULL, NULL, NULL, NULL}
};

PyTypeObject aepy_jpl = {
  PyObject_HEAD_INIT(NULL)
  0,                               // ob_size
  "aephem.jpl",                    // tp_name
  sizeof(struct aepy_jpl_t),       // tp_basicsize
  0,                               // tp_itemsize
  (destructor)aepy_jpl_delete,     // tp_dealloc
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
  AEPY_DS_jpl,                     // tp_doc
  0,                               // tp_traverse
  0,                               // tp_clear
  0,                               // tp_richcompare
  0,                               // tp_weaklistoffset
  0,                               // tp_iter
  0,                               // tp_iternext
  aepy_jpl_methods,                // tp_methods
  0,                               // tp_members
  aepy_jpl_getset,                 // tp_getset
  0,                               // tp_base
  0,                               // tp_dict
  0,                               // tp_descr_get
  0,                               // tp_descr_set
  0,                               // tp_dictoffset
  (initproc)aepy_jpl_init,         // tp_init
  0,                               // tp_alloc
  aepy_jpl_create                  // tp_new
};
