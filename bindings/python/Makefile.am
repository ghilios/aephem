python_LTLIBRARIES   = aephem.la
BUILT_SOURCES        = const.c
AM_CFLAGS            = -Wall -fno-strict-aliasing
AM_CPPFLAGS          = $(CPPFLAGS) $(PYTHON_CPPFLAGS)
aephem_la_SOURCES    = aephem_py.c jpl.c orb.c phys.c star.c
include_HEADERS      = aephem_py.h help.h
aephem_la_LDFLAGS    = -module -avoid-version
aephem_la_LIBADD     = $(PYTHON_LIBS) ../../src/libaephem.la
nodist_aephem_la_SOURCES = const.c

const.c: ../generate_const
	../generate_const "../../src/aephem.h" p > $@

../generate_const: ../generate_const.c ../../src/aephem.h
	cd .. && ${MAKE} generate_const
