lib_LTLIBRARIES		  = libaephem.la
libaephem_la_SOURCES  = annual_ab.c brightness.c calendar.c constel.c coord.c \
                        diurnal_ab.c diurnal_px.c epsilon.c fk4_to_fk5.c \
                        init.c jpl_eph.c geocentric.c gplan.c kepler.c \
                        kfiles.c light_t.c math_util.c nutate.c phys_data.c \
                        physical.c precess.c refrac.c relativity.c time_data.c \
                        topocentric.c v_orbit.c mer404.c ven404.c ear404.c \
                        mar404.c jup404.c sat404.c ura404.c nep404.c plu404.c \
                        mlr404.c mlat404.c
AM_CFLAGS             = -Wall -lm
libaephem_la_LDFLAGS  = -lm -version-info $(AEPHEM_VERSION) \
						-release $(AEPHEM_RELEASE)
include_HEADERS       = aephem.h

noinst_PROGRAMS 	  = parse_iers
parse_iers_SOURCES    = parse_iers.c
parse_iers_LDFLAGS    = -lm
