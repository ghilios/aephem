#
# Makefile for the ephemeris.com ephcom package to use JPL ephemeris files.
#
# This should run as is on Linux, and other Unix systems.  If there are
# any problems with a particular Unix system, email bugs@ephemeris.com.
#
# Set JPLDIR to the location of the JPL ephemeris files, copied from
# the CD from Willman-Bell or from JPL's FTP site:
#
#      ftp://ssd.jpl.nasa.gov/pub/eph/export
#
# You should be able to download and use the binary files from that site
# directly with this software, regardless of whether your host machine is
# Big-Endian or Little-Endian.
#
VERSION=1.0
DATE = `date "+%Y.%m.%d"`
JPLDIR = ../jplcd
CC = gcc
CFLAGS = -g -O2

PROGS = eph2eph asc2eph eph2asc testeph ephcmp headcmp ephcoeff vtransit
MAN1PAGES = $(PROGS)

INSTALL = install
MAN2HTML = man2html
BINDIR = /usr/local/bin
MAN1DIR = /usr/local/man

all: $(PROGS)

eph2eph: eph2eph.c gnulliver.o ephcom.h ephcom.o
	$(CC) $(CFLAGS) eph2eph.c gnulliver.o ephcom.o -o eph2eph

asc2eph: asc2eph.c gnulliver.o ephcom.h ephcom.o
	$(CC) $(CFLAGS) asc2eph.c gnulliver.o ephcom.o -o asc2eph

eph2asc: eph2asc.c gnulliver.o ephcom.h ephcom.o
	$(CC) $(CFLAGS) eph2asc.c gnulliver.o ephcom.o -o eph2asc

testeph: testeph.c gnulliver.o ephcom.h ephcom.o
	$(CC) $(CFLAGS) testeph.c gnulliver.o ephcom.o -o testeph

ephcmp: ephcmp.c gnulliver.o ephcom.h ephcom.o
	$(CC) $(CFLAGS) ephcmp.c gnulliver.o ephcom.o -o ephcmp

headcmp: headcmp.c gnulliver.o ephcom.h ephcom.o
	$(CC) $(CFLAGS) headcmp.c gnulliver.o ephcom.o -o headcmp

ephcoeff: ephcoeff.c gnulliver.o ephcom.h ephcom.o
	$(CC) $(CFLAGS) ephcoeff.c gnulliver.o ephcom.o -o ephcoeff

# vtransit calls sqrt(), and so must link with the math library (-lm)
vtransit: vtransit.c gnulliver.o ephcom.h ephcom.o
	$(CC) $(CFLAGS) vtransit.c gnulliver.o ephcom.o -o vtransit -lm

# html:
# 	for FILE in $(MAN1PAGES) ; do \
# 	   nroff -man $$FILE.1 | $(MAN2HTML) > $$FILE.html ;\
# 	done

html:
	for FILE in $(MAN1PAGES) ; do \
	   $(MAN2HTML) $$FILE.1 > $$FILE.html ;\
	done

install: all
	for FILE in $(PROGS) ; do \
	   $(INSTALL) $$FILE -m 0755 $(BINDIR) ;\
	done
	for FILE in $(MAN1PAGES) ; do \
	   $(INSTALL) $$FILE.1 -m 0755 $(MAN1DIR)/man1 ;\
	done

tar: distclean
	/bin/rm -f ephcom-$(VERSION)
	mkdir ephcom-$(VERSION)
	cp LICENSE.txt README.txt Makefile \
	   ephcom.h ephcom.c gnulliver.c \
	   asc2eph.c eph2asc.c eph2eph.c testeph.c ephcmp.c headcmp.c \
	   ephcoeff.c vtransit.c \
	   eph2eph.1 eph2asc.1 asc2eph.1 testeph.1 headcmp.1 \
	   ephcmp.1 ephcoeff.1 vtransit.1 \
	   eph2eph.html eph2asc.html asc2eph.html testeph.html \
	   headcmp.html ephcmp.html ephcoeff.html vtransit.html \
	   ephcom-$(VERSION)
	tar cf ephcom-$(VERSION)-$(DATE).tar ephcom-$(VERSION)

test: all test200 test405 test406 vtransit
	./vtransit $(JPLDIR)/UNIX.405
	./ephcoeff $(JPLDIR)/UNIX.405 2453164.8477853 >test/vtransit-coeffs

test200: all
	if [ ! -d test ] ; then mkdir test ; \
	else /bin/rm -f test/th.200 test/td.200 \
	   test/teph.200 test/ascout.200 test/thcmp.200 \
	   test/testout.200 test/testunix.200 ; \
	fi
	./eph2asc $(JPLDIR)/UNIX.200 test/th.200 test/td.200
	./asc2eph test/th.200 test/td.200 test/teph.200 >test/ascout.200
	./headcmp test/teph.200 $(JPLDIR)/UNIX.200 test/thcmp.200
	./ephcmp test/teph.200 $(JPLDIR)/UNIX.200
	./testeph test/teph.200 $(JPLDIR)/TESTPO.200 \
	    test/testout.200
	./testeph $(JPLDIR)/UNIX.200 $(JPLDIR)/TESTPO.200 \
	   test/testunix.200
	diff test/testout.200 test/testunix.200

test405: all
	if [ ! -d test ]; then mkdir test; \
	else /bin/rm -f test/th.405 test/td.405 \
	   test/teph.405 test/ascout.405 test/thcmp.405 \
	   test/testout.405 test/testunix.405 ; \
	fi
	./eph2asc $(JPLDIR)/UNIX.405 test/th.405 test/td.405
	./asc2eph test/th.405 test/td.405 test/teph.405 \
	   >test/ascout.405
	./headcmp test/teph.405 $(JPLDIR)/UNIX.405 test/thcmp.405
	./ephcmp test/teph.405 $(JPLDIR)/UNIX.405
	./testeph test/teph.405 $(JPLDIR)/TESTPO.405 \
	    test/testout.405
	./testeph $(JPLDIR)/UNIX.405 $(JPLDIR)/TESTPO.405 \
	   test/testunix.405
	diff test/testout.405 test/testunix.405

test406: all
	if [ ! -d test ]; then mkdir test; \
	else /bin/rm -f test/th.406 test/td.406 \
	   test/teph.406 test/ascout.406 test/thcmp.406 \
	   test/testout.406 test/testunix.406 ; \
	fi
	./eph2asc $(JPLDIR)/UNIX.406 test/th.406 test/td.406
	./asc2eph test/th.406 test/td.406 test/teph.406 \
	   >test/ascout.406
	./headcmp test/teph.406 $(JPLDIR)/UNIX.406 test/thcmp.406
	./ephcmp test/teph.406 $(JPLDIR)/UNIX.406
	./testeph test/teph.406 $(JPLDIR)/TESTPO.406 \
	    test/testout.406
	./testeph $(JPLDIR)/UNIX.406 $(JPLDIR)/TESTPO.406 \
	   test/testunix.406
	diff test/testout.406 test/testunix.406

#
# Note: "make exercise" is way too rigorous for the casual user.
# It is only intended for thorough testing of I/O in a new distribution.
# "make test" suffices for testing functionality of all ephcom routines
# used by all accompanying programs, including interpolation of the
# Chebyshev coefficients.
#
# eph2eph: can't read from stdin or write to stdout
# eph2asc: can't read from stdin or write to stdout
# asc2eph: can read from stdin; can't write to stdout (test stdin read)
# testeph: can't read from stdin; can write to stdout (test stdout write)
# headcmp: can't read from stdin; can write to stdout (test stdout write)
# ephcmp:  can't read from stdin; can write to stdout (test stdout write)
# vtransit: can't read from stdin; always writes to stdout
#
exercise: all test405
	if [ ! -d test ]; then mkdir test; \
	else /bin/rm -f test/teph2000.405 test/tiny0.405 test/tiny1.405 \
	   test/tiny01.405 ; \
	fi
	cat	$(JPLDIR)/HEADER.405     \
		$(JPLDIR)/ASCP2000.405   \
		$(JPLDIR)/ASCP2020.405 | \
		./asc2eph - - test/teph2000.405 >>/dev/null
	./testeph test/teph2000.405 $(JPLDIR)/TESTPO.405
	/bin/rm -f test/teph2000.405
	./headcmp test/teph.405 $(JPLDIR)/UNIX.405 >>/dev/null
	./ephcmp  test/teph.405 $(JPLDIR)/UNIX.405 >>/dev/null
	./eph2eph test/teph.405 test/teph.405 test/tiny0.405 0 2510000
	ls -l test/tiny0.405
	./eph2eph test/teph.405 test/tiny0.405 test/tiny1.405 2510000 2520000
	ls -l test/tiny1.405
	./eph2eph test/tiny0.405 test/tiny1.405 test/tiny01.405
	ls -l test/tiny01.405
	/bin/rm -f test/teph2000.405 test/tiny0.405 test/tiny1.405 \
	   test/tiny01.405
	./eph2eph test/teph.405 test/teph.405 test/tinye.405 2500000 2500001
	ls -l test/tinye.405
	./asc2eph test/th.405 test/td.405 test/tinya.405 2500000 2500001 \
	   >>/dev/null
	ls -l test/tinya.405
	./ephcmp test/tinye.405 test/tinya.405
	/bin/rm -f test/tinya.405 test/tinye.405

clean:
	/bin/rm -rf test
	/bin/rm -f *.o

distclean:
	/bin/rm -rf test
	/bin/rm -f *.o
	for FILE in $(PROGS) ; do \
	   /bin/rm -f $$FILE ;\
	done
