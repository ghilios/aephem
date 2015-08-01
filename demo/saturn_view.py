# This example shows some of the basic features of the python bindings.
# The the ra/dec, distance, solid angle and sub-observer position of Saturn as
# observed from the Atacama Cosmology Telesocpe is calculated, assuming a
# visible-light atmospheric refraction correction.

import aephem as a
from time import *

# Set the latitude, longitude and altitude of the Atacama Cosmology Telescope.
# Also set some atmospheric conditions.
lat = a.decimal(-22, 57, 31)
lon = a.decimal(-67, 47, 15)
alt = 5196
pres = 550
temp = -5.2

# Get the current Julian date in UT1 and TT measures.
jd_utc = a.ctime_to_jd(time())
jd_ut1 = jd_utc + a.dut1(jd_utc) * a.D_PER_S
jd_tt = jd_ut1 + a.delta_t(jd_ut1) * a.D_PER_S
print "Julian Date in UT1 is %.8f; delta T is %.2f sec." % \
      (jd_ut1, a.delta_t(jd_ut1))

# Get the geocentric ra/dec of Saturn, using the built-in ephemerides.
gra, gdec, dist = a.geocentric(jd_tt, a.p_orb["earth"], a.p_orb["saturn"])

# Now get topocentric coordinates, taking into account visible-light refraction.
# First, though, get the geocentric latitude (even though this can be done
# automatically if the optional arguments are not passed to
# aephem.topocentric().
tlat, trho = a.geocentric_lat(lat, alt)
ra, dec = a.topocentric(jd_ut1, lat, lon, gra, gdec, \
                        dist = dist, jd_tt = jd_tt, \
                        refrac = (a.refrac_visible, pres, temp))
print "Ra/dec of Saturn is: %s, %s." % \
      (a.hms_str(ra), a.dms_str(dec))

# Get the local apparent sidereal time and then convert ra/dec to alt/az.
last = a.last(jd_ut1, lon)
alt, az = a.radec_to_altaz(last, lat, ra, dec)
print "Alt/az of Saturn is: %s, %s" % (a.dms_str(alt), a.dms_str(az))

# Now get the subobserver coordinates.
slat, slon, dist = a.subobs_point(jd_ut1, a.p_orb["earth"], a.p_orb["saturn"], \
                                  a.p_phys["saturn"])
print "Sub-earth lat/lon of Saturn is: %s, %s." % \
      (a.dms_str(slat), a.dms_str(slon))

# Finally, get the angular size, correcting the polar radius for the viewing
# angle.
eq = a.p_phys["saturn"].r_eq / dist / a.AU * a.RTD
pol = a.p_phys["saturn"].r_pole / dist / a.AU * a.RTD
pol = a.disc_semiminor(eq, pol, slat)
print "Saturn disc semi-major and -minor angular sizes: %.3f\", %.3f\"" % \
      (eq * 3600.0, pol * 3600.0)
print "Distance to Saturn: %.5f AU." % (dist)
