#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <aephem.h>

#define LATITUDE     -12.345     // Degrees east of the meridian.
#define LONGITUDE    -67.890     // Degrees north of the equator.
#define ALTITUDE     5200        // Metres above sea level.
 
int main(int argc, char *argv[]) {
  double jd_ut, jd_ut1, jd_tt, last, ra, dec, dist, alt, az;
 
  jd_ut = ae_ctime_to_jd(time(NULL));  // Current Julian date in UT measure.
  jd_ut1 = jd_ut + ae_dut1(jd_ut) * AE_D_PER_S;     // UT1 measure.
  jd_tt = jd_ut1 + ae_delta_t(jd_ut1) * AE_D_PER_S; // TT measure.
  last = aes_last(jd_ut1, LONGITUDE);   // Local apparent sidereal time.
 
  // Get the geocentric ra/dec of Mars.
  ae_geocentric_from_orbit(jd_tt, &ae_orb_earth, &ae_orb_mars, &ra, &dec,
                           &dist);
 
  // Get the apparent ra/dec of Mars (no atmospheric model).
  aes_topocentric(jd_ut1, LATITUDE, LONGITUDE, dist, &ra, &dec);
 
  // Convert to alt/az.
  ae_radec_to_altaz(last, LATITUDE, ra, dec, &alt, &az);
 
  printf("JD (UT1) = %.5f, JD (TT) = %.5f, LAST = %.5f deg\n", jd_ut1, jd_tt,
         last);
  printf("Apparent ra/dec of Mars = " ae_hms_fmt ", " ae_dms_fmt "\n",
         ae_hms_arg(ra), ae_dms_arg(dec));
  printf("Az/alt of Mars = " ae_dms_fmt ", " ae_dms_fmt "\n", 
         ae_dms_arg(az), ae_dms_arg(alt));
 
  return 0;
}
