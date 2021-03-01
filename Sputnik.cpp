//
//  Sputnik - estimate the orbit of Sputnik-1 from original measurements
//  D.P. Mitchell  14/05/2011.
//
#include "stdafx.h"
#include "SGP4UNIT.H"
#pragma intrinsic(sqrt, sin, cos)
#define D_DTR               0.01745329251994329576923690768488612713
#define D_RTD               57.295779513082320876798154814105
#define D_2PI               6.2831853071795864769252867665590
#define C_EARTHRADIUS       6378137.0               // Meters (GPS WGS 84)
#define C_EARTHECCENTRIC    0.006694379990141317    // e**2 from WGS 84
#define C_MINUTESPERDAY     1440.0
#define C_MILE          1609.344                    // Meters (5280 feet)
#define J2000       2451545.0                       // 2000 January 1.5 TT, 11:58:55.816 UTC
#define JULIANCENTURY 36525.0                       // days per century
#define SOLARDAY    86400.0                         // seconds per day
//
//  Joel Runes's estimated orbital elements for Sputnik-1: 
//  1 00001E 57  1  A 57277.80437500  .00301870  00000-0  00000-0 0    12
//  2 00001  65.1000 340.3821 0520478  58.0000 306.9536 14.96977024    06
//
//  An approximate BSTAR drag term of 8.2e-4 was calculated, using the
//  old FORTRAN SGP8 propagator.  9.2e-4 reflects the new SGP4 model.
//
double EPOCH = 57277.80437500;    // Oct 4, 19:18:18 (virtual asc. node before the launch)
double BSTAR = 9.27e-4;           // Drag term
double XINCL = 65.1000;           // Inclination of orbital plane
double XNODEO = 340.3821;         // Longitude of ascending node
double EO = 0.0520478;            // Eccentricity
double OMEGAO = 58.0000;          // Argument of perigee
double XMO = 306.9536;            // Mean anomaly
double XNO = 14.96977024;         // Mean Motion (revolutions/day)
double DUMMY;
double *pf = &DUMMY;
double r[3], vel[3], fLat, fLong, rStation[3];
int nVerbose;
//
//  Latitude crossings from radio interferometers and radio doppler measurements
//
struct Crossing {
    int     nYear;
    int     nMonth;
    int     nDay;
    int     nHour;
    int     nMinute;
    double  fSeconds;
    double  fLatitude;
    double  fLongitude;
    double  fAltitude;
    char    chStation;
};

#define RAE(D,H,M,S,A) {1957, 10, D, H, M, S, 51.17916667, 0.0, (1853.184*A), 'R'}
#define URB(D,H,M,S) {1957, 10, D, H, M, S, 40.01833333, 0.0, 0.0, 'U'}
#define CAM(D,H,M,S) {1957, 10, D, H, M, S, 51.71, 0.0, 0.0, 'C'}

Crossing rgc[] = {
    URB( 5, 13, 30,  1),        // University of Illinois, Urbana
    URB( 6, 13, 32, 13),
    URB( 7, 13, 33, 55), 
    URB( 8, 13, 35,  9), 
    URB(10, 13, 36,  1), 
    URB(11, 13, 35, 41), 
    URB(12, 11, 58, 50), 
    URB(12, 13, 34, 48), 
    URB(13, 11, 57, 28), 
    URB(14, 11, 55, 31), 
    URB(15, 10, 17, 11), 
    URB(15, 11, 53,  0), 
    URB(16, 10, 14,  9), 
    URB(16, 11, 49, 55), 
    URB(17, 10, 10, 31), 
    URB(17, 11, 46, 15), 
    URB(18, 10,  6, 18), 
    URB(18, 11, 42,  0), 
    URB(19, 10,  1, 29), 
    URB(20,  9, 56,  4), 
    URB(21,  9, 50,  2), 
    URB(22,  9, 43, 22), 
    URB(23,  8,  0, 34),
    URB(24,  7, 52, 41), 
    URB(24,  9, 28,  8), 
    URB(25,  7, 44,  8), 

    URB( 7,  5, 10, 56), 
    URB(10,  3, 37, 33), 
    URB(11,  3, 37, 25), 
    URB(12,  3, 36, 44), 
    URB(13,  1, 59, 38), 
    URB(13,  3, 35, 32), 
    URB(14,  1, 57, 54), 
    URB(14,  3, 33, 46), 
    URB(15,  1, 55, 36), 
    URB(15,  3, 31, 26), 
    URB(16,  1, 52, 45), 
    URB(17,  1, 49, 19), 
    URB(18,  1, 45, 18), 
    URB(19,  1, 40, 40), 
    URB(19, 23, 59, 50), 
    URB(20,  1, 35, 28), 
    URB(20, 23, 54,  3), 
    URB(21,  1, 29, 38), 
    URB(21, 23, 47, 38), 
    URB(22,  1, 23, 11), 
    URB(22, 23, 40, 36), 
    URB(23,  1, 16,  6), 
    URB(23, 23, 32, 55), 
    URB(24,  1,  8, 23), 
    URB(24, 23, 24, 36), 
    URB(25,  1,  0,  1), 
    URB(25, 23, 15, 37), 

    RAE(11,  5, 32, 15, 265),   // Royal Aircraft Establishment
    RAE(13,  5, 30, 16, 252),   //REVIEW: altitude geodetic or above mean earth?
    RAE(14,  5, 28, 30, 253),   // RAE(14,  5, 29, 22, 253), Royal Society values more accurate
    RAE(14, 21, 11, 30, 138),   // RAE(14, 21, 11, 20, 138),
    RAE(15, 21,  8, 43, 133),
    RAE(16, 21,  5, 25.5, 124),
    RAE(17, 21,  1, 38.5, 125),
    RAE(24, 18, 41, 30, 123),

    CAM(15, 21,  8, 56)         // Mullard Observatory, Cambridge
};
//
//  Optical sightings by Moonwatch observatories
//
#define SYDNEY       151.094722,-33.911944,   15
#define SENDAI       140.865555, 38.256111,   45
#define STATECOLL    282.11,     40.7333,    393
#define USNAVALOBS   282.934167, 38.920556,   86
#define SACRAMENTO   238.723611, 38.6405,    180
#define MILWAUKEE    271.851561, 42.9687583, 294
#define HIROSHIMA    132.469444, 34.368889,    3
#define GREENSBORO   280.132661, 36.0773694, 265
#define MUSASHINO    139.576389, 35.716111,   75
#define KURUMEMACH   139.53,     35.755,      52
#define LASCRUCES    253.152778, 32.328333, 1186
#define ASHIGAWA     142.363889, 43.774167,  113
#define CAMBRIDGE    288.870625, 42.379889,   24
#define SKALNATE      26.233333, 49.195,    1783
#define LOSALAMOS    253.677778, 35.875,    2256
#define AMARILLO     258.057778, 35.191667, 1100
#define SACRAPEAK    254.1797947,32.7874639,2823
#define ORGANPASS    253.4477058,32.4234861,1651
#define MILLBROOK    286.375833, 41.858333,  175
#define BALTIMORE    283.409722, 39.407778,    7
#define DOVER        285.470833, 40.958333,  175
#define PEORIA       270.4025,   40.755278,  152

struct Sighting {
    int     nYear;
    int     nMonth;
    int     nDay;
    int     nHour;
    int     nMinute;
    double  fSeconds;
    double  fLongitude;
    double  fLatitude;
    double  fAltitude;
    double  fRightAscend;
    double  fDeclination;
    double  fAzimuth;
    double  fElevation;
};

#define NOVALUE -1000.0

#define RADEC(MO,D,H,M,S, LATLONALT, RD,RM,RS, DD,DM,DS) \
{1957,MO,D,H,M,S, LATLONALT, RD+RM/60.0+RS/3600.0, DD+DM/60.0+DS/3600.0, NOVALUE, NOVALUE}

#define AZELE(MO,D,H,M,S, LATLONALT, ZD,ZM,ZS, ED,EM,ES) \
{1957,MO,D,H,M,S, LATLONALT, NOVALUE, NOVALUE, ZD+ZM/60.0+ZS/3600.0, ED+EM/60.0+ES/3600.0}

Sighting rgs[] = {
    RADEC(10, 9, 9,39,31,   SYDNEY,     15,18, 0, -44,-48, 0),  // One doc says Oct 10, but error is 10x then
    RADEC(10, 9, 9,41,15.5, SYDNEY,     18,30,16, -10,  0, 0),  // Doc says Oct 10, but only works on Oct 9, a1?
    RADEC(10, 9, 9,41,30.5, SYDNEY,     18,35, 0,  -9,  0, 0),
    RADEC(10,13,10,17,28,   MILLBROOK,  10,54, 0,  31, 26, 0),  // large error
    RADEC(10,13,10,20,23,   BALTIMORE,  11,43, 0,  15, 30, 0),  // large error
    RADEC(10,14,19,56,43,   SENDAI,     13,11, 0,  55, 0, 0),   // large error, maybe a3?
    RADEC(10,14,10,18,12,   STATECOLL,  12,17, 0,  57, 0, 0),
    AZELE(10,14,19,52,4.5,  MUSASHINO,   0, 0, 0,  10.7,0,0),
    AZELE(10,15,10,20,40,   USNAVALOBS, 111, 0, 0, 41, 0, 0),
    RADEC(10,15,10,26, 0,   STATECOLL,   11,27, 0, 54,30, 0),
    AZELE(10,15,19,48,57,   KURUMEMACH,  45, 0, 0, 34, 0, 0),
    AZELE(10,16,10, 8, 0,   DOVER,       10,30, 0, 82, 0, 0),   // source 2
    RADEC(10,16,13,27,27,   SACRAMENTO,  7,25,30,   0,54,36),
    RADEC(10,17,13,24,48,   SACRAMENTO,  7, 7,48, -21,-18, 0),
    RADEC(10,18,13,21,33.6, SACRAMENTO,  6,40,30, -31, 0, 0),
    RADEC(10,18,10, 4,53.5, MILWAUKEE,   9,47, 0,  54, 0, 0),
    RADEC(10,18,19,40,23,   HIROSHIMA,   6,18, 0,  81, 0, 0),
    RADEC(10,19, 9,58, 0,   GREENSBORO,  8,40, 0,  79, 0, 0),   // "a2?" large error
    AZELE(10,22, 9,23,51.7, MUSASHINO,  315, 0, 0, 22,30, 0),
    RADEC(10,22, 9,24,21.8, KURUMEMACH, 14, 7, 0,  55,54, 0),
    AZELE(10,22, 1,20,56,   LASCRUCES,  122, 0, 0, 70, 0, 0),
    AZELE(10,23, 9,13,15,   KURUMEMACH, 315, 0, 0, 21, 0, 0),
    AZELE(10,24, 9, 7, 2,   ASHIGAWA,   275, 0, 0, 26, 0, 0),
    AZELE(10,24,22,48, 7.2, CAMBRIDGE,  299, 0, 0, 85,30, 0),  // "399" typo
    AZELE(10,24,22,48,11,   CAMBRIDGE,  312, 0, 0, 88, 0, 0),
    AZELE(10,25,16,57, 3,   SKALNATE,     0, 0, 0, 90, 0, 0),
    RADEC(10,26,16,47,33.6, SKALNATE,   17,37,18,  55,17, 0),
    AZELE(10,26,16,47,53.4, SKALNATE,     0, 0, 0, 52,42, 0),

    AZELE(11,24,22,48, 7.2, CAMBRIDGE,  339, 0, 0, 85,30, 0),   // source2
    AZELE(11,24,22,48,11.0, CAMBRIDGE,  312, 0, 0, 88, 0, 0),   // source2
    RADEC(11,26, 1,25, 1,   LOSALAMOS,  21,14, 0,  10,48, 0),
    RADEC(11,26, 3, 5, 5.6, SACRAMENTO, 21,15,28, -32,-25,-12),
    AZELE(11,27, 0,53, 6.6, AMARILLO,   265, 0, 0, 65, 0, 0),
    RADEC(11,29, 1,14,59,   SACRAPEAK,   19, 9, 0, -17,-36, 0),
    RADEC(11,29, 1,14,13,   ORGANPASS,  18,25,36,  15, 0, 0),
    RADEC(11,29, 1,14,53,   ORGANPASS,  19,24,48,  -6,-48, 0),
    RADEC(11,29,16,16,18,   SKALNATE,   10,23, 0,  66,20, 0),
    RADEC(11,29,16,17,38,   SKALNATE,    6, 4,18,  49, 8, 0),
    AZELE(11,30,12,23,45,   PEORIA,     90, 0, 0,  63, 0, 0)    // source2
};

double
GeodeticToCentric(double fLatitude, double fAltitude)
{
    double f, r, z, fN;


    f = sin(fLatitude*D_DTR);
    fN = C_EARTHRADIUS/sqrt(1.0 - C_EARTHECCENTRIC*f*f);
    r = (fN + fAltitude)*cos(fLatitude*D_DTR);
    z = (fN*(1.0 - C_EARTHECCENTRIC) + fAltitude)*sin(fLatitude*D_DTR);
    return D_RTD * atan(z/r);
}

double
GeocentricLatitude(double fGeodetic)
{
    return D_RTD*atan((1.0 - C_EARTHECCENTRIC)*tan(fGeodetic*D_DTR));
}

double
GeodeticLatitude(double fGeocentric)
{
    return D_RTD*atan(tan(fGeocentric*D_DTR)/(1.0 - C_EARTHECCENTRIC));
}

double
JulianDay(int nYear, int nMonth, int nDay, int nHour, int nMinute, double fSecond)
{
    int a, y, m, jd;

    a = (14 - nMonth)/12;
    y = nYear + 4800 - a;
    m = nMonth + 12*a - 3;
    jd = nDay + (153*m + 2)/5 + 365*y + y/4 - y/100 + y/400 - 32045;
    return jd + (double(nHour-12) + (double(nMinute) + fSecond/60.0)/60.0)/24.0;
}
//
//  Greenwich Mean Sidereal Time
//
double
GMST(double tJulianUTC)
{
    double Du, Dfrac, T, Theta, fSeconds, fDegrees;

    Du = tJulianUTC - J2000;
    Dfrac = fmod(tJulianUTC, 1.0);
    T = (tJulianUTC - J2000)/JULIANCENTURY;
    Theta = fmod(0.7790572732640 + Dfrac + 0.00273781191135448*Du, 1.0);   // Earth Rotation Angle (cycles)
    if (Theta < 0.0)
        Theta += 1.0;
    fSeconds = SOLARDAY*Theta + (0.014506 + T*(4612.156534 + T*(1.3915817 + T*(-0.00000044
                            + T*(-0.000029956 - T*0.0000000368))))) / 15.0; // Precession correction
    fDegrees = 360.0*fSeconds/SOLARDAY;
    if (fDegrees < 0.0)
        fDegrees += 360.0;
    return fDegrees;
}

void
Geodetic(double fLatitude, double fLongitude, double fAltitude, double tJulian)
{
    double f, r, z, fN;

    f = sin(fLatitude*D_DTR);
    fN = C_EARTHRADIUS/sqrt(1.0 - C_EARTHECCENTRIC*f*f);
    r = (fN + fAltitude)*cos(fLatitude*D_DTR);
    z = (fN*(1.0 - C_EARTHECCENTRIC) + fAltitude)*sin(fLatitude*D_DTR);
    rStation[0] = r*cos(fLongitude*D_DTR)/1000.0;
    rStation[1] = r*sin(fLongitude*D_DTR)/1000.0;
    rStation[2] = z/1000.0;
}
//REVIEW: this only works for 1957
double
Epoch(double tJulian)
{
    return 57277.0 + tJulian - JulianDay(1957, 10, 4, 0, 0, 0.0);   // Julian day into NORAD TLE time
}

double
DistanceFromLine(double vOrg[3], double vDir[3], double vPnt[3])
{
    double v[3], d1, d2, d3, d;
    int i;

    d1 = d2 = d3 = 0.0;
    for (i = 0; i < 3; i++) {
        v[i] = vOrg[i] - vPnt[i];
        d1 += v[i]*v[i];
        d2 += vDir[i]*vDir[i];
        d3 += v[i]*vDir[i];
    }
    d = (d1*d2 - d3*d3);    //REVIEW: numercially unstable, but adequate
    if (d < 0.0)
        d = 0.0;
    return sqrt(d/(d2*d2));
}
//
//  Use NORAD's SGP4 orbital propogator, upgraded in 2006 by Vallado et al, to include
//  numerous improvements made over the years.
//
elsetrec satrec;
double tStart = 0, tStop = 5000000.0;

double
SputnikError(double f)
{
    int i, n;
    double tJulian, fAlt, fError, fTotalErr, tEPOCH, tSince, R, rError[3], vDir[3];
    double RA1950, DEC1950, RA2000, DEC2000;

    //
    //  Process radio locations
    //
    n = sizeof(rgc)/sizeof(Crossing);
    fTotalErr = 0;
    *pf = f;
    tEPOCH = EPOCH - 57277 + 2557 + 277;    // days since Jan 0, 1950, 0 hour
    sgp4init(wgs84, 0, tEPOCH, BSTAR, EO, OMEGAO*D_DTR, XINCL*D_DTR,
             XMO*D_DTR, XNO*D_2PI/C_MINUTESPERDAY, XNODEO*D_DTR, satrec);
    for (i = 0; i < n; i++) {
        tJulian = JulianDay(rgc[i].nYear, rgc[i].nMonth, rgc[i].nDay, rgc[i].nHour,
                            rgc[i].nMinute, rgc[i].fSeconds);
        tSince = (Epoch(tJulian) - EPOCH)*1440.0;
        sgp4(wgs84, satrec, tSince, r, vel);
        R = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
        fLat = D_RTD*asin(r[2]/R);
        fLong = fmod(D_RTD*atan2(r[1], r[0]) - GMST(tJulian), 360.0);
        if (fLong < -180.0)
            fLong += 360.0;
        if (rgc[i].fAltitude) {
            fAlt = rgc[i].fAltitude;
        } else {
            if (i <= 25)
                fAlt = 332.0*C_MILE;
            else if (i <= 52)
                fAlt = 141.7*C_MILE;
            else
                fAlt = 300000.0;
        }
        fError = fLat - GeodeticToCentric(rgc[i].fLatitude, fAlt);
        if (nVerbose) printf("%c %4d %2d:%2d:%4.1f  %7.3f %8.3f  %8.4f\n", rgc[i].chStation,
            rgc[i].nDay, rgc[i].nHour, rgc[i].nMinute, rgc[i].fSeconds,
            fLat, fLong, fError);
        if (tJulian < tStart || tJulian > tStop)
            continue;
        fTotalErr += fError*fError;
        if (i < 53) {
            fError = fLong - rgc[i].fLongitude;
            if (nVerbose) printf("   longitude error %f\n", fError);
            fTotalErr += fError*fError;
        }
    }
    tJulian = JulianDay(1957, 10, 15, 21, 8, 41.0);  // Cambridge zero meridian crossing
    tSince = (Epoch(tJulian) - EPOCH)*1440.0;
    sgp4(wgs84, satrec, tSince, r, vel);
    R = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
    fLat = D_RTD*asin(r[2]/R);
    fLong = fmod(D_RTD*atan2(r[1], r[0]) - GMST(tJulian), 360.0);
    if (fLong < -180.0)
        fLong += 360.0;
    fError = fLong - 0.0;
    if (nVerbose) printf("%c %4d %2d:%2d:%4.1f  %7.3f %8.3f  %8.4f\n", 'C',
    15, 21, 8, 41.0, fLat, fLong, fError);
    fTotalErr += fError*fError;
    //
    //  Process optical sightings
    //
    n = sizeof(rgs)/sizeof(Sighting);
    for (i = 0; i < n; i++) {
        tJulian = JulianDay(rgs[i].nYear, rgs[i].nMonth, rgs[i].nDay, rgs[i].nHour,
                            rgs[i].nMinute, rgs[i].fSeconds);
        tSince = (Epoch(tJulian) - EPOCH)*1440.0;
        sgp4(wgs84, satrec, tSince, r, vel);
        fLong = fmod(rgs[i].fLongitude + GMST(tJulian), 360.0);
        if (fLong < -180.0)
            fLong += 360.0;
        Geodetic(rgs[i].fLatitude, fLong, rgs[i].fAltitude, tJulian);
        if (rgs[i].fDeclination == NOVALUE) {
            rError[0] = r[0] - rStation[0];
            rError[1] = r[1] - rStation[1];
            rError[2] = r[2] - rStation[2];
            if (nVerbose) printf("Station %2d position error: %8.2f %8.2f %8.2f\n", 
                i, rError[0], rError[1], rError[2]);
        } else {
            //
            //  Sputnik sightings were in the old Besselian (B1950) coordinate system
            //
            RA1950 = rgs[i].fRightAscend*15.0;
            DEC1950 = rgs[i].fDeclination;
            RA2000 = RA1950 + 0.640265 + 0.278369 * sin(RA1950*D_DTR) * tan(DEC1950*D_DTR);
	        DEC2000 = DEC1950 + 0.278369 * cos(RA1950*D_DTR);
            RA2000 = RA1950;
            DEC2000 = DEC1950;
            vDir[0] = cos(DEC2000*D_DTR)*cos(RA2000*D_DTR);
            vDir[1] = cos(DEC2000*D_DTR)*sin(RA2000*D_DTR);
            vDir[2] = sin(DEC2000*D_DTR);
            if (nVerbose) printf("Station %2d dist from line of sight: %f\n", i,
                DistanceFromLine(rStation, vDir, r));
        }
    }
    fTotalErr = sqrt(fTotalErr);
    if (nVerbose) printf("Total Error = %f\n", (fTotalErr));

    return fTotalErr;
}
//
//  Brent's localmin algorithm
//
double
ML_LocalMinimum(double a, double b, double (*f)(double))
{
	static double eps = 1.0;
	static double gold;
	double c, d, e;
	double p, q, r, u, v, w;
	double mid, tol1, tol2, eps1;
	double fu, fv, fw, fx, x;

	if (eps > 0.5) {
		gold = 0.5*(3.0 - sqrt(5.0));
		do {
			eps *= 0.5;
			eps1 = eps + 1.0;
		} while (eps1 > 1.0);
		eps = sqrt(eps);
	}
	c = gold;
	v = w = x = a + c*(b - a);
	e = 0.0;
	d = b - a;			/* just to keep lint happy */
	fx = (*f)(x);
	fv = fx;
	fw = fx;
	for (;;) {
		mid = 0.5*(a + b);
		tol1 = eps*fabs(x);
		tol2 = tol1 + tol1;
		if (fabs(x - mid) <= tol2 - 0.5*(b - a))
			return x;
		p = q = r = 0.0;
		/*
		 *	Fit parabola
		 */
		if (fabs(e) > tol1) {
			r = (x - w)*(fx - fv);
			q = (x - v)*(fx - fw);
			p = (x - v)*q - (x - w)*r;
			q = 2.0*(q - r);
			if (q > 0.0)
				p = -p;
			else
				q = -q;
			r = e;
			e = d;
		}
		if (fabs(p) < fabs(0.5*q*r) &&
			    p > q*(a - x) && p < q*(b - x)) {
			/*
			 *	parabolic interpolation step
			 */
			d = p/q;
			u = x + d;
			/*
			 *	f must not be evaluated too close to a or b
			 */
			if ((u - a) < tol2 || (b - u) < tol2)
				d = (x < mid) ? tol1 : -tol1;
		} else {
			/*
			 *	golden-section step
			 */
			if (x >= mid)
				e = a - x;
			else
				e = b - x;
			d = c*e;
		}
		/*
		 *	f must not be evaluted to close t x
		 */
		if (fabs(d) >= tol1)
			u = x + d;
		else
			u = x + ((d > 0) ? tol1 : -tol1);
		fu = (*f)(u);
		/*
		 *	update a, b, v, w, x
		 */
		if (fu <= fx) {
			if (u >= x)
				a = x;
			else
				b = x;
			v = w;
			fv = fw;
			w = x;
			fw = fx;
			x = u;
			fx = fu;
		} else {
			if (u < x)
				a = u;
			else
				b = u;
			if (fu <= fw || w == x) {
				v = w;
				fv = fw;
				w = u;
				fw = fu;
			} else if (fu <= fv || v == x || v == w) {
				v = u;
				fv = fu;
			}
		}
	}
}

void
OptimizeSputnik()
{
    double fMin;
    int i;

    for (i = 0; i < 100; i++) {
        pf = &BSTAR;
        nVerbose = 0;
        fMin = ML_LocalMinimum(5e-4, 10e-4, SputnikError);
        nVerbose = 0;
        SputnikError(fMin);
        printf("Optimal BSTAR = %0.8f\n", BSTAR);

        pf = &EPOCH;
        nVerbose = 0;
        fMin = ML_LocalMinimum(57277.80437500-1.0/720.0, 57277.80437500+1.0/720.0, SputnikError);
        nVerbose = 0;
        SputnikError(fMin);
        printf("Optimal EPOCH = %0.8f\n", EPOCH);

        pf = &XNO;
        nVerbose = 0;
        fMin = ML_LocalMinimum(14.9, 15.0, SputnikError);
        nVerbose = (i == 99);
        SputnikError(fMin);
        printf("Optimal PERIOD = %0.8f\n", 1440.0/XNO);
    }
    for (i = 0; i < 500; i++) {
        pf = &OMEGAO;
        nVerbose = 0;
        fMin = ML_LocalMinimum(45.0, 70.0, SputnikError);
        nVerbose = 1;
        SputnikError(fMin);
        printf("Optimal Arg Perigee = %0.8f\n", OMEGAO);

        pf = &EO;
        nVerbose = 0;
        fMin = ML_LocalMinimum(0.05, 0.055, SputnikError);
        nVerbose = 0;
        SputnikError(fMin);
        printf("Optimal Ecc = %0.8f\n", EO);

        pf = &XMO;
        nVerbose = 0;
        fMin = ML_LocalMinimum(300.0, 315.0, SputnikError);
        nVerbose = 0;
        SputnikError(fMin);
        printf("Optimal Anomaly = %0.8f\n", XMO);

        pf = &BSTAR;
        nVerbose = 0;
        fMin = ML_LocalMinimum(5e-4, 10e-4, SputnikError);
        nVerbose = 0;
        SputnikError(fMin);
        printf("Optimal BSTAR = %0.8f\n", BSTAR);

        pf = &EPOCH;
        nVerbose = 0;
        fMin = ML_LocalMinimum(57277.80437500-1.0/720.0, 57277.80437500+1.0/720.0, SputnikError);
        nVerbose = 0;
        SputnikError(fMin);
        printf("Optimal EPOCH = %0.8f\n", EPOCH);

        pf = &XINCL;
        nVerbose = 0;
        fMin = ML_LocalMinimum(64.0, 66.0, SputnikError);
        nVerbose = 0;
        SputnikError(fMin);
        printf("Optimal INCL = %0.8f\n", XINCL);

        pf = &XNODEO;
        nVerbose = 0;
        fMin = ML_LocalMinimum(330.0, 350.0, SputnikError);
        nVerbose = 0;
        SputnikError(fMin);
        printf("Optimal XNODEO = %0.8f\n", XNODEO);

        pf = &XNO;
        nVerbose = 0;
        fMin = ML_LocalMinimum(14.9, 15.0, SputnikError);
        nVerbose = 0;
        SputnikError(fMin);
        printf("Optimal PERIOD = %0.8f\n", 1440.0/XNO);
    }
}

int
OrbitSputnik(double tJulian)
{
    double R, tEPOCH, tSince;
    int n;

    tEPOCH = EPOCH - 57277 + 2557 + 277;    // days since Jan 0, 1950, 0 hour
    sgp4init(wgs84, 0, tEPOCH, BSTAR, EO, OMEGAO*D_DTR, XINCL*D_DTR,
             XMO*D_DTR, XNO*D_2PI/C_MINUTESPERDAY, XNODEO*D_DTR, satrec);
    tSince = (Epoch(tJulian) - EPOCH)*1440.0;
    n = sgp4(wgs84, satrec, tSince, r, vel);
    R = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
    fLat = D_RTD*asin(r[2]/R);
    fLong = fmod(D_RTD*atan2(r[1], r[0]) - GMST(tJulian), 360.0);
    if (fLong < -180.0)
        fLong += 360.0;
    return n;
}

void
AnalyzeSputnik()
{
    double tJulian, tSince, fLastLat, fLastNode, fPeriod;
    int i, nVal, nOrbit;

    pf = &DUMMY;
    tJulian = JulianDay(1957, 10, 4, 0, 0, 0.0) + EPOCH - 57277.0;  // NORAD TLE time into Julian Day
    tSince = (Epoch(tJulian) - EPOCH)*1440.0;                       // Minutes since TLE epoch
    nVal = OrbitSputnik(tJulian);
    printf("%10.3f minutes  %7.3f lat  %8.3f long %d\n", tSince, fLat, fLong, nVal);
    fLastLat = fLat;
    fLastNode = 0.0;
    nOrbit = 0;
    for (i = 1; i < 8000000; i++) {
        nVal = OrbitSputnik(tJulian + double(i)/SOLARDAY);
        if (nVal)
            break;
        if (fLastLat <= 0.0 && fLat > 0.0) {
            fPeriod = double(i) - fLastNode;
            fLastNode = double(i);
            nOrbit++;
        }
        fLastLat = fLat;
        if ((i % 864000) == 1)
            printf("%4d days  %0.2f period  %4d orbits\n", i/86400, fPeriod/60.0, nOrbit);
    }
    printf("nVal = %d, %d days\n", nVal, i);
}
//
//  Convert between True, Eccentric and Mean anomalies.
//
double
ML_TrueFromEccentricAnomaly(double E, double e)
{
    return 2.0*atan2(sqrt(1.0 + e)*sin(0.5*E), sqrt(1.0 - e)*cos(0.5*E));
}

double
ML_EccentricFromTrueAnomaly(double T, double e)
{
    return 2.0*atan2(sqrt(1.0 - e)*sin(0.5*T), sqrt(1.0 + e)*cos(0.5*T));
}

double
ML_MeanFromEccentricAnomaly(double E, double e)
{
    return E - e*sin(E);
}

double
ML_EccentricFromMeanAnomaly(double M, double e)
{
    double E, Delta, Step;
    int i;

    E = 0.0;
    M = fmod(M, D_2PI);
    for (i = 0; i < 25; i++) {
        Step = M + e*sin(E);
        Delta = (Step - E)/(1.0 - e*cos(E));    // Newton method step
        if (fabs(Delta) > 2.0) {
            E = Step;                           // Simple interation if Newton looks unstable
        } else {
            E = E + Delta;
        }
        if (fabs(Delta) < 5.0e-13)
            break;
    }
    return E;
}

void
TestAnomaly()
{
    double e, EDeg, E, M, EfromM, T, EfromT;

    for (e = 0.0001; e < 1.0; e += 0.005) {
        for (EDeg = 0.0; EDeg <= 360.0; EDeg += 0.005) {
            E = D_DTR*EDeg;
            M = ML_MeanFromEccentricAnomaly(E, e);
            EfromM = ML_EccentricFromMeanAnomaly(M, e);
            if (fabs(E - EfromM) > 1.0e-12)
                printf("%8.5f %6.1f %g\n", e, EDeg, EfromM - E);
            T = ML_TrueFromEccentricAnomaly(E, e);
            EfromT = ML_EccentricFromTrueAnomaly(T, e);
            if (fabs(E - EfromT) > 1.0e-12)
                printf("%8.5f %6.1f %g\n", e, EDeg, EfromT - E);
        }
    }
}

int _tmain(int argc, _TCHAR* argv[])
{
    extern int Urbana();

    TestAnomaly(); return 0;
    Urbana();
    tStart = JulianDay(1957, 10, 4, 0, 0, 0.0);
    tStop = JulianDay(1957, 10, 13, 0, 0, 0.0);
    OptimizeSputnik();
    AnalyzeSputnik();
	return 0;
}

