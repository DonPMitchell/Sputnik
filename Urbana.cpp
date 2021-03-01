//
//  Urbana - reconstruct the longitude values from Urbana interferometer
//
#include "stdafx.h"

extern double GeodeticToCentric(double fLatitude, double fAltitude);
extern double GeodeticLatitude(double fGeocentric);

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

extern Crossing rgc[];

struct Ephemeris {
    int nDay;
    int nHour;
    int nMinute;
    int nSecond;
    double fAlt;
    double fCentralLat;
    double fLong;
    double fGeodeticLat;
    int nOrbit;
};

Ephemeris rge[9998];
//
//  Delta longitudes: Urbana interometer observations - NRL ephemeris calculated values
//
double rgfDeltaLongitude[] = {
     0,
    -0.07,
    -0.13,
    -0.17,
    -0.29,
    -1.83,
    -0.38,
     0,
    -0.04,
    -0.08,
    -0.09,
    -0.10,
    -0.13,
    -0.13,
    +0.01,
    +0.33,
    +0.31,
     0,
     0,
     0,
    -0.02,
    +0.04,
    +0.09,
    //+0.09,  duplicated in Swenson's paper
    +0.15,
    +0.13,
    +0.18,

    +0.54,
    +0.95,
    +0.44,
    +0.29,
    +0.13,
    +0.13,
     0,
     0,
    -0.13,
    -0.15,
    -0.25,
     0,
    -0.11,
    -0.20,
    -0.06,
    -0.03,
    -0.38,
     0,
    -0.05,
    -0.04,
    -0.07,
    -0.10,
    -0.11,
    -0.10,
    -0.16,
    -0.16,
    -0.21
};
//
//  Interpolate longitude values at actual Urbana latitude crossings.  Then
//  apply their published deltas from ephemeris valuyes to reconstruct their
//  original longitude values calculated from interferometer measurments.
//  Altitude is height about geoid, but latitude is geocentric, not geodetic.
//
static int
Load()
{
    FILE *pf;
    int i, nDay, nHour, nMinute, nSecond, nOrbit, nLines;
    double fLat, fLong, fAlt;
    char rgch[128];

    pf = fopen("SputnikEphemeris.txt", "r");
    nLines = 0;
    for (i = 0; i < 1000000; i++) {
        if (fgets(rgch, 128, pf) == 0)
            break;
        if (rgch[0] == 'S' || rgch[0] == 'P')
            continue;
        if (sscanf(rgch, "%d%d%d%d%lf%lf%lf%d", &nDay, &nHour, &nMinute, &nSecond,
                    &fAlt, &fLat, &fLong, &nOrbit) != 8)
                    printf("Bad line: %s\n", rgch);
        rge[nLines].nDay = nDay;
        rge[nLines].nHour = nHour;
        rge[nLines].nMinute = nMinute;
        rge[nLines].nSecond = nSecond;
        rge[nLines].fAlt = fAlt;            // altitude above geoid
        rge[nLines].fCentralLat = fLat;     // sub-satellite geocentric latitude
        rge[nLines].fLong = fLong;
        rge[nLines].nOrbit = nOrbit;
        rge[nLines].fGeodeticLat = GeodeticLatitude(fLat);
        nLines++;
    }
    printf("%d ephemeris lines loaded\n", nLines);
    return i;
}
//
//  Cubic interpolation of longitude from NRL ephermeris
//
double fLat0, fLat1, fLat2, fLat3;

static double
Bisection(double a, double c, double (*f)(double))
{
    double b, fa, fb, fc;
    int i;

    fa = f(a);
    fc = f(c);
    for (i = 0; i < 48; i++) {
        b = a + (c-a)*0.5;
        fb = f(b);
        if (fa*fb < 0.0)
            c = b;
        else
            a = b;
    }
    return b;
}

static double
UrbanaLatitude(double t)
{
    return 0.5*((2.0*fLat1) +                                      // Catmull-Rom spline
               (-fLat0 + fLat2)*t +
               (2.0*fLat0 - 5.0*fLat1 + 4.0*fLat2 - fLat3)*t*t +
               (-fLat0 + 3.0*fLat1 - 3.0*fLat2 + fLat3)*t*t*t) - 40.01833333;
}

double
ProcessCrossing(Crossing &c)
{
    int i, nCrossSec, nEphemSec1, nEphemSec2;
    double fLon0, fLon1, fLon2, fLon3;
    double t, fLong;

    nCrossSec = int(c.fSeconds) + 60*(c.nMinute + 60*(c.nHour + 24*c.nDay));
    nEphemSec1 = 0;
    for (i = 0; i < 9998; i++) {
        nEphemSec2 = rge[i].nSecond + 60*(rge[i].nMinute + 60*(rge[i].nHour + 24*rge[i].nDay));
        if (nEphemSec1 <= nCrossSec && nCrossSec < nEphemSec2)
            break;
        nEphemSec1 = nEphemSec2;
    }
    fLat0 = rge[i-2].fGeodeticLat;
    fLat1 = rge[i-1].fGeodeticLat;
    fLat2 = rge[i+0].fGeodeticLat;
    fLat3 = rge[i+1].fGeodeticLat;
    fLon0 = rge[i-2].fLong;
    fLon1 = rge[i-1].fLong;
    fLon2 = rge[i+0].fLong;
    fLon3 = rge[i+1].fLong;
    t = Bisection(0.0, 1.0, UrbanaLatitude);
    fLong = 0.5*((2.0*fLon1) +                                      // Catmull-Rom spline
                 (-fLon0 + fLon2)*t +
                 (2.0*fLon0 - 5.0*fLon1 + 4.0*fLon2 - fLon3)*t*t +
                 (-fLon0 + 3.0*fLon1 - 3.0*fLon2 + fLon3)*t*t*t);
    //printf("Latitudes: %f  40.01833333  %f  %f\n", fLat1, fLat2, UrbanaLatitude(t));
    //printf("   Longitudes: %f %f %f\n", fLon1, fLon2, fLong);
    return fLong;
}
//
//  Reconstruct the original longitude data for the Urbana interferometer measurements
//
int
Urbana()
{
    int i;

    Load();
    printf("%d longitude corrections\n", sizeof(rgfDeltaLongitude)/sizeof(double));
    rgc[0].fLongitude = rge[0].fLong;
    for (i = 1; i < 53; i++) {
        rgc[i].fLongitude = ProcessCrossing(rgc[i]) + rgfDeltaLongitude[i];

    }
    return 1;
}