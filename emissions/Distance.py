#!/bin/env python

# import numpy to get this to work!

def Distance(LAT1, LON1, LAT2, LON2):
    R = 6.371e6
    G2R = 2*np.pi/360
    RLAT1 = LAT1*G2R
    RLAT2 = LAT2*G2R
    RLON1 = LON1*G2R
    RLON2 = LON2*G2R
    DIST = R * np.arccos(np.sin(RLAT1)*np.sin(RLAT2) + \
                      np.cos(RLAT1)*np.cos(RLAT2)*np.cos(RLON1-RLON2) )
    MINID = np.argmin(DIST)
    MIN = np.min(DIST)
    return(MINID,MIN)
