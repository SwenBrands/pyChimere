#!/bin/env python

'''asdf'''

from netCDF4 import Dataset
from sys import argv
from datetime import timedelta, datetime
import pandas as pd
from numpy import sin, cos, arccos, pi, where
import os
from shutil import copy
from shutil import move

PROD=os.environ['PROD']
PDIR = PROD+"/BDMOD"
INPUT = os.environ['SCRIPTS']

#NAME_FILEIN = argv[1]
DATA = argv[1]
HORA = argv[2]
#PDIR = argv[4]

NAME_FILEIN = "out."+DATA+"_"+HORA+"_gal3.nc"
NAME_FILEOUT = "CHIMERE_"+DATA+HORA+"_gal3.txt"
NAME_FILEEST = INPUT+"/estacions.csv"

ALM = PDIR+"/"+DATA+"/"+HORA

FILEIN = Dataset(NAME_FILEIN, 'r')
FILEOUT = open(NAME_FILEOUT, 'w')


# Id de modelo para el CHIMERE_00Z_GAL3 en
# la base de datos Variables
ID_M = 31

# Ids de las variables en la base de datos Variables
VARS = {
     "NO" : 32,
     "NO2" : 33,
     "NOX" : 36,
     "SO2" : 30,
     "CO" : 31,
     "O3" : 37,
     "PM25" : 34,
     "PM10" : 35
        }

# Leemos los datos de las estaciones
EST = pd.read_csv(NAME_FILEEST, index_col=0)

# Tomamos la lat y lon de las estaciones
LAT = FILEIN.variables["lat"][:]
LON = FILEIN.variables["lon"][:]

# Buscamos el punto mas cercano
def Distance(LAT1, LON1, LAT2, LON2):
    R = 6.371e6
    G2R = 2*pi/360
    RLAT1 = LAT1*G2R
    RLAT2 = LAT2*G2R
    RLON1 = LON1*G2R
    RLON2 = LON2*G2R
    DIST = R * arccos(sin(RLAT1)*sin(RLAT2) + \
                      cos(RLAT1)*cos(RLAT2)*cos(RLON1-RLON2) )
    return DIST

# Tomamos los tiempos del fichero
TIMES_M = FILEIN.variables['Times'][:]
TIMES = []
for ele in TIMES_M:
    TIMES.append(''.join(ele))

# Dimensiones
NT, NZ, NY, NX = FILEIN.variables["NO"].shape

# Recorremos todas las estaciones
for ID_L in EST.index:
    print "Estacion", EST.loc[ID_L, 'NOME']
    # Tomamos la LAT/LON de la estacion
    LAT_EST = EST.loc[ID_L, 'LAT']
    LON_EST = EST.loc[ID_L, 'LON']
    # Calculamos las distancias
    DISTS = Distance(LAT_EST, LON_EST, LAT, LON)
    # Buscamos el punto del minimo
    MIN = where(DISTS == DISTS.min())
    I_MIN = MIN[1][0]
    J_MIN = MIN[0][0]
    # Buscamos os valores no arquivo
    for NVAR, ID_V in VARS.iteritems():
        # Procesando variables
        for t in range(1, NT, 1):
            VALUE = FILEIN.variables[NVAR][t, 0, J_MIN, I_MIN]
            # El CO se mide en mgr/m3, lo convertimos
            if NVAR == "CO":
                VALUE = VALUE/1000
            ID_R = (t -1)/24
            FILEOUT.write("%s %f %d %d %d %d\n"
                          %(TIMES[t].replace("_", " ").replace("-", "/"),
                            VALUE, ID_R, ID_M, ID_V, ID_L))

FILEOUT.close()
        
if not os.path.isdir(ALM):
    os.makedirs(ALM)

print "Copiando %s a %s" %(NAME_FILEOUT,ALM)
copy(NAME_FILEOUT, ALM)

