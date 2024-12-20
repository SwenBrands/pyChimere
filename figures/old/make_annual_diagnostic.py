import sys
import netCDF4
import numpy
import pandas as pd
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, cm
import os
from pylab import *
from datetime import date, timedelta, datetime
from dateutil.rrule import rrule, DAILY
import time
import matplotlib.colors as col
#meto en el path el directorio donde guardo modulos
sys.path.append('modules')
import shapefile
from osgeo import osr, gdal
cpoolSO2 = ['#13A7C7', '#BEE0AE','#F28686','#CE1126']
import subprocess
##CAPA ZONAS
##http:\\www.meteogalicia.es/geoserver/wms?service=WMS&version=1.1.0&request=GetMap&layers=calidadeaire:Zonas

#geoserver de meteogalicia
geoserver="http:\\www.meteogalicia.es/geoserver/ows"
#directorio donde estan los datos del chimere
directorio='/mnt/netapp2/Store_uni/home/xunta/mca/ase/CHIMERE'
#fichero que usaremos como base para saber si un punto es de mar o no
ficherosuelo = '/mnt/netapp2/Store_uni/home/xunta/mca/ase/capas/LANDUSE_GLCF_GAL3.nc'
domain='gal3'
chimere_type='out'
#Definimos la fecha de inicio de la evaluacion. Solo indimamos el inicio ya que el final sera el mismo anho pero a 31 de diciembre
ano=2017
mes=1
dia=1

ndias=0

#Esta variable nos va a servir como base para extraer informacion de las mayas
f=Dataset('/mnt/netapp2/Store_uni/home/xunta/mca/ase/CHIMERE/2017/01/20170101/00/out.20170101_00_gal3.nc','r') 

#Variables que vamos a utilizar para almacenar los distintos datasets
#dataset con las medias de NO2
matrizNO2Media=None
#dataset con la media de NOx
matrizNOxMedia=None
#dataset con las superaciones de NO2
matrizNO2Sup=None
#dataset con las medias de PM10
matrizPM10Media=None
#dataset con las superaciones de PM10
matrizPM10Sup=None
#dataset con los valores de la AOT40 de O3
matrizAOT40=None
#dataset con los maximos octohorarios de O3
matrizMaxO3=None
#dataset con las superacines de maximos octohorarios de O3
matrizSupMaxO3=None
#dataset con las superaciones de CO
matrizCO=None
#dataset con las medias de PM25
matrizPM25Media=None
#dataset con las superaciones de la media de PM25
matrizPM25Sup=None
#dataset con la media de SO2
matrizSO2Med=None
#dataset con las superaciones horarias de SO2
matrizSO2SupH=None
#dataset con las superaciones diarias de SO2
matrizSO2SupD=None
#dataset con las medias invernales de SO2
matrizSO2Inv=None



# Escala de cores
SEISMIC = plt.get_cmap('seismic_r')

def crearStringDia(n):
    salida=str(n)
    if(n<10):       
        salida='0'+salida
    return salida

def crearStringMes(n):
    salida=str(n)
    if(n<10):       
        salida='0'+salida
    return salida


def verVariables():
    f=Dataset('/mnt/netapp2/Store_uni/home/xunta/mca/ase/CHIMERE/2017/01/20170101/00/out.20170101_00_gal3.nc','r')
    for var in f.variables:
        print var

def crearNombreFichero(ano, mes, dia):
    smes = crearStringMes(mes)
    sdia = crearStringDia(dia)
    return "/"+str(ano)+"/"+smes+"/"+str(ano)+smes+sdia+"/00/"+chimere_type+"."+str(ano)+smes+sdia+"_00_"+domain+".nc"

def daterange(start_date, end_date):
    for n in range(int ((end_date - start_date).days)):
        yield start_date + timedelta(n)

def pintarFechas():
    fechaInicio = date(ano, mes, dia)
    fechaFin = datetime.date(ano,12,31)
    for dt in rrule(DAILY, dtstart=fechaInicio, until=fechaFin):
        print dt.strftime("%Y-%m-%d")
    

#la salida de la funcion son dos matrices, Superaciones Horarias, Media Anual
def crearMatricesNO2():
    print "Procedemos a crear las matrices de datos de NOx"
        #para almacenar la media anual
    salidaNO2Media=None

        #para almacenar las superaciones de 200 horario
    salidaNO2Sup=None
    fechaInicio = date(ano, mes, dia)
    fechaFin = date(ano,12,31)
    ndias=0
    for dt in rrule(DAILY, dtstart=fechaInicio, until=fechaFin):
        nombreFichero = crearNombreFichero(dt.year, dt.month, dt.day)
        #print nombreFichero
        nombreFichero=directorio+nombreFichero
        try:
            d = Dataset(nombreFichero)
            #print 'Hacemos los calculos para el dia '+ str(dt.year)+str(dt.month)+str(dt.day)
            #vamos hacer los calculos para el no2 en el dia indicado
            datosNO2 = d.variables['NO2']
            #solo necesito las primeras 24 horas de la prediccion de la capa 0
            datosNO2Nivel = datosNO2[0:24,0,:,:]
                        #nh=numero de horas, ni=numero de posiciones en el eje x, nj=numero de posiciones en el eje y
            nh,ni,nj=shape(datosNO2Nivel)
                #si es la primera iteraccion inicializo las salidas de datos
            if(ndias==0):
                salidaNO2Media=zeros((ni,nj),dtype=float)
                salidaNO2Sup=zeros((ni,nj),dtype=int)
                #creo una variable intermedia para almacenar las medias del dia
            mediaNO2= np.mean(datosNO2Nivel,axis=0)
            #las sumo a las que ya tengo
            salidaNO2Media+=mediaNO2
                #creo una matriz booleana que indica las posicione que superan el valor de 200
            superaciones = greater(datosNO2Nivel,200)
                #sumo las superaciones para el dia en cada posicion
            superacionesDia = np.sum(superaciones,axis=0)
                #las sumo a las que ya tengo
            salidaNO2Sup+=superacionesDia
            #finalizamos los calculos de no2 para el dia
                #como procesamos un fichero y vamos a tener un fichero por dia incrementamos el dia en 1
            ndias+=1
        except:     
            print "no hay fichero: "+nombreFichero      
    #print ndias
        ##hacemos la media de NO2 para el total de dias 
    salidaNO2Media = salidaNO2Media / ndias
    return salidaNO2Media, salidaNO2Sup

#la salida de la funcion es una matriz con la media anual de nox
def crearMatrizNOx():
    print "Procedemos a crear las matrices de datos de NOx"
        #para almacenar la media anual
    salidaNOxMedia=None

    fechaInicio = date(ano, mes, dia)
    fechaFin = date(ano,12,31)
    ndias=0
    for dt in rrule(DAILY, dtstart=fechaInicio, until=fechaFin):
        nombreFichero = crearNombreFichero(dt.year, dt.month, dt.day)
        print nombreFichero
        nombreFichero=directorio+nombreFichero
        try:
            d = Dataset(nombreFichero)
            print 'Hacemos los calculos para el dia '+ str(dt.year)+str(dt.month)+str(dt.day)
            #vamos hacer los calculos para el nox en el dia indicado
            datosNOx = d.variables['NOX']
            #solo necesito las primeras 24 horas de la prediccion de la capa 0
            datosNOxNivel = datosNOx[0:24,0,:,:]
                        #nh=numero de horas, ni=numero de posiciones en el eje x, nj=numero de posiciones en el eje y
            nh,ni,nj=shape(datosNOxNivel)
                #si es la primera iteraccion inicializo las salidas de datos
            if(ndias==0):
                salidaNOxMedia=zeros((ni,nj),dtype=float)
                
                #creo una variable intermedia para almacenar las medias del dia
            mediaNOx= np.mean(datosNOxNivel,axis=0)
            #las sumo a las que ya tengo
            salidaNOxMedia+=mediaNOx
                
            #finalizamos los calculos de nox para el dia
                #como procesamos un fichero y vamos a tener un fichero por dia incrementamos el dia en 1
            ndias+=1
        except:     
            print "no hay fichero: "+nombreFichero      
    #print ndias
        ##hacemos la media de NO2 para el total de dias 
    salidaNOxMedia = salidaNOxMedia / ndias
    return salidaNOxMedia


#la salida de la funcion son cuatro matrices,Media anual, Superaciones Horarias, Superaciones Diarias, Media Invernal
def crearMatricesSO2():
    print "Procedemos a crear las matrices de SO2"
    #definimos los intervalos para la media invernal
    fechaIntervalo1=datetime(ano,1,1)
    fechaIntervalo2=datetime(ano,3,31)
    fechaIntervalo3=datetime(ano,10,1)
    fechaIntervalo4=datetime(ano,12,31)
        #para almacenar las superaciones horarias de 350
    salidaSO2SupH=None
    
        #para almacenar las superaciones diarias de 125
    salidaSO2SupD=None

    #para almacenar la media de ano civil e invierno
    salidaSO2Inv=None
    
    #para almacenar la media anual
    salidaSO2Med=None

    fechaInicio = date(ano, mes, dia)
    fechaFin = date(ano,12,31)
    ndias=0
    ndiasInvierno=0
    for dt in rrule(DAILY, dtstart=fechaInicio, until=fechaFin):
        nombreFichero = crearNombreFichero(dt.year, dt.month, dt.day)
        #print nombreFichero
        nombreFichero=directorio+nombreFichero
        try:
            d = Dataset(nombreFichero)
            #print 'Hacemos los calculos para el dia '+ str(dt.year)+str(dt.month)+str(dt.day)
            #vamos hacer los calculos para el SO2 en el dia indicado
            datosSO2 = d.variables['SO2']
            #solo necesito las primeras 24 horas de la prediccion de la capa 0
            datosSO2Nivel = datosSO2[0:24,0,:,:]
                        #nh=numero de horas, ni=numero de posiciones en el eje x, nj=numero de posiciones en el eje y
            nh,ni,nj=shape(datosSO2Nivel)
                #si es la primera iteraccion inicializo las salidas de datos
            if(ndias==0):
                salidaSO2SupH = zeros((ni,nj),dtype=int)
                salidaSO2SupD = zeros((ni,nj),dtype=int)
                salidaSO2Inv  = zeros((ni,nj),dtype=float)
                salidaSO2Med  = zeros((ni,nj),dtype=float)
                #creo una variable intermedia para almacenar las medias del dia
            mediaSO2= np.mean(datosSO2Nivel,axis=0)
            #sumo las medias
            salidaSO2Med+=mediaSO2
                        #recojo las superaciones diarias de 125
            superacionesD = greater(mediaSO2,125)
            superacionesDia = np.sum(superacionesD)
                        #las sumo a las que ya tengo            
            salidaSO2SupD+=superacionesDia
            #recojo las superaciones horarias de 350
            superacionesH = greater(datosSO2Nivel,350)
            superacionesHora = np.sum(superacionesH,axis=0)
            #las sumo a las que ya tengo
            salidaSO2SupH+=superacionesHora
            
            #tenemos que comprobar si se puede usar los datos para la media invernal
            if (dt >= fechaIntervalo1 and dt<=fechaIntervalo2) or (dt>=fechaIntervalo3 and dt<=fechaIntervalo4):
                salidaSO2Inv+=mediaSO2
                ndiasInvierno+=1
                
            ndias+=1
            
        except Exception as ins: 
            print ins
            print "no hay fichero: "+nombreFichero      
    #print ndias
        ##hacemos la media de SO2 para el total de dias 
    salidaSO2Inv = salidaSO2Inv / ndiasInvierno
    salidaSO2Med = salidaSO2Med / ndias
    return salidaSO2Med, salidaSO2SupH, salidaSO2SupD,salidaSO2Inv

#la salida de la funcion son dos matrices, Superaciones Diarias y Media Anual
def crearMatricesPM10():
    print "Procedemos a crear las matrices de PM10"
        #para almacenar la media anual
    salidaPM10Media=None

        #para almacenar las superaciones de 50 diario
    salidaPM10Sup=None
    fechaInicio = date(ano, mes, dia)
    fechaFin = date(ano,12,31)
    ndias=0
    for dt in rrule(DAILY, dtstart=fechaInicio, until=fechaFin):
        nombreFichero = crearNombreFichero(dt.year, dt.month, dt.day)
        #print nombreFichero
        nombreFichero=directorio+nombreFichero
        try:
            d = Dataset(nombreFichero)
            #print 'Hacemos los calculos para el dia '+ str(dt.year)+str(dt.month)+str(dt.day)
            #vamos hacer los calculos para el PM10 en el dia indicado
            datosPM10 = d.variables['PM10']
            #solo necesito las primeras 24 horas de la prediccion de la capa 0
            datosPM10Nivel = datosPM10[0:24,0,:,:]
                        #nh=numero de horas, ni=numero de posiciones en el eje x, nj=numero de posiciones en el eje y
            nh,ni,nj=shape(datosPM10Nivel)
                #si es la primera iteraccion inicializo las salidas de datos
            if(ndias==0):
                salidaPM10Media=zeros((ni,nj),dtype=float)
                salidaPM10Sup=zeros((ni,nj),dtype=int)
                #creo una variable intermedia para almacenar las medias del dia
            mediaPM10= np.mean(datosPM10Nivel,axis=0)
            #las sumo a las que ya tengo
            salidaPM10Media+=mediaPM10
                #creo una matriz booleana que indica las posiciones que superan el valor de 50
            superaciones = greater(mediaPM10,50)
                #sumo las superaciones para el dia en cada posicion
            #superacionesDia = np.sum(superaciones,axis=0)
                #las sumo a las que ya tengo
            salidaPM10Sup+=superaciones
            #finalizamos los calculos de no2 para el dia
                #como procesamos un fichero y vamos a tener un fichero por dia incrementamos el dia en 1
            ndias+=1
        except Exception as ins: 
            print ins
            print "no hay fichero"      
    print ndias
        ##hacemos la media de PM10 para el total de dias    
    salidaPM10Media = salidaPM10Media / ndias
    #print salidaPM10Media
    return salidaPM10Media, salidaPM10Sup

#la salida de la funcion en una matriz con la Media Anual
def crearMatricesPM25():
    print "Procedemos a crear las Matrices de PM25"
        #para almacenar la media anual
    salidaPM25Media=None
    salidaPM25Sup=None
    fechaInicio = date(ano, mes, dia)
    fechaFin = date(ano,12,31)
    ndias=0
    for dt in rrule(DAILY, dtstart=fechaInicio, until=fechaFin):
        nombreFichero = crearNombreFichero(dt.year, dt.month, dt.day)
        #print nombreFichero
        nombreFichero=directorio+nombreFichero
        try:
            d = Dataset(nombreFichero)
            #print 'Hacemos los calculos para el dia '+ str(dt.year)+str(dt.month)+str(dt.day)
            #vamos hacer los calculos para el PM25 en el dia indicado
            datosPM25 = d.variables['PM25']
            #solo necesito las primeras 24 horas de la prediccion de la capa 0
            datosPM25Nivel = datosPM25[0:24,0,:,:]
                        #nh=numero de horas, ni=numero de posiciones en el eje x, nj=numero de posiciones en el eje y
            nh,ni,nj=shape(datosPM25Nivel)
                #si es la primera iteraccion inicializo las salidas de datos
            if(ndias==0):
                salidaPM25Media=zeros((ni,nj),dtype=float)
                salidaPM25Sup=zeros((ni,nj),dtype=int)
                #creo una variable intermedia para almacenar las medias del dia
            mediaPM25= np.mean(datosPM25Nivel,axis=0)
            #las sumo a las que ya tengo
            salidaPM25Media+=mediaPM25
                
            ndias+=1
        except Exception as ins: 
            print ins
            print "no hay fichero"      
    #print ndias
        ##hacemos la media de PM25 para el total de dias    
    salidaPM25Media = salidaPM25Media / ndias
    salidaPM25Sup = greater(salidaPM25Media,20)
    #print salidaPM25Media
    return salidaPM25Media, salidaPM25Sup



## Seccion para GIS:

def exportar_geotiff(fichero, origen, tam, resolucion, rot, EPSG, datos ):

    ## Funcion para exportar la batimetria generada a GeoTiff
           
    ##        - origen: Es el origen de la malla
    ##        - tam:    Tamanho en pixeles
    ##        - rot:    Matriz de rotacion como se calcula con utils.rotacion
    ##        - EPSG:   Id del sistema de referencia 
    ##        - datos:  Array a almacenar en el fichero.  '''

    format = "GTiff"
    #format="ESRI Shapefile"
    driver = gdal.GetDriverByName( format )

    #fichero = 'prueba_bat_%im_filt.tiff' % resolucion

    dst_ds = driver.Create(fichero, tam[1], tam[0], 1, gdal.GDT_Float32)
    dst_ds.GetRasterBand(1).WriteArray(datos)

    # top left x, w-e pixel resolution, rotation, top left y, rotation, n-s pixel resolution
    rot = rot*resolucion
    transformacion = [origen[0], rot[0,0], rot[0,1], origen[1], rot[1,0], rot[1,1] ]
    dst_ds.SetGeoTransform( transformacion )

    # set the reference info
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(EPSG)

    dst_ds.SetProjection( srs.ExportToWkt() )

    # Cerramos el fichero (!)
    dst_ds = None

#la salida de la funcion en una matriz con el numero de superaciones del valor maximo octohorario
def crearMatrizCO():
    print "Procedemos a crear la matriz de superaciones de CO"
        #para almacenar las superaciones
    salidaCO=None

    fechaInicio = date(ano, mes, dia)
    fechaFin = date(ano,12,31)
    ndias=0
    for dt in rrule(DAILY, dtstart=fechaInicio, until=fechaFin):
        nombreFichero = crearNombreFichero(dt.year, dt.month, dt.day)
        #print nombreFichero
        nombreFichero=directorio+nombreFichero
        try:
            d2=None
            d = Dataset(nombreFichero)
                    #como son medias octoHorarias las que se necesitan para CO necesito el dia anterior         
            #como puede no existir lo meto en un try
            try:
                dt2 = dt - timedelta(days=1)
                nombreFichero2 = directorio+crearNombreFichero(dt2.year, dt2.month, dt2.day)
                #print "buscamos el dia anterior "+ nombreFichero2
                d2 = Dataset(nombreFichero2)
            except Exception as in2:
                print "No hay fichero dia anterior"

            #print 'Hacemos los calculos para el dia '+ str(dt.year)+str(dt.month)+str(dt.day)
            #vamos hacer los calculos para el CO en el dia indicado
            datosCO = d.variables['CO']
            #solo necesito las primeras 24 horas de la prediccion de la capa 0
            datosCONivel = datosCO[0:24,0,:,:]
            octoDia = zeros((24,55,55),dtype=float)
            datosCOTotal = datosCONivel
            #cojo los ultimos 7 datos del dia anterior, como puede que no existan lo compruebo
            i=0
                    #si no tenemos los datos del dia anterior tenemos que empezar en la posicion 6 por criterios de agregacion
            if d2 is None:
                i=6
            else:
                datosCO2 = d2.variables['CO']
                    #los ultimos 7 datos del dia anterior son los que van de la hora 17 a la 24
                    datosCO2Nivel = datosCO2[17:24,0,:,:]
                    #creo una matriz para usarla para calcular las medias octohorarias del dia
                    datosCOTotal = np.concatenate((datosCO2Nivel,datosCONivel))
            #en la matriz las 8 primeras posiciones seran del dia anterior
            #por lo que hay que empezar en la posicion 8
            j=0
            k=0
            while i < 24:
                while j < 55:
                    while k < 55:
                        #calculamos el datos octohorarios
                        if d2 is None:
                            octoDia[(i-6),j,k] = np.mean(datosCOTotal[i:(i+8),j,k])/1000
                        else:
                            octoDia[i,j,k] = np.mean(datosCOTotal[i:(i+8),j,k])/1000
                        if octoDia[i,j,k]>10.5:
                            print "("+str(i)+","+str(j)+","+str(k)+") = "+str(octoDia[i,j,k])
                        k+=1
                    k=0
                    j+=1
                    ##fin while
                j=0
                i+=1
                ##fin while 
            #nh=numero de horas, ni=numero de posiciones en el eje x, nj=numero de posiciones en el eje y
            nh,ni,nj=shape(octoDia)
                #si es la primera iteraccion inicializo las salidas de datos
            if(ndias==0):
                salidaCO=zeros((ni,nj),dtype=float)
                #creo una variable intermedia para almacenar los maximos del dia
            maxCO= np.max(octoDia,axis=0)
            
            #busco los valore superiores a 10
            
            superaciones = greater(maxCO,10.5)
                
            #sumo las superaciones para el dia en cada posicion
            
            salidaCO+=superaciones
            
                
            ndias+=1
        except Exception as ins:
            #print ins 
            print "no hay fichero del dia"      
    #print ndias
        
    return salidaCO

#la salida de la funcion en una matriz con el numero de superaciones del valor maximo octohorario
def crearMatricesO3():
    
        #para almacenar la media anual
    salidaO3=None
    salidaMax=None
    fechaInicio = date(ano, mes, dia)
    fechaFin = date(ano,12,31)
    ndias=0
    for dt in rrule(DAILY, dtstart=fechaInicio, until=fechaFin):
        nombreFichero = crearNombreFichero(dt.year, dt.month, dt.day)
        #print nombreFichero
        nombreFichero=directorio+nombreFichero
        try:
            d2=None
            d = Dataset(nombreFichero)
                    #como son medias octoHorarias las que se necesitan para O3 necesito el dia anterior         
            #como puede no existir lo meto en un try
                
            try:
                dt2 = dt - timedelta(days=1)
                nombreFichero2 = directorio+crearNombreFichero(dt2.year, dt2.month, dt2.day)
                #print "buscamos el dia anterior "+ nombreFichero2
                d2 = Dataset(nombreFichero2)
            except Exception as in2:
                print "No hay fichero dia anterior"

            #print 'Hacemos los calculos para el dia '+ str(dt.year)+str(dt.month)+str(dt.day)
            #vamos hacer los calculos para el O3 en el dia indicado
            datosO3 = d.variables['O3']
            #solo necesito las primeras 24 horas de la prediccion de la capa 0
            datosO3Nivel = datosO3[0:24,0,:,:]
            octoDia = zeros((24,55,55),dtype=float)
            datosO3Total = datosO3Nivel
            #cojo los ultimos 7 datos del dia anterior, como puede que no existan lo compruebo
            i=0
                    #si no tenemos los datos del dia anterior tenemos que empezar en la posicion 6 por criterios de agregacion
            if d2 is None:
                i=6
            else:
                datosO32 = d2.variables['O3']
                    #los ultimos 7 datos del dia anterior son los que van de la hora 17 a la 24
                    datosO32Nivel = datosO32[17:24,0,:,:]
                    #creo una matriz para usarla para calcular las medias octohorarias del dia
                    datosO3Total = np.concatenate((datosO32Nivel,datosO3Nivel))
            #en la matriz las 8 primeras posiciones seran del dia anterior
            #por lo que hay que empezar en la posicion 8
            j=0
            k=0
            while i < 24:
                while j < 55:
                    while k < 55:
                        #calculamos el datos octohorarios
                        if d2 is None:
                            octoDia[(i-6),j,k] = np.mean(datosO3Total[i:(i+8),j,k])
                        else:
                            octoDia[i,j,k] = np.mean(datosO3Total[i:(i+8),j,k])
                        #if octoDia[i,j,k]>120.5:
                            #print "("+str(i)+","+str(j)+","+str(k)+") = "+str(octoDia[i,j,k])
                        k+=1
                    k=0
                    j+=1
                    ##fin while
                j=0
                i+=1
                ##fin while 
            #nh=numero de horas, ni=numero de posiciones en el eje x, nj=numero de posiciones en el eje y
            nh,ni,nj=shape(octoDia)
                #si es la primera iteraccion inicializo las salidas de datos
            if(ndias==0):
                salidaO3=zeros((ni,nj),dtype=float)
                salidaMax=zeros((ni,nj),dtype=float)
                #creo una variable intermedia para almacenar los maximos del dia
            maxO3= np.max(octoDia,axis=0)
            
            #busco los valore superiores a 10
            
            superaciones = greater(maxO3,120.5)
                
            #sumo las superaciones para el dia en cada posicion
            salidaMax=np.maximum(salidaMax,maxO3)
            salidaO3+=superaciones
            
                
            ndias+=1
        except Exception as ins:
            print ins 
            print "no hay fichero del dia"      
    print ndias
        
    return salidaO3,salidaMax

#la salida de la funcion en una matriz con la AOt40 del anho
def crearMatrizAOT40O3():
        #para almacenar el resultado
    salidaAOt40=None

    fechaInicio = date(ano, 5, 1)
    fechaFin = date(ano,7,31)
    ndias=0
    for dt in rrule(DAILY, dtstart=fechaInicio, until=fechaFin):
        nombreFichero = crearNombreFichero(dt.year, dt.month, dt.day)
        print nombreFichero
        nombreFichero=directorio+nombreFichero
        try:
            d = Dataset(nombreFichero)
            print 'Hacemos los calculos para el dia '+ str(dt.year)+str(dt.month)+str(dt.day)
            #vamos hacer los calculos para el O3 en el dia indicado
            datosO3 = d.variables['O3']
            #en la Aot40 solo se necesita de las 08 a las 20 hora local
            #como los datos estan en utc y solo vamos a trabajar con jornadas de verano
            #le restamos 2 horas al intervalo quedaria de 06 a 22 hora utc
            
            datosO3Nivel = datosO3[6:22,0,:,:]
            aotDia = zeros((16,55,55),dtype=float)
            i=0
            j=0
            k=0
            while i < 16:
                while j < 55:
                    while k < 55:
                        #miramos si los valores son mayores a 80
                        if datosO3Nivel[i,j,k] > 80:
                            aotDia[i,j,k] = (datosO3Nivel[i,j,k]) - 80
                        else:
                            aotDia[i,j,k] = 0
                        k+=1
                    k=0
                    j+=1
                    ##fin while
                j=0
                i+=1
                ##fin while

                        #nh=numero de horas, ni=numero de posiciones en el eje x, nj=numero de posiciones en el eje y
            nh,ni,nj=shape(datosO3Nivel)
                #si es la primera iteraccion inicializo las salidas de datos
            if(ndias==0):
                salidaAOt40=zeros((ni,nj),dtype=float)
                
                #creo una variable intermedia para almacenar las sumas del dia
            
            salidaDia= np.sum(aotDia,axis=0)
            #las sumo a las que ya tengo
            salidaAOt40+=salidaDia
                
            ndias+=1
        except Exception as ins: 
            print ins
            print "no hay fichero"      
    #print ndias
        
    return salidaAOt40


#funcion que a partir de un netcdf y una matriz(x,y) grafica la matriz. 
#Esta preparado para medias
def mostrarGrafica(d, datos):
    #cojo las longitudes y latitudes del netcdf
    
    lons = d.variables['lon'][:]
    lats = d.variables['lat'][:]

    #calculo la longitud y latitud centro del mapa
    lon_0 = lons.mean()
    lat_0 = lats.mean()
    #creo el Basemap en projeccion rectangular

    m=Basemap(resolution='h', projection='merc', lat_ts=40, llcrnrlat=lats[0,0],llcrnrlon=lons[0,0], urcrnrlon=lons[-1,-1], urcrnrlat=lats[-1,-1])
    #m=Basemap(resolution='h', projection='merc', lat_ts=40, llcrnrlat=xysila,llcrnrlon=xysilo, urcrnrlon=xyidlo, urcrnrlat=xyidla)
    x,y = m(lons,lats)

    cs = m.pcolor(x,y,np.squeeze(datos))
    m.drawcoastlines()
    
    m.drawstates()
    
    m.drawcountries()   


    cbar = m.colorbar(cs,location='bottom', pad='10%')

    plt.show()

#Funcion que crea un geoTiff para los datos indicados
def crearGraficaDatos(d,datos, nombreFichero):

    #d->Netcdf descriptivo
    #datos->Matriz de 2 dimensiones con los datos a graficar
    #nombreFichero-> Nombre que le queremos dar al fichero de salida

    #Como queremos eliminar de la imagen los puntos del mar cargamos el dataset con los datos de suelo
    landshapeDataset = Dataset(ficherosuelo)
    #la variable ocean contiene el porcentaje de 0 a 1 de la superficie de mar de la celda
    capa = landshapeDataset.variables['Ocean'][:,:] 
    #vamos a recorrer todos los puntos de la malla para eliminar los datos de puntos que sean completamente mar
    i=0
    while i < datos.shape[0] :      
        j=0
        while j < datos.shape[1]:
            if capa[i,j]==1:
                datos[i,j]=-9999
            #str(datos[i,j]).replace(',', '.')
            j+=1
        i+=1
    #Creamos el nombre del fichero  
    fichero = 'salidas/'+nombreFichero+'.tiff'

    #para recortar la malla y centrar el mapa en galicia empezamos en la posicion x=10 e y=10
    #tambien se definen los finales del eje x e y
    iniciox = 10
    finy = -14
    finx=-1
    inicioy=10

    #recogemos las coordenadas del punto de inicio
    lonOrigen = [d.variables["lon"][inicioy,iniciox]]
    latOrigen = [d.variables["lat"][inicioy,iniciox]]
    
    #obtenemos de la matriz de datos los encuadrados entre la posiciones de inicio y final
    datosImagen = datos[inicioy:finy,iniciox:finx]
    origen = np.array([lonOrigen,latOrigen])
    #tamanho de los datos 
    tam        = datosImagen.shape
    #tamanho de la celda
        resolucion = 0.069999695 #tamanho de la celda

    rot        = np.eye(2) #crea una matriz vacia de 2x2
    #creamos el geotiff
    exportar_geotiff(fichero, origen, tam, resolucion, rot, 4326, datosImagen)

#funcion que crea todas las matrices de datos que se usaran para crear las graficas
def cargarMatrices():
#dataset con las medias de NO2
    global matrizNO2Media
#dataset con las superaciones de NO2
    global matrizNO2Sup
#dataset con las medias de NOx
    global matrizNOxMedia
#dataset con las medias de PM10
    global matrizPM10Media
#dataset con las superaciones de PM10
    global matrizPM10Sup
#dataset con los valores de la AOT40 de O3
    global matrizAOT40
#dataset con los maximos octohorarios de O3
    global matrizMaxO3
#dataset con las superaciones de maximos octohorarios de O3
    global matrizSupMaxO3
#dataset con las superaciones de CO
    global matrizCO
#dataset con las medias de PM25
    global matrizPM25Media
#dataset con las superaciones de la media de PM25
    global matrizPM25Sup
#dataset con la media de SO2
    global matrizSO2Med
#dataset con las superaciones horarias de SO2
    global matrizSO2SupH
#dataset con las superaciones diarias de SO2
    global matrizSO2SupD
#dataset con las medias invernales de SO2
    global matrizSO2Inv
    matrizNO2Media,matrizNO2Sup=crearMatricesNO2()
    matrizNOxMedia=crearMatrizNOx()
    matrizPM10Media,matrizPM10Sup=crearMatricesPM10()
    matrizAOT40=crearMatrizAOT40O3()
    matrizSupMaxO3,matrizMaxO3=crearMatricesO3()
    matrizCO=crearMatrizCO()
    matrizPM25Media, matrizPM25Sup=crearMatricesPM25()
    matrizSO2Med,matrizSO2SupH,matrizSO2SupD,matrizSO2Inv=crearMatricesSO2()

def cargarGraficas():
    crearGraficaDatos(f,matrizNO2Media,"NO2Media2017")
    crearGraficaDatos(f,matrizNO2Sup,"NO2Sup2017")
    crearGraficaDatos(f,matrizPM10Media,"PM10Media2017")
    crearGraficaDatos(f,matrizPM10Sup,"PM10Sup2017")
    crearGraficaDatos(f,matrizAOT40,"AOT402017")
    crearGraficaDatos(f,matrizSupMaxO3,"O3SupMaxOH2017")
    crearGraficaDatos(f,matrizMaxO3,"O3MaxOD2017")
    crearGraficaDatos(f,matrizCO,"CO2017")
    crearGraficaDatos(f,matrizPM25Media,"PM25Media2017")
    crearGraficaDatos(f,matrizSO2Med,"SO2Media2017")
    crearGraficaDatos(f,matrizSO2SupH,"SO2SupH2017")
    crearGraficaDatos(f,matrizSO2SupD,"SO2SupD2017")

    crearGraficaDatos(f,matrizSO2Inv,"SO2Inv2017")
    crearGraficaDatos(f,matrizPM25Sup,"PM25Sup2017")
    crearGraficaDatos(f,matrizNOxMedia,"NOxMedia2017")

def cargarGraficasPorPasos():
#dataset con las medias de NO2
    global matrizNO2Media
#dataset con las superaciones de NO2
    global matrizNO2Sup
#dataset con las medias de NOx
    global matrizNOxMedia
#dataset con las medias de PM10
    global matrizPM10Media
#dataset con las superaciones de PM10
    global matrizPM10Sup
#dataset con los valores de la AOT40 de O3
    global matrizAOT40
#dataset con los maximos octohorarios de O3
    global matrizMaxO3
#dataset con las superaciones de maximos octohorarios de O3
    global matrizSupMaxO3
#dataset con las superaciones de CO
    global matrizCO
#dataset con las medias de PM25
    global matrizPM25Media
#dataset con las superaciones de la media de PM25
    global matrizPM25Sup
#dataset con la media de SO2
    global matrizSO2Med
#dataset con las superaciones horarias de SO2
    global matrizSO2SupH
#dataset con las superaciones diarias de SO2
    global matrizSO2SupD
#dataset con las medias invernales de SO2
    global matrizSO2Inv
    matrizNO2Media,matrizNO2Sup=crearMatricesNO2()
    crearGraficaDatos(f,matrizNO2Media,"NO2Media2017")
    crearGraficaDatos(f,matrizNO2Sup,"NO2Sup2017")
    del matrizNO2Media
    del matrizNO2Sup
    matrizNOxMedia=crearMatrizNOx()
    crearGraficaDatos(f,matrizNOxMedia,"NOXMedia2017")
    del matrizNOxMedia
    matrizPM10Media,matrizPM10Sup=crearMatricesPM10()
    crearGraficaDatos(f,matrizPM10Media,"PM10Media2017")
    crearGraficaDatos(f,matrizPM10Sup,"PM10Sup2017")
    del matrizPM10Media
    del matrizPM10Sup
    matrizAOT40=crearMatrizAOT40O3()
    crearGraficaDatos(f,matrizAOT40,"AOT402017")
    del matrizAOT40
    matrizSupMaxO3,matrizMaxO3=crearMatricesO3()
    crearGraficaDatos(f,matrizMaxO3,"O3MaxOD2017")
    crearGraficaDatos(f,matrizSupMaxO3,"O3SupMaxO2017")
    del matrizMaxO3
    del matrizSupMaxO3
    matrizCO=crearMatrizCO()
    crearGraficaDatos(f,matrizCO,"CO2017")
    del matrizCO
    matrizPM25Media, matrizPM25Sup=crearMatricesPM25()
    crearGraficaDatos(f,matrizPM25Media,"PM25Media2017")
    crearGraficaDatos(f,matrizPM25Sup,"PM25Sup2017")
    del matrizPM25Media
    del matrizPM25Sup
    matrizSO2Med,matrizSO2SupH,matrizSO2SupD,matrizSO2Inv=crearMatricesSO2()
    crearGraficaDatos(f,matrizSO2Med,"SO2Media2017")
    crearGraficaDatos(f,matrizSO2SupH,"SO2SupH2017")
    crearGraficaDatos(f,matrizSO2SupD,"SO2SupD2017")
    crearGraficaDatos(f,matrizSO2Inv,"SO2Inv2017")
    del matrizSO2Med
    del matrizSO2SupH
    del matrizSO2SupD
    del matrizSO2Inv
    subprocess.call(['./subirImagenes.sh'])
    
def realizarTodo():
    cargarMatrices()
    cargarGraficas()
    subprocess.call(['./subirImagenes.sh'])
