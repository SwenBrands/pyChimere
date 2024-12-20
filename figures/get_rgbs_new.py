# -*- coding: utf-8 -*-

from utiles import crea_cmap
import numpy as np
execfile('rgbs2string.py')

rgbs_o3 = np.array([[238, 254, 253], [185, 254, 253], [136, 254, 252], [0, 255, 250], \
    [104, 234, 206], [80, 200, 163], \
    [245, 245, 118], [238, 238, 97], [255, 255, 0], \
    [255, 194, 79], [255, 148, 79], [255, 79, 92], \
    [192, 0, 0], [114, 2, 2]])/np.float32(255)
rgbs_o3 = rgbs2string(rgbs_o3)
hex_o3 = ['#EEFEFD', '#B9FEFD', '#88FEFC', '#00FFFA', \
    '#68EACE', '#50C8A3', \
    '#F5F576', '#EEEE61', '#FFFF00', \
    '#FFC24F', '#FF944F', '#FF4F5C', \
    '#C00000', '#720202']
cbounds_o3 = np.array([0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260], np.float32)


rgbs_no2 = np.array([[238, 254, 253], [185, 254, 253], [136, 254, 252], [0, 255, 250], \
    [104, 234, 206], [89, 214, 175], [80, 200, 163], \
    [255, 255, 204], [245, 245, 118], [238, 238, 97], [230, 230, 6], [255, 255, 0], \
    [255, 194, 79], [255, 175, 79], [255, 148, 79], [255, 122, 79], [255, 79, 92], \
    [192, 0, 0], [114, 2, 2]])/np.float32(255)
rgbs_no2 = rgbs2string(rgbs_no2)
hex_no2 = ['#EEFEFD', '#B9FEFD', '#88FEFC', '#00FFFA', \
    '#68EACE', '#59D6AF', '#50C8A3', \
    '#FFFFCC', '#F5F576', '#EEEE61', '#E6E606', '#FFFF00', \
    '#FFC24F', '#FFAF4F', '#FF944F', '#FF7A4F', '#FF4F5C', \
    '#C00000', '#720202']
cbounds_no2 = np.array([0, 10, 20, 30, 40, 60, 80, 100, 120, 140, 160, 180, 200, 240, 280, 320, 360, 400, 440], np.float32)


rgbs_so2 = np.array([[238, 254, 253], [185, 254, 253], [136, 254, 252], [0, 255, 250], \
    [104, 234, 206], [97, 222, 184], [89, 214, 175], [80, 200, 163], \
    [245, 245, 118], [238, 238, 97], [255, 255, 0], \
    [255, 194, 79], [255, 148, 79], [255, 79, 92], \
    [192, 0, 0], [114, 2, 2]])/np.float32(255)
rgbs_so2 = rgbs2string(rgbs_so2)
hex_so2 = ['#EEFEFD', '#B9FEFD', '#88FEFC', '#00FFFA', \
    '#68EACE', '#61DEB8', '#59D6AF', '#50C8A3', \
    '#F5F576', '#EEEE61', '#FFFF00', \
    '#FFC24F', '#FF944F', '#FF4F5C', \
    '#C00000', '#720202']
cbounds_so2 = np.array([0, 25, 50, 75, 100, 125, 150, 175, 200, 250, 300, 350, 400, 450, 500, 550], np.float32)


rgbs_co = np.array([[238, 254, 253], [136, 254, 252], [0, 255, 250], \
    [104, 234, 206], [80, 200, 163], \
    [245, 245, 118], [238, 238, 97], [255, 255, 0], \
    [255, 148, 79], \
    [192, 0, 0], [114, 2, 2]])/np.float32(255)
rgbs_co = rgbs2string(rgbs_co)
hex_co = ['#EEFEFD', '#88FEFC', '#00FFFA', \
    '#68EACE', '#50C8A3', \
    '#F5F576', '#EEEE61', '#FFFF00', \
    '#FF944F', \
    '#C00000', '#720202']
cbounds_co = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11], np.float32)


rgbs_pm10 = np.array([[238, 254, 253], [105, 254, 253], [136, 254, 252], [0, 255, 250], \
    [104, 234, 206], [97, 222, 184], [89, 214, 175], \
    [245, 245, 118], [238, 238, 97], [255, 255, 0], \
    [255, 194, 79], [255, 175, 79], [255, 148, 79], [255, 122, 79], [255, 79, 92], \
    [192, 0, 0], [114, 2, 2]])/np.float32(255)
rgbs_pm10 = rgbs2string(rgbs_pm10)
hex_pm10 = ['#EEFEFD', '#88FEFC', '#69FEFD', '#00FFFA', \
    '#68EACE', '#61DEB8', '#59D6AF', \
    '#F5F576', '#EEEE61', '#FFFF00', \
    '#FFC24F', '#FFAF4F', '#FF944F', '#FF7A4F', '#FF4F5C', \
    '#C00000', '#720202']
cbounds_pm10 = np.array([0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 90, 100, 110], np.float32)

rgbs_pm10ant = rgbs_pm10
rgbs_pm10bio = rgbs_pm10
rgbs_pdust = rgbs_pm10
rgbs_psalt = rgbs_pm10
cbounds_pm10ant = cbounds_pm10
cbounds_pm10bio = cbounds_pm10
cbounds_pdust = cbounds_pm10
cbounds_psalt = cbounds_pm10


rgbs_pm25 = np.array([[238, 254, 253],[105, 254, 253], [136, 254, 252], [0 ,255, 250], \
    [104, 234, 206], [97, 222, 184], [89, 214, 175], [80, 200, 163], \
    [245, 245, 118], [255, 255, 0], \
    [255, 194, 79], [255, 175, 79], [255, 148, 79], [255, 122, 79], [255, 79, 92], \
    [192, 0, 0], [114, 2, 2]])/np.float32(255)
rgbs_pm25 = rgbs2string(rgbs_pm25)
hex_pm25 = ['#EEFEFD', '#88FEFC', '#69FEFD', '#00FFFA', \
    '#68EACE', '#61DEB8', '#59D6AF', '#50C8A3', \
    '#F5F576', '#FFFF00', \
    '#FFC24F', '#FFAF4F', '#FF944F', '#FF7A4F', '#FF4F5C', \
    '#C00000', '#720202']
cbounds_pm25 = np.array([0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20, 22.5, 25, 30, 35, 40, 45, 50, 55], np.float32)

rgbs_pm25ant = rgbs_pm25
rgbs_pm25bio = rgbs_pm25
cbounds_pm25ant = cbounds_pm25
cbounds_pm25bio = cbounds_pm25

rgbs_ica = np.array([[0, 255, 255], [60, 179, 113], [250, 215, 0], [255, 99, 71], [165, 42, 42]])/np.float32(255)
rgbs_ica = rgbs2string(rgbs_ica)
hex_ica = ['#00fffa', '#50c8a3', '#ffff00', '#ff4f5c', '#C00000']
hex_under_ica = '#00fffa'
hex_over_ica = '#C00000'
cbounds_ica = np.array([0, 1, 2, 3, 4, 5], np.float32)

#create cbounds for bcar and hght using color options from pm10 and no2 in plot_results.py
rgbs_pbcar = rgbs_pm10
cbounds_pbcar = np.array([0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8], np.float32)
rgbs_hght = rgbs_no2
cbounds_hght = np.array([0, 20, 40, 60, 80, 100, 150, 200, 250, 300, 400, 500, 600, 700, 800, 900, 1000, 1200, 1400], np.float32)

#rgbs_stand = rgbs_no2

#define cbounds for 8h variables
cbounds_o3_8h = cbounds_o3
cbounds_no2_8h = cbounds_no2
cbounds_so2_8h = cbounds_so2
cbounds_co_8h = cbounds_co
cbounds_pm10_8h = cbounds_pm10
cbounds_pm10ant_8h = cbounds_pm10
cbounds_pm10bio_8h = cbounds_pm10
cbounds_pm25_8h = cbounds_pm25
cbounds_pm25ant_8h = cbounds_pm25
cbounds_pm25bio_8h = cbounds_pm25
cbounds_hght_8h = cbounds_hght
cbounds_pbcar_8h = cbounds_pbcar
cbounds_pdust_8h = cbounds_pdust
cbounds_psalt_8h = cbounds_psalt

#define rgbs for the 8h variables
rgbs_o3_8h = rgbs_o3
rgbs_no2_8h = rgbs_no2
rgbs_so2_8h = rgbs_so2
rgbs_co_8h = rgbs_co
rgbs_pm10_8h = rgbs_pm10
rgbs_pm10ant_8h = rgbs_pm10
rgbs_pm10bio_8h = rgbs_pm10
rgbs_pm25_8h = rgbs_pm25
rgbs_pm25ant_8h = rgbs_pm25
rgbs_pm25bio_8h = rgbs_pm25
rgbs_hght_8h = rgbs_hght
rgbs_pbcar_8h = rgbs_pbcar
rgbs_pdust_8h = rgbs_pdust
rgbs_psalt_8h = rgbs_psalt

### HEREAFTER, colorbars are created, using hex code, implement this in future versions of the script!!

###create the colormaps for hourly variables and ica, cb = colorbar
#cb_o3 = crea_cmap(cbounds_o3, rgbs_o3, rgbs_o3[0], rgbs_o3[-1]) #input arguments are cbounds, under and over
#cb_no2 = crea_cmap(cbounds_no2, rgbs_no2, rgbs_no2[0], rgbs_no2[-1])
#cb_so2 = crea_cmap(cbounds_so2, rgbs_so2, rgbs_so2[0], rgbs_so2[-1])
#cb_co = crea_cmap(cbounds_co, rgbs_co, rgbs_co[0], rgbs_co[-1])
#cb_pm10 = crea_cmap(cbounds_pm10, rgbs_pm10, rgbs_pm10[0], rgbs_pm10[-1])
#cb_pm10ant = crea_cmap(cbounds_pm10ant, rgbs_pm10ant, rgbs_pm10ant[0], rgbs_pm10ant[-1])
#cb_pm10bio = crea_cmap(cbounds_pm10bio, rgbs_pm10bio, rgbs_pm10bio[0], rgbs_pm10bio[-1])
#cb_pm25 = crea_cmap(cbounds_pm25, rgbs_pm25, rgbs_pm25[0], rgbs_pm25[-1])
#cb_pm25ant = crea_cmap(cbounds_pm25ant, rgbs_pm25ant, rgbs_pm25ant[0], rgbs_pm25ant[-1])
#cb_pm25bio = crea_cmap(cbounds_pm25bio, rgbs_pm25bio, rgbs_pm25bio[0], rgbs_pm25bio[-1])
#cb_hght = crea_cmap(cbounds_hght, rgbs_hght, rgbs_hght[0], rgbs_hght[-1])
#cb_pdust = crea_cmap(cbounds_pdust, rgbs_pdust, rgbs_pdust[0], rgbs_pdust[-1])
#cb_psalt = crea_cmap(cbounds_psalt, rgbs_psalt, rgbs_psalt[0], rgbs_psalt[-1])
#cb_pbcar = crea_cmap(cbounds_pbcar, rgbs_pbcar, rgbs_pbcar[0], rgbs_pbcar[-1])
#cb_ica = crea_cmap(cbounds_ica, rgbs_ica, rgbs_ica[0], rgbs_ica[-1])

###create the colormaps for the 8h variables
#cb_o3_8h = cb_o3
#cb_no2_8h = cb_no2
#cb_so2_8h = cb_so2
#cb_co_8h = cb_co
#cb_pm10_8h = cb_pm10
#cb_pm10ant_8h = cb_pm10ant
#cb_pm10bio_8h = cb_pm10bio
#cb_pm25_8h = cb_pm25
#cb_pm25ant_8h = cb_pm25ant
#cb_pm25bio_8h = cb_pm25bio
#cb_hght_8h = cb_hght
#cb_pdust_8h = cb_pdust
#cb_psalt_8h = cb_psalt
#cb_pbcar_8h = cb_pbcar
