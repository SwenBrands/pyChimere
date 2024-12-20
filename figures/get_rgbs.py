# -*- coding: utf-8 -*-

#define the colorbar
BLUE = '#6699cc'
GRAY = '#999999'

rgbs_o3 = [
 '#e0ffff', #azul, casi blanco
 '#afeeee', #azul
 '#00ffff', #azul
 '#00ced1', #azul
 '#32cd32', #verde
 '#7fff00', #verde
 '#ffff00', #amarillo de 200 a 350 en pasos de 50
 '#ffd700', #amarillo
 '#ffa500', #amarillo
 '#ff6347', #rojo de 350 a 500 en pasos de 50
 '#ff4500', #rojo
 '#ff0000', #rojo
 '#a52a2a'] #marrón
under_o3 = '#e0ffff'
over_o3 = '#a52a2a'
cbounds_o3 = [0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260]

rgbs_no2 = [
 '#e0ffff', #azul, casi blanco
 '#afeeee', #azul
 '#00ffff', #azul
 '#00ced1', #azul
 '#66cdaa', #verde de 100 a 200 en pasos de 25
 '#32cd32', #verde
 '#7fff00', #verde
 '#d9fc04', #amarillo
 '#f1fc04', #amarillo
 '#fcea04', #amarillo
 '#fcc604', #amarillo
 '#fca104', #amarillo
 '#ff7f50', #rojo de 50 a 100 en pasos de 10
 '#f08080', #rojo
 '#ff6347', #rojo
 '#ff4500', #rojo
 '#ff0000', #rojo
 '#a52a2a'] #marrón
under_no2 = '#e0ffff'
over_no2 = '#a52a2a'
cbounds_no2 = [0, 10, 20, 30, 40, 60, 80, 100, 120, 140, 160, 180, 200, 240, 280, 320, 360, 400, 440]

rgbs_so2 = [
 '#e0ffff', #azul, casi blanco
 '#afeeee', #azul
 '#00ffff', #azul
 '#00ced1', #azul
 '#66cdaa', #verde de 100 a 200 en pasos de 25
 '#20b2aa', #verde
 '#32cd32', #verde
 '#7fff00', #verde
 '#ffff00', #amarillo de 200 a 350 en pasos de 50
 '#ffd700', #amarillo
 '#ffa500', #amarillo
 '#ff6347', #rojo de 350 a 500 en pasos de 50
 '#ff4500', #rojo
 '#ff0000', #rojo
 '#a52a2a'] #marrón
under_so2 = '#e0ffff'
over_so2 = '#a52a2a'
cbounds_so2 = [0, 25, 50, 75, 100, 125, 150, 175, 200, 250, 300, 350, 400, 450, 500, 550]

rgbs_co = [
 '#e0ffff', #azul, casi blanco
 '#00ffff', #azul
 '#00ced1', #azul
 '#32cd32', #verde
 '#7fff00', #verde
 '#ffff00', #amarillo de 200 a 350 en pasos de 50
 '#ffd700', #amarillo
 '#ffa500', #amarillo
 '#ff6347', #rojo de 350 a 500 en pasos de 50
 '#ff0000', #rojo
 '#a52a2a'] #marrón
under_co = '#e0ffff'
over_co = '#a52a2a'
cbounds_co = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]

rgbs_pm10 = [
 '#e0ffff', #azul, casi blanco
 '#afeeee', #azul
 '#00ffff', #azul
 '#00ced1', #azul
 '#66cdaa', #verde de 100 a 200 en pasos de 25
 '#32cd32', #verde
 '#7fff00', #verde
 '#ffff00', #amarillo de 200 a 350 en pasos de 50
 '#ffd700', #amarillo
 '#ffa500', #amarillo
 '#ff7f50', #rojo de 50 a 100 en pasos de 10
 '#f08080', #rojo
 '#ff6347', #rojo
 '#ff4500', #rojo
 '#ff0000', #rojo
 '#a52a2a'] #marrón
under_pm10 = '#e0ffff'
over_pm10 = '#a52a2a'
cbounds_pm10 = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 90, 100, 110]

rgbs_pm25 = [
 '#e0ffff', #azul, casi blanco
 '#afeeee', #azul
 '#00ffff', #azul
 '#00ced1', #azul
 '#66cdaa', #verde de 100 a 200 en pasos de 25
 '#20b2aa', #verde
 '#32cd32', #verde
 '#7fff00', #verde
 '#ffff00', #amarillo
 '#ffa500', #amarillo
 '#ff7f50', #rojo de 50 a 100 en pasos de 10
 '#f08080', #rojo
 '#ff6347', #rojo
 '#ff4500', #rojo
 '#ff0000', #rojo
 '#a52a2a'] #marrón
under_pm25 = '#e0ffff'
over_pm25 = '#a52a2a'
cbounds_pm25 = [0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20, 22.5, 25, 30, 35, 40, 45, 50, 55]

rgbs_ica = [
 '#00ffff', #azul
 '#3cb371', #verde
 '#ffd700', #amarillo
 '#ff6347', #rojo tomate
 '#a52a2a'] #marrón
under_ica = '#00ffff'
over_ica = '#a52a2a'
#cbounds_ica = [-0.5, 0.5, 1.5, 2.5, 3.5, 4.5]
cbounds_ica = [0, 1, 2, 3, 4, 5]

rgbs_oldica = [
 '#009900', #verde ica
 '#FFFF00', #amarillo ica
 '#FF0000', #rojo ica
 '#9900CC'] #violeta ica
under_oldica = '#009900'
over_oldica = '#9900CC'
cbounds_oldica = [0, 1, 2, 3, 4]

#create cbounds for bcar and hght using color options from pm10 and no2 in plot_results.py
cbounds_pbcar = [0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5]
cbounds_hght = [0, 20, 40, 60, 80, 100, 150, 200, 250, 300, 400, 500, 600, 700, 800, 900, 1000, 1200, 1400]

rgbs_stand = [
 '#3200fe', #azul
 '#0032fe', #azul
 '#0096fe', #azul
 '#00e6fe', #azul
 '#0ef2ee', #azul
 '#00e677', #verde
 '#00e650', #verde 
 '#00fa00', #verde
 '#c5ed12', #amarillo
 '#fee100', #amarillo
 '#feae00', #naranja
 '#e67d00', #naranja
 '#e66400', #naranja
 '#c8321d', #marron
 '#aa001d', #marron
 '#c80064', #violeta
 '#f20c86', #violeta
 '#fa00fe', #violeta
 '#9600fe', #violeta
 '#960096'] #violeta
under_stand = '#320032'
over_stand = '#960096'

 
