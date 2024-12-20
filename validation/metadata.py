# -*- coding: utf-8 -*-

import numpy as np

#station locations
station_loc = [42.876011,-8.559122,
      42.431042,-8.068685,
      42.219002,-8.742058,
      43.480659,-8.241543,
      42.061626,-7.408726,
      42.997922,-7.550728,
      42.353033,-7.877667,
      42.181389,-8.523689,
      42.420431,-8.664424,
      43.367091,-8.420599,
      42.887759,-8.531091,
      43.382786,-8.409211,
      43.668400,-7.380212,
      43.682562,-7.433673,
      43.685662,-7.507667,
      43.335241,-8.472004,
      43.300000,-8.523333,
      43.300000,-8.465278,
      43.259747,-8.554653,
      43.334827,-8.467230,
      43.306398,-8.583548,
      43.324167,-8.501944,
      43.341111,-8.483333,
      43.306944,-8.491667,
      43.302295,-8.505986,
      43.335241,-8.472004,
      43.354167,-8.424722,
      43.354167,-8.424722,
      43.491634,-8.252325,
      43.405399,-7.988619,
      43.536111,-7.740278,
      43.445278,-7.918056,
      43.450000,-7.847222,
      43.470833,-7.840000,
      43.312778,-7.693238,
      43.185833,-8.468611,
      43.173889,-8.188333,
      43.235000,-8.317222,
      43.095000,-8.494722,
      42.726111,-7.451950,
      42.713439,-7.451197,
      42.210278,-8.736944,
      42.202778,-8.746944,
      42.404722,-8.670278,
      42.434167,-8.673611,
      42.95357565,-9.1861455,
      43.008889,-9.118861,
      42.906389,-8.506667,
      43.116389,-8.352222,
      43.256389,-8.511667,
      43.203611,-8.413889,
      43.147500,-8.545278,
      43.354167,-8.424722,
      43.519800,-8.150118,
      42.060486,-7.730241,
      42.656080,-8.107217,
      42.420431,-8.664424,
      42.222428,-8.711082,
      43.377452,-8.436899,
      42.425886,-8.644014,
      42.0,-8.5]
      
#die letzten zwei sind Fake, weil Anthony versehentlich Lanas und Sorrizo in die csv Dateien eingefuegt hat, wo es keine Daten gibt

#seperate lats and lons
latid = range(0,len(station_loc),2)
lonid = range(1,len(station_loc),2)
lat = [station_loc[ii] for ii in latid]
lon = [station_loc[ii] for ii in lonid]

#assign station names
station_names_meta = ['Santiago_Campus','Carballino','Vigo_Coia','Ferrol','Laza',
    'Lugo','Ourense','Ponteareas','Pontevedra_Mollabao','A_Coruna_Riazor','Santiago_San_Caetano',
    'A_Coruna_Torre_Hercules','Burela','Rio_Cobo','Xove','Pastoriza_A','Armenton',
    'Bordeiras','Paiosaco','Plaza_Pastoriza','Sorrizo','Sabon','Suevos','Arteixo',
    'Centro_Civico','Pastoriza','A_Grela_S','A_Grela_C','A_Cabana','Fraga_Redonda','Louseiras',
    'Macineira','Magdalena','Marraxon','Mourence','Cerceda','Paraxon','San_Vicente_de_Vigo',
    'Villagudin','NNW','Sur','Este_Estacion_1','Oeste_Estacion_2','Areeiro','Campelo',
    'Cee','Dumbria','Campo_de_Futbol','Buscas','Cendon','MonteXalo','Rodis','SGL_Carbon','Xubia',
    'Xinzo_de_Limia','Lalin','Pontevedra','Lope_de_Vega','San_Pedro','Pontevedra_Campolongo','Lanas'] #Noia, O_Savinao y Lanas are FAKE

#T=traffic, B=Background, I=Industry, S=sulphur dioxide influence, C=Coast
exposure = ['T','-','T','T','B',
            'T','T','B','-','TS','T',
            'B','B','I','B','IS','-',
            '-','B','-','-','IS','-','-',
            'IS','IS','IS','IS','T','B','B',
            'B','B','B','B','B','B','B',
            'B','IS','IS','T','T','B','T',
            'B','B','I','B','-','-','B','I','T',
            'T','T','-','T','TS','T','-']

##setting exposure entry for these station to '-' will exclude this station in run_validation_seaborn.py, currently centro civico is an outlier for NO2 on the gal0504r domains, check out why!
#exposure[24] = '-'

#bring lat and lon into the same order than the order of the stations in the csv files
#station_names_csv = ['Santiago_Campus','Carballino','Vigo_Coia','Ferrol','Laza','Lugo','Ourense','Ponteareas','Pontevedra','A_Coruna_Riazor','Santiago_San_Caetano','A_Coruna_Torre_Hercules','Burela','Rio_Cobo','Xove','Pastoriza_A','Lanas','Paiosaco','Sorrizo','Centro_Civico','Pastoriza','A_Grela_S','A_Cabana','Fraga_Redonda','Louseiras','Macineira','Magdalena','Marraxon','Mourence','Cerceda','Paraxon','San_Vicente_de_Vigo','Villagudin','NNW','Sur','Este_Estacion_1','Oeste_Estacion_2','Areeiro','Campelo','Cee','Dumbria','Campo_de_Futbol','Buscas','Cendon','MonteXalo','Rodis','SGL_Carbon','Xubia','Pontevedra_Campolongo','Lope_de_Vega','Xinzo_de_Limia','Lalin','San_Pedro','A_Grela_S']
#station_names_csv = ['Santiago_Campus','Carballino','Vigo_Coia','Ferrol','Laza','Lugo','Ourense','Ponteareas','Pontevedra','A_Coruna_Riazor','Santiago_San_Caetano','A_Coruna_Torre_Hercules','Burela','Rio_Cobo','Xove','Pastoriza_A','Armenton','Bordeiras','Lanas','Paiosaco','Plaza_Pastoriza','Sorrizo','Sabon','Suevos','Arteixo','Centro_Civico','Pastoriza','A_Grela_S','A_Cabana','Fraga_Redonda','Louseiras','Macineira','Magdalena','Marraxon','Mourence','Cerceda','Paraxon','San_Vicente_de_Vigo','Villagudin','NNW','Sur','Este_Estacion_1','Oeste_Estacion_2','Areeiro','Campelo','Cee','Dumbria','Campo_de_Futbol','Buscas','Cendon','MonteXalo','Rodis','SGL_Carbon','Xubia','Pontevedra_Campolongo','Lope_de_Vega','Xinzo_de_Limia','Lalin','San_Pedro','A_Grela_C']

#Santiago - Campus;Vigo - Coia;Ferrol;Laza;Lugo;Ourense;Ponteareas;A Coruña - Riazor;Santiago - San Caetano;A Coruña - Torre Hércules;Burela;Río Cobo;Xove;Pastoriza (A);Paiosaco ;Sabón;Centro cívico;Pastoriza ;A Grela (S);A Cabana ;Fraga Redonda;Louseiras;Maciñeira;Magdalena;Marraxón;Mourence;Cerceda;Paraxón;San Vicente de Vigo;Villagudín;NNW;Sur;Este - Estacion 1;Oeste - Estacion 2;Areeiro;Campelo;Cee;Dumbría;Campo de Fútbol ;Buscás;Rodís;SGL Carbon;Xubia;Pontevedra - Campolongo;Vigo - Lope de Vega;Xinzo de Limia;Lalín;San Pedro;A Grela C;
station_names_csv = ['Santiago_Campus','Vigo_Coia','Ferrol','Laza','Lugo','Ourense','Ponteareas','A_Coruna_Riazor','Santiago_San_Caetano','A_Coruna_Torre_Hercules','Burela','Rio_Cobo', \
                    'Xove','Pastoriza_A','Paiosaco','Sabon','Centro_Civico','Pastoriza','A_Grela_S','A_Cabana','Fraga_Redonda','Louseiras','Macineira','Magdalena','Marraxon','Mourence', \
                    'Cerceda','Paraxon','San_Vicente_de_Vigo','Villagudin','NNW','Sur','Este_Estacion_1','Oeste_Estacion_2','Areeiro','Campelo','Cee','Dumbria','Campo_de_Futbol', \
                    'Buscas','Rodis','SGL_Carbon','Xubia','Pontevedra_Campolongo','Lope_de_Vega','Xinzo_de_Limia','Lalin','San_Pedro','A_Grela_C']


tarind = [station_names_meta.index(station_names_csv[ii]) for ii in xrange(len(station_names_csv))]
tarlat = [lat[ii] for ii in tarind]
tarlon = [lon[ii] for ii in tarind]
emclass = [exposure[ii] for ii in tarind]

# A_Grela_A = A_Grela_C, Pastoriza = Pastoriza_A
# former A_Grela coordinates 43.345336,-8.434088
