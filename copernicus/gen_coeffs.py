#!/usr/bin/env python

from netCDF4 import Dataset
import numpy as np

LUSTRE='/mnt/lustre/scratch//home/cesga/orballo'
nc = Dataset(LUSTRE+'/OP/DATOS/CICC/chimere2017r4/MACC/coeffs.nc', 'w',format='NETCDF3_64BIT') 
lev = nc.createDimension('lev', 60)
alpha = nc.createVariable('HYAM',np.float32, ('lev',))
beta = nc.createVariable('HYBM',np.float32, ('lev',))
#alpha[:] = np.array((20.000000,38.425343,63.647804,95.636963,134.483307,180.584351,234.779053,298.495789,373.971924,464.618134,575.651001,713.218079,883.660522,1094.834717,1356.474609,1680.640259,2082.273926,2579.888672,3196.421631,3960.291504,4906.708496,6018.019531,7306.631348,8765.053711,10376.126953,12077.446289,13775.325195,15379.805664,16819.474609,18045.183594,19027.695313,19755.109375,20222.205078,20429.863281,20384.480469,20097.402344,19584.330078,18864.750000,17961.357422,16899.468750,15706.447266,14411.124023,13043.218750,11632.758789,10209.500977,8802.356445,7438.803223,6144.314941,4941.778320,3850.913330,2887.696533,2063.779785,1385.912598,855.361755,467.333588,210.393890,65.889244,7.367743,0.000000,0.000000))
#beta[:] = np.array((0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000076,0.000461,0.001815,0.005081,0.011143,0.020678,0.034121,0.051690,0.073534,0.099675,0.130023,0.164384,0.202476,0.243933,0.288323,0.335155,0.383892,0.433963,0.484772,0.535710,0.586168,0.635547,0.683269,0.728786,0.771597,0.811253,0.847375,0.879657,0.907884,0.931940,0.951822,0.967645,0.979663,0.988270,0.994019,0.997630,1.000000))
alpha[:] = np.array((10, 29.2126693725586, 51.0365676879883, 79.6423797607422, 
    115.060134887695, 157.533828735352, 207.681701660156, 266.637451171875, 
    336.23388671875, 419.295043945312, 520.134643554688, 644.4345703125, 
    798.439208984375, 989.24755859375, 1225.65466308594, 1518.55749511719, 
    1881.45715332031, 2331.08129882812, 2888.15515136719, 3578.35656738281, 
    4433.49926757812, 5462.36328125, 6662.326171875, 8035.84375, 
    9570.58984375, 11226.78515625, 12926.384765625, 14577.564453125, 
    16099.638671875, 17432.328125, 18536.439453125, 19391.40234375, 
    19988.65625, 20326.033203125, 20407.171875, 20240.94140625, 
    19840.865234375, 19224.5390625, 18413.0546875, 17430.4140625, 
    16302.958984375, 15058.787109375, 13727.171875, 12337.98828125, 
    10921.12890625, 9505.927734375, 8120.580078125, 6791.560546875, 
    5543.046875, 4396.34533691406, 3369.30493164062, 2475.73815917969, 
    1724.84619140625, 1120.63720703125, 661.34765625, 338.863693237305, 
    138.141563415527, 36.6284894943237, 3.68387126922607, 0))
beta[:] = np.array((0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    3.79117482225411e-05, 0.000268609197519254, 0.00113827548921108, 
    0.00344813661649823, 0.00811201333999634, 0.01591039262712, 
    0.0273995194584131, 0.0429057851433754, 0.0626121200621128, 
    0.0866042673587799, 0.114848613739014, 0.147203415632248, 
    0.183430105447769, 0.223204523324966, 0.266128063201904, 
    0.311738938093185, 0.359523504972458, 0.408927530050278, 
    0.459367245435715, 0.510240733623505, 0.560939162969589, 
    0.610857933759689, 0.659408032894135, 0.706027209758759, 
    0.750191211700439, 0.79142501950264, 0.829314172267914, 
    0.863515913486481, 0.893770396709442, 0.919912099838257, 
    0.941880911588669, 0.959733366966248, 0.973653972148895, 
    0.983966410160065, 0.991144776344299, 0.995824784040451, 0.998815059661865))
nc.close()
