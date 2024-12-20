# -*- coding: utf-8 -*-

from matplotlib.colors import ListedColormap, LinearSegmentedColormap, BoundaryNorm
from matplotlib.colorbar import ColorbarBase

def hex_to_rgb(value):

    value = value.lstrip('#')
    lv    = len(value)

    return tuple(int(value[i:i+lv/3], 16) for i in range(0, lv, lv/3))

def crea_cmap(clevs, rgbs, under, over):

    rgbs = [under] + rgbs + [over]
    rgbs = [hex_to_rgb(rgb) for rgb in rgbs]

    x0 = float(clevs[0])
    xN = float(clevs[-1])
    y0 = 0.
    yN = 1.

    clevsn = []

    for el in clevs:
        m = ((yN-y0)/(xN-x0))
        b = yN - m*xN
        y = m * float(el) + b
        clevsn.append(y)

    l_r=[]
    l_g=[]
    l_b=[]


    for i in range(len(clevsn)):

        r, g, b = rgbs[i]
        R, G, B = rgbs[i+1]

        l_r.append( (clevsn[i],r/255.,R/255.) )
        l_g.append( (clevsn[i],g/255.,G/255.) )
        l_b.append( (clevsn[i],b/255.,B/255.) )

    cdict={}
    cdict['red'  ] = tuple(l_r)
    cdict['green'] = tuple(l_g)
    cdict['blue' ] = tuple(l_b)

    # Generamos el ListedColormap:
    cmap = LinearSegmentedColormap('mar', cdict)
    cmap.set_over(over)
    cmap.set_under(under)

    return cmap
 
