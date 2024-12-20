# -*- coding: utf-8 -*-

def rgbs2string(rgbs_in):
    decs = 2 #decimals used for rounding
    rgbs_out = "(/"
    for ii in range(len(rgbs_in)):
        rgb = rgbs_in[ii]
        #rgbs_out[ii] = '"(/'+str(round(rgb[0],decs))+','+str(round(rgb[1],decs))+','+str(round(rgb[2],decs))+'/)",'
        #rgbs_step = '"(/'+str(round(rgb[0],decs))+','+str(round(rgb[1],decs))+','+str(round(rgb[2],decs))+'/)",' #this works
        #rgbs_step = '(/'+str(round(rgb[0],decs))+','+str(round(rgb[1],decs))+','+str(round(rgb[2],decs))+'/) '
        #rgbs_step = '(/'+str(round(rgb[0],decs))+','+str(round(rgb[1],decs))+','+str(round(rgb[2],decs))+'/), '
        rgbs_step = '(/'+str(round(rgb[0],decs))+','+str(round(rgb[1],decs))+','+str(round(rgb[2],decs))+'/) ' #this works
        rgbs_out = rgbs_out+rgbs_step

    #str_start = rgbs_out[0]
    #str_start = '(/'+str_start
    #rgbs_out[0] = str_start
    #str_end = rgbs_out[-1]
    #str_end = str_end[0:-1]+' )'
    #rgbs_out[-1] = str_end    
    
    rgbs_out = rgbs_out[0:-1]+'/)'
    #rgbs_out = rgbs_out[0:-2]+' )'
    
    return(rgbs_out)
