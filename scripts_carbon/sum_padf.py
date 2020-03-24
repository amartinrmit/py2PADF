#!/usr/bin/python


import numpy as np
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import os
import struct
import array
import scipy as sp

def read_dbin( fname, swapbyteorder=0 ):

     size = os.path.getsize(fname)
     #print "sizes:", size, size/8
     b = array.array('d')
     f = open(fname, "r")
     b.fromfile( f, size/8 )
     f.close();
     l = b.tolist()
     
     n = len(l)
     nx = int(round(n**(1.0/2.0)))
     print "nx", nx
     output = np.array(l).reshape( nx, nx)
     if swapbyteorder == 1: output = output.newbyteorder()
     return output

def write_dbin( fname, data ):

	f = open( fname, "wb")
	fmt='<'+'d'*data.size
	bin = struct.pack(fmt, *data.flatten()[:] )
	f.write( bin )
	f.close()



def read_correlation( fname, swapbyteorder=0 ):

     size = os.path.getsize(fname)
     #print "sizes:", size, size/8
     b = array.array('d')
     f = open(fname, "r")
     b.fromfile( f, size/8 )
     f.close();
     l = b.tolist()
     output = np.array(l)
     
     if swapbyteorder == 1: output = output.newbyteorder()
     return output

def corr_rescale( plane, gamma ):

     disp = plane*0.0
     ihigh = np.where( plane > 0 )
     ilow = np.where( plane < 0 )
     disp[ihigh] = np.abs( plane[ihigh] )**gamma
     disp[ilow] = - np.abs( plane[ilow] )**gamma

     return disp

#runtag = "_20_20_20nm_20nm_CL_230mm_CA_5um_0_5s_1_npat400_hp55_centred_meancorrected_new_sym" #_hp20_centred"
sample =  "C9"  #"unactivated" 
ldcodes =  ["12.25.40","12.55.43", "13.15.36" ]     # C9 
#ldcodes =  ["11.20.24","12.35.05", "13.04.29"]       # unactivated
#rtag =  ["","_sym", ""]  #C9
#rtag =  ["","", ""]

dtag =  ["_recen2_fz2","_recen2_fz2", "_recen2"]
rtag =  ["_recen2_fz2","_recen2_fz2", "_recen2"]

#dtag =  ["_recen3","_recen3_900", "_recen2"]
#rtag =  ["_recen3","_recen3", "_recen2"]

nr = 256
nth = 402

i=0
dstr = '_x291019_'
for dcode in ldcodes:
     path = "/Users/e38496/Work/Research/Experiments/Monash/2016/Carbon_dec16/Results/"+sample+"/fem/data_"+dcode+dtag[i]+"/"
     tag = "carbon_dec16_"+sample+"_"+dcode+rtag[i]
     cname = path+tag+"_padf2_padf.dbin"
     bname = path+tag+"_padf2blrr_0.dbin"

     
     corr = read_correlation( cname , 0)
     bl0 = read_dbin( bname, 0 )
  #  corr = corr.reshape( nr, nr, nth )

     if i==0:
          corrsum = np.copy(corr)
          blsum = np.copy(bl0)
          corrsumsq = corr*corr
     else:
          corrsum += np.copy(corr)
          blsum += np.copy(bl0)
          corrsumsq += corr*corr
    
     dstr = dstr+"_"+dcode
     i += 1

corrsum *= 1.0/float(i)
blsum *= 1.0/float(i)
corrsumsq *= 1/float(i)
corrsumsq = np.sqrt( corrsumsq - corrsum*corrsum ) 

path = "/Users/e38496/Work/Research/Experiments/Monash/2016/Carbon_dec16/Results/padf_sums/"
tag = sample+dstr
cname = path+tag+"_padf2_padf_averaged.dbin"
sname = path+tag+"_padf2_padf_sigma.dbin"
bname = path+tag+"_padf2_bl0_averaged.dbin"
write_dbin( cname, corrsum )
write_dbin( bname, blsum )
write_dbin( sname, corrsumsq )


