#!/usr/bin/python


import numpy as np
import matplotlib
#matplotlib.use('Agg')
from matplotlib import cm
import matplotlib.pyplot as plt
from matplotlib import rc
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


#rc('font', **{'family':'sans-serif','sans-serif':['Palatino']})
#rc('text', usetex=True)

sample = "C9" #"unactivated" #
#dataset = "13.04.29" # "11.20.24" #"12.35.05" #
dataset = "12.55.43" #  "12.25.40" # 
path = "/Users/e38496/Work/Research/Experiments/Monash/2016/Carbon_dec16/Results/"+sample+"/fem/"+"data_"+dataset+"_recen2/"
tag = "carbon_dec16_"+sample+"_"+dataset+"_recen2"



#cname = path+tag+"_padf2_padf.dbin"
#cname = path+tag+"_padf2_padf_averaged.dbin"
#cname = path+tag+"_padf2_symmetric__padf.dbin"
#cname = path+tag+"_padfXcorr_correlation_bgsum.dbin"
#cname = path+tag+"_padfcorr_correlation_sum_bgsub.dbin"
cname = path+tag+"_padfcorr_correlation_sum_hpfilter.dbin"
#cname = path+tag+"_padfcorr_correlation_sum_symfilter.dbin"
#cname = path+tag+"_padfcorr_correlation_sum_maskcorrected.dbin"
#cname = path+tag+"_mask_correlation.dbin"
#cname = path+tag+"_padfcorr_correlation_sum.dbin"

#cname = path+tag+"_bg_0_correlation.dbin"
#cname = path+tag+"_0_correlation.dbin"


corr = read_correlation( cname , 0) #* -1
nr = 128
nth = 402
npatterns = 2
print corr.size, nr*nr*nth
corr = corr.reshape( nr, nr, nth )
print "test val :", corr[9,10,10:15], np.max(corr), np.min(corr)
corr *= 1.0/float(npatterns)
corr *= 1.


th = np.arange(nth)*360.0/float(nth)
pix = nr #* 1.0
rmax = nr #150.3
lim = ( rmax  ) * (pix/float(nr))

gamma = 0.25


#
# range to set to zero in pixels
#
#r0 = [80, 114]
#r1 = [90, 124]
#th0 = [0, 0]
#th1 = [402, 402]

#
# range for 12.25.40
#
#r0 = [104] #, 74, 87]
#r1 = [128] #, 84, 93]
#th0 = [0] #, 0, 0]
#th1 = [402] #, 402, 402]

#range for 12.55.43
r0 = [74 ]
r1 = [84]
th0 = [0] 
th1 = [402]



#
# filter
#
corr_out = np.copy(corr)
for i in np.arange(len(r0)):
     corr_out[:, r0[i]:r1[i], th0[i]:th1[i]] = 0.0
     corr_out[r0[i]:r1[i], :, th0[i]:th1[i]] = 0.0

#
# write out file
#
outname = cname[:-5]+"_zeromasked_v2.dbin"
write_dbin( outname, corr_out )

#
# plot filtered function
#
r = np.arange( nr )
disp = np.zeros( (nr, nth) )
disp2 = np.zeros( (nr, nth) )
dr = np.zeros( nr )
for i in np.arange( nr ):
     disp[i,:] = corr[i,i,:] *r[i]*r[i] 
     disp2[i,:] = corr_out[i,i,:] *r[i]*r[i] 


disp += - np.outer( np.average( disp, 1), np.ones( nth) )
plt.figure()
disp = corr_rescale( disp, gamma )
plt.imshow( disp[:pix,:nth], extent=[0,360,0,pix], origin='lower' , aspect=1 )
#plt.clim(-(5e57**gamma),5e57**gamma)
plt.clim(np.min(disp)*0.2,np.max(disp)*0.2)

disp2 += - np.outer( np.average( disp2, 1), np.ones( nth) )
plt.figure()
disp2 = corr_rescale( disp2, gamma )
plt.imshow( disp2[:pix,:nth], extent=[0,360,0,pix], origin='lower' , aspect=1 )
#plt.clim(-(5e57**gamma),5e57**gamma)
plt.clim(np.min(disp2)*0.2,np.max(disp2)*0.2)
plt.colorbar()



plt.draw()
plt.show()
