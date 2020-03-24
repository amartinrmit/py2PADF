#!/usr/bin/python


import numpy as np
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


path = "/scratch/amartin1/Work/Results/padf/incoh-150615/"
path2 = "/scratch/amartin1/Work/Results/padf/debug-sampling-090615-2/"
tag = "polyhedra"
tag2 = "polyhedra"


#cname = path+tag+"_padf2_padf.dbin"
cname = path+tag+"_padfcorr_correlation_sum.dbin"
cname2 = path2+tag2+"_padfcorr_correlation_sum.dbin"

nr = 64
nth = 200

corr = read_correlation( cname , 0)
corr = corr.reshape( nr, nr, nth ) / 5.0

corr2 = read_correlation( cname2 , 0)
corr2 = corr2.reshape( nr, nr, nth ) 

cav = np.average( corr, 2 )
for i in np.arange(nth):
     corr[:,:,i] += - cav 


cav2 = np.average( corr2, 2 )
for i in np.arange(nth):
     corr2[:,:,i] += - cav2 


diff  = corr - corr2

print "np.max(diff):", np.max(diff), np.max(corr)

plt.imshow( diff[20,:,:]**1.0 )

plt.figure()
r, r2 = 20, 20
plt.plot( corr[r,r2,:] - corr2[r,r2,:] )
plt.plot( corr[r,r2,:] )
plt.plot( corr2[r,r2,:] )

plt.draw()
plt.show()





#bname = "/home/amartin1/Work/codes/C-code/padfcorr/blar0.dbin"
#br0 = read_dbin( bname, 0 )
#print "test br0 max min :", np.max(br0), np.min(br0)

#bname = "/home/amartin1/Work/codes/C-code/padfcorr/blaq0.dbin"
#bq0 = read_dbin( bname, 0 )
#print "test bq0 max min :", np.max(bq0), np.min(bq0)

#plt.figure()
#plt.imshow( np.abs(br0)**0.3 )


#plt.draw()
#plt.show()
