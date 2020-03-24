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


#path = "/scratch/amartin1/Work/Results/padf/water/onemol/"
path = "/scratch/amartin1/Work/Results/padf/lammps/si_anneal_141215_padf/"
#path = "/scratch/amartin1/Work/Results/padf/lammps/si_anneal_071215_padf/"

tag = "aSi071215_100_noG_"
#tag = "aSi071215_900_"


#cname = path+tag+"_padf2_padf.dbin"
cname = path+tag+"_padfcorr_correlation_sum.dbin"
#cname = path+tag+"_correlation.dbin"
corr = read_correlation( cname , 0)
nr = 64
nth = 402
corr = corr.reshape( nr, nr, nth )
th = np.arange(nth)*360.0/float(nth)


cav = np.average( corr, 2)
cmod = np.copy( corr )
for i in np.arange( corr.shape[2] ):
     cmod[:,:,i] = corr[:,:,i] - cav



plt.imshow( cmod[:,:,0] )
plt.figure()
plt.imshow( cmod[:,42,:] )
plt.draw()
plt.show()
sys.exit()

sp.misc.imsave( cname+'.png', np.abs(corr[:,:,0]) )
plt.figure()
plt.plot( th, corr[8,8,:] )
plt.savefig( cname+'lineplot.png')
sp.misc.imsave( cname+'angle10.png', np.abs(corr[8,:,:])**0.5 )
#plt.draw()
#plt.show()





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
