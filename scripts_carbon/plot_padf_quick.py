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

sample =   "unactivated" # # "C5" # 
dataset =    "12.35.05" #"13.04.29" #  
#dataset = "14.08.14" # "13.52.36" 
#dataset = "11.20.24" # "13.04.29" # 
path = "/Users/e38496/Work/Research/Experiments/Monash/2016/Carbon_dec16/Results/"+sample+"/fem/"+"data_"+dataset+"_recen3_900/"
tag = "carbon_dec16_"+sample+"_"+dataset+"_recen3" # #+"_sym" #

#path = "/Users/e38496/Work/Research/Experiments/Monash/2016/Carbon_dec16/Results/padf_sums/"
#tag = "C9_12.25.40_12.55.43_13.15.36"
#tag = "unactivated_11.20.24_12.35.05_13.04.29"

cname = path+tag+"_padf2_padf.dbin"
#cname = path+tag+"_padf2_padf_averaged.dbin"


usebl0 = False
bname = path+tag+"_padf2_bl0_averaged.dbin"


corr = read_correlation( cname , 0) #* -1
if usebl0 == True:
     bl0 = read_dbin( bname, 0 )

nr = 256
nth = 402
npatterns = 2
print corr.size, nr*nr*nth
corr = corr.reshape( nr, nr, nth )
print "test val :", corr[9,10,10:15], np.max(corr), np.min(corr)
corr *= 1.0/float(npatterns)
corr *= 1.

th = np.arange(nth)*360.0/float(nth)
pix = int(nr * 0.8)
rmax = 29.41 *0.95 # 81.1
lim = ( rmax  ) * (pix/float(nr))

gamma = 1.0
disp = corr_rescale( corr[:pix,:pix,90], gamma )
#plt.imshow( corr[:pix,:pix,67], extent=[0,lim,0,lim], origin='lower' )
plt.imshow( disp, extent=[0,lim,0,lim], origin='lower' )
#plt.imshow( corr[:,:,0], origin='lower' )
plt.xlabel(r'r'+ur'(\u00c5)')
plt.ylabel(r'r$^\prime$'+ur'(\u00c5)')
plt.colorbar()
cmult = 0.1
plt.clim(np.min(disp)*cmult,np.max(disp)*cmult)
#plt.clim(-3e57,3e57) 
#plt.clim(-(8e59**gamma),8e59**gamma) 


print "r samples:", (np.arange(nr) + 0.5 ) * rmax / float(nr)

plt.figure()
#r1, r2 = 50,50
r1, r2 = 1.4, 1.4 #rmax/2, rmax/2 #
ir = int( r1   * (nr/float(rmax)) )
ir2 = int( r2  * (nr/float(rmax)) )
#ir, ir2 = 75,75
print ir, ir2
plt.plot( th, corr[ir,ir2,:] )

plt.figure()
#r1, r2 = 50,50
r1, r2 = 2.9, 2.9 #rmax/2, rmax/2 #
ir = int( r1   * (nr/float(rmax)) )
ir2 = int( r2   * (nr/float(rmax)) )
#ir, ir2 = 75,75
plt.plot( th, corr[ir,ir2,:] )

costh = np.cos(th * 2.0*np.pi / 360.0 )
lp2 = costh**2-0.5
lp2 *= 1.0/np.sqrt(np.sum(lp2**2))
ic = np.where( np.abs(lp2) > 0.05 )
#plt.figure()
#plt.plot( lp2 )


sintharr = np.outer( np.ones(nr), np.abs(np.sin(th * 2.0*np.pi / 360.0 ) ) )

r1 = 1.85
ir = int( r1 * (nr/float(rmax)))
r = np.arange( nr )
disp = np.zeros( (nr, nth) )
dr = np.zeros( nr )
for i in np.arange( nr ):
     if usebl0 == True:
          disp[i,:] = (corr[i,i,:]  +  np.ones(nth)*bl0[i,i]) *r[i]*r[i] #*r[i]*r[i]
     else:
          disp[i,:] = corr[i,i,:] *r[i]*r[i] #*r[i]*r[i] 
#          disp[i,:] = corr[ir,i,:] *r[ir]*r[i] #*r[i]*r[i] 


     #pc = np.sum( lp2*corr[i,i,:] )
     #disp[i,:] = (corr[i,i,:] - pc*lp2) *r[i]*r[i] 

#     disp[i,ic] = corr[i,i,ic] *r[i]*r[i]  / lp2[ic]
     if usebl0 == True:
          dr[i] = (corr[i,i,0] + bl0[i,i] )* r[i] * r[i]
     else:
          dr[i] = corr[i,i,0] #* r[i] * r[i]
plt.figure()
disp = corr_rescale( disp, gamma )
rmin = 0.0
irmin = int( rmin   * (nr/float(rmax)) )
plt.imshow( disp[irmin:pix,:nth/2], extent=[0,180,rmin,lim], origin='lower' , aspect=12 )
plt.ylabel(r'r = r$^\prime$'+ur'(\u00c5)')
plt.xlabel(r'$\theta$ (degrees)')
#plt.clim(-(5e57**gamma),5e57**gamma)
cmult = 0.25
plt.clim(np.min(disp)*cmult,np.max(disp)*cmult)
#plt.clim( -40, -38 )
#plt.clim( 8, 9.6 )
plt.colorbar()

print 1.4, "A is at pixel : ", (1.4/rmax)*nr
print "3.9 A is at pixel : ", (3.9/rmax)*nr

px = np.arange( nr ) 
x = px * rmax / float(nr)

#for i in np.arange( nr ):
#     print px[i], x[i]

x0 = 0
plt.figure()
plt.plot( x[x0:], dr[x0:] )

plt.draw()
plt.show()
