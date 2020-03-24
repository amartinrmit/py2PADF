#!/usr/bin/python


import numpy as np
import matplotlib
#matplotlib.use('Agg')
from matplotlib import cm
import matplotlib.pyplot as plt
from matplotlib import rc
#from mpl_toolkits.axes_grid1 import make_axes_locatable
import sys
import os
import struct
import array
import scipy as sp

plt.rcParams.update({'font.size': 24})

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

sample = "unactivated" # "C9" #"C5"  #
dataset = "12.35.05" #"13.15.36" #"13.04.29" # "11.20.24" # #
#dataset = "13.04.29"  ## "14.08.14" #    "13.15.36" # 
path = "/Users/e38496/Work/Research/Experiments/Monash/2016/Carbon_dec16/Results/"+sample+"/fem/"+"data_"+dataset+"_recen2/"
tag = "carbon_dec16_"+sample+"_"+dataset+"_recen2"



#cname = path+tag+"_padf2_padf.dbin"
#cname = path+tag+"_padf2_padf_averaged.dbin"
#cname = path+tag+"_padf2_symmetric__padf.dbin"
#cname = path+tag+"_padfXcorr_correlation_bgsum.dbin"
#cname = path+tag+"_padfcorr_correlation_sum_bgsub.dbin"
#cname = path+tag+"_padfcorr_correlation_sum_hpfilter_zeromasked_v2.dbin"
#cname = path+tag+"_padfcorr_correlation_sum_hpfilter.dbin"
#cname = path+tag+"_padfcorr_correlation_sum_maskcorrected.dbin"
#cname = path+tag+"_padfcorr_correlation_sum_symfilter.dbin"
#cname = path+tag+"_padfcorr_correlation_sum_maskcorrected.dbin"
#cname = path+tag+"_padfcorr_correlation_sum_zeromasked.dbin"
#cname = path+tag+"_mask_correlation.dbin"
cname = path+tag+"_padfcorr_correlation_sum.dbin"

#cname = path+tag+"_bg_0_correlation.dbin"
#cname = path+tag+"_0_correlation.dbin"


corr = read_correlation( cname , 0) #* -1
nr = 128
nth = 402
npatterns = 1
print corr.size, nr*nr*nth
corr = corr.reshape( nr, nr, nth )
print "test val :", corr[9,10,10:15], np.max(corr), np.min(corr)
#corr *= 1.0/float(npatterns)
corr *= 1. / np.max(corr)

th = np.arange(nth)*360.0/float(nth)
pix = int(nr * 1.0)
rmax = nr #45.65 * 1.21   #2.36/( nbin)
lim = ( rmax  ) * (pix/float(nr))

gamma = 0.25
disp = corr_rescale( corr[:pix,:pix,0], gamma )
#plt.imshow( corr[:pix,:pix,67], extent=[0,lim,0,lim], origin='lower' )
plt.imshow( disp, extent=[0,lim,0,lim], origin='lower' )
#plt.imshow( corr[:,:,0], origin='lower' )
plt.xlabel(r'r'+ur'(\u00c5)')
plt.ylabel(r'r$^\prime$'+ur'(\u00c5)')
plt.colorbar()
#plt.clim(-3e57,3e57) 
#plt.clim(-(8e59**gamma),8e59**gamma) 


plt.figure()
#r1, r2 = 50,50
r1, r2 = 1.4, 1.4 #rmax/2, rmax/2 #
ir = ( r1  ) * (nr/float(rmax))
ir2 = ( r2  ) * (nr/float(rmax))
ir, ir2 = 75,75
plt.plot( th, corr[ir,ir2,:] )

costh = np.cos(th * 2.0*np.pi / 360.0 )
lp2 = costh**2-0.5
lp2 *= 1.0/np.sqrt(np.sum(lp2**2))
ic = np.where( np.abs(lp2) > 0.05 )
#plt.figure()
#plt.plot( lp2 )


sintharr = np.outer( np.ones(nr), np.abs(np.sin(th * 2.0*np.pi / 360.0 ) ) )

r = np.arange( nr )
disp = np.zeros( (nr, nth) )
dr = np.zeros( nr )
for i in np.arange( nr ):
     disp[i,:] = ( corr[i,i,:] - np.average(corr[i,i,:]) ) *r[i]*r[i] # # #*r[i]*r[i] 

     #pc = np.sum( lp2*corr[i,i,:] )
     #disp[i,:] = (corr[i,i,:] - pc*lp2) *r[i]*r[i] 

#     disp[i,ic] = corr[i,i,ic] *r[i]*r[i]  / lp2[ic]
     dr[i] = corr[i,i,0] #* r[i] * r[i]
plt.figure()
disp = corr_rescale( disp, gamma )
qscale = 2.1758000000
x0 = 0
qx0 = qscale * (x0 / float(nr) )
disp *= 1.0/np.max(disp[x0:pix,:nth])
im = plt.imshow( disp[x0:pix,:nth], extent=[0,360,qx0,qscale], origin='lower' , aspect=150 )
#plt.clim(-(5e57**gamma),5e57**gamma)
cscale = 1.0
plt.clim(np.min(disp[x0:pix,:nth])*cscale,np.max(disp[x0:pix,:nth])*cscale)
#plt.clim( -40, -38 )
#plt.clim( 8, 9.6 )
plt.xlabel(r'$\theta$ (degrees)')
#plt.ylabel(r'q'+ur' (x 10$^{10}}$ \u00c5 $^{-1}$)')
plt.ylabel(r'q'+ur' (\u00c5 $^{-1}$)')


plt.colorbar(im,fraction=0.03, pad=0.04)
plt.tight_layout()

#ax = plt.gca()
# create an axes on the right side of ax. The width of cax will be 5%
# of ax and the padding between cax and ax will be fixed at 0.05 inch.
#divider = make_axes_locatable(ax)
#cax = divider.append_axes("right", size="5%", pad=0.05)
#plt.colorbar(im, cax=cax)


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
