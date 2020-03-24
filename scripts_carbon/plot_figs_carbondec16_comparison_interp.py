#!/usr/bin/python


import numpy as np
import matplotlib
#matplotlib.use('Agg')
#from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
from matplotlib import rc
import sys
import os
import struct
import array
import scipy as sp
import scipy.interpolate as spi

fs = 16
plt.rcParams.update({'font.size': fs})

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

   #  disp[ihigh] = np.log(np.abs( plane[ihigh] ))
   #  disp[ilow] = - np.log(np.abs( plane[ilow] ))

     return disp

# shift - a 2D version of numpy's roll
def array_shift(array,xshift=0,yshift=0):
     array = np.roll(array,xshift,0)
     array = np.roll(array,yshift,1)
     return array


## make a 2D array with a gaussian
def make_gaussian(nx, ny, rad=None, rady=-1., cenx=None, ceny=None, invert=0, norm=False ): 
    #set defaults
    if rad is None: rad = np.min(nx,ny)/2
    if cenx is None: cenx = nx/2
    if ceny is None: ceny = ny/2
    radsq = rad**2
    if rady == -1.:
        radysq = radsq
    else:
        radysq = rady**2

    # define the circle
    x = np.outer(np.arange(0-nx/2,nx-nx/2,1),np.ones(ny))
    #print x.size, x.shape
    y = np.outer(np.ones(nx),np.arange(0-ny/2,ny-ny/2,1))
    #print y.size, y.shape
    a = np.zeros([nx,ny])
    a = np.exp(-(x**2)/radsq  - ( y**2)/radysq)
    a[ nx/2, ny/2 ] = 1.0

    a = array_shift(a,cenx-nx/2,ceny-ny/2)

    # normalise if required
    if norm == True: a *= 1./np.sum(a)
    
    return a
#end make_gaussian

def convolve_gaussian( image, rad=3, rady=1):

    c = make_gaussian( image.shape[0], image.shape[1], rad, rady, cenx=0, ceny=0, norm=True )
    fc = np.fft.fft2( c )
    fimage = np.fft.fft2( image )
    output = np.real(np.fft.ifft2( np.conjugate(fc)*fimage ))
    return output

def interpolate2d( data, lims, samplingfactor=2.0):

     s= data.shape
     l = lims
#     x = np.mgrid[l[0]:l[1]:s[0], l[2]:l[3]:s[1]]
     x = np.mgrid[:s[0], :s[1]]
     print "interpolating..."
     dinterp = spi.interp2d( x[0], x[1], data, kind='linear')                                                                                                                                
#     xnew = np.mgrid[l[0]:l[1]:s[0]*samplingfactor, l[2]:l[3]:s[1]*samplingfactor]
     xnew = np.mgrid[:s[0]:s[0]*int(samplingfactor), :s[1]:s[1]*int(samplingfactor)]
     sx = xnew[0].shape
     ynew = dinterp(xnew[0].flatten(),xnew[1].flatten()).reshape(sx)  #.reshape( (s[0]*int(samplingfactor),s[1]*int(samplingfactor)))
     return ynew


#rc('font', **{'family':'sans-serif','sans-serif':['Palatino']})
#rc('text', usetex=True)


path = "/Users/e38496/Work/Research/Experiments/Monash/2016/Carbon_dec16/Results/padf_sums/"
path2 = path

tag = "unactivated291019__11.20.24_12.35.05_13.04.29"
tag2 = "C9291019__12.25.40_12.55.43_13.15.36"




outpath = "/Users/e38496/Work/Research/Experiments/Monash/2016/Carbon_dec16/Results/figs/" 
outtag = "C9_vs_unactivated_110320"
cname = path+tag+"_padf2_padf_averaged.dbin"
cnamen = path2+tag2+"_padf2_padf_averaged.dbin"
#cname5 = path3+tag3+"_padf2_padf_averaged.dbin"

corr = read_correlation( cname , 0)
corrn = read_correlation( cnamen , 0)
#corr5 = read_correlation( cname5 , 0)
nr = 256 
nth = 402
#natoms = 678
npatterns = 1
corr = corr.reshape( nr, nr, nth )
corrn = corrn.reshape( nr, nr, nth )
#corr5 = corr5.reshape( nr, nr, nth )

#fs = 24 # see top of script

cmax = np.max(corr)
corr *= 1.0/cmax
corrn *= 1.0/cmax  # keep both correlation functions on the same scale


th = np.arange(nth)*360.0/float(nth)

pix = int( nr * 0.17 )
rmax = 29.41 * 0.968 # 0.965 1st version # * 0.987
lim = ( rmax  ) * (pix/float(nr))
x = np.arange( nr ) * rmax / float(nr)
x2 = np.outer( x, x )


gamma = 1.0
#disp = corr_rescale( corr[:pix,:pix,0], gamma )
#disp2 = corr_rescale( corrn[:pix,:pix,0], gamma )
#disp3 = corr_rescale( corr5[:pix,:pix,0], gamma )
disp = corr[:pix,:pix,0] * x2[:pix,:pix] #* x2[:pix,:pix]
disp2 = corrn[:pix,:pix,0] * x2[:pix,:pix] #* x2[:pix,:pix]
#disp3 = corr5[:pix,:pix,0] * x2[:pix,:pix] * x2[:pix,:pix]
disp = corr_rescale( disp, gamma )
disp2 = corr_rescale( disp2, gamma )
#disp3 = corr_rescale( disp3, gamma )

clim = np.max(disp2[:pix,:pix]/2.)
climl = np.min(disp2[:pix,:pix]/2.)

plt.imshow( disp[:pix,:pix], extent=[0,lim,0,lim], origin='lower' )
#plt.imshow( corr[:pix,:pix,0], origin='lower' )
plt.xlabel(r'r'+ur'(\u00c5)', fontsize=fs)
plt.ylabel(r'r$^\prime$'+ur'(\u00c5)', fontsize=fs)
plt.colorbar()
plt.clim(climl,clim)
plt.savefig( outpath+outtag+'_theta0_unactivated.png')

plt.figure()
plt.imshow( disp2[:pix,:pix], extent=[0,lim,0,lim], origin='lower' )
#plt.imshow( corr[:pix,:pix,0], origin='lower' )
plt.xlabel(r'r'+ur'(\u00c5)', fontsize=fs)
plt.ylabel(r'r$^\prime$'+ur'(\u00c5)', fontsize=fs)
plt.colorbar()
plt.clim(climl,clim)
plt.savefig( outpath+outtag+'_theta0_C9.png')

#plt.figure()
#plt.imshow( disp3, extent=[0,lim,0,lim], origin='lower' )
##plt.imshow( corr[:pix,:pix,0], origin='lower' )
#plt.xlabel(r'r'+ur'(\u00c5)')
#plt.ylabel(r'r$^\prime$'+ur'(\u00c5)')
#plt.colorbar()
#plt.clim(climl,clim)
#plt.savefig( outpath+outtag+'_theta0_C5.png')


#clim2 = clim
r1, r2 = 1.45, 1.45 #rmax/2, rmax/2 #
ir = int( r1   * (nr/float(rmax)) )
ir2 = int( r2  * (nr/float(rmax)) )
print "ir, ir2", ir, ir2



thpix = 201
thlim = 360.0 * (thpix/float(nth))
asp = 50 # 32
gamma2 = 0.5
disp = np.zeros( (nr, nth) )
disp2 = np.zeros( (nr, nth) )
disp3 = np.zeros( (nr,nth) ) 
disp4 = np.zeros( (nr, nth) )
disp5 = np.zeros( (nr,nth) )
for i in np.arange( nr ):
     disp[i,:] = corr[i,i,:]
     disp2[i,:] = corrn[i,i,:]
#     disp3[i,:] = corr5[i,i,:]
     disp4[i,:] = corr[ir,i,:]
     disp5[i,:] = corrn[ir,i,:]

x2th = np.outer( x*x, np.ones(nth) )
#disp = disp*x2th
#disp2 = disp2*x2th
#disp4 = disp4*x2th
#disp5 = disp5*x2th


disp = corr_rescale( disp, gamma2 )
disp2 = corr_rescale( disp2, gamma2 )
disp4 = corr_rescale( disp4, gamma2 )
disp5 = corr_rescale( disp5, gamma2 )
#disp3 = corr_rescale( disp3*x2th, gamma2 )

#rwid, thwid = 2, 5
#disp = convolve_gaussian(disp, rwid, thwid )
#disp2 = convolve_gaussian(disp2, rwid, thwid )
#disp5 = convolve_gaussian(disp5, rwid, thwid )
#disp4 = convolve_gaussian(disp4, rwid, thwid )

clim2 = np.max(disp2/8.)     # /5. - 1st version
clim2l = np.min(disp2/8.)    # /5. - 1st version
minres = 1.0
x0 = int((minres/rmax)*nr)


sc = 20
x = np.arange(0, nth, 1)
y = np.arange(0, nr, 1)
f = spi.interp2d(x, y, disp, kind='cubic')
xnew = np.arange(0, nth, 1/float(sc))
ynew = np.arange(0, nr, 1/float(sc))
disp = f(xnew,ynew)
print "length xnew", len(xnew), len(x)
print x0, pix, x0*sc, pix*sc

f = spi.interp2d(x, y, disp2, kind='cubic')
xnew = np.arange(0, nth, 1/float(sc))
ynew = np.arange(0, nr, 1/float(sc))
disp2 = f(xnew,ynew)


plt.figure()
sp = 10.0
#disp_interp = interpolate2d( disp, [0,thlim,minres,lim], sp )
plt.imshow( disp[x0*sc:pix*sc,:thpix*sc], extent=[0,thlim,minres,lim], origin='lower', aspect=asp, interpolation='lanczos' )
#plt.imshow( disp_interp[int(x0*sp):int(pix*sp),:int(thpix*sp)], extent=[0,thlim,minres,lim], origin='lower', aspect=asp, interpolation='lanczos' )
plt.ylabel(r'r'+ur'(\u00c5)'+' = '+r'r$^\prime$'+ur'(\u00c5)', fontsize=fs)
plt.xlabel(r'$\theta$ (degrees)', fontsize=fs)
#plt.colorbar()
plt.clim(clim2l,clim2)

atxt = ['G1', 'G2', 'G3', 'G4', 'G5', 'P1', 'P2', 'A1', 'A2', 'D1','D2','B1','B2']
xpnt = [120, 60, 60, 90, 150, 107, 180-146, 77.8, 108.8, 97, 60, 90, 66]
ypnt = [1.45,2.51, 2.9, 2.7, 2.7, 1.45, 2.33, 2.22, 2.22, 1.56, 1.56, 1.725, 2.25]

# diamond
atxt += ['Dia']*7 # ['Dia1','Dia2','Dia3','Dia4','Dia5','Dia6','Dia7']
xpnt += [70.5, 180-60, 90, 35.5, 50.5, 62.5, 180-85]
ypnt += [1.54, 2.53, 2.53, 2.97, 2.97, 2.97, 2.97]
# more graphene
#atxt += ['G']*5
#xpnt += [21.5, 38.5, 60, 82.5, 60]
#ypnt += [3.75]*4 + [4.25]
# more graphite
atxt += ['Gra']*1
xpnt += [0 ] #, 90, 22.6, 45.4, 71, 35, 72]
ypnt += [3.35 ] #, 4.2, 3.6, 3.6, 3.6, 4.3, 4.3]

#heptagons
atxt += ['H1', 'H2', 'H3']
xpnt += [128.6, 77.2, 26]
ypnt += [1.45, 2.61, 2.61]

dltx, dlty  = [1]*len(atxt), [0.05]*len(atxt)
dltx[0], dlty[0] = 1.0, 0.03
dltx[1], dlty[1] = 2.0, -0.05
dltx[2], dlty[2] = 2.2, -0.05
dltx[5], dlty[5] = 2.5, -0.05
dltx[8], dlty[8] = -2, 0.05
dltx[21], dlty[21] = 1, -0.1
dltx[22], dlty[22] = 1, -0.1
dltx[23], dlty[23] = 1, -0.1

#clr = ['black', 'black', 'black', 'black', 'black', 'red', 'red','green','green','black','black','black','black']
#clr = ['black', 'black', 'black', 'black', 'black', 'black','black','black','black','black','black','black','black']
#clrw = ['white','white','white',    'white','white','white','white','white','white','white','white','white','white']
#mkr = ['o', 'o', 'o', 'o', 'o', 'p','p','^','^','+','+','v', 'v']
clr = ['black']*len(atxt)
clrw = ['white']*len(atxt)
mkr = ['H', 'H', 'H', 'H', 'H', 'p','p','^','^','+','+','v', 'v']
mkr += ['D']*7
#mkr += ['H']*5
mkr += ['s']*1
mkr += ['o']*3


afs = 11
symbol_size = 18

for i in np.arange(len(atxt)):
     plt.annotate(atxt[i], xy=(xpnt[i], ypnt[i]), xytext=(xpnt[i]+dltx[i],ypnt[i]+dlty[i]), xycoords='data',fontsize=afs, color=clr[i])
     plt.scatter(xpnt[i], ypnt[i], s=symbol_size, c=clr[i], marker=mkr[i]) #,edgecolors='b')
     plt.scatter(180.0-xpnt[i], ypnt[i], s=14, c=clrw[i], marker=mkr[i]) #,edgecolors='b')

plt.xlim([0,180])
plt.ylim([1,3.5])
plt.savefig( outpath+outtag+'_r=rprime_theta_unactivated.png')


plt.figure()
plt.imshow( disp2[x0*sc:pix*sc,:thpix*sc], extent=[0,thlim,minres,lim], origin='lower', aspect=asp )
plt.ylabel(r'r'+ur'(\u00c5)'+' = '+r'r$^\prime$'+ur'(\u00c5)', fontsize=fs)
plt.xlabel(r'$\theta$ (degrees)', fontsize=fs)
#plt.colorbar()
plt.clim(clim2l,clim2)
for i in np.arange(len(atxt)):
     plt.annotate(atxt[i], xy=(xpnt[i], ypnt[i]), xytext=(xpnt[i]+dltx[i],ypnt[i]+dlty[i]), xycoords='data',fontsize=afs, color=clr[i])
     plt.scatter(xpnt[i], ypnt[i], s=symbol_size, c=clr[i], marker=mkr[i]) #,edgecolors='b')
     plt.scatter(180.0-xpnt[i], ypnt[i], s=14, c=clrw[i], marker=mkr[i]) #,edgecolors='b')
plt.xlim([0,180])
plt.ylim([1,3.5])
plt.savefig( outpath+outtag+'_r=rprime_theta_C9.png')
#exit()

plt.figure()
plt.imshow( disp4[x0:pix,:thpix], extent=[0,thlim,minres,lim], origin='lower', aspect=asp )
plt.ylabel(r'r'+ur'(\u00c5)', fontsize=fs)
plt.xlabel(r'$\theta$ (degrees)', fontsize=fs)
#plt.colorbar()
plt.clim(clim2l,clim2)
plt.savefig( outpath+outtag+'_r='+'{:.1f}'.format(r1)+'_rprime_theta_unactivated.png')

plt.figure()
plt.imshow( disp5[x0:pix,:thpix], extent=[0,thlim,minres,lim], origin='lower', aspect=asp )
plt.ylabel(r'r'+ur'(\u00c5)', fontsize=fs)
plt.xlabel(r'$\theta$ (degrees)', fontsize=fs)
#plt.colorbar()
plt.clim(clim2l,clim2)
plt.savefig( outpath+outtag+'_r='+'{:.1f}'.format(r1)+'_rprime_theta_C9.png')

#plt.figure()
#plt.imshow( disp3[x0:pix,:thpix], extent=[0,thlim,minres,lim], origin='lower', aspect=asp )
#plt.ylabel(r'r'+ur'(\u00c5)'+' = '+r'r$^\prime$'+ur'(\u00c5)')
#plt.xlabel(r'$\theta$ (degrees)')
##plt.colorbar()
#plt.clim(clim2l,clim2)
#plt.savefig( outpath+outtag+'_r=rprime_theta_C5.png')

#zoom in on a region of the r=r' plot
minres = 1.0
x0 = int((minres/rmax)*nr)
maxres = 3.5
x1 = int((maxres/rmax)*nr)
minth = 0.
th0 = int((minth/360.)*nth)
maxth = 180.
th1 = int((maxth/360.)*nth)
plt.figure()
plt.imshow( disp[x0*sc:x1*sc,th0*sc:th1*sc], extent=[minth,maxth,minres,maxres], origin='lower', aspect=asp )
plt.ylabel(r'r'+ur'(\u00c5)'+' = '+r'r$^\prime$'+ur'(\u00c5)', fontsize=fs)
plt.xlabel(r'$\theta$ (degrees)', fontsize=fs)
#plt.colorbar()
plt.clim(clim2l,clim2)
plt.savefig( outpath+outtag+'_r=rprime_theta_unactivated_zoom.png')


plt.figure()
plt.imshow( disp2[x0*sc:x1*sc,th0*sc:th1*sc], extent=[minth,maxth,minres,maxres], origin='lower', aspect=asp )
plt.ylabel(r'r'+ur'(\u00c5)'+' = '+r'r$^\prime$'+ur'(\u00c5)', fontsize=fs)
plt.xlabel(r'$\theta$ (degrees)', fontsize=fs)
#plt.colorbar()
plt.clim(clim2l,clim2)
plt.savefig( outpath+outtag+'_r=rprime_theta_C9_zoom.png')

plt.figure()
plt.imshow( disp4[x0:x1,th0:th1], extent=[minth,maxth,minres,maxres], origin='lower', aspect=asp )
plt.ylabel(r'r'+ur'(\u00c5)', fontsize=fs)
plt.xlabel(r'$\theta$ (degrees)', fontsize=fs)
#plt.colorbar()
plt.clim(clim2l,clim2)
plt.axvline(x=3.36)
plt.savefig( outpath+outtag+'_r='+'{:.1f}'.format(r1)+'_rprime_theta_unactivated_zoom.png')

plt.figure()
plt.imshow( disp5[x0:x1,th0:th1], extent=[minth,maxth,minres,maxres], origin='lower', aspect=asp )
plt.ylabel(r'r'+ur'(\u00c5)', fontsize=fs)
plt.xlabel(r'$\theta$ (degrees)', fontsize=fs)
#plt.colorbar()
plt.clim(clim2l,clim2)
plt.savefig( outpath+outtag+'_r='+'{:.1f}'.format(r1)+'_rprime_theta_C9_zoom.png')

#plt.figure()
#plt.imshow( disp3[x0:x1,th0:th1], extent=[minth,maxth,minres,maxres], origin='lower', aspect=asp )
#plt.ylabel(r'r'+ur'(\u00c5)'+' = '+r'r$^\prime$'+ur'(\u00c5)')
#plt.xlabel(r'$\theta$ (degrees)')
##plt.colorbar()
#plt.clim(clim2l,clim2)
#plt.savefig( outpath+outtag+'_r=rprime_theta_C5_zoom.png')


exit()


#
#
#     LINE PLOTS FROM HERE
#
#



xlines = [1.5, 2.3, 3.3] 
xlabels = [str(x)+ur'\u00c5' for x in xlines]
xlpos = 0.05

dr = np.zeros( nr )
drn = np.zeros( nr )
dr5 = np.zeros( nr )
x = np.arange( nr ) * rmax / float(nr)
for i in np.arange(nr):
     dr[i] = corr[i,i,0]
     drn[i] = corrn[i,i,0]
#     dr5[i] = corr5[i,i,0]
fig, ax = plt.subplots()
#plt.figure()
xlim = 128
x0=0
p,  = plt.plot( x[:xlim], dr[:xlim] ) #*x[x0:xlim]*x[x0:xlim] )
pa, = plt.plot( x[:xlim], drn[:xlim] ) # *x[x0:xlim]*x[x0:xlim] )
#pb, = plt.plot( x[:xlim], dr5[:xlim] ) #*x[x0:xlim]*x[x0:xlim] )
plt.ylabel(r'$\Theta(r,r^\prime,\theta) - G(r)G(r^\prime)$', fontsize=fs)
plt.xlabel(r'$r$ ('+ur'\u00c5'+')', fontsize=fs)
plt.title( 'r = r$^\prime$;'+r'$\theta$'+'= 0', fontsize=fs)
#plt.legend( [p,pa, pb], ['Unactivated FEM 900 frames', 'C9 FEM 900 frames', 'C5 FEM 225 frames' ])
plt.legend( [p,pa], ['Unactivated', 'Activated' ], fontsize=fs)
x_bounds = ax.get_xlim()
for i, xv in enumerate(xlines):
     plt.axvline( xv, label=xlabels[i], color='k')
     ax.annotate(s=xlabels[i], xy =(((xv+xlpos-x_bounds[0])/(x_bounds[1]-x_bounds[0])),0.9), xycoords='axes fraction', verticalalignment='right', horizontalalignment='right bottom' , rotation = 270)

plt.savefig( outpath+outtag+'diagonal_line_plot.png')

#plt.figure()
fig, ax = plt.subplots()
xlim = 128
x0=0
p,  = plt.plot( x[:xlim], dr[:xlim] *x[x0:xlim]*x[x0:xlim]) # *x[x0:xlim]*x[x0:xlim] )
pa, = plt.plot( x[:xlim], drn[:xlim] *x[x0:xlim]*x[x0:xlim] ) #] *x[x0:xlim]*x[x0:xlim] )
#pb, = plt.plot( x[:xlim], dr5[:xlim]*x[x0:xlim]*x[x0:xlim] ) #*x[x0:xlim]*x[x0:xlim] )
plt.ylabel(r'$\Theta(r,r^\prime,\theta) - G(r)G(r^\prime)$', fontsize=fs)
plt.xlabel(r'$r$ ('+ur'\u00c5'+')', fontsize=fs)
plt.title( 'r = r$^\prime$;'+r'$\theta$'+'= 0', fontsize=fs)
#plt.legend( [p,pa, pb], ['Unactivated FEM 900 frames', 'C9 FEM 900 frames', 'C5 FEM 225 frames' ])
for i, xv in enumerate(xlines):
     plt.axvline( xv, label=xlabels[i], color='k')
     ax.annotate(s=xlabels[i], xy =(((xv+xlpos-x_bounds[0])/(x_bounds[1]-x_bounds[0])),0.9), xycoords='axes fraction', verticalalignment='right', horizontalalignment='right bottom' , rotation = 270)
plt.legend( [p,pa], ['Unactivated', 'Activated' ], fontsize=fs)
plt.savefig( outpath+outtag+'diagonal_line_plot2.png')

#plt.figure()
fig, ax = plt.subplots()
x0 = 12
p,  = plt.plot( x[x0:xlim], dr[x0:xlim] )
pa, = plt.plot( x[x0:xlim], drn[x0:xlim] )
#pb, = plt.plot( x[x0:xlim], dr5[x0:xlim] )
plt.ylabel(r'$\Theta(r,r^\prime,\theta) - G(r)G(r^\prime)$', fontsize=fs)
plt.xlabel(r'$r$ ('+ur'\u00c5'+')', fontsize=fs)
plt.title( 'r = r$^\prime$;'+r'$\theta$'+'= 0', fontsize=fs)
#plt.legend( [p,pa, pb], ['Unactivated FEM', 'C9 FEM', 'C5 FEM' ])
plt.xlim([1.4,6])
plt.ylim([0, 0.2])
x_bounds = ax.get_xlim()
for i, xv in enumerate(xlines):
     plt.axvline( xv, label=xlabels[i], color='k')
     ax.annotate(s=xlabels[i], xy =(((xv+xlpos-x_bounds[0])/(x_bounds[1]-x_bounds[0])),0.8), xycoords='axes fraction', verticalalignment='right', horizontalalignment='right bottom' , rotation = 270)
plt.legend( [p,pa], ['Unactivated', 'Activated' ], fontsize=fs)
plt.yticks(np.arange(4)*0.1)
plt.savefig( outpath+outtag+'diagonal_line_plot_zoom.png')

#plt.figure()
fig, ax = plt.subplots()
x0 = 12
p,  = plt.plot( x[x0:xlim], dr[x0:xlim]*x[x0:xlim]*x[x0:xlim] )
pa, = plt.plot( x[x0:xlim], drn[x0:xlim]*x[x0:xlim]*x[x0:xlim] )
#pb, = plt.plot( x[x0:xlim], dr5[x0:xlim]*x[x0:xlim]*x[x0:xlim] )
plt.ylabel(r'$\Theta(r,r^\prime,\theta) - G(r)G(r^\prime)$', fontsize=fs)
plt.xlabel(r'$r$ ('+ur'\u00c5'+')', fontsize=fs)
plt.title( 'r = r$^\prime$;'+r'$\theta$'+'= 0', fontsize=fs)
#plt.legend( [p,pa, pb], ['Unactivated FEM', 'C9 FEM', 'C5 FEM' ])
plt.xlim([1.4,6])
plt.ylim([0, 0.6])
for i, xv in enumerate(xlines):
     plt.axvline( xv, label=xlabels[i], color='k')
     ax.annotate(s=xlabels[i], xy =(((xv+xlpos-x_bounds[0])/(x_bounds[1]-x_bounds[0])),0.8), xycoords='axes fraction', verticalalignment='right', horizontalalignment='right bottom' , rotation = 270)
plt.legend( [p,pa], ['Unactivated', 'Activated' ], fontsize=fs)
plt.savefig( outpath+outtag+'diagonal_line_plot_zoom2.png')



pk = [ 1.4, 2.3, 2.4, 2.9, 1.8, 1.6, 3.35] 
for p in pk:
     print p, "peak at :", p * nr / rmax
#print "peak at :", 2.5 * nr / rmax
#print "peak at :", 2.3 * nr / rmax
#print "peak at :", 2.9 * nr / rmax
#print "peak at :", 1.8 * nr / rmax

#plt.draw()
#plt.show()
#sys.exit()

pr  = np.array([ 13, 15, 17, 13, 13, 13, 13, 13, 13, 21, 22, 27, 13, 13] )
pr2  = np.array([ 13, 15, 17, 17, 20, 21, 22, 23, 15, 21, 22, 27, 30, 33] )
prx  = pr  * rmax / float(nr)
pr2x = pr2 * rmax / float(nr)


u = np.arange( nth/2)*(4/float(nth)) - 1.0 +(1.0/float(nth))
jac = 1.0 / np.sqrt( 1 - u**2 )
lines  = []
lines2 = []
lines3 = []
for i in np.arange( len(pr) ):
     lines.append( corr[ pr[i], pr2[i], :] )
     lines2.append( corrn[ pr[i], pr2[i], :] )
  #   lines3.append( corr5[ pr[i], pr2[i], :] )

for i in np.arange( len(pr) ):
     plt.figure()
     p, = plt.plot( th[:nth/2], lines[i][:nth/2]  )
     pa, = plt.plot( th[:nth/2], lines2[i][:nth/2] )
 #    pb, = plt.plot( th, lines3[i] )
     # plt.axis([0,360,y0a,y1a])
     plt.ylabel(r'$\Theta(r,r^\prime,\theta) - G(r)G(r^\prime)$', fontsize=fs)
     plt.xlabel(r'$\theta$ (degrees)', fontsize=fs)
     leg1 = plt.Rectangle((0, 0), 0, 0, alpha=0.0)
     # plt.legend( [leg1], ['r = 2.5 '+ur'\u00c5'+'\n'  + r'r$^\prime$ = 2.5 '+ur'\u00c5'],
     #            handlelength = 0)
     plt.title( 'r = '+'{:.1f}'.format(prx[i])+' '+ur'\u00c5'+'  '+r'r$^\prime$'+'{:.1f}'.format(pr2x[i])+' '+ur'\u00c5', fontsize=fs)
   #  plt.legend( [p,pa, pb], ['Unactivated FEM 900 frames', 'C9 FEM 900 frames', 'C5 FEM 225 frames' ])
     plt.legend( [p,pa], ['Unactivated', 'Activated' ], fontsize=fs)
     plt.tight_layout()
     plt.savefig( outpath+outtag+'_r'+str(prx[i])+'_'+'r'+str(pr2x[i])+'.png')

plt.draw()
plt.show()
sys.exit()



#
# EXTRA ###################################################################
#
y0, y1 = -2, 6
y0a, y1a = -4, 12
plt.figure()
p, = plt.plot( th, line )
pa, = plt.plot( th, nline )
pb, = plt.plot( th, cline )
#plt.axis([0,360,y0a,y1a])
plt.ylabel(r'$\Theta(r,r^\prime,\theta) - G(r)G(r^\prime)$')
plt.xlabel(r'$\theta$ (degrees)')
leg1 = plt.Rectangle((0, 0), 0, 0, alpha=0.0)
#plt.legend( [leg1], ['r = 2.5 '+ur'\u00c5'+'\n'  + r'r$^\prime$ = 2.5 '+ur'\u00c5'],
#            handlelength = 0)
plt.title( 'r = 2.3 '+ur'\u00c5'+'  '  + r'r$^\prime$ = 2.3 '+ur'\u00c5')
plt.legend( [p,pa,pb], ['paracrystal  Si', 'aSi', 'crystal Si' ])
plt.savefig( outpath+outtag+'first_shell_comp.png')

plt.figure()
p2, = plt.plot( th, line2 )
p2a, = plt.plot( th, nline2 )
p2b, = plt.plot( th, cline2 )
#plt.axis([0,360,y0,y1])
plt.ylabel(r'$\Theta(r,r^\prime,\theta) - G(r)G(r^\prime)$')
plt.xlabel(r'$\theta$ (degrees)')
leg1 = plt.Rectangle((0, 0), 0, 0, alpha=0.0)
#plt.legend( [leg1], ['r = 4.3 '+ur'\u00c5'+'\n'  + r'r$^\prime$ = 4.3 '+ur'\u00c5'],
#            handlelength = 0)
plt.title( 'r = 4.0 '+ur'\u00c5'+'  '  + r'r$^\prime$ = 4.0 '+ur'\u00c5')
#plt.legend( [leg1], ['second neighbour shell'] )
plt.legend( [p2,p2a,p2b], ['paracrystal  Si', 'aSi', 'crystal Si' ])
plt.savefig( outpath+outtag+'second_shell_comp.png')


plt.figure()
p3, = plt.plot( th, line3 )
p3a, = plt.plot( th, nline3 )
p3b, = plt.plot( th, cline3 )
#plt.axis([0,360,y0,y1])
plt.ylabel(r'$\Theta(r,r^\prime,\theta) - G(r)G(r^\prime)$')
plt.xlabel(r'$\theta$ (degrees)')
leg1 = plt.Rectangle((0, 0), 0, 0, alpha=0.0)
#plt.legend( [leg1], ['r = 4.3 '+ur'\u00c5'+'\n'  + r'r$^\prime$ = 2.5 '+ur'\u00c5'],
#            handlelength = 0)
#plt.legend( [leg1], ['third neighbour shell'] )
plt.title( 'r = 6.0 '+ur'\u00c5'+'  '  + r'r$^\prime$ = 6.0 '+ur'\u00c5')
plt.legend( [p3,p3a,p3b], ['paracrystal  Si', 'aSi', 'crystal Si' ])
plt.savefig( outpath+outtag+'third_shell_comp.png')

plt.figure()
p4, = plt.plot( th, line4 )
p4a, = plt.plot( th, nline4 )
p4b, = plt.plot( th, cline4 )
#plt.axis([0,360,y0,y1])
plt.ylabel(r'$\Theta(r,r^\prime,\theta) - G(r)G(r^\prime)$')
plt.xlabel(r'$\theta$ (degrees)')
leg1 = plt.Rectangle((0, 0), 0, 0, alpha=0.0)
#plt.legend( [leg1], ['r = 4.3 '+ur'\u00c5'+'\n'  + r'r$^\prime$ = 2.5 '+ur'\u00c5'],
#            handlelength = 0)
#plt.legend( [leg1], [' 4A / 8A'] )
plt.title( 'r = 4.0 '+ur'\u00c5'+'  '  + r'r$^\prime$ = 8.0 '+ur'\u00c5')
plt.legend( [p4,p4a,p4b], ['paracrystal  Si', 'aSi', 'crystal Si' ])
plt.savefig( outpath+outtag+'4A8A_shell_comp.png')

plt.figure()
p5, = plt.plot( th, line5 )
p5a, = plt.plot( th, nline5 )
p5b, = plt.plot( th, cline5 )
#plt.axis([0,360,y0,y1])
plt.ylabel(r'$\Theta(r,r^\prime,\theta) - G(r)G(r^\prime)$')
plt.xlabel(r'$\theta$ (degrees)')
leg1 = plt.Rectangle((0, 0), 0, 0, alpha=0.0)
#plt.legend( [leg1], ['r = 4.3 '+ur'\u00c5'+'\n'  + r'r$^\prime$ = 2.5 '+ur'\u00c5'],
#            handlelength = 0)
#plt.legend( [leg1], ['fourth neighbour shell'] )
plt.title( 'r = 8.0 '+ur'\u00c5'+'  '  + r'r$^\prime$ = 8.0 '+ur'\u00c5')
plt.legend( [p5,p5a,p5b], ['paracrystal  Si', 'aSi', 'crystal Si' ])
plt.savefig( outpath+outtag+'fourth_shell_comp.png')


plt.figure()
p6, = plt.plot( th, line6 )
p6a, = plt.plot( th, nline6 )
p6b, = plt.plot( th, cline6 )
#plt.axis([0,360,y0,y1])
plt.ylabel(r'$\Theta(r,r^\prime,\theta) - G(r)G(r^\prime)$')
plt.xlabel(r'$\theta$ (degrees)')
leg1 = plt.Rectangle((0, 0), 0, 0, alpha=0.0)
#plt.legend( [leg1], ['r = 4.3 '+ur'\u00c5'+'\n'  + r'r$^\prime$ = 2.5 '+ur'\u00c5'],
#            handlelength = 0)
#plt.legend( [leg1], ['fourth neighbour shell'] )
plt.title( 'r = 8.0 '+ur'\u00c5'+'  '  + r'r$^\prime$ = 8.0 '+ur'\u00c5')
plt.legend( [p6,p6a,p6b], ['paracrystal  Si', 'aSi', 'crystal Si' ])
plt.savefig( outpath+outtag+'2.3A_4A_comp.png')


plt.draw()
plt.show()
sys.exit()



sp.misc.imsave( cname+'.png', np.abs(corr[:,:,0]) )
plt.figure()
plt.plot( th, corr[8,8,:] )
plt.savefig( cname+'lineplot_comp.png')
sp.misc.imsave( cname+'angle10_comp.png', np.abs(corr[8,:,:])**0.5 )
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
