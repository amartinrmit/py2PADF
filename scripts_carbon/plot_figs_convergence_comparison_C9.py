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

plt.rcParams.update({'font.size': 16})

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


#rc('font', **{'family':'sans-serif','sans-serif':['Palatino']})
#rc('text', usetex=True)

path = "/Users/e38496/Work/Research/Experiments/Monash/2016/Carbon_dec16/Results/C9/fem/"
path2 = path+"data_12.55.43_recen2_fz2/"
path3 = path+"data_12.25.40_recen2_fz2/"
path  = path+"data_13.15.36_recen2/"

tag = "carbon_dec16_C9_13.15.36_recen2"
tag2 = "carbon_dec16_C9_12.55.43_recen2_fz2"
tag3 = "carbon_dec16_C9_12.25.40_recen2_fz2"



label1 = "activated region1"
label2 = "activated region2"
label3 = "activated region3"


outpath = "/Users/e38496/Work/Research/Experiments/Monash/2016/Carbon_dec16/Results/figs/" 
outtag = "activated_convergence_281019"
cname = path+tag+"_padf2_padf.dbin"
cnamen = path2+tag2+"_padf2_padf.dbin"
cname5 = path3+tag3+"_padf2_padf.dbin"

corr = read_correlation( cname , 0)
corrn = read_correlation( cnamen , 0)
corr5 = read_correlation( cname5 , 0)
nr = 256 
nth = 402
#natoms = 678
npatterns = 1
corr = corr.reshape( nr, nr, nth )
corrn = corrn.reshape( nr, nr, nth )
corr5 = corr5.reshape( nr, nr, nth )

fs = 16

#cmax = np.max(corr)
cmax = 1.1550020221035933e+63 #2.9426273582118624e+66   #from unactivated case
corr *= 1.0/cmax
corrn *= 1.0/cmax  # keep both correlation functions on the same scale
corr5 *= 1.0/cmax  # keep both correlation functions on the same scale


th = np.arange(nth)*360.0/float(nth)

pix = int( nr * 0.10 )
rmax = 29.41 * 0.95
lim = ( rmax  ) * (pix/float(nr))
x = np.arange( nr ) * rmax / float(nr)
x2 = np.outer( x, x )
minres = 1.0
x0 = int((minres/rmax)*nr)


gamma = 0.5
#disp = corr_rescale( corr[:pix,:pix,0], gamma )
#disp2 = corr_rescale( corrn[:pix,:pix,0], gamma )
#disp3 = corr_rescale( corr5[:pix,:pix,0], gamma )
disp = corr[:pix,:pix,0] * x2[:pix,:pix] #* x2[:pix,:pix]
disp2 = corrn[:pix,:pix,0] * x2[:pix,:pix] #* x2[:pix,:pix]
disp3 = corr5[:pix,:pix,0] * x2[:pix,:pix] * x2[:pix,:pix]
disp = corr_rescale( disp, gamma )
disp2 = corr_rescale( disp2, gamma )
disp3 = corr_rescale( disp3, gamma )

clim = np.max(disp2[:pix,:pix]/1.)
climl = np.min(disp2[:pix,:pix]/1.)

plt.imshow( disp[x0:pix,x0:pix], extent=[minres,lim,minres,lim], origin='lower' )
#plt.imshow( corr[:pix,:pix,0], origin='lower' )
plt.xlabel(r'r'+ur'(\u00c5)', fontsize=fs)
plt.ylabel(r'r$^\prime$'+ur'(\u00c5)', fontsize=fs)
plt.colorbar()
plt.clim(climl,clim)
plt.savefig( outpath+outtag+'_theta0_'+label1+'.png')

plt.figure()
plt.imshow( disp2[:pix,:pix], extent=[0,lim,0,lim], origin='lower' )
#plt.imshow( corr[:pix,:pix,0], origin='lower' )
plt.xlabel(r'r'+ur'(\u00c5)', fontsize=fs)
plt.ylabel(r'r$^\prime$'+ur'(\u00c5)', fontsize=fs)
plt.colorbar()
plt.clim(climl,clim)
plt.savefig( outpath+outtag+'_theta0_'+label2+'.png')

#plt.figure()
plt.imshow( disp3, extent=[0,lim,0,lim], origin='lower' )
#plt.imshow( corr[:pix,:pix,0], origin='lower' )
plt.xlabel(r'r'+ur'(\u00c5)')
plt.ylabel(r'r$^\prime$'+ur'(\u00c5)')
plt.colorbar()
plt.clim(climl,clim)
plt.savefig( outpath+outtag+'_theta0_'+label3+'.png')


#clim2 = clim
r1, r2 = 1.4, 1.4 #rmax/2, rmax/2 #
ir = int( r1   * (nr/float(rmax)) )
ir2 = int( r2  * (nr/float(rmax)) )

thpix = 201
thlim = 360.0 * (thpix/float(nth))
asp = 16
gamma2 = 0.5
disp = np.zeros( (nr, nth) )
disp2 = np.zeros( (nr, nth) )
disp3 = np.zeros( (nr,nth) )
disp4 = np.zeros( (nr, nth) )
disp5 = np.zeros( (nr,nth) )
disp6 = np.zeros( (nr,nth) )
for i in np.arange( nr ):
     disp[i,:] = corr[i,i,:]
     disp2[i,:] = corrn[i,i,:]
     disp3[i,:] = corr5[i,i,:]
     disp4[i,:] = corr[ir,i,:]
     disp5[i,:] = corrn[ir,i,:]
     disp6[i,:] = corr5[ir,i,:]

x2th = np.outer( x*x, np.ones(nth) )
#disp = disp*x2th
#disp2 = disp2*x2th
#disp4 = disp4*x2th
#disp5 = disp5*x2th


disp = corr_rescale( disp, gamma2 )
disp2 = corr_rescale( disp2, gamma2 )
disp3 = corr_rescale( disp3, gamma2 )
disp4 = corr_rescale( disp4, gamma2 )
disp5 = corr_rescale( disp5, gamma2 )
disp6 = corr_rescale( disp6, gamma2 )
#disp3 = corr_rescale( disp3*x2th, gamma2 )
clim2 = np.max(disp2/8.)   #/2. - for July 19 draft
clim2l = np.min(disp2/8.)  #/2. - for July 19 draft
minres = 1.0
x0 = int((minres/rmax)*nr)

plt.figure()
plt.imshow( disp[x0:pix,:thpix], extent=[0,thlim,minres,lim], origin='lower', aspect=asp )
plt.ylabel(r'r'+ur'(\u00c5)'+' = '+r'r$^\prime$'+ur'(\u00c5)', fontsize=fs)
plt.xlabel(r'$\theta$ (degrees)', fontsize=fs)
#plt.colorbar()
plt.clim(clim2l,clim2)
plt.savefig( outpath+outtag+'_r=rprime_theta_'+label1+'.png')

plt.figure()
plt.imshow( disp2[x0:pix,:thpix], extent=[0,thlim,minres,lim], origin='lower', aspect=asp )
plt.ylabel(r'r'+ur'(\u00c5)'+' = '+r'r$^\prime$'+ur'(\u00c5)', fontsize=fs)
plt.xlabel(r'$\theta$ (degrees)', fontsize=fs)
#plt.colorbar()
plt.clim(clim2l,clim2)
plt.savefig( outpath+outtag+'_r=rprime_theta_'+label2+'.png')

plt.figure()
plt.imshow( disp4[x0:pix,:thpix], extent=[0,thlim,minres,lim], origin='lower', aspect=asp )
plt.ylabel(r'r'+ur'(\u00c5)', fontsize=fs)
plt.xlabel(r'$\theta$ (degrees)', fontsize=fs)
#plt.colorbar()
plt.clim(clim2l,clim2)
plt.savefig( outpath+outtag+'_r='+'{:.1f}'.format(r1)+'_rprime_theta_'+label1+'.png')

plt.figure()
plt.imshow( disp5[x0:pix,:thpix], extent=[0,thlim,minres,lim], origin='lower', aspect=asp )
plt.ylabel(r'r'+ur'(\u00c5)', fontsize=fs)
plt.xlabel(r'$\theta$ (degrees)', fontsize=fs)
#plt.colorbar()
plt.clim(clim2l,clim2)
plt.savefig( outpath+outtag+'_r='+'{:.1f}'.format(r1)+'_rprime_theta_'+label2+'.png')

plt.figure()
plt.imshow( disp3[x0:pix,:thpix], extent=[0,thlim,minres,lim], origin='lower', aspect=asp )
plt.ylabel(r'r'+ur'(\u00c5)'+' = '+r'r$^\prime$'+ur'(\u00c5)')
plt.xlabel(r'$\theta$ (degrees)')
#plt.colorbar()
plt.clim(clim2l,clim2)
plt.savefig( outpath+outtag+'_r=rprime_theta_'+label3+'.png')


plt.figure()
plt.imshow( disp6[x0:pix,:thpix], extent=[0,thlim,minres,lim], origin='lower', aspect=asp )
plt.ylabel(r'r'+ur'(\u00c5)', fontsize=fs)
plt.xlabel(r'$\theta$ (degrees)', fontsize=fs)
#plt.colorbar()
plt.clim(clim2l,clim2)
plt.savefig( outpath+outtag+'_r='+'{:.1f}'.format(r1)+'_rprime_theta_'+label3+'.png')

#plt.draw()
#plt.show()
#exit()

#zoom in on a region of the r=r' plot
minres = 1.0
x0 = int((minres/rmax)*nr)
maxres = 3.0
x1 = int((maxres/rmax)*nr)
minth = 0.
th0 = int((minth/360.)*nth)
maxth = 180.
th1 = int((maxth/360.)*nth)
plt.figure()
plt.imshow( disp[x0:x1,th0:th1], extent=[minth,maxth,minres,maxres], origin='lower', aspect=asp )
plt.ylabel(r'r'+ur'(\u00c5)'+' = '+r'r$^\prime$'+ur'(\u00c5)', fontsize=fs)
plt.xlabel(r'$\theta$ (degrees)', fontsize=fs)
#plt.colorbar()
plt.clim(clim2l,clim2)
plt.savefig( outpath+outtag+'_r=rprime_theta_unactivated_zoom.png')


plt.figure()
plt.imshow( disp2[x0:x1,th0:th1], extent=[minth,maxth,minres,maxres], origin='lower', aspect=asp )
plt.ylabel(r'r'+ur'(\u00c5)'+' = '+r'r$^\prime$'+ur'(\u00c5)', fontsize=fs)
plt.xlabel(r'$\theta$ (degrees)', fontsize=fs)
#plt.colorbar()
plt.clim(clim2l,clim2)
plt.savefig( outpath+outtag+'_r=rprime_theta_C9_zoom.png')

plt.figure()
plt.imshow( disp4[x0:x1,th0:th1], extent=[minth,maxth,minres,maxres], origin='lower', aspect=asp )
plt.ylabel(r'r'+ur'(\u00c5)'+' = '+r'r$^\prime$'+ur'(\u00c5)', fontsize=fs)
plt.xlabel(r'$\theta$ (degrees)', fontsize=fs)
#plt.colorbar()
plt.clim(clim2l,clim2)
plt.savefig( outpath+outtag+'_r='+'{:.1f}'.format(r1)+'_rprime_theta_unactivated_zoom.png')

plt.figure()
plt.imshow( disp5[x0:x1,th0:th1], extent=[minth,maxth,minres,maxres], origin='lower', aspect=asp )
plt.ylabel(r'r'+ur'(\u00c5)'+' = '+r'r$^\prime$'+ur'(\u00c5)', fontsize=fs)
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
     dr5[i] = corr5[i,i,0]
#plt.figure()
fig, ax = plt.subplots()
xlim = 128
x0=0
p,  = plt.plot( x[:xlim], dr[:xlim] ) #*x[x0:xlim]*x[x0:xlim] )
pa, = plt.plot( x[:xlim], drn[:xlim] ) # *x[x0:xlim]*x[x0:xlim] )
pb, = plt.plot( x[:xlim], dr5[:xlim] ) #*x[x0:xlim]*x[x0:xlim] )
plt.ylabel(r'$\Theta(r,r^\prime,\theta) - G(r)G(r^\prime)$', fontsize=fs)
plt.xlabel(r'$r$ ('+ur'\u00c5'+')', fontsize=fs)
plt.title( 'r = r$^\prime$;'+r'$\theta$'+'= 0', fontsize=fs)
#plt.legend( [p,pa, pb], ['Unactivated FEM 900 frames', 'C9 FEM 900 frames', 'C5 FEM 225 frames' ])
plt.legend( [p,pa, pb], [label1, label2, label3 ], fontsize=fs)
x_bounds = ax.get_xlim()
for i, xv in enumerate(xlines):
     plt.axvline( xv, label=xlabels[i], color='k')
     ax.annotate(s=xlabels[i], xy =(((xv+xlpos-x_bounds[0])/(x_bounds[1]-x_bounds[0])),0.9), xycoords='axes fraction', verticalalignment='right', horizontalalignment='right bottom' , rotation = 270)
plt.savefig( outpath+outtag+'diagonal_line_plot.png')

#plt.figure()
fig, ax = plt.subplots()
xlim = 128
x0=0
p,  = plt.plot( x[:xlim], dr[:xlim] *x[x0:xlim]*x[x0:xlim] ) # *x[x0:xlim]*x[x0:xlim] )
pa, = plt.plot( x[:xlim], drn[:xlim] *x[x0:xlim]*x[x0:xlim] ) # *x[x0:xlim]*x[x0:xlim] )
pb, = plt.plot( x[:xlim], dr5[:xlim]*x[x0:xlim]*x[x0:xlim] ) #*x[x0:xlim]*x[x0:xlim] )
plt.ylabel(r'$\Theta(r,r^\prime,\theta) - G(r)G(r^\prime)$', fontsize=fs)
plt.xlabel(r'$r$ ('+ur'\u00c5'+')', fontsize=fs)
plt.title( 'r = r$^\prime$;'+r'$\theta$'+'= 0', fontsize=fs)
#plt.legend( [p,pa, pb], ['Unactivated FEM 900 frames', 'C9 FEM 900 frames', 'C5 FEM 225 frames' ]
plt.legend( [p,pa, pb], [label1, label2, label3 ], fontsize=fs)
x_bounds = ax.get_xlim()
for i, xv in enumerate(xlines):
     plt.axvline( xv, label=xlabels[i], color='k')
     ax.annotate(s=xlabels[i], xy =(((xv+xlpos-x_bounds[0])/(x_bounds[1]-x_bounds[0])),0.9), xycoords='axes fraction', verticalalignment='right', horizontalalignment='right bottom' , rotation = 270)
plt.savefig( outpath+outtag+'diagonal_line_plot2.png')

#plt.figure()
fig, ax = plt.subplots()
x0 = 12
p,  = plt.plot( x[x0:xlim], dr[x0:xlim] )
pa, = plt.plot( x[x0:xlim], drn[x0:xlim] )
pb, = plt.plot( x[x0:xlim], dr5[x0:xlim] )
plt.ylabel(r'$\Theta(r,r^\prime,\theta) - G(r)G(r^\prime)$', fontsize=fs)
plt.xlabel(r'$r$ ('+ur'\u00c5'+')', fontsize=fs)
plt.title( 'r = r$^\prime$;'+r'$\theta$'+'= 0', fontsize=fs)
#plt.legend( [p,pa, pb], ['Unactivated FEM', 'C9 FEM', 'C5 FEM' ])
plt.legend( [p,pa, pb], [label1, label2, label3 ], fontsize=fs)
x_bounds = ax.get_xlim()
for i, xv in enumerate(xlines):
     plt.axvline( xv, label=xlabels[i], color='k')
     ax.annotate(s=xlabels[i], xy =(((xv+xlpos-x_bounds[0])/(x_bounds[1]-x_bounds[0])),0.9), xycoords='axes fraction', verticalalignment='right', horizontalalignment='right bottom' , rotation = 270)
plt.savefig( outpath+outtag+'diagonal_line_plot_zoom.png')

plt.figure()
fig, ax = plt.subplots()
x0 = 12
p,  = plt.plot( x[x0:xlim], dr[x0:xlim]*x[x0:xlim]*x[x0:xlim] )
pa, = plt.plot( x[x0:xlim], drn[x0:xlim]*x[x0:xlim]*x[x0:xlim] )
pb, = plt.plot( x[x0:xlim], dr5[x0:xlim]*x[x0:xlim]*x[x0:xlim] )
plt.xlabel(r'$r$ ('+ur'\u00c5'+')', fontsize=fs)
plt.ylabel(r'$\Theta(r,r^\prime,\theta) - G(r)G(r^\prime)$', fontsize=fs)
plt.title( 'r = r$^\prime$;'+r'$\theta$'+'= 0', fontsize=fs)
#plt.legend( [p,pa, pb], ['Unactivated FEM', 'C9 FEM', 'C5 FEM' ])
plt.legend( [p,pa, pb], [label1, label2, label3 ], fontsize=fs)
x_bounds = ax.get_xlim()
for i, xv in enumerate(xlines):
     plt.axvline( xv, label=xlabels[i], color='k')
     ax.annotate(s=xlabels[i], xy =(((xv+xlpos-x_bounds[0])/(x_bounds[1]-x_bounds[0])),0.9), xycoords='axes fraction', verticalalignment='right', horizontalalignment='right bottom' , rotation = 270)
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

pr  = np.array([ 13, 15, 17, 13, 13, 13, 13, 21, 22, 27, 13, 13] )
pr2  = np.array([ 13, 15, 17, 17, 21, 22, 15, 21, 22, 27, 30, 33] )
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
     lines3.append( corr5[ pr[i], pr2[i], :] )

for i in np.arange( len(pr) ):
     plt.figure()
     p, = plt.plot( th[:nth/2], lines[i][:nth/2]  )
     pa, = plt.plot( th[:nth/2], lines2[i][:nth/2] )
     pb, = plt.plot( th[:nth/2], lines3[i][:nth/2] )
     # plt.axis([0,360,y0a,y1a])
     plt.ylabel(r'$\Theta(r,r^\prime,\theta) - G(r)G(r^\prime)$', fontsize=fs)
     plt.xlabel(r'$\theta$ (degrees)', fontsize=fs)
     leg1 = plt.Rectangle((0, 0), 0, 0, alpha=0.0)
     # plt.legend( [leg1], ['r = 2.5 '+ur'\u00c5'+'\n'  + r'r$^\prime$ = 2.5 '+ur'\u00c5'],
     #            handlelength = 0)
     plt.title( 'r = '+'{:.1f}'.format(prx[i])+' '+ur'\u00c5'+'  '+r'r$^\prime$'+'{:.1f}'.format(pr2x[i])+' '+ur'\u00c5', fontsize=fs)
     plt.legend( [p,pa, pb], [label1, label2, label3 ], fontsize=fs)
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
