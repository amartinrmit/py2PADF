
import numpy as np
import matplotlib.pyplot as plt
import os
import struct
import array
import sys



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

def read_bin( fname, swapbyteorder=0 ):

     size = os.path.getsize(fname)
     #print "sizes:", size, size/8
     b = array.array('f')
     f = open(fname, "r")
     b.fromfile( f, size/4 )
     f.close();
     l = b.tolist()
     
     n = len(l)
     nx = int(round(n**(1.0/2.0)))
     print "nx", nx
     output = np.array(l).reshape( nx, nx)
     if swapbyteorder == 1: output = output.newbyteorder()
     return output

# shift - a 2D version of numpy's roll
def array_shift(array,xshift=0,yshift=0):
	array = np.roll(array,xshift,0)
	array = np.roll(array,yshift,1)
	return array

## make an array with a circle set to one
def circle(nx, ny, rad=None, cenx=None, ceny=None, invert=0 ): 
	
    # set defaults
    if rad is None: rad = np.min([nx,ny])/2
    if cenx is None: cenx = nx/2
    if ceny is None: ceny = ny/2

    # define the circle
    x = np.outer(np.arange(nx),np.ones(ny)) - nx/2
    y = np.outer(np.ones(nx),np.arange(ny)) - ny/2
   
    dist = np.sqrt(x**2 + y**2)
    a = np.zeros([nx,ny])
    icirc = np.where(dist <= rad)
    a[icirc] = 1.

    if (invert==1):
	    a = abs(1-a)

    out = array_shift(a, -nx/2, -ny/2)
    out = array_shift(out, cenx, ceny)
    return out
#end circle


def radial_plot( data, mask, rmin, rmax, rbins ):

    nx, ny = data.shape
    x = np.outer( np.arange( nx) - nx/2, np.ones( ny ) )
    y = np.outer( np.ones( nx ), np.arange( ny ) - ny/2 )
    r = np.sqrt( x*x + y*y )
    #r = np.roll( np.roll( r, -nx/2, 0), -ny/2, 1)

    hist = np.zeros( rbins )
    hist2 = np.zeros( rbins )
    rstep = (rmax - rmin ) / float(rbins)

    for i in np.arange( rbins ):
        ri = rmin + i*rstep
        
        pix = np.where( (r >= ri)*(r<(ri+rstep))  )
        
        msum = np.sum( mask[pix] )
        if msum > 0.5:
            hist[i] = np.sum( data[pix] ) / msum
            hist2[i] = np.sum( data[pix]*data[pix] ) / msum

    return hist, hist2


def centrosymmetry_metric( data, mask, rmin, rmax ):

    nx, ny = data.shape
    x = np.outer( np.arange( nx) - nx/2, np.ones( ny ) )
    y = np.outer( np.ones( nx ), np.arange( ny ) - ny/2 )
    r = np.sqrt( x*x + y*y )
    ir = np.where( (r > rmin)*(r<rmax) )

    rotdata = np.rot90( data, 2 )
    rotmask = np.rot90( mask, 2 )
    tmask = mask*rotmask

    #error = np.sum( rotdata[ir] * data[ir] ) / np.sum( mask[ir] * rotmask[ir] )
    #error *= 1.0 / np.sqrt( np.sum( mask[ir]*rotmask[ir]*rotdata[ir])*np.sum( data[ir]*mask[ir]*rotmask[ir]) )
    error = np.sum( np.abs( (rotdata[ir] - data[ir])*tmask[ir])  ) / np.sum( tmask[ir] )

    return error






etype = 'fem'
stype =  'unactivated' #     'C5' #
dcode = '13.04.29'
asize = 512
stag = "_240519_centest" 
#path  = "/Users/e38496/Work/Research/Experiments/Monash/2016/Carbon_dec16/Results/"+stype+"/"+etype+"/data_"+str(asize)+"_"+dcode+stag+"/"
path  = "/Users/e38496/Work/Research/Experiments/Monash/2016/Carbon_dec16/Results/datasum/"

#tag = "c_"+stype+"_"+etype+"_"+dcode+"_"+str(asize)+"_0_1"+stag
#fname = path + tag + "_diffraction_sum.dbin"
#fname = path + tag + "_0_diffraction.dbin"
tag = dcode+"_Spectrum_image_1_0_900_sum.dbin"
fname = path + tag ## + "_diffraction_cropped.dbin"

data = read_dbin( fname ) 

maskpath = "/Users/e38496/Work/Research/Experiments/Monash/2016/Carbon_dec16/masks/"
maskname = etype+"_"+stype+"_"+dcode+"_mask_inverted.bin"
#maskname = path + tag + "_mask_cropped.dbin"
mask = read_bin(maskpath+ maskname )
mask *= 1.0 / np.max(mask)
#cenx = 251
#ceny = 251
n = asize
nx = asize


#mshifted = np.roll( np.roll( mask, cenx-asize, 0), cenx-asize, 1
mshifted = np.roll( np.roll( mask, 0, 0), 0, 1)
maskcrop = np.zeros( (asize,asize) )
maskcrop[:,:] = mask / np.max(mask) #mshifted[n-(nx/2):n+(nx/2),n-(nx/2):n+(nx/2)] / np.max(mshifted)
#maskcrop = 1.0 - maskcrop
print "maskcrop max min", np.max(maskcrop), np.min(maskcrop), np.max(mask)
#plt.imshow( maskcrop )
#plt.draw()
#plt.show()

circ = circle(nx, nx, rad=5, cenx=nx/2, ceny=nx/2, invert=1 )
circ2 = circle(nx, nx, rad=75, cenx=nx/2, ceny=nx/2, invert=0 )
#maskcrop *= circ*circ2


rmin, rmax, rbins = 1, asize/2, asize/2
xstart, xend = -10, 10
ystart, yend = -10, 10

plt.ion()
err = np.zeros( (xend-xstart, yend-ystart)) + 1e10
minerr = 1e10
for ix in np.arange( xend-xstart):
    print ix, ix+xstart
    for iy in np.arange( yend-ystart):

        tmp = np.roll( np.roll( data, ix+xstart,0), iy+ystart, 1)
        tmp_mask = np.roll( np.roll( maskcrop, ix+xstart,0), iy+ystart, 1) * circ*circ2
        rad, var = radial_plot( tmp*tmp_mask, tmp_mask, rmin, rmax, rbins )
        
        # current
        #err[ix,iy] = np.average( var[60:100] - rad[60:100]*rad[60:100] )
       
        err[ix,iy] = centrosymmetry_metric( tmp, tmp_mask, 50, 100)
        if (err[ix,iy] < minerr):
             print "ix, iy, current min :", ix+xstart, iy+ystart, np.min(err)
             minerr = err[ix,iy]

        # err[ix,iy] = np.max( var - rad*rad) #np.average(( var  - rad*rad))
        # err[ix,iy] = np.max( rad) #np.average(( var  - rad*rad))
    
        # print ix, iy, err[ix,iy]
        #plt.gcf().clear()
        #plt.plot( var[40:] - rad[40:]*rad[40:] )
        #plt.imshow( tmp )
        #plt.draw()

plt.ioff()

print "min_error, ixmin, iymin :", np.min(err), np.where( err == np.min(err))
print "max_error, ixmax, iymax :", np.max(err), np.where( err == np.max(err))

imin = np.where( err == np.min(err) )
imax = np.where( err == np.max(err) )
#corrected = np.roll( np.roll( data*maskcrop, imin[0][0]+xstart,0), imin[1][0]+ystart, 1)
shift = imin
print imin
print "final shift :", shift[0][0]+xstart, shift[1][0]+ystart
corrected = np.roll( np.roll( data*maskcrop, (shift[0][0]+xstart),0), (shift[1][0]+ystart), 1)

plt.figure()
plt.imshow( np.abs(data)**0.3, origin='lower' )
plt.title( "Diffraction pattern" )
plt.figure()
plt.imshow( maskcrop, origin='lower' )
plt.title( "Mask" )
plt.figure()
#plt.plot( var )
#plt.figure()
plt.imshow( err, extent=[xstart,xend,ystart,yend], origin='lower' )
plt.title( "Centrosymmetric error" )
plt.xlabel( "x pixel shift" )
plt.ylabel( "y pixel shift" )
plt.figure()
plt.imshow( corrected  , origin='lower' )
plt.draw()
plt.show()
