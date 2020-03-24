#!/usr/bin/python


import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import struct
import array


def read_dbin( fname):
        dlist= []
        with open(fname, "rb") as f:  ##does this file get closed??
            bytes = f.read(8)
            while bytes:
                number, = longBitsToFloat(bytes)
                dlist.append( number )
                bytes = f.read(8)
        array = np.array(dlist)
        return array

def read_dbin_2d( fname):       
        array = read_dbin( fname )
	print "array.size ", array.size, np.sqrt(array.size)
        array = array.reshape(np.sqrt(array.size), \
                                  np.sqrt(array.size))
        return array

def write_dbin( fname, data ):

	f = open( fname, "wb")
	fmt='>'+'d'*data.size
	bin = struct.pack(fmt, *data.flatten()[:] )
	f.write( bin )
	f.close()

def write_dbin2( fname, data ):
     b = array.array('d')
     b.fromlist( data.flatten().tolist() )
     f = open(fname, "w")
     b.tofile( f )
     f.close();


def floatToRawLongBits( value):
	return struct.unpack('Q', struct.pack('d', value))[0]

def longBitsToFloat( bits):
	# fl = struct.unpack('q', bits)
	fl = struct.unpack('<d', bits)
        return fl


# shift - a 2D version of numpy's roll
def array_shift(array,xshift=0,yshift=0):
	array = np.roll(array,xshift,0)
	array = np.roll(array,yshift,1)
	return array

## make a 2D array with a gaussian
def make_gaussian(nx, ny, rad=None, rady=-1., cenx=None, ceny=None, invert=0, norm=False, power=2 ): 
    #set defaults
    if rad is None: rad = np.min(nx,ny)/2
    if cenx is None: cenx = nx/2
    if ceny is None: ceny = ny/2
    radsq = rad**power
    if rady == -1.:
        radysq = radsq
    else:
        radysq = rady**power

    # define the circle
    x = np.outer(np.arange(0-nx/2,nx-nx/2,1),np.ones(ny))
    #print x.size, x.shape
    y = np.outer(np.ones(nx),np.arange(0-ny/2,ny-ny/2,1))
    #print y.size, y.shape

    a = np.zeros([nx,ny])
    #a = np.exp(-(x**2)/radsq  - ( y**2)/radysq)
    a = np.exp(-(x**power)/radsq  - ( y**power)/radysq)
    a[ nx/2, ny/2 ] = 1.0

    a = array_shift(a,cenx-nx/2,ceny-ny/2)

    # normalise if required
    if norm == True: a *= 1./np.sum(a)
    
    return a
#end make_gaussian


path = "/scratch/amartin1/Work/Results/padf/lammps/si_anneal_071215_padf/"
fname = "aSi071215_900__1_diffraction_lowq.dbin"


#path = "/home/amartin1/Work/Results/padf/poly500/"
#path  = "/scratch/amartin1/Work/Results/padf/water/onemol/"
#path = "/home/amartin1/Work/Results/padf/polymm/"
#fname = "water_onemol_scaledtest_1_diffraction_lowq.dbin"
#fname = "blaq2.dbin"
#fname = "mm_bl_l2.dbin"

dp = read_dbin_2d( path+fname )
g = make_gaussian( dp.shape[0], dp.shape[1], dp.shape[0]/4, dp.shape[1]/4, dp.shape[0]/2, dp.shape[1]/2 )
autoc = np.fft.fft2( dp * g )

print "max min", np.max( dp ), np.min(dp)

#plt.imshow( np.abs(autoc)**0.3 )
plt.imshow( np.real(dp)**0.3 )
#sp.misc.imsave( path+fname[:-5]+'.png', np.abs(dp)**0.3 )
plt.draw()
plt.show()
