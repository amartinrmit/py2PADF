

import os
import array
import numpy as np
import scipy.ndimage as sdn
import matplotlib.pyplot as plt


def polar_plot( data, nr, nth, rmin, rmax, thmin, thmax, cenx, ceny, submean=False ):

        # r and theta arrays
        rarr = np.outer( np.arange(nr)*(rmax-rmin)/float(nr) + rmin, np.ones(nth) )
        tharr = np.outer( np.ones(nr), np.arange(nth)*(thmax-thmin)/float(nth) + thmin)
        
        newx = rarr*np.cos( tharr ) + cenx
        newy = rarr*np.sin( tharr ) + ceny
        
        newdata = sdn.map_coordinates( data, [newx.flatten(), newy.flatten()], order=3 )

        out = newdata.reshape( nr, nth )
        if submean == True:
            out = self.polar_plot_subtract_rmean( out  )

        return out


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
     print "nx", nx, n
     output = np.array(l).reshape( nx, nx)
     
     if swapbyteorder == 1: output = output.newbyteorder()
     return output

def image_rescale( image, gamma ):
    output = image*0.0
    iup = np.where( image>0.0)
    idown = np.where( image<0.0)
    output[iup] = image[iup]**gamma
    output[idown] = -np.abs(image[idown])**gamma
    return output

class Formatter(object):
    def __init__(self, im):
        self.im = im
    def __call__(self, x, y):
        z = self.im.get_array()[int(y), int(x)]
        return 'x={:.01g}, y={:.01g}, z={:.01g}'.format(x, y, z)




#path = "/Users/e38496/Work/Research/Experiments/Monash/2016/Carbon_dec16/Results/C9/fem/data_12.25.40/" 
#fname = "carbon_dec16_C9_12.25.40diffraction_1.dbin"

path = "/Users/e38496/Work/Research/Experiments/Monash/2016/Carbon_dec16/Results/datasum/" 
fname = "11.20.24_Spectrum_image_1_0_400_sum.dbin"

data = read_dbin( path+fname )
print data.shape
rescaled = image_rescale( data, 0.1 )

s = data.shape

shiftx = 1
shifty = 2
shifted = np.roll( np.roll( data, shiftx, 0), shifty, 1)

polar = polar_plot( shifted, s[0], 402, 0, s[0]/2, 0, 2*np.pi, s[0]/2, s[1]/2 )

fig, ax = plt.subplots()
im = ax.imshow( polar[50:,:] , interpolation='none', origin='lower') #, vmin=-20, vmax=20)
ax.format_coord = Formatter(im)
plt.show()
