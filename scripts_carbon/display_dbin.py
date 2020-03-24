
import numpy as np
import matplotlib.pyplot as plt
import os
import array



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

def read_bin( fname, swapbyteorder=0, nx=-1, ny=-1 ):

    size = os.path.getsize(fname)
    # print "sizes:", size, size/8
    b = array.array('f')
    f = open(fname, "r")
    b.fromfile( f, size/4 )
    f.close();
    l = b.tolist()
    
    n = len(l)
    if (nx==-1) and (ny==-1):
        nx = int(round(n**(1.0/2.0)))
    if ny==-1:
        ny = nx
    output = np.array(l).reshape( nx, ny)
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


#path = "/Volumes/DataStore1/Data/Monash/Gold2018/Gold_5um_460_20ms/dbin_1000/"
#fname = "Gold_5um_460_20ms0.dbin"
#path = "/Users/e38496/Work/Research/Experiments/Monash/2016/Carbon_dec16/masks/" 
#fname = "fem_C9_12.25.40_mask_inverted.bin"

path = "/Users/e38496/Work/Research/Experiments/Monash/2016/Carbon_dec16/Results/C9/fem/data_12.25.40_recen/" 
fname = "carbon_dec16_C9_12.25.40_recen_diffdiffraction_1.dbin"


data = read_dbin( path+fname )
print data.shape
rescaled = image_rescale( data, 0.25 )

fig, ax = plt.subplots()
im = ax.imshow(rescaled, interpolation='none', origin='lower') #, vmin=-20, vmax=20)
ax.format_coord = Formatter(im)
plt.show()
