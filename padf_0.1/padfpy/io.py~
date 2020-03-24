
import numpy as np
import amser as ser
from scipy.ndimage import imread
import os
import array
import struct
import h5py

def read_dbin( fname, swapbyteorder=0, nx=-1, ny=-1 ):

    size = os.path.getsize(fname)
    # print "sizes:", size, size/8
    b = array.array('d')
    f = open(fname, "r")
    b.fromfile( f, size/8 )
    f.close();
    l = b.tolist()
    
    n = len(l)
    if (nx==-1) and (ny==-1):
        nx = int(round(n**(1.0/2.0)))
    if ny==-1:
        ny = nx
    #print "nx", nx
    output = np.array(l).reshape( nx, ny)
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

def read_correlation( fname, swapbyteorder=0 ):

    size = os.path.getsize(fname)
    # print "sizes:", size, size/8
    b = array.array('d')
    f = open(fname, "r")
    b.fromfile( f, size/8 )
    f.close();
    l = b.tolist()
    output = np.array(l)
    
    if swapbyteorder == 1: output = output.newbyteorder()
    return output

def read_correlation_single( fname, swapbyteorder=0 ):
    
    size = os.path.getsize(fname)
    # print "sizes:", size, size/8
    b = array.array('f')
    f = open(fname, "r")
    b.fromfile( f, size/4 )
    f.close();
    l = b.tolist()
    output = np.array(l)
    
    if swapbyteorder == 1: output = output.newbyteorder()
    return output

def write_dbin( fname, data ):
        
    f = open( fname, "wb")
    fmt='<'+'d'*data.size
    bin = struct.pack(fmt, *data.flatten()[:] )
    f.write( bin )
    f.close()

def h5read(filename,field="data/data1"):
     h5file = h5py.File(filename,"r")
     #print field
     h5data = h5file[field]
     image = h5data[...]
     h5file.close()
     return image


def read_image( fname, swapbyteorder=0, ser_index=0, nx=-1, ny=-1, h5field="/data/data1" ):

    fname = fname.rstrip('\n')
    fname = fname.rstrip('\r\n')
    print "DEBUG <io.py> fname", fname
    if fname[-5:]==".dbin":
        image = read_dbin( fname, swapbyteorder, nx=nx, ny=ny )
    elif fname[-4:]==".bin":
        image = read_bin( fname, swapbyteorder ).astype( np.float64) 
    elif fname[-4:]==".ser":
        image = read_ser( fname, ser_index )
    elif fname[-3:]==".h5":
        image = h5read( fname, h5field )
    else:
        try:
            image = imread( fname )
        except:
            print "error reading image - unknown file location or extension. fname = ", fname
            
    return image.astype( np.float64 )

    
def read_ser( fname, ser_index=0 ):

    data, calp, adl, sh = ser.read_ser( fname, [ser_index] )  
    image = data.dlist[0]
    return image


def read_ser_npatterns( fname ):

    data, calp, adl, sh = ser.read_ser( fname, [1] )
    return sh.ValNumElem


def read_ser_arraysize( fname ):

    data, calp, adl, sh = ser.read_ser( fname, [1] ) 
    #print "DEBUG <io.py; read_ser_arraysize> arraysize", calp.ArraysizeX, data.dlist[0].shape
    return calp.ArraysizeX

def get_array_dimensions_from_file( fname ):

    image = read_image( fname )
    s = image.shape()
    return s
