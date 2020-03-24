#!/usr/bin/python


import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import struct
import array
import time


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

def write_dbin( fname, data ):

	f = open( fname, "wb")
	fmt='<'+'d'*data.size
	bin = struct.pack(fmt, *data.flatten()[:] )
	f.write( bin )
	f.close()


path = "/home/amartin1/Work/Results/padf/crosstest/"
tag = "crosstest"

#
#----------------------------------------------------------------
# SET UP THE OBJECT
#
#

#
# Parameters
#
nx = 256

#
#  generate an array of coordinates. e.g. random polyhedra
#

#
# write a pdb file (at least the atom lines )
#

#
# Or a cross
#
nlength = nx / 4
nwid = 6
object = np.zeros( (nx,nx) )
object[ nx/2-nwid/2:nx/2+nwid/2, nx/2-nlength/2:nx/2+nlength/2] = 1.0
object[ nx/2-nlength/2:nx/2+nlength/2, nx/2-nwid/2:nx/2+nwid/2] = 1.0

#
# display the object
#
#plt.imshow( object )
#plt.draw()
#plt.show()


#
#----------------------------------------------------------------
# Calculate the diffraction pattern
#
#

#
# dp parameters
#
wl = 2e-10
pixel_width = 2.0e-5
detector_z = 1e-2
cx = nx/2
cy = nx/2

#
# Calculate dp
#
fobj = np.fft.fft2( object )
dp = np.abs(fobj)**2
dp = np.roll( np.roll( dp, cx, 0 ), cy, 1 )

#
# output the diffraction pattern
#
dpoutname = path+"dp.bin"
write_dbin( dpoutname, dp )

#
# display the diffraction pattern
#
#plt.imshow( dp**0.3 )
#plt.draw()
#plt.show()


#
# display the autocorrelation function
#
autoc = np.fft.ifft2( dp )
autoc = np.roll(np.roll( autoc, nx/2, 0), nx/2, 1 )
plt.imshow( np.real(autoc) )
plt.draw()
plt.show()

sys.exit()

#
#----------------------------------------------------------------
# Calculate the correlation
#
#

#
# Create a configuration file
#
corrfname = path+"corrconfig.txt"
fcorr = open(corrfname,'w')
fcorr.write("input = "+dpoutname+'\n' )
fcorr.write("outpath = "+path+'\n')
fcorr.write("tag = "+tag+"_padfcorr"+'\n')
fcorr.write("wavelength =  "+str(wl)+'\n'  )
fcorr.write("pixel_width =  "+str(pixel_width)+'\n'  )
fcorr.write("detector_z =  "+str(detector_z)+'\n'  )
fcorr.write("cx =  "+str(cx)+'\n'  )
fcorr.write("cy =  "+str(cy)+'\n'  )
fcorr.close()

#
# run padfcorr
#
print "/home/amartin1/Work/codes/C-code/padfcorr/padfcorr "+corrfname

start = time.clock()
#os.system("/home/amartin1/Work/codes/C-code/padfcorr/padfcorr "+corrfname) 
print "padfcorr took :", time.clock() - start, " seconds"

#
# read corr back in and display
#
cname = path+tag+"_padfcorr_correlation.dbin"
corr = read_correlation( cname , 0)
corr = corr.reshape( 128, 128, 804 )

#plt.imshow( corr[12,:,:]**0.3 )
#plt.draw()
#plt.show()


#sys.exit()

#
#----------------------------------------------------------------
# Calculate the padf
#
#

#
# Parameters
#
nl = 10
nth = 804
nr = 128
nq = 128
qmax = 1.24981e9
rmax = 1.02415e-7

#
# Create a configuration file
#
pfname = path+"padf_config.txt"
fcorr = open(pfname,'w')
fcorr.write("correlationfile = "+path+tag+"_padfcorr_correlation.dbin"+'\n' )
fcorr.write("outpath = "+path+'\n')
fcorr.write("tag = "+tag+"_padf"+'\n')
fcorr.write("wavelength =  "+str(wl)+'\n'  )
fcorr.write("nthq =  "+str(nth)+'\n'  )
fcorr.write("nq =  "+str(nq)+'\n'  )
fcorr.write("nr =  "+str(nr)+'\n'  )
fcorr.write("nl =  "+str(nl)+'\n'  )
fcorr.write("qmax =  "+str(qmax)+'\n'  )
fcorr.write("rmax =  "+str(rmax)+'\n'  )
fcorr.close()

#
# run padf
#
#os.system("gdb --args /home/amartin1/Work/codes/C-code/padf/padf "+pfname) 
start = time.clock()
os.system("/home/amartin1/Work/codes/C-code/padf/padf "+pfname) 
print "padf took :", time.clock() - start, " seconds"
