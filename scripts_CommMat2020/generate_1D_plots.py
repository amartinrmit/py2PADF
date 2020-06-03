
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import imread
import os
import glob
import sys

class rplot:

    def __init__( self, mask=None, rmin=0, rmax=100, rbins=100, cenx=0.0, ceny=0.0,
                  s=None):

        self.mask = mask
        self.rmin = rmin
        self.rmax = rmax
        self.rbins = rbins
        self.cenx = cenx
        self.ceny = ceny
        self.s = s

    def setup( self ):
        
        nx, ny = self.s[0], self.s[1]

        x = np.outer( np.arange( nx), np.ones( ny ) )  - cenx
        y = np.outer( np.ones( nx ), np.arange( ny ) ) - ceny
        r = np.sqrt( x*x + y*y )
        self.hist = np.zeros( rbins )
        self.hist2 = np.zeros( rbins )
        rstep = (rmax - rmin ) / float(rbins)

        self.pixlist = []
        self.msum = []
        for i in np.arange( self.rbins ):
            ri = self.rmin + i*rstep
            
            self.pixlist.append( np.where( (r >= ri)*(r<(ri+rstep))  ) )
            
            self.msum.append( np.sum( self.mask[self.pixlist[i]] ) )

        
    def plot( self, data ):

        for i in np.arange( rbins ):
            if self.msum[i] > 0.5:
                self.hist[i] = np.sum( data[self.pixlist[i]] ) / self.msum[i]
        

def radial_plot( data, mask, rmin, rmax, rbins, cenx=0.0, ceny=0.0 ):

    nx, ny = data.shape

    x = np.outer( np.arange( nx), np.ones( ny ) )  - cenx
    y = np.outer( np.ones( nx ), np.arange( ny ) ) - ceny
    r = np.sqrt( x*x + y*y )
    #r = np.roll( np.roll( r, -nx/2, 0), -ny/2, 1)

#    plt.imshow( r )
#    plt.figure() 

    hist = np.zeros( rbins )
    hist2 = np.zeros( rbins )
    rstep = (rmax - rmin ) / float(rbins)

    for i in np.arange( rbins ):
        ri = rmin + i*rstep
        
        pix = np.where( (r >= ri)*(r<(ri+rstep))  )
        
        msum = np.sum( mask[pix] )
        if msum > 0.5:
            hist[i] = np.sum( data[pix] ) / msum
         #   hist2[i] = np.sum( data[pix]*data[pix] ) / msum

    return hist, hist2

datafolder = "/Volumes/DataStore1/AS_SAXS_Feb18_data/0p6_saxs/"
path = datafolder+"images/"
#tag = "run15_plate18"
tag = "plate4"


datalist = sorted( glob.glob( path+tag+"*.tif" ), key=os.path.getmtime )
data = imread( datalist[0]) 
print datalist[0]

#maskname = "/Users/e38496/Work/Research/Experiments/AS/SAXS/SAXS_feb18/image_2D_sums/plate4mask_dilated.tif"
maskname = "/Users/e38496/Work/Research/Experiments/AS/SAXS/SAXS_feb18/image_2D_sums/run15_plate18_nstart_0_n_300_mask_dilated2.tif"
mask = imread( maskname ).astype(np.float)
print np.max(mask)
mask *= 1.0/np.max(mask)


#x0, y0 = 672.803, 335.895
x0, y0 = 434.191,  335.369

s = data.shape
shiftx = s[0]/2 - y0
shifty = x0 - s[1]/2
cenx = s[0]/2 + shiftx
ceny = s[1]/2 + shifty
print "cenx, ceny:", cenx, ceny, shiftx, shifty
rbins = 880
rmin = 18
rmax = rmin + rbins

#debug
#sys.exit()


data_shift = np.roll( np.roll( data, int(shiftx), 0), int(shifty), 1 )

plt.imshow( data )
#plt.draw()
#plt.show()
plt.figure() 

rp = rplot( mask, rmin, rmax, rbins, cenx, ceny, s )
rp.setup()


#intensity1d, tmp  = radial_plot( data, mask, rmin, rmax, rbins, cenx, ceny )
#rp.plot( data)
#intensity1d = rp.hist

#plt.plot( intensity1d)

#
# Check against ScatterBrain 1d plot
#
#rawpath = datafolder+"raw_dat/"
#rawtag = "plate4"
#raw = np.loadtxt( rawpath+rawtag+"_0010.dat") 
# q = raw[:,0]

#print raw.shape
#plt.figure()
#p2, = plt.plot( raw[:,1] )

#print raw[:,0]

# incorrect scatterbrain calibration
#wavelength = 0.619921 #1.033216 #e-10
#detector_z = 0.962478     #0.959563 #
#pixelwidth = 0.000172 
#rstart = 16.5  

# what I think it should be based on LCP / LCP2 log files
wavelength = 1.033216 
detector_z = 0.959563 
pixelwidth = 0.000172   
rstart = rmin

qindices = rstart + np.arange(rbins )
q = 2.0 * np.pi * (2.0/wavelength) * np.sin( np.arctan( qindices*(pixelwidth/detector_z) ) / 2.0 )
#print raw[:5,0]
#print q[0:5]
#print q[0:5]*raw[0,0]/q[0]

#plt.plot( q )
#plt.plot( raw[:,0] )
#plt.draw()
#plt.show()
#sys.exit()


#
# Generate 1D data for all the run of interest
#
outpath = datafolder + "raw_postexpt_TEST/"
if not os.path.exists( outpath ):
    os.makedirs( outpath )

print "Number of files in run :", len(datalist)
count = 0
for f in datalist:
    
    print "Converting to 1D: ", count, "/", len(datalist), ":", f
    count += 1
    # read in 2d data
    data = imread( f )

    # calculate 1D data
    rp.plot( data)
    
    output = np.array( [q, rp.hist] ).transpose()

    # write out to file
    np.savetxt( outpath+os.path.basename(f)[:-4]+".dat", output )
