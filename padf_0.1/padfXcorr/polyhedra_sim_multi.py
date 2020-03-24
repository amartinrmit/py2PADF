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

#
# Random rotations references
#
# http://planning.cs.uiuc.edu/node198.html
#
# K. Shoemake. Uniform random rotations.
# In D. Kirk, editor, Graphics Gems III, pages 124-132. Academic, New York, 1992. 
#
# random sampling of rotation group quaternion: u1, u2, u3 uniformly distributed in [0,1]
# h = [ sqrt(1-u1)*sin(2pi*u2), sqrt(1-u1)*cos(2pi*u2), sqrt(u1) sin(2pi*u3), sqrt(u1)*cos(2pi*u3) ]
#

def random_rotation():

     u1, u2, u3 = np.random.rand(3)
     h = [     np.sqrt(1-u1)*np.sin(2*np.pi*u2), np.sqrt(1-u1)*np.cos(2*np.pi*u2), \
               np.sqrt(u1)*np.sin(2*np.pi*u3),   np.sqrt(u1)*np.cos(2*np.pi*u3) ]
     th = 2.0*np.arccos( h[0] )
     rx = h[1] / np.sin( th / 2.0 )
     ry = h[2] / np.sin( th / 2.0 )
     rz = h[3] / np.sin( th / 2.0 )

     #print "length check :", np.sqrt( rx*rx + ry*ry + rz*rz )
     #print "length check :", np.sqrt( h[0]*h[0] + h[1]*h[1] + h[2]*h[2] + h[3]*h[3] )

     return rx, ry, rz, th

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


def quaternion_multiply( q1, q2 ):

     qout = np.zeros( 4 )
     qout[0] = q1[0]*q2[0]
     for i in np.arange(3)+1:
          qout[0] += -q1[i]*q2[i]
          
     for i in np.arange(3)+1:
          qout[i] = q1[0]*q2[i] + q2[0]*q1[i]
          
          
     qout[1] +=  q1[2]*q2[3] - q1[3]*q2[2]
     qout[2] += -q1[1]*q2[3] + q1[3]*q2[1] 
     qout[3] +=  q1[1]*q2[2] - q1[2]*q2[1] 

     return qout


def rotate_vector( vect, axis, angle):

     alen = np.sqrt(axis[0]*axis[0] + axis[1]*axis[1]
                    + axis[2]*axis[2])
     for i in np.arange(3):
          axis[i] = axis[i]/alen
          
     #
     #  Define the quaternian and conjugate for rotation
     #
     q = np.zeros( 4 )
     q[0] = np.cos(angle/2.0)
     q[1] = axis[0]*np.sin(angle/2.0)
     q[2] = axis[1]*np.sin(angle/2.0)
     q[3] = axis[2]*np.sin(angle/2.0)

     qstar = np.zeros( 4 )
     qvect = np.zeros( 4 )
     qstar[0] = q[0]
     qvect[0] = 0.0
     for i in np.arange(3)+1:
          qstar[i] = -q[i]
          qvect[i] = vect[i-1]
  

     qtemp = quaternion_multiply( qvect, qstar )
     qout  = quaternion_multiply( q    , qtemp )

     vout = np.zeros( 3 )
     for i in np.arange(3):
          vout[i] = qout[i+1];

     return vout


#
#----------------------------------------------------------------
# MAIN CODE BEGINS HERE
#----------------------------------------------------------------
#
#
path = "/scratch/amartin1/Work/Results/padf/incoh-150615/"
tag = "polyhedra"
np.random.seed()

object_flag = 1
dp_flag = 1
g_flag = 1
mm_flag = 0
corr_flag = 1
padf_flag = 1



#
#----------------------------------------------------------------
# SET UP THE OBJECT
#
#
#
# Parameters
#
nx = 128
nth = 200
a = 1.54
elem = "C"
temp = 1
nmol = 5
npatterns = 50000
ninc = 1         # extra option to sum patterns incoherently
fov_radius = 20.0

if object_flag == 1:

     #
     #  generate an array of coordinates. e.g. random polyhedra
     #
     x = [ a / np.sqrt(3), -a/(2.0*np.sqrt(3)), -a/(2.0*np.sqrt(3)), 0.0                ]
     y = [ 0,              -a/2.0,              a/2.0,               0.0                ]
     z = [ 0,              0,                   0,                   np.sqrt(2.0/3.0)*a ]
     # print "distance check :", np.sqrt( (x[0]-x[2])**2 + (y[0]-y[2])**2 + (z[0]-z[2])**2 )
     # print "distance check :", np.sqrt( (x[2]-x[3])**2 + (y[2]-y[3])**2 + (z[2]-z[3])**2 )


     
     for ipat in np.arange(npatterns*ninc):

          #
          # write a pdb file (at least the atom lines )
          #
          pdbname = path+"pdb_polyhedra_"+str(ipat)+".txt"
          fcorr = open(pdbname,'w')

          for imol in np.arange(nmol):

               rx, ry, rz, th = random_rotation()
               sx, sy, sz, th = random_rotation()

               xn, yn, zn = [], [], []
               for i in np.arange( len(x) ):
               
                    v = np.array([x[i],y[i],z[i]])
                    xt, yt, zt = rotate_vector( v, np.array([rx,ry,rz]), th )
                    xn.append( xt + sx*fov_radius)
                    yn.append( yt + sy*fov_radius)
                    zn.append( zt + sz*fov_radius)
                    
                    #print "ATOM  "+str(i+1).ljust(5)+" "+str(elem).rjust(3)+" "+"NNN A "\
                    #    +str(temp).rjust(4)+"    "+"{:.3f}".format(xn[i]).rjust(8)\
                    #    +"{:.3f}".format(yn[i]).rjust(8)+"{:.3f}".format(zn[i]).rjust(8)\
                    #    +"{:.2f}".format(temp).rjust(6)+"{:.2f}".format(temp).rjust(6)\
                    #    +"          "+elem.rjust(2)
          
                    fcorr.write("ATOM  "+str(i+1).ljust(5)+" "+str(elem).rjust(3)+" "+"NNN A "\
                                     +str(temp).rjust(4)+"    "+"{:.3f}".format(xn[i]).rjust(8)\
                                     +"{:.3f}".format(yn[i]).rjust(8)+"{:.3f}".format(zn[i]).rjust(8)\
                                     +"{:.2f}".format(temp).rjust(6)+"{:.2f}".format(temp).rjust(6)\
                                     +"          "+elem.rjust(2)+'\n' )
          
          fcorr.close()

# sys.exit()


#
#----------------------------------------------------------------
# Calculate the diffraction pattern
#
#

#
# dp parameters
#

class dp_param:
    def __init__(self,path="/home/amartin1",tag="tag",nx=128):
        self.wl = 0.5e-10
        self.pixel_width = 4.0e-5
        self.detector_z = 3.0e-3
        self.nphotons = 1e14
        self.beamArea = 1e-16
        self.rx = 0
        self.ry = 0
        self.rz = 0
        self.th = 0
        self.path = path
        self.tag = tag
        self.nx = nx
        self.nrho = 0
        self.pdbname = self.path+"pdb_polyhedra.txt"

    def write_dp_config(self,dpfname="dpconfig.txt"):

        fcorr = open(dpfname,'w')
        fcorr.write("outpath = "+self.path+'\n')
        fcorr.write("tag = "+self.tag+'\n')
        fcorr.write("pdbname = "+self.pdbname+'\n')
        fcorr.write("sfname = "+"/scratch/amartin1/Work/codes/C-code/diffractionSim/Bwxray.fac"+'\n')
        fcorr.write("wavelength =  "+str(self.wl)+'\n'  )
        fcorr.write("pixel_width =  "+str(self.pixel_width)+'\n'  )
        fcorr.write("detector_z =  "+str(self.detector_z)+'\n'  )
        fcorr.write("nx = "+str(self.nx)+'\n')
        fcorr.write("nphotons = "+str(self.nphotons)+'\n')
        fcorr.write("beamArea = "+str(self.beamArea)+'\n')
        fcorr.write("rotation_axis_x = "+str(self.rx)+'\n')
        fcorr.write("rotation_axis_y = "+str(self.ry)+'\n')
        fcorr.write("rotation_axis_z = "+str(self.rz)+'\n')
        fcorr.write("rotation_angle = "+str(self.th)+'\n')
        fcorr.close()

    def append_Bl_config_params(self,dpfname="dpconfig.txt"):

         fcorr = open(dpfname,'a')
         fcorr.write("nrho = "+str(self.nrho)+'\n')
         fcorr.close()

dpp = dp_param(path,tag,nx)
dpfname = path+"dpconfig.txt"

if dp_flag == 1:
    
     #
     # Create a configuration file
     #
     dpp.write_dp_config(dpfname)
     
     #
     # create Gaussian
     #
     gwid = dpp.nx / 4.0
     g = make_gaussian( dpp.nx, dpp.nx, gwid, gwid )


     #
     # MOVED THE CALCULATION OF DIFFRACTION PATTERNS INTO A COMMON LOOP WITH CORRELATIONS
     #


# sys.exit()

#
#----------------------------------------------------------------
# Calculate the correlation
#
#

#
# Create a configuration file
#
class corr_param:
    def __init__(self,path="/home/amartin1",tag="tag",dpoutname="dp.dbin",
                 nx=128,\
                 wl=1e-10,pw=1.0,dz=1.0,nth=200):
                 
                 self.wl = wl
                 self.pixel_width = pw
                 self.detector_z = dz
                 self.path = path
                 self.tag = tag
                 self.cx = nx/2
                 self.cy = nx/2
                 self.dpoutname = dpoutname
                 self.nth = nth

    def write_corr_config(self,corrfname="correlation_config.txt"):

        fcorr = open(corrfname,'w')
        fcorr.write("input = "+self.dpoutname+'\n' )
        fcorr.write("outpath = "+self.path+'\n')
        fcorr.write("tag = "+self.tag+'\n')
        fcorr.write("wavelength =  "+str(self.wl)+'\n'  )
        fcorr.write("pixel_width =  "+str(self.pixel_width)+'\n'  )
        fcorr.write("detector_z =  "+str(self.detector_z)+'\n'  )
        fcorr.write("cx =  "+str(self.cx)+'\n'  )
        fcorr.write("cy =  "+str(self.cy)+'\n'  )
        fcorr.write("nth =  "+str(self.nth)+'\n'  )
        fcorr.close()


if corr_flag == 1:
     corrfname = path+"corrconfig.txt"
     dpoutname = "None"
     cp = corr_param( path,tag,dpoutname,nx,dpp.wl,dpp.pixel_width,dpp.detector_z,nth)


#
# Calculate diffraction patterns and their correlations
#
for i in np.arange( npatterns*ninc ):
          

     if dp_flag == 1:
          print "\nCalculating diffraction pattern ", i+1
          dpp.tag = tag+"_"+str(i) 
          dpp.rx, dpp.ry, dpp.rz, dpp.th = random_rotation()
          dpp.pdbname = path+"pdb_polyhedra_"+str(ipat)+".txt"
          dpp.write_dp_config(dpfname)
          os.system("/scratch/amartin1/Work/codes/C-code/diffractionSim/diffractionSim "+dpfname) 
    
          if g_flag == 1:
               fname = path+tag+"_"+str(i)+"_diffraction_lowq.dbin"
               dpin = read_dbin( fname )
               dpin *= g
               write_dbin( fname, dpin )
     
          if i==0:
               dpsum = np.copy(dpin) #/ float(ninc)
          else:
               dpsum += dpin #/ float(ninc)


          if ( (i%ninc)==0 ):
               fname = path+tag+"_"+str(i%ninc)+"_diffraction_lowq.dbin"
               write_dbin( fname, dpsum )
               dpsum *= 0.0


     if corr_flag == 1:

          #
          # run padfcorr
          #

          if ( (i%ninc)==0 ):
               print "\nCorrelating pattern ", i%ninc

               cp.dpoutname = path+tag+"_"+str(i%ninc)+"_diffraction_lowq.dbin"
               cp.tag = tag #+"_"+str(i)
               cp.write_corr_config(corrfname)
               
               start = time.time()
               os.system("/scratch/amartin1/Work/codes/C-code/padfcorr/padfcorr "+corrfname) 
               print "padfcorr took :", time.time() - start, " seconds"
          
          # cname = path+tag+"_"+str(i)+"_correlation.dbin"
               cname = path+tag+"_correlation.dbin"
               corr = read_correlation( cname , 0)
               corr = corr.reshape( nx/2, nx/2, nth )
          
               if(i==0):
                    corrsum = np.copy(corr)
               else:
                    corrsum += corr

               os.system("rm "+cp.dpoutname )
          
          pdbname = path+"pdb_polyhedra_"+str(i)+".txt"
          os.system("rm "+pdbname )

          pdbname = path+tag+"_"+str(i)+"_log.txt"
          os.system("rm "+pdbname )

if corr_flag == 1:
     corrsum *= 1.0/float(npatterns)

     #
     # Output the summed correlation function
     #
     write_dbin( path+tag+"_padfcorr_correlation_sum.dbin", corrsum )

     #
     # Look at the summed correlation function
     #
     # plt.imshow( corrsum[30,:,:] )
     # plt.draw()
     # plt.show()


#
#----------------------------------------------------------------
# Calculate Bl(q,q') matrices directly
#
#
if mm_flag == 1:
     dpp.tag = "mm"
     dpp.write_dp_config(dpfname)
     dpp.append_Bl_config_params( dpfname)
     os.system("/scratch/amartin1/Work/codes/C-code/molmoment/molmoment "+dpfname)



#
#----------------------------------------------------------------
# Calculate the padf
#
#

if padf_flag == 1:
     #
     # Parameters
     #
     nl = 156
     nthp = nth
     nr = nx #/2
     nq = nx/2
     qmax = 1.38366e10
     rmax = 2.31271e-9    # 4.265e-9   #1.02415e-7

     #rmax *= 2
     #nr *= 2

     #
     # Create a configuration file
     #
     pfname = path+"padf_config.txt"
     fcorr = open(pfname,'w')
     fcorr.write("correlationfile = "+path+tag+"_padfcorr_correlation_sum.dbin"+'\n' )
     fcorr.write("outpath = "+path+'\n')
     fcorr.write("tag = "+tag+"_padf2"+'\n')
     fcorr.write("wavelength =  "+str(dpp.wl)+'\n'  )
     fcorr.write("nthq =  "+str(nthp)+'\n'  )
     fcorr.write("nq =  "+str(nq)+'\n'  )
     fcorr.write("nr =  "+str(nr)+'\n'  )
     fcorr.write("nl =  "+str(nl)+'\n'  )
     fcorr.write("qmax =  "+str(qmax)+'\n'  )
     fcorr.write("rmax =  "+str(rmax)+'\n'  )
     fcorr.close()
     
     #
     # run padf
     #
     print "\nCalculating the padf"
     start = time.clock()
     #os.system("gdb --args /scratch/amartin1/Work/codes/C-code/padf/padf "+pfname) 
     os.system("/scratch/amartin1/Work/codes/C-code/padf/padf "+pfname) 
     print "padf took :", time.clock() - start, " seconds"


