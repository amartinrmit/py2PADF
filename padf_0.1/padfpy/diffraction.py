
import os
import numpy as np
import multiprocessing as mp
# import padf_params as pp
import glob


#TODO: Add in the different ways of estimating multiple protein scattering
#      Add Poisson noise and other kinds of noise to the diffraction pattern
#      Add a background scattering term



#
#----------------------------------------------------------------
# Diffraction pattern class
#
#

class diffraction:
    def __init__(self,path="~/",tag="tag",nx=128,\
                 wl=1e-10,pw=1e-5,dz=1e-3,nph=1e12,beamArea=1e-14,\
                     rx=0,ry=0,rz=0,th=0,nrho=0,pdbname="molecule.pdb",\
                     scaledflag=1,radflag=0,gwid=10.0,\
                     npat=100,nthreads=1,gflag=0):
        self.wl = wl
        self.pixel_width = pw
        self.detector_z = dz
        self.nphotons = nph
        self.beamArea = beamArea
        self.rx = rx
        self.ry = ry
        self.rz = rz
        self.th = th
        self.path = path
        self.tag = tag
        self.nx = nx
        self.nrho = nrho
        self.pdbname = pdbname
        self.scaledflag = scaledflag
        self.radial_average_flag = radflag
        self.pdbInputType = 'pdb'     # single pdb or directory
        self.order = 'seq'    # 'seq' = sequential; 'random' = random order
        self.noiseflag = 0
        
        # parameters not required by diffractionSim
        self.gwid = gwid
        self.npatterns = npat
        self.nthreads = nthreads
        self.gflag = gflag
        self.flist = []
        self.pdblist = []
        self.pdbpath = path
        self.pdbtag = tag


#   def __init__(self,p=pp.params()):
#
#              self.p = p
#


    def write_dp_config(self,dpfname="dpconfig.txt"):

        fcorr = open(dpfname,'w')
        fcorr.write("outpath = "+self.path+'\n')
        fcorr.write("tag = "+self.tag+'\n')
        fcorr.write("pdbname = "+ self.pdbname+'\n' )
        fcorr.write("sfname = "+"/Users/e38496/Work/Research/codes/code_projects/diffractionsim/diffractionsim_0.1/Bwxray.fac"+'\n')
        fcorr.write("henkepath = "+"/Users/e38496/Work/Research/codes/code_projects/diffractionsim/diffractionsim_0.1/sf/"+'\n')
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
        fcorr.write("scaled_flag = "+str(self.scaledflag)+'\n')
        fcorr.write("radial_average_flag = "+str(self.radial_average_flag)+'\n')
        fcorr.write("#gwid = "+str(self.gwid)+'\n')
        fcorr.write("noise_flag = "+str(self.noiseflag)+'\n')
        if self.latticeflag == 1:
            fcorr.write("lattice_flag = yes\n" )
        else:
            fcorr.write("lattice_flag = no\n" )
        fcorr.write("spotwidth = "+str(self.spotwidth)+"\n")
        fcorr.write("na = "+str(self.na)+"\n")
        fcorr.write("nb = "+str(self.nb)+"\n")
        fcorr.write("nc = "+str(self.nc)+"\n")
        fcorr.write("astar_x = "+str(self.astarx)+"\n")
        fcorr.write("astar_y = "+str(self.astary)+"\n")
        fcorr.write("astar_z = "+str(self.astarz)+"\n")
        fcorr.write("bstar_x = "+str(self.bstarx)+"\n")
        fcorr.write("bstar_y = "+str(self.bstary)+"\n")
        fcorr.write("bstar_z = "+str(self.bstarz)+"\n")
        fcorr.write("cstar_x = "+str(self.cstarx)+"\n")
        fcorr.write("cstar_y = "+str(self.cstary)+"\n")
        fcorr.write("cstar_z = "+str(self.cstarz)+"\n")

        fcorr.close()

    def append_Bl_config_params(self,dpfname="dpconfig.txt"):

         fcorr = open(dpfname,'a')
         fcorr.write("nrho = "+str(self.nrho)+'\n')
         fcorr.close()

    def random_rotation( self ):

        u1, u2, u3 = np.random.rand(3)
        h = [     np.sqrt(1-u1)*np.sin(2*np.pi*u2), np.sqrt(1-u1)*np.cos(2*np.pi*u2), \
                      np.sqrt(u1)*np.sin(2*np.pi*u3),   np.sqrt(u1)*np.cos(2*np.pi*u3) ]
        th = 2.0*np.arccos( h[0] )
        rx = h[1] / np.sin( th / 2.0 )
        ry = h[2] / np.sin( th / 2.0 )
        rz = h[3] / np.sin( th / 2.0 )
        
        return rx, ry, rz, th
    
    def dp_calc( self, qfname ):
        dirpath = os.path.dirname(os.path.realpath(__file__))       
        os.system(dirpath+"/../diffractionSim/diffractionSim "+qfname ) 


    def generate_diffraction_patterns(self, npatterns, nstart=0):

        dpfname = self.path+"dpconfig.txt"
        self.write_dp_config(dpfname)

        #
        # Create a configuration file
        #
        self.write_dp_config(dpfname)
     
         # store the tag
        tag = self.tag

        #
        # Calculate dp
        #
        for i in np.arange( npatterns/self.nthreads ):

          processes = []
          d = []
          for j in np.arange( self.nthreads ):
               m = i*self.nthreads+j+1
               print "\nCalculating diffraction pattern ", m

               # set the input pdb to be the next file name in the list
               print "pdbname :", self.pdbname
               if self.pdbInputType == 'dir':
                   if self.order == 'seq':
                       print "debug diffraction.py len(pdblist)", len(self.pdblist)
                       ir = (m-1)%(len(self.pdblist)-nstart) + nstart
                       self.pdbname = self.pdblist[ir]
                   elif self.order == 'random':
                       ir = int( np.random.rand()*len(self.pdblist))
                       self.pdbname = self.pdblist[ir]
                   else:
                       print "<diffraction.py> unknown order for reading pdb files."
                       exit()


               self.tag = tag+"_"+str(m) 
               self.rx, self.ry, self.rz, self.th = self.random_rotation()
               dpfname = "dpconfig"+str(j)+".txt"
               self.write_dp_config(dpfname)
               d.append(dpfname)


          for j in np.arange( self.nthreads ):
#               print "thread j :", j, d[j]
               p = mp.Process(target=self.dp_calc, args=(d[j],))
               p.start()
               processes.append(p)

          for p in processes:
               p.join()

 
        # reset stored values
        self.tag = tag

    # remove all diffraction patterns
    def remove_diffraction_patterns( self, npatterns ):

        #
        # Calculate dp
        #
        for i in np.arange( npatterns/self.nthreads ):

          for j in np.arange( self.nthreads ):
              m = i*self.nthreads+j+1
              fname = self.path+self.tag+"_"+str(m)+"_diffraction_*.dbin"
              os.system("rm -f "+fname)



    # shift - a 2D version of numpy's roll
    def array_shift(array,xshift=0,yshift=0):
	array = np.roll(array,xshift,0)
	array = np.roll(array,yshift,1)
	return array

    # make a 2D array with a gaussian
    def make_gaussian(nx, ny, rad=None, rady=-1., cenx=None, ceny=None, invert=0, norm=False, power=2 ): 
        # set defaults
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
        # print x.size, x.shape
        y = np.outer(np.ones(nx),np.arange(0-ny/2,ny-ny/2,1))
        # print y.size, y.shape

        a = np.zeros([nx,ny])
        # a = np.exp(-(x**2)/radsq  - ( y**2)/radysq)
        a = np.exp(-(x**power)/radsq  - ( y**power)/radysq)
        a[ nx/2, ny/2 ] = 1.0

        a = self.array_shift(a,cenx-nx/2,ceny-ny/2)

        # normalise if required
        if norm == True: a *= 1./np.sum(a)
    
        return a
#end make_gaussian
    
    def set_dp_inputType( self ):

        fext = os.path.splitext( str(self.pdbname) )[-1]        

        if len(fext) == 0:
            self.pdbInputType = 'dir'
            self.pdblist = self.get_file_list( self.pdbpath, self.pdbtag, ".pdb" )

        elif fext > 1:
            if (fext=='.pdb') or (fext=='.txt'):
                self.pdbInputType = 'pdb'
        else:
            print "error <diffraction.py> unknown structure file type. Expected .pdb"



    def get_file_list( self, path, tag, fext=".*" ):

        try:
            if not os.path.exists(path):
                print "sample path does not exist :", path
            flist = sorted(glob.glob( path + tag+"*"+fext ), key=os.path.getmtime )
            if len(flist) != 0:
                return flist
            else:
                print "<padfpy.py> no diffraction files found. Check the samplepath and sampletag"
                print "samplepath: ", path
                print "sampletag: ", tag
                sys.exit()
        except:
            print "Error (padfpy.py). Files could not be read from the samplepath. Check samplepath and sampletag"
            print path + tag+"*"+fext
