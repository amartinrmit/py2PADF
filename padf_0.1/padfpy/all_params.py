
import sys, getopt
import argparse


class params:

    def __init__(self):

        #
        # 
        #
        self.__version__ = "padf_0.1"
        self.init_flags()
        self.init_common_params()
        self.init_diffraction_params()
        self.init_correlation_params()
        self.init_padf_params()
#        self.init_plotting_params()

    def init_flags(self):
        
        # ADD COMMENTS
        self.dp_flag = 0
        self.corr_flag = 0
        self.padf_flag = 0
        self.units_flag = 0
        self.g_flag = 0
        self.sym_flag = 0
        self.mask_flag = 0
        self.hpfilter_flag = 0
        self.diffraction_average_flag = 0
        self.diffraction_av_correlation_flag = 0
        self.interpolation_noise_filter_flag = 0
        self.corr_shift_flag = 0
        self.corr_bg_calc_flag = 0
        self.corr_bg_sub_flag = 0
        self.samplepath_flag = 0
        self.crop_flag = 0
        self.dp_shift_flag = 0
        self.corr_diff_flag = 0


        # flags to model the noise and the background

    def init_common_params(self):

        self.path = "~"
        self.tag = "tag"
        self.inputpath = "None"
        self.inputtag = "None"
        self.inputtail = "None"
        self.h5field = "None"
        self.wavelength = 0.0
        self.detector_z = 0.0
        self.beamArea = 0.0
        self.pixelwidth = 0.0
        self.gwid = 0.0
        self.npatterns = 1
        self.nthreads = 1
        self.samplepath = ""   #if not set this should be assumed to be equal to path etc.
        self.sampletag = ""
        self.hp_filter = 0
        self.lp_filter = 0
        self.maskname = "mask"
        self.nxcrop = 0
        self.nycrop = 0
        self.shiftx = 0
        self.shifty = 0
        self.rebin = 0
        self.nstart = 0


        # simulation params
        self.nphotons = 0

        # simulation / correlation
        self.nx = -1
        self.ny = -1
        self.pwx = -1
        self.pwy = -1

    
    def init_diffraction_params(self):

        self.dp_path = "Default"
        self.dp_tag  = "Default"
        self.rx = 0.0
        self.ry = 0.0
        self.rz = 0.0
        self.th = 0.0
        self.scaledflag = 0
        self.radial_average_flag = 0
        self.pdbname = None
        self.order = 'seq'
        self.delete_sim_dp = 0
        self.noiseflag = 0
        self.latticeflag = 0
        self.spotwidth = 0.02
        self.na = 1
        self.nb = 1
        self.nc = 1
        self.astarx, self.astary, self.astarz = 1.0, 0.0, 0.0
        self.bstarx, self.bstary, self.bstarz = 0.0, 1.0, 0.0
        self.cstarx, self.cstary, self.cstarz = 0.0, 0.0, 1.0
        # continue with other lattice parameters and below...

    def init_correlation_params(self):
        
        self.cx = 0
        self.cy = 0
        self.nth = 0
        self.nq = 0
        self.sym_filter = 0

    def init_padf_params(self):

        self.nl = 0
        self.nlmin = 0
        self.padf_nth = 0
        self.padf_nr = 0
        self.padf_nq = 0
        self.padf_qmax = 0.0
        self.padf_rmax = 0.0
        # maybe add filter flags??

        # are there any special variables for unit correction?

        # plotting parameters
        self.section1 = 0
        self.section2 = 0
        self.section3 = 0
        self.theta = 0
        self.r = 0.0
        self.r2 = 0.0


#    def init_plotting_params(self):



    def write_all_params(self, name="None"):
        
        if name=="None":
            f = open( self.path+self.tag+"padfpy_log_all_params.log", 'w')


        f.write( "#All parameters Pair Angle Distribution Function (padfpy.py). Version "+self.__version__+"\n")
        f.write( "path = "+self.path+"\n" )
        f.write( "tag = "+self.tag+"\n" )
        f.write( "inputpath = "+self.inputpath+"\n" )
        f.write( "inputtag = "+self.inputtag+"\n" )
        f.write( "inputtail = "+self.inputtail+"\n" )
        if not self.h5field == "None":
            f.write( "h5field = "+self.h5field+"\n" )

        f.write( "#**** FLAGS ****\n")
        f.write( "dp_flag = "+str(self.dp_flag)+"\n" )
        f.write( "g_flag = "+str(self.g_flag)+"\n" )
        f.write( "corr_flag = "+str(self.corr_flag)+"\n" )
        f.write( "units_flag = "+str(self.units_flag)+"\n" )
        f.write( "symmetry_flag = "+str(self.symmetry_flag)+"\n" )
        f.write( "padf_flag = "+str(self.padf_flag)+"\n" )
        f.write( "mask_flag = "+str(self.mask_flag)+"\n" )
        f.write( "hpfilter_flag = "+str(self.hpfilter_flag)+"\n" )
        f.write( "diffraction_average_flag = "+str(self.diffraction_average_flag)+"\n")
        f.write( "diffraction_av_correlation_flag = "+str(self.diffraction_av_correlation_flag)+"\n")
        f.write( "interpolation_noise_filter_flag = "+str(self.interpolation_noise_filter_flag)+"\n")
        f.write( "corr_shift_flag = "+str(self.corr_shift_flag)+"\n")
        f.write( "corr_bg_calc_flag = "+str(self.corr_bg_calc_flag)+"\n")
        f.write( "corr_bg_sub_flag = "+str(self.corr_bg_sub_flag)+"\n")
        f.write( "corr_diff_flag = "+str(self.corr_diff_flag)+"\n")
        f.write( "crop_flag = "+str(self.crop_flag)+"\n" )
        f.write( "dp_shift_flag = "+str(self.dp_shift_flag)+"\n" )
        f.write( "\n" )

        f.write( "#**** Common variables ****\n")
        f.write( "samplepath = "+self.samplepath+"\n" )
        f.write( "sampletag = "+self.sampletag+"\n" )
        f.write( "maskname = "+self.maskname+"\n" )
        f.write( "wavelength = "+str(self.wavelength)+"\n")
        f.write( "detector_z = "+str(self.detector_z)+"\n")
        f.write( "beamArea = "+str(self.beamArea)+"\n")
        f.write( "pixelwidth = "+str(self.pixelwidth)+"\n")
        f.write( "gwid = "+str(self.gwid)+"\n")
        f.write( "npatterns = "+str(self.npatterns)+"\n" )
        f.write( "nthreads = "+str(self.nthreads)+"\n" )
        f.write( "hp_filter = "+str(self.hp_filter)+"\n" )
        f.write( "lp_filter = "+str(self.lp_filter)+"\n" )
        f.write( "nstart = "+str(self.nstart)+"\n" )
        f.write( "nx = "+str(self.nx)+"\n" )
        f.write( "ny = "+str(self.ny)+"\n" )
        f.write( "pwx = "+str(self.pwx)+"\n" )
        f.write( "pwy = "+str(self.pwy)+"\n" )
        f.write( "\n")



        f.write("#**** Simulation Parameters ****\n")
        f.write( "nphotons = "+str(self.nphotons)+"\n" )
        f.write( "\n")

        f.write("#**** Diffraction Parameters ****\n")
        f.write("dp_path = "+self.dp_path+"\n")
        f.write("dp_tag = "+self.dp_tag+"\n")
        f.write( "rx = "+str(self.rx)+"\n" )
        f.write( "ry = "+str(self.ry)+"\n" )
        f.write( "rz = "+str(self.rz)+"\n" )
        f.write( "th = "+str(self.th)+"\n" )
        f.write( "order = "+self.order+"\n" )
        f.write( "scaledflag = "+str(self.scaledflag)+"\n" )
        f.write( "radial_average_flag = "+str(self.radial_average_flag)+"\n" )
        f.write( "pdbname = "+str(self.pdbname)+"\n")
        f.write( "delete_sim_dp = "+str(self.delete_sim_dp)+"\n")
        f.write( "noiseflag = "+str(self.noiseflag)+"\n" )
        f.write( "latticeflag = "+str(self.latticeflag)+"\n" )
        f.write( "spotwidth = "+str(self.spotwidth)+"\n" )
        f.write( "na = "+str(self.na)+"\n" )
        f.write( "nb = "+str(self.nb)+"\n" )
        f.write( "nc = "+str(self.nc)+"\n" )
        f.write( "astarx = "+str(self.astarx)+"\n")
        f.write( "astary = "+str(self.astary)+"\n")
        f.write( "astarz = "+str(self.astarz)+"\n")
        f.write( "bstarx = "+str(self.bstarx)+"\n")
        f.write( "bstary = "+str(self.bstary)+"\n")
        f.write( "bstarz = "+str(self.bstarz)+"\n")
        f.write( "cstarx = "+str(self.cstarx)+"\n")
        f.write( "cstary = "+str(self.cstary)+"\n")
        f.write( "cstarz = "+str(self.cstarz)+"\n")


        f.write( "\n")

        f.write("#**** Correlation Parameters ****\n")
        f.write( "cx = "+str(self.cx)+"\n" )
        f.write( "cy = "+str(self.cy)+"\n" )
        f.write( "nth = "+str(self.nth)+"\n" )
        f.write( "nq = "+str(self.nq)+"\n" )
        f.write( "sym_filter = "+str(self.sym_filter)+"\n" )
        f.write( "nxcrop = "+str(self.nxcrop)+"\n" )
        f.write( "nycrop = "+str(self.nycrop)+"\n" )
        f.write( "shiftx = "+str(self.shiftx)+"\n" )
        f.write( "shifty = "+str(self.shifty)+"\n" )
        f.write( "rebin = "+str(self.rebin)+"\n" )
        f.write( "\n")

        f.write("#**** Padf Parameters ****\n")
        f.write( "nl = "+str(self.nl)+"\n" )
        f.write( "padf_nth = "+str(self.padf_nth)+"\n" )
        f.write( "padf_nr = "+str(self.padf_nr)+"\n" )
        f.write( "padf_nq = "+str(self.padf_nq)+"\n" )
        f.write( "padf_qmax = "+str(self.padf_qmax)+"\n" )
        f.write( "padf_rmax = "+str(self.padf_rmax)+"\n" )
        f.write( "section1 = "+str(self.section1)+"\n" )
        f.write( "section2 = "+str(self.section2)+"\n" )
        f.write( "section3 = "+str(self.section3)+"\n" )
        f.write( "padfplot_theta = "+str(self.theta)+"\n" )
        f.write( "padfplot_r = "+str(self.r)+"\n" )
        f.write( "padfplot_r2 = "+str(self.r2)+"\n" )

        f.close()
    

    def read_config_file( self, fname ):
        
        #print "check fname:", fname
        f = open( fname, 'r' )
        #print fname, f

        for line in f:
            
            if line[0] == '#': continue
            if line[0] == '\n': continue

            bits = line.split()
            if bits[0] == "path":
                self.path = bits[2]
            elif bits[0] == "tag":
                self.tag = bits[2]
            elif bits[0] == "inputpath":
                self.inputpath = bits[2]
            elif bits[0] == "inputtag":
                self.inputtag = bits[2]
            elif bits[0] == "input_tail":
                self.inputtail = bits[2]
            elif bits[0] == "h5field":
                self.h5field = bits[2]
            elif bits[0] == "dp_flag":                            # flags
                self.dp_flag = int(bits[2])
            elif bits[0] == "corr_diff_flag":                            # flags
                self.corr_diff_flag = int(bits[2])
            elif bits[0] == "g_flag":
                self.g_flag = int(bits[2])
            elif bits[0] == "corr_flag":
                self.corr_flag = int(bits[2])
            elif bits[0] == "units_flag":
                self.units_flag = int(bits[2])
            elif bits[0] == "symmetry_flag":
                self.symmetry_flag = int(bits[2])
            elif bits[0] == "padf_flag":
                self.padf_flag = int(bits[2])
            elif bits[0] == "mask_flag":
                self.mask_flag = int(bits[2])
            elif bits[0] == "hpfilter_flag":
                self.hpfilter_flag = int(bits[2])
            elif bits[0] == "diffraction_average_flag":
                self.diffraction_average_flag = int(bits[2])
            elif bits[0] == "interpolation_noise_filter_flag":
                self.interpolation_noise_filter_flag = int(bits[2])
            elif bits[0] == "corr_shift_flag":
                self.corr_shift_flag = int(bits[2])
            elif bits[0] == "corr_bg_sub_flag":
                self.corr_bg_sub_flag = int(bits[2])
            elif bits[0] == "corr_bg_calc_flag":
                self.corr_bg_calc_flag = int(bits[2])
            elif bits[0] == "crop_flag":
                self.crop_flag = int(bits[2])
            elif bits[0] == "dp_shift_flag":
                self.dp_shift_flag = int(bits[2])
            elif bits[0] == "samplepath":                          # Common variables
                self.samplepath = bits[2]
            elif bits[0] == "sampletag":
                self.sampletag = bits[2]
            elif bits[0] == "maskname":
                self.maskname = bits[2]
            elif bits[0] == "wavelength":
                self.wavelength = float(bits[2])
            elif bits[0] == "detector_z":
                self.detector_z = float(bits[2])
            elif bits[0] == "beamArea":
                self.beamArea = float(bits[2])
            elif bits[0] == "pixelwidth":
                self.pixelwidth = float(bits[2])
            elif bits[0] == "gwid":
                self.gwid = float(bits[2])
            elif bits[0] == "npatterns":
                self.npatterns = int(bits[2])
            elif bits[0] == "nthreads":
                self.nthreads = int(bits[2])
            elif bits[0] == "hp_filter":
                self.hp_filter = int(bits[2])
            elif bits[0] == "lp_filter":
                self.lp_filter = int(bits[2])
            elif bits[0] == "nstart":
                self.nstart = int(bits[2])
            elif bits[0] == "nphotons":                             # Simulation parameters
                self.nphotons = int(bits[2])
            elif bits[0] == "nx":
                self.nx = int(bits[2])                            
            elif bits[0] == "ny":
                self.ny = int(bits[2])                            
            elif bits[0] == "pwx":
                self.pwx = float(bits[2])                            
            elif bits[0] == "pwy":
                self.pwy = float(bits[2])                         
            elif bits[0] == "dp_path":                                 # diffraction parameters
                self.dp_path = bits[2]
            elif bits[0] == "dp_tag":
                self.dp_tag = bits[2]
            elif bits[0] == "delete_sim_dp":                                 # diffraction parameters
                self.delete_sim_dp = int(bits[2])
            elif bits[0] == "rx":
                self.rx = float(bits[2])
            elif bits[0] == "ry":
                self.ry = float(bits[2])
            elif bits[0] == "rz":
                self.rz = float(bits[2])
            elif bits[0] == "th":
                self.th = float(bits[2])
            elif bits[0] == "order":
                self.order = bits[2]
            elif bits[0] == "na":
                self.na = int(bits[2])
            elif bits[0] == "nb":
                self.nb = int(bits[2])
            elif bits[0] == "nc":
                self.nc = int(bits[2])
            elif bits[0] == "astarx":
                self.astarx = float(bits[2])
            elif bits[0] == "astary":
                self.astary = float(bits[2])
            elif bits[0] == "astarz":
                self.astarz = float(bits[2])
            elif bits[0] == "bstarx":
                self.bstarx = float(bits[2])
            elif bits[0] == "bstary":
                self.bstary = float(bits[2])
            elif bits[0] == "bstarz":
                self.bstarz = float(bits[2])
            elif bits[0] == "cstarx":
                self.cstarx = float(bits[2])
            elif bits[0] == "cstary":
                self.cstary = float(bits[2])
            elif bits[0] == "cstarz":
                self.cstarz = float(bits[2])
            elif bits[0] == "spotwidth":
                self.spotwidth = float(bits[2])
            elif bits[0] == "latticeflag":
                self.latticeflag = int(bits[2])
            elif bits[0] == "noiseflag":
                self.noiseflag = int(bits[2])
            elif bits[0] == "scaledflag":
                self.scaledflag = int(bits[2])
            elif bits[0] == "radial_average_flag":
                self.radial_average_flag = int(bits[2])
            elif bits[0] == "pdbname":
                self.pdbname = bits[2]
            elif bits[0] == "cx":                                 # correlation flags
                self.cx = float(bits[2])
            elif bits[0] == "cy":
                self.cy = float(bits[2])
            elif bits[0] == "nth":
                self.nth = int(bits[2])
            elif bits[0] == "nq":
                self.nq = int(bits[2])
            elif bits[0] == "sym_filter":
                self.sym_filter = int(bits[2])
            elif bits[0] == "nxcrop":
                self.nxcrop = int(bits[2])
            elif bits[0] == "nycrop":
                self.nycrop = int(bits[2])
            elif bits[0] == "shiftx":
                self.shiftx = float(bits[2])
            elif bits[0] == "shifty":
                self.shifty = float(bits[2])
            elif bits[0] == "rebin":
                self.rebin = int(bits[2])
            elif bits[0] == "nl":                                 # padf flags
                self.nl = int(bits[2])
            elif bits[0] == "padf_nth":
                self.padf_nth = int(bits[2])
            elif bits[0] == "padf_nq":
                self.padf_nq = int(bits[2])
            elif bits[0] == "padf_nr":
                self.padf_nr = int(bits[2])
            elif bits[0] == "padf_qmax":
                self.padf_qmax = float(bits[2])
            elif bits[0] == "padf_rmax":
                self.padf_rmax = float(bits[2])
            elif bits[0] == "section1":
                self.section1 = int(bits[2])
            elif bits[0] == "section2":
                self.section2 = int(bits[2])
            elif bits[0] == "section3":
                self.section3 = int(bits[2])
            elif bits[0] == "padfplot_theta":
                self.theta = float(bits[2])
            elif bits[0] == "padfplot_r":
                self.r = float(bits[2])
            elif bits[0] == "padfplot_r2":
                self.r2 = float(bits[2])


        f.close()



## command line arguments could be a separate structure of arguments. 
#  then a second routine would overwrite the main arguments with teh command line arguments..


    def get_commandline_args( self ):

        
        # Not implemented:
        #  interpolation_noise_filter_flag
        #  corr_shift_flag

        parser = argparse.ArgumentParser( description="padfpy.py. Python library for calculating the pair-angle distribution function." )

        parser.add_argument( "--config", "-c", nargs=1, help="configuration file name" )
        parser.add_argument( "--path", "-p", nargs=1, help="Path to working directory" )
        parser.add_argument( "--tag", "-t", nargs=1, help="A short string (tag) to name all output files." )
        parser.add_argument( "--inputpath", nargs=1, help="Path to directory with input files" )
        parser.add_argument( "--inputtag",  nargs=1, help="A short string (tag) to name all input files." )
        parser.add_argument( "--inputtail",  nargs=1, help="the end of the correlation file name to input." )
        parser.add_argument( "--h5field", nargs=1, help="HDF5 field that contains the image (for .h5 input)." )
        parser.add_argument( "--samplepath", nargs=1, help="Path to working directory with diffraction pattern files" )
        parser.add_argument( "--sampletag", nargs=1, help="A short string (tag) used to name diffraction pattern files.." )
        parser.add_argument( "--maskname", nargs=1, help="A filename (including path) of the mask image (= 0 for excluded pixels)" )
        parser.add_argument( "--dp_flag", nargs=1, help="1 - calculate diffraction patterns; 0 - patterns not calculated", type=int)
        parser.add_argument( "--g_flag", nargs=1, help="1 - use Gaussian filter; 0 - No Gaussian filter", type=int)
        parser.add_argument( "--mask_flag", nargs=1, help="1 - mask applied; 0 - mask not applied", type=int)
        parser.add_argument( "--hpfilter_flag", nargs=1, help="1 - high pass filtered applied; 0 - high pass filter not applied", type=int)
        parser.add_argument( "--corr_flag", nargs=1, help="1 - calculate correlation function; 0 - correlation function not calculated", type=int)
        parser.add_argument( "--padf_flag", nargs=1, help="1 - padf calculated; 0 - padf not calculated", type=int)
        parser.add_argument( "--units_flag", nargs=1, help="1 - padf not corrected for units; 0 - padf corrected for units",type=int)
        parser.add_argument( "--symmetry_flag", nargs=1, help="1 - symmetry filter applied to correlation function; 0 - symmetry filter not applied",type=int)
        parser.add_argument( "--diffraction_average_flag", nargs=1, help="Background estimated from diffraction average? 1 - yes; 0 - no",type=int)
        parser.add_argument( "--corr_bg_calc_flag", nargs=1, help="Background estimated from random cross-correlations? 1 - yes; 0 - no",type=int)
        parser.add_argument( "--corr_bg_sub_flag", nargs=1, help="Background subtracted? 1 - yes; 0 - no",type=int)
        parser.add_argument( "--corr_diff_flag", nargs=1, help="Difference correlation? 1 - yes; 0 - no",type=int)
        parser.add_argument( "--scaledflag", nargs=1, help="Output scaled diffraction pattern (simulation)? 1 - yes; 0 - no",type=int)
        parser.add_argument( "--crop_flag", nargs=1, help="Crop the diffraction pattern?? 1 - yes; 0 - no",type=int)
        parser.add_argument( "--dp_shift_flag", nargs=1, help="Shift the centre of the diffraction pattern? 1 - yes; 0 - no",type=int)
        parser.add_argument( "--radial_average_flag", nargs=1, help="Output radial average of diffraction pattern (simulation)? 1 - yes; 0 - no",type=int)
        parser.add_argument( "--pdbname", nargs=1, help="Name of pdbfile for diffraction calculations")
        parser.add_argument( "--wavelength", "-w", nargs=1, help="Wavelength of radiation in metres",type=float)
        parser.add_argument( "--detector_z", "-dz", nargs=1, help="Sample-detector distance in metres",type=float)
        parser.add_argument( "--beamArea", nargs=1, help="Beam area in metres squared",type=float)
        parser.add_argument( "--pixelwidth", nargs=1, help="Pixel width in metres",type=float)
        parser.add_argument( "--gwid", nargs=1, help="Width of Gaussian filter in pixels",type=float)
        parser.add_argument( "--npatterns", "-n", nargs=1, help="Number of diffraction patterns to process",type=int)
        parser.add_argument( "--nthreads", nargs=1, help="Number of threads to use",type=int)
        parser.add_argument( "--hp_filter", nargs=1, help="width of the high pass filter (default 0 - no filter)",type=float)
        parser.add_argument( "--lp_filter", nargs=1, help="width of the low pass (Gaussian) filter (default 0 - no filter)",type=float)
        parser.add_argument( "--nstart", nargs=1, help="Starting diffraction pattern number",type=float)
        parser.add_argument( "--nphotons", nargs=1, help="Number of incident photons in simulation",type=int)
        parser.add_argument( "--nx", nargs=1, help="Size of array x dimension in data/simulation",type=int)
        parser.add_argument( "--ny", nargs=1, help="Size of array y dimension in data/simulation",type=int)
        parser.add_argument( "--pwx", nargs=1, help="X length of pixel (if different from pixel_width)",type=int)
        parser.add_argument( "--pwy", nargs=1, help="Y length of pixel (if different from pixel_width)",type=int)
        parser.add_argument( "--dp_path", nargs=1, help="Path to working directory with diffraction pattern files" )
        parser.add_argument( "--dp_tag", nargs=1, help="A short string (tag) used to name diffraction pattern files.." )
        parser.add_argument( "--rx", nargs=1, help="x coordinate of rotation axis in simulation",type=float)
        parser.add_argument( "--ry", nargs=1, help="y coordinate of rotation axis in simulation",type=float)
        parser.add_argument( "--rz", nargs=1, help="z coordinate of rotation axis in simulation",type=float)
        parser.add_argument( "--th", nargs=1, help="angle of rotation axis in simulation",type=float)
        parser.add_argument( "--na", nargs=1, help="number of reciprocal cells along astar direction",type=int)
        parser.add_argument( "--nb", nargs=1, help="number of reciprocal cells along bstar direction",type=int)
        parser.add_argument( "--nc", nargs=1, help="number of reciprocal cells along cstar direction",type=int)
        parser.add_argument( "--astarx", nargs=1, help="x coordiante of astar vector",type=float)
        parser.add_argument( "--astary", nargs=1, help="y coordiante of astar vector",type=float)
        parser.add_argument( "--astarz", nargs=1, help="z coordiante of astar vector",type=float)
        parser.add_argument( "--bstarx", nargs=1, help="x coordiante of bstar vector",type=float)
        parser.add_argument( "--bstary", nargs=1, help="y coordiante of bstar vector",type=float)
        parser.add_argument( "--bstarz", nargs=1, help="z coordiante of bstar vector",type=float)
        parser.add_argument( "--cstarx", nargs=1, help="x coordiante of cstar vector",type=float)
        parser.add_argument( "--cstary", nargs=1, help="y coordiante of cstar vector",type=float)
        parser.add_argument( "--cstarz", nargs=1, help="z coordiante of cstar vector",type=float)
        parser.add_argument( "--spotwidth", nargs=1, help="width of Bragg peaks in lattice calculation",type=float)
        parser.add_argument( "--order", nargs=1, help="order to read pdb files 'seq' = sequence; 'random' = random")
        parser.add_argument( "--noiseflag", nargs=1, help="0-no noise; 1-Poisson noise", type=int)
        parser.add_argument( "--latticeflag", nargs=1, help="0-pdb atoms only; 1-simulate crystal diffraction", type=int)
        parser.add_argument( "--delete_sim_dp", nargs=1, help="delete simulated diffraction patterns 0-no, 1-yes ", type=int)
        parser.add_argument( "--cx", nargs=1, help="x-coordinate - centre of diffraction pattern (for correlation calc)",type=float)
        parser.add_argument( "--cy", nargs=1, help="y-coordinate - centre of diffraction pattern (for correlation calc)",type=float)
        parser.add_argument( "--nth", nargs=1, help="number of angular samples in correlation function",type=int)
        parser.add_argument( "--nq", nargs=1, help="number of q-samples in correlation function",type=float)
        parser.add_argument( "--sym_filter", nargs=1, help="width of the symmetry filter (default 0 - no filter)",type=float)
        parser.add_argument( "--nxcrop", nargs=1, help="number of pixels in x-direction in cropped array",type=int)
        parser.add_argument( "--nycrop", nargs=1, help="number of pixels in y-direction in cropped array",type=int)
        parser.add_argument( "--shiftx", nargs=1, help="number of pixel shift in x direction",type=float)
        parser.add_argument( "--shifty", nargs=1, help="Number of pixel shift in y direction",type=float)
        parser.add_argument( "--rebin", nargs=1, help="dividing factor to rebin array. A power of two.",type=int)
        parser.add_argument( "--nl", nargs=1, help="maximum order of spherical harmonic co-efficients",type=int)
        parser.add_argument( "--padf_nth", nargs=1, help="number of angular samples in padf calculation (default=nth)",type=int)
        parser.add_argument( "--padf_nq", nargs=1, help="number of q-samples in padf function (default = nq)",type=int)
        parser.add_argument( "--padf_nr", nargs=1, help="number of real-space radial samples in padf function (default = nq)",type=int)
        parser.add_argument( "--padf_qmax", nargs=1, help="maximum qvector (default=derived from other geometric params)",type=float)
        parser.add_argument( "--padf_rmax", nargs=1, help="maximum r-space distance",type=float)
        parser.add_argument( "--section1", nargs=1, help="real-space cross-section to plot",type=float)
        parser.add_argument( "--section2", nargs=1, help="real-space cross-section to plot",type=float)
        parser.add_argument( "--section3", nargs=1, help="real-space cross-section to plot",type=float)
        parser.add_argument( "--padfplot_theta", nargs=1, help="theta value of cross-section",type=float)
        parser.add_argument( "--padfplot_r", nargs=1, help="real-space distance for cross-section",type=float)
        parser.add_argument( "--padfplot_r2", nargs=1, help="2nd real-space distance for cross-section",type=float)
        
        self.args = parser.parse_args()



    def update_with_commandline_args( self ):

        args = self.args


        if args.path != None:
            self.path = args.path

        if args.tag != None:
            self.tag = args.tag

        if args.h5field != None:
            self.h5field = args.h5field

        if args.samplepath != None:
            self.samplepath = args.samplepath

        if args.sampletag != None:
            self.sampletag = args.sampletag

        if args.inputpath != None:
            self.inputpath = args.inputpath
        
        if args.inputtag != None:
            self.inputtag = args.inputtag

        if args.inputtail != None:
            self.inputtail = args.inputtail

        if args.maskname != None:
            self.sampletag = args.maskname

        if args.dp_flag != None:
            self.dp_flag = args.dp_flag
            
        if args.g_flag != None:
            self.g_flag = args.gflag

        if args.corr_flag != None:
            self.corr_flag = args.corr_flag

        if args.units_flag != None:
            self.units_flag = args.units_flag

        if args.symmetry_flag != None:
            self.symmetry_flag = args.symmetry_flag

        if args.padf_flag != None:
            self.padf_flag = args.padf_flag

        if args.mask_flag != None:
            self.mask_flag = args.mask_flag

        if args.hpfilter_flag != None:
            self.hpfilter_flag = args.hpfilter_flag

        if args.diffraction_average_flag != None:
            self.diffraction_average_flag = args.diffraction_average_flag

        if args.corr_bg_calc_flag != None:
            self.corr_bg_calc_flag = args.corr_bg_calc_flag

        if args.corr_bg_sub_flag != None:
            self.corr_bg_sub_flag = args.corr_bg_sub_flag

        if args.corr_diff_flag != None:
            self.corr_diff_flag = args.corr_diff_flag

        if args.scaledflag != None:
            self.scaledflag = scaledflag

        if args.radial_average_flag != None:
            self.radial_average_flag = args.radial_average_flag

        if args.crop_flag != None:
            self.crop_flag = args.crop_flag

        if args.dp_shift_flag != None:
            self.dp_shift_flag = args.dp_shift_flag

        if args.pdbname != None:
            self.pdbname = args.pdbname

        if args.wavelength != None:
            self.wavelength = args.wavelength

        if args.detector_z != None:
            self.detector_z = args.detector_z

        if args.beamArea != None:
            self.beamArea = args.beamArea

        if args.pixelwidth != None:
            self.pixelwidth = args.pixelwidth

        if args.gwid != None:
            self.gwid = args.gwid

        if args.npatterns != None:
            self.gwid = args.npatterns

        if args.nthreads != None:
            self.nthreads = args.nthreads

        if args.hp_filter != None:
            self.hp_filter = args.hp_filter

        if args.lp_filter != None:
            self.lp_filter = args.lp_filter

        if args.nstart != None:
            self.nstart = args.nstart

        if args.nphotons != None:
            self.nphotons = args.nphotons

        if args.nx != None:
            self.nx = args.nx

        if args.ny != None:
            self.ny = args.ny

        if args.pwx != None:
            self.pwx = args.pwx

        if args.pwy != None:
            self.pwy = args.pwy

        if args.dp_path != None:
            self.dp_path = args.dp_path

        if args.dp_tag != None:
            self.dp_tag = args.dp_tag

        if args.rx != None:
            self.rx = args.rx

        if args.ry != None:
            self.ry = args.ry

        if args.rz != None:
            self.rz = args.rz

        if args.th != None:
            self.th = args.th

        if args.na != None:
            self.na = args.na

        if args.nc != None:
            self.nc = args.nc

        if args.nb != None:
            self.nb = args.nb

        if args.astarx != None:
            self.astarx = args.astarx

        if args.astary != None:
            self.astary = args.astary

        if args.astarz != None:
            self.astarz = args.astarz

        if args.bstarx != None:
            self.bstarx = args.bstarx

        if args.bstary != None:
            self.bstary = args.bstary

        if args.bstarz != None:
            self.bstarz = args.bstarz

        if args.cstarx != None:
            self.cstarx = args.cstarx

        if args.cstary != None:
            self.cstary = args.cstary

        if args.cstarz != None:
            self.cstarz = args.cstarz

        if args.spotwidth != None:
            self.spotwidth = args.spotwidth


        if args.delete_sim_dp != None:
            self.delete_sim_dp = args.delete_sim_dp
            
        if args.order != None:
            self.order = args.order

        if args.noiseflag != None:
            self.noiseflag = args.noiseflag

        if args.latticeflag != None:
            self.latticeflag = args.latticeflag

        if args.cx != None:
            self.cx = args.cx

        if args.cy != None:
            self.cy = args.cy

        if args.nth != None:
            self.nth = args.nth

        if args.nq != None:
            self.nq = args.nq

        if args.sym_filter != None:
            self.sym_filter = args.sym_filter

        if args.nxcrop != None:
            self.nxcrop = args.nxcrop

        if args.nycrop != None:
            self.nycrop = args.nycrop

        if args.shiftx != None:
            self.shiftx = args.shiftx

        if args.shifty != None:
            self.shifty = args.shifty

        if args.rebin != None:
            self.rebin = args.rebin

        if args.nl != None:
            self.nl = args.nl

        if args.padf_nth != None:
            self.padf_nth = args.padf_nth

        if args.padf_nq != None:
            self.padf_nq = args.padf_nq

        if args.padf_nr != None:
            self.padf_nr = args.padf_nr

        if args.padf_qmax != None:
            self.padf_qmax = args.padf_qmax

        if args.padf_rmax != None:
            self.padf_rmax = args.padf_rmax

        if args.section1 != None:
            self.section1 = args.section1

        if args.section2 != None:
            self.section2 = args.section2

        if args.section3 != None:
            self.section3 = args.section3
        
        if args.padfplot_theta != None:
            self.theta = args.padfplot_theta

        if args.padfplot_r != None:
            self.r = args.padfplot_r

        if args.padfplot_r2 != None:
            self.r2 = args.padfplot_r2



        #
        # modify defaults
        #
        if self.dp_path == "Default":
            self.dp_path = self.path

        if self.dp_tag == "Default":
            self.dp_tag = self.tag

        if (self.corr_diff_flag == 1) and ( (self.corr_bg_sub_flag == 1) or (self.corr_bg_calc_flag == 1) ):
            print "corr_diff_flag not compatible with subtracting a background correlation. corr_bg flags turned off (=0)"
            self.corr_bg_sub_flag = 0
            self.corr_bg_calc_flag = 0
