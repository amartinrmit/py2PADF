

import numpy as np
import all_params as ap
import padflib
import correlation
import diffraction
import io
import filters
import sys

import matplotlib.pyplot as plt

import os
# import plotting package
import glob



class padfpy:

    def __init__(self, configname="config.txt"):

        self.configname = configname

        self.all_params = ap.params()

        # test that the config file exists !

        #self.all_params.read_config_file( configname )

        self.diffraction = diffraction.diffraction()
        self.correlation = correlation.correlation()
        self.padf = padflib.padf()


        # create an instance of the diffraction class, the correlation class, io class, plotting class etc.


    def run(self):

        # read command line parameters
        self.all_params.get_commandline_args()
           
        # read config file
        self.all_params.read_config_file( self.all_params.args.config[0] )
        print "DEBUG <padfpy.py> all_params nx, ny", self.all_params.nx, self.all_params.ny        

        # read command line parameters and resolve conflicts (with warnings)
        self.all_params.update_with_commandline_args()

        # udpate parameters in all classes
        self.update_diffraction_parameters()
        self.update_correlation_parameters()
        self.update_padf_parameters()
        self.update_rebinned_array_dims()
        self.update_rebinned_pixel_width()

        # set filename
        self.set_inputType()

        self.correlation.set_diffraction_parameters( self.diffraction )

        # get array dimension from file
##        print "DEBUG <padfpy.py; run> ", self.flist[0]
        if (self.all_params.nx == -1) and (self.all_params.ny==-1):
            s = self.get_array_dimensions_from_file()
            self.set_nx( s[0], s[1])
        else:
            self.set_nx( self.all_params.nx, self.all_params.ny)

        print "DEBUG <padfpy.py> nxorig, nyorig", self.correlation.nxorig, self.correlation.nyorig

        # check path and make it if not exists
        if not os.path.exists( self.all_params.path ):
            os.makedirs( self.all_params.path )


        #
        # Write all parameters
        #
        self.all_params.write_all_params()


        # create high pass filter and gaussian (filters.py ???)
        if self.all_params.hpfilter_flag == 1:
            if self.all_params.hp_filter <= 0:
                print "<padfpy.py> hp_filter can not be applied. Check that hp_filter is a positive value."
                self.all_params.hpfilter_flag = 0
            else:
                self.hpf = filters.high_pass_circle_filter( self.correlation.nx, self.all_params.hp_filter )
#                plt.plot( self.hpf)
#                plt.draw()
#                plt.show()
                self.correlation.corr_hp = self.correlation.radial_profile_to_correlation( rad=self.hpf )
                io.write_dbin( self.correlation.path+self.correlation.tag+"_hpfilter.dbin", self.correlation.corr_hp )
#                print "DEBUG <padfpy.py> corr_hp min max", np.min(self.correlation.corr_hp), np.max(self.correlation.corr_hp)
#                sys.exit()


        if self.all_params.g_flag == 1:
            if self.all_params.lp_filter <= 0:
                print "<padfpy.py> Low pass (Gaussian) filter can not be applied. Check that low passfilter is a positive value."
                self.all_params.g_flag = 0
            else:
                self.hpf = filters.gaussian_filter( self.correlation.nx, self.all_params.lp_filter )
                self.correlation.corr_lp = self.correlation.radial_profile_to_correlation( rad=self.hpf )
                io.write_dbin( self.correlation.path+self.correlation.tag+"_lpfilter.dbin", self.correlation.corr_lp )
        
        # calculate mask correlation
        if self.all_params.mask_flag == 1:
            self.correlation.calculate_mask_correlation()
        elif (self.all_params.mask_flag==0) and (self.all_params.corr_shift_flag==1):
            # (if there is a shift, but no imported mask, we may need to correct for zero padding at edge of array)
            self.correlation.calculate_mask_correlation( dontreadfile=True )   
            self.all_params.mask_flag = 1

        # calculate diffraction patterns if required in advance (simulation)
        if self.all_params.dp_flag > 0:
            self.diffraction.set_dp_inputType()
            print "debug diffraction input type:", self.diffraction.pdbInputType

        if self.all_params.dp_flag == 1:
            self.diffraction.generate_diffraction_patterns( self.all_params.npatterns )
            if self.diffraction.noiseflag == 1:
                self.correlation.flist = self.diffraction.get_file_list( self.diffraction.path,\
                                                                             self.diffraction.tag, "noise.dbin")
                self.flist = self.correlation.flist
            else:
                self.correlation.flist = self.diffraction.get_file_list( self.diffraction.path,\
                                                                             self.diffraction.tag, "diffraction.dbin")
                self.flist = self.correlation.flist


        # set the input path, if required
        self.correlation.set_inputpath()
        if self.all_params.inputtail != "None":
            corrtail_current = self.all_params.inputtail
        else: 
            corrtail_current = "_padfcorr_correlation_sum.dbin"

        # calculate the correlation function
        if self.all_params.corr_flag == 1:        
            self.correlation.calculate_correlation()
            self.correlation.set_inputpath()
            corrtail_current = "_padfcorr_correlation_sum.dbin"
            

        # calculate the correlation background if required
        print "DEBUG <padfpy.py> corr_bg_calc_flag", self.all_params.corr_bg_calc_flag
        if self.all_params.corr_bg_calc_flag == 1:        
            self.correlation.bg_estimate = True
            self.correlation.calculate_correlation( "_padfXcorr_correlation_bgsum.dbin")


        # make correlations to correlation function
        if self.all_params.corr_bg_sub_flag == 1:
            bgname = self.correlation.inputpath + self.correlation.inputtag +  "_padfXcorr_correlation_bgsum.dbin"
            corrtail_in = corrtail_current
            corrtail_current = "_padfcorr_correlation_sum_bgsub.dbin"
            self.correlation.subtract_correlation( bgname, ftailin=corrtail_in,\
                             ftailout=corrtail_current)

        if self.all_params.mask_flag ==1:
            corrtail_in = corrtail_current
            corrtail_current = "_padfcorr_correlation_sum_maskcorrected.dbin"
            self.correlation.mask_correction(ftailin=corrtail_in, ftailout=corrtail_current)

        if self.all_params.g_flag == 1:
            corrtail_in = corrtail_current
            corrtail_current = "_padfcorr_correlation_sum_lpfilter.dbin"
            self.correlation.filter( self.correlation.corr_lp, ftailin=corrtail_in, ftailout=corrtail_current)

        if self.all_params.hpfilter_flag == 1:
            corrtail_in = corrtail_current
            corrtail_current = "_padfcorr_correlation_sum_hpfilter.dbin"
            self.correlation.filter( self.correlation.corr_hp, ftailin=corrtail_in, ftailout=corrtail_current)

        if self.all_params.symmetry_flag == 1:
            corrtail_in = corrtail_current
            corrtail_current = "_padfcorr_correlation_sum_symfilter.dbin"
            self.correlation.symmetry_filter( self.correlation.sym_filter, ftailin=corrtail_in, ftailout=corrtail_current)

        print "debug about to do padf"

        # calculate the padf (with units corrections if required)
        if self.all_params.padf_flag == 1:
#            if not self.correlation.correlation_sum_name == None:
            self.padf.corrfile = self.correlation.inputpath + self.correlation.inputtag + corrtail_current
 #           else:
 #               self.padf.corrfile_default()
            self.padf.calculate_padf()

        # make plots
            if (self.all_params.section1>0) and (self.all_params.section1<6):
                self.padf.section = self.all_params.section1
                self.padf.calculate_padfplot( "_padfplot_section1_config.txt" )
            if (self.all_params.section2>0) and (self.all_params.section2<6):
                self.padf.section = self.all_params.section2
                self.padf.calculate_padfplot( "_padfplot_section2_config.txt" )
            if (self.all_params.section3>0) and (self.all_params.section3<6):
                self.padf.section = self.all_params.section3
                self.padf.calculate_padfplot( "_padfplot_section3_config.txt" )


        if (self.all_params.dp_flag == 1) and (self.all_params.delete_sim_dp==1):
            self.diffraction.remove_diffraction_patterns( self.all_params.npatterns )


# read in all the parameters from a configuration file
        
# write relevant parameters into the parameters for specific sections
        
    def update_correlation_parameters( self):

        self.correlation.wl = self.all_params.wavelength
        self.correlation.pixel_width = self.all_params.pixelwidth
        self.correlation.detector_z  = self.all_params.detector_z
        self.correlation.path = self.all_params.path
        self.correlation.tag = self.all_params.tag
        self.correlation.inputpath = self.all_params.inputpath
        self.correlation.inputtag = self.all_params.inputtag
        self.correlation.nth = self.all_params.nth
        self.correlation.nx = self.all_params.nx
        self.correlation.nthreads = self.all_params.nthreads
        self.correlation.npatterns = self.all_params.npatterns
        self.correlation.mask_flag = self.all_params.mask_flag
        self.correlation.maskname = self.all_params.maskname                 
        self.correlation.crop_flag = self.all_params.crop_flag               ### commented lines missing from all_params
#        self.correlation.nxcrop = self.all_params.nxcrop
#        self.correlation.nycrop = self.all_params.nycrop
        self.correlation.set_nxcrop_and_nycrop( nxcrop=self.all_params.nxcrop, nycrop=self.all_params.nycrop)
        self.correlation.cx = self.all_params.cx
        self.correlation.cy = self.all_params.cy
        self.correlation.dp_shift_flag = self.all_params.dp_shift_flag
        self.correlation.shiftx = self.all_params.shiftx
        self.correlation.shifty = self.all_params.shifty
        self.correlation.rebin = self.all_params.rebin
        self.correlation.sym_filter = self.all_params.sym_filter
        self.correlation.nstart = self.all_params.nstart
        self.correlation.h5field = self.all_params.h5field
        self.correlation.pwx = self.all_params.pwx
        self.correlation.pwy = self.all_params.pwy
        self.correlation.set_rebin_flag()
        self.correlation.nxorig = self.all_params.nx
        self.correlation.nyorig = self.all_params.ny
        self.correlation.dp_flag = self.all_params.dp_flag
        if self.all_params.corr_diff_flag == 1:
            self.correlation.diffcorrflag = True

    def update_padf_parameters( self ):
        
        self.padf.path = self.all_params.path
        self.padf.tag = self.all_params.tag
        self.padf.nl = self.all_params.nl
        self.padf.nr = self.all_params.padf_nr
        self.padf.nth = self.all_params.padf_nth
        self.padf.nq = self.all_params.padf_nq
        self.padf.qmax = self.all_params.padf_qmax
        self.padf.rmax = self.all_params.padf_rmax
        self.padf.wl = self.all_params.wavelength
        self.padf.section1 = self.all_params.section1
        self.padf.section2 = self.all_params.section2
        self.padf.section3 = self.all_params.section3
        self.padf.theta = self.all_params.theta
        self.padf.r = self.all_params.r
        self.padf.r2 = self.all_params.r2

    def update_diffraction_parameters( self ):

        self.diffraction.wl = self.all_params.wavelength
        self.diffraction.pixel_width = self.all_params.pixelwidth
        self.diffraction.detector_z  = self.all_params.detector_z
        self.diffraction.path = self.all_params.path
        self.diffraction.tag = self.all_params.tag
        self.diffraction.pdbpath = self.all_params.samplepath
        self.diffraction.pdbtag = self.all_params.sampletag
        self.diffraction.beamArea = self.all_params.beamArea
        self.diffraction.nph = self.all_params.nphotons
        self.diffraction.rx = self.all_params.rx
        self.diffraction.ry = self.all_params.ry
        self.diffraction.rz = self.all_params.rz
        self.diffraction.ry = self.all_params.ry
        self.diffraction.nx = self.all_params.nx
        self.diffraction.nrho = self.all_params.nx/2
        self.diffraction.pdbname = self.all_params.samplepath + self.all_params.sampletag 
        self.diffraction.scaledflag = self.all_params.scaledflag
        self.diffraction.radial_average_flag = self.all_params.radial_average_flag
        self.diffraction.npatterns = self.all_params.npatterns
        self.diffraction.nthreads = self.all_params.nthreads
        self.diffraction.gflag = self.all_params.g_flag
        self.diffraction.gwid = self.all_params.gwid
        self.diffraction.delete_sim_dp = self.all_params.delete_sim_dp
        self.diffraction.order = self.all_params.order
        self.diffraction.noiseflag = self.all_params.noiseflag
        self.diffraction.latticeflag = self.all_params.latticeflag
        self.diffraction.spotwidth = self.all_params.spotwidth
        self.diffraction.na = self.all_params.na
        self.diffraction.nb = self.all_params.nb
        self.diffraction.nc = self.all_params.nc
        self.diffraction.astarx = self.all_params.astarx
        self.diffraction.astary = self.all_params.astary
        self.diffraction.astarz = self.all_params.astarz
        self.diffraction.bstarx = self.all_params.bstarx
        self.diffraction.bstary = self.all_params.bstary
        self.diffraction.bstarz = self.all_params.bstarz
        self.diffraction.cstarx = self.all_params.cstarx
        self.diffraction.cstary = self.all_params.cstary
        self.diffraction.cstarz = self.all_params.cstarz



    def get_file_list( self, path, tag ):

        try:
            #                print "DEBUG (set_inputType) sample string: ", self.all_params.samplepath + self.all_params.sampletag+"*" 
            self.flist = sorted(glob.glob( path + tag+"*" ), key=os.path.getmtime )
            if len(self.flist) != 0:
                self.correlation.flist = self.flist
            else:
                print "<padfpy.py> no diffraction files found. Check the samplepath and sampletag"
                print "samplepath: ", path
                print "sampletag: ", tag
                sys.exit()
        except:
            print "Error (padfpy.py). Files could not be read from the samplepath. Check samplepath and sampletag"
        

    def set_inputType( self ):

        fext = os.path.splitext( str(self.all_params.samplepath + self.all_params.sampletag) )[-1]
        print "DEBUG <padfpy.py; set_inputType> fext", fext

        if self.all_params.dp_flag > 0:
            self.inputType = 'dir'
        elif len(fext) == 0:
            self.inputType = 'dir'
            self.get_file_list( self.all_params.samplepath, self.all_params.sampletag )

        elif fext > 1:
            if (fext=='.ser') or (fext=='.tif'):
                self.inputType = 'multifile'
                self.correlation.datafilename = self.all_params.samplepath + self.all_params.sampletag
            elif (fext=='.list') or (fext=='.txt'):
                self.inputType = 'dir'
                f = open( self.all_params.samplepath + self.all_params.sampletag, 'r')
                self.flist = []
                self.xshift_list = []
                self.yshift_list = []
                for line in f:
                    bits = line.split()
                    self.flist.append( bits[0] )
                    if len(bits)>1: self.xshift_list.append( int(bits[1]) )
                    if len(bits)>2: self.yshift_list.append( int(bits[2]) )
                    if (len(self.xshift_list)>0) and not (len(bits)>1):
                        print "Error: line in file list in missing xshift value (expected 2nd entry on line)"
                        sys.exit()
                    if (len(self.yshift_list)>0) and not (len(bits)>2):
                        print "Error: line in file list in missing yshift value (expected 3nd entry on line)"
                        sys.exit()
                f.close()
                self.correlation.flist = self.flist
                self.correlation.xshift_list = self.xshift_list
                self.correlation.yshift_list = self.yshift_list

 #       print "DEBUG <padfpy.py; set_inputType> inputType", self.inputType
        self.correlation.inputType = self.inputType
            

    def set_wavelength( self, wl ):

        self.all_params.wavelength = wl
        self.diffraction.wl = wl
        self.correlation.wl = wl
        self.padf.wl = wl
        
    def get_wavelength( self ):
        return self.all_params.wavelength


    def set_pixel_width( self, pw ):

        self.all_params.pixel_width = pw
        self.diffraction.pixel_width = pw
        self.correlation.pixel_width = pw

    def set_nx( self, nx, ny=-1 ):
        
        print "DEBUG <padfpy.py> set_nx nx, ny", nx, ny
        if ny==-1:
            ny = nx

        if self.all_params.crop_flag == 0:
            self.all_params.nx = nx 
            self.diffraction.nx = nx
            self.correlation.nx = nx
            # self.all_params.ny = ny
            # self.diffraction.ny = ny
            self.correlation.ny = ny
            self.correlation.nxorig = nx
            self.correlation.nyorig = ny

        elif self.all_params.crop_flag==1:
            self.all_params.nx = self.all_params.nxcrop
            self.diffraction.nx = self.all_params.nxcrop
            self.correlation.nx = self.all_params.nxcrop
            self.correlation.ny = self.all_params.nycrop
            self.all_params.ny = self.all_params.nycrop
            self.diffraction.ny = self.all_params.nycrop
 
            
        self.update_rebinned_array_dims()

            

    def update_rebinned_array_dims( self ):

        if self.all_params.rebin > 1:
            self.all_params.nx = self.all_params.nx / self.all_params.rebin
            self.diffraction.nx = self.diffraction.nx / self.all_params.rebin
            self.correlation.nx = self.correlation.nx / self.all_params.rebin
            # self.all_params.ny = self.all_params.ny / self.all_params.rebin
            # self.diffraction.ny = self.diffraction.ny / self.all_params.rebin
            self.correlation.ny = self.correlation.ny / self.all_params.rebin
        
    def update_rebinned_pixel_width( self ):
        if self.all_params.rebin > 1:
            self.diffraction.pixel_width = self.diffraction.pixel_width * self.all_params.rebin
            self.correlation.pixel_width = self.correlation.pixel_width * self.all_params.rebin
            self.correlation.pwx = self.correlation.pwx * self.all_params.rebin
            self.correlation.pwy = self.correlation.pwy * self.all_params.rebin

        
    def get_pixel_width( self ):
        return self.all_params.pixel_width


    def get_array_dimensions_from_file( self ):

        if self.inputType == 'dir':
            image = io.read_image( self.flist[0] )
        elif self.inputType == 'multifile':
            image = io.read_image( self.all_params.samplepath+self.all_params.sampletag, 0 )
        elif self.inputType == 'simulation':
            print "<padfpy.py: get_array_dimension_from_file(): Simulation diffraction patterns not implemented yet."
            sys.exit()

        return image.shape
