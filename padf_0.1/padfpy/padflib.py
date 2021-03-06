
import os
import time
import io

class padf:

    def __init__(self, path="None", tag="None", nl=2, nr=100, nth=100, nq=100, qmax=0.0, rmax=0.0, corrfile="None",
                 wl=0.0,units_flag=0,section=1,theta=0,r=0,r2=0):


        self.path = path
        self.tag = tag
        self.nl = nl
        self.nlmin = nlmin
        self.nr = nr
        self.nth = nth
        self.nq = nq
        self.qmax = qmax
        self.rmax = rmax
        self.corrfile = corrfile
        self.wl = wl
        self.units_flag = units_flag

        if self.corrfile == "None":
            self.corrfile = self.path+self.tag+"_padfcorr_correlation_sum.dbin"

        # padf plot parameters
        self.section = section
        self.theta   = theta
        self.r       = r
        self.r2      = r2
        

    def write_padf_config(self,pfname):

        fcorr = open(pfname,'w')
        fcorr.write("correlationfile = "+self.corrfile+'\n' )
        # fcorr.write("correlation_sigma_file = "+path+tag+csig_name+'\n' )
        fcorr.write("outpath = "+self.path+'\n')
        fcorr.write("tag = "+self.tag+"_padf2"+'\n')
        fcorr.write("wavelength =  "+str(self.wl)+'\n'  )
        fcorr.write("nthq =  "+str(self.nth)+'\n'  )
        fcorr.write("nq =  "+str(self.nq)+'\n'  )
        fcorr.write("nr =  "+str(self.nr)+'\n'  )
        fcorr.write("nl =  "+str(self.nl)+'\n'  )
        fcorr.write("nlmin =  "+str(self.nlmin)+'\n'  )
        fcorr.write("qmax =  "+str(self.qmax)+'\n'  )
        fcorr.write("rmax =  "+str(self.rmax)+'\n'  )
        fcorr.close()

    def write_padfplot_config(self,pfname):

        fcorr = open(pfname,'w')
        fcorr.write("correlationfile = "+self.corrfile+'\n' )
        # fcorr.write("correlation_sigma_file = "+path+tag+csig_name+'\n' )
        fcorr.write("outpath = "+self.path+'\n')
        fcorr.write("tag = "+self.tag+"_padf2"+'\n')
        fcorr.write("wavelength =  "+str(self.wl)+'\n'  )
        fcorr.write("nthq =  "+str(self.nth)+'\n'  )
        fcorr.write("nq =  "+str(self.nq)+'\n'  )
        fcorr.write("nr =  "+str(self.nr)+'\n'  )
        fcorr.write("nl =  "+str(self.nl)+'\n'  )
        fcorr.write("nlmin =  "+str(self.nlmin)+'\n'  )
        fcorr.write("qmax =  "+str(self.qmax)+'\n'  )
        fcorr.write("rmax =  "+str(self.rmax)+'\n'  )
        fcorr.write("section = "+str(self.section)+'\n')
        if (self.section==1) or (self.section==4):
            fcorr.write("theta = "+str(self.theta)+'\n')
        if (self.section==2) or (self.section==4) or (self.section==5):
            fcorr.write("r = "+str(self.r)+'\n')
        if (self.section==4):
            fcorr.write("r2 = "+str(self.r2)+'\n')
        fcorr.close()


    def calculate_padf(self):

        print "\nCalculating the padf"
        start = time.clock()
        pfname = self.path+self.tag+"_padf_config.txt"
        print "debug <padflib.py> corrfile", self.corrfile
        self.write_padf_config(pfname)          
        dirpath = os.path.dirname(os.path.realpath(__file__))       
        os.system(dirpath+"/../padf/padf "+pfname ) 
#        os.system("lldb "+dirpath+"/../padf/padf "+pfname ) 
        print "padf took :", time.clock() - start, " seconds"

    def calculate_padfplot(self, configtag="_padfplot_config.txt"):

        print "\nCalculating the padfplot"
        start = time.clock()
        pfname = self.path+self.tag+configtag
        self.write_padfplot_config(pfname)          
        dirpath = os.path.dirname(os.path.realpath(__file__))       
        os.system(dirpath+"/../padfplot/padfplot "+pfname ) 
#        os.system("lldb "+dirpath+"/../padf/padf "+pfname ) 
        print "padf took :", time.clock() - start, " seconds"


    def corrfile_default( self) :
        self.corrfile = self.path+self.tag+"_padfcorr_correlation_sum.dbin"

#
#----------------------------------------------------------------
# Corrections for units and interpretation
#
#
    def correct_units(self):
        re = 2.8179403267*1e-15
        av = 6.0221409e23
        if self.units_flag == 1:

            constants = (dpp.nphotons / dpp.beamArea)*re*re
            # mass_density = 2.33                                # g per cm^3
            # molmass = 58.6934					# g per mol    
            # number_density = (mass_density*1e6*av / molmass)
     
            number_density = (5234.0 / 60927.0)*1e30        # direction from simulation
            number_density *= (rmax*1e10 / float(nx))**(-3)  # a correction for the size of the radial bin
            factor = 1.0/(number_density*constants)
            total_correction = nth * factor*factor
     
            print number_density, factor, total_correction

# needs correct io routines...
            cname = path+tag+"_padf2_padf.dbin"
            padf = read_correlation( cname , 0)
            padf *= total_correction
            io.write_dbin( path+tag+"_padf2_padf_units_scaled.dbin", padf )
