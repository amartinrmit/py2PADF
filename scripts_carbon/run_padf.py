
import numpy as np
import matplotlib.pyplot as plt
import all_params as ap
import carbonbg
import shutil
import os

class padf_sample:

    def __init__( self, sample="None" ):

        self.name = sample
        self.runs = []
        self.xshift =[]
        self.yshift =[]



#
# setup the config file
#
configfile = "config_monash_carbon_dec17.txt"
app = ap.params()
app.read_config_file( configfile )


#
# set up location variables
#
outpath = "/Users/e38496/Work/Research/Experiments/Monash/2017/Carbon_dec17/Results/"
outtag = "meancorrected_new"
tag_end = "_20_20_20nm_20nm_CL_230mm_CA_5um_0_5s_1_npat400_hp55_centred_"+outtag

#
# Set up some variables 
# lists of samples, runs, beam centres
#
samplist = []
i = 0

#samplist.append( padf_sample("C25"))
#samplist[i].runs = ["one", "two", "three", "four","five"]
#samplist[i].xshift = [1, 0,  0,  0 ,1]
#samplist[i].yshift = [11, 22, 21, 21 ,21]
#i += 1

#samplist.append( padf_sample("C45"))
#samplist[i].runs = ["one", "two", "three","four","five"]
#samplist[i].xshift = [3,0,0,0,0]
#samplist[i].yshift = [1,2,0,1,0]
#i += 1

#samplist.append( padf_sample("c0_2nd_plasma"))
#samplist[i].runs = ["one", "two", "three", "four","five"]
#samplist[i].xshift = [3,2,1,2,6]
#samplist[i].yshift = [21,16,7,12,29]
#i += 1

samplist.append( padf_sample("c0_2nd_plasma"))
samplist[i].runs = ["six", "seven", "eight", "nine","ten"]
samplist[i].xshift = [2,2,2,2,2]
samplist[i].yshift = [16,17,17,18,18]
i += 1

#samplist.append( padf_sample("C5"))
#samplist[i].runs = ["one", "two", "three", "four", "five"] #, "six", "seven", "eight", "nine","ten"]
#samplist[i].xshift = [1,0,1,1,1] #,1,1,1,0,0]
#samplist[i].yshift = [0,4,2,0,-2] #,-2, 1,-1,-1,-1]
#i += 1

samplist.append( padf_sample("glassy_carbon") )
samplist[i].runs = ["three", "four", "five"] #["one", "two", "three", "four", "five"] #, "six", "seven", "eight", "nine","ten"]
samplist[i].xshift = [3,3,3] #[2,5,3,3,3] #,2,4,4,4,3] 
samplist[i].yshift = [-2, -1, -1] #[21,-6,-2,-1,-1] #,8,1,-1,-2,-1] 
i += 1

samplist.append( padf_sample("C45"))
samplist[i].runs = ["one", "two", "three","four","five"]
samplist[i].xshift = [3,0,0,0,0]
samplist[i].yshift = [1,2,0,1,0]
i += 1

samplist.append( padf_sample("glassy_carbon") )
samplist[i].runs = ["six", "seven", "eight", "nine","ten"]
samplist[i].xshift = [2,4,4,4,3] 
samplist[i].yshift = [8,1,-1,-2,-1] 
i += 1



#
# Set up background class
#
datapath = "/Volumes/DataStore1/Data/Monash/Carbon_dec17/dbin_data/"
n = 399
bg = carbonbg.bgsub( root=datapath, npatterns=n )


#
# loop over all the samples/runs
#
for sample in samplist:

    for i in np.arange( len(sample.runs) ):
   
        #
        # Generate the bg
        #
        tag = sample.runs[i]+"_20_20_20nm_20nm_CL_230mm_CA_5um_0_5s_1"
        bg.update_dataset( sample.name, outtag="_"+outtag )
        bg.update_tag( tag, outtag="_"+outtag )
            
        if outtag[:2] == "bg":
            print "<run_padf.py> generate bg corrected patterns"
            bg.generate_bg_corrected_patterns()
            app.corr_bg_calc_flag = 1
            app.corr_bg_sub_flag = 1

        if outtag[:4] == "mean":
            print "<run_padf.py> generate mean corrected patterns"
            bg.data_set_mean()
            bg.generate_mean_corrected_patterns()
            app.corr_bg_calc_flag = 0
            app.corr_bg_sub_flag = 0


        #
        # Set up the config file
        #
        configname = outpath+tag+"_"+sample.name+"_"+sample.runs[i]+"config.txt"
        app.path = outpath + sample.name +"/FEM/"+sample.runs[i]+"_20_20_20nm_20nm_CL_230mm_CA_5um_0_5s_1_"+outtag+"/"
        app.tag = sample.runs[i]+tag_end 
        app.samplepath = bg.outpath
        print "DEBUG <run_padf.py> samplepath", app.samplepath
        app.sampletag  = sample.runs[i]+"_20_20_20nm_20nm_CL_230mm_CA_5um_0_5s_1"
        app.shiftx = sample.xshift[i]
        app.shifty = sample.yshift[i]
        app.npatterns = n
        app.nthreads = 3
        app.write_all_params( configname )

        #
        # run padf
        #
        os.system( "padf_0.1 -c "+configname )

        #
        # delete the files
        #
        shutil.rmtree( bg.outpath )
