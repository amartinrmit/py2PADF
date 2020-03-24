


import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import struct
import array
import time
import multiprocessing as mp
import amser as ser


def write_dbin( fname, data ):

	f = open( fname, "wb")
	fmt='<'+'d'*data.size
	bin = struct.pack(fmt, *data.flatten()[:] )
	f.write( bin )
	f.close()


etype = 'fem'
stype = 'Unactivated' # 'C5' #    'unactivated' #     
dcode =  '13.04.29' #'13.52.36' #  '12.35.05' #        '12.55.43' # 
sample_path = "/Volumes/DataStore1/Data/Monash/Carbon_dec16/"+stype+"/"+etype+"/"
sample_fname = dcode+"_Spectrum_image_1.ser"

outpath  = "/Users/e38496/Work/Research/Experiments/Monash/2016/Carbon_dec16/Results/datasum/"
#plt.ion()

nstart = 0
npatterns = 900 - nstart
freq = 300
for i in np.arange( npatterns ) + nstart:

    data, calp, adl, sh = ser.read_ser( sample_path+sample_fname, [i] )

    if i == nstart:
        datasum = data.dlist[0].copy()
    else:
        datasum += data.dlist[0].copy()

    if i%freq == 0:
        fname = sample_fname[:-4]+"_"+str(i)+"_sum.dbin"
        write_dbin( outpath+fname, datasum )


fname = sample_fname[:-4]+"_"+str(nstart)+"_"+str(npatterns)+"_sum.dbin"
write_dbin( outpath+fname, datasum )


#    plt.imshow( np.abs(data.dlist[0])**0.3 )
#    plt.draw()
#    plt.show()
    #time.sleep(0.1)

#plt.ioff()
    
    
