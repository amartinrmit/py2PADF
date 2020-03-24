


import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import struct
import array
import time
import multiprocessing as mp
import amser as ser


etype = 'fem'
stype = 'C9' # 'C5' #    'unactivated' #     
dcode =  '12.25.40' #'13.52.36' #  '12.35.05' #        '12.55.43' # 
sample_path = "/Volumes/DataStore1/Data/Monash/Carbon_dec16/"+stype+"/"+etype+"/"
sample_fname = dcode+"_Spectrum_image_1.ser"

#plt.ion()

nstart = 260
npatterns = 900 - nstart
for i in np.arange( npatterns ) + nstart:

    data, calp, adl, sh = ser.read_ser( sample_path+sample_fname, [i] )

    plt.imshow( np.abs(data.dlist[0])**0.3 )
    plt.draw()
    plt.show()
    #time.sleep(0.1)

#plt.ioff()
    
    
