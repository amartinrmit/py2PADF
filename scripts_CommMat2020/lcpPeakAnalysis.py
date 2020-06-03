"""

lcpPeakAnalysis.py

- reads 1D intensity profile from txt file
- compares intensity to predicted peaks from cubic, hexagonal, laminar phases 

Andrew Martin Aug 2019

"""



import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import scipy.interpolate as spi


def namestr(obj, namespace):
    return [name for name in namespace if namespace[name] is obj]

def scatterbrain_q_calibration( rbins ):
    # scatterbrain calibration
    #wavelength = 0.619921 #1.033216 #e-10
    #detector_z = 0.962478     #0.959563 #
    #pixelwidth = 0.000172 
    rstart = 16.5  

    # what I think it should be based on LCP / LCP2 log files
    wavelength = 1.033216 
    detector_z = 0.959563 
    pixelwidth = 0.000172   
    # rstart = rmin

    qindices = rstart + np.arange(rbins )
    q = 2.0 * np.pi * (2.0/wavelength) * np.sin( np.arctan( qindices*(pixelwidth/detector_z) ) / 2.0 )
    return q


##################################################################################
#
# DEFINE THE LOCATION AND NAME OF THE INPUT FILE
#
path = "/Users/e38496/Work/Research/Experiments/AS/SAXS/SAXS_feb18/Results/radial_1D_sums/"
#fname = "scan27_MO_DL2_buffer_nstart_0_n_1000_average.txt"
fname = "scan19_MO_mQ_nstart_0_n_4488_average.txt"
#fname = "run9_plate22_nstart_0_n_3072_average.txt"


# Load the data from file
data = np.loadtxt( path+fname )

# Extract the q-values
dataqflag = False
if dataqflag == True:
    qcalib = data[:,0]
else:
    qcalib = scatterbrain_q_calibration( data.shape[0] )


#
# Specify range of q-values for peak search
#
qmin = 0.05
qmax = 0.35



# indices of q values in search range 
iq = np.where( (qcalib[:] > qmin)*(qcalib[:]<qmax))

# find the peaks
peakind = signal.find_peaks_cwt(data[iq,1][0], np.arange(0.1,1), min_snr=1.1)


print "\n"
print "INPUT PARAMEToERS:"
print "Path :", path
print "Input filename:", fname
print "Peak search range :", qmin, " to", qmax


#
# Find the peaks again, but interpolate to improve the accuracy
#
samplingfactor = 100.0
y = data[iq,1]
x = np.linspace(0,10,num=y.size)
dinterp = spi.interp1d( x, y, kind='cubic')
xnew = np.linspace(0,10,num=y.size*samplingfactor)
ynew = np.squeeze(dinterp(xnew))
peakind_interp = signal.find_peaks_cwt(ynew, np.arange(0.1,1), min_snr=1.1)

y = qcalib[:]
x = np.linspace(0,10,num=y.size)
dinterp = spi.interp1d( x, y, kind='cubic')
xnew = np.linspace(0,10,num=y.size*samplingfactor)
ynew = np.squeeze(dinterp(xnew))
qcalib_interp = ynew
iq_interp = np.where( (qcalib_interp[:] > qmin)*(qcalib_interp[:]<qmax))

#
# Plot the original data
#
y0, y1 = 50, 200
print "Plot q-range (pixels):", y0, y1
plt.plot( qcalib[y0:y1], data[y0:y1,1] )

#
# Add lines to plot where peaks were found
#
for i in np.arange(len(peakind)):
    plt.axvline( x=qcalib_interp[iq_interp[0][peakind_interp[i]]]+0.0001, color='g' )
plt.xlabel( r'q (nm$^{-1}$)' )
plt.ylabel( r'Intensity (aribrary units)' )

#plt.draw()
#plt.show()

#
# Choose which reference peak to compare phase
#
peaklist = [0, 1]


#
# Set lattice peak corresponding to reference peak
#
latRefPeak = [0, 0]

clist = ['r','k']
lslist = ["--"]*2

# find q value of reference peak
q0list = []
qlist = []
for peak in peaklist: 
    q0 = iq_interp[0][peakind_interp[peak]]
    q = qcalib_interp[q0]
    q0list.append(q0)
    qlist.append(q)


#
# Define the lattice ratios
#
pn3m = np.array([np.sqrt(2), np.sqrt(3), np.sqrt(4), np.sqrt(6),\
                     np.sqrt(8),np.sqrt(9), np.sqrt(10),np.sqrt(11)] )
###                     np.sqrt(12), np.sqrt(14), np.sqrt(16), np.sqrt(17) ])

hex = np.array( [1, np.sqrt(3), 2, np.sqrt(12), 3, np.sqrt(10) ] ) 

cubic = np.array( [1, np.sqrt(2), np.sqrt(3), 2] )

ia3d = np.array([np.sqrt(6), np.sqrt(8), np.sqrt(14), np.sqrt(16), np.sqrt(20) ] )

laminar = np.array([1,2,3,4,5,6])

#
# Select the lattice(s) to be compared
#
lattice = [pn3m, cubic]

lattice_names = [ namestr(l, globals())[0] for l in lattice ]

print "\nLattice types:", lattice_names
print "Reference peaks:", peaklist
print "Lattice peak to match to reference:", latRefPeak
print "Lattice line colors:", clist
print "Lattice line style:", lslist



# calculate q values of 

for ip in np.arange(len(peaklist)):
    
    qlat = qlist[ip]/(lattice[ip][latRefPeak[ip]])
    for i in np.arange(lattice[ip].size):
        plt.axvline( x=qlat*lattice[ip][i], color=clist[ip], linestyle=lslist[ip] )



print "\nOUTPUT:"
print "peakind, qcalib", peakind, qcalib[iq[0][peakind]]

#calculate the lattice parameters
for ip in np.arange(len(peaklist)):
    qlat = qlist[ip]/(lattice[ip][latRefPeak[ip]]) / (2.0*np.pi)
    latparam = 1.0/qlat
    print "lattice no. , name, spacing :", ip, lattice_names[ip], latparam
    print qlist[ip]


plt.draw()
plt.show()
