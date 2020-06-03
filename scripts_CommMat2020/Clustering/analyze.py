from os import listdir
from os.path import isfile, join
import os
from scipy import signal
from scipy import linalg
import numpy as np
from Laplace import GetDW
from Laplace import mknn
import matplotlib.pyplot as plt
from math import sqrt
from sklearn.cluster import KMeans
import glob
import sys
import matplotlib
matplotlib.rcParams.update({'font.size': 16})


class Data:
    path = ""
    names = []
    Peaks = []
    Label = []

    def histo(self):
        Histo = {}
        for elem in self.Peaks:
            l = len(elem)
            if l in Histo:        
                Histo[l] += 1
            else:
                Histo.update({l : 1})
        
        return Histo

    def SortByLabel(self):
        # Sorts data by Cluster Label, and withing that label by number of Peaks
        TMP = []
        NumClust = {}
        for i, elem in enumerate(self.Label):
            TMP.append([self.names[i], self.Peaks[i], self.Label[i]])
            if elem in NumClust:
                NumClust[elem][0] += 1
            else:
                NumClust.update({elem : [1, 0]})
        # NumClust[elem][1] indicates the postion of the first element of that cluster.
        # NumClust[elem][0] indicates how many elements that cluster has.
        for elem in NumClust:
            if elem > 0:
                NumClust[elem][1] = NumClust[elem-1][0] + NumClust[elem-1][1]

        TMP.sort(key = lambda x : x[2])
        for elem in NumClust:
            if elem == 0:
                TMP[:NumClust[elem][1]] = sorted(TMP[:NumClust[elem][1]], key = lambda x: len(x[1]))
            else:
                TMP[NumClust[elem-1][1]:NumClust[elem][1]] = sorted(TMP[NumClust[elem-1][1]:NumClust[elem][1]], key = lambda x: len(x[1]))

        for i, elem in enumerate(TMP):
            self.names[i] = elem[0]
            self.Peaks[i] = elem[1]
            self.Label[i] = elem[2]

        return NumClust

class cluster_params:

    def __init__( self ):

        self.datapath = "None"
        self.datatag =  "None"
        self.outputpath = "None"
        self.outtag = "None"
        self.sing_val_threshold = 1e-9
        self.sgm = 20.0
        self.sim_metric = "sum"     #other options "prod"
        self.A = 0.1
        self.B = 0.5
        self.nmin = 0
        self.nmax = 100
        self.ncluster_thresh = 10



    def write_cluster_params_to_file( self ):

        f = open( self.outputpath+self.outtag + "_cluster_params_log.txt", 'w' )
        f.write("datapath = "+self.datapath+'\n')
        f.write("datatag  = "+self.datatag+'\n')
        f.write("outputpath  = "+self.outputpath+'\n')
        f.write("outtag  = "+self.outtag+'\n')
        f.write("sing_val_threshold  = "+str(self.sing_val_threshold)+'\n')
        f.write("sgm  = "+str(self.sgm)+'\n')
        f.write("sim_metric  = "+str(self.sim_metric)+'\n')
        f.write("A  = "+str(self.A)+'\n')
        f.write("B  = "+str(self.B)+'\n')
        f.write("nmin  = "+str(self.nmin)+'\n')
        f.write("nmax  = "+str(self.nmax)+'\n')
        f.write("ncluster_thresh  = "+str(self.ncluster_thresh)+'\n')
        f.close()


    def read_clusters_params_from_file( self, fname ):

        f = open( fname, 'r' )

        for line in f:
            
            if line[0] == '#': continue
            if line[0] == '\n': continue

            bits = line.split()
            if bits[0] == "datapath":
                self.datapath = bits[2]

        f.close()


def write_cluster_names_to_file( path, Input, LabelHisto, outpath=None, outtag=None ):


    for elem in LabelHisto:
        print "[ ", elem, " : ", LabelHisto[elem][1], "-", LabelHisto[elem][1] + LabelHisto[elem][0] - 1," ]"
        nstart = LabelHisto[elem][1]
        nstop = LabelHisto[elem][1] + LabelHisto[elem][0] - 1

        if (outpath != None) and (outtag != None):
            fout = open( outpath+outtag+"_"+str(elem)+"_filelist.txt", "w")
            for id in np.arange( nstop-nstart)+nstart:
                
                fout.write( path + Input.names[id]+"\n" )
               
            fout.close()
        else:
            print "Warning! Cluster file names not written to file, because outpath and/or outtag not specified."


def plot_cluster_sums( LabelHisto, cl_sums, threshold, display=False, outpath=None, outtag=None, q=np.zeros(1), tail="_cluster_sums_plot.png",\
                           offset=0):

        num_plots = 0
        for elem in LabelHisto:
            if LabelHisto[elem][0] > threshold:
                num_plots += 1
        if q.size == 1:
            q = np.arange( cl_sums.size)

        # Have a look at the colormaps here and decide which one you'd like:
        # http://matplotlib.org/1.2.1/examples/pylab_examples/show_colormaps.html
        colormap = plt.cm.gist_ncar
        plt.gca().set_color_cycle([colormap(i) for i in np.linspace(0, 0.9, num_plots)])
        plist, labellist = [], []
        fig = plt.figure(figsize=(16,16),dpi=80)
#        for cl in cl_sums:
        i = 0
        for elem in LabelHisto:
            if LabelHisto[elem][0] > threshold:
                p, = plt.plot( q, cl_sums[elem] + i*offset )
                plist.append( p )
                labellist.append( elem )
                i += 1
        plt.legend( tuple(plist), tuple(labellist), ncol=1, loc='center left', bbox_to_anchor=(1,0.5) ) #, prop={'size',10} )
        plt.ylabel( "Intensity (counts)" )
        plt.xlabel( r'q (nm$^{-1})$' ) 
        fig.tight_layout( rect=[0.,0.,0.5,1.0])

        plt.draw()
        if (outpath != None) and (outtag != None):
            plt.savefig( outpath+outtag+tail )

        if display==True:
            plt.show()


def Plot(FileName, A, B):
    # Get X and Y for plotting.
    Xdata = []
    Ydata = []
    ReadData(FileName, Xdata, Ydata, A, B)

    peakind = signal.find_peaks_cwt(Ydata, np.arange(0.1,1,0.2), min_snr=1.1)
    # Plot the data.
    fig = plt.figure(figsize=(15, 10))
    ax = fig.add_subplot(111)
   # ay = fig.add_subplot(111)
    ax.plot(Xdata, Ydata, lw = 2)
   # ay.plot(Xdata, YdataF, lw = 2, color = 'orange')
    ax.set_xlabel('Q, 1/A')
    ax.set_ylabel('Intensity')
    ax.set_title('Diffraction data')
    for elem in peakind:
        plt.axvline(x = Xdata[elem], color='red', linestyle = '--')

    plt.show()

def ReadData(FileName, Qdata, Idata, A, B):
    File = open(FileName)
        # read the content into a list "Rates"
    data = []
    for line in File:
        data.append(line.split())
    
    #print len(data)
    for elem in data:
        tmp_Q = float(elem[0])
        if tmp_Q > A and tmp_Q < B: 
            Qdata.append(tmp_Q)
            tmp_I = float(elem[1])
            Idata.append(tmp_I)

def SumClusters( path, Input, LabelHisto, A, B, average=False ):

    cl_sums = []
    for elem in LabelHisto:
        print "[ ", elem, " : ", LabelHisto[elem][1], "-", LabelHisto[elem][1] + LabelHisto[elem][0] - 1," ]"
        nstart = LabelHisto[elem][1]
        nstop = LabelHisto[elem][1] + LabelHisto[elem][0] - 1

        for id in np.arange( nstop-nstart)+nstart:
            Qdata, Idata = [], []
            ReadData(path + Input.names[id], Qdata, Idata, A, B)
            if id==nstart:
                datasum = np.array( Idata )
            else:
                #print "nstart", nstart, nstop, id, np.array(Idata).size, datasum.size
                datasum += np.array( Idata )

        if average==True:
            datasum /= float(nstop-nstart)
            print "dividing"
        print datasum[50]

        cl_sums.append( datasum ) #/float(nstop-nstart) )
    return cl_sums


def main( cp = cluster_params() ):
    # A is min(Q), B is max(Q) in the scattering data (Q is on x axis).
#    A = 0.1
#    B = 0.5
    A = cp.A
    B = cp.B

    Input = Data()

#    Input.path = "../Clustering/raw_dat(copy)"
#    Input.path = "/Volumes//DataStore1/AS_SAXS_feb18_data/0p6_saxs/raw_dat/"
#    Input.path = "/Volumes/DataStore1/AS_SAXS_feb18_data/0p6_saxs/raw_postexpt"
#    mypath = Input.path + '/'
    #tag = "run15_plate18"
#    tag = "plate4"
    mypath = cp.path + '/'
    tag = cp.tag

    # open a file to log information to
    flog = open( cp.outputpath+cp.outtag+"_clustering_log.txt", "w")

    # Get all files in a directory.
    #    Input.names = [f for f in listdir(mypath) if isfile(join(mypath, f))]
    nmin, nmax = cp.nmin, cp.nmax
#    flist = sorted(glob.glob( mypath+tag+"*.dat" ), key=os.path.getmtime )
#    print "total images:", len(flist)
    Input.names = [ os.path.basename(x) for x in sorted(glob.glob( mypath+tag+"*.dat" )[nmin:nmax], key=os.path.getmtime ) ]
    print len(Input.names)
    if len(Input.names)==0:
        print "No files found. Check the path and tag inputs"
        sys.exit()

    print "Reading the files in the folder..."
    Qlen = 0
    for i, f in enumerate(Input.names):
        #print f
        # Get Intensity and Q values.
        Qdata = []
        Idata = []
        peaks = []
        ReadData(mypath + f, Qdata, Idata, A, B)
        # Find maxima in Idata.
        peakind = signal.find_peaks_cwt(Idata, np.arange(0.1,2), min_snr=1.1)
        for elem in peakind:
            peaks.append(Qdata[elem])
        Input.Peaks.append(peaks)
        Qlen = len(Qdata)      
            
    print "================================================="
    print "[num peaks : num images]:"
    flog.write( "=================================================\n" )
    flog.write( "[num peaks : num images]:\n" )
    Histo = Input.histo()
    for key in sorted(Histo.iterkeys()):
        print '[', key, ':', Histo[key], ']'
        flog.write('['+str(key)+':'+str(Histo[key])+']'+'\n')

    print "================================================="   
    print "Constructing graph Laplacian for spectral clustering..."
    sgm = cp.sgm
    W = GetDW(Input.Peaks, sgm)
    Lgth = len(W)
    #mknn(W, 40)
    mknn(W, int(sqrt(Lgth)))
    D = np.zeros((Lgth, Lgth))
    L = np.zeros((Lgth, Lgth))
    for i, row in enumerate(W):
        tmp = np.sum(row)
        if tmp != 0:
            D[i][i] = tmp
        else:
            D[i][i] = 1
        for j, col in enumerate(row):
            if i != j:
                L[i][j] = -W[i][j]
            else:
                L[i][i] = D[i][i] - W[i][i]

    # Solve generalized eigenvalue problem for Laplacian.      
    Vals, Vecs = linalg.eigh(L, D)
    Vecs = Vecs.transpose()
    
    print "First few Eigenvalues: "
    flog.write("First few Eigenvalues: \n")
    for i in range(0, 20, 1):
        if i >= len(Vals):
            break
        print Vals[i]
        flog.write(str(Vals[i])+"\n")

    print "================================================="
    print "Number of clusters identified: "
    flog.write("=================================================\n")
    flog.write("Number of clusters identified: \n")
    EigZero = 0
    for elem in Vals:
        if elem < cp.sing_val_threshold:
            EigZero += 1
        else:
            break
    print EigZero

    X = np.zeros((EigZero, len(Vals)))
    for i in range(0, EigZero):
        X[i] = Vecs[i]

    X = np.transpose(X)
        
    Input.Label = KMeans(n_clusters=EigZero,random_state=0).fit_predict(X)

    LabelHisto = Input.SortByLabel()

    print "[cluster : first img - last img]:"
    flog.write( "[cluster : first img - last img]:\n")
    for elem in LabelHisto:
        print "[ ", elem, " : ", LabelHisto[elem][1], "-", LabelHisto[elem][1] + LabelHisto[elem][0] - 1," ]"
        flog.write("[ "+str(elem)+" : "+str(LabelHisto[elem][1])+"-"+str(LabelHisto[elem][1] + LabelHisto[elem][0] - 1)+" ]\n")

    cl_sums = SumClusters( mypath, Input, LabelHisto, A, B )
    cl_avs = SumClusters( mypath, Input, LabelHisto, A, B, average=True )

    write_cluster_names_to_file( mypath, Input, LabelHisto, outpath=cp.outputpath, outtag=cp.outtag)
    plot_cluster_sums( LabelHisto, cl_sums, cp.ncluster_thresh, display=False, outpath=cp.outputpath, outtag=cp.outtag, q=np.array(Qdata) )
    plot_cluster_sums( LabelHisto, cl_avs, cp.ncluster_thresh, display=False, outpath=cp.outputpath, outtag=cp.outtag, q=np.array(Qdata),\
                           tail="_cluster_avs_plot.png", offset=1000 )

    print "Images are sorted by Cluster label and by number of Peaks within the Cluster."
    flog.write("Images are sorted by Cluster label and by number of Peaks within the Cluster.\n")
    while True:
        for elem in LabelHisto:
            print "[ ", elem, " : ", LabelHisto[elem][1], "-", LabelHisto[elem][1] + LabelHisto[elem][0] - 1," ]", LabelHisto[elem][0]
            flog.write("[ "+str(elem)+" : "+str(LabelHisto[elem][1])+"-"+str(LabelHisto[elem][1] + LabelHisto[elem][0] - 1)+" ] "+str(LabelHisto[elem][0])+"\n")
        print "Type q to exit the program. Please provide:"
#        Img = raw_input("Image number: ")
    


        Img = raw_input("Cluster number: ")
        if Img == 'q':
            break
        try:
#            print "Displaying file:", mypath + Input.names[int(Img)]
#            Plot(mypath + Input.names[int(Img)], A, B)
            plt.plot( cl_sums[int(Img)])
            plt.draw()
            plt.show()
        except TypeError:
            print "Incorrect value, please try again."

    flog.close()

if __name__ == "__main__":

    cp = cluster_params()

#    cp.path = "../Clustering/raw_dat(copy)"
    cp.path = "/Volumes/DataStore1/AS_SAXS_feb18_data/0p6_saxs/raw_dat/"
#    cp.path = "/Volumes/DataStore1/AS_SAXS_feb18_data/0p6_saxs/raw_postexpt"
#    cp.path = "/Volumes/DataStore1/AS_SAXS_feb18_data/LCP2/raw_dat/"
#    cp.path = "/Volumes/DataStore1/AS_SAXS_feb18_data/LCP2/raw_postexpt/"
    cp.tag = "run5_plate19" #"scan19_MO_mQ" #"run9_plate22" # "scan27_MO_DL2_buffer_2" #   "run15_plate18" #
#"run7_plate21" #"scan19_MO_mQ" # "run9_plate22" # "scan25_MO_mQ_lys2min" #"scan18_syringe2_MO_milliQ" #"scan27_MO_DL2_buffer_2" # #"run10_plate03" #
    
    cp.outputpath = "/Users/e38496/Work/Research/Experiments/AS/SAXS/SAXS_feb18/Results/clustering/"+cp.tag+"/"
    cp.outtag = cp.tag+"_svt2_sgm5_all_PLOT"
    cp.outputpath = "/Users/e38496/Work/Research/Experiments/AS/SAXS/SAXS_feb18/Results/clustering/"+cp.outtag+"/"
    if not os.path.exists( cp.outputpath):
        os.makedirs( cp.outputpath) 

    cp.A, cp.B = 0.06, 0.5
    cp.sgm = abs(cp.B-cp.A)/5
#   cp.sgm = 20 * abs(B-A)/float(Qlen)
    
    cp.sing_val_threshold = 1e-2

    cp.nmin = 0
    cp.nmax = 3000000000
    cp.ncluster_thresh = 0

    cp.write_cluster_params_to_file()
    main( cp )
