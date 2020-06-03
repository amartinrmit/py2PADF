
import matplotlib
import numpy as np
import scipy as sp
import scipy.ndimage as spn
from scipy.misc import imsave
import matplotlib.pyplot as plt
import glob
import struct
import os

matplotlib.rcParams.update({'font.size': 16})



def write_dbin( fname, data ):

	f = open( fname, "wb")
	fmt='<'+'d'*data.size
	bin = struct.pack(fmt, *data.flatten()[:] )
	f.write( bin )
	f.close()

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


dataqflag = False
#datapath = "/Volumes/DataStore1/AS_SAXS_feb18_data/0p6_SAXS/images/"
#run = "run9_plate22_svt2_sgm5_all"
#run = "run8_plate16_svt2_sgm5_all"
#run = "run2_plate15_svt2_sgm5_10_all"
#run = "plate4_svt2_sgm20_10_all"
run = "run15_plate18_svt2_sgm20_10_all"
#run = "run7_plate21_svt2_sgm5_all"
#run = "run5_plate19_svt2_sgm5_all"
#run = "scan26_MO_mQ_lys60min"
#run = "scan27_MO_DL2_buffer_2_svt2_sgm5_3000"
tag = ""
samplepath = "/Users/e38496/Work/Research/Experiments/AS/SAXS/SAXS_feb18/Results/clustering/"+run+"/"
#sampletag = "run15_plate18_test_9_filelist.txt"
#sampletag = run +"_"+tag+"_"+cluster+"_filelist.txt"
outputtag = '141119'
outpath = "/Users/e38496/Work/Research/Experiments/AS/SAXS/SAXS_feb18/Results/clustering/"+run+"/maunual_recluster/"+outputtag+"/"
if not os.path.exists( outpath):
	os.makedirs( outpath )

groups = []
#fontize = 16


#  scan27_MO_DL2_buffer_2_svt2_sgm5_3000
#groups.append( ['0',  '4',  '6',  '8', \
#			'12',  '14', '15', '16', '17',  \
#			'22', '23', '25', '28', '36'\
#			,'38','39','41','47','49','50','51','52', \
#			'54','55','57'] )
#groups.append( ['1', '2', '3', '7', '10', '11', '13', '18', '19', '21',  '26', '27', '29', '30', '31', '32','33','34','35',\
#			'37', '40','42','43','44','45', '48', '53'] )  # diffuse bump
#groups.append( ['5', '24' ] )  # blank
#groups.append( ['9', '56' ] )  # nothing
#groups.append( ['20' ] )  # weird

# run 9 plate 22
#groups.append( ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10',\
#			'11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21',\
#			'22', '23', '24', '25', '26', '27', '28', '29', '30', '31', '32','33','34','35','36'\
#			,'37','38','39' ] )
#groups.append( ['0', '19', '20'] ) # blank
#groups.append( ['15', '23', '33'] ) # almost blank
#groups.append( ['10'] ) # low q
#groups.append( ['0'] ) # odd

# run 15 plate 18
groups.append( ['0', '1', '9'] )
groups.append( ['2', '6', '7', '8'] ) 
groups.append( ['5', '10'] ) 
groups.append( ['4'] ) 


# run8 plate 16
#groups.append( ['1', '2', '3', '4', '5', '6', '8', '9', '10',\
#			'11', '13', '14', '15', '17', '18', '19', '20', \
#			'24', '25', '26', '27', '28', '29', '30', '31' ] )
#
#groups.append( ['22'] ) # blank
#groups.append( ['7', '12', '21', '23'] ) # almost blank
#groups.append( ['16'] ) # low q
#groups.append( ['0'] ) # odd


# run 7 plate 21
#groups.append( ['0', '2', '3', '4', '5', '6', '7', '8', '9', '10',\
#			'12', '13', '15', '16', '17', '19', '20', '21',\
#			'22', '23', '24', '26', '27', '28', '29', '31', '32','33','34','35','36'\
#			,'37','38','39','40','41' ] ) 
#groups.append( ['30'] ) # blank
#groups.append( ['14', '18', '25'] ) # almost blank
#groups.append( ['11'] ) # odd
#groups.append( ['1'] ) # odd


# run 6 plate 17
#groups.append( ['0', '2', '3', '5', '6', '7', '8', '10',\
#			'11', '12', '13', '15', '16', '18', '19', '20', '21',\
#			'22', '24', '26', '27', '28', '29', '31', '32' ] )
#groups.append( ['4', '30'] ) # blank
#groups.append( ['1', '9', '17', '23'] ) # low q
#groups.append( ['25'] ) # odd
#groups.append( ['14'] ) # almost blank


# run 5 plate 19
#groups.append( ['0', '13', '14'] ) # low res
#groups.append( ['30'] ) # blank
#groups.append( ['9', '21', '22', '25'] ) # almost blank
#groups.append( ['1', '2', '3', '4', '5', '6', '7', '8', '10',\
#			'11', '12', '15', '16', '17', '18', '19', '20', \
#			'23', '24', '26', '27', '28', '29', '31', '32','33','34','35'] )


# run 2 plate 15
#groups.append( ['0', '2', '3', '6', '8', '9', '10'] ) # hexagonal
#groups.append( ['1'] ) # low res scat 1
#groups.append( ['4'] ) # low res scat 2
#groups.append( ['5'] ) # odd structure
#groups.append( ['7'] ) # low res scat 3
#groups.append( ['11'] ) # blank


#plate 4
#groups.append( ['0', '2', '3', '4', '6', '7', '8', '9', '10',\
#			'11', '12', '13', '14', '16', '17', '18', '19', '21',\
#			'22', '23', '24', '26', '27', '28', '29', '31', '32','33','34','35','36'\
#			,'37','38','39','40','41','42','43','44','45','47','48','49','50','51','52'] )#normal
#groups.append( ['20'] )  # v short cell
#groups.append( ['25'] )  # odd smaller cell
#groups.append( ['15', '30', '46'] )# extra second peak
#groups.append( ['1', '5'] ) # shoulder


# reference list
#groups.append( ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10',\
#			'11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21',\
#			'22', '23', '24', '25', '26', '27', '28', '29', '30', '31', '32','33','34','35','36'\
#			,'37','38','39','40','41','42','43','44','45','47','48','49','50','51','52'] )#normal
#			,'53','54','55','56','57','58','59','60','61','62','63','64','65','66','67'] )#normal
#			,'68','69','70','71','72','73','74','75','76','77','78','79','80','81','82'] )#normal
#			,'83','84','85','86','87','88','89','90','91','92','93','94','95','96','97'] )#normal

#
# make the file lists for each group
#
for j, g in enumerate(groups):

	inputstring = ""
	clusterstring = ""
	for k, cluster in enumerate(g):

		inputstring += samplepath+run+"_"+cluster+"_filelist.txt "
		clusterstring += cluster+"_"

	outputstring = outpath+run+"_g"+str(j)+"_clusters_"+clusterstring+"filelist.txt"
	os.system( "cat "+inputstring+" > "+outputstring)

#
#
#
group_intensity = []
group_npatterns = []
group_cnpat = []
group_cint = []
for j, g in enumerate(groups):
	gnpat = 0
	cintensity = []
	cnpat = []
	for k, cluster in enumerate(g):


		sampletag = run+"_"+cluster+"_filelist.txt"
		print sampletag

		i = 0
		f = open( samplepath+sampletag, "r")
		for fname in f:
			i += 1
	#	fname = files[i]
			if fname[-1] == '\n':
				fname = fname[:-1]
				print i, fname
	
			data = np.loadtxt( fname ).astype(np.float)
			print data.shape, np.sum(data)
	
			if i > 1:
				datasum += data[:,1]
			elif i == 1:
				datasum = np.copy(data[:,1])

		if k == 0:
			csum  = np.copy(datasum)
		else:
			csum += datasum
	
		cnpat.append( i )
		cintensity.append( datasum )

		gnpat += i 
	group_intensity.append( csum )
	group_npatterns.append( gnpat )
	group_cnpat.append( cnpat )
	group_cint.append( cintensity )

#	if i>5:
#		break
#    plt.imshow( np.abs(image)**0.1 )
#    plt.draw()
#    plt.show()

offset = 1000 #25

# write out the groups to a log file
f = open( outpath + outputtag+"_recluster_log.txt", 'w' )
f.write( "recluster.py output" )
f.write( "run = "+run+"\n" )
f.write( "tag = "+tag+"\n" )
f.write( "samplepath = "+samplepath+"\n" )
f.write( "outpath = "+outpath+"\n" )
f.write( "outputtag = "+outputtag+"\n" )
f.write( "dataqflag = "+str(dataqflag)+"\n" )
f.write( "offset = "+str(offset)+"\n" )
for i, g in enumerate(groups):
	strng = ""
	for cluster in g:
		strng += " "+cluster
	strng += "\n"
	f.write( "group "+str(i)+" : "+strng )
f.close()

#
# output each individual cluster intensity and comparison
#
for i, cint in enumerate( group_cint ):
	f = open( outpath + outputtag+"_recluster_log.txt", 'a' )
	f.write( "group_"+str(i)+" npatterns"+str(group_npatterns[i])+"\n")
	f.close()


	for j, inten in enumerate(cint): 
		datasum = inten/float(group_cnpat[i][j])
		output = np.transpose(np.array( [data[:,0], datasum] ))
		np.savetxt( outpath+sampletag[:-4]+"_"+outputtag+"_group_"+str(i)+"_cluster_"+str(j)+"_average.txt", output )

		if dataqflag == True:
			qcalib = data[:,0]
		else:
			qcalib = scatterbrain_q_calibration( data.shape[0] )
		
		plt.figure()
		p, =  plt.plot( qcalib[:], datasum )
		plt.ylabel( "Intensity (counts)" )
		plt.xlabel( r'q (nm$^{-1}$' )
		plt.savefig( outpath+sampletag[:-4]+"_"+outputtag+"_group_"+str(i)+"_cluster_"+str(j)+"_average_plot.png" )
		plt.close()

		f = open( outpath + outputtag+"_recluster_log.txt", 'a' )
		f.write( "group_"+str(i)+"_cluster_"+str(j)+" npatterns"+str(group_cnpat[i][j])+"\n")
		f.close()

# comparison plot
for i, cint in enumerate( group_cint ):
	plt.figure()
	plist  =  []
	llist = []
	for j, inten in enumerate(cint): 
		datasum = inten/float(group_cnpat[i][j])
		output = np.transpose(np.array( [data[:,0], datasum] ))

		if dataqflag == True:
			qcalib = data[:,0]
		else:
			qcalib = scatterbrain_q_calibration( data.shape[0] )
		
		p, =  plt.plot( qcalib[:], datasum + j*offset )
		plist.append( p )
		llist.append( groups[i][j] )

	plt.legend( plist, llist )
	plt.ylabel( "Intensity (counts)" )
	plt.xlabel( r'q (nm$^{-1}$' )
	plt.savefig( outpath+sampletag[:-4]+"_"+outputtag+"_group_"+str(i)+"_cluster_comparison_plot.png" )
	plt.close()

#
# output each group intensity and comparison
#
for i, inten in enumerate( group_intensity ):
	datasum = inten/float(group_npatterns[i])
	output = np.transpose(np.array( [data[:,0], datasum] ))
	np.savetxt( outpath+sampletag[:-4]+"_"+outputtag+"_group_"+str(i)+"_average.txt", output )

	plt.figure()
	p, =  plt.plot( qcalib[:], datasum + i*offset )
	plt.ylabel( "Intensity (counts)" )
	plt.xlabel( r'q (nm$^{-1}$' )
	plt.savefig( outpath+sampletag[:-4]+"_"+outputtag+"_group_"+str(i)+"_average_plot.png" )
	plt.close()

# group comparison plot
plt.figure()
plist  =  []
llist = []
for i, inten in enumerate( group_intensity ):
	datasum = inten/float(group_npatterns[i])
	output = np.transpose(np.array( [data[:,0], datasum] ))
	np.savetxt( outpath+sampletag[:-4]+"_"+outputtag+"_group_"+str(i)+"_average.txt", output )

	p, =  plt.plot( qcalib[:], datasum + i*offset )
	plist.append( p )
	llist.append( i )

plt.legend( plist, llist )
plt.ylabel( "Intensity (counts)" )
plt.xlabel( r'q (nm$^{-1})$' )
plt.savefig( outpath+sampletag[:-4]+"_"+outputtag+"_group_comparison_plot.png" )



# write out the intensities for each group (txt and plot )
# write out the number of patterns to a text file ( gnpat and cnpat )
# write out each individual cluster intensity to plot and txt file

# a comparison plot of the group intensities and the cluster intensities in each group.
# reclustering needs a particular output tag and probably its own directory





plt.plot( datasum )
plt.draw()
plt.show()
