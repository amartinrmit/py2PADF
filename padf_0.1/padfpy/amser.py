
import numpy as np
import matplotlib.pyplot as plt
import array
import os
import struct
import sys


class serheader:

      def __init__(self):

          self.placeholder = 0
          self.Head_DataTypeID = 0
          self.Head_TagTypeID = 0
          self.TagTypeID = 0
          self.Head_TotalNumberElements = 0
          self.Head_ValidNumberElements = 0
          self.ValNumElem = 0
          self.offarroff = 0
          self.Head_OffsetArrayOffset = 0
          self.NumDim = 0 
          self.Head_NumberDimensions = 0
          self.DataDim = 2

      def read_ser_header(self, fc, verbose=0):


          # first three short integers to evaluate file   
          tmp = fc[0:6]
        
          head  = struct.unpack("hhh",tmp)
        
          # validity checks
          if head[0] != 18761:
              print 'Error <amser.py> read_ser_header()'
              print 'Error in Head_Byteorder expecting 18761'
              sys.exit()

          if head[1] != 407:
              print 'Error <amser.py> read_ser_header()'
              print 'Error in Head_SeriesID expecting 407'
              sys.exit()

          if head[2] != 528:
              print 'Error <amser.py> read_ser_header()'
              print 'Error in Head_SeriesVersion expecting 528'
              sys.exit()

          # DATATYPEID CONTAINS DIMENSIONALITY OF INDIVIDUAL DATA ELEMENTS
          self.Head_DataTypeID, = struct.unpack("i", fc[6:10])
          self.DataDim = 1
          if self.Head_DataTypeID == 16674:
              self.DataDim = 2

          # print "DEBUG head_datatypeID :", Head_DataTypeID
          if verbose == 1: print "Data elements are "+str(self.DataDim)+"-dimensional arrays"

          # TAGTYPEID CONTAINES INFOMRATION ON TIME/POSITION CONTENT OF TAB PROCEEDING EACH DATA ELEMENT

          self.Head_TagTypeID, = struct.unpack( "i", fc[10:14] )

          self.TagTypeID=0
          if self.Head_TagTypeID == 16722:
              print 'Tag contains time only'
          if self.Head_TagTypeID == 16706:
              self.TagTypeID=1
          if verbose == 1: print 'Tag contains 2D position with time'

          # TOTALNUMBERELEMENTS CONTAINS TOTAL NUMBER OF DATA ELEMENTS (AND TAG ARRAYS)

          self.Head_TotalNumberElements, = struct.unpack( "i", fc[14:18] )
    
          self.TotNumElem = self.Head_TotalNumberElements

          if verbose == 1: print 'Total number of data elements (inividual spectra/images) is ', self.TotNumElem

          # VALIDNUMBERELEMENTS INDICATES HOW MANY ACTUALLY HAVE DATA (USUALLY THE SAME AS TOTNUMELEM)

          self.Head_ValidNumberElements, = struct.unpack( "i", fc[18:22] )

          self.ValNumElem = self.Head_ValidNumberElements

          if verbose == 1: print 'Total number of valid data elements (inividual spectra/images) is ', self.ValNumElem

          # OFFSETARRAYOFFSET INDICATES THE ABSOLUTE NUMBER OF BYTE OFFSET OF THE START OF THE DATA OFFSET ARRAY

          self.Head_OffsetArrayOffset, = struct.unpack( "i", fc[22:26] ) 

          self.offarroff = self.Head_OffsetArrayOffset

          if verbose == 1: print 'Offset position of offset array (containing data locations) is at ', self.offarroff

          # NUMBER DIMENSIONS DESCRIBES THE NUMBER OF DIMENSION ARRAYS IN THE FILE NOT OF THE DATA

          self.Head_NumberDimensions, = struct.unpack( "i", fc[26:30] )  

          self.NumDim = self.Head_NumberDimensions

          if verbose == 1: print 'There are '+str(self.NumDim)+' dimensions in the data structure (NOT THE DATA)'

    # END EVALUATE HEADER

          

class arrayDims:

    def __init__(self):

        self.Dim1_DimensionSize = 0
        self.Dim1_CalibrationOffset = 0.
        self.Dim1_CalibrationDelta = 0.
        self.Dim1_CalibrationElement = 0
        self.Dim1_DescriptionLength = 0
        self.Dim1_Description = ""
        self.Dim1_UnitsLength = 0
        self.Dim1_Units = ""


    def evaluate_array( self, fc, offset, verbose=0 ):

        # DIMENSIONS SIZE INDICATES NUMBER OF ELEMENTS IN THIS DIMENSION

        self.Dim1_DimensionSize, = struct.unpack( "h", fc[offset:offset+2] )
        if verbose==1: print 'Dimension1 Size is ', self.Dim1_DimensionSize
        offset=offset+4

        # CALIBRATIONOFFSET ?
        self.Dim1_CalibrationOffset, = struct.unpack(  "d", fc[offset:offset+8])
        if verbose==1: print 'Calibration offset is ', self.Dim1_CalibrationOffset
        offset=offset+8

        # CALIBRATIONDELTA ?
        self.Dim1_CalibrationDelta, = struct.unpack(  "d", fc[offset:offset+8])
        if verbose==1: print 'Calibration delta is ', self.Dim1_CalibrationDelta
        offset=offset+8

        # CALIBRATIONELEMENT ?
        self.Dim1_CalibrationElement, =struct.unpack(  "i", fc[offset:offset+4])
        if verbose==1: print 'Calibration element is ', self.Dim1_CalibrationElement
        offset=offset+4

        # DESCRIPTIONLENGTH IS THE LENGTH OF THE DESCRIPTION STRING THAT FOLLOWS
        self.Dim1_DescriptionLength, = struct.unpack(  "i", fc[offset:offset+4])
        if verbose==1: print 'Description length is ', self.Dim1_DescriptionLength
        offset=offset+4

        # ;NEED TO START HAVING OFFSETS TO READ VARIABLE SIZE DIMENSION ARRAY
        #offset=58+self.Dim1_DescriptionLength
        #if verbose==1: print offset

        # ;DESCRIPTION IS A STRING THAT DESCRIBES THE DIMENSION
        self.Dim1_Description, = struct.unpack(  str(self.Dim1_DescriptionLength)+"s",\
                                                     fc[offset:offset+self.Dim1_DescriptionLength] )
        if verbose==1: print 'Dimension description is ', self.Dim1_Description
        offset=offset+self.Dim1_DescriptionLength

        # ;UNITSLENGTH DESCRIBES THE LENGTH OF THE UNITS STRING THAT FOLLOWS
        self.Dim1_UnitsLength, = struct.unpack(  "i", fc[offset:offset+4])
        if verbose==1: print 'Units length is ', self.Dim1_UnitsLength
    
        # ;UPDATE OFFSET
        offset=offset+4
        offset2=offset+self.Dim1_UnitsLength

        # ;UNITS IS A STRING THAT DESCRIBES THE UNIT OF THE DIMENSION
        if offset2 != offset:
            self.Dim1_Units, =struct.unpack(  str(self.Dim1_UnitsLength)+"s", fc[offset:offset2] )
        if offset2 == offset:
            self.Dim1_Units='n/a'
        if verbose==1: print 'Unit description is ', self.Dim1_Units

        offset=offset2

        return  offset


class serdata:

  def __init__(self):

      self.dlist = []
      self.ndata = 0


  def read_image_and_append_to_list( self, fc, cp, datoffarr, imgnum, verbose=0 ):


      if cp.Datatype == 2:
          print 'Error <amser.py> read_image_and_append_to_list()'
          print 'Datatype is Unsigned 2-byte integer'
          print 'This datatype is not currently support by read_ser.py - exiting'
          sys.exit()
      elif cp.Datatype == 6:
          if verbose == 1:
                print 'Datatype is Signed 4-byte integer'
     
      
          # read in the image and append it to the data list (self.dlist)
          offset= datoffarr[imgnum]+50
          TotNumElem = cp.ArraysizeX*cp.ArraysizeY
          tmparr = np.zeros( (cp.ArraysizeX*cp.ArraysizeY) , dtype=np.int)
          for i in range(TotNumElem):
              ind = offset + i*4
              tmparr[i], = struct.unpack( "i", fc[ind:ind+4] ) 
 
          self.dlist.append( tmparr.reshape( cp.ArraysizeX, cp.ArraysizeY) ) 
          self.ndata += 1
      else:
          print 'Error <amser.py> read_image_and_append_to_list()'
          print "Datatype ", cp.Datatype, " is not known - error - exiting"
          sys.exit()

#
# Calibration parameters
#
class calparams:

    def __init__(self):

        self.Caloffsetx = 0.
        self.Caldeltax = 0.
        self.Calelementx = 0
        self.Caloffsety = 0.
        self.Caldeltay = 0.
        self.Calelementy = 0
        self.Datatype = 0
        self.ArraysizeX = 0
        self.ArraysizeY = 0

    def read_cal_params( self, fc, offset, verbose=0 ):

	# ;EVALUATE CALIBRATION OFFSET X - ABSOLUTE START OF X
	self.Caloffsetx, =struct.unpack(  "d", fc[offset:offset+8]) #double(data(offset:offset+7),0)
        offset += 8

	#;EVALUATE CALIBRATION DELTA X - X PIXEL SIZE
	self.Caldeltax, =struct.unpack(  "d", fc[offset:offset+8]) #double(data(offset+8:offset+15),0)
        offset += 8

	# ;EVALUATE CALIBRATION ELEMENT X - POSITION OF CALIBRATION OFFSET
	self.Calelementx, =struct.unpack(  "i", fc[offset:offset+4]) #long(data(offset+16:offset+19),0)
        offset += 4

	# ;EVALUATE CALIBRATION OFFSET Y - ABSOLUTE START OFY
	self.Caloffsety, =struct.unpack(  "d", fc[offset:offset+8]) #double(data(offset+20:offset+27),0)
        offset += 8

	# ;EVALUATE CALIBRATION DELTA Y - Y PIXEL SIZE
	self.Caldeltay, =struct.unpack(  "d", fc[offset:offset+8]) #double(data(offset+28:offset+35),0)
        offset += 8

	# ;EVALUATE CALIBRATION ELEMENT Y - POSITION OF CALIBRATION OFFSET
	self.Calelementy, =struct.unpack(  "i", fc[offset:offset+4]) #long(data(offset+36:offset+39),0)
        offset += 4

	# ;DETERMINE DATATYPE OF DATA
	self.Datatype, =struct.unpack(  "h", fc[offset:offset+2]) #fix(data(offset+40:offset+41),0)
        offset += 2
        if verbose==1: print 'Datatype is ', self.Datatype

	# ;DETERIMINE ARRAY SIZES
	self.ArraysizeX, =struct.unpack(  "i", fc[offset:offset+4]) #long(data(offset+42:offset+45),0)
        offset += 4
	self.ArraysizeY, =struct.unpack(  "i", fc[offset:offset+4]) #long(data(offset+46:offset+49),0)
        offset += 4

	if verbose==1: print 'Image is '+str(self.ArraysizeX)+' by '+str(self.ArraysizeY)+' pixels'

        return offset


#
# Function to read data from a particular file
#
# Frame - specifies which image to read. Default -1: read all frames into a list
#
def read_ser( fname, frame=[], verbose=False ):

    if verbose==True:
          v=1
    else:
          v=0

    size = os.path.getsize(fname)
    f = open( fname, 'rb')
    fc = file.read(f)
#    tmp = fc[2:0:-1]

    # READ HEADER
    sh = serheader()
    sh.read_ser_header( fc, v )
  
    # BEGIN EVALUATE DIMENSION ARRAY(S)----------------------------
    # EVALUATE FIRST ARRAY (ALWAYS TRUE)
    if verbose==True:
          print 'Evaluating Dimension Array 1-------------------------------'

    arrayDimsList = []
    arrayDimsList.append( arrayDims())
    offset = 30
    offset = arrayDimsList[0].evaluate_array(  fc, offset, verbose=v )
    if sh.NumDim > 1:
          if verbose==True:
                print 'Evaluating Dimension Array 1-------------------------------'
          arrayDimsList.append( arrayDims())
          offset = arrayDimsList[1].evaluate_array( fc, offset, verbose=v )

    # END EVALUATE DIMENSION ARRAY (S)----------------------------------------
    # BEGIN EVALUATE DATA/TAG OFFSET ARRAYS ----------------------------------
    if verbose==True:
          print 'Evaluating Data/Tag array offsets-------------------------------'

    # ;DATA OFFSET ARRAY IS AN ARRAY CONTIANING THE LOCATIONS OF ALL THE BYTE OFFSETS OF ALL THE
    # ;DATA ELEMENTS

    offset=sh.offarroff
    offset2=sh.offarroff+(sh.TotNumElem*4)

    if verbose==True:
          print offset, offset2
    datoffarr = np.zeros(sh.TotNumElem, dtype=np.int)
    for i in range(sh.TotNumElem):
        ind = offset + i*4
        datoffarr[i], = struct.unpack( "i", fc[ind:ind+4] ) 
    #    print datoffarr

    # ;TAG OFFSET ARRAY IS AN ARRAY CONTAINING THE LOCATIONS OF ALL THE BYTE OFFSET OF THE TAG ARRAYS
    # ;THAT CONTAIN THE TIME/POSITION INFORMATION OF EACH TYPE DATASET

    offset=offset2
    offset2=offset+(sh.TotNumElem*4)

    tagoffarr = np.zeros(sh.TotNumElem, dtype=np.int)
    for i in range(sh.TotNumElem):
        ind = offset + i*4
        tagoffarr[i], = struct.unpack( "i", fc[ind:ind+4] ) 
    # print tagoffarr
    
    # ;END EVALUATE DATA/TAG OFFSET ARRAYS----------------------------------

    if sh.DataDim == 1:
        print "<amser.py> 1D arrays not supported by the python .ser reader - quitting"
        sys.exit()

    if sh.DataDim == 2:

        # ;SET OFFSET TO FIRST DATAOFFSET
	offset=datoffarr[0]

        # READ CALIBRATION PARAMETERS/ARRAY DIMENSIONS FROM FIRST IMAGE
        cp = calparams()
        offset = cp.read_cal_params( fc, offset, v )

        sd = serdata()
        if verbose==True:
              print "DEBUG reached data reading section", frame

        if len(frame) != 0:
            for find in frame:
                sd.read_image_and_append_to_list( fc, cp, datoffarr, find )
        else:
            for i in np.arange(sh.ValNumElem):
                sd.read_image_and_append_to_list( fc, cp, datoffarr, i )

    f.close()
   
    return sd, cp, arrayDimsList, sh




#path = "/scratch/amartin1/Work/Experiments/Monash/Carbon_dec16/C9/fem/"
#fname = "12.25.40 Spectrum image_1.ser"

#data, cp, adl, sh = read_ser( path+fname, [50, 150, 250] )
#print len(data.dlist), data.ndata

#plt.imshow( data.dlist[2]**0.1 )
#plt.figure()
#plt.imshow( data.dlist[1]**0.1 )



#plt.draw()
#plt.show()
