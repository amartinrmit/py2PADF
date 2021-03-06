
import os
import numpy as np
import multiprocessing as mp
import io
import sys
import diffraction

# for debugging
import matplotlib.pyplot as plt

#
#----------------------------------------------------------------
# Calculate the correlation
#
#

class correlation:
    def __init__(self,path="/scratch/amartin1",tag="tag",dpname="dp.dbin",dpname2="None",\
                     nx=128,ny=0,wl=1e-10,pw=1.0,dz=1.0,nth=200,\
                     nthreads=1,npatterns=100, bg_estimate=False,\
                     mask_flag = 0, crop_flag=0,nxcrop=0,nycrop=0,\
                     dp_shift_flag = 0, shiftx = 0, shifty = 0,\
                     maskname="None", rebin=-1, nstart=0,\
                     diff=diffraction.diffraction(), diffcorr=False):
                
                 self.wl = wl
                 self.pixel_width = pw
                 self.detector_z = dz
                 self.path = path
                 self.tag = tag
                 self.h5field = "None"
                 self.cx = nx/2
                 self.cy = ny/2
                 self.dpname = dpname
                 self.dpname2 = dpname2
                 self.nth = nth
                 self.nx = nx
                 self.bg_estimate = bg_estimate
                 self.diffcorrflag = diffcorr
                 
                 self.ny = -1
                 self.pwx = pw
                 self.pwy = pw

                 # other parameters
                 self.nthreads = nthreads
                 self.npatterns = npatterns
                 self.nstart = nstart

                 # flags
                 self.mask_flag = mask_flag
                 self.maskname = maskname
                 self.crop_flag = crop_flag
                 self.nxcrop = nxcrop
                 self.nycrop = nycrop
                 self.dp_shift_flag = dp_shift_flag
                 self.shiftx = shiftx
                 self.shifty = shifty
                 self.rebin = rebin

                 if self.crop_flag == 1:
                     self.nxorig = nx
                     self.nyorig = ny
                     self.nx = int(nxcrop)
                     self.ny = int(nycrop)
                     self.cx = self.nx/2
                     self.cy = self.ny/2


                 self.set_rebin_flag()

                 self.inputType = 'dir'
                 self.flist = []
                 self.datafilename = None
                 self.sym_filter = 0
                 self.correlation_sum_name = None
                 
                 self.xshift_list = []
                 self.yshift_list = []

                 self.inputpath = "None"
                 self.inputtag  = "None"

                 self.dp_flag = 0
                 self.diffraction = diff


    def write_corr_config(self,corrfname="correlation_config.txt"):

        fcorr = open(corrfname,'w')
        fcorr.write("input = "+self.dpname+'\n' )
        if self.bg_estimate == True: fcorr.write("inputb = "+self.dpname2+'\n' )
        if self.bg_estimate == True: fcorr.write("Xflag = 1\n" )
        fcorr.write("outpath = "+self.path+'\n')
        fcorr.write("tag = "+self.tag+'\n')
        fcorr.write("wavelength =  "+str(self.wl)+'\n'  )
        fcorr.write("pixel_width =  "+str(self.pixel_width)+'\n'  )
        fcorr.write("detector_z =  "+str(self.detector_z)+'\n'  )
        fcorr.write("cx =  "+str(self.cx)+'\n'  )
        fcorr.write("cy =  "+str(self.cy)+'\n'  )
        fcorr.write("nth =  "+str(self.nth)+'\n'  )
        fcorr.write("nx =  "+str(self.nx)+'\n' )
        fcorr.write("ny =  "+str(self.ny)+'\n' )
        fcorr.write("pwx = "+str(self.pwx)+'\n' )
        fcorr.write("pwy = "+str(self.pwy)+'\n' )
        fcorr.close()


    def calculate_correlation(self, ftailout="_padfcorr_correlation_sum.dbin"):

        if self.mask_flag == 1:
            maskname = self.maskname ##self.path+self.tag+"_mask_processed.dbin"
            mask = io.read_image( maskname, nx=self.nxorig, ny=self.nyorig ).astype( np.float64 )
            mask *= 1.0 / np.max(mask)
            #print "DEBUG <correlation.py calculate_correlation()> mask.shape", mask.shape
            mask_scb = self.shift_crop_bin( mask, True)

        print "DEBUG calculate correlation npatterns", self.npatterns, len(self.flist), self.bg_estimate

        corrsum_def = False
        corrsqsum_def = False
        tag_save = self.tag
        tag = self.tag
        if self.bg_estimate == True:
            tag += "_bg"
        if self.diffcorrflag == True:
            tag += "_diff"
        

        for i in np.arange(self.npatterns/self.nthreads):
            processes = []
            d = []
            if self.dp_flag == 2:
                self.diffraction.generate_diffraction_patterns( self.nthreads, nstart=i*self.nthreads )
                if self.diffraction.noiseflag == 1:
                    self.flist = self.diffraction.get_file_list( self.diffraction.path,\
                                                                                 self.diffraction.tag, "noise.dbin")
                else:
                    self.flist = self.diffraction.get_file_list( self.diffraction.path,\
                                                                     self.diffraction.tag, "diffraction.dbin")


            for j in np.arange( self.nthreads ):
                m = i*self.nthreads+j+1 + self.nstart
                if self.bg_estimate == False:
                    print "Calculating correlation ", m
                elif self.bg_estimate == True:
                    print "Calculating background cross-correlation ", m

                    

                print "debug correlation.py inputType", self.inputType
                if self.inputType == 'dir':
                   
                    if (self.bg_estimate == True) or (self.diffcorrflag == True):
                        m= np.int( np.random.rand()*self.npatterns )
                        image = io.read_image( self.flist[m-1], nx=self.nxorig, ny=self.nyorig )
                        m2 = np.int( np.random.rand()*self.npatterns )
                        image2 = io.read_image( self.flist[m2-1], nx=self.nxorig, ny=self.nyorig )
                    else:
                        image = io.read_image( self.flist[m-1], nx=self.nxorig, ny=self.nyorig )                        

                elif self.inputType == 'multifile':
                                                        
                    if (self.bg_estimate == True) or (self.diffcorrflag == True):
                        m= np.int( np.random.rand()*self.npatterns )
                        image = io.read_image( self.datafilename, ser_index=m-1 )
                        m2 = np.int( np.random.rand()*self.npatterns )
                        print "DEBUG <correlation.py> m2 :", m2
                        image2 = io.read_image( self.datafilename, ser_index=m2-1 )
                    else:
                        image = io.read_image( self.datafilename, ser_index=m-1 )

                elif self.inputType == 'sim':

                    print "<correlation.py: calculate_correlation(): Simulation diffraction patterns not implemented yet"
                    sys.exit()

                #print "<correlation.py: calculate_correlation()> before image corrections: self.nx image.dtype", self.nx, image.shape, image.dtype 

                # mask diffraction pattern
#                if self.mask_flag == 1:
#                    image *= mask
#                    if self.bg_estimate == True:
#                        image2 *= mask

                # shift diffraction pattern
                if self.dp_shift_flag == 1:
                    if not len(self.xshift_list) == 0:
                        self.shiftx = self.xshift_list[m-1]
                    if not len(self.yshift_list) == 0:
                        self.shifty = self.yshift_list[m-1]

                    print "DEBUG <correlation.py> shiftx shifty", self.shiftx, self.shifty
                    image = self.array_shift(image, self.shiftx, self.shifty)

                    if (self.bg_estimate == True) or (self.diffcorrflag == True):
                        if not len(self.xshift_list) == 0:
                            self.shiftx = self.xshift_list[m2-1]
                        if not len(self.yshift_list) == 0:
                            self.shifty = self.yshift_list[m2-1]

                        image2 = self.array_shift(image2, self.shiftx, self.shifty)

                # crop diffraction pattern
                if self.crop_flag == 1:
                    image = self.crop_image(image, self.nxcrop)
                    if (self.bg_estimate == True) or (self.diffcorrflag == True):
                        image2 = self.crop_image(image2, self.nxcrop)
                    

                # rebin image
                if self.rebin_flag == 1:
                    image = self.rebin_pattern( image, self.rebin )
                    if (self.bg_estimate == True) or (self.diffcorrflag == True):
                        image2 = self.rebin_pattern(image2, self.rebin )

                # mask diffraction pattern
                if self.mask_flag == 1:
                    image *= mask_scb
                    if (self.bg_estimate == True) or (self.diffcorrflag == True):
                        image2 *= mask_scb


                self.nx = image.shape[0]
                #print "<correlation.py: calculate_correlation(): self.nx image.dtype", self.nx, image.shape, image.dtype

                self.dpname = self.path+tag+"diffraction_"+str(j)+".dbin"
                if (self.diffcorrflag == False):
                    io.write_dbin( self.dpname, image )
                else:
                    io.write_dbin( self.dpname, image - image2 )
                print "DEBUG <correlation.py> sum(image) =", np.sum(image)
                if (self.bg_estimate == True):
                    self.dpname2 = self.path+tag+"diffraction2_"+str(j)+".dbin"
                    io.write_dbin( self.dpname2, image2 )
                    print "DEBUG <correlation.py> image diff=", np.sum(image2), np.sum(image - image2)

                self.tag = tag+"_"+str(j)
                self.corrfname = self.path+tag+"corr_config"+str(j)+".txt"
                self.write_corr_config(self.corrfname)
                print "DEBUG <correlation.py> pixel_width ", self.pixel_width
                d.append( self.corrfname )

            for j in np.arange( self.nthreads ):
                #print "thread j :", j, d[j]

                p = mp.Process(target=self.corr_calc, args=(d[j],))
#                if self.bg_estimate == False:
#                    p = mp.Process(target=self.corr_calc, args=(d[j],))
#                elif self.bg_estimate == True:
#                    p = mp.Process(target=self.corr_X_calc, args=(d[j],))
                p.start()
                processes.append(p)

                for p in processes:
                    p.join()

            for j in np.arange( self.nthreads ):
                self.tag = tag+"_"+str(j)
                
                cname = self.path+self.tag+"_correlation.dbin"
                corr = io.read_correlation( cname , 0)
                print "DEBUG <correlation.py> cname ", cname, np.max(corr)
                #print "DEBUG (calculate correlation) corr.size, corr.shape", corr.size, corr.shape, self.nx, self.nth,\
                #    (self.nx/2)*(self.nx/2)*self.nth
                corr = corr.reshape( self.nx/2, self.nx/2, self.nth )

                if corrsum_def == False:
                    corrsum = np.copy(corr)
                    corrsum_def = True
                else:
                    corrsum += corr

                if corrsqsum_def == False:
                    corrsqsum = corr*corr
                    corrsqsum_def = True
                else:
                    corrsqsum += corr*corr


            # remove the original simulated diffraction patterns:
            if self.dp_flag == 2:
                self.diffraction.remove_diffraction_patterns( self.nthreads )


        self.tag = tag_save  # reset the tag
        corrsum *= 1.0/float(self.npatterns)
        print "DEBUG <correlation.py> correlation sum output name :", self.path+self.tag+ftailout, np.max(corrsum)
        io.write_dbin( self.path+self.tag+ftailout, corrsum )
        if self.bg_estimate == False:
            self.correlation_sum_name = self.path+self.tag+ftailout

    def corr_calc( self, qfname ):
        dirpath = os.path.dirname(os.path.realpath(__file__))       
        #print "DEBUG (corr_calc) dirpath:", dirpath
        os.system(dirpath+"/../padfcorr/padfcorr "+qfname ) 

    def corr_X_calc( self, qfname ):
        dirpath = os.path.dirname(os.path.realpath(__file__))
#        os.system("lldb "+dirpath+"/../padfXcorr/padfXcorr "+qfname) 
        os.system(dirpath+"/../padfXcorr/padfXcorr "+qfname) 


    def correlate_by_chunk(self):
        
        # loop over the number of chunks

        #   calculate the diffraction patterns in chunk
        #   or define a list of file names to be included in the chunk

        #   correlate the chunk
        #   add the chunk to sum or output the correlation of the chunk


        csum_name = "_padfcorr_correlation_sum.dbin"
        csig_name = "_padfcorr_correlation_sigma.dbin"


    def radial_profile_to_correlation( self, rad=np.ones(1) ):

        if rad.size==1: rad = np.ones(self.nx/2)

        corr_out = np.zeros( (self.nx/2, self.nx/2, self.nth) )
        for i in np.arange( self.nth ):
            corr_out[:,:,i] = np.outer( rad, rad )
                
        return corr_out

    #
    # this should be sufficient for both hpfilter and the gaussian filter
    #
    def filter( self, filter=np.ones(1),\
                      ftailin="_padfcorr_correlation_sum.dbin",\
                      ftailout="_padfcorr_correlation_sum.dbin"):

        if filter.size == 1 : 
            filter = np.ones(self.nx/2,self.nx/2,self.nth)
        # routine to modify correlation with a hp filter
        cs = io.read_correlation( self.path+self.tag+ftailin )
        cs = cs.reshape( self.nx/2, self.nx/2, self.nth )
        cs *= filter
        io.write_dbin( self.path+self.tag+ftailout, cs )
        # add option not to overwrite??



    # routine to symmeterize the correlation function
    #
    # symmetry noise filter
    #
    def symmetry_filter( self, width=10,\
                             ftailin="_padfcorr_correlation_sum.dbin",\
                             ftailout="_padfcorr_correlation_symmetric.dbin"):

        cs = io.read_correlation( self.path+self.tag+ftailin)
        cs = cs.reshape( self.nx/2, self.nx/2, self.nth )
      
        filter = cs*0.0
        filter[:,:,width:self.nth-width] = 1.0    
        csmod = (cs*filter + np.roll(cs*filter,self.nth/2,2)) \
            / (filter + np.roll(filter,self.nth/2,2))
        io.write_dbin( self.path+self.tag+ftailout, csmod )
        self.correlation_sum_name = self.path+self.tag+ftailout
     
    #
    # Divide by mask function correlation
    #
    def mask_correction( self,\
                             ftailin="_padfcorr_correlation_sum.dbin",\
                             ftailout="_padfcorr_correlation_sum.dbin"):

        maskcorr = io.read_correlation( self.path+self.tag+"_mask_correlation.dbin" )
#        print "DEBUG <correlation.py; mask_correlation()> maskcorr.shape", maskcorr.shape
        maskcorr = maskcorr.reshape( self.nx/2, self.nx/2, self.nth )
        imaskcorr = np.where( maskcorr > 0.05*np.max(maskcorr) )
        cs = io.read_correlation( self.path+self.tag+ftailin )
        cs = cs.reshape( self.nx/2, self.nx/2, self.nth )
        corr_masked = cs * 0.0
#        print "DEBUG <correlation.py; mask_correlation()> cs.shape, maskcorr.shape", cs.shape, maskcorr.shape
        corr_masked[imaskcorr] = cs[imaskcorr] / maskcorr[imaskcorr]
        cs = corr_masked
        io.write_dbin( self.path+self.tag+ftailout, cs )
        self.correlation_sum_name = self.path+self.tag+ftailout
        
    #
    # Subtract background correlation function 
    #
    def subtract_correlation( self, bgname="bg.dbin",\
                             ftailin="_padfcorr_correlation_sum.dbin",\
                             ftailout="_padfcorr_correlation_sum.dbin"):
        
        cs = io.read_correlation( self.path+self.tag+ftailin )
        cs = cs.reshape( self.nx/2, self.nx/2, self.nth )
        print "DEBUG <correlation.py> subtract_correlation() max(cs)", np.max(cs), self.path+self.tag+ftailin

        csbg = io.read_correlation(  bgname )
        csbg = csbg.reshape( self.nx/2, self.nx/2, self.nth )
        print "DEBUG <correlation.py> subtract_correlation() max(csbg)", np.max(csbg), bgname

        cs += -csbg
        io.write_dbin( self.path+self.tag+ftailout, cs )
        self.correlation_sum_name = self.path+self.tag+ftailout
        print "DEBUG <correlation.py> subtract_correlation() max(cs) (corrected)", np.max(cs)

    # shift - a 2D version of numpy's roll
#    def array_shift(self,array,xshift=0,yshift=0):
#        
#        nx, ny = array.shape[0], array.shape[1]
#        nx0 = nx - nx/2
#        ny0 = ny - ny/2
#        tmp = np.zeros( (nx*2, ny*2) ) 
#        tmp[ nx-nx/2:nx+nx/2+1, ny-ny/2:ny+ny/2+1] = array
#        tmp[ nx0:nx0+nx, ny0:ny0+ny] = array
#        
#        tmp = np.roll(tmp,xshift,0)
#        tmp = np.roll(tmp,yshift,1)
#
#        array = tmp[ nx0:nx0+nx, ny0:ny0+ny]
#        
#	array = np.roll(array,xshift,0)
#	array = np.roll(array,yshift,1)
#	return array

    # shift - a 2D version of numpy's roll
    # with sub-pixel shifts
    # interpolation by nearest neighbour interpolation
    def array_shift(self,array,xshiftin=0.0,yshiftin=0.0):
        
        xshift = int(xshiftin)
        yshift = int(yshiftin)
        xrem = xshiftin - xshift
        yrem = yshiftin - yshift

        if xrem >= 0.0:
            mx = 0
            alpha = xrem
        else:
            mx = -1
            alpha = 1 + xrem

        if yrem >= 0.0:
            my = 0
            beta = yrem
        else:
            my = -1
            beta = 1 + yrem

        print xshift, yshift, xrem, yrem
        print mx, my, alpha, beta


        nx, ny = array.shape[0], array.shape[1]
        nx0 = nx - nx/2
        ny0 = ny - ny/2
        tmp = np.zeros( (nx*2, ny*2) ) 
#        tmp[ nx-nx/2:nx+nx/2+1, ny-ny/2:ny+ny/2+1] = array
        tmp[ nx0:nx0+nx, ny0:ny0+ny] = array
        
        # integer shift
        tmp = np.roll(tmp,xshift,0)
        tmp = np.roll(tmp,yshift,1)

        # sub-pixel shift
        tmpsum = (1- alpha) * (1 - beta) * np.roll( np.roll( tmp, mx, 0 ), my, 1)
        tmpsum += alpha * (1- beta) * np.roll( np.roll( tmp, (mx+1), 0 ), my, 1)
        tmpsum += (1 - alpha) * beta * np.roll( np.roll( tmp, mx, 0 ), (my+1), 1)
        tmpsum += alpha * beta * np.roll( np.roll( tmp, (mx+1), 0 ), (my+1), 1)

        tmp = tmpsum

        array = tmp[ nx0:nx0+nx, ny0:ny0+ny]
        
#	array = np.roll(array,xshift,0)
#	array = np.roll(array,yshift,1)
	return array



    def crop_image(self,image,nxcrop,nycrop=-1):
 
 #       print "DEBUG <correlation.py; crop_image> nxcrop", nxcrop
        if nycrop==-1:
            nycrop = nxcrop
        
        dcrop = np.zeros( (nxcrop,nycrop) ) 
        nx, ny = image.shape[0], image.shape[1]
  #      print "DEBUG <correlation.py; crop_image> nx ny", nx, ny, nxcrop, nycrop

        xoffset = 0
        xlow = nx/2-(nxcrop/2)
        if xlow <0 : 
            xlow = 0
            xoffset = (nxcrop-nx)/2

        yoffset = 0
        ylow = ny/2-(nycrop/2)
        if ylow <0 : 
            ylow = 0
            yoffset = (nycrop-ny)/2

        xhigh = nx/2+(nxcrop/2)
        if xhigh > nx: xhigh = nx

        yhigh = ny/2+(nycrop/2)
        if yhigh > ny: yhigh = ny

#        print "DEBUG <correlation.py: crop_image> x/yoffset, x/ylow, x/yhigh", xoffset, yoffset,\
#            xlow, ylow, xhigh, yhigh, dcrop.shape, image.shape
        dcrop[:,:] = image[xoffset+xlow:xoffset+xhigh,yoffset+ylow:yoffset+yhigh]
#        dcrop[xoffset+xlow:xoffset+xhigh,yoffset+ylow:yoffset+yhigh] = image[xlow:xhigh,ylow:yhigh]
#        dcrop[:,:] = image[nx/2-(nxcrop/2):nx/2+(nxcrop/2),ny/2-(nycrop/2):ny/2+(nycrop/2)]
#        plt.imshow( image )
#        plt.figure()
#        plt.imshow( dcrop )
#        plt.draw()
#        plt.show()
        return dcrop

#    def rebin_pattern(self, image, nbin ):
#
#        imagesum = image * 0.0
#        for i in np.arange( nbin)-nbin/2:
#            for j in np.arange( nbin)-nbin/2:
#                imagesum += np.roll(np.roll( image, i, 0), j, 1) 
#        imagesum *= 1.0/float(nbin*nbin)
#        out = imagesum[::nbin,::nbin]
#        return out

    def rebin_pattern(self, imagein, nbin ):

        image = np.copy(imagein)
        
        # nearest neighbour interpolation for even rebin values
        if( (nbin%2)==0):
            imagesum = image * 0.0
            for i in np.arange( nbin)-nbin/2:
                for j in np.arange( nbin)-nbin/2:
                    imagesum += np.roll(np.roll( image, i, 0), j, 1) 
                image *= imagesum/float(nbin*nbin)


        # rebinning
        imagesum = image * 0.0
        for i in np.arange( nbin)-nbin/2:
            for j in np.arange( nbin)-nbin/2:
                imagesum += np.roll(np.roll( image, i, 0), j, 1) 
        imagesum *= 1.0/float(nbin*nbin)
        out = imagesum[::nbin,::nbin]
        return out



    def set_nxcrop_and_nycrop( self, nxcrop, nycrop=-1 ):

        self.nxcrop = nxcrop
        self.nx = nxcrop
        self.cx = self.nx/2

        if nycrop==-1:
            self.nycrop = nxcrop
            self.ny = nxcrop
            self.cy = self.nx/2
        else:
            self.nycrop = nycrop
            self.ny = nycrop
            self.cy = self.ny/2

    def set_rebin_flag( self ):

         if (self.rebin != -1) \
                 and (self.rebin < self.nx) and (self.rebin < self.ny) \
                 and (self.rebin > 0) :
                 
             self.rebin_flag = 1
         else:
             self.rebin_flag = 0

    def binary_mask( self, mask, tol=0.99 ):

        mask *= 1.0/np.max(mask)
        ilow = np.where( mask < tol*np.max(mask) )
        mask[ilow] = 0.0
        return mask

    def shift_crop_bin( self, image, binarize=False ):

        tmp = np.copy(image)

        # shift diffraction pattern
        if self.dp_shift_flag == 1:
            tmp = self.array_shift(image, self.shiftx, self.shifty)
        

         # crop diffraction pattern
        if self.crop_flag == 1:
            tmp = self.crop_image(tmp, self.nxcrop)
                    

         # rebin image
        if self.rebin_flag == 1:
            tmp = self.rebin_pattern( tmp, self.rebin )

        if binarize==True:
            tmp = self.binary_mask( tmp )

        return tmp
        

    def calculate_mask_correlation( self, dontreadfile=False ):

        
        if (dontreadfile==True):
            mask = np.ones( (self.nxorig, self.nyorig) )
        else:
            mask = io.read_image( self.maskname, nx=self.nxorig, ny=self.nyorig ).astype( np.float64 )
            mask *= 1.0/np.max(mask)
            # print "DEBUG <correlation.py; calculate_mask_correlation> mask.shape", mask.shape, mask.dtype
 
        mask = self.shift_crop_bin( mask, True)
       
        # shift diffraction pattern
#        if self.dp_shift_flag == 1:
#            mask = self.array_shift(mask, self.shiftx, self.shifty)
        

#         # crop diffraction pattern
#        if self.crop_flag == 1:
#            mask = self.crop_image(mask, self.nxcrop)
                    

         # rebin image
#        if self.rebin_flag == 1:
#            mask = self.rebin_pattern( mask, self.rebin )
#            mask = self.binary_mask( mask )

        tag = self.tag
        self.tag = tag + "_mask"
        maskname = self.path+self.tag+"_processed.dbin"
        io.write_dbin( maskname, mask )
    #     print "DEBUG <correlation.py; calculate_mask_correlation> processsed mask.shape", mask.shape, mask.dtype
        self.dpname = maskname
        self.corrfname = self.path+self.tag+"_corr_config.txt"
        self.write_corr_config(self.corrfname)
        self.corr_calc( self.corrfname )
        self.tag = tag


    def set_diffraction_parameters( self, diff ):
        self.diffraction = diff

    def set_inputpath( self ):

        if self.inputpath == "None":
            self.inputpath = self.path
        if self.inputtag == "None":
            self.inputtag = self.tag

