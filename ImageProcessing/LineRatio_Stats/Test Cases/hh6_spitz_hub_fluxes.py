#purpose is to generate images in the steps below, first collecting some files
'''
Paths and file needs:
*imglams and spitzer_conversions are excel files, right now I have it so you need to put it as same directory as your code (but could later maybe just give it a path to go to - would be smarter)
*paths to images and data in general
'''
#now the steps
'''
1) read in all the data by noting all the paths to given spitzer and hubble images
2) loop through all the data, read it in, convert units
3) cutout all the data as appropriate
3) create a loop or otherwise hardcode going through all the combinations of images by hand...
4) regrid all the images
5) de-extinct all the images
6) create apertures as appropriate for all the knots
7) perform relevant analyses: e.g. taking ratio and then finding EDFs, summing up the intensities of each knot for noting and saving
'''

from platform import python_version
print(python_version())


# In[3]:

#importing libraries
from astropy.io import fits
from astropy.convolution import convolve, Gaussian2DKernel, Box2DKernel
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
from reproject import reproject_exact  #a package that can be added to astropy using anaconda or pip (see their docs pg)
from reproject import reproject_interp

import glob
import matplotlib 
matplotlib.use('Agg') #invokved b/c just plain matplotlib was insufficient
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys


#switches for the three different parts of this code
switch1 = 'on' #convolving images [needed to put it on for switch 3 at min...need to figure out other solution, eh]
switch1b = 'on' #regridding...
switch2 = 'on' #solving equations
switch3 = 'on' #plotting / graphics of solutions


if switch1 == 'on':

    # In[4]:

    # #finding the path to every fits images in a directory
    def im_name_finder(path, file_type):
        #Using glob (it's a unix command similar to ls)
        #WARNING: using recursive=True...depending how many images you use this could be very slow, it's recommended not to have too many subfolders
        #if needed, some example code is commented towards the latter half of this code that could help make an alternative
        all_names = glob.glob(path, recursive=True)

        #IMPORTANT: Using "fit" here because it is inclusive of both fits and FIT...some files end in "FIT" and need to be included
        #using s.lower() include uppercase names
        im_names = [s for s in all_names if 'fit' in s.lower()]

        return im_names


    # In[5]:

    '''now convolve my image with a PSF of the image we're projecting ONTO
    an approx PSF can be found by assuming a 2D Gaussian func with a width (a FWHM) of the diffrac limit
    that is the st dev of the Gaussian is about the st dev is about = lambda/D
    a list of PSFs are found on https://docs.astropy.org/en/stable/convolution/kernels.html

    Notes:
    FIRST: always must convert hdu1_pixtorad to radians! It's inconsistent otherwise, and lambda/D is generally in radians

    what we're using for the gaussian width is the FWHM, not the radius of the first ring of the diffraction pattern,
    so it's 1.2 not 1.22 times lambda/D

    D is 85 cm for spitzer
    D is 2.4 m for hubble
    '''

    def im_conv(low_res_name, D, hdu_pix_torad, hdu_dat, kern):
        #unfortuantely no good way to find wavelength from header right now. can enter it manually, but I tried to automate it

        #reading in excel file of wavelengths...right now needs to be in same directory as this code
        #first col is a substring of the fits image file name, the second col is the wavelengths in microns
        df = pd.read_excel('imglams.xlsx')
        cols = df.columns
        cols_str = [str(i) for i in df[cols[0]]]
        #some test cases I was using

        #gaussian kernel
        if kern == 'gauss':
            #this finds the loc in the excel file where the image substring matches our image name
            #it then finds the wavelength value corresponding to that loc
            lam =  df.loc[np.where([i in low_res_name for i in cols_str])[0][0]].values[1] #lambda in microns
            
            #finding angular resolution...the FWHM of our Gaussian PSF
            res = 1.2 * lam / D         #resolution in radians
            res = res / hdu_pix_torad        #so converting to pixels

            #finding PSF and then calculating the convolution of our image and the PSF of the image we're projecting onto
            kernel = Gaussian2DKernel(res)
        
        #box kernel
        if kern == 'box':
            kernel = Box2DKernel(16.)

        hdu_conv = convolve(hdu_dat, kernel)
        return hdu_conv



    # In[27]:

    #setting up a new fits file to be saved and viewed in DS9
    #primarily to save the image we reprojected, but can also be used to save the convolved images
    def fits_saver(array, wcs_header, name, save_path):
        '''
        array is a 2d array of data - could be from reprojecting one image onto another or from convolution
        wcs_header is a header containing the wcs coords of the image that we projected onto or of the orig image (if from the convolution)
        name is the path to some image you're using. It will get string split at the / character, and the func only takes the last element of that splitting
        save_path is the folder you want to save to...recommended to also add something to the start of the images names to make it clear what you did to them (e.g. 'Regridded/regrid_')
        '''

        #creating a new file and adding the reprojected array of data as well as the WCS that we projected onto
        hdu_new = fits.PrimaryHDU(array, header=wcs_header)
        hdul = fits.HDUList([hdu_new])

        #saving the file
        new_filename = name.split('/')[-1]  #grabs the file name we were using from before
        hdul.writeto(save_path+new_filename, overwrite=True)

        return (save_path+new_filename)



    # In[21]:

    #EX: grabbing all the fits image paths in a directory, so they can be looped through and their data opened
    #set your path to some directory with images (the images can be in subdirectories)
    #using ** will grab all files even in subdirectories...WARNING this will take longer
    # path = '../../../ngc1333_fits/'
    # im_names_hub_dash = im_name_finder(path+'*', 'fit')
    # im_names_hub_dash = [i.replace('\\', '/') for i in im_names_hub_dash]
    # im_names_hub = [path+'126_image.fits', path+'128_image.fits', path+'164_image.fits', path+'672_image.fits',
    #                 path+'halph_hart_image.fits']

    # In[28]:
    #this time setting up the file names by hand since I've found that easier...
    #order: halpha or .656 mic, 0.672 mic, 1.26, 1.28, 1.64
    # files_units = ['../../../../ngc1333_fits/unregridded/656_image.fits', 
    #                '../../../../ngc1333_fits/unregridded/0301_flt.fits', 
    #                '../../../../ngc1333_fits/unregridded/0501_flt.fits', 
    #                '../../../../ngc1333_fits/126build_shift_2_drz.fits', 
    #                '../../../../ngc1333_fits/128build_shift_2_drz.fits', 
    #                '../../../../ngc1333_fits/164build_shift_2_drz.fits']
    # hdu_list_units = [fits.open(i) for i in files_units]
    # files_data = ['../../../../ngc1333_fits/656_hareproject.fits', 
    #               '../../../../ngc1333_fits/0301_oIreproject2.fits', 
    #               '../../../../ngc1333_fits/672_sIIreproject.fits', 
    #               '../../../../ngc1333_fits/Background_corr/background_corr_126_aligned.fits', 
    #               '../../../../ngc1333_fits/Background_corr/background_corr_128_aligned.fits', 
    #               '../../../../ngc1333_fits/Background_corr/background_corr_164_aligned.fits']
    # hdu_list = [fits.open(i) for i in files_data]
    files_units = ['../../../ngc1333_fits/126build_shift_2_drz.fits', 
                    '../../../ngc1333_fits/128build_shift_2_drz.fits', 
                    '../../../ngc1333_fits/164build_shift_2_drz.fits']
    hdu_list_units = [fits.open(i) for i in files_units]
    files_data = ['../../../ngc1333_fits/Background_corr/background_corr_126_aligned.fits', 
                    '../../../ngc1333_fits/Background_corr/background_corr_128_aligned.fits', 
                    '../../../ngc1333_fits/Background_corr/background_corr_164_aligned.fits']

    path = '../../../n1333_photometry_ds9.bck.dir/**'
    im_names_spitz = im_name_finder(path, 'fit')
    im_names_spitz = [i.replace('\\', '/') for i in im_names_spitz]

    files_units = files_units + im_names_spitz
    files_data = files_data + im_names_spitz

    hdu_list = [fits.open(i) for i in files_data]

    hdu_pix_list = []
    hdu_pixtorad_list = []
    hdu_fnu_list = []
    hdu_flam_list = []
    hdu_bw_list = []
    hdu_data_list = []
    hdu_header_list = []
    # throughput_list = [1., 1., 1., 1., 1., 1.] # [0.242, 1., 0.246, 0.496, 0.521, 0.470] #also has to be done by hand, not in the headers?


    #I'm using count here just to point to specific indices that I've set up...unfortunately some have different headers...
    #the only diff between the if and else cases are the indexing of the hdu's, some need 1 and some need 0
    #I've tried to group it for convience, so the the first two have the same headers, the last 3 have the same headers
    count = 0
    for (hdu_units,hdu_data) in zip(hdu_list_units, hdu_list):
        if count <= 2: 
            #reading in conversions
            hdu_pix_list.append(hdu_units[0].header['D001SCAL'])  #D001SCAL is the keyword for Hubble images
            hdu_pixtorad_list.append(hdu_pix_list[count] / 206265.)
            # hdu_fnu_list.append(hdu_units[0].header['PHOTFNU'])
            hdu_flam_list.append(hdu_units[0].header['PHOTFLAM'])
            hdu_bw_list.append(hdu_units[0].header['PHOTBW'])

            #reading in datafor general use  and header for wcs
            hdu_data_list.append(hdu_data[0].data)
            hdu_header_list.append(hdu_data[0].header)

        else:
            #reading in conversions
            hdu_pix_list.append(hdu_units[0].header['CDELT2'])  #D001SCAL is the keyword for Hubble images
            hdu_pixtorad_list.append(hdu_pix_list[count] * np.pi / 180.)
            # hdu_fnu_list.append(hdu_units[0].header['PHOTFNU'])
            # hdu_flam_list.append(hdu_units[0].header['PHOTFLAM'])
            hdu_bw_list.append(hdu_units[0].header['PHOTBW'])

            #reading in datafor general use  and header for wcs
            hdu_data_list.append(hdu_data[0].data * 1e6 * 1e-23)  #the spiter data is in MJy / sr, so let's convert out the MJy to Flam units
            hdu_header_list.append(hdu_data[0].header)

        count += 1


    #can update later...but basically the sulfur II image header isn't avail...I think this was for hh 7-11 only
    #header info taken from https://www.stsci.edu/hst/instrumentation/wfc3/data-analysis/photometric-calibration/uvis-photometric-calibration/quad-filter-photometry
    # hdu_flam_list[2] = 1.3699e-17
    # hdu_bw_list[2] = 69.98

    #converting spitzer units
    #reading in excel file of bandwidths
    #first col is a substring of the fits image file name, the second col is the wavelengths in microns
    df = pd.read_excel('spitzer_conversions.xlsx')
    cols = df.columns

    for image_name in files_data:
	    cols_str = [str(image_name) for image_name in df[cols[0]]]
	    lam =  df.loc[np.where([image_name in low_res_name for image_name in cols_str])[0][0]].values[1] #lambda in microns


    print('loaded data!')
    sys.exit()










    '''
    need to do mutual convolutions...it doesn't look like a loop would help much
    I really couldn't find a loop, so I just do this mostly by hand...
    the basic idea is to run all combinations of convolutions so each image is convolved with each other's psf

    here I'm using the function im_conv, which requires the data file name to find the wavelength, the D for the resolution to convolve to, the pixel conversion, and the data that we are convolving
    '''

    res_str = 'flam' #used to label what we're saving, usually related to units or whether we're doing a gaussian or box convolution, etc
    # resize = 60. #if trying to adjust size of gaussian convolution
    D = 2.4 #/ resize #that of Hubble, in m
    D *= 1e6 #converting to microns since x m / 1 m * 1E6 microns gets microns, the unit of our wavelength file

    #format: conv(data convolving with, Diameter, pix size convolving with, image data to be convolved, convolution method)
    #do each with its own psf
    hdu_conv_list = [im_conv(i, D, j, k, 'gauss') for i,j,k in zip(files_data, hdu_pixtorad_list, hdu_data_list)]

    # #656 and 631
    hdu656_conv = im_conv(files_data[1], D, hdu_pixtorad_list[1], hdu_conv_list[0], 'gauss')

    # #656 and 672
    hdu656_conv = im_conv(files_data[2], D, hdu_pixtorad_list[2], hdu_conv_list[0], 'gauss')
    hdu672_conv = im_conv(files_data[0], D, hdu_pixtorad_list[0], hdu_conv_list[2], 'gauss')

    # #126 and 128
    hdu126_conv = im_conv(files_data[4], D, hdu_pixtorad_list[4], hdu_conv_list[3], 'gauss')
    hdu128_conv = im_conv(files_data[3], D, hdu_pixtorad_list[3], hdu_conv_list[4], 'gauss')

    # #656 and 126
    hdu656_conv = im_conv(files_data[3], D, hdu_pixtorad_list[3], hdu656_conv, 'gauss')
    hdu126_conv = im_conv(files_data[0], D, hdu_pixtorad_list[0], hdu126_conv, 'gauss')

    # #672 and 128
    hdu672_conv = im_conv(files_data[4], D, hdu_pixtorad_list[4], hdu672_conv, 'gauss')
    hdu128_conv = im_conv(files_data[2], D, hdu_pixtorad_list[2], hdu128_conv, 'gauss')

    #164 convolved with all the others and just one of the others with 164, 656 in this case
    hdu672_conv = im_conv(files_data[5], D, hdu_pixtorad_list[5], hdu672_conv, 'gauss')
    hdu126_conv = im_conv(files_data[5], D, hdu_pixtorad_list[5], hdu126_conv, 'gauss')
    hdu128_conv = im_conv(files_data[5], D, hdu_pixtorad_list[5], hdu128_conv, 'gauss')
    hdu164_conv = im_conv(files_data[0], D, hdu_pixtorad_list[0], hdu_conv_list[5], 'gauss') #164 with 656 independently
    hdu656_conv = im_conv(files_data[5], D, hdu_pixtorad_list[5], hdu656_conv, 'gauss') #656 with 164

    #need to do the OI line with everything...
    hdu631_conv = im_conv(files_data[0], D, hdu_pixtorad_list[0], hdu_conv_list[1], 'gauss')

    #recompiling our list
    hdu_conv_list = [hdu656_conv, hdu631_conv, hdu672_conv, hdu126_conv, hdu128_conv, hdu164_conv]


    print('convolved!')
        
    '''
    onto converting units and regridding...
    '''
    #converting the convolved image to correct units and saving it so we can reproject it
    #conversion needed for hubble case since units are not in terms of surface brightness
    hdu_conv_scaled_list = []

    #for each convolved image, we need to convert from e-/s to flambda units which are in erg/s/cm^2/Angstrom...multiply by bw to get rid of angstrom
    #...also images are divided through by throughput, so multiplying by throughput
    #and lastly dividing by the arcsec^2 to get the image in surface brightness units, to be used in regridding
    for count, i in enumerate(hdu_conv_list):

        #a condition added since the background corrected images have units corrected?
        if  count > 2: #the background subtracted images
            hdu_conv_scaled_list.append(i) #/ hdu_pixtorad_list[count]**2. ) #note if doing regridding should also add in , regrid only takes surface brightness
        elif count == 2: #sulfur image or OI
            hdu_conv_scaled_list.append(i * hdu_flam_list[count] * hdu_bw_list[count]) # / hdu_pixtorad_list[count]**2.)# * hdu_flam_list[count] * hdu_bw_list[count] * throughput_list[count]) #note if doing regridding should also add in  / hdu_pixtorad_list[count]**2., regrid only takes surface brightness
        elif count == 0 or count == 1: #Halpha image, no dividing through by pix^2 b/c it's the basis for regridding
            hdu_conv_scaled_list.append(i * hdu_flam_list[count] * hdu_bw_list[count] )




    #you'll need to set the WCS to be that of the header you're basing this off of...ie the header
    file_start = 'Convolved_Images_Hub/conv_'+res_str+'_'
    conv_path_list = [] #list of paths to the convolved images, can be useful...


    for count, i in enumerate(hdu_header_list):
        #finding wcs for a given image
        w = WCS(i)
        wcs_header = w.to_header()

        #saving each file to some path, conv_path is the path to that file
        conv_path = fits_saver(hdu_conv_scaled_list[count], wcs_header, files_data[count], file_start)
        conv_path_list.append(conv_path)

    print('saved convolved images!')










    #need a wcs standard for regridding and plots, fits files...
    w = WCS(hdu_header_list[0]) #I picked 0 arbitrarily, it shouldn't really matter
    wcs_header = w.to_header()

    # ##############for regridding just in case...
    
    # if switch1b == 'on':
    #     # I put some options below here just in case

    #     #I only put this here in case later someone wants to regrid the images...
    #     #reprojection of one hdu using the header (coords and pixels) of another
    #     #The first input is the path to the file we're reprojecting. The second input is the header of the image we're projecting ONTO
    #     #para is False for large images (like these hubble ones)
    #     #output is array (a 2D array of data) and footprint (the footprint from the analysis)

    #     #an example of reproject by hand here, but you could loop this, makes things more readable!
    #     hdu_compiled_list = []
    #     regrid_path_list = []
    #     regrid_foot_path_list = []

    #     para = False
    #     file_start = 'Regridded_Hub/regrid_'+res_str+'_'

    #     for count, i in enumerate(hdu_header_list):
    #         if count > 0:
    #             array, footprint = reproject_exact(conv_path_list[count], w, shape_out=hdu_conv_scaled_list[0].shape, parallel=para)
    #             hdu_compiled_list.append(array) #* hdu_pixtorad_list[count]**2.)
    #             #here's what you'd have to do if you were saving from reproject...

    #             #saving a new fits file from the reprojected image
    #             #first, grabbing the WCS coords of the appropriate image to be set as the header of the new image
    #             #in this case, the wcs is set to always be the same
    #             #then we multiply by the pix^2 to get the right units...
    #             # remember to have the right header with the wcs below and that it matches the one we're projecting ONTO
    #             regrid_path = fits_saver(array, wcs_header, files_data[count], file_start) #* hdu_pixtorad_list[count]**2., wcs_header, files_data[count], file_start)
    #             regrid_path_list.append(regrid_path)
    #             regrid_foot_path = fits_saver(footprint * hdu_pixtorad_list[count]**2., wcs_header, files_data[count], file_start+'footprint_')
    #             regrid_foot_path_list.append(regrid_foot_path)

    #         else:
    #             hdu_compiled_list.append(hdu_conv_scaled_list[0])


    #     print('saved regridded images!')
    





    #for de-extincting...should do at an earlier stage for the other lines or just not delete them....
	# s_ii_k = 7.279e-1
	# o_i_k = 7.877e-1

	# s_ii_tau = s_ii_k * A_V_arr
	# o_i_tau = o_i_k * A_V_arr

	# s_ii_extinc = np.exp(s_ii_tau)
	# o_i_extinc = np.exp(o_i_tau)

	# file_start_extinc = nonlin_folder + '/nonlin' + res_str + 'deextinc_'
	# conv_extinc_paths = ['fitting_withOI_rms/nonlinflam0p1e-170_631_mic.fits', 'fitting_withOI_rms/nonlinflam0p1e-170_672_mic.fits']
	# hdu_extinc = [fits.open(extinc_path) for extinc_path in conv_extinc_paths]
	# hdu_extinc = [fits.open(extinc_path)[0].data for extinc_path in conv_extinc_paths]
	# hdu_extinc_data = [fits.open(extinc_path)[0].data for extinc_path in conv_extinc_paths]
	# hdu_dextinc = [hdu_extinc_data[0]*o_i_extinc, hdu_extinc_data[1]*s_ii_extinc]

	# fits_saver(hdu_dextinc, wcs_header, 'ex/OI.fits', file_start_extinc)
	# fits_saver(hdu_dextinc[0], wcs_header, 'ex/OI.fits', file_start_extinc)
	# fits_saver(hdu_dextinc[1], wcs_header, 'ex/SII.fits', file_start_extinc)











	    #if you want to divide two of the data files, can do so like
    #data_ratio = np.divide(array, hdu1_conv_scaled, out=np.zeros_like(array), where=hdu1_conv_scaled!=0) #need to do np.divide to guarantee we get no divide by zero issue...

        
