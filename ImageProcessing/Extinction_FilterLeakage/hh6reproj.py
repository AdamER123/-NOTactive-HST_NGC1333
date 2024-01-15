from platform import python_version
print(python_version())


# In[3]:

#importing libraries
from astropy.io import fits
from astropy.convolution import convolve, Gaussian2DKernel
from astropy.wcs import WCS
from reproject import reproject_exact  #a package that can be added to astropy using anaconda or pip (see their docs pg)
from reproject import reproject_interp

import glob
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys


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

def im_conv(low_res_name, D, hdu1_pix_torad, hdu2_dat):
    #unfortuantely no good way to find wavelength from header right now. can enter it manually, but I tried to automate it

    #reading in excel file of wavelengths...right now needs to be in same directory as this code
    #first col is a substring of the fits image file name, the second col is the wavelengths in microns
    df = pd.read_excel('imglams.xlsx')
    cols = df.columns
    cols_str = [str(i) for i in df[cols[0]]]
    #some test cases I was using
    #print(low_res_name)
    #print(cols_str)
    #print([i in low_res_name for i in cols_str])
    #print(np.where([i in low_res_name for i in cols_str]))
    #sys.exit()

    #this finds the loc in the excel file where the image substring matches our image name
    #it then finds the wavelength value corresponding to that loc
    lam =  df.loc[np.where([i in low_res_name for i in cols_str])[0][0]].values[1] #lambda in microns

    #finding angular resolution...the FWHM of our Gaussian PSF
    res = 1.2 * lam / D         #resolution in radians
    res = res / hdu1_pix_torad        #so converting to pixels

    #finding PSF and then calculating the convolution of our image and the PSF of the image we're projecting onto
    gauss_kernel = Gaussian2DKernel(res)
    hdu_conv = convolve(hdu2_dat, gauss_kernel)

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
path = '../../../HH12_Hubble/'
im_names_hub_dash = im_name_finder(path+'*', 'fit')
im_names_hub_dash = [i.replace('\\', '/') for i in im_names_hub_dash]
im_names_hub = [path+'126_image.fits', path+'128_image.fits', path+'164_image.fits']
im_names_hh6 = [path+'HH6_126.fits', path+'HH6_128.fits', path+'HH6_164.fits']

#The output from im names hub ['../HH12_Hubble/HH6_126.fits', '../HH12_Hubble/128_image.fits', '../HH12_Hubble/126build_shift_2_drz.fits', '../HH12_Hubble/HH6_164.fits', '../HH12_Hubble/126_image.fits', '../HH12_Hubble/164_image.fits', '../HH12_Hubble/164build_shift_2_drz.fits', '../HH12_Hubble/128build_shift_2_drz.fits', '../HH12_Hubble/HH6_128.fits']
#print(im_names_hub_dash)
#sys.exit()

# In[28]:

#Minimal loop through all the images, including a try/except in case an image is faulty for whatever reason
#IMPORTANT: A more detailed example between two single images is at the end of this code with many more comments

#First, we need to setup the image we're projecting ONTO
#This will be the same no matter the loop, so only need to do this once
#low_res = [x for x in im_names_spitz if 'n1333_lh_3_SiII_' in x][0]  #finding the lowest res image - LH and long lambda, so [SiII]
#hdu1 = fits.open(low_res)[0]
low_res = im_names_hub_dash[2]   #Could also let this be a hubble image...if so, see the hubble loop below for how to setup grabbing the data, header, pixel conversion from the hdu
hdu1 = fits.open(low_res)

#hdu1_pix = hdu1.header['CDELT2'] #the pixel size in degrees, CDELT is the keyword for Spitzer images
#hdu1_pix_torad = hdu1_pix * np.pi / 180.
hdu1_pix = hdu1[0].header['D001SCAL'] #same as above line, but D001SCAL is the keyword for Hubble images
hdu1_pix_torad = hdu1_pix / 206265.
hdu1_fnu = hdu1[0].header['PHOTFNU']

low_res = im_names_hh6[0]
hdu1 = fits.open(low_res)
hdu1_data = hdu1[0].data
hdu1_header = hdu1[0].header



for name, name2 in [(im_names_hh6[2],im_names_hub_dash[6]), (im_names_hh6[1],im_names_hub_dash[7])]:
    hdu2 = fits.open(name)  #Image we're reprojecting...no index yet since some useful data split between headers
    hdu2 = fits.open(name2)
    try:
        #reading in data
        hdu2_pix = hdu2[0].header['D001SCAL'] #same as above line, but D001SCAL is the keyword for Hubble images
        hdu2_pix_torad = hdu2_pix / 206265.
        hdu2_fnu = hdu2[0].header['PHOTFNU']

        hdu2 = fits.open(name)
        hdu2_data = hdu2[0].data
        hdu2_header = hdu2[0].header


        #convolving images
        D = 2.4 #that of Spitzer, in m
        D *= 1e6 #converting to microns since x m / 1 m * 1E6 microns gets microns, the unit of our wavelength file

        hdu1_conv = im_conv(name, D, hdu2_pix_torad, hdu1_data)
        hdu2_conv = im_conv(low_res, D, hdu1_pix_torad, hdu2_data)

        #converting the convolved image to correct units and saving it so we can reproject it
        #conversion needed for hubble case since units are not in terms of surface brightness
        hdu1_conv_scaled = hdu1_conv
        hdu1_conv_scaled = hdu1_conv_scaled * hdu1_fnu / 1e6 #converting to MJy
        hdu1_conv_scaled = hdu1_conv_scaled / hdu1_pix**2. * 4.25e10 #dividing out sr, D001SCAL is key for pixel size in arcsec

        hdu2_conv_scaled = hdu2_conv
        hdu2_conv_scaled = hdu2_conv_scaled * hdu2_fnu / 1e6 #converting to MJy
        hdu2_conv_scaled = hdu2_conv_scaled / hdu2_pix**2. * 4.25e10 #dividing out sr, D001SCAL is key for pixel size in arcsec


        #you'll need to set the WCS to be that of the header you're basing this off of...ie the header
        w = WCS(hdu2_header)
        wcs_header = w.to_header()
        file_start = 'Convolved_Images_Hub/conv_'
        conv2_path = fits_saver(hdu2_conv_scaled, wcs_header, name, file_start)


        w = WCS(hdu1_header)
        wcs_header = w.to_header()
        file_start = 'Convolved_Images_Hub/conv_'
        conv1_path = fits_saver(hdu1_conv_scaled, wcs_header, low_res, file_start)




        #reprojection of one hdu using the header (coords and pixels) of another
        #The first input is the path to the file we're reprojecting. The second input is the header of the image we're projecting ONTO
        #para is False for large images (like these hubble ones)
        #output is array (a 2D array of data) and footprint (the footprint from the analysis)
        para = False
        print('test1')
        array, footprint = reproject_exact(conv2_path, w, shape_out=hdu1_conv_scaled, parallel=para)
        print('test2')
        sys.exit()

        w = WCS(hdu1_header)
        wcs_header = w.to_header()
        file_start = 'Regridded_Hub/regrid_'
        regrid_2_path = fits_saver(array, wcs_header, name, file_start)
        regrid_2_foot_path = fits_saver(footprint, wcs_header, name, file_start+'footprint_')


        print('test3')

        print(array)
        print(hdu1_conv_scaled)

        #now that we have a reprojected hubble image for hdu2 and both are convolved, need to
        #divide one by the other...then can use the same wcs header that we projected onto (hdu1's)!
        #data_ratio = np.divide(array, hdu1_conv_scaled, out=np.zeros_like(array), where=hdu1_conv_scaled!=0) #need to do np.divide to guarantee we get no divide by zero issue...


        #saving a new fits file from the reprojected image
        #first, grabbing the WCS coords of the appropriate image to be set as the header of the new image
        # remember to have the right header with the wcs below and that it matches the one we're projecting ONTO
        #w = WCS(hdu1_header)
        #wcs_header = w.to_header()
        #save_path = 'Regridded_Hub/regrid_'  #See fits_saver's "save_path" description for explanation
        #fits_saver(data_ratio, wcs_header, 'hub_fe_ratio.fits', save_path)  #saving the reprojected image

    except Exception as e: print(e, name)

sys.exit()

                                                                                                                     