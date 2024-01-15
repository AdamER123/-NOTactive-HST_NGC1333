
# coding: utf-8

# In[2]:

'''
Last Updated at the end of August, 2019
Contact: Adam Rubinstein (arubinst@ur.rochester.edu)
For: Dan Watson's Group

The purpose of this code is to reproject/regrid an image. 
That means changing the WCS of one fits to that of another
and also changing the pixel scale or resolution to that of another (likely lower-res)

Primarily based on the reproject module: https://github.com/astropy/reproject and https://reproject.readthedocs.io/en/stable/
You can get it through anaconda by: https://anaconda.org/astropy/reproject
'''

'''
The general structure of this code is that the first half runs through a bulk list of Spitzer and Hubble images
and reprojects them onto a low-res Spitzer image. Before running the code, you will need to see im_conv first 
and set up a similar method of reading in wavelengths.

This can be modified using an individual case in the latter half of this code (the latter half has various 
examples I was testing).
'''

'''
The detailed steps of how this works as Dan explained...
Convolution theorem: we have an image and a PSF (e.g. a Gaussian or a sinc func)
2D Fourier transf them, multiply them, 2D FT the result back. 
Can use FFT routines built into libraries (EX: astropy). 
May have to alias power down to an order that will reappear. 

The PSF can either be from lambda / D b/c we’re not diffraction limited. Can also measure it from background stars 

Often people first figure out the spatial resolution of each image. 
Say you have 0.2” and 9” for two images. 
Take the low res image, convolve with the high res PSF. 
Take the high res image, convolve with the low res PSF. 

Then have to regrid the images to conserve flux and energy on the same grid. 
Usually the final grid you want to have is usually that of the final image. 
For each point on the sky, you normally need to generate a set of artificial pixels with the orientation you want (from the WCS)
and reorient all your pixels to the same spacing and pixel width. 
'''

#check what version of python you're using - I'm using 3.7.3
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
path = '../n1333_photometry_ds9.bck.dir/**'
im_names_spitz = im_name_finder(path, 'fit')
im_names_spitz = [i.replace('\\', '/') for i in im_names_spitz]

path = '../HH12_Hubble/*'
im_names_hub = im_name_finder(path, 'fit')
im_names_hub = [i.replace('\\', '/') for i in im_names_hub]


# In[28]:

#Minimal loop through all the images, including a try/except in case an image is faulty for whatever reason
#IMPORTANT: A more detailed example between two single images is at the end of this code with many more comments

#First, we need to setup the image we're projecting ONTO
#This will be the same no matter the loop, so only need to do this once
low_res = [x for x in im_names_spitz if 'n1333_lh_3_SiII_' in x][0]  #finding the lowest res image - LH and long lambda, so [SiII]
hdu1 = fits.open(low_res)[0]
#low_res = im_names_hub[ind_proj_to_hub]   #Could also let this be a hubble image...if so, see the hubble loop below for how to setup grabbing the data, header, pixel conversion from the hdu
#hdu1 = fits.open(im_names_hub[ind_proj_to_hub])
hdu1_data = hdu1.data
hdu1_header = hdu1.header
hdu1_pix = hdu1.header['CDELT2'] #the pixel size in degrees, CDELT is the keyword for Spitzer images
hdu1_pix_torad = hdu1_pix * np.pi / 180.
#hdu1_pix = hdu1[0].header['D001SCAL'] #same as above line, but D001SCAL is the keyword for Hubble images

#reprojecting images in loop
for name in im_names_spitz:
    hdu2 = fits.open(name)[0]  #Image we're reprojecting...indexed by 0 since only header
    try:
        #reading in data
        hdu2_data = hdu2.data
        hdu2_header = hdu2.header        
        
        #convolving images
        D = 85 #that of Spitzer, in cm
        D *= 1e4 #converting to microns since x cm / 1 cm * 1E4 microns gets microns, the unit of our wavelength file
        
        hdu_conv = im_conv(low_res, D, hdu1_pix, hdu2_data)
        
        
        #converting the convolved image to correct units and saving it so we can reproject it
        #you'll need to set the WCS to be that of the header you're basing this off of...ie the header 
        w = WCS(hdu2_header)
        wcs_header = w.to_header()
        file_start = 'Convolved_Images/conv_' 
        conv_path = fits_saver(hdu_conv, wcs_header, name, file_start)
        
        
        #reprojection of one hdu using the header (coords and pixels) of another
        #The first input is the path to the file we're reprojecting. The second input is the header of the image we're projecting ONTO
        #para is False for large images (like these hubble ones)
        #output is array (a 2D array of data) and footprint (the footprint from the analysis)
        para = True
        array, footprint = reproject_exact((hdu_conv, hdu2_header), hdu1_header, parallel=para)
        
        
        #saving a new fits file from the reprojected image
        #first, grabbing the WCS coords of the appropriate image to be set as the header of the new image
        # remember to have the right header with the wcs below and that it matches the one we're projecting ONTO
        w = WCS(hdu1_header)
        wcs_header = w.to_header()
        save_path = 'Regridded/regrid_'  #See fits_saver's "save_path" description for explanation
        fits_saver(array, wcs_header, name, save_path)  #saving the reprojected image
        
    except: print('fail', name)


for name in im_names_hub:
    hdu2 = fits.open(name)  #Image we're reprojecting...no index yet since some useful data split between headers
    try: 
        #reading in data
        hdu2_data = hdu2[1].data
        hdu2_header = hdu2[1].header
        hdu2_pix = hdu2[0].header['D001SCAL'] #same as above line, but D001SCAL is the keyword for Hubble images
        hdu2_fnu = hdu2[0].header['PHOTFNU']
        
        #convolving images
        D = 2.4 #that of Spitzer, in m
        D *= 1e6 #converting to microns since x m / 1 m * 1E6 microns gets microns, the unit of our wavelength file
        
        hdu_conv = im_conv(low_res, D, hdu1_pix, hdu2_data)
        
        #converting the convolved image to correct units and saving it so we can reproject it
        #conversion needed for hubble case since units are not in terms of surface brightness
        hdu_conv_scaled = hdu_conv
        hdu_conv_scaled = hdu_conv_scaled * hdu2_fnu / 1e6 #converting to MJy
        hdu_conv_scaled = hdu_conv_scaled / hdu2_pix**2. * 4.25e10 #dividing out sr, D001SCAL is key for pixel size in arcsec
        
        #you'll need to set the WCS to be that of the header you're basing this off of...ie the header 
        w = WCS(hdu2_header) 
        wcs_header = w.to_header()
        file_start = 'Convolved_Images/conv_'  
        conv_path = fits_saver(hdu_conv_scaled, wcs_header, name, file_start)

        
        #reprojection of one hdu using the header (coords and pixels) of another
        #The first input is the path to the file we're reprojecting. The second input is the header of the image we're projecting ONTO
        #para is False for large images (like these hubble ones)
        #output is array (a 2D array of data) and footprint (the footprint from the analysis)
        para = False
        array, footprint = reproject_exact(conv_path, hdu1_header, parallel=para)
        
        
        #saving a new fits file from the reprojected image
        #first, grabbing the WCS coords of the appropriate image to be set as the header of the new image
        # remember to have the right header with the wcs below and that it matches the one we're projecting ONTO
        w = WCS(hdu1_header)
        wcs_header = w.to_header()
        save_path = 'Regridded/regrid_'  #See fits_saver's "save_path" description for explanation
        fits_saver(array, wcs_header, name, save_path)  #saving the reprojected image

    except: print('fail', name)
        
sys.exit()
























# In[ ]:

#alternative method to find all the paths in a directory I was working on if glob(...recursive) doesn't work
# def im_name_finder(path, im_name_end):
#     #Using glob (it's a unix command similar to ls)
#     #* is a wild card, meaning "all", so this grabs all files
#     frames = glob.glob(path)
        
#     #looping through images to grab the image files
#     im_names = []
#     for i in frames:
#         i = i.replace('\\', '/')                                          #replacing a substring b/c glob replaces some /'s with // (maybe a windows thing?)
#         im_names.append(glob.glob(i + '/n1333*')[0].replace('\\', '/'))   #in this case, n1333* means any file beginning with n1333
    
#     return im_names

# #set your path to some directory with images (the images can be in subdirectories)
# path = '../n1333_photometry_ds9.bck.dir/*'

# #set im_name_start...this should be the start of an image file name
# im_name_start = 'n1333'

# im_names = im_name_finder(path, im_name_start)


# In[ ]:

#More detailed example/explanations of the steps for an individual case

#EX: Reading in data from our list of file paths
#first setting up indices to grab file names for the list of image names 
ind_proj = 2            #the index of the image that's being modified, regridded, or reprojected
ind_proj_to = -6        #the index of the new image coordinates we're projecting onto
ind_proj_to_hub = 1     #if wanting to project onto a hubble image (e.g. reproject the 1.2 micron image onto the 1.6 micron one)

#Now opening two example fits files
#For spitzer or hubble images, gave examples below of whether to use index 0 or 1 for the hdu when grabbing data or header info
#Since this header info is in inconsistent locations (two headers),
#probably fast enough to just grab this outside of a function

#Image we're projecting onto
low_res = [x for x in im_names_spitz if 'n1333_lh_3_SiII_' in x][0]
# low_res = im_names_hub[ind_proj_to_hub]
hdu1 = fits.open(low_res)[0]
#hdu1 = fits.open(im_names_hub[ind_proj_to_hub])
hdu1_data = hdu1.data
hdu1_header = hdu1.header
hdu1_pix = hdu1.header['CDELT2'] #the pixel size in arcsecs, CDELT is the keyword for Spitzer images
hdu1_pix_torad = hdu1_pix * np.pi / 180.
#hdu1_pix = hdu1[0].header['D001SCAL'] #same as above line, but D001SCAL is the keyword for Hubble images and is in arcsecs!

#Image we're reprojecting
#hdu2 = fits.open(im_names_spitz[ind_proj])[0]
#hdu2_data = hdu2.data
#hdu2_header = hdu2.header
#hdu2_pix = hdu2.header['D001SCAL']
hdu2 = fits.open(im_names_hub[ind_proj])
hdu2_data = hdu2[1].data
hdu2_header = hdu2[1].header
# hdu1_pix = hdu1[0].header['CDELT2'] #the pixel size in degs
hdu2_pix = hdu2[0].header['D001SCAL'] #same as above line, but D001SCAL is the keyword for Hubble images and is in arcsecs!


#Convolving images
'''
one more function or loop to make is after convolving all images to the same resolution as a stack
(e.g. the longest lambda, SiII image's PSF) 
if we're interested in high precision (e.g. extinction with the two Hubble iron lines)
we should then convolve both images with each other's psf
'''

#Need to set the telescope's D for whichever telescope we're working with
#This is because the convolution convolves a given image with the PSF of the image we're projecting ONTO, which depends on 
#the FWHM...i.e. theta = lambda/D
#The wavelength used to create the PSF is set by the function matching to a file of image wavelengths. The idea is so that
#order won't matter too much, but that can be adjusted if manual input is preferred

#for spitzer cases
D = 85 #that of Spitzer, in cm
D *= 1e4 #converting to microns since x cm / 1 cm * 1E4 microns gets microns, the unit of our wavelength file

#for hubble cases
# D = 2.4 #that of Spitzer, in m
# D *= 1e6 #converting to microns since x m / 1 m * 1E6 microns gets microns, the unit of our wavelength file

hdu_conv = im_conv(low_res, D, hdu1_pix_torad, hdu2_data)


# In[ ]:

#first part of running reproject: rescaling and saving the convolved image
#important! before using reproject, one may have to convert to surface brightness units!
#see note at top of: https://reproject.readthedocs.io/en/stable/celestial.html#spherical-polygon-intersection
#see explanation of surface brightness units in  http://coolwiki.ipac.caltech.edu/index.php/Units#Units_of_Spitzer_Images

#Spitzer is in MJy/sr which is fine BUT Hubble is in units of elecs/sec
#hubble ex:
hdu_conv_scaled = hdu_conv
hdu_conv_scaled = hdu_conv_scaled * hdu2[0].header['PHOTFNU'] / 1e6 #converting to MJy
hdu_conv_scaled = hdu_conv_scaled / hdu2[0].header['D001SCAL']**2. * 4.25e10 #dividing out sr, D001SCAL is key for pixel size in arcsec

#you'll need to set the WCS to be that of the header you're basing this off of...ie the header 
w = WCS(hdu2_header)
wcs_header = w.to_header()
file_start = 'Convolved_Images/conv_' 
conv_path = fits_saver(hdu_conv_scaled, wcs_header, im_names_hub[ind_proj], file_start)


# In[ ]:

#reprojection of one hdu using the header (coords and pixels) of another
#The first input is the path to the file we're reprojecting
#The second input is the header of the image we're projecting ONTO
#Set parallel equal to False if working on very large images (i.e. Hubble imgs)
#output is array (a 2D array of data) and footprint (the footprint from the analysis)

para = True
array, footprint = reproject_exact(conv_path, hdu1_header, parallel=para)


# In[ ]:

#saving a new fits file from the reprojected image
#first, grabbing the WCS coords of the appropriate image to be set as the header of the new image
# remember to have the right header with the wcs below and that it matches the one we're projecting ONTO
w = WCS(hdu1_header)
wcs_header = w.to_header()
save_path = 'Regridded/regrid_'  #See fits_saver's "save_path" description for explanation
fits_saver(array, wcs_header, im_names_hub[ind_proj], save_path)  #saving the reprojected image


# In[ ]:

#how to read in WCS objects
#w2 = WCS(hdu2_header)
#wcs_header2 = w2.to_header()

#w1 = WCS(hdu1.header)
#wcs_header1 = w1.to_header()

#how to replace nan values
#where_are_NaNs = np.isnan(hdu_conv_scaled)
#hdu_conv_scaled[where_are_NaNs] = 0

#print(np.amax(hdu_conv_scaled))
#print(hdu1.data.shape)
#print(hdu_conv_scaled.shape)
#array, footprint = reproject_exact((hdu_conv, hdu2.header), hdu1.header)
#array, footprint = reproject_exact((hdu_conv_scaled, hdu2_header), hdu1[1].header)


# In[ ]:

#example plot 1 to check how things look
#plotting the original two images we're comparing

#plotting hdu1
ax2 = plt.subplot(1,3,1, projection=WCS(hdu1.header))
ax2.imshow(hdu1.data, origin='lower')
ax2.coords.grid(color='white')
# ax1.coords['ra'].set_axislabel('Right Ascension')
# ax1.coords['dec'].set_axislabel('Declination')
# ax1.set_title('2MASS K-band')

#plotting hdu2
ax3 = plt.subplot(1,3,2, projection=WCS(hdu2[1].header))
ax3.imshow(hdu2_data, origin='lower')
ax3.coords.grid(color='white')
# ax2.coords['glon'].set_axislabel('Galactic Longitude')
# ax2.coords['glat'].set_axislabel('Galactic Latitude')
# ax2.coords['glat'].set_axislabel_position('r')
# ax2.coords['glat'].set_ticklabel_position('r')
# ax2.set_title('MSX band E')

#plotting convolved hdu
ax1 = plt.subplot(1,3,3, projection=WCS(hdu2[1].header))
ax1.imshow(hdu_conv, origin='lower')
ax1.coords.grid(color='white')


# In[ ]:

#example plot 2 from the docs: trying to check the look of the image

#plotting hdu1 for comparison
ax2 = plt.subplot(1,3,1, projection=WCS(hdu1.header))
ax2.imshow(hdu1.data, origin='lower')
ax2.coords.grid(color='white')
# ax1.coords['ra'].set_axislabel('Right Ascension')
# ax1.coords['dec'].set_axislabel('Declination')
# ax1.set_title('2MASS K-band')

#plotting the convolved and reprojected image
ax1 = plt.subplot(1,3,2, projection=WCS(hdu1.header)) 
ax1.imshow(array, origin='lower')#, vmin=1, vmax=100
ax1.coords.grid(color='white')
# ax1.coords['ra'].set_axislabel('Right Ascension')
# ax1.coords['dec'].set_axislabel('Declination')
# ax1.set_title('Reprojected MSX band E image')

ax2 = plt.subplot(1,3,3, projection=WCS(hdu1.header)) 
ax2.imshow(footprint, origin='lower')#, vmin=1, vmax=100
ax2.coords.grid(color='white')
# ax1.coords['ra'].set_axislabel('Right Ascension')
# ax1.coords['dec'].set_axislabel('Declination')
# ax2.coords['dec'].set_axislabel_position('r')
# ax2.coords['dec'].set_ticklabel_position('r')
# ax2.set_title('MSX band E image footprint')

