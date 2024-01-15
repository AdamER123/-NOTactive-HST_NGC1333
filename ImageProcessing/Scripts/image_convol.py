
# coding: utf-8

# In[1]:

'''
The purpose of this code is to reproject/regrid an image. 
That means changing the WCS of one fits to that of another
and also changing the pixel scale or resolution to that of another (likely lower-res)

Primarily based on the reproject module: https://github.com/astropy/reproject and https://reproject.readthedocs.io/en/stable/
You can get it through anaconda by: https://anaconda.org/astropy/reproject
'''

'''
The general steps of how this works...
Convolution theorem...we have an image and a PSF (e.g. a Gaussian or a sinc func)… 
2D Fourier transf them, multiply them, 2D FT the result back. 
Need to use FFT routines built into libraries (EX: astropy). 
Just have to alias power down to an order that will reappear. 

The PSF can either be from 1.22 lambda / D b/c we’re not diffraction limited … you can also measure it from background stars 

Often people first figure out the spatial resolution of each image. 
Say you have 0.2” and 9” for two images. 
Take the low res image, convolve with the high res PSF. 
Take the high res image, convolve with the low res PSF. 

Then have to regrid the images to conserve flux and energy on the same grid. Usually the final grid you want it to be is usually that of the final image. For each point on the sky, you normally need to generate a set of artificial pixels with the orientation you want (from the WCS) and reorient all your pixels to the same spacing and pixel width. 
'''

#check what version of python you're using - I'm using 3.7.3
from platform import python_version
print(python_version())


# In[2]:

#importing libraries
from astropy.io import fits
from astropy.convolution import convolve, Gaussian2DKernel
from astropy.wcs import WCS
from reproject import reproject_exact  #a package that can be added to astropy using anaconda
from reproject import reproject_interp

import glob
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys


# In[3]:

# #finding the path to every fits images in a directory
def im_name_finder(path, file_type):
    #Using glob (it's a unix command similar to ls)
    #WARNING: using recursive=True...depending how many images you use this could be very slow, it's recommended not to have too many subfolders
    #if needed, some example code is commented below that could help make an alternative
    all_names = glob.glob(path, recursive=True)
    
    #IMPORTANT: Using "fit" here because it is inclusive of both fits and FIT...some files end in "FIT" and need to be included
    #using s.lower() include uppercase names
    im_names = [s for s in all_names if 'fit' in s.lower()]
    
    return im_names

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


# In[4]:

'''now need to convolve my image with a PSF
an approx PSF can be found by assuming a 2D Gaussian func with a width (a FWHM) of the diffrac limit
that is the st dev of the Gaussian is abt the st dev is abt lambda/D
a list of PSFs are found on https://docs.astropy.org/en/stable/convolution/kernels.html

Notes:
what we're using for the gaussian width is the FWHM not the radius of the first ring of the diffraction pattern, 
so it's 1.2 not 1.22

D is 85 cm for spitzer
D is 2.4 m for hubble
'''

def im_conv(im_names, ind_proj_to, cols_str, D, hdu1_pixtoarc, hdu2_dat):
    #unfortuantely no good way to find wavelength from header
    #this finds the loc in the dataframe checking every string 
    lam =  df.loc[np.where([i in im_names[ind_proj_to] for i in cols_str])[0][0]].values[1] # in microns

    #finding res...the FWHM of our Gaussian PSF
    res = 1.2 * lam / D #resolution in radians
    res *= 206265 #to convert to arcsec
    res /= hdu1_pixtoarc*3600. #CDELT is how many degs are in a pixel, so converting to pixels, the 3600 is to convert to arcsecs
    
    #finding PSF and then calculating the convolution of our image and the PSF of the image we're projecting onto
    gauss_kernel = Gaussian2DKernel(res) 
    hdu_conv = convolve(hdu2_dat, gauss_kernel)
    
    return hdu_conv


# In[5]:

#set your path to some directory with images (the images can be in subdirectories)
# recommended to use * or **, which is a wild card, meaning "all", so this grabs all files under a certain directory

#using ** will grab all files even in subdirectories...WARNING this will take longer
path = '../n1333_photometry_ds9.bck.dir/**'
im_names_spitz = im_name_finder(path, 'fit')
im_names_spitz = [i.replace('\\', '/') for i in im_names_spitz]

path = '../HH12_Hubble/*'
im_names_hub = im_name_finder(path, 'fit')
im_names_hub = [i.replace('\\', '/') for i in im_names_hub]


#EX
#Reading in data

#to loop would just start with a for loop lke for im in im_names: ...
ind_proj = 2 #the index of the image that's being modified, regridded, reprojected
ind_proj_to = -6 #the index of the new image coordinates we're projecting onto

#opening two example fits files, reading in its data from the 0th index, then closing it
hdu1 = fits.open(im_names_spitz[ind_proj_to])[0]
hdu2 = fits.open(im_names_hub[ind_proj])
hdu2_data = hdu2[1].data

#Convolving images
'''
one more function or loop to make is after convolving all images to the same resolution as a stack
(e.g. the longest lambda, SiII image's PSF) 
then if we're interested in high precision (e.g. extinction with the two Hubble iron lines)
we should then convolve both images with each other's psf
'''

#reading in pandas file of wavelengths...right now needs to be in same directory as this code
#first col is at least a substring of the image file name, the second col is the wavelengths in microns
df = pd.read_excel('imglams.xlsx')
cols = df.columns
cols_str = [i for i in df[cols[0]]]

#for spitzer cases
# D = 85 #that of Spitzer, in cm
# D *= 1e4 #converting to microns since x cm / 1 cm * 1E4 microns gets microns, the unit of our wavelength file

# hdu_conv = im_conv(D, im_names, ind_proj_to, hdu1.header['CDELT2'])

#for hubble cases
D = 2.4 #that of Spitzer, in m
D *= 1e6 #converting to microns since x m / 1 m * 1E6 microns gets microns, the unit of our wavelength file

hdu_conv = im_conv(im_names_spitz, ind_proj_to, cols_str, D, hdu1.header['CDELT2'], hdu2_data)


# In[6]:

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


# In[7]:

#reprojection of one hdu using the header (coords and pixels) of another

#important! before using reproject, one must convert to surface brightness units!
#see note at top of: https://reproject.readthedocs.io/en/stable/celestial.html#spherical-polygon-intersection
#see explanation of surface brightness units in  http://coolwiki.ipac.caltech.edu/index.php/Units#Units_of_Spitzer_Images
#important to note that Spitzer is MJy/sr which is fine but Hubble is in units of elecs/sec
#hubble ex:
hdu_conv_scaled = hdu2_data
hdu_conv_scaled *= hdu2[0].header['PHOTFNU'] / 1e6 #converting to MJy
hdu_conv_scaled /= hdu2[0].header['D001SCAL']**2. * 4.25e10 #dividing out sr, D001SCAL is key for pixel size in arcsec

#reprojecting
#output is array (a 2D array of data) and footprint (the footprint from the analysis)
#you'll need to set the WCS to be that of the header you're basing this off of
array, footprint = reproject_exact((hdu_conv_scaled, hdu2[1].header), hdu1.header)
# array, footprint = reproject_interp((hdu_conv_scaled, hdu2[1].header), hdu1.header)


# In[8]:

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


# In[9]:

#setting up a new fits file to be saved and viewed in DS9
def fits_saver(array, wcs_header, name, subtrac=False):
    '''
    array is the 2d array of data from reprojecting one image onto another
    wcs_header is a header containing the wcs coords of the image that we projected onto
    name is the path to some image you're using. It will get string split at the / character, and the func only takes the last element of that splitting
    subtrac is if you want to look at the footprint. The other inputs don't need to change, but the name does
    '''
    
    #creating a new file and adding the reprojected array of data as well as the WCS that we projected onto
    hdu_new = fits.PrimaryHDU(array, header=wcs_header)
    hdul = fits.HDUList([hdu_new])

    #saving the file
    if subtrac == False:
        new_filename = name.split('/')[-1]  #grabs the file name we were using from before
    else:
        new_filename = 'footprint_'+name.split('/')[-1]  #grabs the file name we were using from before
        
    hdul.writeto('Regridded/regrid_'+new_filename+'.fits', overwrite=True)
    
    return 0


# In[10]:

#saving a new fits file from the reprojected image

#first, grabbing the WCS coords of the appropriate image to be set as the header of the new image
# header = fits.getheader(hdu1)

'''
remember to have the right header with the wcs below and that it matches the one we're projecting ONTO
'''
w = WCS(hdu1.header)
wcs_header = w.to_header()

fits_saver(array, wcs_header, im_names_hub[ind_proj], subtrac=False)


# In[ ]:




# In[ ]:



