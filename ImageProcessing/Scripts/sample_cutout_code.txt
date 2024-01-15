#NOTE PLS READ: I put functions at the end of the file instead of beginning to give main idea of how to use Cutout2D, which is what I use
#this is a mix of code from me and Amanda (this was adapted by Amanda from my code, then I edited it again)
#the code obviously won't run, but one thing I tend to do is use my "implot" function or print statements to check the outputs
#EX: implot(data, w, False, np.mean(data)), often also mean, median, max, min of data, etc
#for extinction, I usually need a file with images of the 1.26 to 1.64 ratio, so this example makes a cutout of each of the images and divides them

#some key modules to import (mostly just astropy / Cutout2D is what to search)
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy import wcs
#if you want to try this with various images, glob can help
import glob #EX:  filenames = glob.glob('../../Convolved_Images_Hub/*126_image*')
import numpy as np
import sys

# file = '../../Convolved_Images_Hub/conv_126_image.fits'
# file = 'regrid_hub_dash_noleakage_ston_0.8.fits'  #this was the one I used before...
path = 'southernoutflows_amanda/Frame3-selected/'
# file = path+'regrid_lamflam_126_dash_filled_final.fits' # AB: either 126_ or 164_
filelist = [path+'regrid_lamflam_126_dash_filled_final.fits', path+'regrid_lamflam_164_dash_filled_final.fits'] #doing it by hand instead of glob bc 2...glob or a function would be nicer

for file in filelist:
    hdu1 = fits.open(file)        #open fits image
    w = wcs.WCS(hdu1[0].header)   #get wcs coords from header
    data = hdu1[0].data           #get data from hdu
    #hdu1.close()                  #closing just in case? # AB: moved the location of the close

    #plotting to review what data looks like, need to send wcs regardless
    #probably need to make wcs an optional param...
    implot(data, w, False, 0.3)  # AB: change the last number value so plot was easier to see
    
    #initial guesses for cutout coords of scattered(?) light
    #NOTE: you can actually get the exact values in ds9 if you put a box region around your star or HH object (Region -> Shape -> Box)
    #NOTE: Then you have to double click the region or click region info, then switch fk5 and degrees to image 
    #NOTE: I switch to use units of image pixels to locate my regions as opposed to something "more astronomically sensible" like RA/DEC...you can use RA/DEC which is probably more precise, I just found it more difficult to locate precise locations with Cutout2D
    #NOTE: One thing I do convolutedly sometimes is also put the RA/DEC in fk5 and region size in degrees or arcsec, then still switch to image (pix) units...
    #guessing 
    # AB: updated to values for lower region
    
    #For each object or star, I noted down a set of coords and box region size to loop through (Cutout2D is always a box shape I believe)
    #format: [(x, y), (delta_y, delta_x)] 
    #don't ask me why the box width and height are flipped:)
    coords_list = [[(1995-6, 1460+1), (100, 150)],
                  [(2120-6, 1371+1), (100, 100)],
                  [(1200-6, 2000+1), (150, 300)], 
                  [(1480-6, 1950+1), (100,120)],
                  [(2010-6, 1920+1), (180, 200)],
                  [(2480-6, 1960+1), (200, 200)],
                  ] #(Units: Image Pix)

    #plotting cutouts, probably better loops for this
    for i in range(len(coords_list)):
        #unpacking coords
        position = coords_list[i][0] #x, y
        size = coords_list[i][1] #delta_y, delta_x

        #cutting out coordinates using Cutout2D
        cutout = Cutout2D(data, position, size, wcs = w)
        datacut = cutout.data
        wcscut = cutout.wcs

        #plotting *cutout*
        implot(datacut, wcscut, False, 0.3) 
        #plt.savefig('knot'+str(i+1)) #+6 b/c of indexing, would have to adjust that and path

        # AB: add on to try and save the cutouts as a fits file to do the division
        #NOTE: I typically don't do my saves like this, instead I actually have a fits_saver function
        #NOTE: This is because this overwrites the hdu
        #The outline works something like this:
        #data_ratio = np.divide(flux01, flux02, out=np.zeros_like(flux02), where=flux02!=0.) #need to do np.divide to guarantee we get no divide by zero issue...
        ## remember to have the right header with the wcs below and that it matches the one we're projecting ONTO
        #save_path = 'hh_div/regrid_'  #See fits_saver's "save_path" description for explanation
        #fits_saver(data_ratio, hdu_header1, name.split('/')[1].split('_')[0]+'_ston_'+str(perc)+'.fits', save_path)  #saving the reprojected image
        hdu1[0].data = datacut
        hdu1[0].header.update(wcscut.to_header()) #this is the part that often messes up for me...
        hdu1[0].writeto('hh_cutouts/amanda_'+file.split('_')[-4]+'_knot'+str(i+1)+'_shift.fits', overwrite=True)
        #NOTE: For example, here I could do fits_saver(datacut, wcscut.to_header(), filename, save_path) ... would have to look at fits_saver to see how the filename and save_path should be formatted

    hdu1.close()



#some ideas for functions, most used ones I use are definitely implot and fits_saver
#...I meant to work on file_opening and cutout_saving since I do it so often, but I found file and path differences confuse me, needed to get stuff done...ah well

#a plotting code to review what we're analyzing while in python
#data is the input data
#w is the wcs
#wcscond is True or False, True means the axes will be in RA/Dec, and False means axes in pixels
import matplotlib.pyplot as plt
def implot(data, w, wcscond, vmax_p):
    fig = plt.figure()

    if  wcscond == True:
        fig.add_subplot(111, projection=w)
    else:
        fig.add_subplot(111)

    #for christmas turn on GnRd
    #plt.cm.get_cmap('Blues', 6) is another option
    #can also use RdBu...
    #otherwise just use plt.cm.viridis b/c it's not bad
    plt.imshow(data, origin='lower', cmap=plt.cm.viridis, vmin =0, vmax=vmax_p)
    plt.xlabel('RA')
    plt.ylabel('Dec')
    
    

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


#determining cutout coordinates for each knot...so far just do it by eye and in pixel space
#could probably do it in wcs coords, RA/DEC, but I couldn't get that working

def file_open(file):
    hdu1 = fits.open(file)  #import image
    w = wcs.WCS(hdu1[0].header)   #get wcs coords
    data = hdu1[0].data  #getting data from hdu
    hdu1.close()

    return w, data

#only saving the fits files
def cutout_saver(filenames, pos, size, name, save=False):
    for file in filenames:
        hdu1 = fits.open(file)  #import image
        w = wcs.WCS(hdu1[1].header)   #get wcs coords
    #     print(w.array_shape)
    #     w = wcs.utils.wcs_to_celestial_frame(w)

        #cuting out data and wcs
        data = hdu1[1].data    
        cutout = Cutout2D(data, position, size, wcs = w.celestial)
        datacut = cutout.data
        wcscut = cutout.wcs
    #     print(wcscut.is_celestial)

        #updating header with WCS info
        newhead = hdu1[0].header.update(wcscut.to_header())
        hdu1.close()

        #plotting
        implot(datacut, wcscut, False, np.mean(datacut))  
    #     implot(data, new_wcs)     #plot
    #     plt.savefig('datacut.png')
    #     sys.exit()


        #saving full fits file...
        if save == True:
            lamnum = file[file.index('build')-3:file.index('build')]
            fits.writeto('hh_cutouts/'+name+lamnum+".fits", datacut, wcscut.to_header(), overwrite=True)
        #     fits.writeto('HH6_'+lamnum+".fits", datacut, newhead, overwrite=True)

        #     output_hdul = new_wcs.to_fits()
        #     output_hdul[0].data = data
        #     output_hdul.writeto('HH6_'+file[:3]+".fits", overwrite=True)
        #     sys.exit()