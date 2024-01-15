from astropy.io import fits
import sys

#based on https://stackoverflow.com/questions/32643206/how-do-i-open-a-fits-file-from-a-url-in-astropy
#URL that has fits files
url = 'http://www.pas.rochester.edu/~dmw/Documents/Natalie-HST/Build_2/'
start = ['126', '128', '164'] 	#making a list to loop through and grab each file
ending = 'build_shift_2_drz.fits'

for file in start:
	full_url = url + file + ending	#combining parts of file name to have full path
	dat = fits.open(full_url)		#opening fits data at url using astropy fits module

	print dat
	sys.exit()


