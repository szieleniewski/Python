'''Short scripts to handle QFitsView
extracted spectra for use with IRAF
telluric functions.

Written by Simon Zieleniewski

Last updated 11-03-14

'''

import astropy.io.fits as p
import sys
import os

filename = sys.argv[1]

dat, head = p.getdata(filename, header=True)


head.update('WCSDIM', 1)
head.update('CTYPE1', 'LINEAR')
head.update('CUNIT1', 'microns')
head.update('CDELT1', 0.0001)
head.update('CD1_1', 0.0001)
head.update('CRVAL1', 0.63)

print head['WCSDIM']
print head['CTYPE1']
print head['CUNIT1']
print head['CDELT1']
print head['CRVAL1']
print head['CD1_1']


try:
    del head['LTM1_2']
    del head['LTM3_1']
    del head['LTM2_3']
except:
    print 'No weird 3D headers present'
try:
    del head['LTM1_1']
    del head['LTM2_2']
    del head['CD2_2']
    del head['CTYPE2']
    del head['CDELT2']
except:
    print 'No post-telluric correction weird headers present'

os.remove(filename)
p.writeto(filename, dat, header=head)
print 'Created new version of ', filename




