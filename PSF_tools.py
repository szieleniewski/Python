'''Functions to measure PSFs

Author: Simon zieleniewski

Last updated: 15-09-16

'''

import numpy as n
import astropy.io.fits as p
from scipy.special import j1


#------------#
def makearray(dim, spaxel):
    '''Function that creates an array of angular positions from
    the centre in radians per spaxel.

    Inputs:
        - dim: array side length [spaxels]
        - spaxel: spaxel size [mas]

    Outputs:

        - arr: Array of angular distances from centre pixel [radians]

    '''

    #Array
    dim = int(dim)
    #detector array
    arr = n.ones((dim,dim))
    rad_conv = 1.E-3/206265. #factor for converting from mas to radians
    spaxel *= rad_conv #in [rad/spaxel]


    for i in xrange(dim):
        for j in xrange(dim):

            arr[j,i] = n.sqrt(abs( ((dim/2-j))**2 + ((dim/2-i))**2 ))
            if (arr[j,i] == 0.):
                arr[j,i] = 1.E-50

    arr *= spaxel
    return arr
#------------#


#------------#
def r_values(spaxels, spans):
    '''Generates an array of size [spaxels,spaxels] containing the distance
       of each element from the centre point of the array.

       =========================== Input ==============================
       spaxels       - List, number of spaxels in [x, y] directions.
                       x-y grid in the format: datacube[Intensity, y, x].
                       (Should be even)
       span          - List, the width of the edges of the image in arc
                       seconds [x, y].
       centre        - BOOL, if True, image is central, if False then
                       image is offset to centre on a pixel at just off
                       the centre of the spaxel array. (currently removed)
       ============================ Out ===============================
       r             - 2D numpy array. Contains the distance of each
                       element from the centre point of the array.'''

    spans = n.array(spans, dtype=float)

    i=(spaxels[0]/2)-1
    j=(spaxels[1]/2)-1
    incrx=spans[0]/(2*spaxels[0])#Half a spaxel width in arcsec
    incry=spans[1]/(2*spaxels[1])#Half a spaxel width in arcsec
    x=n.linspace(-(spans[0]/2.),i*2*incrx,spaxels[0])
    y=n.linspace(-(spans[1]/2.),i*2*incry,spaxels[1])

    cart=n.meshgrid(x,y)
    r=n.sqrt(cart[0]**2+cart[1]**2)        #Pythagoras
    r /= 206265.                           #Values in radians
    r = n.where(r < 1.E-15, 1.E-30, r)     #Quick and dirty fix to ensure central value is acceptable for Airy function
                                           #1.E-30 seems to give good result for Airy peak

    return r
#------------#


#------------#
def EE(PSF, img_spax, spaxel=(4,4), num_spax=2):
    '''Function to calculate the ensquared energy of a PSF over
    a square of given box size in mas.

    Inputs:

        - PSF: 2D array of PSF
        - img_spax: Tuple (x,y) of input PSF spaxel scale [mas, mas]
        - spaxel: Tuple (x_new,y_new) of box size to measure EE [mas, mas]
        - num_spax: Side length of box of spaxels to sum over. E.g. entering
                  2 means summing over 4 spaxels (2x2 box)

    Outputs:

        - ee: Ensquared Energy fraction

    '''

    loc = n.argwhere(PSF == PSF.max())[0,:]

    img_spax_x = img_spax[0]
    img_spax_y = img_spax[1]

    #No. of spaxels from edge of peak spaxel in x direction
    pix_num_x = (num_spax*spaxel[0]/2.)/float(img_spax_x) - 0.5
    #No. of spaxels from edge of peak spaxel in y direction
    pix_num_y = (num_spax*spaxel[1]/2.)/float(img_spax_y) - 0.5

    box_x = int(pix_num_x)
    box_y = int(pix_num_y)
    rem_x = round(pix_num_x - box_x, 6)
    rem_y = round(pix_num_y - box_y, 6)

    if num_spax*spaxel[0] <= img_spax_x and num_spax*spaxel[1] <= img_spax_y:
        ee = PSF[loc[0],loc[1]]*(spaxel[0]*spaxel[1]/float(img_spax_x*img_spax_y))

    if num_spax*spaxel[0] > img_spax_x and num_spax*spaxel[1] > img_spax_y:
        #Sum of spaxels fully within square
        sum1 = PSF[loc[0]-box_y:loc[0]+box_y+1,loc[1]-box_x:loc[1]+box_x+1].sum()
        #Sum of horizontal edge spaxels
        sum_2_1 = PSF[loc[0]+box_y+1,loc[1]-box_x:loc[1]+box_x+1].sum()*rem_y
        sum_2_2 = PSF[loc[0]-box_y-1,loc[1]-box_x:loc[1]+box_x+1].sum()*rem_y
        sum_2_3 = PSF[loc[0]-box_y:loc[0]+box_y+1,loc[1]+box_x+1].sum()*rem_x
        sum_2_4 = PSF[loc[0]-box_y:loc[0]+box_y+1,loc[1]-box_x-1].sum()*rem_x
        sum2s = sum_2_1 + sum_2_2 + sum_2_3 + sum_2_4
        #Sum of corner spaxels
        sum_3_1 = PSF[loc[0]+box_y+1,loc[1]+box_x+1]*rem_x*rem_y
        sum_3_2 = PSF[loc[0]+box_y+1,loc[1]-box_x-1]*rem_x*rem_y
        sum_3_3 = PSF[loc[0]-box_y-1,loc[1]+box_x+1]*rem_x*rem_y
        sum_3_4 = PSF[loc[0]-box_y-1,loc[1]-box_x-1]*rem_x*rem_y
        sum3s = sum_3_1 + sum_3_2 + sum_3_3 + sum_3_4

        ee = sum1 + sum2s + sum3s

    ee = ee/PSF.sum()

    return ee
#------------#


#------------#
def SR(psf, spaxel, D, eps, wavelength):
    '''Function to calculate the Strehl ratio for a PSF.
    Function creates a normalised diffraction limited PSF
    and divides the two peak intensities by oneanother.

    Inputs:

        - psf: 2D numpy array
        - spaxel: PSF spaxel scale [mas]
        - D: telescope diameter [m]
        - eps: telescope obscuration ratio
        - wavelength: PSF wavelength [um]

    Outputs:

        - sr: Strehl ratio fraction

    '''

    #PSF size
    y, x = psf.shape
    k = 2.*n.pi/float(wavelength*1.E-6)
    #Create diff. limited PSF
    # arr = makearray(x, spaxel)
    arr = r_values([x,y], [x*spaxel/1000., y*spaxel/1000.])
    arr = n.array(arr, dtype=n.float64)
    arr = n.sin(arr)*(D/2.)*k
    array = abs((1/(1-eps**2)**2)*(2*j1(arr)/arr - eps*2*j1(eps*arr)/arr)**2)
    array /= n.sum(array)
    #Calculate Strehl
    sr = n.max(psf)/n.max(array)
    return sr
#------------#


#------------#
def calc_FWHM(obj, samp):
    '''Function that computes the FWHM of a (Gaussian) source

    Inputs:

        - obj: Image (object at centre of image)
        - samp: Spatial sampling [mas]

    Outputs:

        - FWHM: Source FWHMs (x, y) ([mas], [mas])

    '''

    cmaxx = n.max(obj[obj.shape[0]/2,:])
    cmaxy = n.max(obj[:,obj.shape[1]/2])

    cslicex = obj[obj.shape[0]/2,:]
    cslicey = obj[:,obj.shape[1]/2]

    xpixels = cslicex[n.where(cslicex > cmaxx/2.)]
    ypixels = cslicey[n.where(cslicey > cmaxy/2.)]

    xfwhm = len(xpixels)*float(samp)
    yfwhm = len(ypixels)*float(samp)

    #print 'FWHM X = ', xfwhm, 'mas'
    #print 'FWHM Y = ', yfwhm, 'mas'

    return (xfwhm, yfwhm)
#------------#
