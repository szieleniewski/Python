''' Mask function for kinematic plots to mask out random noisey patches

Written by Simon Zieleniewski

Last updated: 15-01-15

'''

import numpy as n


def app_mask(img, thres='default', num_pix=5):
    #image shape
    y, x = img.shape
    #magnitude of image
    img = n.abs(img)

    print n.min(img)
    print n.median(img)
    print n.max(img)
    
    mask = n.ones((y,x), dtype=float)

    #standard threshold value
    if thres == 'default':
        print 'Default threshold'
        thr = n.median(img)/10.
    #custom threshold value
    elif thres != 'default':
        thr = float(thres)
    #loop over x and y dimensions
    for i in xrange(x):
        for j in xrange(y):
            if img[j,i] >= thr:
                npix = n.where(img[j-1:j+2,i-1:i+2] >= thr, 1., 0.)
                if n.sum(npix) < num_pix:
                    mask[j,i] = 0.

    return mask
