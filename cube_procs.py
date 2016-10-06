''' Code to perofrm simple SWIFT datacube procedures

Author: Simon Zieleniewski

Created: 23-10-15

Last updated: 23-10-15

'''

import astropy.io.fits as p
import numpy as n


def div_cube(sci, div, string, var=None):

    cube, head = p.getdata(sci, header=True)
    div, dhead = p.getdata(div, header=True)

    if dhead['NAXIS'] == 1:
        div.shape = dhead['NAXIS1'], 1, 1

    ncube = cube/div

    if var:
        vcube, vhead = p.getdata(var, header=True)
        nvcube = vcube/(div**2)

    print 'Saving new cubes into current directory'


    p.writeto('./'+sci.split('.fits')[0]+'.'+string+'.fits', ncube, head)
    if var:
        p.writeto('./'+var.split('.var.fits')[0]+'.'+string+'.var.fits', nvcube, vhead)




if __name__=="__main__":
    import sys
    print len(sys.argv)
    if len(sys.argv) < 2 or sys.argv[1] not in ['divide']:
        print ''
        print '1. divide: [SCI cube, Dividing cube (e.g. telluric), output suffix string, [VAR cube]]'
        print ''
    elif sys.argv[1] == 'divide':
        if len(sys.argv) == 6:
            div_cube(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
        else:
            div_cube(sys.argv[2], sys.argv[3], sys.argv[4])


                
