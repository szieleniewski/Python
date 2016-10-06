'''Code to create CSPs from MIUSCAt SSPs

Written by Simon Zieleniewski


Last updated 12-01-16

'''


import numpy as n
import pylab as P
import scipy.interpolate as si
import astropy.io.fits as p
import V03tools as v03
import CvD12tools as cvd
import os
import spectools as s
import measure_real_indices as mri



def csp(s_age, e_age, mod='V12', outsig=0., v12imf=1, cvdimf=1):
    '''Function to generate CSPs from V12 and CvD12 models

    Inputs:

        s_age: starting age of SF [Gyr]
        e_age: end age of SF [Gyr]
        mod: model, 'V12' or 'CvD12'
        outsig: common resolution to convolve model to [km/s]
        v12imf: IMf of V12 models [1: x=3, 2: Salp]
        cvdimf: IMF of CvD12 models [1: x=3, 2: Salp, 3: Chab]

    Outputs:

        spectrum class of CSP

    '''

    if mod == 'V12':
        if v12imf == 2:
            v03path = '/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_un/MIUSCAT_un_1.30_fits_safe/'
            print 'V12 IMF: x=2.3'
        if v12imf == 1:
            print 'V12 IMF: x=3'
            v03path = '/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_un/MIUSCAT_un_2.00_fits_safe/'
        files = os.listdir(v03path)
    ###Z = 0 ssps
        salfiles = files[100:150]
    ###Z=0.22
        #salfiles = files[150:200]
        v03ssps = []
        v03ages = []

        #Firstly for V12 models
        for i in salfiles:
            age = float(i.split('T')[1][0:7])
            v03ages.append(age)
            if outsig:
                print 'Convolving to common sigma = %.2f km/s' % outsig
                v03ssps.append(mri.measure_convolved_spec_index(v03path+i, None, 'V12', outsig, aorv='air'))
            elif not outsig:
                v03ssps.append(v03.loadV03spec(v03path+i))

        v03ages = n.array(v03ages)
        print v03ages
        start = n.where(v03ages >= s_age)[0][0]
        end = n.where(v03ages >= e_age)[0][0]
        print start, end
        agestosum = v03ages[start:end+1]
        print agestosum
        cspec = n.zeros(len(v03ssps[0].lam),dtype=float)
        for i in xrange(len(agestosum)-1):
            print agestosum[i], agestosum[i+1]
            print n.where(v03ages == agestosum[i])[0][0]
            print n.where(v03ages == agestosum[i+1])[0][0]
            cspec += ((agestosum[i+1] - agestosum[i])/2.)*\
                     (v03ssps[n.where(v03ages == agestosum[i])[0][0]].flam[0]+v03ssps[n.where(v03ages == agestosum[i+1])[0][0]].flam[0])

        return s.spectrum(lam=v03ssps[0].lam, lamspec=cspec)

    elif mod == 'CvD12':
        cvdssps =[]
        cvdages = []
        cvdpath = '/Users/zieleniewski/Data/stelpopmods/CvD12/'
        cvdfiles = [f for f in os.listdir(cvdpath) if not f.startswith('.')]#Ignore .DS_Store
        del cvdfiles[0:4]#Remove non CvD12 model files
        del cvdfiles[5:8]#Remove [alpha/Fe] model files
        del cvdfiles[-1]#Remove varelem model file
        print cvdfiles
        #Gives list of CvD files of Z = 0 and [alpha/Fe] = 0.
        cvdages = []
        for i in cvdfiles:
            cvdages.append(float(i[1:5]))
            if outsig:
                print 'Convolving to common sigma = %.2f km/s' % outsig
                cvdssps.append(mri.measure_convolved_spec_index(cvdpath+i, None, 'CvD12', outsig))
            elif not outsig:
                cvdssps.append(cvd.loadCvD12spec(cvdpath+i))

        cvdages = n.array(cvdages)
        if s_age < 3.0:
            print 'Lowest CvD12 age = 3 Gyr'
            start = n.where(cvdages == 3.)[0]
        else:
            start = n.where(cvdages >= s_age)[0][0]
        if e_age > 13.5:
            print 'Largest CvD12 age = 13.5 Gyr'
            end = n.where(cvdages == 13.5)[0]
        else:
            end = n.where(cvdages >= e_age)[0][0]
        agestosum = cvdages[start:end+1]
        print agestosum
        cspec = n.zeros(len(cvdssps[0].lam),dtype=float)
        for i in xrange(len(agestosum)-1):
            print agestosum[i], agestosum[i+1]
            print n.where(cvdages == agestosum[i])[0][0]
            print n.where(cvdages == agestosum[i+1])[0][0]
            cspec += ((agestosum[i+1] - agestosum[i])/2.)*\
                     (cvdssps[n.where(cvdages == agestosum[i])[0][0]].flam[cvdimf]+cvdssps[n.where(cvdages == agestosum[i+1])[0][0]].flam[cvdimf])

        return s.spectrum(lam=cvdssps[0].lam, lamspec=cspec)
