'''Function to test whether convolving spectra up to
common velocity dispersion has a bad effect on index
measurements. Also, compare to using model spectra
to compute 'correcting function' for index measurements
at different sigmas

Author: Simon Zieleniewski

Last updated: 17-08-16

'''

import numpy as n
import pylab as P
import spectools as s
import spectools_SZ as ssz
import scipy.interpolate as si
import measure_real_indices as mri
import CvD12tools as cvd
import V03tools as v03
import scipy.constants as sc
import sys


####FONTS
tfsize = 28
lfsize = 24
gensize = 18
majtsize = 14
mintsize = 6
twidth = 2
msize = 10
####



def convolve_spec(spec, sig, losvd, var=None):
    '''Convolve spectrum up to common velocity dispersion

    Inputs:

        - spec: spectrum class
        - sig: sigma [km/s]
        - losvd: Line of sight velocity dispersion (array)
        - var: variance spectrum class

    Outputs:

        - nspec
        - nvar (if var given)

    '''

    spec.gaussConvolve(0.0, sig, losvd=losvd[:,1])
    linlambs = n.linspace(spec.loglam[0], spec.loglam[-1], len(spec.lam))
    intspec = si.interp1d(spec.loglam, spec.conflam)
    newflux = intspec(linlambs)
    if var:
        var.gaussConvolve(0.0, sig, losvd=losvd[:,1])
        vlinlambs = n.linspace(var.loglam[0], var.loglam[-1], len(var.lam))
        intspec = si.interp1d(var.loglam, var.conflam)
        newvar = intspec(vlinlambs)

    nspec = s.spectrum(lam=linlambs, lamspec=newflux)
    if not var:
        return nspec
    else:
        nvar = s.spectrum(lam=vlinlambs, lamspec=newvar)
        return nspec, nvar
#------------#



def test_convolve(ind):
    '''Function to test whether convolving to a common
    velocity dispersion has a bad effect on index measurements.

    Inputs:

        -

    Outputs:

        - Plot showing

    '''

    cvd13 = cvd.loadCvD12spec('/Users/zieleniewski/Data/stelpopmods/CvD12/t13.5_solar.ssp')

    chab = n.column_stack((cvd13.lam, cvd13.flam[3]))
    bh = n.column_stack((cvd13.lam, cvd13.flam[1]))

    # chabind =  cvd13.irindex(cvd13.lam[1]-cvd13.lam[0], index=ind)[0][3]
    # bhind = cvd13.irindex(cvd13.lam[1]-cvd13.lam[0], index=ind)[0][1]

    if ind in ['NaI', 'NaIsdss']:
        start = n.where(chab[:,0] > 8050.)[0][0]
        end = n.where(chab[:,0] > 8400.)[0][0]
    #start = n.where(chab[:,0] > 8300.)[0][0]
    #end = n.where(chab[:,0] > 8900.)[0][0]
    if ind == 'FeH':
        start = n.where(chab[:,0] > 9800.)[0][0]
        end = n.where(chab[:,0] > 10000.)[0][0]

    chab = chab[start:end,:]
    bh = bh[start:end,:]

    chab[:,1] /= n.median(chab[:,1])
    bh[:,1] /= n.median(bh[:,1])

    #Want Error as a function of sigma for three different random noise values: 1, 2, 5%.
    #Loop over sigma array
    #For each sigma generate random noise 50 times and compute index at given sigma then store error results_v6
    #Print errors as function of sigma for three noise values (and two IMFs)
    num=50
    sigs = n.arange(180., 430., 24.)
    meance1 = n.zeros_like(sigs)
    rance1 = n.zeros_like(sigs)
    sysce1 = n.zeros_like(sigs)
    meance2 = n.zeros_like(sigs)
    rance2 = n.zeros_like(sigs)
    sysce2 = n.zeros_like(sigs)
    meance3 = n.zeros_like(sigs)
    rance3 = n.zeros_like(sigs)
    sysce3 = n.zeros_like(sigs)

    meanbe1 = n.zeros_like(sigs)
    ranbe1 = n.zeros_like(sigs)
    sysbe1 = n.zeros_like(sigs)
    meanbe2 = n.zeros_like(sigs)
    ranbe2 = n.zeros_like(sigs)
    sysbe2 = n.zeros_like(sigs)
    meanbe3 = n.zeros_like(sigs)
    ranbe3 = n.zeros_like(sigs)
    sysbe3 = n.zeros_like(sigs)


    for i in xrange(len(sigs)):

        cinds1 = n.zeros(num)
        cinds2 = n.zeros(num)
        cinds3 = n.zeros(num)
        bhinds1 = n.zeros(num)
        bhinds2 = n.zeros(num)
        bhinds3 = n.zeros(num)

        specres = sc.c/2000.
        specsig = specres/(2.*n.sqrt(2.*n.log(2.)))/1000.
        convsig = n.sqrt(sigs[i]**2 - specsig**2)
        loglambs = n.linspace(n.log10(chab[0,0]), n.log10(chab[-1,0]), len(chab[:,0]))
        veldisp = (10**(loglambs[1]-loglambs[0]) - 1)*sc.c/1000.
        gauss = ssz.Gauss(convsig, veldisp)

        schab = s.spectrum(chab[:,0], chab[:,1])
        sbh = s.spectrum(bh[:,0], bh[:,1])
        schab = convolve_spec(schab, sigs[i], gauss)
        sbh = convolve_spec(sbh, sigs[i], gauss)
        disp = cvd13.lam[1]-cvd13.lam[0]
        chabind = schab.irindex(disp, ind)[0]
        bhind = sbh.irindex(disp, ind)[0]

        for nn in xrange(num):

            #1, 2, 5% random noise (SNR = 100, 50, 20)
            sig1 = 0.01
            sig2 = 0.02
            sig3 = 0.05
            n1 = n.random.normal(n.zeros(len(bh[:,0])), sig1)
            n2 = n.random.normal(n.zeros(len(bh[:,0])), sig2)
            n3 = n.random.normal(n.zeros(len(bh[:,0])), sig3)
            chab1 = n.column_stack((chab[:,0], chab[:,1]+n1))
            chab2 = n.column_stack((chab[:,0], chab[:,1]+n2))
            chab3 = n.column_stack((chab[:,0], chab[:,1]+n3))
            bh1 = n.column_stack((bh[:,0], bh[:,1]+n1))
            bh2 = n.column_stack((bh[:,0], bh[:,1]+n2))
            bh3 = n.column_stack((bh[:,0], bh[:,1]+n3))

            schab1 = s.spectrum(chab1[:,0], chab1[:,1])
            schab2 = s.spectrum(chab2[:,0], chab2[:,1])
            schab3 = s.spectrum(chab3[:,0], chab3[:,1])
            sbh1 = s.spectrum(bh1[:,0], bh1[:,1])
            sbh2 = s.spectrum(bh2[:,0], bh2[:,1])
            sbh3 = s.spectrum(bh3[:,0], bh3[:,1])

            schab1 = convolve_spec(schab1, sigs[i], gauss)
            schab2 = convolve_spec(schab2, sigs[i], gauss)
            schab3 = convolve_spec(schab3, sigs[i], gauss)
            sbh1 = convolve_spec(sbh1, sigs[i], gauss)
            sbh2 = convolve_spec(sbh2, sigs[i], gauss)
            sbh3 = convolve_spec(sbh3, sigs[i], gauss)

            cinds1[nn] = schab1.irindex(disp, ind)[0]
            cinds2[nn] = schab2.irindex(disp, ind)[0]
            cinds3[nn] = schab3.irindex(disp, ind)[0]
            bhinds1[nn] = sbh1.irindex(disp, ind)[0]
            bhinds2[nn] = sbh2.irindex(disp, ind)[0]
            bhinds3[nn] = sbh3.irindex(disp, ind)[0]

        meance1[i] = n.mean(abs(cinds1-chabind)/chabind)
        meance2[i] = n.mean(abs(cinds2-chabind)/chabind)
        meance3[i] = n.mean(abs(cinds3-chabind)/chabind)

        rance1[i] = n.std(cinds1)
        rance2[i] = n.std(cinds2)
        rance3[i] = n.std(cinds3)

        sysce1[i] = n.mean(cinds1-chabind)
        sysce2[i] = n.mean(cinds2-chabind)
        sysce3[i] = n.mean(cinds3-chabind)

        meanbe1[i] = n.mean(abs(bhinds1-bhind)/bhind)
        meanbe2[i] = n.mean(abs(bhinds2-bhind)/bhind)
        meanbe3[i] = n.mean(abs(bhinds3-bhind)/bhind)

        ranbe1[i] = n.std(bhinds1)
        ranbe2[i] = n.std(bhinds2)
        ranbe3[i] = n.std(bhinds3)

        sysbe1[i] = n.mean(bhinds1-bhind)
        sysbe2[i] = n.mean(bhinds2-bhind)
        sysbe3[i] = n.mean(bhinds3-bhind)

    print 'SIGMAS [km/s] = ', sigs
    print 'MEAN CHAB ERR (1% noise) = ', meance1
    print 'RAN CHAB ERR (1%) = ', rance1
    print 'SYS CHAB ERR (1%) = ', sysce1

    print 'MEAN CHAB ERR (2% noise) = ', meance2
    print 'RAN CHAB ERR (2%) = ', rance2
    print 'SYS CHAB ERR (2%) = ', sysce2

    print 'MEAN CHAB ERR (5% noise) = ', meance3
    print 'RAN CHAB ERR (5%) = ', rance3
    print 'SYS CHAB ERR (5%) = ', sysce3

    print 'MEAN BH ERR (1% noise) = ', meanbe1
    print 'RAN BH ERR (1%) = ', ranbe1
    print 'SYS BH ERR (1%) = ', sysbe1

    print 'MEAN BH ERR (2% noise) = ', meanbe2
    print 'RAN BH ERR (2%) = ', ranbe2
    print 'SYS BH ERR (2%) = ', sysbe2

    print 'MEAN BH ERR (5% noise) = ', meanbe3
    print 'RAN BH ERR (5%) = ', ranbe3
    print 'SYS BH ERR (5%) = ', sysbe3

    #Plot MEAN error v sigma
    P.plot(sigs, meance1*100, 'k-',
           sigs, meance2*100, 'r-',
           sigs, meance3*100, 'b-', lw=2.0)
    P.xlabel('$\sigma$ [km/s]', fontsize=tfsize)
    P.ylabel('Mean Error [%]', fontsize=tfsize)
    P.xticks(size=lfsize)
    P.yticks(size=lfsize)
    P.minorticks_on()
    P.tick_params(axis='both', direction='in', length=majtsize,
                  width=1.5, which='major')
    P.tick_params(axis='both', direction='in', length=mintsize,
                  width=1.5, which='minor')
    leg = P.legend(['1% noise (S/N=100)', '2% noise (S/N=50)', '5% noise (S/N=20)'],
    loc='lower right', frameon=False, borderaxespad=0.3, labelspacing=0.2, numpoints=1, scatterpoints=1)
    P.tight_layout()
    P.show()

    #Plot RANDOM error v sigma
    P.plot(sigs, rance1*100, 'k-',
           sigs, rance2*100, 'r-',
           sigs, rance3*100, 'b-', lw=2.0)
    P.xlabel('$\sigma$ [km/s]', fontsize=tfsize)
    P.ylabel('Random Error [%]', fontsize=tfsize)
    P.xticks(size=lfsize)
    P.yticks(size=lfsize)
    P.minorticks_on()
    P.tick_params(axis='both', direction='in', length=majtsize,
                  width=1.5, which='major')
    P.tick_params(axis='both', direction='in', length=mintsize,
                  width=1.5, which='minor')
    leg = P.legend(['1% noise (S/N=100)', '2% noise (S/N=50)', '5% noise (S/N=20)'],
    loc='lower right', frameon=False, borderaxespad=0.3, labelspacing=0.2, numpoints=1, scatterpoints=1)
    P.tight_layout()
    P.show()

    #Plot SYSTEMATIC error v sigma
    P.plot(sigs, sysce1*100, 'k-',
           sigs, sysce2*100, 'r-',
           sigs, sysce3*100, 'b-', lw=2.0)
    P.xlabel('$\sigma$ [km/s]', fontsize=tfsize)
    P.ylabel('Systematic Error [%]', fontsize=tfsize)
    P.xticks(size=lfsize)
    P.yticks(size=lfsize)
    P.minorticks_on()
    P.tick_params(axis='both', direction='in', length=majtsize,
                  width=1.5, which='major')
    P.tick_params(axis='both', direction='in', length=mintsize,
                  width=1.5, which='minor')
    leg = P.legend(['1% noise (S/N=100)', '2% noise (S/N=50)', '5% noise (S/N=20)'],
    loc='lower right', frameon=False, borderaxespad=0.3, labelspacing=0.2, numpoints=1, scatterpoints=1)
    P.tight_layout()
    P.show()








############
def CvD12_indexvsig(ind, insig, outsig, savevals=False, plot=False):
    '''Code to create correction factor as a function of sigma
    using CvD12 models

    Inputs:

        - ind: index
        - insig: input sigma of spectrum [km/s]
        - outsig: output sigma of spectrum [km/s]

    Outputs:

        - correction factor to multiply index by

    '''

    insig = float(insig)
    outsig = float(outsig)

    delsig = 30.0
    if savevals == 'True':
        sigs = n.arange(80., 430., delsig)
        inds9c = n.zeros_like(sigs)
        inds9x3 = n.zeros_like(sigs)
        inds13c = n.zeros_like(sigs)
        inds13x3 = n.zeros_like(sigs)
        inds13afc = n.zeros_like(sigs)
        inds13afx3 = n.zeros_like(sigs)

        # vinds14c = n.zeros_like(sigs)
        # vinds14x3 = n.zeros_like(sigs)

        chab = 3 #Chabrier IMF
        x3 = 1 #x = 3 IMF
        file1 = '/Users/zieleniewski/Data/stelpopmods/CvD12/t09.0_solar.ssp'
        file2 = '/Users/zieleniewski/Data/stelpopmods/CvD12/t13.5_solar.ssp'
        file3 = '/Users/zieleniewski/Data/stelpopmods/CvD12/t13.5_afe+0.2.ssp'

        #Calculate indices at intrinsic model resolution
        spec9 = cvd.loadCvD12spec(file1)
        spec13 = cvd.loadCvD12spec(file2)
        spec13af = cvd.loadCvD12spec(file3)
        if ind in ['NaD', 'NaI', 'NaIsdss', 'MgI', 'FeH', 'CaT', 'PaT', 'MgFe', 'Fe']:
            ind9c =  spec9.irindex(spec9.lam[1]-spec9.lam[0], index=ind)[0][chab]
            ind9x3 = spec9.irindex(spec9.lam[1]-spec9.lam[0], index=ind)[0][x3]
            indc =  spec13.irindex(spec13.lam[1]-spec13.lam[0], index=ind)[0][chab]
            indx3 = spec13.irindex(spec13.lam[1]-spec13.lam[0], index=ind)[0][x3]
            ind13afc =  spec13af.irindex(spec13af.lam[1]-spec13af.lam[0], index=ind)[0][chab]
            ind13afx3 = spec13af.irindex(spec13af.lam[1]-spec13af.lam[0], index=ind)[0][x3]
        elif ind == 'TiO':
            ind9c =  spec9.TiOindex(spec9.lam[1]-spec9.lam[0])[0][chab]
            ind9x3 = spec9.TiOindex(spec9.lam[1]-spec9.lam[0])[0][x3]
            indc =  spec13.TiOindex(spec13.lam[1]-spec13.lam[0])[0][chab]
            indx3 = spec13.TiOindex(spec13.lam[1]-spec13.lam[0])[0][x3]
            ind13afc =  spec13af.TiOindex(spec13af.lam[1]-spec13af.lam[0])[0][chab]
            ind13afx3 = spec13af.TiOindex(spec13af.lam[1]-spec13af.lam[0])[0][x3]

        # spec14c = v03.loadV15spec('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MIUSCAT_Padova00_bi_1.30_fits/Ibi1.30Zp0.00T14.1254_iPp0.00_baseFe.fits')
        # spec14x3 = v03.loadV15spec('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MIUSCAT_Padova00_un_2.00_fits/Iun2.00Zp0.00T14.1254_iPp0.00_baseFe.fits')
        # vindc =  spec14c.irindex(spec14c.lam[1]-spec14c.lam[0], index=ind)[0]
        # vindx3 = spec14x3.irindex(spec14x3.lam[1]-spec14x3.lam[0], index=ind)[0]

        for i in xrange(len(sigs)):
            cvd9 = mri.measure_convolved_spec_index(file1, varspecf=None, kins='CvD12', out_sig=sigs[i], index=ind)
            cvd13 = mri.measure_convolved_spec_index(file2, varspecf=None, kins='CvD12', out_sig=sigs[i], index=ind)
            cvd13af = mri.measure_convolved_spec_index(file3, varspecf=None, kins='CvD12', out_sig=sigs[i], index=ind)

            # v14c = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MIUSCAT_Padova00_bi_1.30_fits/Ibi1.30Zp0.00T14.1254_iPp0.00_baseFe.fits',
            #                                         varspecf=None, kins='V15', out_sig=sigs[i], index=ind)
            # v14x3 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MIUSCAT_Padova00_un_2.00_fits/Iun2.00Zp0.00T14.1254_iPp0.00_baseFe.fits',
            #                                         varspecf=None, kins='V15', out_sig=sigs[i], index=ind)

            if ind in ['NaD', 'NaI', 'NaIsdss', 'MgI', 'FeH', 'CaT', 'PaT', 'MgFe', 'Fe']:
                inds9c[i] = ind9c/float(cvd9.irindex(cvd9.lam[1]-cvd9.lam[0], index=ind)[0][chab])
                inds9x3[i] = ind9x3/float(cvd9.irindex(cvd9.lam[1]-cvd9.lam[0], index=ind)[0][x3])
                inds13c[i] = indc/float(cvd13.irindex(cvd13.lam[1]-cvd13.lam[0], index=ind)[0][chab])
                inds13x3[i] = indx3/float(cvd13.irindex(cvd13.lam[1]-cvd13.lam[0], index=ind)[0][x3])
                inds13afc[i] = ind13afc/float(cvd13af.irindex(cvd13af.lam[1]-cvd13af.lam[0], index=ind)[0][chab])
                inds13afx3[i] = ind13afx3/float(cvd13af.irindex(cvd13af.lam[1]-cvd13af.lam[0], index=ind)[0][x3])

                # vinds14c[i] = vindc/float(v14c.irindex(v14c.lam[1]-v14c.lam[0], index=ind, vac_or_air='vac')[0])
                # vinds14x3[i] = vindx3/float(v14x3.irindex(v14x3.lam[1]-v14x3.lam[0], index=ind, vac_or_air='vac')[0])
            elif ind == 'TiO':
                inds9c[i] = ind9c/float(cvd9.TiOindex(cvd9.lam[1]-cvd9.lam[0])[0][chab])
                inds9x3[i] = ind9x3/float(cvd9.TiOindex(cvd9.lam[1]-cvd9.lam[0])[0][x3])
                inds13c[i] = indc/float(cvd13.TiOindex(cvd13.lam[1]-cvd13.lam[0])[0][chab])
                inds13x3[i] = indx3/float(cvd13.TiOindex(cvd13.lam[1]-cvd13.lam[0])[0][x3])
                inds13afc[i] = ind13afc/float(cvd13af.TiOindex(cvd13af.lam[1]-cvd13af.lam[0])[0][chab])
                inds13afx3[i] = ind13afx3/float(cvd13af.TiOindex(cvd13af.lam[1]-cvd13af.lam[0])[0][x3])

                # vinds14c[i] = vindc/float(v14c.TiOindex(v14c.lam[1]-v14c.lam[0], vac_or_air='vac')[0])
                # vinds14x3[i] = vindx3/float(v14x3.TiOindex(v14x3.lam[1]-v14x3.lam[0], vac_or_air='vac')[0])

        n.savetxt('./CvD12_Model_'+ind+'_v_sigma_delsig'+str(int(delsig))+'.txt', n.column_stack((sigs, inds9c, inds9x3, inds13c, inds13x3, inds13afc, inds13afx3)), fmt='%.6f', delimiter=',')
        print 'Values Saved!'
        sys.exit()

    else:
        dat = n.genfromtxt('/Users/zieleniewski/Documents/Oxford/D.Phil/stelpops/CvD12_Model_'+ind+'_v_sigma_delsig'+str(int(delsig))+'.txt', delimiter=',')
        sigs = dat[:,0]
        inds9c = dat[:,1]
        inds9x3 = dat[:,2]
        inds13c = dat[:,3]
        inds13x3 = dat[:,4]
        inds13afc = dat[:,5]
        inds13afx3 = dat[:,6]

    avg = (inds13c+inds13x3+inds9c+inds9x3+inds13afc+inds13afx3)/6.
    avgint = si.interp1d(sigs, avg, 1)
    inds9cint = si.interp1d(sigs, inds9c, 1)
    inds9x3int = si.interp1d(sigs, inds9x3, 1)
    inds13cint = si.interp1d(sigs, inds13c, 1)
    inds13x3int = si.interp1d(sigs, inds13x3, 1)
    inds13afcint = si.interp1d(sigs, inds13afc, 1)
    inds13afx3int = si.interp1d(sigs, inds13afx3, 1)

    inpoint = float(avgint(insig))
    inpoint9c = float(inds9cint(insig))
    inpoint9x3 = float(inds9x3int(insig))
    inpoint13c = float(inds13cint(insig))
    inpoint13x3 = float(inds13x3int(insig))
    inpoint13afc = float(inds13afcint(insig))
    inpoint13afx3 = float(inds13afx3int(insig))

    inds9c = inds9c / inpoint9c
    inds9x3 = inds9x3 / inpoint9x3
    inds13c = inds13c / inpoint13c
    inds13x3 = inds13x3 / inpoint13x3
    inds13afc = inds13afc / inpoint13afc
    inds13afx3 = inds13afx3 / inpoint13afx3
    avg = avg / inpoint

    #Interpolate index values as a function of sigma_0
    avgint = si.interp1d(sigs, avg, 1)
    # inds9cint = si.interp1d(sigs, inds9c, 1)
    # inds9x3int = si.interp1d(sigs, inds9x3, 1)
    # inds13cint = si.interp1d(sigs, inds13c, 1)
    # inds13x3int = si.interp1d(sigs, inds13x3, 1)
    # inds13afcint = si.interp1d(sigs, inds13afc, 1)
    # inds13afx3int = si.interp1d(sigs, inds13afx3, 1)

    outpoint = avgint(outsig)
    # outpoint9c = inds9cint(outsig)
    # outpoint9x3 = inds9x3int(outsig)
    # outpoint13c = inds13cint(outsig)
    # outpoint13x3 = inds13x3int(outsig)
    # outpoint13afc = inds13afcint(outsig)
    # outpoint13afx3 = inds13afx3int(outsig)

    outfac = n.round(outpoint, 4)
    outerr = n.round(n.std([inds13c,inds13x3,inds9c,inds9x3,inds13afc,inds13afx3]), 4)
    print 'OUTSIG = ', n.round(outsig, 1), ' km/s'
    print 'OUTFAC = ', outfac
    print 'OUTERR = ', outerr

    if plot=='True':
        #Plot
        P.rcParams['font.family'] = 'Times New Roman'
        P.rcParams.update({'font.size': gensize})
        P.figure(figsize=(13,10), dpi=72)
        P.subplots_adjust(left=0.14, bottom=0.12)
        #Plot lines for all
        P.plot(sigs, inds9c, 'k--',
               sigs, inds9x3, 'r--',
            #    sigs, (inds9c+inds9x3)/2., 'b--',
               sigs, inds13c, 'k-',
               sigs, inds13x3, 'r-',
            #    sigs, (inds13c+inds13x3)/2., 'b-',
               sigs, inds13afc, 'k:',
               sigs, inds13afx3, 'r:',
            #    sigs, (inds13afc+inds13afx3)/2., 'b:',
               sigs, (inds13c+inds13x3+inds9c+inds9x3+inds13afc+inds13afx3)/6., 'g-',
            #    sigs, vinds14c, 'k--',
            #    sigs, vinds14x3, 'r--',
            #    sigs, (vinds14c+vinds14x3)/2., 'b--',
               lw=2.0)
               #sigs, inds13afe[IMF], 'g-', lw=2.0)

        P.xlabel('$\sigma$ [km/s]', fontsize=tfsize)
        if ind != 'TiO':
            P.ylabel('C($\sigma$) = I($\sigma_0$)/I($\sigma$)', fontsize=tfsize)
        P.xticks(size=lfsize)
        P.yticks(size=lfsize)
        P.minorticks_on()
        P.tick_params(axis='both', direction='in', length=majtsize,
                      width=1.5, which='major')
        P.tick_params(axis='both', direction='in', length=mintsize,
                      width=1.5, which='minor')

        leg = P.legend(['CvD12 9Gyr Chab', '9Gyr $x=3$', '13Gyr Chab', '13Gyr $x=3$',
                        '13Gyr [$\\alpha$/Fe]=+0.2 Chab', '13Gyr [$\\alpha$/Fe]=+0.2 $x=3$', 'Average'],
        loc='lower right', frameon=False, borderaxespad=0.3, labelspacing=0.2, numpoints=1, scatterpoints=1)
        # for i in xrange(len(leg.legendHandles)):
        #     leg.legendHandles[i].set_color(r_colour[0])

        P.tight_layout()
        P.show()

    return outfac, outerr


if __name__=="__main__":

    import sys

    if len(sys.argv)<2 or sys.argv[1] not in ['corrfac']:
        print ''
        print '1. corrfac [index, insigma, outsigma, savevals(bool), plot(bool)]'
        print ''
        sys.exit()
    elif sys.argv[1] == 'corrfac':
        CvD12_indexvsig(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6])
