'''Code to measure index and error of spectra
from SWIFT IFU for several near-IR indices:
CaT*, NaI and FeH.

Written by Simon Zieleniewski

Last updated: 03-08-16

'''


import numpy as n
import pylab as P
import scipy.interpolate as si
import scipy.constants as sc
import spectools as s
import spectools_SZ as ssz
import CvD12tools as cvd
import V03tools as v03
from os.path import expanduser as eu
import pdb



def measure_convolved_spec_index(specf, varspecf, kins='CvD12', out_sig=200.0, index='CaT_star', aorv='air'):
    '''Function to convolve spectrum up to a common resolution before
    measuring absorption feature strength.

    Inputs:

        specf - FITS file path for science spectrum
        varspecf - FITS file path for variance spectrum
        kins - kinematics file path for real data. Use 'CvD12', 'V12' or 'V15' if using SPS models
        out_sig - final common resolution [km/s]
        index - Absorption feature [Hb, Mgb, Fe, MgFe, NaI, CaT, PaT, CaT_star, MgI, TiO, FeH]
        aorv: air or vacuum wavelength definitions

    Outputs:

        ind - index value [A]
        ind_sig - index sigma [A]
    '''

    if kins == 'CvD12':
        spec = cvd.loadCvD12spec(specf)
        #Convolving sigma [km/s]
        if index in ['Hb', 'Mgb', 'Fe', 'MgFe', 'Fe52', 'Fe53', 'NaD']:
            print 'Optical index'
            specres = 2.5#A
            specsig = specres*sc.c/5200./(2.*n.sqrt(2.*n.log(2.)))/1000.
            print 'sigma = %.2f km/s' % specsig
        else:
            print 'NIR index'
            specres = sc.c/2000.
            specsig = specres/(2.*n.sqrt(2.*n.log(2.)))/1000.
            print 'sigma = %.2f km/s' % specsig
        convsig = n.sqrt(out_sig**2 - specsig**2)

        loglambs = n.linspace(n.log10(spec.lam[0]), n.log10(spec.lam[-1]), len(spec.lam))
        veldisp = (10**(loglambs[1]-loglambs[0]) - 1)*sc.c/1000. #in km/s

        gauss = ssz.Gauss(convsig, veldisp)
        spec.gaussConvolve(0.0, convsig, losvd=gauss[:,1])

        linlambs = n.linspace(spec.loglam[0], spec.loglam[-1], len(spec.lam))
        intspec = si.interp1d(spec.loglam, spec.conflam)
        newflux = intspec(linlambs)
        conspec = s.spectrum(lam=linlambs, lamspec=newflux, mass=spec.mass)

        return conspec


    if kins in ['V12', 'V15']:
        if kins == 'V12':
            spec = v03.loadV03spec(specf)
        elif kins == 'V15':
            spec = v03.loadV15spec(specf)
        #Convolving sigma [km/s]
        lamres = 2.51 #A FWHM
        lamdisp = 0.9 #A/pix

        loglambs = n.linspace(n.log10(spec.lam[0]), n.log10(spec.lam[-1]), len(spec.lam))
        veldisp = (10**(loglambs[1]-loglambs[0]) - 1)*sc.c/1000. #in km/s
        specsig = veldisp*(lamres/lamdisp)/(2.*n.sqrt(2.*n.log(2.)))


        #For V12 models, set sigma for each specific wavelength.
        #sigma = 2.51 [A] * c [km/s] *(1/(2.*n.sqrt(2.*n.log(2.))))/ wavelength_of_feature [A]
        if index in ['NaI', 'NaIsdss']:
            specsig = 39.0 #km/s
        elif index in ['CaT_star', 'CaT', 'PaT']:
            specsig = 37.1
        elif index == 'FeH':
            specsig = 32.2
        elif index == 'MgI':
            specsig = 35.9
        elif index == 'TiO':
            specsig = 36.3
        elif index in ['Hb', 'Mgb', 'Fe', 'MgFe', 'Fe52', 'Fe53', 'NaD']:
            print 'Optical index'
            specsig = 61.5
        print 'specsig = ', specsig

        convsig = n.sqrt(out_sig**2 - specsig**2)

        #gauss = ssz.Gauss(convsig, veldisp)
        gauss = ssz.Gauss(convsig, veldisp)
        spec.gaussConvolve(0.0, convsig, losvd=gauss[:,1])

        linlambs = n.linspace(spec.loglam[0], spec.loglam[-1], len(spec.lam))
        intspec = si.interp1d(spec.loglam, spec.conflam)
        newflux = intspec(linlambs)
        conspec = s.spectrum(lam=linlambs, lamspec=newflux)

        return conspec


    if kins != 'CvD12' or kins != 'V12':
        cat = ssz.cattoolset()
        kinfile = n.genfromtxt(kins)

        #redshift
        z = kinfile[0,0]*1000./sc.c
        #convolving sigma [km/s] (now with instrument resolution = 50.4 km/s for NaCaT; 63.3 km/s for FeH)
        if index == 'FeH':
            instsig = 63.3
        else:
            instsig = 50.4
        print instsig
        convsig = n.sqrt(out_sig**2 - (kinfile[1,0]**2 + instsig**2))
        #Load spectra
        spec = cat.getspec(specf)
        varspec = cat.getspec(varspecf)
        #Redshift correct
        specz = ssz.z_correct(spec, z)
        varspecz = ssz.z_correct(varspec, z)

        specc = s.spectrum(lam=specz[:,0], lamspec=specz[:,1])
        vspecc = s.spectrum(lam=varspecz[:,0], lamspec=varspecz[:,1])


##        ###
##        ###For calculating indices at intrinsic velocity dispersion
##        if index == 'CaT_star':
##            ind, indvar = specc.CaT_star(specc.lam[1]-specc.lam[0], vspecc, aorv)
##            return ind, indvar
##        elif index == 'TiO':
##            ind, indvar = specc.TiOindex(specc.lam[1]-specc.lam[0], vspecc, aorv)
##            return ind, indvar
##        else:
##            ind, indvar = specc.irindex(specc.lam[1]-specc.lam[0], index, vspecc, aorv)
##            return ind, indvar
##        ###
##        ###


        loglambs = n.linspace(n.log10(specc.lam[0]), n.log10(specc.lam[-1]), len(specc.lam))
        veldisp = (10**(loglambs[1]-loglambs[0]) - 1)*sc.c/1000. #in km/s

        gauss = ssz.Gauss(convsig, veldisp)

        specc.gaussConvolve(0.0, convsig, losvd=gauss[:,1])
        vspecc.gaussConvolve(0.0, convsig, losvd=gauss[:,1])

        linlambs = n.linspace(specc.loglam[0], specc.loglam[-1], len(specc.lam))
        intspec = si.interp1d(specc.loglam, specc.conflam)
        intvspec = si.interp1d(vspecc.loglam, vspecc.conflam)
        newflux = intspec(linlambs)
        newvflux = intvspec(linlambs)
        conspec = s.spectrum(lam=linlambs, lamspec=newflux)
        convspec = s.spectrum(lam=linlambs, lamspec=newvflux)

##        P.plot(conspec.lam, conspec.flam[0], 'k-')
##        P.xlabel('$\lambda$',fontsize=18)
##        P.ylabel('F$_\lambda$',fontsize=18)
##        P.xlim([8130, 8280])
##        P.show()


        if index == 'CaT_star':
            print 'CaT*'
            ind, indvar = conspec.CaT_star(conspec.lam[1]-conspec.lam[0], convspec, aorv)
        elif index == 'TiO':
            print 'TiO'
            ind, indvar = conspec.TiOindex(conspec.lam[1]-conspec.lam[0], convspec, aorv)
        elif index == 'Fe':
            print 'Fe'
            ind, indvar = conspec.Fe(conspec.lam[1]-conspec.lam[0], convspec, aorv)
        elif index == 'MgFe':
            print "[MgFe]'"
            ind, indvar = conspec.MgFe(conspec.lam[1]-conspec.lam[0], convspec, aorv)
        else:
            print 'Else: '+index
            ind, indvar = conspec.irindex(conspec.lam[1]-conspec.lam[0], index, convspec, aorv)

        #pdb.set_trace()
        return ind, indvar#, conspec









if __name__=="__main__":

    import sys

###Code to measure index
##    if len(sys.argv) == 9:
##        sed = sys.argv[1]
##        vsed = sys.argv[2]
##        v = float(sys.argv[3])
##        vsig = float(sys.argv[4])
##        vdisp = float(sys.argv[5])
##        vdispsig = float(sys.argv[6])
##        mod = str(sys.argv[7])
##        ind = sys.argv[8]
##
##        print sed
##        print vsed
##        print v
##        print vsig
##        print vdisp
##        print vdispsig
##        print mod
##        print ind
##
##        measure_index(sed, vsed, v, vsig, vdisp, vdispsig, mod, ind)
##
##    elif len(sys.argv) != 9:
##        print""
##        print '    --------     '
##        print 'Valculate index value and error'
##        print '    --------     '
##        print 'USAGE'
##        print 'Enter command line arguments in following order:'
##        print '1. Science spectrum FITS file'
##        print '2. Variance spectrum FITS file'
##        print '3. Velocity [km/s]'
##        print '4. Velocity sigma [km/s]'
##        print '5. Velocity dispersion [km/s]'
##        print '6. Velocity dispersion sigma [km/s]'
##        print '7. Model choice [CvD12, V03]'
##        print '8. Index [CaT_star, NaI, FeH]'
##        print ""


#Code to measure index convolved up to a common resolution

    if len(sys.argv) == 7:
        sed = sys.argv[1]
        vsed = sys.argv[2]
        kins = sys.argv[3]
        outsig = float(sys.argv[4])
        ind = sys.argv[5]
        aorv = sys.argv[6]

        print sed
        print vsed
        print kins
        print outsig
        print ind
        print aorv

        indval, indvar = measure_convolved_spec_index(sed, vsed, kins, outsig, ind, aorv)
        print indval
        print n.sqrt(indvar)

    elif len(sys.argv) != 7:
        print""
        print '    --------     '
        print 'Calculate index value and error'
        print '    --------     '
        print 'USAGE'
        print 'Enter command line arguments in following order:'
        print '1. Science spectrum FITS file'
        print '2. Variance spectrum FITS file'
        print '3. Kinematics file'
        print '4. output resolution sigma [km/s]'
        print '5. Index [CaT_star, NaI, FeH, MgI, TiO, CaT, PaT]'
        print '6. air or vacuum wavelengths [air, vac]'
        print ""


##########
#Code to create text files of index v vel/sigma/sigmacorrfac

##    mod = sys.argv[1]
##    ind = sys.argv[2]
##    var = sys.argv[3]
##    fname = sys.argv[4]
##
##    print mod
##    print ind
##    print var
##    print fname
##
##    CvD12_index_v_data(mod, ind, var, fname)
##
##

#'/Users/zieleniewski/Data/swift/m31/results_v6/data_v6/m31_1_science_spec.T2.fits'
#'/Users/zieleniewski/Data/swift/m31/2012-09-22/m31_1_combined_median_divided_spec_opt0p82_ms065_067.var.fits'
#'/Users/zieleniewski/Data/swift/m31/results_v6/data_v6/m31_1_combined_median_divided_spec_opt0p82_ms065_067.kin.T1.txt'
