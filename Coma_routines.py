'''Code to generate various plots for
SWIFT Coma data.


Written by Simon Zieleniewski

Last updated 15-08-16

'''

import numpy as n
import pylab as P
import time
import astropy.io.fits as p
from scipy.io.idl import readsav
import sys
sys.path.append('/Users/zieleniewski/Documents/Oxford/D.Phil/stelpops')
import spectools_SZ as ssz
import spectools as s
import measure_real_indices as mri
import scipy.constants as sc
import scipy.interpolate as si
from CvD12_Index_v_sigma import CvD12_indexvsig
sys.path.append('/Users/zieleniewski/Data/SDSS_Coma_spectra')
import sdss_handle_spec as shs
import NWCC
from matplotlib.ticker import AutoMinorLocator
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle



####FONTS
tfsize = 28
lfsize = 24
gensize = 20
majtsize = 14
mintsize = 6
twidth = 2
msize = 10
####

#------------#
def get_data(gal):

    if gal == 'GMP2921':
        # nacat = readsav('/Users/zieleniewski/Data/swift/Coma_stellar_pops/'+gal+'/NaCaT_v2/Conroy_0.79_0.905/ppxfSECT_SINFONI.idl')
        nacat = readsav('/Users/zieleniewski/Data/swift/Coma_stellar_pops/'+gal+'/NaCaT_molecfit/Conroy_0.79_0.905/ppxfSECT_SINFONI.idl')
        na = readsav('/Users/zieleniewski/Data/swift/Coma_stellar_pops/'+gal+'/Na/Conroy_0.8_0.86/ppxfSECT_SINFONI.idl')
        fehc = readsav('/Users/zieleniewski/Data/swift/Coma_stellar_pops/'+gal+'/FeH_cont_v2_LegCont_v5/Conroy_0.98_1.03/ppxfSECT_SINFONI.idl')
        # fehc = readsav('/Users/zieleniewski/Data/swift/Coma_stellar_pops/'+gal+'/FeH_cont_v2_LegCont_v6/Conroy_0.98_1.03/ppxfSECT_SINFONI.idl')
        # fehos = readsav('/Users/zieleniewski/Data/swift/Coma_stellar_pops/'+gal+'/FeH_cont_v2_LegCont_O-S_fitted/Conroy_0.98_1.03/ppxfSECT_SINFONI.idl')
        fehosos = readsav('/Users/zieleniewski/Data/swift/Coma_stellar_pops/'+gal+'/FeH_cont_v2_LegCont_O-S_fitted_oldskycube/Conroy_0.98_1.03/ppxfSECT_SINFONI.idl')
        feho = readsav('/Users/zieleniewski/Data/swift/Coma_stellar_pops/'+gal+'/FeH_opt_v2/Conroy_0.98_1.03/ppxfSECT_SINFONI.idl')

    if gal == 'GMP3329':
        # nacat = readsav('/Users/zieleniewski/Data/swift/Coma_stellar_pops/'+gal+'/NaCaT_v3/Conroy_0.80_0.92/ppxfSECT_SINFONI.idl')
        nacat = readsav('/Users/zieleniewski/Data/swift/Coma_stellar_pops/'+gal+'/NaCaT_molecfit/Conroy_0.80_0.92/ppxfSECT_SINFONI.idl')
        na = readsav('/Users/zieleniewski/Data/swift/Coma_stellar_pops/'+gal+'/Na_v2/Conroy_0.8_0.86/ppxfSECT_SINFONI.idl')
        fehc = readsav('/Users/zieleniewski/Data/swift/Coma_stellar_pops/'+gal+'/FeH_cont_v2_LegCont_v6/Conroy_0.98_1.04/ppxfSECT_SINFONI.idl')
        # fehc = readsav('/Users/zieleniewski/Data/swift/Coma_stellar_pops/'+gal+'/FeH_cont_v2_LegCont_v9/Conroy_0.98_1.04/ppxfSECT_SINFONI.idl')
        # fehos = readsav('/Users/zieleniewski/Data/swift/Coma_stellar_pops/'+gal+'/FeH_cont_v2_LegCont_O-S_fitted/Conroy_0.98_1.04/ppxfSECT_SINFONI.idl')
        fehosos = readsav('/Users/zieleniewski/Data/swift/Coma_stellar_pops/'+gal+'/FeH_cont_v2_LegCont_O-S_fitted_noskycube_nomedsub/Conroy_0.98_1.04/ppxfSECT_SINFONI.idl')
        feho = readsav('/Users/zieleniewski/Data/swift/Coma_stellar_pops/'+gal+'/FeH_opt_v2/Conroy_0.98_1.04/ppxfSECT_SINFONI.idl')

    if gal == 'GMP4928':
        # nacat = readsav('/Users/zieleniewski/Data/swift/Coma_stellar_pops/'+gal+'/NaCaT_v5/Conroy_0.80_0.92/ppxfSECT_SINFONI.idl')
        # nacat = readsav('/Users/zieleniewski/Data/swift/Coma_stellar_pops/'+gal+'/NaCaT_ms059_060/Conroy_0.80_0.92/ppxfSECT_SINFONI.idl')
        # nacat = readsav('/Users/zieleniewski/Data/swift/Coma_stellar_pops/'+gal+'/NaCaT_ms064_065/Conroy_0.80_0.92/ppxfSECT_SINFONI.idl')
        nacat = readsav('/Users/zieleniewski/Data/swift/Coma_stellar_pops/'+gal+'/NaCaT_molecfit/Conroy_0.80_0.92/ppxfSECT_SINFONI.idl')
        na = readsav('/Users/zieleniewski/Data/swift/Coma_stellar_pops/'+gal+'/Na_v3/Conroy_0.8_0.86/ppxfSECT_SINFONI.idl')
        fehc = readsav('/Users/zieleniewski/Data/swift/Coma_stellar_pops/'+gal+'/FeH_cont_v2_LegCont_v6_newbpm_v2/Conroy_0.98_1.04/ppxfSECT_SINFONI.idl')
        # fehc = readsav('/Users/zieleniewski/Data/swift/Coma_stellar_pops/'+gal+'/FeH_cont_v2_LegCont_ms059_060/Conroy_0.98_1.04/ppxfSECT_SINFONI.idl')
        # fehc = readsav('/Users/zieleniewski/Data/swift/Coma_stellar_pops/'+gal+'/FeH_cont_v2_LegCont_ms064_065/Conroy_0.98_1.04/ppxfSECT_SINFONI.idl')
        # fehc = readsav('/Users/zieleniewski/Data/swift/Coma_stellar_pops/'+gal+'/FeH_cont_v2_LegCont_v7_newbpm_v2/Conroy_0.98_1.04/ppxfSECT_SINFONI.idl') #New cube with sky continuum subtracted by Sam Vaughan
        # fehc = readsav('/Users/zieleniewski/Data/swift/Coma_stellar_pops/'+gal+'/FeH_cont_v2_LegCont_v8_newbpm_v2/Conroy_0.98_1.04/ppxfSECT_SINFONI.idl') #New cube with sky continuum subtracted by Sam Vaughan AND removing SSPLIT definitions that fall over FeH feature
        # fehc = readsav('/Users/zieleniewski/Data/swift/Coma_stellar_pops/'+gal+'/FeH_cont_v2_LegCont_v9_newbpm_v2/Conroy_0.98_1.04/ppxfSECT_SINFONI.idl') #Original cube with removing SSPLIT definitions that fall over FeH feature
        # fehos = readsav('/Users/zieleniewski/Data/swift/Coma_stellar_pops/'+gal+'/FeH_cont_v2_LegCont_O-S_fitted_oldskycube_newbpm/Conroy_0.98_1.04/ppxfSECT_SINFONI.idl')
        fehosos = readsav('/Users/zieleniewski/Data/swift/Coma_stellar_pops/'+gal+'/FeH_cont_v2_LegCont_O-S_fitted_noskycube_newbpm_nomedsub/Conroy_0.98_1.04/ppxfSECT_SINFONI.idl')
        feho = readsav('/Users/zieleniewski/Data/swift/Coma_stellar_pops/'+gal+'/FeH_opt_v2_newbpm/Conroy_0.98_1.04/ppxfSECT_SINFONI.idl')

    if gal == 'GMP3367':
        nacat = readsav('/Users/zieleniewski/Data/swift/Coma_stellar_pops/'+gal+'/NaCaT_v2/Conroy_0.80_0.92/ppxfSECT_SINFONI.idl')
        # nacat = readsav('/Users/zieleniewski/Data/swift/Coma_stellar_pops/'+gal+'/NaCaT_v3/Conroy_0.80_0.92/ppxfSECT_SINFONI.idl')
        na = readsav('/Users/zieleniewski/Data/swift/Coma_stellar_pops/'+gal+'/Na/Conroy_0.8_0.86/ppxfSECT_SINFONI.idl')
        fehc = readsav('/Users/zieleniewski/Data/swift/Coma_stellar_pops/'+gal+'/FeH_cont_LegCont_v3/Conroy_0.98_1.03/ppxfSECT_SINFONI.idl')
        # fehos = readsav('/Users/zieleniewski/Data/swift/Coma_stellar_pops/'+gal+'/FeH_cont_LegCont_v2/Conroy_0.98_1.03/ppxfSECT_SINFONI.idl')
        fehosos = readsav('/Users/zieleniewski/Data/swift/Coma_stellar_pops/'+gal+'/FeH_cont_LegCont_v2/Conroy_0.98_1.03/ppxfSECT_SINFONI.idl')
        feho = readsav('/Users/zieleniewski/Data/swift/Coma_stellar_pops/'+gal+'/FeH_opt/Conroy_0.98_1.04/ppxfSECT_SINFONI.idl')
    skylines_nai = n.genfromtxt('/Users/zieleniewski/Data/swift/Coma_stellar_pops/'+gal+'/skylines_NaI.txt', delimiter=',')
    skylines_cat = n.genfromtxt('/Users/zieleniewski/Data/swift/Coma_stellar_pops/'+gal+'/skylines_CaT.txt', delimiter=',')
    skylines_mgi = n.genfromtxt('/Users/zieleniewski/Data/swift/Coma_stellar_pops/'+gal+'/skylines_MgI.txt', delimiter=',')
    skylines_tio = n.genfromtxt('/Users/zieleniewski/Data/swift/Coma_stellar_pops/'+gal+'/skylines_TiO.txt', delimiter=',')
    skylines_feh = n.genfromtxt('/Users/zieleniewski/Data/swift/Coma_stellar_pops/'+gal+'/skylines_FeH.txt', delimiter=',')

    # return na, nacat, fehc, fehos, fehosos, feho, skylines_nai, skylines_cat, skylines_mgi, skylines_tio, skylines_feh
    return na, nacat, fehc, fehosos, feho, skylines_nai, skylines_cat, skylines_mgi, skylines_tio, skylines_feh
#------------#



#------------#
def indices_v_radius(gal, out_sig, save='True', plot='True', skymask='True', fehopt='False', osfitted='False', gspecs='False', corrfac='False', corrval=0.0):
    #Galaxy [GMP2921, GMP3329, GMP4928, GMP3367
    #Common velocity dispersion [km/s]
    #out_sig = 400.0#For GMP2921
    #out_sig = 270.0#For GMP3329
    #out_sig = 250.0#for GMP4928
    #out_sig = ###.0#for GMP3367_FR
    print gal, out_sig, save, plot, skymask, fehopt, osfitted, gspecs, corrfac, corrval

    # na, nacat, fehc, fehcos, fehcosos, feho, skylines_nai, skylines_cat, skylines_mgi, skylines_tio, skylines_feh = get_data(gal)
    na, nacat, fehc, fehcosos, feho, skylines_nai, skylines_cat, skylines_mgi, skylines_tio, skylines_feh = get_data(gal)

    v03 = ssz.cattoolset()

    lams = n.round(nacat["wave"]*10000., 3)
    rs = nacat['rbar']*0.235
    rs[0] = 0.0
    print 'RS = ', rs

    naiinds = n.zeros((len(nacat['binflux'][0,:])),dtype=float)
    naiinds_var = n.zeros((len(nacat['binflux'][0,:])),dtype=float)
    naiminds = n.zeros((len(nacat['binflux'][0,:])),dtype=float)
    catinds = n.zeros((len(nacat['binflux'][0,:])),dtype=float)
    catinds_var = n.zeros((len(nacat['binflux'][0,:])),dtype=float)
    catminds = n.zeros((len(nacat['binflux'][0,:])),dtype=float)
    mgiinds = n.zeros((len(nacat['binflux'][0,:])),dtype=float)
    mgiinds_var = n.zeros((len(nacat['binflux'][0,:])),dtype=float)
    mgiminds = n.zeros((len(nacat['binflux'][0,:])),dtype=float)
    tioinds = n.zeros((len(nacat['binflux'][0,:])),dtype=float)
    tioinds_var = n.zeros((len(nacat['binflux'][0,:])),dtype=float)
    tiominds = n.zeros((len(nacat['binflux'][0,:])),dtype=float)
    fehinds = n.zeros((len(nacat['binflux'][0,:])),dtype=float)
    fehinds_var = n.zeros((len(nacat['binflux'][0,:])),dtype=float)
    fehminds = n.zeros((len(nacat['binflux'][0,:])),dtype=float)
    fehosinds = n.zeros((len(nacat['binflux'][0,:])),dtype=float)
    fehosinds_var = n.zeros((len(nacat['binflux'][0,:])),dtype=float)
    fehosminds = n.zeros((len(nacat['binflux'][0,:])),dtype=float)
    fehososinds = n.zeros((len(nacat['binflux'][0,:])),dtype=float)
    fehososinds_var = n.zeros((len(nacat['binflux'][0,:])),dtype=float)
    fehososminds = n.zeros((len(nacat['binflux'][0,:])),dtype=float)
    fehindso = n.zeros((len(nacat['binflux'][0,:])),dtype=float)
    fehindso_var = n.zeros((len(nacat['binflux'][0,:])),dtype=float)
    fehmindso = n.zeros((len(nacat['binflux'][0,:])),dtype=float)

    if gspecs == 'True':
        naistart = n.where(lams > 8250.)[0][0]
        naiend = n.where(lams > 8500.)[0][0]
        nacatstart = n.where(lams > 8400.)[0][0]
        nacatend = n.where(lams > 9150.)[0][0]
        start = n.where(lams > 9700.)[0][0]
        end = n.where(lams > 10400.)[0][0]
        gz = nacat['gsol'][0]*1000./sc.c
        nacatgspec = ssz.z_correct(n.column_stack((lams[nacatstart:nacatend],nacat['skysubtdggal'][nacatstart:nacatend]/n.median(nacat['skysubtdggal'][nacatstart:nacatend]))), gz)
        nacatgvar = ssz.z_correct(n.column_stack((lams[nacatstart:nacatend],nacat['skysubtdggerr'][nacatstart:nacatend]**2/n.median(nacat['skysubtdggal'][nacatstart:nacatend]**2))), gz)
        nagspec = ssz.z_correct(n.column_stack((lams[naistart:naiend],na['skysubtdggal'][naistart:naiend]/n.median(na['skysubtdggal'][naistart:naiend]))), gz)
        nagvar = ssz.z_correct(n.column_stack((lams[naistart:naiend],na['skysubtdggerr'][naistart:naiend]**2/n.median(na['skysubtdggal'][naistart:naiend]**2))), gz)
        fehcgspec = ssz.z_correct(n.column_stack((lams[start:end],fehc['skysubtdggal'][start:end]/n.median(fehc['skysubtdggal'][start:end]))), gz)
        fehcgvar = ssz.z_correct(n.column_stack((lams[start:end],fehc['skysubtdggerr'][start:end]**2/n.median(fehc['skysubtdggal'][start:end]**2))), gz)
        fehcososgspec = ssz.z_correct(n.column_stack((lams[start:end],fehcosos['skysubtdggal'][start:end]/n.median(fehcosos['skysubtdggal'][start:end]))), gz)
        fehcososgvar = ssz.z_correct(n.column_stack((lams[start:end],fehcosos['skysubtdggerr'][start:end]**2/n.median(fehcosos['skysubtdggal'][start:end]**2))), gz)
        fehogspec = ssz.z_correct(n.column_stack((lams[start:end],feho['skysubtdggal'][start:end]/n.median(feho['skysubtdggal'][start:end]))), gz)
        fehogvar = ssz.z_correct(n.column_stack((lams[start:end],feho['skysubtdggerr'][start:end]**2/n.median(feho['skysubtdggal'][start:end]**2))), gz)
        if gal == 'GMP3329' and skymask == 'True':
            nacatgmspec = ssz.z_correct(n.column_stack((lams[nacatstart:nacatend],nacat['bestfitgMODEL'][nacatstart:nacatend]/n.median(nacat['skysubtdggal'][nacatstart:nacatend]))), gz)
            nacatgspec, nacatglinepos = mask_spec(nacatgspec, 'glob', skylines_cat, nacatgmspec)
            nacatgspec, nacatglinepos = mask_spec(nacatgspec, 'glob', skylines_mgi, nacatgmspec)
        if gal == 'GMP3367' and skymask == 'True':
            fehcgmspec = ssz.z_correct(n.column_stack((lams[start:end],fehc['bestfitgMODEL'][start:end]/n.median(fehc['skysubtdggal'][start:end]))), gz)
            fehcgspec, feholinepos = mask_spec(fehcgspec, 'glob', skylines_feh, fehcgmspec)
            fehcososgmspec = ssz.z_correct(n.column_stack((lams[start:end],fehcosos['bestfitgMODEL'][start:end]/n.median(fehcosos['skysubtdggal'][start:end]))), gz)
            fehcososgspec, feholinepos = mask_spec(fehcososgspec, 'glob', skylines_feh, fehcososgmspec)
            fehogmspec = ssz.z_correct(n.column_stack((lams[start:end],feho['bestfitgMODEL'][start:end]/n.median(feho['skysubtdggal'][start:end]))), gz)
            fehogspec, feholinepos = mask_spec(fehogspec, 'glob', skylines_feh, fehogmspec)
        nacatgspec = s.spectrum(lam=nacatgspec[:,0], lamspec=nacatgspec[:,1])
        nacatgvar = s.spectrum(lam=nacatgvar[:,0], lamspec=nacatgvar[:,1])
        nagspec = s.spectrum(lam=nagspec[:,0], lamspec=nagspec[:,1])
        nagvar = s.spectrum(lam=nagvar[:,0], lamspec=nagvar[:,1])
        fehcgspec = s.spectrum(lam=fehcgspec[:,0], lamspec=fehcgspec[:,1])
        fehcgvar = s.spectrum(lam=fehcgvar[:,0], lamspec=fehcgvar[:,1])
        fehcososgspec = s.spectrum(lam=fehcososgspec[:,0], lamspec=fehcososgspec[:,1])
        fehcososgvar = s.spectrum(lam=fehcososgvar[:,0], lamspec=fehcososgvar[:,1])
        fehogspec = s.spectrum(lam=fehogspec[:,0], lamspec=fehogspec[:,1])
        fehogvar = s.spectrum(lam=fehogvar[:,0], lamspec=fehogvar[:,1])

    for i in xrange(len(nacat['binflux'][0,:])):
        print i
        #redshift and velocity dispersion
        z = nacat['param'][0,i]*1000./sc.c
        sig = nacat['param'][1,i]
        print 'z = ', z
        print 'sig = ', sig, 'km/s'

        #Load spectra
        nacatspec = ssz.z_correct(n.column_stack((lams,nacat['skysubtdgal'][:,i])), z)#/n.median(nacat['skysubtdgal'][:,i]))), z)
        nacatmspec = ssz.z_correct(n.column_stack((lams,nacat['bestfitMODEL'][:,i])), z)#/n.median(nacat['bestfitMODEL'][:,i]))), z)
        nacatvar = ssz.z_correct(n.column_stack((lams,nacat['skysubtdgerr'][:,i]**2)), z)#/n.median(nacat['skysubtdgal'][:,i]**2))), z)

        naspec = ssz.z_correct(n.column_stack((lams,na['skysubtdgal'][:,i])), z)#/n.median(na['skysubtdgal'][:,i]))), z)
        namspec = ssz.z_correct(n.column_stack((lams,na['bestfitMODEL'][:,i])), z)#/n.median(na['bestfitMODEL'][:,i]))), z)
        navar = ssz.z_correct(n.column_stack((lams,na['skysubtdgerr'][:,i]**2)), z)#/n.median(na['skysubtdgal'][:,i]**2))), z)

        fehcspec = ssz.z_correct(n.column_stack((lams,fehc['skysubtdgal'][:,i])), z)#/n.median(fehc['skysubtdgal'][:,i]))), z)
        fehcmspec = ssz.z_correct(n.column_stack((lams,fehc['bestfitMODEL'][:,i])), z)#/n.median(fehc['bestfitMODEL'][:,i]))), z)
        fehcvar = ssz.z_correct(n.column_stack((lams,fehc['skysubtdgerr'][:,i]**2)), z)#/n.median(fehc['skysubtdgal'][:,i]**2))), z)

        # fehcosspec = ssz.z_correct(n.column_stack((lams,fehcos['skysubtdgal'][:,i])), z)#/n.median(fehc['skysubtdgal'][:,i]))), z)
        # fehcosmspec = ssz.z_correct(n.column_stack((lams,fehcos['bestfitMODEL'][:,i])), z)#/n.median(fehc['bestfitMODEL'][:,i]))), z)
        # fehcosvar = ssz.z_correct(n.column_stack((lams,fehcos['skysubtdgerr'][:,i]**2)), z)#/n.median(fehc['skysubtdgal'][:,i]**2))), z)
        fehcososspec = ssz.z_correct(n.column_stack((lams,fehcosos['skysubtdgal'][:,i])), z)#/n.median(fehc['skysubtdgal'][:,i]))), z)
        fehcososmspec = ssz.z_correct(n.column_stack((lams,fehcosos['bestfitMODEL'][:,i])), z)#/n.median(fehc['bestfitMODEL'][:,i]))), z)
        fehcososvar = ssz.z_correct(n.column_stack((lams,fehcosos['skysubtdgerr'][:,i]**2)), z)#/n.median(fehc['skysubtdgal'][:,i]**2))), z)
        fehospec = ssz.z_correct(n.column_stack((lams,feho['skysubtdgal'][:,i])), z)#/n.median(feho['skysubtdgal'][:,i]))), z)
        fehomspec = ssz.z_correct(n.column_stack((lams,feho['bestfitMODEL'][:,i])), z)#/n.median(feho['bestfitMODEL'][:,i]))), z)
        fehovar = ssz.z_correct(n.column_stack((lams,feho['skysubtdgerr'][:,i]**2)), z)#/n.median(feho['skysubtdgal'][:,i]**2))), z)

        #####Iterate over sky lines, mask, smooth, replace masked region with smoothed pixels x5
        if skymask=='True':
            print 'Finding and masking skylines'
            naspec, nalinepos = mask_spec(naspec, i, skylines_nai, namspec)
            nacatspec, nacatlinepos = mask_spec(nacatspec, i, skylines_cat, nacatmspec)
            nacatspec, nacatlinepos = mask_spec(nacatspec, i, skylines_mgi, nacatmspec)
            nacatspec, nacatlinepos = mask_spec(nacatspec, i, skylines_tio, nacatmspec)
            fehcspec, fehclinepos = mask_spec(fehcspec, i, skylines_feh, fehcmspec)
            # fehcosspec, fehcoslinepos = mask_spec(fehcosspec, i, skylines_feh, fehcosmspec)
            fehcososspec, fehcososlinepos = mask_spec(fehcososspec, i, skylines_feh, fehcososmspec)
            fehospec, feholinepos = mask_spec(fehospec, i, skylines_feh, fehomspec)

        #Load as spectrum class
        nacatspec = s.spectrum(lam=nacatspec[:,0], lamspec=nacatspec[:,1])
        nacatmspec = s.spectrum(lam=nacatmspec[:,0], lamspec=nacatmspec[:,1])
        nacatvar = s.spectrum(lam=nacatvar[:,0], lamspec=nacatvar[:,1])
        naspec = s.spectrum(lam=naspec[:,0], lamspec=naspec[:,1])
        namspec = s.spectrum(lam=namspec[:,0], lamspec=namspec[:,1])
        navar = s.spectrum(lam=navar[:,0], lamspec=navar[:,1])
        fehcspec = s.spectrum(lam=fehcspec[:,0], lamspec=fehcspec[:,1])
        fehcmspec = s.spectrum(lam=fehcmspec[:,0], lamspec=fehcmspec[:,1])
        fehcvar = s.spectrum(lam=fehcvar[:,0], lamspec=fehcvar[:,1])

        # fehcosspec = s.spectrum(lam=fehcosspec[:,0], lamspec=fehcosspec[:,1])
        # fehcosmspec = s.spectrum(lam=fehcosmspec[:,0], lamspec=fehcosmspec[:,1])
        # fehcosvar = s.spectrum(lam=fehcosvar[:,0], lamspec=fehcosvar[:,1])
        fehcososspec = s.spectrum(lam=fehcososspec[:,0], lamspec=fehcososspec[:,1])
        fehcososmspec = s.spectrum(lam=fehcososmspec[:,0], lamspec=fehcososmspec[:,1])
        fehcososvar = s.spectrum(lam=fehcososvar[:,0], lamspec=fehcososvar[:,1])
        fehospec = s.spectrum(lam=fehospec[:,0], lamspec=fehospec[:,1])
        fehomspec = s.spectrum(lam=fehomspec[:,0], lamspec=fehomspec[:,1])
        fehovar = s.spectrum(lam=fehovar[:,0], lamspec=fehovar[:,1])

        if out_sig != 'intrinsic':
            out_sig = float(out_sig)
            #Instrument sigma = 50.4 km/s for NaCaT; 63.3 km/s for FeH)
            nainstsig = 0.0#50.4
            fehinstsig = 0.0#63.3
            naconvsig = n.sqrt(out_sig**2 - (nacat['param'][1,i]**2+nainstsig**2))
            fehconvsig = n.sqrt(out_sig**2 - (nacat['param'][1,i]**2+fehinstsig**2))
            #Convolve to common vlocity dispersion
            loglambs = n.linspace(n.log10(naspec.lam[0]), n.log10(naspec.lam[-1]), len(naspec.lam))
            veldisp = (10**(loglambs[1]-loglambs[0]) - 1)*sc.c/1000. #in km/s
            nagauss = ssz.Gauss(naconvsig, veldisp)
            fehgauss = ssz.Gauss(fehconvsig, veldisp)

            nacatspec, nacatvar = convolve_spec(nacatspec, naconvsig, nagauss, nacatvar)
            naspec, navar = convolve_spec(naspec, naconvsig, nagauss, navar)
            fehcspec, fehcvar = convolve_spec(fehcspec, fehconvsig, fehgauss, fehcvar)
            # fehcosspec, fehcosvar = convolve_spec(fehcosspec, fehconvsig, fehgauss, fehcosvar)
            fehcososspec, fehcososvar = convolve_spec(fehcososspec, fehconvsig, fehgauss, fehcososvar)
            fehospec, fehovar = convolve_spec(fehospec, fehconvsig, fehgauss, fehovar)
        else:
            pass

        #Calculate indices
        naiind, naiindvar = naspec.irindex(naspec.lam[1]-naspec.lam[0], 'NaIsdss', navar, 'vac')
        naiinds[i] = naiind; naiinds_var[i] = naiindvar#+(0.01**2)
        naimind, naimindvar = namspec.irindex(namspec.lam[1]-namspec.lam[0], 'NaIsdss', None, 'vac')
        naiminds[i] = naimind

        catind, catindvar = nacatspec.irindex(nacatspec.lam[1]-nacatspec.lam[0], 'CaT', nacatvar, 'vac')
        catinds[i] = catind; catinds_var[i] = catindvar#+(0.01**2)
        catmind, catmindvar = nacatmspec.irindex(nacatmspec.lam[1]-nacatmspec.lam[0], 'CaT', None, 'vac')
        catminds[i] = catmind

        mgiind, mgiindvar = nacatspec.irindex(nacatspec.lam[1]-nacatspec.lam[0], 'MgI', nacatvar, 'vac')
        mgiinds[i] = mgiind; mgiinds_var[i] = mgiindvar#+(0.01**2)
        mgimind, mgimindvar = nacatmspec.irindex(nacatmspec.lam[1]-nacatmspec.lam[0], 'MgI', None, 'vac')
        mgiminds[i] = mgimind

        tioind, tioindvar = nacatspec.TiOindex(nacatspec.lam[1]-nacatspec.lam[0], nacatvar, 'vac')
        tioinds[i] = tioind; tioinds_var[i] = tioindvar#+(0.01**2)
        tiomind, tiomindvar = nacatmspec.TiOindex(nacatmspec.lam[1]-nacatmspec.lam[0], None, 'vac')
        tiominds[i] = tiomind

        fehind, fehindvar = fehcspec.irindex(fehcspec.lam[1]-fehcspec.lam[0], 'FeH', fehcvar, 'vac')
        fehinds[i] = fehind; fehinds_var[i] = fehindvar
        fehmind, fehmindvar = fehcmspec.irindex(fehcmspec.lam[1]-fehcmspec.lam[0], 'FeH', None, 'vac')
        fehminds[i] = fehmind

        # fehind, fehindvar = fehcosspec.irindex(fehcosspec.lam[1]-fehcosspec.lam[0], 'FeH', fehcosvar, 'vac')
        # fehosinds[i] = fehind; fehosinds_var[i] = fehindvar
        # fehmind, fehmindvar = fehcosmspec.irindex(fehcosmspec.lam[1]-fehcosmspec.lam[0], 'FeH', None, 'vac')
        # fehosminds[i] = fehmind
        fehind, fehindvar = fehcososspec.irindex(fehcososspec.lam[1]-fehcososspec.lam[0], 'FeH', fehcososvar, 'vac')
        fehososinds[i] = fehind; fehososinds_var[i] = fehindvar
        fehmind, fehmindvar = fehcososmspec.irindex(fehcososmspec.lam[1]-fehcososmspec.lam[0], 'FeH', None, 'vac')
        fehososminds[i] = fehmind
        fehind, fehindvar = fehospec.irindex(fehospec.lam[1]-fehospec.lam[0], 'FeH', fehovar, 'vac')
        fehindso[i] = fehind; fehindso_var[i] = fehindvar
        fehmind, fehmindvar = fehomspec.irindex(fehomspec.lam[1]-fehomspec.lam[0], 'FeH', None, 'vac')
        fehmindso[i] = fehmind

    if gspecs == 'True':
        if out_sig != 'intrinsic':
            naconvgsig = n.sqrt(out_sig**2 - (nacat['gsol'][1]**2+nainstsig**2))
            fehconvgsig = n.sqrt(out_sig**2 - (nacat['gsol'][1]**2+fehinstsig**2))
            fehloglambs = n.linspace(n.log10(fehcgspec.lam[0]), n.log10(fehcgspec.lam[-1]), len(fehcgspec.lam))
            fehveldisp = (10**(fehloglambs[1]-fehloglambs[0]) - 1)*sc.c/1000. #in km/s
            naggauss = ssz.Gauss(naconvgsig, veldisp)
            fehggauss = ssz.Gauss(fehconvgsig, fehveldisp)
            nacatgspec, nacatgvar = convolve_spec(nacatgspec, naconvgsig, naggauss, nacatgvar)
            nagspec, nagvar = convolve_spec(nagspec, naconvgsig, naggauss, nagvar)
            fehcgspec, fehcgvar = convolve_spec(fehcgspec, fehconvgsig, fehggauss, fehcgvar)
            fehcososgspec, fehcososgvar = convolve_spec(fehcososgspec, fehconvgsig, fehggauss, fehcososgvar)
            fehogspec, fehogvar = convolve_spec(fehogspec, fehconvgsig, fehggauss, fehogvar)
        else:
            pass

        # fcont = nacatgspec.fehcont
        # feat = nacatgspec.feh
        # print 'FEAT = ', feat
        # print 'FCONT = ', fcont
        # indbmid = (fcont[0]+fcont[1])/2.
        # indrmid = (fcont[-1]+fcont[-1])/2.
        # fehcgspec, fehcgvar = normalise_spec(fehcgspec, fcont, [indbmid, indrmid], fehcgvar)
        # fehcososgspec, fehcososgvar = normalise_spec(fehcososgspec, fcont, [indbmid, indrmid], fehcososgvar)
        # fehogspec, fehogvar = normalise_spec(fehogspec, fcont, [indbmid, indrmid], fehogvar)
        gnai, gnaivar = nagspec.irindex(nagspec.lam[1]-nagspec.lam[0], 'NaIsdss', nagvar, 'vac')
        gcat, gcatvar = nacatgspec.irindex(nacatgspec.lam[1]-nacatgspec.lam[0], 'CaT', nacatgvar, 'vac')
        gmgi, gmgivar = nacatgspec.irindex(nacatgspec.lam[1]-nacatgspec.lam[0], 'MgI', nacatgvar, 'vac')
        gtio, gtiovar = nacatgspec.TiOindex(nacatgspec.lam[1]-nacatgspec.lam[0], nacatgvar, 'vac')
        gfeh, gfehvar = fehcgspec.irindex(fehcgspec.lam[1]-fehcgspec.lam[0], 'FeH', fehcgvar, 'vac')
        gfehosos, gfehososvar = fehcososgspec.irindex(fehcososgspec.lam[1]-fehcososgspec.lam[0], 'FeH', fehcososgvar, 'vac')
        gfeho, gfehovar = fehogspec.irindex(fehogspec.lam[1]-fehogspec.lam[0], 'FeH', fehogvar, 'vac')

    if out_sig == 'intrinsic' and corrfac == 'True':
        #Use correction factor to scale indices to central values
        if corrval == 'central':
            corrsig = nacat['param'][1,0]
        else:
            corrsig = n.float(corrval)
        naisyserr = n.zeros(len(naiinds))
        catsyserr = n.zeros_like(naisyserr)
        mgisyserr = n.zeros_like(naisyserr)
        tiosyserr = n.zeros_like(naisyserr)
        fehsyserr = n.zeros_like(naisyserr)
        for i in xrange(len(nacat['binflux'][0,:])):
            naifac, naisys = CvD12_indexvsig('NaIsdss', nacat['param'][1,i], corrsig, False, False)
            catfac, catsys = CvD12_indexvsig('CaT', nacat['param'][1,i], corrsig, False, False)
            mgifac, mgisys = CvD12_indexvsig('MgI', nacat['param'][1,i], corrsig, False, False)
            tiofac, tiosys = CvD12_indexvsig('TiO', nacat['param'][1,i], corrsig, False, False)
            fehfac, fehsys = CvD12_indexvsig('FeH', nacat['param'][1,i], corrsig, False, False)
            naiinds[i] /= naifac; naisyserr[i] = naisys
            catinds[i] /= catfac; catsyserr[i] = catsys
            mgiinds[i] /= mgifac; mgisyserr[i] = mgisys
            tioinds[i] /= tiofac; tiosyserr[i] = tiosys
            fehinds[i] /= fehfac; fehsyserr[i] = fehsys
            fehososinds[i] /= fehfac
            fehindso[i] /= fehfac
        if gspecs == 'True':
            gnaifac, gnaisys = CvD12_indexvsig('NaIsdss', nacat['gsol'][1], corrsig, False, False)
            gcatfac, gcatsys = CvD12_indexvsig('CaT', nacat['gsol'][1], corrsig, False, False)
            gmgifac, gmgisys = CvD12_indexvsig('MgI', nacat['gsol'][1], corrsig, False, False)
            gtiofac, gtiosys = CvD12_indexvsig('TiO', nacat['gsol'][1], corrsig, False, False)
            gfehfac, gfehsys = CvD12_indexvsig('FeH', nacat['gsol'][1], corrsig, False, False)
            gnai /= gnaifac
            gcat /= gcatfac
            gmgi /= gmgifac
            gtio /= gtiofac
            gfeh /= gfehfac
            gfehosos /= gfehfac
            gfeho /= gfehfac



    if plot == 'True':
        #Generate plot
        P.rcParams['font.family'] = 'Times New Roman'
        P.rcParams.update({'font.size': 20})
        f, (ax1, ax2, ax3, ax4, ax5) = P.subplots(5, sharex=True, figsize=(16,24))
        #NaI
        p1a = ax1.errorbar(rs, naiinds, yerr=n.sqrt(naiinds_var), fmt='ko', ms=msize,
                     capthick=2, elinewidth=1.5, alpha=1)
        # ax1.plot(rs, naiminds, c='r', marker='s', ms=msize-2, linestyle='None', alpha=1)
        #CaT*
        ax2.errorbar(rs, catinds, yerr=n.sqrt(catinds_var), fmt='ko', ms=msize,
                     capthick=2, elinewidth=1.5, alpha=1)
        # ax2.plot(rs, catminds, c='r', marker='s', ms=msize-2, linestyle='None', alpha=1)
        #MgI
        ax3.errorbar(rs, mgiinds, yerr=n.sqrt(mgiinds_var), fmt='ko', ms=msize,
                     capthick=2, elinewidth=1.5, alpha=1)
        # ax3.plot(rs, mgiminds, c='r', marker='s', ms=msize-2, linestyle='None', alpha=1)
        #TiO
        ax4.errorbar(rs, tioinds, yerr=n.sqrt(tioinds_var), fmt='ko', ms=msize,
                     capthick=2, elinewidth=1.5, alpha=1)
        # ax4.plot(rs, tiominds, c='r', marker='s', ms=msize-2, linestyle='None', alpha=1)
        #FeH
        ax5.errorbar(rs, fehinds, yerr=n.sqrt(fehinds_var), fmt='ko', ms=msize,
                     capthick=2, elinewidth=1.5, alpha=1, label='Fit Div.')
        if fehopt == 'True':
            ax5.errorbar(rs, fehindso, yerr=n.sqrt(fehindso_var), fmt='r*', ms=msize+5,
                        capthick=2, elinewidth=1.5, alpha=1, label='Tell. Corr.')
            if gal == 'GMP2921':
                ax5.legend(loc='upper right', frameon=False, borderaxespad=0.2, labelspacing=0.18, numpoints=1)
        if osfitted == 'True':
            ax5.errorbar(rs, fehososinds, yerr=n.sqrt(fehososinds_var), fmt='bs', ms=msize,
                        capthick=2, elinewidth=1.5, alpha=1, label='O-S fitted')
            if gal == 'GMP2921':
                ax5.legend(loc='upper right', frameon=False, borderaxespad=0.2, labelspacing=0.18, numpoints=1)

        if gspecs == 'True':
            if gal in ['GMP2921', 'GMP3329', 'GMP4928']:
                grval = 4.0
            elif gal == 'GMP3367':
                grval = 4.0
                # handles.append(P.Line2D((),(), color=colours[gals.index('3367')], marker=markers[gals.index('3367')], linestyle='', markersize=msize))
            fc = 'w'
            #NaI
            ax1.errorbar(grval, gnai, yerr=n.sqrt(gnaivar), fmt='o', ms=msize+1, mec='k',
            mew=1.3, mfc=fc, ls='', ecolor='k', capthick=2, elinewidth=1.5, zorder=10)
            #CaT*
            ax2.errorbar(grval, gcat, yerr=n.sqrt(gcatvar), fmt='o', ms=msize+1, mec='k',
            mew=1.3, mfc=fc, ls='', ecolor='k', capthick=2, elinewidth=1.5, zorder=10)
            #MgI
            ax3.errorbar(grval, gmgi, yerr=n.sqrt(gmgivar), fmt='o', ms=msize+1, mec='k',
            mew=1.3, mfc=fc, ls='', ecolor='k', capthick=2, elinewidth=1.5, zorder=10)
            #TiO
            ax4.errorbar(grval, gtio, yerr=n.sqrt(gtiovar), fmt='o', ms=msize+1, mec='k',
            mew=1.3, mfc=fc, ls='', ecolor='k', capthick=2, elinewidth=1.5, zorder=10)
            #Fe
            ax5.errorbar(grval, gfeh, yerr=n.sqrt(gfehvar), fmt='o', ms=msize+1, mec='k',
            mew=1.3, mfc=fc, ls='', ecolor='k', capthick=2, elinewidth=1.5, zorder=10)
            if fehopt == 'True':
                ax5.errorbar(grval, gfeho, yerr=n.sqrt(gfehovar), fmt='*', ms=msize+5, mec='r',
                            capthick=2, mfc=fc, mew=1.3, ecolor='r', elinewidth=1.5, alpha=1)
            if osfitted == 'True':
                ax5.errorbar(grval, gfehosos, yerr=n.sqrt(gfehososvar), fmt='s', ms=msize, mec='b',
                            capthick=2, mfc=fc, mew=1.3, ecolor='b', elinewidth=1.5, alpha=1)

        if gal=='GMP4928':
            ax5.set_ylim([0.0, 1.0])
        else:
            ax5.set_ylim([0.0,0.8])

        # Fine-tune figure; make subplots close to each other and hide x ticks for all but bottom plot.
        f.subplots_adjust(hspace=0.1)
        #P.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
        P.xlabel('r [arcsec]', fontsize=tfsize+2)
        ax1.set_ylabel('NaI$_{SDSS}$ [$\AA$]', fontsize=tfsize+2)
        ax2.set_ylabel('CaT [$\AA$]', fontsize=tfsize+2)
        ax3.set_ylabel('MgI0.88 [$\AA$]', fontsize=tfsize+2)
        ax4.set_ylabel('TiO0.89', fontsize=tfsize+2)
        ax5.set_ylabel('FeH0.99 [$\AA$]', fontsize=tfsize+2)
        ax1.yaxis.label.set_fontsize(lfsize+2)
        ax1.yaxis.label.set_fontsize(lfsize+2)
        ax1.tick_params(axis='both', which='major', direction='in', length=majtsize, width=twidth)
        ax1.tick_params(axis='both', which='minor', direction='in', length=mintsize, width=twidth)
        ax2.tick_params(axis='both', which='major', direction='in', length=majtsize, width=twidth)
        ax2.tick_params(axis='both', which='minor', direction='in', length=mintsize, width=twidth)
        ax3.tick_params(axis='both', which='major', direction='in', length=majtsize, width=twidth)
        ax3.tick_params(axis='both', which='minor', direction='in', length=mintsize, width=twidth)
        ax4.tick_params(axis='both', which='major', direction='in', length=majtsize, width=twidth)
        ax4.tick_params(axis='both', which='minor', direction='in', length=mintsize, width=twidth)
        ax5.tick_params(axis='both', which='major', direction='in', length=majtsize, width=twidth)
        ax5.tick_params(axis='both', which='minor', direction='in', length=mintsize, width=twidth)
        P.minorticks_on()
    #    leg = ax5.legend(['FeH_cont', 'FeH_opt'], loc='upper right',
    #                frameon=False, borderaxespad=0.5, labelspacing=0.5, numpoints=1)
    #    ltext = leg.get_texts()
    #    P.setp(ltext[0], fontsize = lfsize, color = 'k')
    #    P.setp(ltext[1], fontsize = lfsize, color = 'r')

        P.xlim([-0.5,10.])
        P.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
        f.tight_layout()
        #P.savefig('/Users/zieleniewski/Data/swift/Coma/Coma_stellar_pops/'+gal+'/Indices_v_radius_'+gal+'.png')
        if corrfac=='False':
            P.savefig('/Users/zieleniewski/Data/swift/Coma_stellar_pops/'+gal+'/Indices_v_radius_'+gal+'_'+str(out_sig)[:-2]+'kms_data_'+str(time.gmtime()[0])+'-'+str(time.gmtime()[1])+'-'+str(time.gmtime()[2])+'.png')
        else:
            P.savefig('/Users/zieleniewski/Data/swift/Coma_stellar_pops/'+gal+'/Indices_v_radius_'+gal+'_'+str(corrsig)[:-2]+'kms_data_corrfac_'+str(time.gmtime()[0])+'-'+str(time.gmtime()[1])+'-'+str(time.gmtime()[2])+'.png')
        P.show()

    if save == 'True':
        if corrfac=='False':
            sfilename = '/Users/zieleniewski/Data/swift/Coma_stellar_pops/'+gal+\
                      '/Indices_v_radius_'+gal+'_'+str(out_sig)[:-2]+'kms_data_'+\
                      str(time.gmtime()[0])+'-'+str(time.gmtime()[1])+'-'+str(time.gmtime()[2])+'.txt'
        else:
            sfilename = '/Users/zieleniewski/Data/swift/Coma_stellar_pops/'+gal+\
                      '/Indices_v_radius_'+gal+'_'+str(corrsig)[:-2]+'kms_data_corrfac_'+\
                      str(time.gmtime()[0])+'-'+str(time.gmtime()[1])+'-'+str(time.gmtime()[2])+'.txt'
        n.savetxt(sfilename,
                  n.column_stack((rs, naiinds, naiinds_var, naiminds,
                                  catinds, catinds_var, catminds,
                                  mgiinds, mgiinds_var, mgiminds,
                                  tioinds, tioinds_var, tiominds,
                                  fehinds, fehinds_var, fehminds,
                                  #fehosinds, fehosinds_var, fehosminds,
                                  fehososinds, fehososinds_var, fehososminds,
                                  fehindso, fehindso_var, fehmindso)),
                  delimiter=',', fmt='%.7f',
                  header='r [arcsec], NaI [A], NaIvar, NaImodel, CaT [A], CaTvar, CaTmodel, MgI [A], MgIvar, MgImodel, TiO, TiOvar, TiOmodel, FeH [A], FeHvar, FeHmodel, FeH O-S [A], FeHvar O-S, FeHmodel O-S, FeH opt [A], FeHvar opt , FeHmodel opt')
    if save == 'True' and gspecs == 'True':
        if corrfac=='False':
            gsfilename = '/Users/zieleniewski/Data/swift/Coma_stellar_pops/'+gal+\
                      '/'+gal+'_global_spec_index_values_'+str(out_sig)[:-2]+'kms_data_'+\
                      str(time.gmtime()[0])+'-'+str(time.gmtime()[1])+'-'+str(time.gmtime()[2])+'.txt'
        else:
            gsfilename = '/Users/zieleniewski/Data/swift/Coma_stellar_pops/'+gal+\
                      '/'+gal+'_global_spec_index_values_'+str(corrsig)[:-2]+'kms_data_corrfac_'+\
                      str(time.gmtime()[0])+'-'+str(time.gmtime()[1])+'-'+str(time.gmtime()[2])+'.txt'
        n.savetxt(gsfilename,
                  n.column_stack((gnai, n.sqrt(gnaivar),
                                  gcat, n.sqrt(gcatvar),
                                  gmgi, n.sqrt(gmgivar),
                                  gtio, n.sqrt(gtiovar),
                                  gfeh, n.sqrt(gfehvar),
                                  gfehosos, n.sqrt(gfehososvar),
                                  gfeho, n.sqrt(gfehovar))),
                  delimiter=',', fmt='%.3f',
                  header='\n'.join(['#'+gal+' global index values at '+str(out_sig)[:-2]+' km/s',
                  '#NaIsdss [A], NaIsdss_err, CaT [A], CaT_err, MgI [A], MgI_err, TiO, TiO_err, FeH [A], FeH_err, FeH (O-S), FeH (O-S) err, FeH (opt), FeH (opt) err']))
#------------#




#------------#
def indices_v_radius_allgals(dir2921='None', dir3329='None', dir4928='None', dir3367='None', save='False', model='False', gvals=False):

    #GMP2921
    gals = []
    galsdata = []
    sigs = []
    if dir2921 != 'None':
        data2921 = n.genfromtxt(dir2921, delimiter=',')
        gals.append('2921')
        galsdata.append(data2921)
        sigs.append(dir2921.split('kms')[0][-3:])
    if dir3329 != 'None':
        data3329 = n.genfromtxt(dir3329, delimiter=',')
        gals.append('3329')
        galsdata.append(data3329)
        sigs.append(dir3329.split('kms')[0][-3:])
    if dir4928 != 'None':
        data4928 = n.genfromtxt(dir4928, delimiter=',')
        gals.append('4928')
        galsdata.append(data4928)
        sigs.append(dir4928.split('kms')[0][-3:])
    if dir3367 != 'None':
        data3367 = n.genfromtxt(dir3367, delimiter=',')
        gals.append('3367')
        # galsdata.append(data3367)
        sigs.append(dir3367.split('kms')[0][-3:])

    lgals = ['NGC4889', 'NGC4874', 'NGC4839', 'NGC4873']
    print 'GALS = ', gals
    print 'GALS = ', lgals

    handles = []
    #Generate plot
    P.rcParams['font.family'] = 'Times New Roman'
    P.rcParams.update({'font.size': 20})
    f, (ax1, ax2, ax3, ax4, ax5) = P.subplots(5, sharex=True, figsize=(16,24))
    markers = ['o', 's', '^', 'D']
    colours = ['k', 'r', 'b', 'g']
    offsets = [-0.1,-0.05,0.05,0.1]

    for i in xrange(len(galsdata)):
        galsdata[i][0,0] += offsets[i]
        #NaI
        if gals[i] == '4928':#Outermost point of GMP4928 affected by sky
            p1a = ax1.errorbar(galsdata[i][:-1,0], galsdata[i][:-1,1], yerr=n.sqrt(galsdata[i][:-1,2]), fmt=colours[i]+markers[i], ms=msize,
                         capthick=2, elinewidth=1.5, alpha=1)
            p1b = ax1.errorbar(galsdata[i][-1,0], galsdata[i][-1,1], yerr=n.sqrt(galsdata[i][-1,2]), fmt=colours[i]+markers[i], ms=msize,
                         capthick=2, elinewidth=1.5, alpha=1)
            p1b[-1][0].set_linestyle('--')
        elif gals[i] == '3329':#Outermost point of GMP3329 affected by sky
            p1a = ax1.errorbar(galsdata[i][:-1,0], galsdata[i][:-1,1], yerr=n.sqrt(galsdata[i][:-1,2]), fmt=colours[i]+markers[i], ms=msize,
                         capthick=2, elinewidth=1.5, alpha=1)
            p1b = ax1.errorbar(galsdata[i][-1,0], galsdata[i][-1,1], yerr=n.sqrt(galsdata[i][-1,2]), fmt=colours[i]+markers[i], ms=msize,
                         capthick=2, elinewidth=1.5, alpha=1)
            p1b[-1][0].set_linestyle('--')
        else:
            p1a = ax1.errorbar(galsdata[i][:,0], galsdata[i][:,1], yerr=n.sqrt(galsdata[i][:,2]), fmt=colours[i]+markers[i], ms=msize,
                         capthick=2, elinewidth=1.5, alpha=1)
        #CaT
        if gals[i] == '4928':#Outermost point of GMP4928 affected by sky
            p2a = ax2.errorbar(galsdata[i][:-1,0], galsdata[i][:-1,4], yerr=n.sqrt(galsdata[i][:-1,5]), fmt=colours[i]+markers[i], ms=msize,
                         capthick=2, elinewidth=1.5, alpha=1)
            p2b = ax2.errorbar(galsdata[i][-1,0], galsdata[i][-1,4], yerr=n.sqrt(galsdata[i][-1,5]), fmt=colours[i]+markers[i], ms=msize,
                         capthick=2, elinewidth=1.5, alpha=1)
            p2b[-1][0].set_linestyle('--')
        else:
            p2a = ax2.errorbar(galsdata[i][:,0], galsdata[i][:,4], yerr=n.sqrt(galsdata[i][:,5]), fmt=colours[i]+markers[i], ms=msize,
                         capthick=2, elinewidth=1.5, alpha=1)
        #MgI
        p3a = ax3.errorbar(galsdata[i][:,0], galsdata[i][:,7], yerr=n.sqrt(galsdata[i][:,8]), fmt=colours[i]+markers[i], ms=msize,
                     capthick=2, elinewidth=1.5, alpha=1)
        #TiO
        p4a = ax4.errorbar(galsdata[i][:,0], galsdata[i][:,10], yerr=n.sqrt(galsdata[i][:,11]), fmt=colours[i]+markers[i], ms=msize,
                     capthick=2, elinewidth=1.5, alpha=1)
        #FeH
        if gals[i] == '3329':#Second from last point of GMP3329 affected by sky
            p5a = ax5.errorbar(galsdata[i][:-2,0], galsdata[i][:-2,13], yerr=n.sqrt(galsdata[i][:-2,14]), fmt=colours[i]+markers[i], ms=msize,
                         capthick=2, elinewidth=1.5, alpha=1)
            p5c = ax5.errorbar(galsdata[i][-1,0], galsdata[i][-1,13], yerr=n.sqrt(galsdata[i][-1,14]), fmt=colours[i]+markers[i], ms=msize,
                         capthick=2, elinewidth=1.5, alpha=1)
            p5b = ax5.errorbar(galsdata[i][-2,0], galsdata[i][-2,13], yerr=n.sqrt(galsdata[i][-2,14]), fmt=colours[i]+markers[i], ms=msize,
                         capthick=2, elinewidth=1.5, alpha=1)
            p5b[-1][0].set_linestyle('--')
        else:
            p5a = ax5.errorbar(galsdata[i][:,0], galsdata[i][:,13], yerr=n.sqrt(galsdata[i][:,14]), fmt=colours[i]+markers[i], ms=msize,
                         capthick=2, elinewidth=1.5, alpha=1)
        if model == 'True':
            plb = ax1.plot(galsdata[i][:,0], galsdata[i][:,3], c=colours[i], marker='x', ms=msize, linestyle='None', alpha=1)
            ax2.plot(galsdata[i][:,0], galsdata[i][:,6], c=colours[i], marker='x', ms=msize, linestyle='None', alpha=1)
            ax3.plot(galsdata[i][:,0], galsdata[i][:,9], c=colours[i], marker='x', ms=msize, linestyle='None', alpha=1)
            ax4.plot(galsdata[i][:,0], galsdata[i][:,12], c=colours[i], marker='x', ms=msize, linestyle='None', alpha=1)
            ax5.plot(galsdata[i][:,0], galsdata[i][:,15], c=colours[i], marker='x', ms=msize, linestyle='None', alpha=1)
        handles.append(P.Line2D((),(), color=colours[i], marker=markers[i], linestyle='', markersize=msize))
    ax1.set_ylim([0.2, 1.1])
    ax4.set_ylim([1.045,1.085])
    ax5.set_ylim([0.,1.0])

    gdat = []
    if gvals:
        for i in xrange(len(gvals)):
            dat = n.genfromtxt(gvals[i], delimiter=',')
            gdat.append(dat)
            print gals[i]
            grval = [3.9, 3.95, 4.05, 4.1]
            if gals[i] in ['2921', '3329', '4928']:
                fc = 'w'
                zord = 10.0
            if gals[i] == '3367':
                fc = 'w'
                zord= 5.0
                handles.append(P.Line2D((),(), color=colours[gals.index('3367')], marker=markers[gals.index('3367')], linestyle='', markersize=msize))
            #NaI
            ax1.errorbar(grval[i], dat[0], yerr=dat[1], fmt=markers[i], ms=msize+1, mec=colours[i],
            mew=1.3, mfc=fc, ls='', ecolor=colours[i], capthick=2, elinewidth=1.5, zorder=zord)
            #CaT*
            ax2.errorbar(grval[i], dat[2], yerr=dat[3], fmt=markers[i], ms=msize+1, mec=colours[i],
            mew=1.3, mfc=fc, ls='', ecolor=colours[i], capthick=2, elinewidth=1.5, zorder=zord)
            #MgI
            ax3.errorbar(grval[i], dat[4], yerr=dat[5], fmt=markers[i], ms=msize+1, mec=colours[i],
            mew=1.3, mfc=fc, ls='', ecolor=colours[i], capthick=2, elinewidth=1.5, zorder=zord)
            #TiO
            ax4.errorbar(grval[i], dat[6], yerr=dat[7], fmt=markers[i], ms=msize+1, mec=colours[i],
            mew=1.3, mfc=fc, ls='', ecolor=colours[i], capthick=2, elinewidth=1.5, zorder=zord)
            #Fe
            ax5.errorbar(grval[i], dat[8], yerr=dat[9], fmt=markers[i], ms=msize+1, mec=colours[i],
            mew=1.3, mfc=fc, ls='', ecolor=colours[i], capthick=2, elinewidth=1.5, zorder=zord)

    # Fine-tune figure; make subplots close to each other and hide x ticks for all but bottom plot.
    f.subplots_adjust(hspace=0.08)
    P.xlabel('r [arcsec]', fontsize=tfsize+2)
    ax1.minorticks_on()
    ax2.minorticks_on()
    ax3.minorticks_on()
    ax4.minorticks_on()
    ax5.minorticks_on()
    ax1.set_ylabel('NaI$_{\\rm{SDSS}}$0.82 [$\AA$]', fontsize=tfsize+2)
    ax2.set_ylabel('CaT [$\AA$]', fontsize=tfsize+2)
    ax3.set_ylabel('MgI0.88 [$\AA$]', fontsize=tfsize+2)
    ax4.set_ylabel('TiO0.89', fontsize=tfsize+2)
    ax5.set_ylabel('FeH0.99 [$\AA$]', fontsize=tfsize+2)
    ax1.yaxis.label.set_fontsize(tfsize+2)
    ax2.yaxis.label.set_fontsize(tfsize+2)
    ax3.yaxis.label.set_fontsize(tfsize+2)
    ax4.yaxis.label.set_fontsize(tfsize+2)
    ax5.yaxis.label.set_fontsize(tfsize+2)
    ax5.xaxis.label.set_fontsize(tfsize+2)
    ax1.tick_params(axis='both', which='major', direction='in', length=majtsize, width=twidth)
    ax1.tick_params(axis='both', which='minor', direction='in', length=mintsize, width=twidth)
    ax2.tick_params(axis='both', which='major', direction='in', length=majtsize, width=twidth)
    ax2.tick_params(axis='both', which='minor', direction='in', length=mintsize, width=twidth)
    ax3.tick_params(axis='both', which='major', direction='in', length=majtsize, width=twidth)
    ax3.tick_params(axis='both', which='minor', direction='in', length=mintsize, width=twidth)
    ax4.tick_params(axis='both', which='major', direction='in', length=majtsize, width=twidth)
    ax4.tick_params(axis='both', which='minor', direction='in', length=mintsize, width=twidth)
    ax5.tick_params(axis='both', which='major', direction='in', length=majtsize, width=twidth)
    ax5.tick_params(axis='both', which='minor', direction='in', length=mintsize, width=twidth)
    ax1.set_xlim([-0.5,10.])

    leg = ax2.legend(handles, [lgals[i]+' ('+sigs[i]+' km s$^{-1}$)' for i in xrange(len(gals))], loc='upper right',
                frameon=False, borderaxespad=0.5, labelspacing=0.2, numpoints=1)
    ltext = leg.get_texts()
    for j in xrange(len(ltext)):
        P.setp(ltext[j], fontsize=20, color=colours[j])

    ###Add kpc x-axis along top
    ax11 = ax1.twiny()
    new_tick_locations = n.array([0, 2, 4, 6, 8, 10])
    kpc_per_arcsec = NWCC.calc(0.024, 68.0, 0.3, 0.7, 'general')[2]
    new_ticks = n.round(n.array([0,2,4,6,8,10])*kpc_per_arcsec, 0)
    ax11.set_xlim(ax1.get_xlim())
    ax11.set_xticks(new_tick_locations)
    ax11.set_xticklabels(new_ticks.astype(n.int))
    ax11.set_xlabel('r [kpc]', fontsize=tfsize+2)
    ax11.xaxis.label.set_fontsize(tfsize+2)
    ax11.tick_params(axis='both', which='major', direction='in', length=majtsize, width=twidth)
    ax11.tick_params(axis='both', which='minor', direction='in', length=mintsize, width=twidth)
    ax11.xaxis.set_minor_locator(AutoMinorLocator())
    ###

    # ###
    # #Add min and max systematic error bars!
    # dat2921 = n.genfromtxt('/Users/zieleniewski/Data/swift/Coma_stellar_pops/GMP2921_FeH_min_max_dat.txt', delimiter=',')
    # dat3329 = n.genfromtxt('/Users/zieleniewski/Data/swift/Coma_stellar_pops/GMP3329_FeH_min_max_dat.txt', delimiter=',')
    # dat4928 = n.genfromtxt('/Users/zieleniewski/Data/swift/Coma_stellar_pops/GMP4928_FeH_min_max_dat.txt', delimiter=',')
    # dats = [dat2921, dat3329, dat4928]
    # for i in xrange(len(dats)):
    #     print dats[i][:,0]
    #     print galsdata[i][:,13]
    #     print [gdat[i][8]]
    #     alldats = n.concatenate((galsdata[i][:,13],[gdat[i][8]]))
    #     print n.row_stack((alldats-dats[i][:,1],dats[i][:,2]-alldats))

    #     dats[i][0,0] += offsets[i]
    #     dats[i][-1,0] += offsets[i]
    #     plsys = ax5.errorbar(dats[i][:,0], alldats,
    #     yerr=n.row_stack((alldats-dats[i][:,1],dats[i][:,2]-alldats)),
    #     marker='', ls='', ecolor=colours[i], capthick=2, elinewidth=2.5, zorder=zord)
    #     plsys[-1][0].set_linestyle('--')  
    # ###

    ###
    #Add min and max systematic range as rectangle
    dat2921 = n.genfromtxt('/Users/zieleniewski/Data/swift/Coma_stellar_pops/GMP2921_FeH_min_max_dat.txt', delimiter=',')
    dat3329 = n.genfromtxt('/Users/zieleniewski/Data/swift/Coma_stellar_pops/GMP3329_FeH_min_max_dat.txt', delimiter=',')
    dat4928 = n.genfromtxt('/Users/zieleniewski/Data/swift/Coma_stellar_pops/GMP4928_FeH_min_max_dat.txt', delimiter=',')
    dats = [dat2921, dat3329, dat4928]
    #
    def makeErrorBoxes(xdata,ydata,xerror,yerror,fc='r',ec='None',alpha=0.25):

        # Create list for all the error patches
        errorboxes = []

        # Loop over data points; create box from errors at each point
        for xc,yc,xe,ye in zip(xdata,ydata,xerror.T,yerror.T):
            print xc, yc, xe, ye
            rect = Rectangle((xc-xe[0],yc-ye[0]),xe.sum(),ye.sum())
            errorboxes.append(rect)

        # Create patch collection with specified colour/alpha
        pc = PatchCollection(errorboxes,facecolor=fc,alpha=alpha,edgecolor=ec)

        # Add collection to axes
        ax5.add_collection(pc)
    #
    for i in xrange(len(dats)):
        alldats = n.concatenate((galsdata[i][:,13],[gdat[i][8]]))
        dats[i][0,0] += offsets[i]
        dats[i][-1,0] += offsets[i]
        makeErrorBoxes(dats[i][:,0], alldats, xerror=n.row_stack((n.array([0.05]*len(dats[i][:,0])),n.array([0.05]*len(dats[i][:,0])))),
        yerror=n.row_stack((alldats-dats[i][:,1],dats[i][:,2]-alldats)), fc=colours[i], ec=colours[i])
    ###

    # P.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
    f.tight_layout()
    f.subplots_adjust(hspace=0.08)
    if save != 'False':
        P.savefig('/Users/zieleniewski/Data/swift/Coma_stellar_pops/'+str(save)+'.png')
    P.show()
#------------#



#------------#
def indices_v_radius_compare_indices(file1, gfile1, file2, gfile2, save='False'):

    #GMP2921
    gals = []
    galsdata = []
    ggalsdata = []
    galsdatacf = []
    ggalsdatacf = []
    sigs = []

    data1 = n.genfromtxt(file1, delimiter=',')
    gdata1 = n.genfromtxt(gfile1, delimiter=',')
    gals.append(file1.split('/')[-1].split('_')[3])
    galsdata.append(data1)
    ggalsdata.append(gdata1)
    sigs.append(file1.split('kms')[0][-3:])

    data2 = n.genfromtxt(file2, delimiter=',')
    gdata2 = n.genfromtxt(gfile2, delimiter=',')
    galsdatacf.append(data2)
    ggalsdatacf.append(gdata2)

    lgals = ['NGC4889', 'NGC4874', 'NGC4839', 'NGC4873']
    print 'GALS = ', gals
    print 'GALS = ', lgals

    handles = []
    #Generate plot
    P.rcParams['font.family'] = 'Times New Roman'
    P.rcParams.update({'font.size': 20})
    f, (ax1, ax2, ax3, ax4, ax5) = P.subplots(5, sharex=True, figsize=(16,24))
    markers = ['o', 's', '^', 'D']
    colours = ['k', 'r', 'b', 'g']

    for i in xrange(len(galsdata)):
        #NaI
        if gals[i] == '4928':#Outermost point of GMP4928 affected by sky
            p1a = ax1.errorbar(galsdata[i][:-1,0], galsdata[i][:-1,1], yerr=n.sqrt(galsdata[i][:-1,2]), fmt=colours[i]+markers[i], ms=msize,
                         capthick=2, elinewidth=1.5, alpha=1)
            p1b = ax1.errorbar(galsdata[i][-1,0], galsdata[i][-1,1], yerr=n.sqrt(galsdata[i][-1,2]), fmt=colours[i]+markers[i], ms=msize,
                         capthick=2, elinewidth=1.5, alpha=1)
            p1b[-1][0].set_linestyle('--')
        else:
            p1a = ax1.errorbar(galsdata[i][:,0], galsdata[i][:,1], yerr=n.sqrt(galsdata[i][:,2]), fmt=colours[i]+markers[i], ms=msize,
                         capthick=2, elinewidth=1.5, alpha=1)
            p1aa = ax1.errorbar(galsdatacf[i][:,0], galsdatacf[i][:,1], yerr=n.sqrt(galsdatacf[i][:,2]), fmt=colours[i+1]+markers[i+1], ms=msize,
                         capthick=2, elinewidth=1.5, alpha=0.6)
        #CaT
        if gals[i] == '4928':#Outermost point of GMP4928 affected by sky
            p2a = ax2.errorbar(galsdata[i][:-1,0], galsdata[i][:-1,4], yerr=n.sqrt(galsdata[i][:-1,5]), fmt=colours[i]+markers[i], ms=msize,
                         capthick=2, elinewidth=1.5, alpha=1)
            p2b = ax2.errorbar(galsdata[i][-1,0], galsdata[i][-1,4], yerr=n.sqrt(galsdata[i][-1,5]), fmt=colours[i]+markers[i], ms=msize,
                         capthick=2, elinewidth=1.5, alpha=1)
            p2b[-1][0].set_linestyle('--')
        else:
            p2a = ax2.errorbar(galsdata[i][:,0], galsdata[i][:,4], yerr=n.sqrt(galsdata[i][:,5]), fmt=colours[i]+markers[i], ms=msize,
                         capthick=2, elinewidth=1.5, alpha=1)
            p2a = ax2.errorbar(galsdatacf[i][:,0], galsdatacf[i][:,4], yerr=n.sqrt(galsdatacf[i][:,5]), fmt=colours[i+1]+markers[i+1], ms=msize,
                         capthick=2, elinewidth=1.5, alpha=0.6)
        #MgI
        p3a = ax3.errorbar(galsdata[i][:,0], galsdata[i][:,7], yerr=n.sqrt(galsdata[i][:,8]), fmt=colours[i]+markers[i], ms=msize,
                     capthick=2, elinewidth=1.5, alpha=1)
        p3a = ax3.errorbar(galsdatacf[i][:,0], galsdatacf[i][:,7], yerr=n.sqrt(galsdatacf[i][:,8]), fmt=colours[i+1]+markers[i+1], ms=msize,
                     capthick=2, elinewidth=1.5, alpha=0.6)
        #TiO
        p4a = ax4.errorbar(galsdata[i][:,0], galsdata[i][:,10], yerr=n.sqrt(galsdata[i][:,11]), fmt=colours[i]+markers[i], ms=msize,
                     capthick=2, elinewidth=1.5, alpha=1)
        p4a = ax4.errorbar(galsdatacf[i][:,0], galsdatacf[i][:,10], yerr=n.sqrt(galsdatacf[i][:,11]), fmt=colours[i+1]+markers[i+1], ms=msize,
                     capthick=2, elinewidth=1.5, alpha=0.6)
        #FeH
        if gals[i] == '3329':#Second from last point of GMP3329 affected by sky
            p5a = ax5.errorbar(galsdata[i][:-2,0], galsdata[i][:-2,13], yerr=n.sqrt(galsdata[i][:-2,14]), fmt=colours[i]+markers[i], ms=msize,
                         capthick=2, elinewidth=1.5, alpha=1)
            p5c = ax5.errorbar(galsdata[i][-1,0], galsdata[i][-1,13], yerr=n.sqrt(galsdata[i][-1,14]), fmt=colours[i]+markers[i], ms=msize,
                         capthick=2, elinewidth=1.5, alpha=1)
            p5b = ax5.errorbar(galsdata[i][-2,0], galsdata[i][-2,13], yerr=n.sqrt(galsdata[i][-2,14]), fmt=colours[i]+markers[i], ms=msize,
                         capthick=2, elinewidth=1.5, alpha=1)
            p5b[-1][0].set_linestyle('--')
        else:
            p5a = ax5.errorbar(galsdata[i][:,0], galsdata[i][:,13], yerr=n.sqrt(galsdata[i][:,14]), fmt=colours[i]+markers[i], ms=msize,
                         capthick=2, elinewidth=1.5, alpha=1)
            p5a = ax5.errorbar(galsdatacf[i][:,0], galsdatacf[i][:,13], yerr=n.sqrt(galsdatacf[i][:,14]), fmt=colours[i+1]+markers[i+1], ms=msize,
                         capthick=2, elinewidth=1.5, alpha=0.6)
        # if model == 'True':
        #     plb = ax1.plot(galsdata[i][:,0], galsdata[i][:,3], c=colours[i], marker='x', ms=msize, linestyle='None', alpha=1)
        #     ax2.plot(galsdata[i][:,0], galsdata[i][:,6], c=colours[i], marker='x', ms=msize, linestyle='None', alpha=1)
        #     ax3.plot(galsdata[i][:,0], galsdata[i][:,9], c=colours[i], marker='x', ms=msize, linestyle='None', alpha=1)
        #     ax4.plot(galsdata[i][:,0], galsdata[i][:,12], c=colours[i], marker='x', ms=msize, linestyle='None', alpha=1)
        #     ax5.plot(galsdata[i][:,0], galsdata[i][:,15], c=colours[i], marker='x', ms=msize, linestyle='None', alpha=1)
        handles.append(P.Line2D((),(), color=colours[i], marker=markers[i], linestyle='', markersize=msize))
    ax1.set_ylim([0.2, 0.9])
    ax5.set_ylim([0.,1.0])

    for i in xrange(len(ggalsdata)):
        dat = ggalsdata[i]
        datcf = ggalsdatacf[i]
        grval = 4.0
        fc = 'w'
        # handles.append(P.Line2D((),(), color=colours[gals.index('3367')], marker=markers[gals.index('3367')], linestyle='', markersize=msize))
        #NaI
        ax1.errorbar(grval, dat[0], yerr=dat[1], fmt=markers[i], ms=msize+1, mec=colours[i],
        mew=1.3, mfc=fc, ls='', ecolor=colours[i], capthick=2, elinewidth=1.5, zorder=10)
        ax1.errorbar(grval, datcf[0], yerr=datcf[1], fmt=markers[i+1], ms=msize+1, mec=colours[i+1],
        mew=1.3, mfc=fc, ls='', ecolor=colours[i+1], capthick=2, elinewidth=1.5, zorder=10, alpha=0.6)
        #CaT*
        ax2.errorbar(grval, dat[2], yerr=dat[3], fmt=markers[i], ms=msize+1, mec=colours[i],
        mew=1.3, mfc=fc, ls='', ecolor=colours[i], capthick=2, elinewidth=1.5, zorder=10)
        ax2.errorbar(grval, datcf[2], yerr=datcf[3], fmt=markers[i+1], ms=msize+1, mec=colours[i+1],
        mew=1.3, mfc=fc, ls='', ecolor=colours[i+1], capthick=2, elinewidth=1.5, zorder=10, alpha=0.6)
        #MgI
        ax3.errorbar(grval, dat[4], yerr=dat[5], fmt=markers[i], ms=msize+1, mec=colours[i],
        mew=1.3, mfc=fc, ls='', ecolor=colours[i], capthick=2, elinewidth=1.5, zorder=10)
        ax3.errorbar(grval, datcf[4], yerr=datcf[5], fmt=markers[i+1], ms=msize+1, mec=colours[i+1],
        mew=1.3, mfc=fc, ls='', ecolor=colours[i+1], capthick=2, elinewidth=1.5, zorder=10, alpha=0.6)
        #TiO
        ax4.errorbar(grval, dat[6], yerr=dat[7], fmt=markers[i], ms=msize+1, mec=colours[i],
        mew=1.3, mfc=fc, ls='', ecolor=colours[i], capthick=2, elinewidth=1.5, zorder=10)
        ax4.errorbar(grval, datcf[6], yerr=datcf[7], fmt=markers[i+1], ms=msize+1, mec=colours[i+1],
        mew=1.3, mfc=fc, ls='', ecolor=colours[i+1], capthick=2, elinewidth=1.5, zorder=10, alpha=0.6)
        #Fe
        ax5.errorbar(grval, dat[8], yerr=dat[9], fmt=markers[i], ms=msize+1, mec=colours[i],
        mew=1.3, mfc=fc, ls='', ecolor=colours[i], capthick=2, elinewidth=1.5, zorder=10)
        ax5.errorbar(grval, datcf[8], yerr=datcf[9], fmt=markers[i+1], ms=msize+1, mec=colours[i+1],
        mew=1.3, mfc=fc, ls='', ecolor=colours[i+1], capthick=2, elinewidth=1.5, zorder=10, alpha=0.6)

    # Fine-tune figure; make subplots close to each other and hide x ticks for all but bottom plot.
    f.subplots_adjust(hspace=0.08)
    P.xlabel('r [arcsec]', fontsize=tfsize+2)
    ax1.minorticks_on()
    ax2.minorticks_on()
    ax3.minorticks_on()
    ax4.minorticks_on()
    ax5.minorticks_on()
    ax1.set_ylabel('NaI$_{\\rm{SDSS}}$0.82 [$\AA$]', fontsize=tfsize+2)
    ax2.set_ylabel('CaT [$\AA$]', fontsize=tfsize+2)
    ax3.set_ylabel('MgI0.88 [$\AA$]', fontsize=tfsize+2)
    ax4.set_ylabel('TiO0.89', fontsize=tfsize+2)
    ax5.set_ylabel('FeH0.99 [$\AA$]', fontsize=tfsize+2)
    ax1.yaxis.label.set_fontsize(tfsize+2)
    ax2.yaxis.label.set_fontsize(tfsize+2)
    ax3.yaxis.label.set_fontsize(tfsize+2)
    ax4.yaxis.label.set_fontsize(tfsize+2)
    ax5.yaxis.label.set_fontsize(tfsize+2)
    ax5.xaxis.label.set_fontsize(tfsize+2)
    ax1.tick_params(axis='both', which='major', direction='in', length=majtsize, width=twidth)
    ax1.tick_params(axis='both', which='minor', direction='in', length=mintsize, width=twidth)
    ax2.tick_params(axis='both', which='major', direction='in', length=majtsize, width=twidth)
    ax2.tick_params(axis='both', which='minor', direction='in', length=mintsize, width=twidth)
    ax3.tick_params(axis='both', which='major', direction='in', length=majtsize, width=twidth)
    ax3.tick_params(axis='both', which='minor', direction='in', length=mintsize, width=twidth)
    ax4.tick_params(axis='both', which='major', direction='in', length=majtsize, width=twidth)
    ax4.tick_params(axis='both', which='minor', direction='in', length=mintsize, width=twidth)
    ax5.tick_params(axis='both', which='major', direction='in', length=majtsize, width=twidth)
    ax5.tick_params(axis='both', which='minor', direction='in', length=mintsize, width=twidth)
    ax1.set_xlim([-0.5,10.])

    # leg = ax2.legend(handles, [lgals[i]+' ('+sigs[i]+' km s$^{-1}$)' for i in xrange(len(gals))], loc='upper right',
    #             frameon=False, borderaxespad=0.5, labelspacing=0.2, numpoints=1)
    leg = ax2.legend([p1a, p1aa], [lgals[i]+' ('+sigs[i]+' km s$^{-1}$)', lgals[i]+' (w. Correction Factor)'], loc='upper right',
                frameon=False, borderaxespad=0.5, labelspacing=0.2, numpoints=1)
    ltext = leg.get_texts()
    for j in xrange(len(ltext)):
        P.setp(ltext[j], fontsize=20, color=colours[j])

    ###Add kpc x-axis along top
    ax11 = ax1.twiny()
    new_tick_locations = n.array([0, 2, 4, 6, 8, 10])
    kpc_per_arcsec = NWCC.calc(0.024, 68.0, 0.3, 0.7, 'general')[2]
    new_ticks = n.round(n.array([0,2,4,6,8,10])*kpc_per_arcsec, 0)
    ax11.set_xlim(ax1.get_xlim())
    ax11.set_xticks(new_tick_locations)
    ax11.set_xticklabels(new_ticks.astype(n.int))
    ax11.set_xlabel('r [kpc]', fontsize=tfsize+2)
    ax11.xaxis.label.set_fontsize(tfsize+2)
    ax11.tick_params(axis='both', which='major', direction='in', length=majtsize, width=twidth)
    ax11.tick_params(axis='both', which='minor', direction='in', length=mintsize, width=twidth)
    ax11.xaxis.set_minor_locator(AutoMinorLocator())
    ###

    # P.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
    f.tight_layout()
    f.subplots_adjust(hspace=0.08)
    if save != 'False':
        P.savefig('/Users/zieleniewski/Data/swift/Coma_stellar_pops/'+str(save)+'.png')
    P.show()
#------------#



#------------#
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




#------------#
def normalise_spec(spec, fcont, mids, vspec=None):
    '''Normalise spectrum by psdueo-continuum

    Inputs:

        - spec: spectrum class
        - fcont: array of continuum positions [A]
        - mids: midpoints of the pseudo-continuum bands (bmid, rmid) (A, A)
        - vspec: variance spectrum class. Divided by (continuum fit)^2 if given

    Outputs:

        - divspec

    '''

    #Normalise spectrum by pseudo-continuum
    indbmid = mids[0]; indrmid = mids[1]
    spec_bcont = n.mean(spec.flam[0][n.where(spec.lam > fcont[0])[0][0]:n.where(spec.lam > fcont[1])[0][0]])
    spec_rcont = n.mean(spec.flam[0][n.where(spec.lam > fcont[-2])[0][0]:n.where(spec.lam > fcont[-1])[0][0]])
    spec_contfit = n.polyfit([indbmid, indrmid], [spec_bcont, spec_rcont], 1)
    contspec = n.column_stack((spec.lam, spec.lam*spec_contfit[0] + spec_contfit[1]))
    divspec = s.spectrum(lam=spec.lam, lamspec=(spec.flam[0]/contspec[:,1]))
    if vspec:
        divvspec = s.spectrum(lam=vspec.lam, lamspec=(vspec.flam[0]/(contspec[:,1])**2))
        return divspec, divvspec

    return divspec
#------------#



#------------#
def mask_spec(inspec, it, lines, mspec):
    '''Function to mask residual bad regions of spectrum

    Inputs:

        - spec: spectrum [wave, flux]
        - it: spectrum iteration (number)
        - lines: numpy array of bad regions [s_wave, e_wave]
        - mspec: model spec for replacing bad regions [wave, flux]

    Outputs:

        - maskspec

    '''

    #####Mask remaining skylines

    spec = inspec.copy()
    skyline_pos = []
    #Find and mask sky line regions
    for j in xrange(lines.shape[0]):
        if it == 'glob' or int(lines[j,2]) == int(it) or int(lines[j,2]) == -1:
            sline_start = n.where(spec[:,0]>lines[j,0])[0][0]
            sline_end = n.where(spec[:,0]>lines[j,1])[0][0]
            if lines[j,3] == 1:
                sline_interp = si.interp1d([spec[sline_start,0],spec[sline_end,0]],[spec[sline_start,1],spec[sline_end,1]])
                maskspec = sline_interp(spec[sline_start:sline_end,0])
                spec[sline_start:sline_end,1] = maskspec
            elif lines[j,3] == 0:
                spec[sline_start:sline_end,1] = mspec[sline_start:sline_end,1]
            # skyline_positions = n.arange(sline_start,sline_end+1,1)
            skyline_pos.append(n.column_stack((spec[sline_start:sline_end,0],spec[sline_start:sline_end,1])))

    return spec, skyline_pos
#------------#



#------------#
def plot_all_specs(gal, ind, out_sig='intrinsic', skymask='False', savevals='False', fehopt='False', osfitted='False', gspecs='False', snrs='False'):
    #gal = sys.argv[1]
    #ind = sys.argv[3]

    # na, nacat, fehc, fehcos, fehcosos, feho, skylines_nai, skylines_cat, skylines_mgi, skylines_tio, skylines_feh = get_data(gal)
    na, nacat, fehc, fehcosos, feho, skylines_nai, skylines_cat, skylines_mgi, skylines_tio, skylines_feh = get_data(gal)

    #Wavelengths
    lams = nacat['wave']*10000.

    #Sausage spectrum to get index definitions
    sspec = s.spectrum(lam=lams, lamspec=na['skysubtdgal'][:,0]/n.median(na['skysubtdgal'][:,0]))

    if ind == 'NaI':
        start = n.where(lams > 8250.)[0][0]
        end = n.where(lams > 8500.)[0][0]
        fcont = sspec.nacont
        feat = sspec.na
        skylines = skylines_nai

    if ind == 'NaIsdss':
        start = n.where(lams > 8250.)[0][0]
        end = n.where(lams > 8500.)[0][0]
        fcont = sspec.nacont_sdss
        feat = sspec.na_sdss
        skylines = skylines_nai

    if ind == 'CaT':
        start = n.where(lams > 8400.)[0][0]
        end = n.where(lams > 9100.)[0][0]
        fcont = sspec.cacont
        feat = sspec.ca
        skylines = skylines_cat

    if ind == 'MgI':
        start = n.where(lams > 8900.)[0][0]
        end = n.where(lams > 9150.)[0][0]
        fcont = sspec.mgcont
        feat = sspec.mg
        skylines = skylines_mgi

    if ind == 'TiO':
        start = n.where(lams > 8950.)[0][0]
        end = n.where(lams > 9150.)[0][0]
        fcont = sspec.tiocont
        feat = sspec.tio
        skylines = skylines_tio

    if ind == 'FeH':
        # start = n.where(lams > 9800.)[0][0]
        # end = n.where(lams > 10300.)[0][0]
        start = n.where(lams > 9700.)[0][0]
        end = n.where(lams > 10400.)[0][0]
        fcont = sspec.fehcont
        feat = sspec.feh
        skylines = skylines_feh

    print 'ind = ', ind
    print 'Feature = ', feat
    print 'Continuum = ', fcont

    #Cut wavelength range
    indlams = lams[start:end]

    nacatspecs = []
    naspecs = []
    fehcspecs = []
    fehcosspecs = []
    fehcososspecs = []
    fehospecs = []
    nacatvars = []
    navars = []
    fehcvars = []
    fehcosvars = []
    fehcososvars = []
    fehovars = []

##    if ind == 'CaT':
##        pass
##
##    else:
    print 'FEAT = ', feat
    print 'FCONT = ', fcont
    indbmid = (fcont[0]+fcont[1])/2.
    indrmid = (fcont[-1]+fcont[-1])/2.

    rs = nacat['rbar']*0.235
    rs[0] = 0.0
    print 'RS = ', rs

    #Instrument sigma = 50.4 km/s for NaCaT; 63.3 km/s for FeH
    fehinstsig = 0.0#63.3
    nainstsig = 0.0#50.4

    if gspecs == 'True':
        gz = nacat['gsol'][0]*1000./sc.c
        nacatgspec = ssz.z_correct(n.column_stack((indlams,nacat['skysubtdggal'][start:end]/n.median(nacat['skysubtdggal'][start:end]))), gz)
        nacatgvar = ssz.z_correct(n.column_stack((indlams,nacat['skysubtdggerr'][start:end]**2/n.median(nacat['skysubtdggal'][start:end]**2))), gz)
        nagspec = ssz.z_correct(n.column_stack((indlams,na['skysubtdggal'][start:end]/n.median(na['skysubtdggal'][start:end]))), gz)
        nagvar = ssz.z_correct(n.column_stack((indlams,na['skysubtdggerr'][start:end]**2/n.median(na['skysubtdggal'][start:end]**2))), gz)
        fehcgspec = ssz.z_correct(n.column_stack((indlams,fehc['skysubtdggal'][start:end]/n.median(fehc['skysubtdggal'][start:end]))), gz)
        fehcgvar = ssz.z_correct(n.column_stack((indlams,fehc['skysubtdggerr'][start:end]**2/n.median(fehc['skysubtdggal'][start:end]**2))), gz)
        if gal == 'GMP3329' and ind in ['CaT', 'MgI'] and skymask == 'True':
            nacatgmspec = ssz.z_correct(n.column_stack((indlams,nacat['bestfitgMODEL'][start:end]/n.median(nacat['skysubtdggal'][start:end]))), gz)
            nacatgspec, nacatglinepos = mask_spec(nacatgspec, 'glob', skylines, nacatgmspec)
        if gal == 'GMP3367' and ind =='FeH' and skymask == 'True':
            fehcgmspec = ssz.z_correct(n.column_stack((indlams,fehc['bestfitgMODEL'][start:end]/n.median(fehc['skysubtdggal'][start:end]))), gz)
            fehcgspec, feholinepos = mask_spec(fehcgspec, 'glob', skylines, fehcgmspec)
        nacatgspec = s.spectrum(lam=nacatgspec[:,0], lamspec=nacatgspec[:,1])
        nacatgvar = s.spectrum(lam=nacatgvar[:,0], lamspec=nacatgvar[:,1])
        nagspec = s.spectrum(lam=nagspec[:,0], lamspec=nagspec[:,1])
        nagvar = s.spectrum(lam=nagvar[:,0], lamspec=nagvar[:,1])
        fehcgspec = s.spectrum(lam=fehcgspec[:,0], lamspec=fehcgspec[:,1])
        fehcgvar = s.spectrum(lam=fehcgvar[:,0], lamspec=fehcgvar[:,1])

    for i in xrange(len(nacat['binflux'][0,:])):
        #Redshift correct spectra
        z = nacat['param'][0,i]*1000./sc.c
        sig = nacat['param'][1,i]
        print 'z = ', z
        print 'sig = ', sig, 'km/s'
        nacatspec = ssz.z_correct(n.column_stack((indlams,nacat['skysubtdgal'][start:end,i]/n.median(nacat['skysubtdgal'][start:end,i]))), z)
        nacatmspec = ssz.z_correct(n.column_stack((indlams,nacat['bestfitMODEL'][start:end,i]/n.median(nacat['bestfitMODEL'][start:end,i]))), z)
        nacatvar = ssz.z_correct(n.column_stack((indlams,nacat['skysubtdgerr'][start:end,i]**2/n.median(nacat['skysubtdgal'][start:end,i]**2))), z)

        naspec = ssz.z_correct(n.column_stack((indlams,na['skysubtdgal'][start:end,i]/n.median(na['skysubtdgal'][start:end,i]))), z)
        namspec = ssz.z_correct(n.column_stack((indlams,na['bestfitMODEL'][start:end,i]/n.median(na['bestfitMODEL'][start:end,i]))), z)
        navar = ssz.z_correct(n.column_stack((indlams,na['skysubtdgerr'][start:end,i]**2/n.median(na['skysubtdgal'][start:end,i]**2))), z)

        fehcspec = ssz.z_correct(n.column_stack((indlams,fehc['skysubtdgal'][start:end,i]/n.median(fehc['skysubtdgal'][start:end,i]))), z)
        fehcmspec = ssz.z_correct(n.column_stack((indlams,fehc['bestfitMODEL'][start:end,i]/n.median(fehc['bestfitMODEL'][start:end,i]))), z)
        fehcvar = ssz.z_correct(n.column_stack((indlams,fehc['skysubtdgerr'][start:end,i]**2/n.median(fehc['skysubtdgal'][start:end,i]**2))), z)

        # fehcosspec = ssz.z_correct(n.column_stack((indlams,fehcos['skysubtdgal'][start:end,i]/n.median(fehcos['skysubtdgal'][start:end,i]))), z)
        # fehcosmspec = ssz.z_correct(n.column_stack((indlams,fehcos['bestfitMODEL'][start:end,i]/n.median(fehcos['bestfitMODEL'][start:end,i]))), z)
        # fehcosvar = ssz.z_correct(n.column_stack((indlams,fehcos['skysubtdgerr'][start:end,i]**2/n.median(fehcos['skysubtdgal'][start:end,i]**2))), z)

        fehcososspec = ssz.z_correct(n.column_stack((indlams,fehcosos['skysubtdgal'][start:end,i]/n.median(fehcosos['skysubtdgal'][start:end,i]))), z)
        fehcososmspec = ssz.z_correct(n.column_stack((indlams,fehcosos['bestfitMODEL'][start:end,i]/n.median(fehcosos['bestfitMODEL'][start:end,i]))), z)
        fehcososvar = ssz.z_correct(n.column_stack((indlams,fehcosos['skysubtdgerr'][start:end,i]**2/n.median(fehcosos['skysubtdgal'][start:end,i]**2))), z)

        fehospec = ssz.z_correct(n.column_stack((indlams,feho['skysubtdgal'][start:end,i]/n.median(feho['skysubtdgal'][start:end,i]))), z)
        fehomspec = ssz.z_correct(n.column_stack((indlams,feho['bestfitMODEL'][start:end,i]/n.median(feho['bestfitMODEL'][start:end,i]))), z)
        fehovar = ssz.z_correct(n.column_stack((indlams,feho['skysubtdgerr'][start:end,i]**2/n.median(feho['skysubtdgal'][start:end,i]**2))), z)

        #####Mask remaining skylines
        if ind in ['NaI', 'NaIsdss', 'CaT', 'MgI', 'TiO', 'FeH'] and skymask=='True':
            print 'Masking skylines'
            naspec, nalinepos = mask_spec(naspec, i, skylines, namspec)
            nacatspec, nacatlinepos = mask_spec(nacatspec, i, skylines, nacatmspec)
            fehcspec, fehclinepos = mask_spec(fehcspec, i, skylines, fehcmspec)
            # fehcosspec, fehcoslinepos = mask_spec(fehcosspec, i, skylines, fehcosmspec)
            fehcososspec, fehcososlinepos = mask_spec(fehcososspec, i, skylines, fehcososmspec)
            fehospec, feholinepos = mask_spec(fehospec, i, skylines, fehomspec)

        #Load as spectrum class
        nacatspec = s.spectrum(lam=nacatspec[:,0], lamspec=nacatspec[:,1])
        nacatvar = s.spectrum(lam=nacatvar[:,0], lamspec=nacatvar[:,1])
        naspec = s.spectrum(lam=naspec[:,0], lamspec=naspec[:,1])
        navar = s.spectrum(lam=navar[:,0], lamspec=navar[:,1])
        fehcspec = s.spectrum(lam=fehcspec[:,0], lamspec=fehcspec[:,1])
        fehcvar = s.spectrum(lam=fehcvar[:,0], lamspec=fehcvar[:,1])
        # fehcosspec = s.spectrum(lam=fehcosspec[:,0], lamspec=fehcosspec[:,1])
        # fehcosvar = s.spectrum(lam=fehcosvar[:,0], lamspec=fehcosvar[:,1])
        fehcososspec = s.spectrum(lam=fehcososspec[:,0], lamspec=fehcososspec[:,1])
        fehcososvar = s.spectrum(lam=fehcososvar[:,0], lamspec=fehcososvar[:,1])
        fehospec = s.spectrum(lam=fehospec[:,0], lamspec=fehospec[:,1])
        fehovar = s.spectrum(lam=fehovar[:,0], lamspec=fehovar[:,1])

        ###Convolve to common velocity dispersion
        if out_sig != 'intrinsic':
            out_sig = float(out_sig)

            naconvsig = n.sqrt(out_sig**2 - (nacat['param'][1,i]**2+nainstsig**2))
            fehconvsig = n.sqrt(out_sig**2 - (nacat['param'][1,i]**2+fehinstsig**2))

            #convsig = n.sqrt(out_sig**2 - sig**2)
            loglambs = n.linspace(n.log10(naspec.lam[0]), n.log10(naspec.lam[-1]), len(naspec.lam))
            veldisp = (10**(loglambs[1]-loglambs[0]) - 1)*sc.c/1000. #in km/s
            nagauss = ssz.Gauss(naconvsig, veldisp)
            fehgauss = ssz.Gauss(fehconvsig, veldisp)

            nacatspec, nacatvar = convolve_spec(nacatspec, naconvsig, nagauss, nacatvar)

            naspec, navar = convolve_spec(naspec, naconvsig, nagauss, navar)

            fehcspec, fehcvar = convolve_spec(fehcspec, fehconvsig, fehgauss, fehcvar)

            # fehcosspec, fehcosvar = convolve_spec(fehcosspec, fehconvsig, fehgauss, fehcosvar)

            fehcososspec, fehcososvar = convolve_spec(fehcososspec, fehconvsig, fehgauss, fehcososvar)

            fehospec, fehovar = convolve_spec(fehospec, fehconvsig, fehgauss, fehovar)


        else:
            pass

        nacatspecs.append(nacatspec)
        naspecs.append(naspec)
        nacatvars.append(nacatvar)
        navars.append(navar)

        fehcdivspec, fehcdivvar = normalise_spec(fehcspec, fcont, [indbmid, indrmid], fehcvar)
        # fehcosdivspec, fehcosdivvar = normalise_spec(fehcosspec, fcont, [indbmid, indrmid], fehcosvar)
        fehcososdivspec, fehcososdivvar = normalise_spec(fehcososspec, fcont, [indbmid, indrmid], fehcososvar)
        fehodivspec, fehodivvar = normalise_spec(fehospec, fcont, [indbmid, indrmid], fehovar)
        fehcspecs.append(fehcdivspec)
        # fehcosspecs.append(fehcosdivspec)
        fehcososspecs.append(fehcososdivspec)
        fehospecs.append(fehodivspec)
        fehcvars.append(fehcdivvar)
        # fehcosvars.append(fehcosdivvar)
        fehcososvars.append(fehcososdivvar)
        fehovars.append(fehodivvar)

    if gspecs == 'True':
        if out_sig != 'intrinsic':
            naconvgsig = n.sqrt(out_sig**2 - (nacat['gsol'][1]**2+nainstsig**2))
            fehconvgsig = n.sqrt(out_sig**2 - (nacat['gsol'][1]**2+fehinstsig**2))
            naggauss = ssz.Gauss(naconvgsig, veldisp)
            fehggauss = ssz.Gauss(fehconvgsig, veldisp)
            nacatgspec = convolve_spec(nacatgspec, naconvgsig, naggauss)
            nacatgvar = convolve_spec(nacatgvar, naconvgsig, naggauss)
            nagspec = convolve_spec(nagspec, naconvgsig, naggauss)
            nagvar = convolve_spec(nagvar, naconvgsig, naggauss)
            fehcgspec = convolve_spec(fehcgspec, fehconvgsig, fehggauss)
            fehcgvar = convolve_spec(fehcgvar, fehconvgsig, fehggauss)
        else:
            pass
        fehcdivgspec, fehcdivgvar = normalise_spec(fehcgspec, fcont, [indbmid, indrmid], fehcgvar)

    #PLOT SCIENCE SPECTRA
    P.rcParams['font.family'] = 'Times New Roman'
    P.rcParams.update({'font.size': 20})
    P.figure(figsize=(10,8), dpi=72)

    number = len(naspecs)
    mycmap = P.get_cmap('copper')
    #colours = [cmap(i) for i in n.linspace(0, 1, number)]
    colours = [mycmap(cols) for cols in (rs/float(len(rs)))]
    #print 'COLOURS = ', colours

    if gal in ['GMP2921', 'GMP3329', 'GMP4928']:
        if ind in ['CaT', 'MgI', 'TiO']:
            #if gal == 'GMP4928':
            for ii in xrange(len(nacatspecs)):
                P.plot(nacatspecs[ii].lam, nacatspecs[ii].flam[0], color=colours[ii], ls='-',
                       lw=2.5, label=str(n.round(rs[ii],1))+"''")

        elif ind in ['NaI', 'NaIsdss']:
            for ii in xrange(len(naspecs)):
                P.plot(naspecs[ii].lam, naspecs[ii].flam[0], color=colours[ii], ls='-',
                       lw=2.5, label=str(n.round(rs[ii],1))+"''")

        elif ind == 'FeH':
            P.plot(fehcspecs[0].lam, fehcspecs[0].flam[0], color=colours[0], ls='-',
                   lw=2.5, label='Fit Div.')
            if osfitted == 'True':
                    P.plot(fehcososspecs[0].lam, fehcososspecs[0].flam[0], color=colours[0], ls='-.',
                           lw=2.5, label='O-S fitted')
                    P.legend(loc='best', frameon=False, borderaxespad=0.5, labelspacing=0.18, handlelength=2.2)
            if fehopt == 'True':
                    P.plot(fehospecs[0].lam, fehospecs[0].flam[0], color=colours[0], ls=':',
                        lw=2.5, label='Tell. Corr.')
                    P.legend(loc='best', frameon=False, borderaxespad=0.5, labelspacing=0.18, handlelength=2.2)
            for ii in xrange(len(naspecs)):
                P.plot(fehcspecs[ii].lam, fehcspecs[ii].flam[0], color=colours[ii], ls='-',
                       lw=2.5)
                if osfitted=='True':
                    P.plot(fehcososspecs[ii].lam, fehcososspecs[ii].flam[0], color=colours[ii], ls='-.',
                           lw=2.5)
                if fehopt == 'True':
                    P.plot(fehospecs[ii].lam, fehospecs[ii].flam[0], color=colours[ii], ls=':',
                        lw=2.5)

    if gspecs == 'True':
        if ind in ['NaI', 'NaIsdss']:
            P.plot(nagspec.lam, nagspec.flam[0], color='k', ls='--',
                   lw=2.5, label='Global spec')
        if ind in ['CaT', 'MgI', 'TiO']:
            P.plot(nacatgspec.lam, nacatgspec.flam[0], color='k', ls='--',
                   lw=2.5, label='Global spec')
        if ind == 'FeH':
            P.plot(fehcdivgspec.lam, fehcdivgspec.flam[0], color='k', ls='--',
                   lw=2.5, label='Global spec')
        # P.legend(loc='upper right', frameon=False, borderaxespad=0.2, labelspacing=0.18, handlelength=2.2)

    if ind == 'NaI' or ind == 'NaIsdss':
        leg = P.legend(loc='lower left', frameon=False, borderaxespad=0.2, labelspacing=0.18, handlelength=2.2)
        ltext = P.gca().get_legend().get_texts()
        for legi in xrange(len(ltext)):
            P.setp(ltext[legi], fontsize=tfsize)
        for legobj in leg.legendHandles:
            legobj.set_linewidth(3.0)

    if ind in ['NaI', 'NaIsdss', 'MgI', 'TiO', 'FeH']:
        # P.plot([fcont[0],fcont[0]], [0.5,1.5], 'k-',
        #        [fcont[1],fcont[1]], [0.5,1.5], 'k-',
        #        [fcont[-2],fcont[-2]], [0.5,1.5], 'k-',
        #        [fcont[-1],fcont[-1]], [0.5,1.5], 'k-', alpha=0.5, lw=2.2)
        # P.plot([feat[0],feat[0]], [0.5,1.5], 'r-', alpha=0.5, zorder=10, lw=2.2)
        # P.plot([feat[1],feat[1]], [0.5,1.5], 'r-', alpha=0.5, zorder=10, lw=2.2)
        P.axvspan(feat[0], feat[1], color='r', alpha=0.15, lw=0.0)
        P.axvspan(fcont[0], fcont[1], color='k', alpha=0.15, lw=0.0)
        P.axvspan(fcont[2], fcont[3], color='k', alpha=0.15, lw=0.0)
    elif ind == 'CaT':
        # P.plot([fcont[0],fcont[0]], [0.5,1.5], 'k-',
        #        [fcont[1],fcont[1]], [0.5,1.5], 'k-',
        #        [fcont[2],fcont[2]], [0.5,1.5], 'k-',
        #        [fcont[3],fcont[3]], [0.5,1.5], 'k-',
        #        [fcont[4],fcont[4]], [0.5,1.5], 'k-',
        #        [fcont[5],fcont[5]], [0.5,1.5], 'k-',
        #        [fcont[6],fcont[6]], [0.5,1.5], 'k-',
        #        [fcont[7],fcont[7]], [0.5,1.5], 'k-',
        #        [fcont[-2],fcont[-2]], [0.5,1.5], 'k-',
        #        [fcont[-1],fcont[-1]], [0.5,1.5], 'k-', alpha=0.5, lw=2.2)
        # P.plot([feat[0],feat[0]], [0.5,1.5], 'r--', alpha=0.5, zorder=10, lw=2.2)
        # P.plot([feat[1],feat[1]], [0.5,1.5], 'r--', alpha=0.5, zorder=10, lw=2.2)
        # P.plot([feat[2],feat[2]], [0.5,1.5], 'r--', alpha=0.5, zorder=10, lw=2.2)
        # P.plot([feat[3],feat[3]], [0.5,1.5], 'r--', alpha=0.5, zorder=10, lw=2.2)
        # P.plot([feat[4],feat[4]], [0.5,1.5], 'r--', alpha=0.5, zorder=10, lw=2.2)
        # P.plot([feat[5],feat[5]], [0.5,1.5], 'r--', alpha=0.5, zorder=10, lw=2.2)
        P.axvspan(feat[0], feat[1], color='r', alpha=0.15, lw=0.0)
        P.axvspan(feat[2], feat[3], color='r', alpha=0.15, lw=0.0)
        P.axvspan(feat[4], feat[5], color='r', alpha=0.15, lw=0.0)
        P.axvspan(fcont[0], fcont[1], color='k', alpha=0.15, lw=0.0)
        P.axvspan(fcont[2], fcont[3], color='k', alpha=0.15, lw=0.0)
        P.axvspan(fcont[4], fcont[5], color='k', alpha=0.15, lw=0.0)
        P.axvspan(fcont[6], fcont[7], color='k', alpha=0.15, lw=0.0)
        P.axvspan(fcont[8], fcont[9], color='k', alpha=0.15, lw=0.0)

    if skymask == 'True' and out_sig=='intrinsic':
        if ind in ['NaI', 'NaIsdss']:
            for j in xrange(len(nalinepos)):
                P.plot(nalinepos[j][:,0], nalinepos[j][:,1], color='r', ls='--', lw=2.5)
        if ind in ['CaT', 'MgI', 'TiO']:
            for j in xrange(len(nacatlinepos)):
                P.plot(nacatlinepos[j][:,0], nacatlinepos[j][:,1], color='r', ls='--', lw=2.5)
        if ind == 'FeH':
            for j in xrange(len(fehclinepos)):
                P.plot(fehclinepos[j][:,0], fehclinepos[j][:,1], color='r', ls='--', lw=2.5)

    #For NaI
    if ind == 'NaI' or ind == 'NaIsdss':
        P.xlim([8130., 8250.])
        P.ylim([0.94,1.06])
        P.ylabel('Normalised flux', fontsize=tfsize+3)

    #For CaT
    elif ind == 'CaT':
        P.xlim([8400., 8800.])
        P.ylim([0.85,1.05])

    #For FeH
    elif ind == 'FeH':
        P.xlim([9840., 9980.])
        P.ylim([0.96,1.04])

    #For MgI
    elif ind == 'MgI':
        P.xlim([8765., 8865.])
        P.ylim([0.93,1.13])

    #For TiO
    elif ind == 'TiO':
        P.xlim([8810., 8910.])
        P.ylim([0.93,1.13])

    P.xlabel('$\lambda$ [$\AA$]', fontsize=tfsize+3)

    P.xticks(size=lfsize+3)
    P.yticks(size=lfsize+3)
    P.minorticks_on()
    P.tick_params(axis='both', which='major', direction='in', length=majtsize, width=twidth)
    P.tick_params(axis='both', which='major', direction='in', length=mintsize, width=twidth)

    # if out_sig != 'intrinsic':
    #     P.annotate('$\sigma$ = '+str(n.round(out_sig,1))+'km s$^{-1}$', xy=(feat[0]+2,1.02), textcoords='axes fraction',
    #                xytext=(0.35,0.8), fontsize=lfsize, color='k')

    P.tight_layout()
    P.show()
    #END SCIENCE SPEC PLOTS

    #PLOT S/N SPECTRA
    if snrs=='True':
        P.figure(figsize=(10,8), dpi=72)
        if ind in ['CaT', 'MgI', 'TiO']:
            #if gal == 'GMP4928':
            for ii in xrange(len(nacatspecs)):
                P.plot(nacatspecs[ii].lam, nacatspecs[ii].flam[0]/n.sqrt(nacatvars[ii].flam[0]), color=colours[ii], ls='-',
                       lw=2.5, label=str(n.round(rs[ii],1))+"''")

        elif ind in ['NaI', 'NaIsdss']:
            for ii in xrange(len(naspecs)):
                P.plot(naspecs[ii].lam, naspecs[ii].flam[0]/n.sqrt(navars[ii].flam[0]), color=colours[ii], ls='-',
                       lw=2.5, label=str(n.round(rs[ii],1))+"''")

        elif ind == 'FeH':
            P.plot(fehcspecs[0].lam, fehcspecs[0].flam[0]/n.sqrt(fehcvars[0].flam[0]), color=colours[0], ls='-',
                   lw=2.5, label='Fit Div.')
            if osfitted == 'True':
                    P.plot(fehcososspecs[0].lam, fehcososspecs[0].flam[0]/n.sqrt(fehcososvars[0].flam[0]), color=colours[0], ls='-.',
                           lw=2.5, label='O-S fitted')
                    P.legend(loc='best', frameon=False, borderaxespad=0.5, labelspacing=0.18, handlelength=2.2)
            if fehopt == 'True':
                    P.plot(fehospecs[0].lam, fehospecs[0].flam[0]/n.sqrt(fehovars[0].flam[0]), color=colours[0], ls=':',
                        lw=2.5, label='Tell. Corr.')
                    P.legend(loc='best', frameon=False, borderaxespad=0.5, labelspacing=0.18, handlelength=2.2)
            for ii in xrange(len(naspecs)):
                P.plot(fehcspecs[ii].lam, fehcspecs[ii].flam[0]/n.sqrt(fehcvars[ii].flam[0]), color=colours[ii], ls='-',
                       lw=2.5)
                if osfitted=='True':
                    P.plot(fehcososspecs[ii].lam, fehcososspecs[ii].flam[0]/n.sqrt(fehcososvars[ii].flam[0]), color=colours[ii], ls='-.',
                           lw=2.5)
                if fehopt == 'True':
                    P.plot(fehospecs[ii].lam, fehospecs[ii].flam[0]/n.sqrt(fehovars[ii].flam[0]), color=colours[ii], ls=':',
                        lw=2.5)

        if gspecs == 'True':
            if ind in ['NaI', 'NaIsdss']:
                P.plot(nagspec.lam, nagspec.flam[0]/n.sqrt(nagvar.flam[0]), color='k', ls='--',
                       lw=2.5, label='Global spec')
            if ind in ['CaT', 'MgI', 'TiO']:
                P.plot(nacatgspec.lam, nacatgspec.flam[0]/n.sqrt(nacatgvar.flam[0]), color='k', ls='--',
                       lw=2.5, label='Global spec')
            if ind == 'FeH':
                P.plot(fehcdivgspec.lam, fehcdivgspec.flam[0]/n.sqrt(fehcdivgvar.flam[0]), color='k', ls='--',
                       lw=2.5, label='Global spec')
            # P.legend(loc='upper right', frameon=False, borderaxespad=0.2, labelspacing=0.18, handlelength=2.2)

        if ind == 'NaI' or ind == 'NaIsdss':
            P.legend(loc='lower left', frameon=False, borderaxespad=0.2, labelspacing=0.18, handlelength=2.2)

        #For NaI
        if ind == 'NaI' or ind == 'NaIsdss':
            P.xlim([8130., 8250.])
            P.ylabel('S/N', fontsize=tfsize)

        #For CaT
        elif ind == 'CaT':
            P.xlim([8400., 8800.])

        #For FeH
        elif ind == 'FeH':
            P.xlim([9840., 9980.])

        #For MgI
        elif ind == 'MgI':
            P.xlim([8765., 8865.])

        #For TiO
        elif ind == 'TiO':
            P.xlim([8810., 8910.])

        P.xlabel('$\lambda$ [$\AA$]', fontsize=tfsize)

        P.xticks(size=lfsize)
        P.yticks(size=lfsize)

        P.tick_params(axis='both', which='major', direction='in', length=majtsize, width=twidth)
        P.tick_params(axis='both', which='major', direction='in', length=mintsize, width=twidth)

        P.tight_layout()
        P.show()
        #END S/N SPEC PLOTS

        #CALCULATE MEDIAN S/N
        print 'S/N values: '
        snrs = n.zeros(len(naspecs)+1)
        newrs = n.zeros(len(rs)+1)
        newrs[:-1] = rs; newrs[-1] = 0.0
        for i in xrange(len(naspecs)):
            if ind in ['NaI', 'NaIsdss']:
                snrs[i] = n.median(naspecs[i].flam[0]/n.sqrt(navars[i].flam[0]))
            elif ind in ['CaT', 'MgI', 'TiO']:
                snrs[i] = n.median(nacatspecs[i].flam[0]/n.sqrt(nacatvars[i].flam[0]))
            elif ind == 'FeH':
                snrs[i] = n.median(fehcspecs[i].flam[0]/n.sqrt(fehcvars[i].flam[0]))
                # if osfitted == 'True':
                #     snrs[i] = n.median(fehcososspecs[ii].flam[0]/n.sqrt(fehcososvars[ii].flam[0]))
                # if fehopt == 'True':
                #     snrs[i] = n.median(fehospecs[ii].flam[0]/n.sqrt(fehovars[ii].flam[0]))
        if gspecs == 'True':
            'Global spec S/Ns:'
            if ind in ['NaI', 'NaIsdss']:
                gsnr = n.median(nagspec.flam[0]/n.sqrt(nagvar.flam[0]))
                print 'Global S/N = ', gsnr
            elif ind in ['CaT', 'MgI', 'TiO']:
                gsnr = n.median(nacatgspec.flam[0]/n.sqrt(nacatgvar.flam[0]))
                print 'Global S/N = ', gsnr
            elif ind == 'FeH':
                gsnr = n.median(fehcgspec.flam[0]/n.sqrt(fehcgvar.flam[0]))
                print 'Global S/N = ', gsnr
            snrs[-1] = gsnr
        print 'Rs = ', newrs
        print 'S/Ns: = ', snrs
        n.savetxt('/Users/zieleniewski/Data/swift/Coma_stellar_pops/'+gal+'/'+gal+'_'+ind+'_'+str(out_sig)[:-2]+'_kms_SNRs.txt',
                  n.column_stack((newrs, snrs)), delimiter=',', fmt='%.3f', header=gal+' '+ind+' '+str(out_sig)[:-2]+'kms SNRs')
        #END CALCULATE S/N VALUES

    #CALCULATE INDEX VALUES
    print 'Index values: '
    inds = n.zeros(len(naspecs))
    indvs = n.zeros(len(naspecs))
    for i in xrange(len(naspecs)):
        if ind in ['NaI', 'NaIsdss']:
            inds[i], indvs[i] = naspecs[i].irindex(naspecs[i].lam[1]-naspecs[i].lam[0], ind, navar, 'vac')
        elif ind in ['CaT', 'MgI']:
            inds[i], indvs[i] = nacatspecs[i].irindex(nacatspecs[i].lam[1]-nacatspecs[i].lam[0], ind, nacatvar, 'vac')
        elif ind == 'TiO':
            inds[i], indvs[i] = nacatspecs[i].TiOindex(nacatspecs[i].lam[1]-nacatspecs[i].lam[0], nacatvar, 'vac')
        elif ind == 'FeH':
            inds[i], indvs[i] = fehcspecs[i].irindex(fehcspecs[i].lam[1]-fehcspecs[i].lam[0], 'FeH', fehcvar, 'vac')
            if osfitted == 'True':
                fehcososspecs[i].irindex(fehcososspecs[i].lam[1]-fehcososspecs[i].lam[0], 'FeH', fehcvar, 'vac')
            if fehopt == 'True':
                fehospecs[i].irindex(fehospecs[i].lam[1]-fehospecs[i].lam[0], 'FeH', fehovar, 'vac')
    if gspecs == 'True':
        if ind in ['NaI', 'NaIsdss']:
            nagspec.irindex(nagspec.lam[1]-nagspec.lam[0], ind, nagvar, 'vac')
        elif ind in ['CaT', 'MgI']:
            nacatgspec.irindex(nacatgspec.lam[1]-nacatgspec.lam[0], ind, nacatgvar, 'vac')
        elif ind == 'TiO':
            nacatgspec.TiOindex(nacatgspec.lam[1]-nacatgspec.lam[0], nacatgvar, 'vac')
        elif ind == 'FeH':
            fehcgspec.irindex(fehcgspec.lam[1]-fehcgspec.lam[0], 'FeH', fehcgvar, 'vac')


    if savevals == 'True':
        n.savetxt('/Users/zieleniewski/Data/swift/Coma_stellar_pops/'+gal+'/'+gal+'_'+ind+'_'+str(out_sig)[:-2]+'_kms_skymask='+skymask+'_index_vals.txt',
                  n.column_stack((inds, indvs)), delimiter=',', fmt='%.5f', header=gal+' '+ind+' '+str(out_sig)[:-2]+'km/s skymask='+skymask)
#------------#



#------------#
def plot_global_spec(gal, ind, out_sig='intrinsic', skymask='False', saveval='True', iterate='False'):

    na, nacat, fehc, fehos, fehosos, feho, skylines_nai, skylines_cat, skylines_mgi, skylines_tio, skylines_feh = get_data(gal)

    #Wavelengths
    lams = nacat['wave']*10000.

    #Sausage spectrum to get index definitions
    sspec = s.spectrum(lam=lams, lamspec=na['skysubtdgal'][:,0]/n.median(na['skysubtdgal'][:,0]))

    if ind == 'NaI':
        start = n.where(lams > 8250.)[0][0]
        end = n.where(lams > 8500.)[0][0]
        fcont = sspec.nacont
        feat = sspec.na

    if ind == 'NaIsdss':
        start = n.where(lams > 8250.)[0][0]
        end = n.where(lams > 8500.)[0][0]
        fcont = sspec.nacont_sdss
        feat = sspec.na_sdss

    if ind == 'CaT':
        start = n.where(lams > 8400.)[0][0]
        end = n.where(lams > 9100.)[0][0]
        fcont = sspec.cacont
        feat = sspec.ca

    if ind == 'MgI':
        start = n.where(lams > 8900.)[0][0]
        end = n.where(lams > 9150.)[0][0]
        fcont = sspec.mgcont
        feat = sspec.mg

    if ind == 'TiO':
        start = n.where(lams > 8950.)[0][0]
        end = n.where(lams > 9150.)[0][0]
        fcont = sspec.tiocont
        feat = sspec.tio

    if ind == 'FeH':
        start = n.where(lams > 9800.)[0][0]
        end = n.where(lams > 10300.)[0][0]
        fcont = sspec.fehcont
        feat = sspec.feh

    print 'ind = ', ind
    print 'Feature = ', feat
    print 'Continuum = ', fcont

    #Cut wavelength range
    indlams = lams[start:end]

    print 'FEAT = ', feat
    print 'FCONT = ', fcont
    indbmid = (fcont[0]+fcont[1])/2.
    indrmid = (fcont[-1]+fcont[-1])/2.

    rs = nacat['rbar']*0.235
    rs[0] = 0.0
    print 'RS = ', rs


    #Redshift correct spectra
    z = nacat['param'][0,0]*1000./sc.c
    sig = nacat['param'][1,0]
    print 'z = ', z
    print 'sig = ', sig, 'km/s'
    nacatspec = ssz.z_correct(n.column_stack((indlams,nacat['skysubtdggal'][start:end]/n.median(nacat['skysubtdggal'][start:end]))), z)
    nacatmspec = ssz.z_correct(n.column_stack((indlams,nacat['bestfitgMODEL'][start:end]/n.median(nacat['bestfitgMODEL'][start:end]))), z)
    nacatvar = ssz.z_correct(n.column_stack((indlams,nacat['skysubtdggerr'][start:end]**2/n.median(nacat['skysubtdggal'][start:end]**2))), z)

    naspec = ssz.z_correct(n.column_stack((indlams,na['skysubtdggal'][start:end]/n.median(na['skysubtdggal'][start:end]))), z)
    namspec = ssz.z_correct(n.column_stack((indlams,na['bestfitgMODEL'][start:end]/n.median(na['bestfitgMODEL'][start:end]))), z)
    navar = ssz.z_correct(n.column_stack((indlams,na['skysubtdggerr'][start:end]**2/n.median(na['skysubtdggal'][start:end]**2))), z)

    fehcspec = ssz.z_correct(n.column_stack((indlams,fehc['skysubtdggal'][start:end]/n.median(fehc['skysubtdggal'][start:end]))), z)
    fehcmspec = ssz.z_correct(n.column_stack((indlams,fehc['bestfitgMODEL'][start:end]/n.median(fehc['bestfitgMODEL'][start:end]))), z)
    fehcvar = ssz.z_correct(n.column_stack((indlams,fehc['skysubtdggerr'][start:end]**2/n.median(fehc['skysubtdggal'][start:end]**2))), z)

    fehospec = ssz.z_correct(n.column_stack((indlams,feho['skysubtdggal'][start:end]/n.median(feho['skysubtdggal'][start:end]))), z)
    fehomspec = ssz.z_correct(n.column_stack((indlams,feho['bestfitgMODEL'][start:end]/n.median(feho['bestfitgMODEL'][start:end]))), z)
    fehovar = ssz.z_correct(n.column_stack((indlams,feho['skysubtdggerr'][start:end]**2/n.median(feho['skysubtdggal'][start:end]**2))), z)

    #####Iterate over sky lines, mask, smooth, replace masked region with smoothed pixels x5
    if ind in ['NaI', 'NaIsdss', 'CaT', 'MgI', 'FeH'] and skymask=='True':
        print 'Finding and masking skylines'

        #Find and mask sky line regions
        skyline_pos = []
        if ind == 'NaI' or ind == 'NaIsdss':
            for j in xrange(skylines_nai.shape[0]):
                sline_start = n.where(naspec[:,0]>skylines_nai[j,0])[0][0]
                sline_end = n.where(naspec[:,0]>skylines_nai[j,1])[0][0]
                sline_interp = si.interp1d([naspec[sline_start,0],naspec[sline_end,0]],[naspec[sline_start,1],naspec[sline_end,1]])
                maskspec = sline_interp(naspec[sline_start:sline_end,0])
                if gal=='GMP3329' and j==1:#replace with modelspec for red continuum region
                    naspec[sline_start:sline_end,1] = namspec[sline_start:sline_end,1]
                if gal=='GMP3329' and j==0 and i==(len(na['binflux'][0,:])-1):#interpolate over bad region for outermost spectrum
                    naspec[sline_start:sline_end,1] = maskspec
                elif gal!='GMP3329':#replace with modelspec for all other galaxies
                    naspec[sline_start:sline_end,1] = namspec[sline_start:sline_end,1]
                skyline_positions = n.arange(sline_start,sline_end+1,1)
                for jj in xrange(len(skyline_positions)):
                    skyline_pos.append(skyline_positions[jj])
        if ind == 'CaT':
            for j in xrange(skylines_cat.shape[0]):
                sline_start = n.where(nacatspec[:,0]>skylines_cat[j,0])[0][0]
                sline_end = n.where(nacatspec[:,0]>skylines_cat[j,1])[0][0]
                sline_interp = si.interp1d([nacatspec[sline_start,0],nacatspec[sline_end,0]],[nacatspec[sline_start,1],nacatspec[sline_end,1]])
                maskspec = sline_interp(nacatspec[sline_start:sline_end,0])
                #nacatspec[sline_start:sline_end,1] = maskspec
                nacatspec[sline_start:sline_end,1] = nacatmspec[sline_start:sline_end,1]
                skyline_positions = n.arange(sline_start,sline_end+1,1)
                for jj in xrange(len(skyline_positions)):
                    skyline_pos.append(skyline_positions[jj])
        elif ind == 'MgI':
            for j in xrange(skylines_mgi.shape[0]):
                sline_start = n.where(nacatspec[:,0]>skylines_mgi[j,0])[0][0]
                sline_end = n.where(nacatspec[:,0]>skylines_mgi[j,1])[0][0]
                sline_interp = si.interp1d([nacatspec[sline_start,0],nacatspec[sline_end,0]],[nacatspec[sline_start,1],nacatspec[sline_end,1]])
                maskspec = sline_interp(nacatspec[sline_start:sline_end,0])
                #nacatspec[sline_start:sline_end,1] = maskspec
                nacatspec[sline_start:sline_end,1] = nacatmspec[sline_start:sline_end,1]
                skyline_positions = n.arange(sline_start,sline_end+1,1)
                for jj in xrange(len(skyline_positions)):
                    skyline_pos.append(skyline_positions[jj])
        elif ind == 'TiO':
            for j in xrange(skylines_tio.shape[0]):
                sline_start = n.where(nacatspec[:,0]>skylines_tio[j,0])[0][0]
                sline_end = n.where(nacatspec[:,0]>skylines_tio[j,1])[0][0]
                sline_interp = si.interp1d([nacatspec[sline_start,0],nacatspec[sline_end,0]],[nacatspec[sline_start,1],nacatspec[sline_end,1]])
                maskspec = sline_interp(nacatspec[sline_start:sline_end,0])
                #nacatspec[sline_start:sline_end,1] = maskspec
                nacatspec[sline_start:sline_end,1] = nacatmspec[sline_start:sline_end,1]
                skyline_positions = n.arange(sline_start,sline_end+1,1)
                for jj in xrange(len(skyline_positions)):
                    skyline_pos.append(skyline_positions[jj])
        elif ind == 'FeH':
            for j in xrange(skylines_feh.shape[0]):
                sline_start = n.where(fehcspec[:,0]>skylines_feh[j,0])[0][0]
                sline_end = n.where(fehcspec[:,0]>skylines_feh[j,1])[0][0]
                sline_interp = si.interp1d([fehcspec[sline_start,0],fehcspec[sline_end,0]],[fehcspec[sline_start,1],fehcspec[sline_end,1]])
                maskspec = sline_interp(fehcspec[sline_start:sline_end,0])
                #fehcspec[sline_start:sline_end,1] = maskspec
                fehcspec[sline_start:sline_end,1] = fehcmspec[sline_start:sline_end,1]
                fehospec[sline_start:sline_end,1] = fehomspec[sline_start:sline_end,1]
                skyline_positions = n.arange(sline_start,sline_end+1,1)
                for jj in xrange(len(skyline_positions)):
                    skyline_pos.append(skyline_positions[jj])

        if iterate=='True':
            print 'Entering sky line iteration loop'
            for k in xrange(5):
                print 'Loop = ', k

                #Load as spectrum class
                nacatsubspec = s.spectrum(lam=nacatspec[:,0], lamspec=nacatspec[:,1])
                nasubspec = s.spectrum(lam=naspec[:,0], lamspec=naspec[:,1])
                fehcsubspec = s.spectrum(lam=fehcspec[:,0], lamspec=fehcspec[:,1])
                fehosubspec = s.spectrum(lam=fehospec[:,0], lamspec=fehospec[:,1])

                #Convolve to common velocity dispersion
                if out_sig != 'intrinsic':
                    out_sig = float(out_sig)
                    convsig = n.sqrt(out_sig**2 - sig**2)
                    loglambs = n.linspace(n.log10(nasubspec.lam[0]), n.log10(nasubspec.lam[-1]), len(nasubspec.lam))
                    veldisp = (10**(loglambs[1]-loglambs[0]) - 1)*sc.c/1000. #in km/s
                    gauss = ssz.Gauss(convsig, veldisp)
                    nacatsubspec.gaussConvolve(0.0, convsig, losvd=gauss[:,1])
                    linlambs = n.linspace(nacatsubspec.loglam[0], nacatsubspec.loglam[-1], len(nacatsubspec.lam))
                    intspec = si.interp1d(nacatsubspec.loglam, nacatsubspec.conflam)
                    nanewflux = intspec(linlambs)

                    nasubspec.gaussConvolve(0.0, convsig, losvd=gauss[:,1])
                    linlambs = n.linspace(nasubspec.loglam[0], nasubspec.loglam[-1], len(nasubspec.lam))
                    intspec = si.interp1d(nasubspec.loglam, nasubspec.conflam)
                    nanewflux = intspec(linlambs)

                    fehcsubspec.gaussConvolve(0.0, convsig, losvd=gauss[:,1])
                    linlambs = n.linspace(fehcsubspec.loglam[0], fehcsubspec.loglam[-1], len(fehcsubspec.lam))
                    intspec = si.interp1d(fehcsubspec.loglam, fehcsubspec.conflam)
                    fehcnewflux = intspec(linlambs)

                    fehosubspec.gaussConvolve(0.0, convsig, losvd=gauss[:,1])
                    linlambs = n.linspace(fehosubspec.loglam[0], fehosubspec.loglam[-1], len(fehosubspec.lam))
                    intspec = si.interp1d(fehosubspec.loglam, fehosubspec.conflam)
                    fehonewflux = intspec(linlambs)

                    #Replace masked pixels in original spectrum
                    nacatspec[skyline_pos,1] = nacatnewflux[0][skyline_pos]
                    naspec[skyline_pos,1] = nanewflux[0][skyline_pos]
                    fehcspec[skyline_pos,1] = fehcnewflux[0][skyline_pos]
                    fehospec[skyline_pos,1] = fehonewflux[0][skyline_pos]
                else:
                    pass

    #####End iteration

    #Load as spectrum class
    nacatspec = s.spectrum(lam=nacatspec[:,0], lamspec=nacatspec[:,1])
    nacatvar = s.spectrum(lam=nacatvar[:,0], lamspec=nacatvar[:,1])
    naspec = s.spectrum(lam=naspec[:,0], lamspec=naspec[:,1])
    navar = s.spectrum(lam=navar[:,0], lamspec=navar[:,1])
    fehcspec = s.spectrum(lam=fehcspec[:,0], lamspec=fehcspec[:,1])
    fehcvar = s.spectrum(lam=fehcvar[:,0], lamspec=fehcvar[:,1])
    fehospec = s.spectrum(lam=fehospec[:,0], lamspec=fehospec[:,1])
    fehovar = s.spectrum(lam=fehovar[:,0], lamspec=fehovar[:,1])

    ###Convolve to common velocity dispersion
    if out_sig != 'intrinsic':
        out_sig = float(out_sig)

        #Instrument sigma = 50.4 km/s for NaCaT; 63.3 km/s for FeH)
        if ind == 'FeH':
            instsig = 63.3
        else:
            instsig = 50.4
        convsig = n.sqrt(out_sig**2 - (na['param'][1,0]**2+instsig**2))

        convsig = n.sqrt(out_sig**2 - sig**2)
        loglambs = n.linspace(n.log10(naspec.lam[0]), n.log10(naspec.lam[-1]), len(naspec.lam))
        veldisp = (10**(loglambs[1]-loglambs[0]) - 1)*sc.c/1000. #in km/s
        gauss = ssz.Gauss(convsig, veldisp)

        nacatspec.gaussConvolve(0.0, convsig, losvd=gauss[:,1])
        nacatlinlambs = n.linspace(nacatspec.loglam[0], nacatspec.loglam[-1], len(nacatspec.lam))
        intspec = si.interp1d(nacatspec.loglam, nacatspec.conflam)
        nacatnewflux = intspec(nacatlinlambs)
        nacatvar.gaussConvolve(0.0, convsig, losvd=gauss[:,1])
        nacatvlinlambs = n.linspace(nacatvar.loglam[0], nacatvar.loglam[-1], len(nacatvar.lam))
        intspec = si.interp1d(nacatvar.loglam, nacatvar.conflam)
        nacatnewvar = intspec(nacatvlinlambs)

        naspec.gaussConvolve(0.0, convsig, losvd=gauss[:,1])
        nalinlambs = n.linspace(naspec.loglam[0], naspec.loglam[-1], len(naspec.lam))
        intspec = si.interp1d(naspec.loglam, naspec.conflam)
        nanewflux = intspec(nalinlambs)
        navar.gaussConvolve(0.0, convsig, losvd=gauss[:,1])
        navlinlambs = n.linspace(navar.loglam[0], navar.loglam[-1], len(navar.lam))
        intspec = si.interp1d(navar.loglam, navar.conflam)
        nanewvar = intspec(navlinlambs)

        fehcspec.gaussConvolve(0.0, convsig, losvd=gauss[:,1])
        fehclinlambs = n.linspace(fehcspec.loglam[0], fehcspec.loglam[-1], len(fehcspec.lam))
        intspec = si.interp1d(fehcspec.loglam, fehcspec.conflam)
        fehcnewflux = intspec(fehclinlambs)
        fehcvar.gaussConvolve(0.0, convsig, losvd=gauss[:,1])
        fehcvlinlambs = n.linspace(fehcvar.loglam[0], fehcvar.loglam[-1], len(fehcvar.lam))
        intspec = si.interp1d(fehcvar.loglam, fehcvar.conflam)
        fehcnewvar = intspec(fehcvlinlambs)

        fehospec.gaussConvolve(0.0, convsig, losvd=gauss[:,1])
        feholinlambs = n.linspace(fehospec.loglam[0], fehospec.loglam[-1], len(fehospec.lam))
        intspec = si.interp1d(fehospec.loglam, fehospec.conflam)
        fehonewflux = intspec(feholinlambs)
        fehovar.gaussConvolve(0.0, convsig, losvd=gauss[:,1])
        fehovlinlambs = n.linspace(fehovar.loglam[0], fehovar.loglam[-1], len(fehovar.lam))
        intspec = si.interp1d(fehovar.loglam, fehovar.conflam)
        fehonewvar = intspec(fehovlinlambs)

        nacatspec = s.spectrum(lam=nacatlinlambs, lamspec=nacatnewflux)
        nacatvar = s.spectrum(lam=nacatvlinlambs, lamspec=nacatnewvar)
        naspec = s.spectrum(lam=nalinlambs, lamspec=nanewflux)
        navar = s.spectrum(lam=navlinlambs, lamspec=nanewvar)
        fehcspec = s.spectrum(lam=fehclinlambs, lamspec=fehcnewflux)
        fehcvar = s.spectrum(lam=fehcvlinlambs, lamspec=fehcnewvar)
        fehospec = s.spectrum(lam=feholinlambs, lamspec=fehonewflux)
        fehovar = s.spectrum(lam=fehovlinlambs, lamspec=fehonewvar)

    else:
        pass
    ###

    #Normalise spectrum by pseudo-continuum
    nacatspec_bcont = n.mean(nacatspec.flam[0][n.where(nacatspec.lam > fcont[0])[0][0]:n.where(nacatspec.lam > fcont[1])[0][0]])
    nacatspec_rcont = n.mean(nacatspec.flam[0][n.where(nacatspec.lam > fcont[-2])[0][0]:n.where(nacatspec.lam > fcont[-1])[0][0]])
    nacatspec_contfit = n.polyfit([indbmid, indrmid], [nacatspec_bcont, nacatspec_rcont], 1)
    nacatspec_contspec = n.column_stack((nacatspec.lam, nacatspec.lam*nacatspec_contfit[0] + nacatspec_contfit[1]))
    #nacatdivspec = s.spectrum(lam=naspec.lam, lamspec=(naspec.flam[0]/naspec_contspec[:,1]))
    nacatdivspec = nacatspec

    naspec_bcont = n.mean(naspec.flam[0][n.where(naspec.lam > fcont[0])[0][0]:n.where(naspec.lam > fcont[1])[0][0]])
    naspec_rcont = n.mean(naspec.flam[0][n.where(naspec.lam > fcont[-2])[0][0]:n.where(naspec.lam > fcont[-1])[0][0]])
    naspec_contfit = n.polyfit([indbmid, indrmid], [naspec_bcont, naspec_rcont], 1)
    naspec_contspec = n.column_stack((naspec.lam, naspec.lam*naspec_contfit[0] + naspec_contfit[1]))
    #nadivspec = s.spectrum(lam=naspec.lam, lamspec=(naspec.flam[0]/naspec_contspec[:,1]))
    nadivspec = naspec

    fehcspec_bcont = n.mean(fehcspec.flam[0][n.where(fehcspec.lam > fcont[0])[0][0]:n.where(fehcspec.lam > fcont[1])[0][0]])
    fehcspec_rcont = n.mean(fehcspec.flam[0][n.where(fehcspec.lam > fcont[-2])[0][0]:n.where(fehcspec.lam > fcont[-1])[0][0]])
    fehcspec_contfit = n.polyfit([indbmid, indrmid], [fehcspec_bcont, fehcspec_rcont], 1)
    fehcspec_contspec = n.column_stack((fehcspec.lam, fehcspec.lam*fehcspec_contfit[0] + fehcspec_contfit[1]))
    fehcdivspec = s.spectrum(lam=fehcspec.lam, lamspec=(fehcspec.flam[0]/fehcspec_contspec[:,1]))
    fehcdivvar = s.spectrum(lam=fehcvar.lam, lamspec=(fehcvar.flam[0]/(fehcspec_contspec[:,1])**2))
    #fehcdivspec = fehcspec

    fehospec_bcont = n.mean(fehospec.flam[0][n.where(fehospec.lam > fcont[0])[0][0]:n.where(fehospec.lam > fcont[1])[0][0]])
    fehospec_rcont = n.mean(fehospec.flam[0][n.where(fehospec.lam > fcont[-2])[0][0]:n.where(fehospec.lam > fcont[-1])[0][0]])
    fehospec_contfit = n.polyfit([indbmid, indrmid], [fehospec_bcont, fehospec_rcont], 1)
    fehospec_contspec = n.column_stack((fehospec.lam, fehospec.lam*fehospec_contfit[0] + fehospec_contfit[1]))
    fehodivspec = s.spectrum(lam=fehospec.lam, lamspec=(fehospec.flam[0]/fehospec_contspec[:,1]))
    fehodivvar = s.spectrum(lam=fehovar.lam, lamspec=(fehovar.flam[0]/(fehospec_contspec[:,1])**2))
    #fehodivspec = fehospec


    P.rcParams['font.family'] = 'Times New Roman'
    P.rcParams.update({'font.size': 20})
    P.figure(figsize=(10,8), dpi=72)

    mycmap = P.get_cmap('copper')

    if ind in ['CaT', 'MgI', 'TiO']:
        P.plot(nacatdivspec.lam, nacatdivspec.flam[0], color='k', ls='-',
               lw=2.5, label='Global Spec')

    elif ind in ['NaI', 'NaIsdss']:
        P.plot(nadivspec.lam, nadivspec.flam[0], color='k', ls='-',
               lw=2.5, label='Global Spec')

    elif ind == 'FeH':
        P.plot(fehcdivspec.lam, fehcdivspec.flam[0], color='k', ls='-',
               lw=2.5, label='Global Spec')
        P.plot(fehodivspec.lam, fehodivspec.flam[0], color='r', ls='-',
               lw=2.5, label='Global Spec')

    #sm = P.cm.ScalarMappable(cmap=mycmap, norm=P.Normalize(vmin=rs[0], vmax=rs[-1]))
    # fake up the array of the scalar mappable. Urgh...
    #sm._A = []
    #sm.set_array(rs)
    #clb = P.colorbar(sm)
    #clb = P.colorbar(plot[0])
    #clb.set_label('r [arcsec]', fontsize=tfsize)
    if ind == 'NaI' or ind == 'NaIsdss':
        P.legend(loc='best', frameon=False, borderaxespad=0.5, labelspacing=0.18, handlelength=2.2)

    P.plot([fcont[0],fcont[0]], [0.5,1.5], 'k-',
           [fcont[1],fcont[1]], [0.5,1.5], 'k-',
           [fcont[-2],fcont[-2]], [0.5,1.5], 'k-',
           [fcont[-1],fcont[-1]], [0.5,1.5], 'k-', alpha=0.5, lw=2.2)
    P.plot([feat[0],feat[0]], [0.5,1.5], 'r-', alpha=0.5, zorder=10, lw=2.2)
    P.plot([feat[1],feat[1]], [0.5,1.5], 'r-', alpha=0.5, zorder=10, lw=2.2)


#    P.xlim([fcont[0]-10,fcont[3]+10])
#    P.ylim([0.97,1.03])

    #For NaI
    if ind == 'NaI' or ind == 'NaIsdss':
        P.xlim([8130., 8250.])
        P.ylim([0.94,1.06])

    #For NaI
    elif ind == 'CaT':
        P.xlim([8400., 8800.])
        P.ylim([0.85,1.05])

    #For FeH
    elif ind == 'FeH':
        P.xlim([9840., 9980.])
        P.ylim([0.93,1.05])

    #For MgI
    elif ind == 'MgI':
        P.xlim([8765., 8865.])
        P.ylim([0.93,1.13])

    #For TiO
    elif ind == 'TiO':
        P.xlim([8810., 8910.])
        P.ylim([0.93,1.13])

    P.xlabel('$\lambda$ [$\AA$]', fontsize=tfsize)
    P.ylabel('Normalised flux', fontsize=tfsize)

    P.xticks(size=lfsize)
    P.yticks(size=lfsize)

    P.tick_params(axis='both', which='major', direction='in', length=majtsize, width=twidth)
    P.tick_params(axis='both', which='major', direction='in', length=mintsize, width=twidth)

    if out_sig != 'intrinsic':
        P.annotate('$\sigma$ = '+str(n.round(out_sig,1))+'km s$^{-1}$', xy=(feat[0]+2,1.02), textcoords='axes fraction',
                   xytext=(0.35,0.8), fontsize=lfsize, color='k')

    P.tight_layout()
    P.show()

    if ind in ['NaI', 'NaIsdss']:
        inds, indvs = nadivspec.irindex(nadivspec.lam[1]-nadivspec.lam[0], ind, navar, 'vac')
    if ind in ['CaT', 'MgI']:
        inds, indvs = nacatdivspec.irindex(nacatdivspec.lam[1]-nacatdivspec.lam[0], ind, navar, 'vac')
    elif ind == 'TiO':
        inds, indvs = nadivspec.TiOindex(nadivspec.lam[1]-nadivspec.lam[0], navar, 'vac')
    elif ind == 'FeH':
        inds, indvs = fehcdivspec.irindex(fehcdivspec.lam[1]-fehcdivspec.lam[0], 'FeH', fehcdivvar, 'vac')
        fehodivspec.irindex(fehodivspec.lam[1]-fehodivspec.lam[0], 'FeH', fehodivvar, 'vac')

    if saveval == 'True':
        n.savetxt('/Users/zieleniewski/Data/swift/Coma/Coma_stellar_pops/'+gal+'/'+gal+'_'+ind+'_'+str(out_sig)[:-2]+'_kms_skymask='+skymask+'_global_index_val.txt',
                  n.column_stack((inds, indvs)), delimiter=',', fmt='%.5f', header=gal+' '+ind+' '+str(out_sig)[:-2]+'km/s skymask='+skymask+'_global_spec')
#------------#




#------------#
def colour_map(gal):
    '''Function that creates SDSS g-r colour map
    using images from the Montage software. Maps are
    generated to match the binning procedure
    by pPXF on the galaxy datacubes.

    '''

    #Galaxy [GMP2921, GMP3329, GMP4928, GMP3367]

    #Load galaxy IDL file
    if gal == 'GMP2921':
        data = readsav('/Users/zieleniewski/Data/swift/Coma_stellar_pops/'+gal+'/NaCaT_v2/Conroy_0.79_0.905/ppxfSECT_SINFONI.idl')
        cent = n.array([24,44])

    if gal == 'GMP3329':
        data = readsav('/Users/zieleniewski/Data/swift/Coma_stellar_pops/'+gal+'/NaCaT_v3/Conroy_0.80_0.92/ppxfSECT_SINFONI.idl')
        cent = n.array([22,44])

    if gal == 'GMP4928':
        data = readsav('/Users/zieleniewski/Data/swift/Coma_stellar_pops/'+gal+'/NaCaT_v4/Conroy_0.80_0.92/ppxfSECT_SINFONI.idl')
        cent = n.array([22,55])

    if gal == 'GMP3367':
        data = readsav('/Users/zieleniewski/Data/swift/Coma_stellar_pops/'+gal+'/NaCaT/Conroy_0.79_0.905/ppxfSECT_SINFONI.idl')
        #cent = n.array([22,44])

    #Create map of bins
    image = n.copy(data["fluxav"])*0.0
    for q in range(data["numnotnans"]):
        image[data["one2two"][1,q], data["one2two"][0,q]] = data["binnum"][q]

    print 'MIN = ', n.min(image)
    print 'MAX = ', n.max(image)

    #Load up SDSS images
    sdss_g, ghead = p.getdata('/Users/zieleniewski/Data/Coma_colour_data/'+gal+'_g.fits', header=True)
    sdss_r, rhead = p.getdata('/Users/zieleniewski/Data/Coma_colour_data/'+gal+'_r.fits', header=True)

    #Image spaxel scale in mas
    spax = ghead['CDELT1']
    print 'SPAX = ', spax, ' mas'

    import frebin as fre
    sys.path.append('/Users/zieleniewski/Documents/Oxford/D.Phil/Simulator/hsim/src/modules/')
    import Gaussians as gs
    import psf_convolution as pc

    #Frebin images up to 235mas
    sdss_g_235, ghead = fre.spaxel_scale(sdss_g, ghead, (235.,235.))
    sdss_r_235, rhead = fre.spaxel_scale(sdss_r, rhead, (235.,235.))


    #Gaussian representing 1 arcsec FWHM seeing.
    see = 1000. #[mas]
    psf = gs.Gauss2D(sdss_g_235.shape[0], fwhm=(see/235.))

    sdss_g_235_conv = pc.psf_convolve(sdss_g_235, psf)
    sdss_r_235_conv = pc.psf_convolve(sdss_r_235, psf)

    #Cut down to SWIFT FoV
    y, x = sdss_g_235_conv.shape
    print 'IMAGE SIZE = ', y, x

    yy, xx = image.shape
    print 'SWIFT SIZE = ', yy, xx

    sdss_g_235_conv_swift = sdss_g_235_conv[y/2-cent[0]:y/2+(yy-cent[0]), x/2-cent[1]:x/2+(xx-cent[1])]
    sdss_r_235_conv_swift = sdss_r_235_conv[y/2-cent[0]:y/2+(yy-cent[0]), x/2-cent[1]:x/2+(xx-cent[1])]

    print 'NEW GIMG SIZE = ', sdss_g_235_conv_swift.shape
    print 'NEW RIMG SIZE = ', sdss_r_235_conv_swift.shape


    #Arrange colour arrays
    colimg = n.copy(data['fluxav'])*0.0
    for q in range(data['numnotnans']):
        colimg[data['one2two'][1,q], data['one2two'][0,q]] = data['binnum'][q]
    gim = n.zeros_like(colimg)
    rim = n.zeros_like(colimg)
    for i in xrange(data['binflux'].shape[1]):
        args = n.where(colimg == i)
        num_g = n.sum(sdss_g_235_conv_swift[args])
        mag_g = -2.5*n.log10(num_g)+ghead['MAGZP']
        gim[args] = mag_g
        num_r = n.sum(sdss_r_235_conv_swift[args])
        mag_r = -2.5*n.log10(num_r)+rhead['MAGZP']
        rim[args] = mag_r

    colmap = gim-rim

    P.rcParams.update({'font.size': 20})
    P.figure(figsize=(12,6),dpi=72)
    gext = n.array([-((image.shape[1]/2.)*235./1000.),((image.shape[1]/2.)*235./1000.),
                    -((image.shape[0]/2.)*235./1000.),((image.shape[0]/2.)*235./1000.)])
    P.imshow(image, interpolation='nearest', aspect='equal', extent=gext)
    P.xlabel('kpc', fontsize=lfsize)
    P.ylabel('kpc', fontsize=lfsize)
    print data['binnum'].min()
    print data['binnum'].max()
    P.colorbar()
    P.tight_layout()
    P.show()

    P.rcParams.update({'font.size': 20})
    P.figure(figsize=(12,6),dpi=72)
    gext = n.array([-((image.shape[1]/2.)*235./1000.),((image.shape[1]/2.)*235./1000.),
                    -((image.shape[0]/2.)*235./1000.),((image.shape[0]/2.)*235./1000.)])
    P.imshow(colmap, interpolation='nearest', aspect='equal', extent=gext, vmin=0.45, vmax=0.7)
    P.xlabel('kpc', fontsize=lfsize)
    P.ylabel('kpc', fontsize=lfsize)
    P.colorbar()
    P.tight_layout()
    P.show()
#------------#



#------------#
def sdss_specs_ind(gal):
    '''Function to measure indices from SDSS
    Coma spectra. Currently measures index at sigma
    of central resolved bin for each galaxy
    (approximately 3 arcsec diameter)
    From Price+2011 sigmas [km/s]:
    GMP2921 - 386.9
    GMP3329 - 278.4
    GMP4928 - 273.5
    GMP3367 - 173.4

    Inputs:

        - gal: galaxy
        - ind: index
        - insig: input velocity dispersion [km/s]

    Outputs:

        - indval
        - indvar

    '''
    fpath = '/Users/zieleniewski/Data/SDSS_Coma_spectra/'
    if gal == 'GMP2921':
        sdss = shs.load_spec(fpath+'spec-2241-54169-0542.fits')[0]
        na, nacat, fehc, fehcosos, feho, skylines_nai, skylines_cat, skylines_mgi, skylines_tio, skylines_feh = get_data(gal)
        insig = 386.9
        sig = 400.0


    if gal == 'GMP3329':
        sdss = shs.load_spec(fpath+'spec-2240-53823-0585.fits')[0]
        na, nacat, fehc, fehcosos, feho, skylines_nai, skylines_cat, skylines_mgi, skylines_tio, skylines_feh = get_data(gal)
        insig = 278.4
        sig = 270.0

    if gal == 'GMP4928':
        sdss = shs.load_spec(fpath+'spec-2240-53823-0157.fits')[0]
        na, nacat, fehc, fehcosos, feho, skylines_nai, skylines_cat, skylines_mgi, skylines_tio, skylines_feh = get_data(gal)
        insig = 273.5
        sig = 270.0

    if gal == 'GMP3367':
        sdss = shs.load_spec(fpath+'spec-2240-53823-0562.fits')[0]
        na, nacat, fehc, fehcosos, feho, skylines_nai, skylines_cat, skylines_mgi, skylines_tio, skylines_feh = get_data(gal)
        insig = 173.4
        sig = 200.0

    z = nacat['param'][0,0]*1000./sc.c
    # sig = nacat['param'][1,0]
    print 'SWIFT sig = ', n.round(sig, 1), ' km/s'

    if sig > insig:
        'SWIFT SIG > SDSS SIG: convolving up to SWIFT SIG (%.1f)' % sig
        convsig = n.sqrt(sig**2 - insig**2)
        fsig = sig
    elif insig > sig:
        'SDSS SIG > SWIFT SIG: convolving up to SDSS SIG (%.1f)' % insig
        convsig = n.sqrt(insig**2 - sig**2)
        fsig = insig

    zspec = ssz.z_correct(n.column_stack((sdss[:,0],sdss[:,1]/n.median(sdss[:,1]))), z)
    zvar = ssz.z_correct(n.column_stack((sdss[:,0],sdss[:,2]/(n.median(sdss[:,1]))**2)), z)
    sspec = s.spectrum(lam=zspec[:,0], lamspec=zspec[:,1])
    svar = s.spectrum(lam=zvar[:,0], lamspec=zvar[:,1])

    # loglambs = n.linspace(n.log10(sspec.lam[0]), n.log10(sspec.lam[-1]), len(sspec.lam))
    # veldisp = (10**(loglambs[1]-loglambs[0]) - 1)*sc.c/1000. #in km/s
    # gauss = ssz.Gauss(convsig, veldisp)
    # sspec, svar = convolve_spec(sspec, convsig, gauss, svar)

    nadind, nadindvar = sspec.irindex(sspec.lam[1]-sspec.lam[0], 'NaD', svar, 'vac')
    naiind, naiindvar = sspec.irindex(sspec.lam[1]-sspec.lam[0], 'NaIsdss', svar, 'vac')
    catind, catindvar = sspec.irindex(sspec.lam[1]-sspec.lam[0], 'CaT', svar, 'vac')
    mgiind, mgiindvar = sspec.irindex(sspec.lam[1]-sspec.lam[0], 'MgI', svar, 'vac')
    tioind, tioindvar = sspec.TiOindex(sspec.lam[1]-sspec.lam[0], svar, 'vac')

    print ''
    print 'INDICES at %.1f' % fsig
    print ''
    print 'NaD = %.3f pm %.3f Angstroms' % (nadind, n.sqrt(nadindvar))
    print 'NaIsdss = %.3f pm %.3f Angstroms' % (naiind, n.sqrt(naiindvar))
    print 'CaT = %.3f pm %.3f Angstroms' % (catind, n.sqrt(catindvar))
    print 'MgI = %.3f pm %.3f Angstroms' % (mgiind, n.sqrt(mgiindvar))
    print 'TiO = %.4f pm %.4f ' % (tioind, n.sqrt(tioindvar))


#------------#









if __name__=="__main__":
    import sys
    print len(sys.argv)
    if len(sys.argv) < 2:
        print 'Choose from:'
        print ' '
        print 'indices [galaxy, sigma, [save], [plot], [skymask], [FeH_opt], [osfitted], [gspecs], [corrfac] [cfacval]]'
        print 'allspecs [galaxy, index, sigma, [mask_skylines], [savevals], [fehopt], [osfitted], [gspecs], [S/N specs]]'
        print 'globalspec [galaxy, index, sigma, [mask_skylines], [saveval]]'
        print 'colmap [galaxy]'
        print 'allgals [dir2921, dir3329, dir4928, dir3367, [savefig], [models] [global_spec_files]]'
        print 'corrfac [galfile1, galgfile1, galfile2, galgfile2, [save]]'
        print 'sdss [galaxy]'
        print ' '
        sys.exit()

    if sys.argv[1] == 'indices':
        if len(sys.argv)==4:
            indices_v_radius(sys.argv[2], sys.argv[3])
        elif len(sys.argv)==5:
            indices_v_radius(sys.argv[2], sys.argv[3], sys.argv[4])
        elif len(sys.argv)==6:
            indices_v_radius(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
        elif len(sys.argv)==7:
            indices_v_radius(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6])
        elif len(sys.argv)==8:
            indices_v_radius(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7])
        elif len(sys.argv)==9:
            indices_v_radius(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7], sys.argv[8])
        elif len(sys.argv)==10:
            indices_v_radius(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7], sys.argv[8], sys.argv[9])
        elif len(sys.argv)==11:
            indices_v_radius(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7], sys.argv[8], sys.argv[9], sys.argv[10])
        elif len(sys.argv)==12:
            indices_v_radius(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7], sys.argv[8], sys.argv[9], sys.argv[10], sys.argv[11])

    elif sys.argv[1] == 'globalspec':
        if len(sys.argv)==5:
            plot_global_spec(sys.argv[2], sys.argv[3], sys.argv[4])
        elif len(sys.argv)==6:
            plot_global_spec(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
        elif len(sys.argv)==7:
            plot_global_spec(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6])

    elif sys.argv[1] == 'colmap':
        colour_map(sys.argv[2])

    elif sys.argv[1] == 'allgals':
        if len(sys.argv)==6:
            indices_v_radius_allgals(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
        elif len(sys.argv)==7:
            indices_v_radius_allgals(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6])
        elif len(sys.argv)==8:
            indices_v_radius_allgals(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7])
        elif len(sys.argv)==9:
            indices_v_radius_allgals(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7], [sys.argv[8]])
        elif len(sys.argv)>9:
            indices_v_radius_allgals(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7], sys.argv[8:])

    elif sys.argv[1] == 'corrfac':
        if len(sys.argv)==6:
            indices_v_radius_compare_indices(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
        elif len(sys.argv)==7:
            indices_v_radius_compare_indices(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6])

    elif sys.argv[1] == 'allspecs':
        if len(sys.argv)==5:
            plot_all_specs(sys.argv[2], sys.argv[3], sys.argv[4])
        elif len(sys.argv)==6:
            plot_all_specs(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
        elif len(sys.argv)==7:
            plot_all_specs(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6])
        elif len(sys.argv)==8:
            plot_all_specs(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7])
        elif len(sys.argv)==9:
            plot_all_specs(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7], sys.argv[8])
        elif len(sys.argv)==10:
            plot_all_specs(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7], sys.argv[8], sys.argv[9])
        elif len(sys.argv)==11:
            plot_all_specs(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7], sys.argv[8], sys.argv[9], sys.argv[10])

    elif sys.argv[1] == 'sdss':
        sdss_specs_ind(sys.argv[2])
