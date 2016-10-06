'''Code to calculate and plot Index v Colour plots
for V12 model SEDs.

Written by Simon Zieleniewski

Created 12-06-14

Last updated 01-01-16
'''

import sys
import spectools as s
import spectools_SZ as ssz
import CvD12tools as c12
import numpy as n
import pylab as P
import matplotlib
import measure_real_indices as mri


####FONTS
tfsize = 28
lfsize = 24
gensize = 20
majtsize = 14
mintsize = 6
twidth = 2
msize = 10
####


def indexindexplot(ind1, ind2, imf, comdisp, data=None):
    comdisp = float(comdisp)
    #Index 1
    #Varying IMF at 14 Gyr
    spec14_0p3i1 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MIUSCAT_Padova00_'+imf+'_0.30_fits/I'+imf+'0.30Zp0.00T14.1254_iPp0.00_baseFe.fits',
                                            varspecf=None, kins='V15', out_sig=comdisp, index=ind1)
    spec14_0p8i1 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MIUSCAT_Padova00_'+imf+'_0.80_fits/I'+imf+'0.80Zp0.00T14.1254_iPp0.00_baseFe.fits',
                                            varspecf=None, kins='V15', out_sig=comdisp, index=ind1)
    spec14_1p3i1 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MIUSCAT_Padova00_'+imf+'_1.30_fits/I'+imf+'1.30Zp0.00T14.1254_iPp0.00_baseFe.fits',
                                            varspecf=None, kins='V15', out_sig=comdisp, index=ind1)
    spec14_1p8i1 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MIUSCAT_Padova00_'+imf+'_1.80_fits/I'+imf+'1.80Zp0.00T14.1254_iPp0.00_baseFe.fits',
                                            varspecf=None, kins='V15', out_sig=comdisp, index=ind1)
    spec14_2p0i1 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MIUSCAT_Padova00_'+imf+'_2.00_fits/I'+imf+'2.00Zp0.00T14.1254_iPp0.00_baseFe.fits',
                                            varspecf=None, kins='V15', out_sig=comdisp, index=ind1)

    #Varying age at IMF mu = 1.3
    spec3_1p3i1 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MIUSCAT_Padova00_'+imf+'_1.30_fits/I'+imf+'1.30Zp0.00T03.1623_iPp0.00_baseFe.fits',
                                            varspecf=None, kins='V15', out_sig=comdisp, index=ind1)
    spec5_1p3i1 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MIUSCAT_Padova00_'+imf+'_1.30_fits/I'+imf+'1.30Zp0.00T05.0119_iPp0.00_baseFe.fits',
                                            varspecf=None, kins='V15', out_sig=comdisp, index=ind1)
    spec7_1p3i1 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MIUSCAT_Padova00_'+imf+'_1.30_fits/I'+imf+'1.30Zp0.00T07.0795_iPp0.00_baseFe.fits',
                                            varspecf=None, kins='V15', out_sig=comdisp, index=ind1)
    spec10_1p3i1 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MIUSCAT_Padova00_'+imf+'_1.30_fits/I'+imf+'1.30Zp0.00T10.0000_iPp0.00_baseFe.fits',
                                            varspecf=None, kins='V15', out_sig=comdisp, index=ind1)
    spec11_1p3i1 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MIUSCAT_Padova00_'+imf+'_1.30_fits/I'+imf+'1.30Zp0.00T11.2202_iPp0.00_baseFe.fits',
                                            varspecf=None, kins='V15', out_sig=comdisp, index=ind1)

    #Varying age at IMF mu = 2.0
    spec3_2p0i1 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MIUSCAT_Padova00_'+imf+'_2.00_fits/I'+imf+'2.00Zp0.00T03.1623_iPp0.00_baseFe.fits',
                                            varspecf=None, kins='V15', out_sig=comdisp, index=ind1)
    spec5_2p0i1 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MIUSCAT_Padova00_'+imf+'_2.00_fits/I'+imf+'2.00Zp0.00T05.0119_iPp0.00_baseFe.fits',
                                            varspecf=None, kins='V15', out_sig=comdisp, index=ind1)
    spec7_2p0i1 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MIUSCAT_Padova00_'+imf+'_2.00_fits/I'+imf+'2.00Zp0.00T07.0795_iPp0.00_baseFe.fits',
                                            varspecf=None, kins='V15', out_sig=comdisp, index=ind1)
    spec10_2p0i1 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MIUSCAT_Padova00_'+imf+'_2.00_fits/I'+imf+'2.00Zp0.00T10.0000_iPp0.00_baseFe.fits',
                                            varspecf=None, kins='V15', out_sig=comdisp, index=ind1)
    spec11_2p0i1 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MIUSCAT_Padova00_'+imf+'_2.00_fits/I'+imf+'2.00Zp0.00T11.2202_iPp0.00_baseFe.fits',
                                            varspecf=None, kins='V15', out_sig=comdisp, index=ind1)
    #Varying Z at IMF mu = 1.3, age = 14 Gyr
    spec14_1p3_mp7i1 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MIUSCAT_Padova00_'+imf+'_1.30_fits/I'+imf+'1.30Zm0.71T14.1254_iPp0.00_baseFe.fits',
                                            varspecf=None, kins='V15', out_sig=comdisp, index=ind1)
    spec14_1p3_mp4i1 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MIUSCAT_Padova00_'+imf+'_1.30_fits/I'+imf+'1.30Zm0.40T14.1254_iPp0.00_baseFe.fits',
                                            varspecf=None, kins='V15', out_sig=comdisp, index=ind1)
    spec14_1p3_p2i1 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MIUSCAT_Padova00_'+imf+'_1.30_fits/I'+imf+'1.30Zp0.22T14.1254_iPp0.00_baseFe.fits',
                                            varspecf=None, kins='V15', out_sig=comdisp, index=ind1)
#    spec14_1p8_p4i1 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MIUSCAT_Padova00_'+imf+'_1.30_fits/I'+imf+'1.30Zp0.44T14.1254_iPp0.00_baseFe.fits',
#                                            varspecf=None, kins='V15', out_sig=comdisp, index=ind1)

    #Varying Z at IMF mu = 1.3, age = 7 Gyr
    spec7_1p3_mp7i1 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MIUSCAT_Padova00_'+imf+'_1.30_fits/I'+imf+'1.30Zm0.71T07.0795_iPp0.00_baseFe.fits',
                                            varspecf=None, kins='V15', out_sig=comdisp, index=ind1)
    spec7_1p3_mp4i1 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MIUSCAT_Padova00_'+imf+'_1.30_fits/I'+imf+'1.30Zm0.40T07.0795_iPp0.00_baseFe.fits',
                                            varspecf=None, kins='V15', out_sig=comdisp, index=ind1)
    spec7_1p3_p2i1 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MIUSCAT_Padova00_'+imf+'_1.30_fits/I'+imf+'1.30Zp0.22T07.0795_iPp0.00_baseFe.fits',
                                            varspecf=None, kins='V15', out_sig=comdisp, index=ind1)
#    spec7_1p8_p4i1 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MIUSCAT_Padova00_'+imf+'_1.30_fits/I'+imf+'1.30Zp0.44T07.0795_iPp0.00_baseFe.fits',
#                                            varspecf=None, kins='V15', out_sig=comdisp, index=ind1)
    #Index 2
    #Varying IMF at 14 Gyr
    spec14_0p3i2 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MIUSCAT_Padova00_'+imf+'_0.30_fits/I'+imf+'0.30Zp0.00T14.1254_iPp0.00_baseFe.fits',
                                            varspecf=None, kins='V15', out_sig=comdisp, index=ind2)
    spec14_0p8i2 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MIUSCAT_Padova00_'+imf+'_0.80_fits/I'+imf+'0.80Zp0.00T14.1254_iPp0.00_baseFe.fits',
                                            varspecf=None, kins='V15', out_sig=comdisp, index=ind2)
    spec14_1p3i2 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MIUSCAT_Padova00_'+imf+'_1.30_fits/I'+imf+'1.30Zp0.00T14.1254_iPp0.00_baseFe.fits',
                                            varspecf=None, kins='V15', out_sig=comdisp, index=ind2)
    spec14_1p8i2 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MIUSCAT_Padova00_'+imf+'_1.80_fits/I'+imf+'1.80Zp0.00T14.1254_iPp0.00_baseFe.fits',
                                            varspecf=None, kins='V15', out_sig=comdisp, index=ind2)
    spec14_2p0i2 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MIUSCAT_Padova00_'+imf+'_2.00_fits/I'+imf+'2.00Zp0.00T14.1254_iPp0.00_baseFe.fits',
                                            varspecf=None, kins='V15', out_sig=comdisp, index=ind2)

    #Varying age at IMF mu = 1.3
    spec3_1p3i2 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MIUSCAT_Padova00_'+imf+'_1.30_fits/I'+imf+'1.30Zp0.00T03.1623_iPp0.00_baseFe.fits',
                                            varspecf=None, kins='V15', out_sig=comdisp, index=ind2)
    spec5_1p3i2 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MIUSCAT_Padova00_'+imf+'_1.30_fits/I'+imf+'1.30Zp0.00T05.0119_iPp0.00_baseFe.fits',
                                            varspecf=None, kins='V15', out_sig=comdisp, index=ind2)
    spec7_1p3i2 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MIUSCAT_Padova00_'+imf+'_1.30_fits/I'+imf+'1.30Zp0.00T07.0795_iPp0.00_baseFe.fits',
                                            varspecf=None, kins='V15', out_sig=comdisp, index=ind2)
    spec10_1p3i2 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MIUSCAT_Padova00_'+imf+'_1.30_fits/I'+imf+'1.30Zp0.00T10.0000_iPp0.00_baseFe.fits',
                                            varspecf=None, kins='V15', out_sig=comdisp, index=ind2)
    spec11_1p3i2 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MIUSCAT_Padova00_'+imf+'_1.30_fits/I'+imf+'1.30Zp0.00T11.2202_iPp0.00_baseFe.fits',
                                            varspecf=None, kins='V15', out_sig=comdisp, index=ind2)

    #Varying age at IMF mu = 2.0
    spec3_2p0i2 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MIUSCAT_Padova00_'+imf+'_2.00_fits/I'+imf+'2.00Zp0.00T03.1623_iPp0.00_baseFe.fits',
                                            varspecf=None, kins='V15', out_sig=comdisp, index=ind2)
    spec5_2p0i2 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MIUSCAT_Padova00_'+imf+'_2.00_fits/I'+imf+'2.00Zp0.00T05.0119_iPp0.00_baseFe.fits',
                                            varspecf=None, kins='V15', out_sig=comdisp, index=ind2)
    spec7_2p0i2 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MIUSCAT_Padova00_'+imf+'_2.00_fits/I'+imf+'2.00Zp0.00T07.0795_iPp0.00_baseFe.fits',
                                            varspecf=None, kins='V15', out_sig=comdisp, index=ind2)
    spec10_2p0i2 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MIUSCAT_Padova00_'+imf+'_2.00_fits/I'+imf+'2.00Zp0.00T10.0000_iPp0.00_baseFe.fits',
                                            varspecf=None, kins='V15', out_sig=comdisp, index=ind2)
    spec11_2p0i2 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MIUSCAT_Padova00_'+imf+'_2.00_fits/I'+imf+'2.00Zp0.00T11.2202_iPp0.00_baseFe.fits',
                                            varspecf=None, kins='V15', out_sig=comdisp, index=ind2)
    #Varying Z at IMF mu = 1.3, age = 14 Gyr
    spec14_1p3_mp7i2 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MIUSCAT_Padova00_'+imf+'_1.30_fits/I'+imf+'1.30Zm0.71T14.1254_iPp0.00_baseFe.fits',
                                            varspecf=None, kins='V15', out_sig=comdisp, index=ind2)
    spec14_1p3_mp4i2 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MIUSCAT_Padova00_'+imf+'_1.30_fits/I'+imf+'1.30Zm0.40T14.1254_iPp0.00_baseFe.fits',
                                            varspecf=None, kins='V15', out_sig=comdisp, index=ind2)
    spec14_1p3_p2i2 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MIUSCAT_Padova00_'+imf+'_1.30_fits/I'+imf+'1.30Zp0.22T14.1254_iPp0.00_baseFe.fits',
                                            varspecf=None, kins='V15', out_sig=comdisp, index=ind2)
#    spec14_1p8_p4i2 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MIUSCAT_Padova00_'+imf+'_1.30_fits/I'+imf+'1.30Zp0.44T14.1254_iPp0.00_baseFe.fits',
#                                            varspecf=None, kins='V15', out_sig=comdisp, index=ind2)

    #Varying Z at IMF mu = 1.3, age = 7 Gyr
    spec7_1p3_mp7i2 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MIUSCAT_Padova00_'+imf+'_1.30_fits/I'+imf+'1.30Zm0.71T07.0795_iPp0.00_baseFe.fits',
                                            varspecf=None, kins='V15', out_sig=comdisp, index=ind2)
    spec7_1p3_mp4i2 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MIUSCAT_Padova00_'+imf+'_1.30_fits/I'+imf+'1.30Zm0.40T07.0795_iPp0.00_baseFe.fits',
                                            varspecf=None, kins='V15', out_sig=comdisp, index=ind2)
    spec7_1p3_p2i2 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MIUSCAT_Padova00_'+imf+'_1.30_fits/I'+imf+'1.30Zp0.22T07.0795_iPp0.00_baseFe.fits',
                                            varspecf=None, kins='V15', out_sig=comdisp, index=ind2)
#    spec7_1p8_p4i2 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MIUSCAT_Padova00_'+imf+'_1.30_fits/I'+imf+'1.30Zp0.44T07.0795_iPp0.00_baseFe.fits',
#                                            varspecf=None, kins='V15', out_sig=comdisp, index=ind2)

    if ind1 in ['NaD', 'NaI', 'NaIsdss', 'MgI', 'CaT', 'PaT', 'Mgb']:
        un1 = ' [$\AA$]'
        #Varying age at IMF mu = 1.3, Z = 0
        ind1_3_1p3 = spec3_1p3i1.irindex(spec3_1p3i1.lam[1]-spec3_1p3i1.lam[0], index=ind1, vac_or_air='air')[0]
        ind1_5_1p3 = spec5_1p3i1.irindex(spec5_1p3i1.lam[1]-spec5_1p3i1.lam[0], index=ind1, vac_or_air='air')[0]
        ind1_7_1p3 = spec7_1p3i1.irindex(spec7_1p3i1.lam[1]-spec7_1p3i1.lam[0], index=ind1, vac_or_air='air')[0]
        ind1_10_1p3 = spec10_1p3i1.irindex(spec10_1p3i1.lam[1]-spec10_1p3i1.lam[0], index=ind1, vac_or_air='air')[0]
        ind1_11_1p3 = spec11_1p3i1.irindex(spec11_1p3i1.lam[1]-spec11_1p3i1.lam[0], index=ind1, vac_or_air='air')[0]
        ind1_14_1p3 = spec14_1p3i1.irindex(spec14_1p3i1.lam[1]-spec14_1p3i1.lam[0], index=ind1, vac_or_air='air')[0]
        #Varying age at IMF mu = 2.3, Z = 0
        ind1_3_2p0 = spec3_2p0i1.irindex(spec3_2p0i1.lam[1]-spec3_2p0i1.lam[0], index=ind1, vac_or_air='air')[0]
        ind1_5_2p0 = spec5_2p0i1.irindex(spec5_2p0i1.lam[1]-spec5_2p0i1.lam[0], index=ind1, vac_or_air='air')[0]
        ind1_7_2p0 = spec7_2p0i1.irindex(spec7_2p0i1.lam[1]-spec7_2p0i1.lam[0], index=ind1, vac_or_air='air')[0]
        ind1_10_2p0 = spec10_2p0i1.irindex(spec10_2p0i1.lam[1]-spec10_2p0i1.lam[0], index=ind1, vac_or_air='air')[0]
        ind1_11_2p0 = spec11_2p0i1.irindex(spec11_2p0i1.lam[1]-spec11_2p0i1.lam[0], index=ind1, vac_or_air='air')[0]
        ind1_14_2p0 = spec14_2p0i1.irindex(spec14_2p0i1.lam[1]-spec14_2p0i1.lam[0], index=ind1, vac_or_air='air')[0]
        #Varying IMF at 14 Gyr and Z = 0
        ind1_14_0p3 = spec14_0p3i1.irindex(spec14_0p3i1.lam[1]-spec14_0p3i1.lam[0], index=ind1, vac_or_air='air')[0]
        ind1_14_0p8 = spec14_0p8i1.irindex(spec14_0p8i1.lam[1]-spec14_0p8i1.lam[0], index=ind1, vac_or_air='air')[0]
        ind1_14_1p3 = spec14_1p3i1.irindex(spec14_1p3i1.lam[1]-spec14_1p3i1.lam[0], index=ind1, vac_or_air='air')[0]
        ind1_14_1p8 = spec14_1p8i1.irindex(spec14_1p8i1.lam[1]-spec14_1p8i1.lam[0], index=ind1, vac_or_air='air')[0]
        ind1_14_2p0 = spec14_2p0i1.irindex(spec14_2p0i1.lam[1]-spec14_2p0i1.lam[0], index=ind1, vac_or_air='air')[0]
        #Varying Z at 14 Gyr and IMF mu = 1.3
        ind1_14_1p3_mp4 = spec14_1p3_mp4i1.irindex(spec14_1p3_mp4i1.lam[1]-spec14_1p3_mp4i1.lam[0], index=ind1, vac_or_air='air')[0]
        ind1_14_1p3_p2 = spec14_1p3_p2i1.irindex(spec14_1p3_p2i1.lam[1]-spec14_1p3_p2i1.lam[0], index=ind1, vac_or_air='air')[0]
        #Varying Z at 7 Gyr and IMF mu = 1.3
        ind1_7_1p3_mp4 = spec7_1p3_mp4i1.irindex(spec7_1p3_mp4i1.lam[1]-spec7_1p3_mp4i1.lam[0], index=ind1, vac_or_air='air')[0]
        ind1_7_1p3_p2 = spec7_1p3_p2i1.irindex(spec7_1p3_p2i1.lam[1]-spec7_1p3_p2i1.lam[0], index=ind1, vac_or_air='air')[0]

    elif ind1 == 'TiO':
        un1 = ''
        #Varying age at IMF mu = 1.3, Z = 0
        ind1_3_1p3 = spec3_1p3i1.TiOindex(spec3_1p3i1.lam[1]-spec3_1p3i1.lam[0], vac_or_air='air')[0]
        ind1_5_1p3 = spec5_1p3i1.TiOindex(spec5_1p3i1.lam[1]-spec5_1p3i1.lam[0], vac_or_air='air')[0]
        ind1_7_1p3 = spec7_1p3i1.TiOindex(spec7_1p3i1.lam[1]-spec7_1p3i1.lam[0], vac_or_air='air')[0]
        ind1_10_1p3 = spec10_1p3i1.TiOindex(spec10_1p3i1.lam[1]-spec10_1p3i1.lam[0], vac_or_air='air')[0]
        ind1_11_1p3 = spec11_1p3i1.TiOindex(spec11_1p3i1.lam[1]-spec11_1p3i1.lam[0], vac_or_air='air')[0]
        ind1_14_1p3 = spec14_1p3i1.TiOindex(spec14_1p3i1.lam[1]-spec14_1p3i1.lam[0], vac_or_air='air')[0]
        #Varying age at IMF mu = 2.3, Z = 0
        ind1_3_2p0 = spec3_2p0i1.TiOindex(spec3_2p0i1.lam[1]-spec3_2p0i1.lam[0], vac_or_air='air')[0]
        ind1_5_2p0 = spec5_2p0i1.TiOindex(spec5_2p0i1.lam[1]-spec5_2p0i1.lam[0], vac_or_air='air')[0]
        ind1_7_2p0 = spec7_2p0i1.TiOindex(spec7_2p0i1.lam[1]-spec7_2p0i1.lam[0], vac_or_air='air')[0]
        ind1_10_2p0 = spec10_2p0i1.TiOindex(spec10_2p0i1.lam[1]-spec10_2p0i1.lam[0], vac_or_air='air')[0]
        ind1_11_2p0 = spec11_2p0i1.TiOindex(spec11_2p0i1.lam[1]-spec11_2p0i1.lam[0], vac_or_air='air')[0]
        ind1_14_2p0 = spec14_2p0i1.TiOindex(spec14_2p0i1.lam[1]-spec14_2p0i1.lam[0], vac_or_air='air')[0]
        #Varying IMF at 14 Gyr and Z = 0
        ind1_14_0p3 = spec14_0p3i1.TiOindex(spec14_0p3i1.lam[1]-spec14_0p3i1.lam[0], vac_or_air='air')[0]
        ind1_14_0p8 = spec14_0p8i1.TiOindex(spec14_0p8i1.lam[1]-spec14_0p8i1.lam[0], vac_or_air='air')[0]
        ind1_14_1p3 = spec14_1p3i1.TiOindex(spec14_1p3i1.lam[1]-spec14_1p3i1.lam[0], vac_or_air='air')[0]
        ind1_14_1p8 = spec14_1p8i1.TiOindex(spec14_1p8i1.lam[1]-spec14_1p8i1.lam[0], vac_or_air='air')[0]
        ind1_14_2p0 = spec14_2p0i1.TiOindex(spec14_2p0i1.lam[1]-spec14_2p0i1.lam[0], vac_or_air='air')[0]
        #Varying Z at 14 Gyr and IMF mu = 1.3
        ind1_14_1p3_mp4 = spec14_1p3_mp4i1.TiOindex(spec14_1p3_mp4i1.lam[1]-spec14_1p3_mp4i1.lam[0], vac_or_air='air')[0]
        ind1_14_1p3_p2 = spec14_1p3_p2i1.TiOindex(spec14_1p3_p2i1.lam[1]-spec14_1p3_p2i1.lam[0], vac_or_air='air')[0]
        #Varying Z at 7 Gyr and IMF mu = 1.3
        ind1_7_1p3_mp4 = spec7_1p3_mp4i1.TiOindex(spec7_1p3_mp4i1.lam[1]-spec7_1p3_mp4i1.lam[0], vac_or_air='air')[0]
        ind1_7_1p3_p2 = spec7_1p3_p2i1.TiOindex(spec7_1p3_p2i1.lam[1]-spec7_1p3_p2i1.lam[0], vac_or_air='air')[0]

    elif ind1 == 'Fe':
        un1 = ' [$\AA$]'
        #Varying age at IMF mu = 1.3, Z = 0
        ind1_3_1p3 = spec3_1p3i1.Fe(spec3_1p3i1.lam[1]-spec3_1p3i1.lam[0], vac_or_air='air')[0]
        ind1_5_1p3 = spec5_1p3i1.Fe(spec5_1p3i1.lam[1]-spec5_1p3i1.lam[0], vac_or_air='air')[0]
        ind1_7_1p3 = spec7_1p3i1.Fe(spec7_1p3i1.lam[1]-spec7_1p3i1.lam[0], vac_or_air='air')[0]
        ind1_10_1p3 = spec10_1p3i1.Fe(spec10_1p3i1.lam[1]-spec10_1p3i1.lam[0], vac_or_air='air')[0]
        ind1_11_1p3 = spec11_1p3i1.Fe(spec11_1p3i1.lam[1]-spec11_1p3i1.lam[0], vac_or_air='air')[0]
        ind1_14_1p3 = spec14_1p3i1.Fe(spec14_1p3i1.lam[1]-spec14_1p3i1.lam[0], vac_or_air='air')[0]
        #Varying age at IMF mu = 2.3, Z = 0
        ind1_3_2p0 = spec3_2p0i1.Fe(spec3_2p0i1.lam[1]-spec3_2p0i1.lam[0], vac_or_air='air')[0]
        ind1_5_2p0 = spec5_2p0i1.Fe(spec5_2p0i1.lam[1]-spec5_2p0i1.lam[0], vac_or_air='air')[0]
        ind1_7_2p0 = spec7_2p0i1.Fe(spec7_2p0i1.lam[1]-spec7_2p0i1.lam[0], vac_or_air='air')[0]
        ind1_10_2p0 = spec10_2p0i1.Fe(spec10_2p0i1.lam[1]-spec10_2p0i1.lam[0], vac_or_air='air')[0]
        ind1_11_2p0 = spec11_2p0i1.Fe(spec11_2p0i1.lam[1]-spec11_2p0i1.lam[0], vac_or_air='air')[0]
        ind1_14_2p0 = spec14_2p0i1.Fe(spec14_2p0i1.lam[1]-spec14_2p0i1.lam[0], vac_or_air='air')[0]
        #Varying IMF at 14 Gyr and Z = 0
        ind1_14_0p3 = spec14_0p3i1.Fe(spec14_0p3i1.lam[1]-spec14_0p3i1.lam[0], vac_or_air='air')[0]
        ind1_14_0p8 = spec14_0p8i1.Fe(spec14_0p8i1.lam[1]-spec14_0p8i1.lam[0], vac_or_air='air')[0]
        ind1_14_1p3 = spec14_1p3i1.Fe(spec14_1p3i1.lam[1]-spec14_1p3i1.lam[0], vac_or_air='air')[0]
        ind1_14_1p8 = spec14_1p8i1.Fe(spec14_1p8i1.lam[1]-spec14_1p8i1.lam[0], vac_or_air='air')[0]
        ind1_14_2p0 = spec14_2p0i1.Fe(spec14_2p0i1.lam[1]-spec14_2p0i1.lam[0], vac_or_air='air')[0]
        #Varying Z at 14 Gyr and IMF mu = 1.3
        ind1_14_1p3_mp4 = spec14_1p3_mp4i1.Fe(spec14_1p3_mp4i1.lam[1]-spec14_1p3_mp4i1.lam[0], vac_or_air='air')[0]
        ind1_14_1p3_p2 = spec14_1p3_p2i1.Fe(spec14_1p3_p2i1.lam[1]-spec14_1p3_p2i1.lam[0], vac_or_air='air')[0]
        #Varying Z at 7 Gyr and IMF mu = 1.3
        ind1_7_1p3_mp4 = spec7_1p3_mp4i1.Fe(spec7_1p3_mp4i1.lam[1]-spec7_1p3_mp4i1.lam[0], vac_or_air='air')[0]
        ind1_7_1p3_p2 = spec7_1p3_p2i1.Fe(spec7_1p3_p2i1.lam[1]-spec7_1p3_p2i1.lam[0], vac_or_air='air')[0]

    if ind2 in ['NaD', 'NaI', 'NaIsdss', 'MgI', 'CaT', 'PaT', 'Mgb']:
        un2 = ' [$\AA$]'
        #Varying age at IMF mu = 1.3, Z = 0
        ind2_3_1p3 = spec3_1p3i2.irindex(spec3_1p3i2.lam[1]-spec3_1p3i2.lam[0], index=ind2, vac_or_air='air')[0]
        ind2_5_1p3 = spec5_1p3i2.irindex(spec5_1p3i2.lam[1]-spec5_1p3i2.lam[0], index=ind2, vac_or_air='air')[0]
        ind2_7_1p3 = spec7_1p3i2.irindex(spec7_1p3i2.lam[1]-spec7_1p3i2.lam[0], index=ind2, vac_or_air='air')[0]
        ind2_10_1p3 = spec10_1p3i2.irindex(spec10_1p3i2.lam[1]-spec10_1p3i2.lam[0], index=ind2, vac_or_air='air')[0]
        ind2_11_1p3 = spec11_1p3i2.irindex(spec11_1p3i2.lam[1]-spec11_1p3i2.lam[0], index=ind2, vac_or_air='air')[0]
        ind2_14_1p3 = spec14_1p3i2.irindex(spec14_1p3i2.lam[1]-spec14_1p3i2.lam[0], index=ind2, vac_or_air='air')[0]
        #Varying age at IMF mu = 2.3, Z = 0
        ind2_3_2p0 = spec3_2p0i2.irindex(spec3_2p0i2.lam[1]-spec3_2p0i2.lam[0], index=ind2, vac_or_air='air')[0]
        ind2_5_2p0 = spec5_2p0i2.irindex(spec5_2p0i2.lam[1]-spec5_2p0i2.lam[0], index=ind2, vac_or_air='air')[0]
        ind2_7_2p0 = spec7_2p0i2.irindex(spec7_2p0i2.lam[1]-spec7_2p0i2.lam[0], index=ind2, vac_or_air='air')[0]
        ind2_10_2p0 = spec10_2p0i2.irindex(spec10_2p0i2.lam[1]-spec10_2p0i2.lam[0], index=ind2, vac_or_air='air')[0]
        ind2_11_2p0 = spec11_2p0i2.irindex(spec11_2p0i2.lam[1]-spec11_2p0i2.lam[0], index=ind2, vac_or_air='air')[0]
        ind2_14_2p0 = spec14_2p0i2.irindex(spec14_2p0i2.lam[1]-spec14_2p0i2.lam[0], index=ind2, vac_or_air='air')[0]
        #Varying IMF at 14 Gyr and Z = 0
        ind2_14_0p3 = spec14_0p3i2.irindex(spec14_0p3i2.lam[1]-spec14_0p3i2.lam[0], index=ind2, vac_or_air='air')[0]
        ind2_14_0p8 = spec14_0p8i2.irindex(spec14_0p8i2.lam[1]-spec14_0p8i2.lam[0], index=ind2, vac_or_air='air')[0]
        ind2_14_1p3 = spec14_1p3i2.irindex(spec14_1p3i2.lam[1]-spec14_1p3i2.lam[0], index=ind2, vac_or_air='air')[0]
        ind2_14_1p8 = spec14_1p8i2.irindex(spec14_1p8i2.lam[1]-spec14_1p8i2.lam[0], index=ind2, vac_or_air='air')[0]
        ind2_14_2p0 = spec14_2p0i2.irindex(spec14_2p0i2.lam[1]-spec14_2p0i2.lam[0], index=ind2, vac_or_air='air')[0]
        #Varying Z at 14 Gyr and IMF mu = 1.3
        ind2_14_1p3_mp4 = spec14_1p3_mp4i2.irindex(spec14_1p3_mp4i2.lam[1]-spec14_1p3_mp4i2.lam[0], index=ind2, vac_or_air='air')[0]
        ind2_14_1p3_p2 = spec14_1p3_p2i2.irindex(spec14_1p3_p2i2.lam[1]-spec14_1p3_p2i2.lam[0], index=ind2, vac_or_air='air')[0]
        #Varying Z at 7 Gyr and IMF mu = 1.3
        ind2_7_1p3_mp4 = spec7_1p3_mp4i2.irindex(spec7_1p3_mp4i2.lam[1]-spec7_1p3_mp4i2.lam[0], index=ind2, vac_or_air='air')[0]
        ind2_7_1p3_p2 = spec7_1p3_p2i2.irindex(spec7_1p3_p2i2.lam[1]-spec7_1p3_p2i2.lam[0], index=ind2, vac_or_air='air')[0]

    elif ind2 == 'TiO':
        un2 = ''
        #Varying age at IMF mu = 1.3, Z = 0
        ind2_3_1p3 = spec3_1p3i2.TiOindex(spec3_1p3i2.lam[1]-spec3_1p3i2.lam[0], vac_or_air='air')[0]
        ind2_5_1p3 = spec5_1p3i2.TiOindex(spec5_1p3i2.lam[1]-spec5_1p3i2.lam[0], vac_or_air='air')[0]
        ind2_7_1p3 = spec7_1p3i2.TiOindex(spec7_1p3i2.lam[1]-spec7_1p3i2.lam[0], vac_or_air='air')[0]
        ind2_10_1p3 = spec10_1p3i2.TiOindex(spec10_1p3i2.lam[1]-spec10_1p3i2.lam[0], vac_or_air='air')[0]
        ind2_11_1p3 = spec11_1p3i2.TiOindex(spec11_1p3i2.lam[1]-spec11_1p3i2.lam[0], vac_or_air='air')[0]
        ind2_14_1p3 = spec14_1p3i2.TiOindex(spec14_1p3i2.lam[1]-spec14_1p3i2.lam[0], vac_or_air='air')[0]
        #Varying age at IMF mu = 2.3, Z = 0
        ind2_3_2p0 = spec3_2p0i2.TiOindex(spec3_2p0i2.lam[1]-spec3_2p0i2.lam[0], vac_or_air='air')[0]
        ind2_5_2p0 = spec5_2p0i2.TiOindex(spec5_2p0i2.lam[1]-spec5_2p0i2.lam[0], vac_or_air='air')[0]
        ind2_7_2p0 = spec7_2p0i2.TiOindex(spec7_2p0i2.lam[1]-spec7_2p0i2.lam[0], vac_or_air='air')[0]
        ind2_10_2p0 = spec10_2p0i2.TiOindex(spec10_2p0i2.lam[1]-spec10_2p0i2.lam[0], vac_or_air='air')[0]
        ind2_11_2p0 = spec11_2p0i2.TiOindex(spec11_2p0i2.lam[1]-spec11_2p0i2.lam[0], vac_or_air='air')[0]
        ind2_14_2p0 = spec14_2p0i2.TiOindex(spec14_2p0i2.lam[1]-spec14_2p0i2.lam[0], vac_or_air='air')[0]
        #Varying IMF at 14 Gyr and Z = 0
        ind2_14_0p3 = spec14_0p3i2.TiOindex(spec14_0p3i2.lam[1]-spec14_0p3i2.lam[0], vac_or_air='air')[0]
        ind2_14_0p8 = spec14_0p8i2.TiOindex(spec14_0p8i2.lam[1]-spec14_0p8i2.lam[0], vac_or_air='air')[0]
        ind2_14_1p3 = spec14_1p3i2.TiOindex(spec14_1p3i2.lam[1]-spec14_1p3i2.lam[0], vac_or_air='air')[0]
        ind2_14_1p8 = spec14_1p8i2.TiOindex(spec14_1p8i2.lam[1]-spec14_1p8i2.lam[0], vac_or_air='air')[0]
        ind2_14_2p0 = spec14_2p0i2.TiOindex(spec14_2p0i2.lam[1]-spec14_2p0i2.lam[0], vac_or_air='air')[0]
        #Varying Z at 14 Gyr and IMF mu = 1.3
        ind2_14_1p3_mp4 = spec14_1p3_mp4i2.TiOindex(spec14_1p3_mp4i2.lam[1]-spec14_1p3_mp4i2.lam[0], vac_or_air='air')[0]
        ind2_14_1p3_p2 = spec14_1p3_p2i2.TiOindex(spec14_1p3_p2i2.lam[1]-spec14_1p3_p2i2.lam[0], vac_or_air='air')[0]
        #Varying Z at 7 Gyr and IMF mu = 1.3
        ind2_7_1p3_mp4 = spec7_1p3_mp4i2.TiOindex(spec7_1p3_mp4i2.lam[1]-spec7_1p3_mp4i2.lam[0], vac_or_air='air')[0]
        ind2_7_1p3_p2 = spec7_1p3_p2i2.TiOindex(spec7_1p3_p2i2.lam[1]-spec7_1p3_p2i2.lam[0], vac_or_air='air')[0]

    elif ind2 == 'Fe':
        un2 = ' [$\AA$]'
        #Varying age at IMF mu = 1.3, Z = 0
        ind2_3_1p3 = spec3_1p3i2.Fe(spec3_1p3i2.lam[1]-spec3_1p3i2.lam[0], vac_or_air='air')[0]
        ind2_5_1p3 = spec5_1p3i2.Fe(spec5_1p3i2.lam[1]-spec5_1p3i2.lam[0], vac_or_air='air')[0]
        ind2_7_1p3 = spec7_1p3i2.Fe(spec7_1p3i2.lam[1]-spec7_1p3i2.lam[0], vac_or_air='air')[0]
        ind2_10_1p3 = spec10_1p3i2.Fe(spec10_1p3i2.lam[1]-spec10_1p3i2.lam[0], vac_or_air='air')[0]
        ind2_11_1p3 = spec11_1p3i2.Fe(spec11_1p3i2.lam[1]-spec11_1p3i2.lam[0], vac_or_air='air')[0]
        ind2_14_1p3 = spec14_1p3i2.Fe(spec14_1p3i2.lam[1]-spec14_1p3i2.lam[0], vac_or_air='air')[0]
        #Varying age at IMF mu = 2.3, Z = 0
        ind2_3_2p0 = spec3_2p0i2.Fe(spec3_2p0i2.lam[1]-spec3_2p0i2.lam[0], vac_or_air='air')[0]
        ind2_5_2p0 = spec5_2p0i2.Fe(spec5_2p0i2.lam[1]-spec5_2p0i2.lam[0], vac_or_air='air')[0]
        ind2_7_2p0 = spec7_2p0i2.Fe(spec7_2p0i2.lam[1]-spec7_2p0i2.lam[0], vac_or_air='air')[0]
        ind2_10_2p0 = spec10_2p0i2.Fe(spec10_2p0i2.lam[1]-spec10_2p0i2.lam[0], vac_or_air='air')[0]
        ind2_11_2p0 = spec11_2p0i2.Fe(spec11_2p0i2.lam[1]-spec11_2p0i2.lam[0], vac_or_air='air')[0]
        ind2_14_2p0 = spec14_2p0i2.Fe(spec14_2p0i2.lam[1]-spec14_2p0i2.lam[0], vac_or_air='air')[0]
        #Varying IMF at 14 Gyr and Z = 0
        ind2_14_0p3 = spec14_0p3i2.Fe(spec14_0p3i2.lam[1]-spec14_0p3i2.lam[0], vac_or_air='air')[0]
        ind2_14_0p8 = spec14_0p8i2.Fe(spec14_0p8i2.lam[1]-spec14_0p8i2.lam[0], vac_or_air='air')[0]
        ind2_14_1p3 = spec14_1p3i2.Fe(spec14_1p3i2.lam[1]-spec14_1p3i2.lam[0], vac_or_air='air')[0]
        ind2_14_1p8 = spec14_1p8i2.Fe(spec14_1p8i2.lam[1]-spec14_1p8i2.lam[0], vac_or_air='air')[0]
        ind2_14_2p0 = spec14_2p0i2.Fe(spec14_2p0i2.lam[1]-spec14_2p0i2.lam[0], vac_or_air='air')[0]
        #Varying Z at 14 Gyr and IMF mu = 1.3
        ind2_14_1p3_mp4 = spec14_1p3_mp4i2.Fe(spec14_1p3_mp4i2.lam[1]-spec14_1p3_mp4i2.lam[0], vac_or_air='air')[0]
        ind2_14_1p3_p2 = spec14_1p3_p2i2.Fe(spec14_1p3_p2i2.lam[1]-spec14_1p3_p2i2.lam[0], vac_or_air='air')[0]
        #Varying Z at 7 Gyr and IMF mu = 1.3
        ind2_7_1p3_mp4 = spec7_1p3_mp4i2.Fe(spec7_1p3_mp4i2.lam[1]-spec7_1p3_mp4i2.lam[0], vac_or_air='air')[0]
        ind2_7_1p3_p2 = spec7_1p3_p2i2.Fe(spec7_1p3_p2i2.lam[1]-spec7_1p3_p2i2.lam[0], vac_or_air='air')[0]

    #Arrays
    ind1_2p0 = n.array([ind1_3_2p0, ind1_5_2p0, ind1_7_2p0, ind1_10_2p0, ind1_11_2p0, ind1_14_2p0])
    ind1_1p3 = n.array([ind1_3_1p3, ind1_5_1p3, ind1_7_1p3, ind1_10_1p3, ind1_11_1p3, ind1_14_1p3])
    ind2_2p0 = n.array([ind2_3_2p0, ind2_5_2p0, ind2_7_2p0, ind2_10_2p0, ind2_11_2p0, ind2_14_2p0])
    ind2_1p3 = n.array([ind2_3_1p3, ind2_5_1p3, ind2_7_1p3, ind2_10_1p3, ind2_11_1p3, ind2_14_1p3])

    ind1_imf_Z0p0_14gyr = n.array([ind1_14_0p3, ind1_14_0p8, ind1_14_1p3, ind1_14_1p8, ind1_14_2p0])
    ind1_Z_imf1p3_14gyr = n.array([ind1_14_1p3_mp4, ind1_14_1p3, ind1_14_1p3_p2])#, ind1_14_1p3_p4])
    ind1_Z_imf1p3_7gyr = n.array([ind1_7_1p3_mp4, ind1_7_1p3, ind1_7_1p3_p2])#, ind1_7_1p3_p4])
    ind2_imf_Z0p0_14gyr = n.array([ind2_14_0p3, ind2_14_0p8, ind2_14_1p3, ind2_14_1p8, ind2_14_2p0])
    ind2_Z_imf1p3_14gyr = n.array([ind2_14_1p3_mp4, ind2_14_1p3, ind2_14_1p3_p2])#, ind2_14_1p3_p4])
    ind2_Z_imf1p3_7gyr = n.array([ind2_7_1p3_mp4, ind2_7_1p3, ind2_7_1p3_p2])#, ind2_7_1p3_p4])

    #Plot
    P.rcParams['font.family'] = 'Times New Roman'
    P.rcParams.update({'font.size': gensize})
    P.figure(figsize=(13,10), dpi=72)
    P.subplots_adjust(left=0.14, bottom=0.12)

    p2, = P.plot(ind1_2p0, ind2_2p0, 'r--', lw=1.4)
    p3, = P.plot(ind1_1p3, ind2_1p3, 'r--', lw=1.4)
    p4, = P.plot(ind1_imf_Z0p0_14gyr, ind2_imf_Z0p0_14gyr, 'k-', lw=1.4)
    p5, = P.plot(ind1_Z_imf1p3_14gyr[1:], ind2_Z_imf1p3_14gyr[1:], 'b-', lw=1.4)
    #p6, = P.plot(ind1_Z_imf1p3_7gyr, ind2_Z_imf1p3_7gyr, 'b--', lw=1.4)
    P.scatter(ind1_imf_Z0p0_14gyr, ind2_imf_Z0p0_14gyr, marker='o',c='k',s=[30,50,70,90,110])
    P.scatter(ind1_2p0, ind2_2p0,marker='s',c='r', s=[110,90,70,50,30,10])
    P.scatter(ind1_1p3, ind2_1p3, marker='s',c='r', s=[110,90,70,50,30,10])
    # P.scatter(ind1_Z_imf1p3_14gyr, ind2_Z_imf1p3_14gyr, marker='D',c='b', s=[10,30,50,70,90])
    P.scatter(ind1_Z_imf1p3_14gyr[1:], ind2_Z_imf1p3_14gyr[1:], marker='D',c='b', s=[30,50])
    P.xlabel(ind1+un1, fontsize=tfsize)
    P.ylabel(ind2+un2, fontsize=tfsize)
    P.annotate('$x$ = 3.0', xy=(ind1_imf_Z0p0_14gyr[4], ind2_imf_Z0p0_14gyr[4]),xycoords='data',xytext=(5,5),textcoords='offset points', color='k', fontsize=gensize)
    P.annotate('$x$ = 2.8', xy=(ind1_imf_Z0p0_14gyr[3], ind2_imf_Z0p0_14gyr[3]),xycoords='data',xytext=(5,5),textcoords='offset points', color='k', fontsize=gensize)
    P.annotate('$x$ = 2.3', xy=(ind1_imf_Z0p0_14gyr[2], ind2_imf_Z0p0_14gyr[2]),xycoords='data',xytext=(5,5),textcoords='offset points', color='k', fontsize=gensize)
    P.annotate('$x$ = 1.8', xy=(ind1_imf_Z0p0_14gyr[1], ind2_imf_Z0p0_14gyr[1]),xycoords='data',xytext=(5,5),textcoords='offset points', color='k', fontsize=gensize)
    P.annotate('$x$ = 1.3', xy=(ind1_imf_Z0p0_14gyr[0], ind2_imf_Z0p0_14gyr[0]),xycoords='data',xytext=(5,5),textcoords='offset points', color='k', fontsize=gensize)
    P.annotate('3 Gyr', xy=(ind1_2p0[0],ind2_2p0[0]),xycoords='data',xytext=(10,-10),textcoords='offset points', color='r', fontsize=gensize)
    P.annotate('5', xy=(ind1_2p0[1],ind2_2p0[1]),xycoords='data',xytext=(10,-10),textcoords='offset points', color='r', fontsize=gensize)
    P.annotate('7', xy=(ind1_2p0[2],ind2_2p0[2]),xycoords='data',xytext=(10,-20),textcoords='offset points', color='r', fontsize=gensize)
    P.annotate('10', xy=(ind1_2p0[3],ind2_2p0[3]),xycoords='data',xytext=(5,-15),textcoords='offset points', color='r', fontsize=gensize)
    P.annotate('11', xy=(ind1_2p0[4],ind2_2p0[4]),xycoords='data',xytext=(10,-10),textcoords='offset points', color='r', fontsize=gensize)
    P.annotate('14', xy=(ind1_2p0[5],ind2_2p0[5]),xycoords='data',xytext=(0,-40),textcoords='offset points', color='r', fontsize=gensize)
    # P.annotate('Z = -0.40', xy=(ind1_Z_imf1p3_14gyr[0], ind2_Z_imf1p3_14gyr[0]),xycoords='data',xytext=(10,-15),textcoords='offset points', color='b', fontsize=gensize)
    P.annotate('Z = 0.0', xy=(ind1_Z_imf1p3_14gyr[1], ind2_Z_imf1p3_14gyr[1]),xycoords='data',xytext=(0,-25),textcoords='offset points', color='b', fontsize=gensize)
    P.annotate('Z = 0.22', xy=(ind1_Z_imf1p3_14gyr[2], ind2_Z_imf1p3_14gyr[2]),xycoords='data',xytext=(-30,15),textcoords='offset points', color='b', fontsize=gensize)
    P.xticks(size=lfsize)
    P.yticks(size=lfsize)
    P.minorticks_on()
    P.tick_params(axis='both', direction='in', length=majtsize,
                  width=1.5, which='major')
    P.tick_params(axis='both', direction='in', length=mintsize,
                  width=1.5, which='minor')
    #Plot data from file if given
    if data:
        dat = n.genfromtxt(data, delimiter=',')
        rs = dat[:,0]
        indpos = {'NaI':1, 'NaIsdss':1, 'CaT':4, 'MgI':7, 'TiO':10, 'FeH':13}
        ind1vals = dat[:,indpos[ind1]]
        ind1vars = dat[:,indpos[ind1]+1]
        ind2vals = dat[:,indpos[ind2]]
        ind2vars = dat[:,indpos[ind2]+1]
        markers = [220 -i*20. for i in xrange(len(ind1vals))]

        plot = P.scatter(ind1vals, ind2vals, marker='o', s=[220,200,180,160,140,120,100,80,60,40], c=rs, linewidths=0, cmap=P.cm.copper_r, zorder=10)
        clb = P.colorbar(plot)
        clb.set_label('r [arcsec]', fontsize=tfsize)
        #convert time to a color tuple using the colormap used for scatter
        r_colour = clb.to_rgba(rs)
        #pm0.05 RANDOM ERROR ADDED IN QUADRATURE - 28-05-15 - COMPUTE THIS FOR COMA GALAXIES!!! - 12-11-15
        ind1errfac = 0.00
        ind2errfac = 0.00
        ind1vars = ind1vars+ind1errfac**2
        ind2vars = ind2vars+ind2errfac**2
        a,b,c = P.errorbar(ind1vals, ind2vals, xerr=n.sqrt(ind1vars), yerr=n.sqrt(ind2vars),
                           marker='', ls='', capsize=0, barsabove=False, capthick=1.0,
                           elinewidth=1.8, zorder=10)

        #adjust the color of c[0], which is a LineCollection, to the colormap
        c[0].set_color(r_colour)
        c[1].set_color(r_colour)

    P.show()
#------------#


if __name__=="__main__":
    import sys
    print len(sys.argv)
    if len(sys.argv) not in [5,6]:
        print ''
        print '1. Index 1'
        print '2. Index 2'
        print '3. IMF: [un, bi]'
        print '4. Velocity disp [km/s]'
        print '[5. datafile]'
        print ''
    elif len(sys.argv) == 5:
        indexindexplot(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
    elif len(sys.argv) == 6:
        indexindexplot(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])



#------------#
##if len(sys.argv) != 3:
##    print ' '
##    print '1. Index [NaD, NaI, CaT, MgI, TiO]'
##    print '2. Velocity dispersion [km/s]'
##    print ' '
##    sys.exit()
##
###Common velocity dispersion [km/s]
##comdisp = float(sys.argv[2])
##
##ind = 'NaI'; ax = 'NaI'
##
##sdss = s.loadSDSSFilters()
##vega = s.loadHSTVegaSpec()
##filt1 = sdss['g']
##filt2 = sdss['r']
##
##
###------------#
############
###Plot V12 model Index v colour maps convolved to common dispersion
###Varying IMF at 14 Gyr
##cvd14_0p3 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_un/MIUSCAT_un_0.30_fits_safe/Iun0.30Zp0.00T14.1254.fits',
##                                        varspecf=None, kins='V12', out_sig=comdisp, index=ind)
##cvd14_1p3 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_un/MIUSCAT_un_1.30_fits_safe/Iun1.30Zp0.00T14.1254.fits',
##                                        varspecf=None, kins='V12', out_sig=comdisp, index=ind)
##cvd14_2p3 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_un/MIUSCAT_un_2.30_fits_safe/Iun2.30Zp0.00T14.1254.fits',
##                                        varspecf=None, kins='V12', out_sig=comdisp, index=ind)
##cvd14_3p3 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_un/MIUSCAT_un_3.30_fits_all/Iun3.30Zp0.00T14.1254.fits',
##                                        varspecf=None, kins='V12', out_sig=comdisp, index=ind)
###Varying age at IMF mu = 1.3
##cvd3_1p3 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_un/MIUSCAT_un_1.30_fits_safe/Iun1.30Zp0.00T03.1623.fits',
##                                        varspecf=None, kins='V12', out_sig=comdisp, index=ind)
##cvd5_1p3 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_un/MIUSCAT_un_1.30_fits_safe/Iun1.30Zp0.00T05.0119.fits',
##                                        varspecf=None, kins='V12', out_sig=comdisp, index=ind)
##cvd7_1p3 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_un/MIUSCAT_un_1.30_fits_safe/Iun1.30Zp0.00T07.0795.fits',
##                                        varspecf=None, kins='V12', out_sig=comdisp, index=ind)
##cvd10_1p3 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_un/MIUSCAT_un_1.30_fits_safe/Iun1.30Zp0.00T10.0000.fits',
##                                        varspecf=None, kins='V12', out_sig=comdisp, index=ind)
##cvd11_1p3 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_un/MIUSCAT_un_1.30_fits_safe/Iun1.30Zp0.00T11.2202.fits',
##                                        varspecf=None, kins='V12', out_sig=comdisp, index=ind)
###Varying age at IMF mu = 2.3
##cvd3_2p3 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_un/MIUSCAT_un_2.30_fits_safe/Iun2.30Zp0.00T03.1623.fits',
##                                        varspecf=None, kins='V12', out_sig=comdisp, index=ind)
##cvd5_2p3 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_un/MIUSCAT_un_2.30_fits_safe/Iun2.30Zp0.00T05.0119.fits',
##                                        varspecf=None, kins='V12', out_sig=comdisp, index=ind)
##cvd7_2p3 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_un/MIUSCAT_un_2.30_fits_safe/Iun2.30Zp0.00T07.0795.fits',
##                                        varspecf=None, kins='V12', out_sig=comdisp, index=ind)
##cvd10_2p3 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_un/MIUSCAT_un_2.30_fits_safe/Iun2.30Zp0.00T10.0000.fits',
##                                        varspecf=None, kins='V12', out_sig=comdisp, index=ind)
##cvd11_2p3 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_un/MIUSCAT_un_2.30_fits_safe/Iun2.30Zp0.00T11.2202.fits',
##                                        varspecf=None, kins='V12', out_sig=comdisp, index=ind)
###Varying Z at IMF mu = 1.3, age = 14 Gyr
##cvd14_1p3_mp7 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_un/MIUSCAT_un_1.30_fits_safe/Iun1.30Zm0.71T14.1254.fits',
##                                        varspecf=None, kins='V12', out_sig=comdisp, index=ind)
##cvd14_1p3_mp4 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_un/MIUSCAT_un_1.30_fits_safe/Iun1.30Zm0.40T14.1254.fits',
##                                        varspecf=None, kins='V12', out_sig=comdisp, index=ind)
##cvd14_1p3_p2 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_un/MIUSCAT_un_1.30_fits_safe/Iun1.30Zp0.22T14.1254.fits',
##                                        varspecf=None, kins='V12', out_sig=comdisp, index=ind)
##cvd14_1p3_p4 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_un/MIUSCAT_un_1.30_fits_safe/Iun1.30Zp0.44T14.1254.fits',
##                                        varspecf=None, kins='V12', out_sig=comdisp, index=ind)
##
###Varying Z at IMF mu = 1.3, age = 7 Gyr
##cvd7_1p3_mp7 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_un/MIUSCAT_un_1.30_fits_safe/Iun1.30Zm0.71T07.0795.fits',
##                                        varspecf=None, kins='V12', out_sig=comdisp, index=ind)
##cvd7_1p3_mp4 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_un/MIUSCAT_un_1.30_fits_safe/Iun1.30Zm0.40T07.0795.fits',
##                                        varspecf=None, kins='V12', out_sig=comdisp, index=ind)
##cvd7_1p3_p2 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_un/MIUSCAT_un_1.30_fits_safe/Iun1.30Zp0.22T07.0795.fits',
##                                        varspecf=None, kins='V12', out_sig=comdisp, index=ind)
##cvd7_1p3_p4 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_un/MIUSCAT_un_1.30_fits_safe/Iun1.30Zp0.44T07.0795.fits',
##                                        varspecf=None, kins='V12', out_sig=comdisp, index=ind)
##
###NaI
###Varying age at IMF mu = 1.3, Z = 0
##nai3_1p3 = cvd3_1p3.irindex(cvd3_1p3.lam[1]-cvd3_1p3.lam[0], index='NaI', vac_or_air='air')[0]
##nai5_1p3 = cvd5_1p3.irindex(cvd5_1p3.lam[1]-cvd5_1p3.lam[0], index='NaI', vac_or_air='air')[0]
##nai7_1p3 = cvd7_1p3.irindex(cvd7_1p3.lam[1]-cvd7_1p3.lam[0], index='NaI', vac_or_air='air')[0]
##nai10_1p3 = cvd10_1p3.irindex(cvd10_1p3.lam[1]-cvd10_1p3.lam[0], index='NaI', vac_or_air='air')[0]
##nai11_1p3 = cvd11_1p3.irindex(cvd11_1p3.lam[1]-cvd11_1p3.lam[0], index='NaI', vac_or_air='air')[0]
##nai14_1p3 = cvd14_1p3.irindex(cvd14_1p3.lam[1]-cvd14_1p3.lam[0], index='NaI', vac_or_air='air')[0]
###Varying age at IMF mu = 2.3, Z = 0
##nai3_2p3 = cvd3_2p3.irindex(cvd3_2p3.lam[1]-cvd3_2p3.lam[0], index='NaI', vac_or_air='air')[0]
##nai5_2p3 = cvd5_2p3.irindex(cvd5_2p3.lam[1]-cvd5_2p3.lam[0], index='NaI', vac_or_air='air')[0]
##nai7_2p3 = cvd7_2p3.irindex(cvd7_2p3.lam[1]-cvd7_2p3.lam[0], index='NaI', vac_or_air='air')[0]
##nai10_2p3 = cvd10_2p3.irindex(cvd10_2p3.lam[1]-cvd10_2p3.lam[0], index='NaI', vac_or_air='air')[0]
##nai11_2p3 = cvd11_2p3.irindex(cvd11_2p3.lam[1]-cvd11_2p3.lam[0], index='NaI', vac_or_air='air')[0]
##nai14_2p3 = cvd14_2p3.irindex(cvd14_2p3.lam[1]-cvd14_2p3.lam[0], index='NaI', vac_or_air='air')[0]
###Varying IMF at 14 Gyr and Z = 0
##nai14_0p3 = cvd14_0p3.irindex(cvd14_0p3.lam[1]-cvd14_0p3.lam[0], index='NaI', vac_or_air='air')[0]
##nai14_3p3 = cvd14_3p3.irindex(cvd14_3p3.lam[1]-cvd14_3p3.lam[0], index='NaI', vac_or_air='air')[0]
###Varying Z at 14 Gyr and IMF mu = 1.3
##nai14_1p3_mp7 = cvd14_1p3_mp7.irindex(cvd14_1p3_mp7.lam[1]-cvd14_1p3_mp7.lam[0], index='NaI', vac_or_air='air')[0]
##nai14_1p3_mp4 = cvd14_1p3_mp4.irindex(cvd14_1p3_mp4.lam[1]-cvd14_1p3_mp4.lam[0], index='NaI', vac_or_air='air')[0]
##nai14_1p3_p2 = cvd14_1p3_p2.irindex(cvd14_1p3_p2.lam[1]-cvd14_1p3_p2.lam[0], index='NaI', vac_or_air='air')[0]
##nai14_1p3_p4 = cvd14_1p3_p4.irindex(cvd14_1p3_p4.lam[1]-cvd14_1p3_p4.lam[0], index='NaI', vac_or_air='air')[0]
###Varying Z at 7 Gyr and IMF mu = 1.3
##nai7_1p3_mp7 = cvd7_1p3_mp7.irindex(cvd7_1p3_mp7.lam[1]-cvd7_1p3_mp7.lam[0], index='NaI', vac_or_air='air')[0]
##nai7_1p3_mp4 = cvd7_1p3_mp4.irindex(cvd7_1p3_mp4.lam[1]-cvd7_1p3_mp4.lam[0], index='NaI', vac_or_air='air')[0]
##nai7_1p3_p2 = cvd7_1p3_p2.irindex(cvd7_1p3_p2.lam[1]-cvd7_1p3_p2.lam[0], index='NaI', vac_or_air='air')[0]
##nai7_1p3_p4 = cvd7_1p3_p4.irindex(cvd7_1p3_p4.lam[1]-cvd7_1p3_p4.lam[0], index='NaI', vac_or_air='air')[0]
##
###CaT star
###Varying age at IMF mu = 1.3, Z = 0
##cat3_1p3 = cvd3_1p3.CaT_star(cvd3_1p3.lam[1]-cvd3_1p3.lam[0], vac_or_air='air')[0]
##cat5_1p3 = cvd5_1p3.CaT_star(cvd5_1p3.lam[1]-cvd5_1p3.lam[0], vac_or_air='air')[0]
##cat7_1p3 = cvd7_1p3.CaT_star(cvd7_1p3.lam[1]-cvd7_1p3.lam[0], vac_or_air='air')[0]
##cat10_1p3 = cvd10_1p3.CaT_star(cvd10_1p3.lam[1]-cvd10_1p3.lam[0], vac_or_air='air')[0]
##cat11_1p3 = cvd11_1p3.CaT_star(cvd11_1p3.lam[1]-cvd11_1p3.lam[0], vac_or_air='air')[0]
##cat14_1p3 = cvd14_1p3.CaT_star(cvd14_1p3.lam[1]-cvd14_1p3.lam[0], vac_or_air='air')[0]
###Varying age at IMF mu = 2.3, Z = 0
##cat3_2p3 = cvd3_2p3.CaT_star(cvd3_2p3.lam[1]-cvd3_2p3.lam[0], vac_or_air='air')[0]
##cat5_2p3 = cvd5_2p3.CaT_star(cvd5_2p3.lam[1]-cvd5_2p3.lam[0], vac_or_air='air')[0]
##cat7_2p3 = cvd7_2p3.CaT_star(cvd7_2p3.lam[1]-cvd7_2p3.lam[0], vac_or_air='air')[0]
##cat10_2p3 = cvd10_2p3.CaT_star(cvd10_2p3.lam[1]-cvd10_2p3.lam[0], vac_or_air='air')[0]
##cat11_2p3 = cvd11_2p3.CaT_star(cvd11_2p3.lam[1]-cvd11_2p3.lam[0], vac_or_air='air')[0]
##cat14_2p3 = cvd14_2p3.CaT_star(cvd14_2p3.lam[1]-cvd14_2p3.lam[0], vac_or_air='air')[0]
###Varying IMF at 14 Gyr and Z = 0
##cat14_0p3 = cvd14_0p3.CaT_star(cvd14_0p3.lam[1]-cvd14_0p3.lam[0], vac_or_air='air')[0]
##cat14_3p3 = cvd14_3p3.CaT_star(cvd14_3p3.lam[1]-cvd14_3p3.lam[0], vac_or_air='air')[0]
###Varying Z at 14 Gyr and IMF mu = 1.3
##cat14_1p3_mp7 = cvd14_1p3_mp7.CaT_star(cvd14_1p3_mp7.lam[1]-cvd14_1p3_mp7.lam[0], vac_or_air='air')[0]
##cat14_1p3_mp4 = cvd14_1p3_mp4.CaT_star(cvd14_1p3_mp4.lam[1]-cvd14_1p3_mp4.lam[0], vac_or_air='air')[0]
##cat14_1p3_p2 = cvd14_1p3_p2.CaT_star(cvd14_1p3_p2.lam[1]-cvd14_1p3_p2.lam[0], vac_or_air='air')[0]
##cat14_1p3_p4 = cvd14_1p3_p4.CaT_star(cvd14_1p3_p4.lam[1]-cvd14_1p3_p4.lam[0], vac_or_air='air')[0]
###Varying Z at 7 Gyr and IMF mu = 1.3
##cat7_1p3_mp7 = cvd7_1p3_mp7.CaT_star(cvd7_1p3_mp7.lam[1]-cvd7_1p3_mp7.lam[0], vac_or_air='air')[0]
##cat7_1p3_mp4 = cvd7_1p3_mp4.CaT_star(cvd7_1p3_mp4.lam[1]-cvd7_1p3_mp4.lam[0], vac_or_air='air')[0]
##cat7_1p3_p2 = cvd7_1p3_p2.CaT_star(cvd7_1p3_p2.lam[1]-cvd7_1p3_p2.lam[0], vac_or_air='air')[0]
##cat7_1p3_p4 = cvd7_1p3_p4.CaT_star(cvd7_1p3_p4.lam[1]-cvd7_1p3_p4.lam[0], vac_or_air='air')[0]
##
###CaT
###Varying age at IMF mu = 1.3, Z = 0
##cat3_1p3 = cvd3_1p3.irindex(cvd3_1p3.lam[1]-cvd3_1p3.lam[0], index='CaT', vac_or_air='air')[0]
##cat5_1p3 = cvd5_1p3.irindex(cvd5_1p3.lam[1]-cvd5_1p3.lam[0], index='CaT', vac_or_air='air')[0]
##cat7_1p3 = cvd7_1p3.irindex(cvd7_1p3.lam[1]-cvd7_1p3.lam[0], index='CaT', vac_or_air='air')[0]
##cat10_1p3 = cvd10_1p3.irindex(cvd10_1p3.lam[1]-cvd10_1p3.lam[0], index='CaT', vac_or_air='air')[0]
##cat11_1p3 = cvd11_1p3.irindex(cvd11_1p3.lam[1]-cvd11_1p3.lam[0], index='CaT', vac_or_air='air')[0]
##cat14_1p3 = cvd14_1p3.irindex(cvd14_1p3.lam[1]-cvd14_1p3.lam[0], index='CaT', vac_or_air='air')[0]
###Varying age at IMF mu = 2.3, Z = 0
##cat3_2p3 = cvd3_2p3.irindex(cvd3_2p3.lam[1]-cvd3_2p3.lam[0], index='CaT', vac_or_air='air')[0]
##cat5_2p3 = cvd5_2p3.irindex(cvd5_2p3.lam[1]-cvd5_2p3.lam[0], index='CaT', vac_or_air='air')[0]
##cat7_2p3 = cvd7_2p3.irindex(cvd7_2p3.lam[1]-cvd7_2p3.lam[0], index='CaT', vac_or_air='air')[0]
##cat10_2p3 = cvd10_2p3.irindex(cvd10_2p3.lam[1]-cvd10_2p3.lam[0], index='CaT', vac_or_air='air')[0]
##cat11_2p3 = cvd11_2p3.irindex(cvd11_2p3.lam[1]-cvd11_2p3.lam[0], index='CaT', vac_or_air='air')[0]
##cat14_2p3 = cvd14_2p3.irindex(cvd14_2p3.lam[1]-cvd14_2p3.lam[0], index='CaT', vac_or_air='air')[0]
###Varying IMF at 14 Gyr and Z = 0
##cat14_0p3 = cvd14_0p3.irindex(cvd14_0p3.lam[1]-cvd14_0p3.lam[0], index='CaT', vac_or_air='air')[0]
##cat14_3p3 = cvd14_3p3.irindex(cvd14_3p3.lam[1]-cvd14_3p3.lam[0], index='CaT', vac_or_air='air')[0]
###Varying Z at 14 Gyr and IMF mu = 1.3
##cat14_1p3_mp7 = cvd14_1p3_mp7.irindex(cvd14_1p3_mp7.lam[1]-cvd14_1p3_mp7.lam[0], index='CaT', vac_or_air='air')[0]
##cat14_1p3_mp4 = cvd14_1p3_mp4.irindex(cvd14_1p3_mp4.lam[1]-cvd14_1p3_mp4.lam[0], index='CaT', vac_or_air='air')[0]
##cat14_1p3_p2 = cvd14_1p3_p2.irindex(cvd14_1p3_p2.lam[1]-cvd14_1p3_p2.lam[0], index='CaT', vac_or_air='air')[0]
##cat14_1p3_p4 = cvd14_1p3_p4.irindex(cvd14_1p3_p4.lam[1]-cvd14_1p3_p4.lam[0], index='CaT', vac_or_air='air')[0]
###Varying Z at 7 Gyr and IMF mu = 1.3
##cat7_1p3_mp7 = cvd7_1p3_mp7.irindex(cvd7_1p3_mp7.lam[1]-cvd7_1p3_mp7.lam[0], index='CaT', vac_or_air='air')[0]
##cat7_1p3_mp4 = cvd7_1p3_mp4.irindex(cvd7_1p3_mp4.lam[1]-cvd7_1p3_mp4.lam[0], index='CaT', vac_or_air='air')[0]
##cat7_1p3_p2 = cvd7_1p3_p2.irindex(cvd7_1p3_p2.lam[1]-cvd7_1p3_p2.lam[0], index='CaT', vac_or_air='air')[0]
##cat7_1p3_p4 = cvd7_1p3_p4.irindex(cvd7_1p3_p4.lam[1]-cvd7_1p3_p4.lam[0], index='CaT', vac_or_air='air')[0]
##
###MgI
###Varying age at IMF mu = 1.3, Z = 0
##mgi3_1p3 = cvd3_1p3.irindex(cvd3_1p3.lam[1]-cvd3_1p3.lam[0], index='MgI', vac_or_air='air')[0]
##mgi5_1p3 = cvd5_1p3.irindex(cvd5_1p3.lam[1]-cvd5_1p3.lam[0], index='MgI', vac_or_air='air')[0]
##mgi7_1p3 = cvd7_1p3.irindex(cvd7_1p3.lam[1]-cvd7_1p3.lam[0], index='MgI', vac_or_air='air')[0]
##mgi10_1p3 = cvd10_1p3.irindex(cvd10_1p3.lam[1]-cvd10_1p3.lam[0], index='MgI', vac_or_air='air')[0]
##mgi11_1p3 = cvd11_1p3.irindex(cvd11_1p3.lam[1]-cvd11_1p3.lam[0], index='MgI', vac_or_air='air')[0]
##mgi14_1p3 = cvd14_1p3.irindex(cvd14_1p3.lam[1]-cvd14_1p3.lam[0], index='MgI', vac_or_air='air')[0]
###Varying age at IMF mu = 2.3, Z = 0
##mgi3_2p3 = cvd3_2p3.irindex(cvd3_2p3.lam[1]-cvd3_2p3.lam[0], index='MgI', vac_or_air='air')[0]
##mgi5_2p3 = cvd5_2p3.irindex(cvd5_2p3.lam[1]-cvd5_2p3.lam[0], index='MgI', vac_or_air='air')[0]
##mgi7_2p3 = cvd7_2p3.irindex(cvd7_2p3.lam[1]-cvd7_2p3.lam[0], index='MgI', vac_or_air='air')[0]
##mgi10_2p3 = cvd10_2p3.irindex(cvd10_2p3.lam[1]-cvd10_2p3.lam[0], index='MgI', vac_or_air='air')[0]
##mgi11_2p3 = cvd11_2p3.irindex(cvd11_2p3.lam[1]-cvd11_2p3.lam[0], index='MgI', vac_or_air='air')[0]
##mgi14_2p3 = cvd14_2p3.irindex(cvd14_2p3.lam[1]-cvd14_2p3.lam[0], index='MgI', vac_or_air='air')[0]
###Varying IMF at 14 Gyr and Z = 0
##mgi14_0p3 = cvd14_0p3.irindex(cvd14_0p3.lam[1]-cvd14_0p3.lam[0], index='MgI', vac_or_air='air')[0]
##mgi14_3p3 = cvd14_3p3.irindex(cvd14_3p3.lam[1]-cvd14_3p3.lam[0], index='MgI', vac_or_air='air')[0]
###Varying Z at 14 Gyr and IMF mu = 1.3
##mgi14_1p3_mp7 = cvd14_1p3_mp7.irindex(cvd14_1p3_mp7.lam[1]-cvd14_1p3_mp7.lam[0], index='MgI', vac_or_air='air')[0]
##mgi14_1p3_mp4 = cvd14_1p3_mp4.irindex(cvd14_1p3_mp4.lam[1]-cvd14_1p3_mp4.lam[0], index='MgI', vac_or_air='air')[0]
##mgi14_1p3_p2 = cvd14_1p3_p2.irindex(cvd14_1p3_p2.lam[1]-cvd14_1p3_p2.lam[0], index='MgI', vac_or_air='air')[0]
##mgi14_1p3_p4 = cvd14_1p3_p4.irindex(cvd14_1p3_p4.lam[1]-cvd14_1p3_p4.lam[0], index='MgI', vac_or_air='air')[0]
###Varying Z at 7 Gyr and IMF mu = 1.3
##mgi7_1p3_mp7 = cvd7_1p3_mp7.irindex(cvd7_1p3_mp7.lam[1]-cvd7_1p3_mp7.lam[0], index='MgI', vac_or_air='air')[0]
##mgi7_1p3_mp4 = cvd7_1p3_mp4.irindex(cvd7_1p3_mp4.lam[1]-cvd7_1p3_mp4.lam[0], index='MgI', vac_or_air='air')[0]
##mgi7_1p3_p2 = cvd7_1p3_p2.irindex(cvd7_1p3_p2.lam[1]-cvd7_1p3_p2.lam[0], index='MgI', vac_or_air='air')[0]
##mgi7_1p3_p4 = cvd7_1p3_p4.irindex(cvd7_1p3_p4.lam[1]-cvd7_1p3_p4.lam[0], index='MgI', vac_or_air='air')[0]
##
###TiO
###Varying age at IMF mu = 1.3, Z = 0
##tio3_1p3 = cvd3_1p3.TiOindex(cvd3_1p3.lam[1]-cvd3_1p3.lam[0], vac_or_air='air')[0]
##tio5_1p3 = cvd5_1p3.TiOindex(cvd5_1p3.lam[1]-cvd5_1p3.lam[0], vac_or_air='air')[0]
##tio7_1p3 = cvd7_1p3.TiOindex(cvd7_1p3.lam[1]-cvd7_1p3.lam[0], vac_or_air='air')[0]
##tio10_1p3 = cvd10_1p3.TiOindex(cvd10_1p3.lam[1]-cvd10_1p3.lam[0], vac_or_air='air')[0]
##tio11_1p3 = cvd11_1p3.TiOindex(cvd11_1p3.lam[1]-cvd11_1p3.lam[0], vac_or_air='air')[0]
##tio14_1p3 = cvd14_1p3.TiOindex(cvd14_1p3.lam[1]-cvd14_1p3.lam[0], vac_or_air='air')[0]
###Varying age at IMF mu = 2.3, Z = 0
##tio3_2p3 = cvd3_2p3.TiOindex(cvd3_2p3.lam[1]-cvd3_2p3.lam[0], vac_or_air='air')[0]
##tio5_2p3 = cvd5_2p3.TiOindex(cvd5_2p3.lam[1]-cvd5_2p3.lam[0], vac_or_air='air')[0]
##tio7_2p3 = cvd7_2p3.TiOindex(cvd7_2p3.lam[1]-cvd7_2p3.lam[0], vac_or_air='air')[0]
##tio10_2p3 = cvd10_2p3.TiOindex(cvd10_2p3.lam[1]-cvd10_2p3.lam[0], vac_or_air='air')[0]
##tio11_2p3 = cvd11_2p3.TiOindex(cvd11_2p3.lam[1]-cvd11_2p3.lam[0], vac_or_air='air')[0]
##tio14_2p3 = cvd14_2p3.TiOindex(cvd14_2p3.lam[1]-cvd14_2p3.lam[0], vac_or_air='air')[0]
###Varying IMF at 14 Gyr and Z = 0
##tio14_0p3 = cvd14_0p3.TiOindex(cvd14_0p3.lam[1]-cvd14_0p3.lam[0], vac_or_air='air')[0]
##tio14_3p3 = cvd14_3p3.TiOindex(cvd14_3p3.lam[1]-cvd14_3p3.lam[0], vac_or_air='air')[0]
###Varying Z at 14 Gyr and IMF mu = 1.3
##tio14_1p3_mp7 = cvd14_1p3_mp7.TiOindex(cvd14_1p3_mp7.lam[1]-cvd14_1p3_mp7.lam[0], vac_or_air='air')[0]
##tio14_1p3_mp4 = cvd14_1p3_mp4.TiOindex(cvd14_1p3_mp4.lam[1]-cvd14_1p3_mp4.lam[0], vac_or_air='air')[0]
##tio14_1p3_p2 = cvd14_1p3_p2.TiOindex(cvd14_1p3_p2.lam[1]-cvd14_1p3_p2.lam[0], vac_or_air='air')[0]
##tio14_1p3_p4 = cvd14_1p3_p4.TiOindex(cvd14_1p3_p4.lam[1]-cvd14_1p3_p4.lam[0], vac_or_air='air')[0]
###Varying Z at 7 Gyr and IMF mu = 1.3
##tio7_1p3_mp7 = cvd7_1p3_mp7.TiOindex(cvd7_1p3_mp7.lam[1]-cvd7_1p3_mp7.lam[0], vac_or_air='air')[0]
##tio7_1p3_mp4 = cvd7_1p3_mp4.TiOindex(cvd7_1p3_mp4.lam[1]-cvd7_1p3_mp4.lam[0], vac_or_air='air')[0]
##tio7_1p3_p2 = cvd7_1p3_p2.TiOindex(cvd7_1p3_p2.lam[1]-cvd7_1p3_p2.lam[0], vac_or_air='air')[0]
##tio7_1p3_p4 = cvd7_1p3_p4.TiOindex(cvd7_1p3_p4.lam[1]-cvd7_1p3_p4.lam[0], vac_or_air='air')[0]
##
###NaD
###Varying age at IMF mu = 1.3, Z = 0
##nad3_1p3 = cvd3_1p3.irindex(cvd3_1p3.lam[1]-cvd3_1p3.lam[0], index='NaD', vac_or_air='air')[0]
##nad5_1p3 = cvd5_1p3.irindex(cvd5_1p3.lam[1]-cvd5_1p3.lam[0], index='NaD', vac_or_air='air')[0]
##nad7_1p3 = cvd7_1p3.irindex(cvd7_1p3.lam[1]-cvd7_1p3.lam[0], index='NaD', vac_or_air='air')[0]
##nad10_1p3 = cvd10_1p3.irindex(cvd10_1p3.lam[1]-cvd10_1p3.lam[0], index='NaD', vac_or_air='air')[0]
##nad11_1p3 = cvd11_1p3.irindex(cvd11_1p3.lam[1]-cvd11_1p3.lam[0], index='NaD', vac_or_air='air')[0]
##nad14_1p3 = cvd14_1p3.irindex(cvd14_1p3.lam[1]-cvd14_1p3.lam[0], index='NaD', vac_or_air='air')[0]
###Varying age at IMF mu = 2.3, Z = 0
##nad3_2p3 = cvd3_2p3.irindex(cvd3_2p3.lam[1]-cvd3_2p3.lam[0], index='NaD', vac_or_air='air')[0]
##nad5_2p3 = cvd5_2p3.irindex(cvd5_2p3.lam[1]-cvd5_2p3.lam[0], index='NaD', vac_or_air='air')[0]
##nad7_2p3 = cvd7_2p3.irindex(cvd7_2p3.lam[1]-cvd7_2p3.lam[0], index='NaD', vac_or_air='air')[0]
##nad10_2p3 = cvd10_2p3.irindex(cvd10_2p3.lam[1]-cvd10_2p3.lam[0], index='NaD', vac_or_air='air')[0]
##nad11_2p3 = cvd11_2p3.irindex(cvd11_2p3.lam[1]-cvd11_2p3.lam[0], index='NaD', vac_or_air='air')[0]
##nad14_2p3 = cvd14_2p3.irindex(cvd14_2p3.lam[1]-cvd14_2p3.lam[0], index='NaD', vac_or_air='air')[0]
###Varying IMF at 14 Gyr and Z = 0
##nad14_0p3 = cvd14_0p3.irindex(cvd14_0p3.lam[1]-cvd14_0p3.lam[0], index='NaD', vac_or_air='air')[0]
##nad14_3p3 = cvd14_3p3.irindex(cvd14_3p3.lam[1]-cvd14_3p3.lam[0], index='NaD', vac_or_air='air')[0]
###Varying Z at 14 Gyr and IMF mu = 1.3
##nad14_1p3_mp7 = cvd14_1p3_mp7.irindex(cvd14_1p3_mp7.lam[1]-cvd14_1p3_mp7.lam[0], index='NaD', vac_or_air='air')[0]
##nad14_1p3_mp4 = cvd14_1p3_mp4.irindex(cvd14_1p3_mp4.lam[1]-cvd14_1p3_mp4.lam[0], index='NaD', vac_or_air='air')[0]
##nad14_1p3_p2 = cvd14_1p3_p2.irindex(cvd14_1p3_p2.lam[1]-cvd14_1p3_p2.lam[0], index='NaD', vac_or_air='air')[0]
##nad14_1p3_p4 = cvd14_1p3_p4.irindex(cvd14_1p3_p4.lam[1]-cvd14_1p3_p4.lam[0], index='NaD', vac_or_air='air')[0]
###Varying Z at 7 Gyr and IMF mu = 1.3
##nad7_1p3_mp7 = cvd7_1p3_mp7.irindex(cvd7_1p3_mp7.lam[1]-cvd7_1p3_mp7.lam[0], index='NaD', vac_or_air='air')[0]
##nad7_1p3_mp4 = cvd7_1p3_mp4.irindex(cvd7_1p3_mp4.lam[1]-cvd7_1p3_mp4.lam[0], index='NaD', vac_or_air='air')[0]
##nad7_1p3_p2 = cvd7_1p3_p2.irindex(cvd7_1p3_p2.lam[1]-cvd7_1p3_p2.lam[0], index='NaD', vac_or_air='air')[0]
##nad7_1p3_p4 = cvd7_1p3_p4.irindex(cvd7_1p3_p4.lam[1]-cvd7_1p3_p4.lam[0], index='NaD', vac_or_air='air')[0]
##
###G-R SDSS colour
##col3_1p3 = cvd3_1p3.calcABmag(filt1) - cvd3_1p3.calcABmag(filt2)
##col5_1p3 = cvd5_1p3.calcABmag(filt1) - cvd5_1p3.calcABmag(filt2)
##col7_1p3 = cvd7_1p3.calcABmag(filt1) - cvd7_1p3.calcABmag(filt2)
##col10_1p3 = cvd10_1p3.calcABmag(filt1) - cvd10_1p3.calcABmag(filt2)
##col11_1p3 = cvd11_1p3.calcABmag(filt1) - cvd11_1p3.calcABmag(filt2)
##col14_1p3 = cvd14_1p3.calcABmag(filt1) - cvd14_1p3.calcABmag(filt2)
##col3_2p3 = cvd3_2p3.calcABmag(filt1) - cvd3_2p3.calcABmag(filt2)
##col5_2p3 = cvd5_2p3.calcABmag(filt1) - cvd5_2p3.calcABmag(filt2)
##col7_2p3 = cvd7_2p3.calcABmag(filt1) - cvd7_2p3.calcABmag(filt2)
##col10_2p3 = cvd10_2p3.calcABmag(filt1) - cvd10_2p3.calcABmag(filt2)
##col11_2p3 = cvd11_2p3.calcABmag(filt1) - cvd11_2p3.calcABmag(filt2)
##col14_2p3 = cvd14_2p3.calcABmag(filt1) - cvd14_2p3.calcABmag(filt2)
##col14_0p3 = cvd14_0p3.calcABmag(filt1) - cvd14_0p3.calcABmag(filt2)
##col14_3p3 = cvd14_3p3.calcABmag(filt1) - cvd14_3p3.calcABmag(filt2)
##col14_1p3_mp7 = cvd14_1p3_mp7.calcABmag(filt1) - cvd14_1p3_mp7.calcABmag(filt2)
##col14_1p3_mp4 = cvd14_1p3_mp4.calcABmag(filt1) - cvd14_1p3_mp4.calcABmag(filt2)
##col14_1p3_p2 = cvd14_1p3_p2.calcABmag(filt1) - cvd14_1p3_p2.calcABmag(filt2)
##col14_1p3_p4 = cvd14_1p3_p4.calcABmag(filt1) - cvd14_1p3_p4.calcABmag(filt2)
##col7_1p3_mp7 = cvd7_1p3_mp7.calcABmag(filt1) - cvd7_1p3_mp7.calcABmag(filt2)
##col7_1p3_mp4 = cvd7_1p3_mp4.calcABmag(filt1) - cvd7_1p3_mp4.calcABmag(filt2)
##col7_1p3_p2 = cvd7_1p3_p2.calcABmag(filt1) - cvd7_1p3_p2.calcABmag(filt2)
##col7_1p3_p4 = cvd7_1p3_p4.calcABmag(filt1) - cvd7_1p3_p4.calcABmag(filt2)
##
##
###Arrays
##nai2p3 = n.array([nai3_2p3, nai5_2p3, nai7_2p3, nai10_2p3, nai11_2p3, nai14_2p3])
##nai1p3 = n.array([nai3_1p3, nai5_1p3, nai7_1p3, nai10_1p3, nai11_1p3, nai14_1p3])
##cat2p3 = n.array([cat3_2p3, cat5_2p3, cat7_2p3, cat10_2p3, cat11_2p3, cat14_2p3])
##cat1p3 = n.array([cat3_1p3, cat5_1p3, cat7_1p3, cat10_1p3, cat11_1p3, cat14_1p3])
##mgi2p3 = n.array([mgi3_2p3, mgi5_2p3, mgi7_2p3, mgi10_2p3, mgi11_2p3, mgi14_2p3])
##mgi1p3 = n.array([mgi3_1p3, mgi5_1p3, mgi7_1p3, mgi10_1p3, mgi11_1p3, mgi14_1p3])
##tio2p3 = n.array([tio3_2p3, tio5_2p3, tio7_2p3, tio10_2p3, tio11_2p3, tio14_2p3])
##tio1p3 = n.array([tio3_1p3, tio5_1p3, tio7_1p3, tio10_1p3, tio11_1p3, tio14_1p3])
##nad2p3 = n.array([nad3_2p3, nad5_2p3, nad7_2p3, nad10_2p3, nad11_2p3, nad14_2p3])
##nad1p3 = n.array([nad3_1p3, nad5_1p3, nad7_1p3, nad10_1p3, nad11_1p3, nad14_1p3])
##
##nai_imf_Z0p0_14gyr = n.array([nai14_0p3, nai14_1p3, nai14_2p3])#, nai14_3p3])
##nai_Z_imf1p3_14gyr = n.array([nai14_1p3_mp4, nai14_1p3, nai14_1p3_p2, nai14_1p3_p4])
##nai_Z_imf1p3_7gyr = n.array([nai7_1p3_mp4, nai7_1p3, nai7_1p3_p2, nai7_1p3_p4])
##cat_imf_Z0p0_14gyr = n.array([cat14_0p3, cat14_1p3, cat14_2p3])#, cat14_3p3])
##cat_Z_imf1p3_14gyr = n.array([cat14_1p3_mp4, cat14_1p3, cat14_1p3_p2, cat14_1p3_p4])
##cat_Z_imf1p3_7gyr = n.array([cat7_1p3_mp4, cat7_1p3, cat7_1p3_p2, cat7_1p3_p4])
##mgi_imf_Z0p0_14gyr = n.array([mgi14_0p3, mgi14_1p3, mgi14_2p3])#, mgi14_3p3])
##mgi_Z_imf1p3_14gyr = n.array([mgi14_1p3_mp4, mgi14_1p3, mgi14_1p3_p2, mgi14_1p3_p4])
##mgi_Z_imf1p3_7gyr = n.array([mgi7_1p3_mp4, mgi7_1p3, mgi7_1p3_p2, mgi7_1p3_p4])
##tio_imf_Z0p0_14gyr = n.array([tio14_0p3, tio14_1p3, tio14_2p3])#, tio14_3p3])
##tio_Z_imf1p3_14gyr = n.array([tio14_1p3_mp4, tio14_1p3, tio14_1p3_p2, tio14_1p3_p4])
##tio_Z_imf1p3_7gyr = n.array([tio7_1p3_mp4, tio7_1p3, tio7_1p3_p2, tio7_1p3_p4])
##nad_imf_Z0p0_14gyr = n.array([nad14_0p3, nad14_1p3, nad14_2p3])#, nad14_3p3])
##nad_Z_imf1p3_14gyr = n.array([nad14_1p3_mp4, nad14_1p3, nad14_1p3_p2, nad14_1p3_p4])
##nad_Z_imf1p3_7gyr = n.array([nad7_1p3_mp4, nad7_1p3, nad7_1p3_p2, nad7_1p3_p4])
##
##col2p3 = n.array([col3_2p3, col5_2p3, col7_2p3, col10_2p3, col11_2p3, col14_2p3])
##col1p3 = n.array([col3_1p3, col5_1p3, col7_1p3, col10_1p3, col11_1p3, col14_1p3])
##col_imf_Z0p0_14gyr = n.array([col14_0p3, col14_1p3, col14_2p3])#, col14_3p3])
##col_Z_imf1p3_14gyr = n.array([col14_1p3_mp4, col14_1p3, col14_1p3_p2, col14_1p3_p4])
##col_Z_imf1p3_7gyr = n.array([col7_1p3_mp4, col7_1p3, col7_1p3_p2, col7_1p3_p4])
##
######FONTS
##tfsize = 28
##lfsize = 24
##majtsize = 16
##mintsize = 10
##twidth = 2
##msize = 10
######
##
###Plot
##P.rcParams['font.family'] = 'Times New Roman'
##P.figure(figsize=(13,10), dpi=72)
##P.subplots_adjust(left=0.14, bottom=0.12)
##
#####NaI v NaD
####p2, = P.plot(nai2p3, nad2p3, 'r--', lw=1.4)
####p3, = P.plot(nai1p3, nad1p3, 'r--', lw=1.4)
####p4, = P.plot(nai_imf_Z0p0_14gyr, nad_imf_Z0p0_14gyr, 'k-', lw=1.4)
####p5, = P.plot(nai_Z_imf1p3_14gyr, nad_Z_imf1p3_14gyr, 'b-', lw=1.4)
#####p6, = P.plot(nai_Z_imf1p3_7gyr, nad_Z_imf1p3_7gyr, 'b--', lw=1.4)
####P.scatter(nai_imf_Z0p0_14gyr, nad_imf_Z0p0_14gyr, marker='o',c='k',s=[90,70,50,30])
####P.scatter(nai2p3, nad2p3,marker='s',c='r', s=[110,90,70,50,30,10])
####P.scatter(nai1p3, nad1p3, marker='s',c='r', s=[110,90,70,50,30,10])
####P.scatter(nai_Z_imf1p3_14gyr, nad_Z_imf1p3_14gyr, marker='D',c='b', s=[10,30,50,70,90])
####P.xlabel('NaI0.82 [$\AA$]', fontsize=tfsize)
####P.ylabel('NaD0.59 [$\AA$]', fontsize=tfsize)
##
#######NaD v g-r
##if sys.argv[1] == 'NaD':
##    p2, = P.plot(nad2p3, col2p3, 'r--', lw=1.4)
##    p3, = P.plot(nad1p3, col1p3, 'r--', lw=1.4)
##    p4, = P.plot(nad_imf_Z0p0_14gyr, col_imf_Z0p0_14gyr, 'k-', lw=1.4)
##    p5, = P.plot(nad_Z_imf1p3_14gyr, col_Z_imf1p3_14gyr, 'b-', lw=1.4)
##    P.scatter(nad_imf_Z0p0_14gyr, col_imf_Z0p0_14gyr, marker='o',c='k',s=[90,70,50,30])
##    P.scatter(nad2p3, col2p3,marker='s',c='r', s=[110,90,70,50,30,10])
##    P.scatter(nad1p3, col1p3, marker='s',c='r', s=[110,90,70,50,30,10])
##    P.scatter(nad_Z_imf1p3_14gyr, col_Z_imf1p3_14gyr, marker='D',c='b', s=[10,30,50,70,90])
##    P.xlabel('NaD0.59 [$\AA$]', fontsize=tfsize)
##    P.ylabel('$g - r$', fontsize=tfsize)
##    P.annotate('$x$ = 3.3', xy=(nad_imf_Z0p0_14gyr[2], col_imf_Z0p0_14gyr[2]),xycoords='data',xytext=(-50,10),textcoords='offset points', fontsize=lfsize)
##    P.annotate('$x$ = 2.3', xy=(nad_imf_Z0p0_14gyr[1], col_imf_Z0p0_14gyr[1]),xycoords='data',xytext=(-40,15),textcoords='offset points', fontsize=lfsize)
##    P.annotate('$x$ = 1.3', xy=(nad_imf_Z0p0_14gyr[0], col_imf_Z0p0_14gyr[0]),xycoords='data',xytext=(-40,5),textcoords='offset points', fontsize=lfsize)
##    P.annotate('3 Gyr', xy=(nad2p3[0],col2p3[0]),xycoords='data',xytext=(10,-10),textcoords='offset points', fontsize=lfsize, color='r')
##    P.annotate('5', xy=(nad2p3[1],col2p3[1]),xycoords='data',xytext=(10,-10),textcoords='offset points', fontsize=lfsize, color='r')
##    P.annotate('7', xy=(nad2p3[2],col2p3[2]),xycoords='data',xytext=(10,-5),textcoords='offset points', fontsize=lfsize, color='r')
##    P.annotate('10', xy=(nad2p3[3],col2p3[3]),xycoords='data',xytext=(5,-15),textcoords='offset points', fontsize=lfsize, color='r')
##    P.annotate('11', xy=(nad2p3[4],col2p3[4]),xycoords='data',xytext=(10,-10),textcoords='offset points', fontsize=lfsize, color='r')
##    P.annotate('14', xy=(nad2p3[5],col2p3[5]),xycoords='data',xytext=(5,-20),textcoords='offset points', fontsize=lfsize, color='r')
##    P.annotate('[Z/H] = -0.40', xy=(nad_Z_imf1p3_14gyr[0], col_Z_imf1p3_14gyr[0]),xycoords='data',xytext=(10,-10),textcoords='offset points', fontsize=lfsize, color='b')
##    P.annotate('[Z/H] = 0.0', xy=(nad_Z_imf1p3_14gyr[1], col_Z_imf1p3_14gyr[1]),xycoords='data',xytext=(5,-20),textcoords='offset points', fontsize=lfsize, color='b')
##    P.annotate('[Z/H] = 0.22', xy=(nad_Z_imf1p3_14gyr[2], col_Z_imf1p3_14gyr[2]),xycoords='data',xytext=(-50,20),textcoords='offset points', fontsize=lfsize, color='b')
##    P.annotate('[Z/H] = 0.44', xy=(nad_Z_imf1p3_14gyr[3], col_Z_imf1p3_14gyr[3]),xycoords='data',xytext=(-40,10),textcoords='offset points', fontsize=lfsize, color='b')
##
##
#######NaI v g-r
##elif sys.argv[1] == 'NaI':
##    p2, = P.plot(nai2p3, col2p3, 'r--', lw=1.4)
##    p3, = P.plot(nai1p3, col1p3, 'r--', lw=1.4)
##    p4, = P.plot(nai_imf_Z0p0_14gyr, col_imf_Z0p0_14gyr, 'k-', lw=1.4)
##    p5, = P.plot(nai_Z_imf1p3_14gyr, col_Z_imf1p3_14gyr, 'b-', lw=1.4)
##    P.scatter(nai_imf_Z0p0_14gyr, col_imf_Z0p0_14gyr, marker='o',c='k',s=[90,70,50,30])
##    P.scatter(nai2p3, col2p3,marker='s',c='r', s=[110,90,70,50,30,10])
##    P.scatter(nai1p3, col1p3, marker='s',c='r', s=[110,90,70,50,30,10])
##    P.scatter(nai_Z_imf1p3_14gyr, col_Z_imf1p3_14gyr, marker='D',c='b', s=[10,30,50,70,90])
##    P.xlabel('NaI0.82 [$\AA$]', fontsize=tfsize)
##    P.ylabel('$g - r$', fontsize=tfsize)
##    P.annotate('$x$ = 3.3', xy=(nai_imf_Z0p0_14gyr[2], col_imf_Z0p0_14gyr[2]),xycoords='data',xytext=(-50,10),textcoords='offset points', fontsize=lfsize)
##    P.annotate('$x$ = 2.3', xy=(nai_imf_Z0p0_14gyr[1], col_imf_Z0p0_14gyr[1]),xycoords='data',xytext=(-40,15),textcoords='offset points', fontsize=lfsize)
##    P.annotate('$x$ = 1.3', xy=(nai_imf_Z0p0_14gyr[0], col_imf_Z0p0_14gyr[0]),xycoords='data',xytext=(-40,5),textcoords='offset points', fontsize=lfsize)
##    P.annotate('3 Gyr', xy=(nai2p3[0],col2p3[0]),xycoords='data',xytext=(10,-10),textcoords='offset points', fontsize=lfsize, color='r')
##    P.annotate('5', xy=(nai2p3[1],col2p3[1]),xycoords='data',xytext=(10,-10),textcoords='offset points', fontsize=lfsize, color='r')
##    P.annotate('7', xy=(nai2p3[2],col2p3[2]),xycoords='data',xytext=(10,-5),textcoords='offset points', fontsize=lfsize, color='r')
##    P.annotate('10', xy=(nai2p3[3],col2p3[3]),xycoords='data',xytext=(5,-15),textcoords='offset points', fontsize=lfsize, color='r')
##    P.annotate('11', xy=(nai2p3[4],col2p3[4]),xycoords='data',xytext=(10,-10),textcoords='offset points', fontsize=lfsize, color='r')
##    P.annotate('14', xy=(nai2p3[5],col2p3[5]),xycoords='data',xytext=(5,-20),textcoords='offset points', fontsize=lfsize, color='r')
##    P.annotate('[Z/H] = -0.40', xy=(nai_Z_imf1p3_14gyr[0], col_Z_imf1p3_14gyr[0]),xycoords='data',xytext=(10,-10),textcoords='offset points', fontsize=lfsize, color='b')
##    P.annotate('[Z/H] = 0.0', xy=(nai_Z_imf1p3_14gyr[1], col_Z_imf1p3_14gyr[1]),xycoords='data',xytext=(5,-20),textcoords='offset points', fontsize=lfsize, color='b')
##    P.annotate('[Z/H] = 0.22', xy=(nai_Z_imf1p3_14gyr[2], col_Z_imf1p3_14gyr[2]),xycoords='data',xytext=(-50,20),textcoords='offset points', fontsize=lfsize, color='b')
##    P.annotate('[Z/H] = 0.44', xy=(nai_Z_imf1p3_14gyr[3], col_Z_imf1p3_14gyr[3]),xycoords='data',xytext=(-40,10),textcoords='offset points', fontsize=lfsize, color='b')
##
###CaT v g-r
##elif sys.argv[1] == 'CaT':
##    p2, = P.plot(cat2p3, col2p3, 'r--', lw=1.4)
##    p3, = P.plot(cat1p3, col1p3, 'r--', lw=1.4)
##    p4, = P.plot(cat_imf_Z0p0_14gyr, col_imf_Z0p0_14gyr, 'k-', lw=1.4)
##    p5, = P.plot(cat_Z_imf1p3_14gyr, col_Z_imf1p3_14gyr, 'b-', lw=1.4)
##    P.scatter(cat_imf_Z0p0_14gyr, col_imf_Z0p0_14gyr, marker='o',c='k',s=[90,70,50,30])
##    P.scatter(cat2p3, col2p3,marker='s',c='r', s=[110,90,70,50,30,10])
##    P.scatter(cat1p3, col1p3, marker='s',c='r', s=[110,90,70,50,30,10])
##    P.scatter(cat_Z_imf1p3_14gyr, col_Z_imf1p3_14gyr, marker='D',c='b', s=[10,30,50,70,90])
##    P.xlabel('CaT [$\AA$]', fontsize=tfsize)
##    P.ylabel('$g - r$', fontsize=tfsize)
##    P.annotate('$x$ = 3.3', xy=(cat_imf_Z0p0_14gyr[2], col_imf_Z0p0_14gyr[2]),xycoords='data',xytext=(-50,10),textcoords='offset points', fontsize=lfsize)
##    P.annotate('$x$ = 2.3', xy=(cat_imf_Z0p0_14gyr[1], col_imf_Z0p0_14gyr[1]),xycoords='data',xytext=(-50,20),textcoords='offset points', fontsize=lfsize)
##    P.annotate('$x$ = 1.3', xy=(cat_imf_Z0p0_14gyr[0], col_imf_Z0p0_14gyr[0]),xycoords='data',xytext=(-40,10),textcoords='offset points', fontsize=lfsize)
##    P.annotate('3 Gyr', xy=(cat2p3[0],col2p3[0]),xycoords='data',xytext=(10,-10),textcoords='offset points', fontsize=lfsize, color='r')
##    P.annotate('5', xy=(cat2p3[1],col2p3[1]),xycoords='data',xytext=(10,-10),textcoords='offset points', fontsize=lfsize, color='r')
##    P.annotate('7', xy=(cat2p3[2],col2p3[2]),xycoords='data',xytext=(10,-5),textcoords='offset points', fontsize=lfsize, color='r')
##    P.annotate('10', xy=(cat2p3[3],col2p3[3]),xycoords='data',xytext=(5,-15),textcoords='offset points', fontsize=lfsize, color='r')
##    P.annotate('11', xy=(cat2p3[4],col2p3[4]),xycoords='data',xytext=(10,-10),textcoords='offset points', fontsize=lfsize, color='r')
##    P.annotate('14', xy=(cat2p3[5],col2p3[5]),xycoords='data',xytext=(10,-25),textcoords='offset points', fontsize=lfsize, color='r')
##    P.annotate('[Z/H] = -0.40', xy=(cat_Z_imf1p3_14gyr[0], col_Z_imf1p3_14gyr[0]),xycoords='data',xytext=(30,-5),textcoords='offset points', fontsize=lfsize, color='b')
##    P.annotate('[Z/H] = 0.0', xy=(cat_Z_imf1p3_14gyr[1], col_Z_imf1p3_14gyr[1]),xycoords='data',xytext=(5,-30),textcoords='offset points', fontsize=lfsize, color='b')
##    P.annotate('[Z/H] = 0.22', xy=(cat_Z_imf1p3_14gyr[2], col_Z_imf1p3_14gyr[2]),xycoords='data',xytext=(-60,20),textcoords='offset points', fontsize=lfsize, color='b')
##    P.annotate('[Z/H] = 0.44', xy=(cat_Z_imf1p3_14gyr[3], col_Z_imf1p3_14gyr[3]),xycoords='data',xytext=(-40,10),textcoords='offset points', fontsize=lfsize, color='b')
##
#####CaT star v g-r
##elif sys.argv[1] == 'CaTstar':
##    p2, = P.plot(cat2p3, col2p3, 'r--', lw=1.4)
##    p3, = P.plot(cat1p3, col1p3, 'r--', lw=1.4)
##    p4, = P.plot(cat_imf_Z0p0_14gyr, col_imf_Z0p0_14gyr, 'k-', lw=1.4)
##    p5, = P.plot(cat_Z_imf1p3_14gyr, col_Z_imf1p3_14gyr, 'b-', lw=1.4)
##    P.scatter(cat_imf_Z0p0_14gyr, col_imf_Z0p0_14gyr, marker='o',c='k',s=[90,70,50,30])
##    P.scatter(cat2p3, col2p3,marker='s',c='r', s=[110,90,70,50,30,10])
##    P.scatter(cat1p3, col1p3, marker='s',c='r', s=[110,90,70,50,30,10])
##    P.scatter(cat_Z_imf1p3_14gyr, col_Z_imf1p3_14gyr, marker='D',c='b', s=[10,30,50,70,90])
##    P.xlabel('CaT* [$\AA$]', fontsize=tfsize)
##    P.ylabel('g - r', fontsize=tfsize)
##    P.annotate('$x$ = 3.3', xy=(cat_imf_Z0p0_14gyr[2], col_imf_Z0p0_14gyr[2]),xycoords='data',xytext=(-50,10),textcoords='offset points', fontsize=lfsize)
##    P.annotate('$x$ = 2.3', xy=(cat_imf_Z0p0_14gyr[1], col_imf_Z0p0_14gyr[1]),xycoords='data',xytext=(-50,20),textcoords='offset points', fontsize=lfsize)
##    P.annotate('$x$ = 1.3', xy=(cat_imf_Z0p0_14gyr[0], col_imf_Z0p0_14gyr[0]),xycoords='data',xytext=(-40,5),textcoords='offset points', fontsize=lfsize)
##    P.annotate('3 Gyr', xy=(cat2p3[0],col2p3[0]),xycoords='data',xytext=(10,-10),textcoords='offset points', fontsize=lfsize)
##    P.annotate('5', xy=(cat2p3[1],col2p3[1]),xycoords='data',xytext=(10,-10),textcoords='offset points', fontsize=lfsize)
##    P.annotate('7', xy=(cat2p3[2],col2p3[2]),xycoords='data',xytext=(10,-5),textcoords='offset points', fontsize=lfsize)
##    P.annotate('10', xy=(cat2p3[3],col2p3[3]),xycoords='data',xytext=(5,-15),textcoords='offset points', fontsize=lfsize)
##    P.annotate('11', xy=(cat2p3[4],col2p3[4]),xycoords='data',xytext=(10,-10),textcoords='offset points', fontsize=lfsize)
##    P.annotate('14', xy=(cat2p3[5],col2p3[5]),xycoords='data',xytext=(10,-25),textcoords='offset points', fontsize=lfsize)
##    P.annotate('Z = -0.40', xy=(cat_Z_imf1p3_14gyr[0], col_Z_imf1p3_14gyr[0]),xycoords='data',xytext=(30,-5),textcoords='offset points', fontsize=lfsize)
##    P.annotate('Z = 0.0', xy=(cat_Z_imf1p3_14gyr[1], col_Z_imf1p3_14gyr[1]),xycoords='data',xytext=(5,-30),textcoords='offset points', fontsize=lfsize)
##    P.annotate('Z = 0.22', xy=(cat_Z_imf1p3_14gyr[2], col_Z_imf1p3_14gyr[2]),xycoords='data',xytext=(-60,20),textcoords='offset points', fontsize=lfsize)
##    P.annotate('Z = 0.44', xy=(cat_Z_imf1p3_14gyr[3], col_Z_imf1p3_14gyr[3]),xycoords='data',xytext=(-30,10),textcoords='offset points', fontsize=lfsize)
##
#######MgI v g-r
##elif sys.argv[1] == 'MgI':
##    p2, = P.plot(mgi2p3, col2p3, 'r--', lw=1.4)
##    p3, = P.plot(mgi1p3, col1p3, 'r--', lw=1.4)
##    p4, = P.plot(mgi_imf_Z0p0_14gyr, col_imf_Z0p0_14gyr, 'k-', lw=1.4)
##    p5, = P.plot(mgi_Z_imf1p3_14gyr, col_Z_imf1p3_14gyr, 'b-', lw=1.4)
##    P.scatter(mgi_imf_Z0p0_14gyr, col_imf_Z0p0_14gyr, marker='o',c='k',s=[90,70,50,30])
##    P.scatter(mgi2p3, col2p3,marker='s',c='r', s=[110,90,70,50,30,10])
##    P.scatter(mgi1p3, col1p3, marker='s',c='r', s=[110,90,70,50,30,10])
##    P.scatter(mgi_Z_imf1p3_14gyr, col_Z_imf1p3_14gyr, marker='D',c='b', s=[10,30,50,70,90])
##    P.xlabel('MgI0.88 [$\AA$]', fontsize=tfsize)
##    P.ylabel('g - r', fontsize=tfsize)
##    P.annotate('$x$ = 3.3', xy=(mgi_imf_Z0p0_14gyr[2], col_imf_Z0p0_14gyr[2]),xycoords='data',xytext=(-50,10),textcoords='offset points', fontsize=lfsize)
##    P.annotate('$x$ = 2.3', xy=(mgi_imf_Z0p0_14gyr[1], col_imf_Z0p0_14gyr[1]),xycoords='data',xytext=(-40,10),textcoords='offset points', fontsize=lfsize)
##    P.annotate('$x$ = 1.3', xy=(mgi_imf_Z0p0_14gyr[0], col_imf_Z0p0_14gyr[0]),xycoords='data',xytext=(10,-20),textcoords='offset points', fontsize=lfsize)
##    P.annotate('3 Gyr', xy=(mgi2p3[0],col2p3[0]),xycoords='data',xytext=(10,-10),textcoords='offset points', fontsize=lfsize)
##    P.annotate('5', xy=(mgi2p3[1],col2p3[1]),xycoords='data',xytext=(10,-15),textcoords='offset points', fontsize=lfsize)
##    P.annotate('7', xy=(mgi2p3[2],col2p3[2]),xycoords='data',xytext=(10,-10),textcoords='offset points', fontsize=lfsize)
##    P.annotate('10', xy=(mgi2p3[3],col2p3[3]),xycoords='data',xytext=(5,-15),textcoords='offset points', fontsize=lfsize)
##    P.annotate('11', xy=(mgi2p3[4],col2p3[4]),xycoords='data',xytext=(10,-10),textcoords='offset points', fontsize=lfsize)
##    P.annotate('14', xy=(mgi2p3[5],col2p3[5]),xycoords='data',xytext=(-30,-10),textcoords='offset points', fontsize=lfsize)
##    P.annotate('Z = -0.40', xy=(mgi_Z_imf1p3_14gyr[0], col_Z_imf1p3_14gyr[0]),xycoords='data',xytext=(-80,10),textcoords='offset points', fontsize=lfsize)
##    P.annotate('Z = 0.0', xy=(mgi_Z_imf1p3_14gyr[1], col_Z_imf1p3_14gyr[1]),xycoords='data',xytext=(10,-15),textcoords='offset points', fontsize=lfsize)
##    P.annotate('Z = 0.22', xy=(mgi_Z_imf1p3_14gyr[2], col_Z_imf1p3_14gyr[2]),xycoords='data',xytext=(-50,20),textcoords='offset points', fontsize=lfsize)
##    P.annotate('Z = 0.44', xy=(mgi_Z_imf1p3_14gyr[3], col_Z_imf1p3_14gyr[3]),xycoords='data',xytext=(-30,10),textcoords='offset points', fontsize=lfsize)
##
#########TiO v g-r
##elif sys.argv[1] == 'TiO':
##    p2, = P.plot(tio2p3, col2p3, 'r--', lw=1.4)
##    p3, = P.plot(tio1p3, col1p3, 'r--', lw=1.4)
##    p4, = P.plot(tio_imf_Z0p0_14gyr, col_imf_Z0p0_14gyr, 'k-', lw=1.4)
##    p5, = P.plot(tio_Z_imf1p3_14gyr, col_Z_imf1p3_14gyr, 'b-', lw=1.4)
##    P.scatter(tio_imf_Z0p0_14gyr, col_imf_Z0p0_14gyr, marker='o',c='k',s=[90,70,50,30])
##    P.scatter(tio2p3, col2p3,marker='s',c='r', s=[110,90,70,50,30,10])
##    P.scatter(tio1p3, col1p3, marker='s',c='r', s=[110,90,70,50,30,10])
##    P.scatter(tio_Z_imf1p3_14gyr, col_Z_imf1p3_14gyr, marker='D',c='b', s=[10,30,50,70,90])
##    P.xlabel('TiO0.89', fontsize=tfsize)
##    P.ylabel('g - r', fontsize=tfsize)
##    P.annotate('$x$ = 3.3', xy=(tio_imf_Z0p0_14gyr[2], col_imf_Z0p0_14gyr[2]),xycoords='data',xytext=(5,0),textcoords='offset points', fontsize=lfsize)
##    P.annotate('$x$ = 2.3', xy=(tio_imf_Z0p0_14gyr[1], col_imf_Z0p0_14gyr[1]),xycoords='data',xytext=(-25,10),textcoords='offset points', fontsize=lfsize)
##    P.annotate('$x$ = 1.3', xy=(tio_imf_Z0p0_14gyr[0], col_imf_Z0p0_14gyr[0]),xycoords='data',xytext=(-40,5),textcoords='offset points', fontsize=lfsize)
##    P.annotate('3 Gyr', xy=(tio2p3[0],col2p3[0]),xycoords='data',xytext=(10,-10),textcoords='offset points', fontsize=lfsize)
##    P.annotate('5', xy=(tio2p3[1],col2p3[1]),xycoords='data',xytext=(10,-10),textcoords='offset points', fontsize=lfsize)
##    P.annotate('7', xy=(tio2p3[2],col2p3[2]),xycoords='data',xytext=(10,-20),textcoords='offset points', fontsize=lfsize)
##    P.annotate('10', xy=(tio2p3[3],col2p3[3]),xycoords='data',xytext=(5,-15),textcoords='offset points', fontsize=lfsize)
##    P.annotate('11', xy=(tio2p3[4],col2p3[4]),xycoords='data',xytext=(10,-10),textcoords='offset points', fontsize=lfsize)
##    P.annotate('14', xy=(tio2p3[5],col2p3[5]),xycoords='data',xytext=(0,-40),textcoords='offset points', fontsize=lfsize)
##    P.annotate('Z = -0.40', xy=(tio_Z_imf1p3_14gyr[0], col_Z_imf1p3_14gyr[0]),xycoords='data',xytext=(10,-15),textcoords='offset points', fontsize=lfsize)
##    P.annotate('Z = 0.0', xy=(tio_Z_imf1p3_14gyr[1], col_Z_imf1p3_14gyr[1]),xycoords='data',xytext=(0,-25),textcoords='offset points', fontsize=lfsize)
##    P.annotate('Z = 0.22', xy=(tio_Z_imf1p3_14gyr[2], col_Z_imf1p3_14gyr[2]),xycoords='data',xytext=(-30,15),textcoords='offset points', fontsize=lfsize)
##    P.annotate('Z = 0.44', xy=(tio_Z_imf1p3_14gyr[3], col_Z_imf1p3_14gyr[3]),xycoords='data',xytext=(-30,10),textcoords='offset points', fontsize=lfsize)
##
##P.ylim([0.65,0.95])
##P.xticks(size=lfsize)
##P.yticks(size=lfsize)
##P.tick_params(axis='both', direction='in', length=majtsize,
##              width=1.5, which='major')
##P.tick_params(axis='both', direction='in', length=mintsize,
##              width=1.5, which='minor')
##P.show()
###------------#
