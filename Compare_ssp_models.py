'''Code to compare V15 and CvD12 SPS models

Author: Simon Zieleniewski

Last updated: 02-08-16

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

ind1='Mgb'
ind2 = 'Fe'
comdisp = 200.0
imf='un'
#V15 models
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
spec14_1p3i1 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MIUSCAT_Padova00_'+imf+'_1.30_fits/I'+imf+'1.30Zp0.00T14.1254_iPp0.00_baseFe.fits',
                                        varspecf=None, kins='V15', out_sig=comdisp, index=ind1)

ind1_3_1p3 = spec3_1p3i1.irindex(spec3_1p3i1.lam[1]-spec3_1p3i1.lam[0], index=ind1, vac_or_air='air')[0]
ind1_5_1p3 = spec5_1p3i1.irindex(spec5_1p3i1.lam[1]-spec5_1p3i1.lam[0], index=ind1, vac_or_air='air')[0]
ind1_7_1p3 = spec7_1p3i1.irindex(spec7_1p3i1.lam[1]-spec7_1p3i1.lam[0], index=ind1, vac_or_air='air')[0]
ind1_10_1p3 = spec10_1p3i1.irindex(spec10_1p3i1.lam[1]-spec10_1p3i1.lam[0], index=ind1, vac_or_air='air')[0]
ind1_11_1p3 = spec11_1p3i1.irindex(spec11_1p3i1.lam[1]-spec11_1p3i1.lam[0], index=ind1, vac_or_air='air')[0]
ind1_14_1p3 = spec14_1p3i1.irindex(spec14_1p3i1.lam[1]-spec14_1p3i1.lam[0], index=ind1, vac_or_air='air')[0]

ind2_3_1p3 = spec3_1p3i1.Fe(spec3_1p3i1.lam[1]-spec3_1p3i1.lam[0], vac_or_air='air')[0]
ind2_5_1p3 = spec5_1p3i1.Fe(spec5_1p3i1.lam[1]-spec5_1p3i1.lam[0], vac_or_air='air')[0]
ind2_7_1p3 = spec7_1p3i1.Fe(spec7_1p3i1.lam[1]-spec7_1p3i1.lam[0], vac_or_air='air')[0]
ind2_10_1p3 = spec10_1p3i1.Fe(spec10_1p3i1.lam[1]-spec10_1p3i1.lam[0], vac_or_air='air')[0]
ind2_11_1p3 = spec11_1p3i1.Fe(spec11_1p3i1.lam[1]-spec11_1p3i1.lam[0], vac_or_air='air')[0]
ind2_14_1p3 = spec14_1p3i1.Fe(spec14_1p3i1.lam[1]-spec14_1p3i1.lam[0], vac_or_air='air')[0]

v15ind1 = [ind1_3_1p3, ind1_5_1p3, ind1_7_1p3, ind1_10_1p3, ind1_11_1p3, ind1_14_1p3]
v15ind2 = [ind2_3_1p3, ind2_5_1p3, ind2_7_1p3, ind2_10_1p3, ind2_11_1p3, ind2_14_1p3]
v15ages = [3.0, 5.0, 7.0, 10.0, 11.0, 14.0]

#CvD12 models
cvd3 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/CvD12/t03.0_solar.ssp', varspecf=None, kins='CvD12', out_sig=comdisp, index=ind1)
cvd5 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/CvD12/t05.0_solar.ssp', varspecf=None, kins='CvD12', out_sig=comdisp, index=ind1)
cvd7 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/CvD12/t07.0_solar.ssp', varspecf=None, kins='CvD12', out_sig=comdisp, index=ind1)
cvd9 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/CvD12/t09.0_solar.ssp', varspecf=None, kins='CvD12', out_sig=comdisp, index=ind1)
cvd11 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/CvD12/t11.0_solar.ssp', varspecf=None, kins='CvD12', out_sig=comdisp, index=ind1)
cvd13 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/CvD12/t13.5_solar.ssp', varspecf=None, kins='CvD12', out_sig=comdisp, index=ind1)

cvd3ind1 = cvd3.irindex(cvd3.lam[1]-cvd3.lam[0], index=ind1)[0]
cvd5ind1 = cvd5.irindex(cvd5.lam[1]-cvd5.lam[0], index=ind1)[0]
cvd7ind1 = cvd7.irindex(cvd7.lam[1]-cvd7.lam[0], index=ind1)[0]
cvd9ind1 = cvd9.irindex(cvd9.lam[1]-cvd9.lam[0], index=ind1)[0]
cvd11ind1 = cvd11.irindex(cvd11.lam[1]-cvd11.lam[0], index=ind1)[0]
cvd13ind1 = cvd13.irindex(cvd13.lam[1]-cvd13.lam[0], index=ind1)[0]

cvd3ind2 = cvd3.Fe(cvd3.lam[1]-cvd3.lam[0])[0]
cvd5ind2 = cvd5.Fe(cvd5.lam[1]-cvd5.lam[0])[0]
cvd7ind2 = cvd7.Fe(cvd7.lam[1]-cvd7.lam[0])[0]
cvd9ind2 = cvd9.Fe(cvd9.lam[1]-cvd9.lam[0])[0]
cvd11ind2 = cvd11.Fe(cvd11.lam[1]-cvd11.lam[0])[0]
cvd13ind2 = cvd13.Fe(cvd13.lam[1]-cvd13.lam[0])[0]

cvdind1 = [cvd3ind1[2], cvd5ind1[2], cvd7ind1[2], cvd9ind1[2], cvd11ind1[2], cvd13ind1[2]]
cvdind2 = [cvd3ind2[2], cvd5ind2[2], cvd7ind2[2], cvd9ind2[2], cvd11ind2[2], cvd13ind2[2]]
cvdages = [3.0, 5.0, 7.0, 9.0, 11.0, 13.5]

#Plot Mgb v age
P.rcParams['font.family'] = 'Times New Roman'
P.rcParams.update({'font.size': gensize})
P.figure(figsize=(13,10), dpi=72)
P.subplots_adjust(left=0.14, bottom=0.12)

P.plot(cvdages, cvdind1, 'ko-', lw=2.0, ms=msize)
P.plot(v15ages, v15ind1, 'rs-', lw=2.0, ms=msize)
P.xlabel('Age [Gyr]', fontsize=tfsize)
P.ylabel('Mgb [$\AA$]', fontsize=tfsize)
P.xticks(size=lfsize)
P.yticks(size=lfsize)
P.minorticks_on()
P.tick_params(axis='both', direction='in', length=majtsize,
                width=1.5, which='major')
P.tick_params(axis='both', direction='in', length=mintsize,
                width=1.5, which='minor')
P.legend(['CvD12 models', 'V15 models'],
loc='lower right', frameon=False, borderaxespad=0.3, labelspacing=0.3, handlelength=2.2)
P.tight_layout()
P.show()

#Plot <Fe> v age
P.rcParams['font.family'] = 'Times New Roman'
P.rcParams.update({'font.size': gensize})
P.figure(figsize=(13,10), dpi=72)
P.subplots_adjust(left=0.14, bottom=0.12)

P.plot(cvdages, cvdind2, 'ko-', lw=2.0, ms=msize)
P.plot(v15ages, v15ind2, 'rs-', lw=2.0, ms=msize)
P.xlabel('Age [Gyr]', fontsize=tfsize)
P.ylabel('<Fe> [$\AA$]', fontsize=tfsize)
P.xticks(size=lfsize)
P.yticks(size=lfsize)
P.minorticks_on()
P.tick_params(axis='both', direction='in', length=majtsize,
                width=1.5, which='major')
P.tick_params(axis='both', direction='in', length=mintsize,
                width=1.5, which='minor')
P.legend(['CvD12 models', 'V15 models'],
loc='lower right', frameon=False, borderaxespad=0.3, labelspacing=0.3, handlelength=2.2)
P.tight_layout()
P.show()