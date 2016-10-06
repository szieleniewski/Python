import numpy as n
import astropy.io.fits as p
import sys

#iTp0.00_Ep0.00 or iTp0.40_Ep0.40
#iTp0.00_Ep0.00
afe = 'iTp0.00_Ep0.00'
imf = 'bi'
#Z=0.06 specs
v3,v3head = p.getdata('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MILES_BaSTI_'+imf+'_1.30_fits/M'+imf+'1.30Zp0.06T03.0000_'+afe+'.fits', header=True)
v5,v5head = p.getdata('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MILES_BaSTI_'+imf+'_1.30_fits/M'+imf+'1.30Zp0.06T05.0000_'+afe+'.fits', header=True)
v7,v7head = p.getdata('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MILES_BaSTI_'+imf+'_1.30_fits/M'+imf+'1.30Zp0.06T07.0000_'+afe+'.fits', header=True)
v9,v9head = p.getdata('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MILES_BaSTI_'+imf+'_1.30_fits/M'+imf+'1.30Zp0.06T09.0000_'+afe+'.fits', header=True)
v11,v11head = p.getdata('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MILES_BaSTI_'+imf+'_1.30_fits/M'+imf+'1.30Zp0.06T11.0000_'+afe+'.fits', header=True)
v13,v13head = p.getdata('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MILES_BaSTI_'+imf+'_1.30_fits/M'+imf+'1.30Zp0.06T13.5000_'+afe+'.fits', header=True)

#Z=0.26 specs
v3_0p26,v3head_0p26 = p.getdata('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MILES_BaSTI_'+imf+'_1.30_fits/M'+imf+'1.30Zp0.26T03.0000_'+afe+'.fits', header=True)
v5_0p26,v5head_0p26 = p.getdata('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MILES_BaSTI_'+imf+'_1.30_fits/M'+imf+'1.30Zp0.26T05.0000_'+afe+'.fits', header=True)
v7_0p26,v7head_0p26 = p.getdata('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MILES_BaSTI_'+imf+'_1.30_fits/M'+imf+'1.30Zp0.26T07.0000_'+afe+'.fits', header=True)
v9_0p26,v9head_0p26 = p.getdata('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MILES_BaSTI_'+imf+'_1.30_fits/M'+imf+'1.30Zp0.26T09.0000_'+afe+'.fits', header=True)
v11_0p26,v11head_0p26 = p.getdata('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MILES_BaSTI_'+imf+'_1.30_fits/M'+imf+'1.30Zp0.26T11.0000_'+afe+'.fits', header=True)
v13_0p26,v13head_0p26 = p.getdata('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MILES_BaSTI_'+imf+'_1.30_fits/M'+imf+'1.30Zp0.26T13.5000_'+afe+'.fits', header=True)

#Z=0.40 specs
v3_0p4,v3head_0p4 = p.getdata('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MILES_BaSTI_'+imf+'_1.30_fits/M'+imf+'1.30Zp0.40T03.0000_'+afe+'.fits', header=True)
v5_0p4,v5head_0p4 = p.getdata('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MILES_BaSTI_'+imf+'_1.30_fits/M'+imf+'1.30Zp0.40T05.0000_'+afe+'.fits', header=True)
v7_0p4,v7head_0p4 = p.getdata('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MILES_BaSTI_'+imf+'_1.30_fits/M'+imf+'1.30Zp0.40T07.0000_'+afe+'.fits', header=True)
v9_0p4,v9head_0p4 = p.getdata('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MILES_BaSTI_'+imf+'_1.30_fits/M'+imf+'1.30Zp0.40T09.0000_'+afe+'.fits', header=True)
v11_0p4,v11head_0p4 = p.getdata('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MILES_BaSTI_'+imf+'_1.30_fits/M'+imf+'1.30Zp0.40T11.0000_'+afe+'.fits', header=True)
v13_0p4,v13head_0p4 = p.getdata('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MILES_BaSTI_'+imf+'_1.30_fits/M'+imf+'1.30Zp0.40T13.5000_'+afe+'.fits', header=True)

#Z=-0.25 specs
v3_m0p25,v3head_m0p25 = p.getdata('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MILES_BaSTI_'+imf+'_1.30_fits/M'+imf+'1.30Zm0.25T03.0000_'+afe+'.fits', header=True)
v5_m0p25,v5head_m0p25 = p.getdata('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MILES_BaSTI_'+imf+'_1.30_fits/M'+imf+'1.30Zm0.25T05.0000_'+afe+'.fits', header=True)
v7_m0p25,v7head_m0p25 = p.getdata('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MILES_BaSTI_'+imf+'_1.30_fits/M'+imf+'1.30Zm0.25T07.0000_'+afe+'.fits', header=True)
v9_m0p25,v9head_m0p25 = p.getdata('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MILES_BaSTI_'+imf+'_1.30_fits/M'+imf+'1.30Zm0.25T09.0000_'+afe+'.fits', header=True)
v11_m0p25,v11head_m0p25 = p.getdata('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MILES_BaSTI_'+imf+'_1.30_fits/M'+imf+'1.30Zm0.25T11.0000_'+afe+'.fits', header=True)
v13_m0p25,v13head_m0p25 = p.getdata('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MILES_BaSTI_'+imf+'_1.30_fits/M'+imf+'1.30Zm0.25T13.5000_'+afe+'.fits', header=True)

#iTp0.40_Ep0.40
afe = 'iTp0.40_Ep0.40'
#Z=0.06 specs
v3_afe4,v3head_afe4 = p.getdata('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MILES_BaSTI_'+imf+'_1.30_fits/M'+imf+'1.30Zp0.06T03.0000_'+afe+'.fits', header=True)
v5_afe4,v5head_afe4 = p.getdata('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MILES_BaSTI_'+imf+'_1.30_fits/M'+imf+'1.30Zp0.06T05.0000_'+afe+'.fits', header=True)
v7_afe4,v7head_afe4 = p.getdata('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MILES_BaSTI_'+imf+'_1.30_fits/M'+imf+'1.30Zp0.06T07.0000_'+afe+'.fits', header=True)
v9_afe4,v9head_afe4 = p.getdata('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MILES_BaSTI_'+imf+'_1.30_fits/M'+imf+'1.30Zp0.06T09.0000_'+afe+'.fits', header=True)
v11_afe4,v11head_afe4 = p.getdata('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MILES_BaSTI_'+imf+'_1.30_fits/M'+imf+'1.30Zp0.06T11.0000_'+afe+'.fits', header=True)
v13_afe4,v13head_afe4 = p.getdata('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MILES_BaSTI_'+imf+'_1.30_fits/M'+imf+'1.30Zp0.06T13.5000_'+afe+'.fits', header=True)

#Z=0.26 specs
v3_0p26_afe4,v3head_0p26_afe4 = p.getdata('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MILES_BaSTI_'+imf+'_1.30_fits/M'+imf+'1.30Zp0.26T03.0000_'+afe+'.fits', header=True)
v5_0p26_afe4,v5head_0p26_afe4 = p.getdata('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MILES_BaSTI_'+imf+'_1.30_fits/M'+imf+'1.30Zp0.26T05.0000_'+afe+'.fits', header=True)
v7_0p26_afe4,v7head_0p26_afe4 = p.getdata('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MILES_BaSTI_'+imf+'_1.30_fits/M'+imf+'1.30Zp0.26T07.0000_'+afe+'.fits', header=True)
v9_0p26_afe4,v9head_0p26_afe4 = p.getdata('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MILES_BaSTI_'+imf+'_1.30_fits/M'+imf+'1.30Zp0.26T09.0000_'+afe+'.fits', header=True)
v11_0p26_afe4,v11head_0p26_afe4 = p.getdata('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MILES_BaSTI_'+imf+'_1.30_fits/M'+imf+'1.30Zp0.26T11.0000_'+afe+'.fits', header=True)
v13_0p26_afe4,v13head_0p26_afe4 = p.getdata('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MILES_BaSTI_'+imf+'_1.30_fits/M'+imf+'1.30Zp0.26T13.5000_'+afe+'.fits', header=True)

#Z=0.40 specs
v3_0p4_afe4,v3head_0p4_afe4 = p.getdata('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MILES_BaSTI_'+imf+'_1.30_fits/M'+imf+'1.30Zp0.40T03.0000_'+afe+'.fits', header=True)
v5_0p4_afe4,v5head_0p4_afe4 = p.getdata('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MILES_BaSTI_'+imf+'_1.30_fits/M'+imf+'1.30Zp0.40T05.0000_'+afe+'.fits', header=True)
v7_0p4_afe4,v7head_0p4_afe4 = p.getdata('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MILES_BaSTI_'+imf+'_1.30_fits/M'+imf+'1.30Zp0.40T07.0000_'+afe+'.fits', header=True)
v9_0p4_afe4,v9head_0p4_afe4 = p.getdata('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MILES_BaSTI_'+imf+'_1.30_fits/M'+imf+'1.30Zp0.40T09.0000_'+afe+'.fits', header=True)
v11_0p4_afe4,v11head_0p4_afe4 = p.getdata('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MILES_BaSTI_'+imf+'_1.30_fits/M'+imf+'1.30Zp0.40T11.0000_'+afe+'.fits', header=True)
v13_0p4_afe4,v13head_0p4_afe4 = p.getdata('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MILES_BaSTI_'+imf+'_1.30_fits/M'+imf+'1.30Zp0.40T13.5000_'+afe+'.fits', header=True)

#Z=-0.25 specs
v3_m0p25_afe4,v3head_m0p25_afe4 = p.getdata('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MILES_BaSTI_'+imf+'_1.30_fits/M'+imf+'1.30Zm0.25T03.0000_'+afe+'.fits', header=True)
v5_m0p25_afe4,v5head_m0p25_afe4 = p.getdata('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MILES_BaSTI_'+imf+'_1.30_fits/M'+imf+'1.30Zm0.25T05.0000_'+afe+'.fits', header=True)
v7_m0p25_afe4,v7head_m0p25_afe4 = p.getdata('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MILES_BaSTI_'+imf+'_1.30_fits/M'+imf+'1.30Zm0.25T07.0000_'+afe+'.fits', header=True)
v9_m0p25_afe4,v9head_m0p25_afe4 = p.getdata('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MILES_BaSTI_'+imf+'_1.30_fits/M'+imf+'1.30Zm0.25T09.0000_'+afe+'.fits', header=True)
v11_m0p25_afe4,v11head_m0p25_afe4 = p.getdata('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MILES_BaSTI_'+imf+'_1.30_fits/M'+imf+'1.30Zm0.25T11.0000_'+afe+'.fits', header=True)
v13_m0p25_afe4,v13head_m0p25_afe4 = p.getdata('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MILES_BaSTI_'+imf+'_1.30_fits/M'+imf+'1.30Zm0.25T13.5000_'+afe+'.fits', header=True)

#Create new SSPs with *arbitrary* amount of [a/Fe]
afeval = 0.20
v3_nafe = v3 * ((v3_afe4/v3 - 1)*(afeval/0.4) + 1)
v5_nafe = v5 * ((v5_afe4/v5 - 1)*(afeval/0.4) + 1)
v7_nafe = v7 * ((v7_afe4/v7 - 1)*(afeval/0.4) + 1)
v9_nafe = v9 * ((v9_afe4/v9 - 1)*(afeval/0.4) + 1)
v11_nafe = v11 * ((v11_afe4/v11 - 1)*(afeval/0.4) + 1)
v13_nafe = v13 * ((v13_afe4/v13 - 1)*(afeval/0.4) + 1)

v3_0p26_nafe = v3_0p26 * ((v3_0p26_afe4/v3_0p26 - 1)*(afeval/0.4) + 1)
v5_0p26_nafe = v5_0p26 * ((v5_0p26_afe4/v5_0p26 - 1)*(afeval/0.4) + 1)
v7_0p26_nafe = v7_0p26 * ((v7_0p26_afe4/v7_0p26 - 1)*(afeval/0.4) + 1)
v9_0p26_nafe = v9_0p26 * ((v9_0p26_afe4/v9_0p26 - 1)*(afeval/0.4) + 1)
v11_0p26_nafe = v11_0p26 * ((v11_0p26_afe4/v11_0p26 - 1)*(afeval/0.4) + 1)
v13_0p26_nafe = v13_0p26 * ((v13_0p26_afe4/v13_0p26 - 1)*(afeval/0.4) + 1)

v3_0p4_nafe = v3_0p4 * ((v3_0p4_afe4/v3_0p4 - 1)*(afeval/0.4) + 1)
v5_0p4_nafe = v5_0p4 * ((v5_0p4_afe4/v5_0p4 - 1)*(afeval/0.4) + 1)
v7_0p4_nafe = v7_0p4 * ((v7_0p4_afe4/v7_0p4 - 1)*(afeval/0.4) + 1)
v9_0p4_nafe = v9_0p4 * ((v9_0p4_afe4/v9_0p4 - 1)*(afeval/0.4) + 1)
v11_0p4_nafe = v11_0p4 * ((v11_0p4_afe4/v11_0p4 - 1)*(afeval/0.4) + 1)
v13_0p4_nafe = v13_0p4 * ((v13_0p4_afe4/v13_0p4 - 1)*(afeval/0.4) + 1)

v3_m0p25_nafe = v3_m0p25 * ((v3_m0p25_afe4/v3_m0p25 - 1)*(afeval/0.4) + 1)
v5_m0p25_nafe = v5_m0p25 * ((v5_m0p25_afe4/v5_m0p25 - 1)*(afeval/0.4) + 1)
v7_m0p25_nafe = v7_m0p25 * ((v7_m0p25_afe4/v7_m0p25 - 1)*(afeval/0.4) + 1)
v9_m0p25_nafe = v9_m0p25 * ((v9_m0p25_afe4/v9_m0p25 - 1)*(afeval/0.4) + 1)
v11_m0p25_nafe = v11_m0p25 * ((v11_m0p25_afe4/v11_m0p25 - 1)*(afeval/0.4) + 1)
v13_m0p25_nafe = v13_m0p25 * ((v13_m0p25_afe4/v13_m0p25 - 1)*(afeval/0.4) + 1)

p.writeto('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MILES_BaSTI_'+imf+'_1.30_fits/M'+imf+'1.30Zp0.06T03.0000_iTp'+str(afeval)+'_Ep'+str(afeval)+'.fits', v3_nafe, v3head, clobber=True)
p.writeto('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MILES_BaSTI_'+imf+'_1.30_fits/M'+imf+'1.30Zp0.06T05.0000_iTp'+str(afeval)+'_Ep'+str(afeval)+'.fits', v5_nafe, v5head, clobber=True)
p.writeto('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MILES_BaSTI_'+imf+'_1.30_fits/M'+imf+'1.30Zp0.06T07.0000_iTp'+str(afeval)+'_Ep'+str(afeval)+'.fits', v7_nafe, v7head, clobber=True)
p.writeto('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MILES_BaSTI_'+imf+'_1.30_fits/M'+imf+'1.30Zp0.06T09.0000_iTp'+str(afeval)+'_Ep'+str(afeval)+'.fits', v9_nafe, v9head, clobber=True)
p.writeto('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MILES_BaSTI_'+imf+'_1.30_fits/M'+imf+'1.30Zp0.06T11.0000_iTp'+str(afeval)+'_Ep'+str(afeval)+'.fits', v11_nafe, v11head, clobber=True)
p.writeto('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MILES_BaSTI_'+imf+'_1.30_fits/M'+imf+'1.30Zp0.06T13.5000_iTp'+str(afeval)+'_Ep'+str(afeval)+'.fits', v13_nafe, v13head, clobber=True)

p.writeto('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MILES_BaSTI_'+imf+'_1.30_fits/M'+imf+'1.30Zp0.26T03.0000_iTp'+str(afeval)+'_Ep'+str(afeval)+'.fits', v3_0p26_nafe, v3head_0p26, clobber=True)
p.writeto('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MILES_BaSTI_'+imf+'_1.30_fits/M'+imf+'1.30Zp0.26T05.0000_iTp'+str(afeval)+'_Ep'+str(afeval)+'.fits', v5_0p26_nafe, v5head_0p26, clobber=True)
p.writeto('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MILES_BaSTI_'+imf+'_1.30_fits/M'+imf+'1.30Zp0.26T07.0000_iTp'+str(afeval)+'_Ep'+str(afeval)+'.fits', v7_0p26_nafe, v7head_0p26, clobber=True)
p.writeto('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MILES_BaSTI_'+imf+'_1.30_fits/M'+imf+'1.30Zp0.26T09.0000_iTp'+str(afeval)+'_Ep'+str(afeval)+'.fits', v9_0p26_nafe, v9head_0p26, clobber=True)
p.writeto('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MILES_BaSTI_'+imf+'_1.30_fits/M'+imf+'1.30Zp0.26T11.0000_iTp'+str(afeval)+'_Ep'+str(afeval)+'.fits', v11_0p26_nafe, v11head_0p26, clobber=True)
p.writeto('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MILES_BaSTI_'+imf+'_1.30_fits/M'+imf+'1.30Zp0.26T13.5000_iTp'+str(afeval)+'_Ep'+str(afeval)+'.fits', v13_0p26_nafe, v13head_0p26, clobber=True)

p.writeto('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MILES_BaSTI_'+imf+'_1.30_fits/M'+imf+'1.30Zp0.40T03.0000_iTp'+str(afeval)+'_Ep'+str(afeval)+'.fits', v3_0p4_nafe, v3head_0p4, clobber=True)
p.writeto('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MILES_BaSTI_'+imf+'_1.30_fits/M'+imf+'1.30Zp0.40T05.0000_iTp'+str(afeval)+'_Ep'+str(afeval)+'.fits', v5_0p4_nafe, v5head_0p4, clobber=True)
p.writeto('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MILES_BaSTI_'+imf+'_1.30_fits/M'+imf+'1.30Zp0.40T07.0000_iTp'+str(afeval)+'_Ep'+str(afeval)+'.fits', v7_0p4_nafe, v7head_0p4, clobber=True)
p.writeto('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MILES_BaSTI_'+imf+'_1.30_fits/M'+imf+'1.30Zp0.40T09.0000_iTp'+str(afeval)+'_Ep'+str(afeval)+'.fits', v9_0p4_nafe, v9head_0p4, clobber=True)
p.writeto('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MILES_BaSTI_'+imf+'_1.30_fits/M'+imf+'1.30Zp0.40T11.0000_iTp'+str(afeval)+'_Ep'+str(afeval)+'.fits', v11_0p4_nafe, v11head_0p4, clobber=True)
p.writeto('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MILES_BaSTI_'+imf+'_1.30_fits/M'+imf+'1.30Zp0.40T13.5000_iTp'+str(afeval)+'_Ep'+str(afeval)+'.fits', v13_0p4_nafe, v13head_0p4, clobber=True)

p.writeto('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MILES_BaSTI_'+imf+'_1.30_fits/M'+imf+'1.30Zm0.25T03.0000_iTp'+str(afeval)+'_Ep'+str(afeval)+'.fits', v3_m0p25_nafe, v3head_m0p25, clobber=True)
p.writeto('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MILES_BaSTI_'+imf+'_1.30_fits/M'+imf+'1.30Zm0.25T05.0000_iTp'+str(afeval)+'_Ep'+str(afeval)+'.fits', v5_m0p25_nafe, v5head_m0p25, clobber=True)
p.writeto('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MILES_BaSTI_'+imf+'_1.30_fits/M'+imf+'1.30Zm0.25T07.0000_iTp'+str(afeval)+'_Ep'+str(afeval)+'.fits', v7_m0p25_nafe, v7head_m0p25, clobber=True)
p.writeto('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MILES_BaSTI_'+imf+'_1.30_fits/M'+imf+'1.30Zm0.25T09.0000_iTp'+str(afeval)+'_Ep'+str(afeval)+'.fits', v9_m0p25_nafe, v9head_m0p25, clobber=True)
p.writeto('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MILES_BaSTI_'+imf+'_1.30_fits/M'+imf+'1.30Zm0.25T11.0000_iTp'+str(afeval)+'_Ep'+str(afeval)+'.fits', v11_m0p25_nafe, v11head_m0p25, clobber=True)
p.writeto('/Users/zieleniewski/Data/stelpopmods/V03/MIUSCAT_2015/MILES_BaSTI_'+imf+'_1.30_fits/M'+imf+'1.30Zm0.25T13.5000_iTp'+str(afeval)+'_Ep'+str(afeval)+'.fits', v13_m0p25_nafe, v13head_m0p25, clobber=True)

print 'New SSPs created!'