'''Code to plot variation of Index and SDSS g-r Colour
for CvD12 models, as functions of age, IMF, alpha enhancement.

By Simon Zieleniewski

Last updated 31-08-16

'''


import spectools as s
import CvD12tools as cvd
import pylab as P
import sys
import numpy as n
import scipy.interpolate as si
import measure_real_indices as mri
import matplotlib.collections as mcoll
import matplotlib.path as mpath
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
def colorline(
    x, y, z=None, cmap=P.cm.copper, norm=P.Normalize(0.0, 1.0),
        linewidth=1.5, alpha=1.0):
    """
    http://nbviewer.ipython.org/github/dpsanders/matplotlib-examples/blob/master/colorline.ipynb
    http://matplotlib.org/examples/pylab_examples/multicolored_line.html
    Plot a colored line with coordinates x and y
    Optionally specify colors in the array z
    Optionally specify a colormap, a norm function and a line width
    """

    # Default colors equally spaced on [0,1]:
    if z is None:
        z = n.linspace(0.0, 1.0, len(x))

    # Special case if a single number:
    if not hasattr(z, "__iter__"):  # to check for numerical input -- this is a hack
        z = n.array([z])

    z = n.asarray(z)

    segments = make_segments(x, y)
    lc = mcoll.LineCollection(segments, array=z, cmap=cmap,
                              linewidth=linewidth, alpha=alpha)

    ax = P.gca()
    ax.add_collection(lc)

    return lc

def make_segments(x, y):
    """
    Create list of line segments from x and y coordinates, in the correct format
    for LineCollection: an array of the form numlines x (points per line) x 2 (x
    and y) array
    """

    points = n.array([x, y]).T.reshape(-1, 1, 2)
    segments = n.concatenate([points[:-1], points[1:]], axis=1)
    return segments

def makeErrorBoxes(xdata,ydata,xerror,yerror,fc='k',ec='None',alpha=0.25):

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
    P.gca().add_collection(pc)
#------------#



#------------#
def varelems(ind):
    if ind in ['NaD', 'NaI', 'NaIsdss']:
        vposm = 1; vposp = 0; labm = '[Na/Fe] = -0.3'; labp = '[Na/Fe] = +0.3'; unit = ' [$\AA$]'
        #Code to plot GMP2921 NaIsdss v CaT
        # vposm = 0; vposp = 0; labm = '[Na/Fe] = +0.3'; labp = '[Na/Fe] = +0.3'; unit = ' [$\AA$]'
    elif ind == 'CaT':
        vposm = 3; vposp = 2; labm = '[Ca/Fe] = -0.15'; labp = '[Ca/Fe] = +0.15'; unit = ' [$\AA$]'
        #Code to plot GMP2921 NaIsdss v CaT
        # vposm = 3; vposp = 2; labm = '[Ca/Fe] = -0.15'; labp = ''; unit = ' [$\AA$]'
    elif ind in ['MgI', 'Mgb', 'Hb']:
        vposm = 15; vposp = 14; labm = '[Mg/Fe] = -0.3'; labp = '[Mg/Fe] = +0.3'; unit = ' [$\AA$]'
    elif ind == 'TiO':
        vposm = 13; vposp = 12; labm = '[Ti/Fe] = -0.3'; labp = '[Ti/Fe] = +0.3'; unit = ''
    elif ind in ['FeH', 'Fe', 'MgFe']:
        vposm = 5; vposp = 4; labm = '[Fe/H] = -0.3'; labp = '[Fe/H] = +0.3'; unit = ' [$\AA$]'

    return vposm, vposp, labm, labp, unit


def indexindexplot(ind1, ind2, comdisp, comadata, fehsys, labels, data='None', gdata='None', data2='None', gdata2='None', data3='None', gdata3='None', data4='None', gdata4='None'):
    comdisp = float(comdisp)

    cvd3 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/CvD12/t03.0_solar.ssp', varspecf=None, kins='CvD12', out_sig=comdisp, index=ind1)
    cvd5 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/CvD12/t05.0_solar.ssp', varspecf=None, kins='CvD12', out_sig=comdisp, index=ind1)
    cvd7 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/CvD12/t07.0_solar.ssp', varspecf=None, kins='CvD12', out_sig=comdisp, index=ind1)
    cvd9 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/CvD12/t09.0_solar.ssp', varspecf=None, kins='CvD12', out_sig=comdisp, index=ind1)
    cvd11 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/CvD12/t11.0_solar.ssp', varspecf=None, kins='CvD12', out_sig=comdisp, index=ind1)
    cvd13 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/CvD12/t13.5_solar.ssp', varspecf=None, kins='CvD12', out_sig=comdisp, index=ind1)
    cvdafe2 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/CvD12/t13.5_afe+0.2.ssp', varspecf=None, kins='CvD12', out_sig=comdisp, index=ind1)
    cvdafe3 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/CvD12/t13.5_afe+0.3.ssp', varspecf=None, kins='CvD12', out_sig=comdisp, index=ind1)
    cvdvar = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/CvD12/t13.5_varelem.ssp', varspecf=None, kins='CvD12', out_sig=comdisp, index=ind1)

    if ind1 in ['NaD', 'NaI', 'NaIsdss', 'MgI', 'FeH', 'CaT', 'PaT', 'Mgb', 'Hb']:
        cvd3ind1 = cvd3.irindex(cvd3.lam[1]-cvd3.lam[0], index=ind1)
        cvd5ind1 = cvd5.irindex(cvd5.lam[1]-cvd5.lam[0], index=ind1)
        cvd7ind1 = cvd7.irindex(cvd7.lam[1]-cvd7.lam[0], index=ind1)
        cvd9ind1 = cvd9.irindex(cvd9.lam[1]-cvd9.lam[0], index=ind1)
        cvd11ind1 = cvd11.irindex(cvd11.lam[1]-cvd11.lam[0], index=ind1)
        cvd13ind1 = cvd13.irindex(cvd13.lam[1]-cvd13.lam[0], index=ind1)
        cvdafe2ind1 = cvdafe2.irindex(cvdafe2.lam[1]-cvdafe2.lam[0], index=ind1)
        cvdafe3ind1 = cvdafe3.irindex(cvdafe3.lam[1]-cvdafe3.lam[0], index=ind1)
        cvdvarind1 = cvdvar.irindex(cvdvar.lam[1]-cvdvar.lam[0], index=ind1)
    elif ind1 == 'TiO':
        cvd3ind1 = cvd3.TiOindex(cvd3.lam[1]-cvd3.lam[0])
        cvd5ind1 = cvd5.TiOindex(cvd5.lam[1]-cvd5.lam[0])
        cvd7ind1 = cvd7.TiOindex(cvd7.lam[1]-cvd7.lam[0])
        cvd9ind1 = cvd9.TiOindex(cvd9.lam[1]-cvd9.lam[0])
        cvd11ind1 = cvd11.TiOindex(cvd11.lam[1]-cvd11.lam[0])
        cvd13ind1 = cvd13.TiOindex(cvd13.lam[1]-cvd13.lam[0])
        cvdafe2ind1 = cvdafe2.TiOindex(cvdafe2.lam[1]-cvdafe2.lam[0])
        cvdafe3ind1 = cvdafe3.TiOindex(cvdafe3.lam[1]-cvdafe3.lam[0])
        cvdvarind1 = cvdvar.TiOindex(cvdvar.lam[1]-cvdvar.lam[0])
    elif ind1 == 'Fe':
        cvd3ind1 = cvd3.Fe(cvd3.lam[1]-cvd3.lam[0])
        cvd5ind1 = cvd5.Fe(cvd5.lam[1]-cvd5.lam[0])
        cvd7ind1 = cvd7.Fe(cvd7.lam[1]-cvd7.lam[0])
        cvd9ind1 = cvd9.Fe(cvd9.lam[1]-cvd9.lam[0])
        cvd11ind1 = cvd11.Fe(cvd11.lam[1]-cvd11.lam[0])
        cvd13ind1 = cvd13.Fe(cvd13.lam[1]-cvd13.lam[0])
        cvdafe2ind1 = cvdafe2.Fe(cvdafe2.lam[1]-cvdafe2.lam[0])
        cvdafe3ind1 = cvdafe3.Fe(cvdafe3.lam[1]-cvdafe3.lam[0])
        cvdvarind1 = cvdvar.Fe(cvdvar.lam[1]-cvdvar.lam[0])
    elif ind1 == 'MgFe':
        cvd3ind1 = cvd3.MgFe(cvd3.lam[1]-cvd3.lam[0])
        cvd5ind1 = cvd5.MgFe(cvd5.lam[1]-cvd5.lam[0])
        cvd7ind1 = cvd7.MgFe(cvd7.lam[1]-cvd7.lam[0])
        cvd9ind1 = cvd9.MgFe(cvd9.lam[1]-cvd9.lam[0])
        cvd11ind1 = cvd11.MgFe(cvd11.lam[1]-cvd11.lam[0])
        cvd13ind1 = cvd13.MgFe(cvd13.lam[1]-cvd13.lam[0])
        cvdafe2ind1 = cvdafe2.MgFe(cvdafe2.lam[1]-cvdafe2.lam[0])
        cvdafe3ind1 = cvdafe3.MgFe(cvdafe3.lam[1]-cvdafe3.lam[0])
        cvdvarind1 = cvdvar.MgFe(cvdvar.lam[1]-cvdvar.lam[0])

    if ind2 in ['NaD', 'NaI', 'NaIsdss', 'MgI', 'FeH', 'CaT', 'PaT', 'Mgb', 'Hb']:
        cvd3ind2 = cvd3.irindex(cvd3.lam[1]-cvd3.lam[0], index=ind2)#/cvd3.Fe(cvd3.lam[1]-cvd3.lam[0])[0]
        cvd5ind2 = cvd5.irindex(cvd5.lam[1]-cvd5.lam[0], index=ind2)#/cvd5.Fe(cvd3.lam[1]-cvd3.lam[0])[0]
        cvd7ind2 = cvd7.irindex(cvd7.lam[1]-cvd7.lam[0], index=ind2)#/cvd7.Fe(cvd3.lam[1]-cvd3.lam[0])[0]
        cvd9ind2 = cvd9.irindex(cvd9.lam[1]-cvd9.lam[0], index=ind2)#/cvd9.Fe(cvd3.lam[1]-cvd3.lam[0])[0]
        cvd11ind2 = cvd11.irindex(cvd11.lam[1]-cvd11.lam[0], index=ind2)#/cvd11.Fe(cvd3.lam[1]-cvd3.lam[0])[0]
        cvd13ind2 = cvd13.irindex(cvd13.lam[1]-cvd13.lam[0], index=ind2)#/cvd13.Fe(cvd3.lam[1]-cvd3.lam[0])[0]
        cvdafe2ind2 = cvdafe2.irindex(cvdafe2.lam[1]-cvdafe2.lam[0], index=ind2)#/cvdafe2.Fe(cvd3.lam[1]-cvd3.lam[0])[0]
        cvdafe3ind2 = cvdafe3.irindex(cvdafe3.lam[1]-cvdafe3.lam[0], index=ind2)#/cvdafe3.Fe(cvd3.lam[1]-cvd3.lam[0])[0]
        cvdvarind2 = cvdvar.irindex(cvdvar.lam[1]-cvdvar.lam[0], index=ind2)#/cvdvar.Fe(cvd3.lam[1]-cvd3.lam[0])[0]
    elif ind2 == 'TiO':
        cvd3ind2 = cvd3.TiOindex(cvd3.lam[1]-cvd3.lam[0])
        cvd5ind2 = cvd5.TiOindex(cvd5.lam[1]-cvd5.lam[0])
        cvd7ind2 = cvd7.TiOindex(cvd7.lam[1]-cvd7.lam[0])
        cvd9ind2 = cvd9.TiOindex(cvd9.lam[1]-cvd9.lam[0])
        cvd11ind2 = cvd11.TiOindex(cvd11.lam[1]-cvd11.lam[0])
        cvd13ind2 = cvd13.TiOindex(cvd13.lam[1]-cvd13.lam[0])
        cvdafe2ind2 = cvdafe2.TiOindex(cvdafe2.lam[1]-cvdafe2.lam[0])
        cvdafe3ind2 = cvdafe3.TiOindex(cvdafe3.lam[1]-cvdafe3.lam[0])
        cvdvarind2 = cvdvar.TiOindex(cvdvar.lam[1]-cvdvar.lam[0])
    elif ind2 == 'Fe':
        cvd3ind2 = cvd3.Fe(cvd3.lam[1]-cvd3.lam[0])
        cvd5ind2 = cvd5.Fe(cvd5.lam[1]-cvd5.lam[0])
        cvd7ind2 = cvd7.Fe(cvd7.lam[1]-cvd7.lam[0])
        cvd9ind2 = cvd9.Fe(cvd9.lam[1]-cvd9.lam[0])
        cvd11ind2 = cvd11.Fe(cvd11.lam[1]-cvd11.lam[0])
        cvd13ind2 = cvd13.Fe(cvd13.lam[1]-cvd13.lam[0])
        cvdafe2ind2 = cvdafe2.Fe(cvdafe2.lam[1]-cvdafe2.lam[0])
        cvdafe3ind2 = cvdafe3.Fe(cvdafe3.lam[1]-cvdafe3.lam[0])
        cvdvarind2 = cvdvar.Fe(cvdvar.lam[1]-cvdvar.lam[0])
    elif ind2 == 'MgFe':
        cvd3ind2 = cvd3.MgFe(cvd3.lam[1]-cvd3.lam[0])
        cvd5ind2 = cvd5.MgFe(cvd5.lam[1]-cvd5.lam[0])
        cvd7ind2 = cvd7.MgFe(cvd7.lam[1]-cvd7.lam[0])
        cvd9ind2 = cvd9.MgFe(cvd9.lam[1]-cvd9.lam[0])
        cvd11ind2 = cvd11.MgFe(cvd11.lam[1]-cvd11.lam[0])
        cvd13ind2 = cvd13.MgFe(cvd13.lam[1]-cvd13.lam[0])
        cvdafe2ind2 = cvdafe2.MgFe(cvdafe2.lam[1]-cvdafe2.lam[0])
        cvdafe3ind2 = cvdafe3.MgFe(cvdafe3.lam[1]-cvdafe3.lam[0])
        cvdvarind2 = cvdvar.MgFe(cvdvar.lam[1]-cvdvar.lam[0])

    #Determine varying element spectra
    vposm, vposp, labm, labp, un1 = varelems(ind1)
    vposm2, vposp2, labm2, labp2, un2 = varelems(ind2)

    #Plot
    P.rcParams['font.family'] = 'Times New Roman'
    P.rcParams.update({'font.size': gensize-3})
    P.figure(figsize=(13,10), dpi=72)
    P.subplots_adjust(left=0.14, bottom=0.12)
    #Varying IMF
    P.scatter(cvd13ind1[0][1:], cvd13ind2[0][1:], marker='o',c='k',s=[90,70,50,30])
    #Varying age
    # P.scatter([cvd3ind1[0][1],cvd5ind1[0][1],cvd7ind1[0][1],cvd9ind1[0][1],cvd11ind1[0][1],cvd13ind1[0][1]],
    #           [cvd3ind2[0][1],cvd5ind2[0][1],cvd7ind2[0][1],cvd9ind2[0][1],cvd11ind2[0][1],cvd13ind2[0][1]],
    #           marker='s',c='r', s=[110,90,70,50,30,10])
    P.scatter([cvd3ind1[0][3],cvd5ind1[0][3],cvd7ind1[0][3],cvd9ind1[0][3],cvd11ind1[0][3],cvd13ind1[0][3]],
              [cvd3ind2[0][3],cvd5ind2[0][3],cvd7ind2[0][3],cvd9ind2[0][3],cvd11ind2[0][3],cvd13ind2[0][3]],
              marker='s',c='r', s=[110,90,70,50,30,10])
    # #Varying [alpha/Fe]
    # P.scatter([cvd13ind1[0][3],cvdafe2ind1[0][3],cvdafe3ind1[0][3]],
    #           [cvd13ind2[0][3],cvdafe2ind2[0][3],cvdafe3ind2[0][3]],
    #           marker='D',c='b',s=[10,30,50])
    #Varying [alpha/Fe]
    P.scatter([cvd13ind1[0][3],cvdafe3ind1[0][3]],
              [cvd13ind2[0][3],cvdafe3ind2[0][3]],
              marker='D',c='b',s=[10,50])
    #Varying element abundance
    P.scatter([cvdvarind1[0][vposm],cvd13ind1[0][3],cvdvarind1[0][vposp]],
              [cvdvarind2[0][vposm],cvd13ind2[0][3],cvdvarind2[0][vposp]],
              marker='D',c='g',s=[50,5,50])
    # P.scatter([cvdvarind1[0][vposm2],cvd13ind1[0][3],cvdvarind1[0][vposp2]],
    #           [cvdvarind2[0][vposm2],cvd13ind2[0][3],cvdvarind2[0][vposp2]],
    #           marker='D',c='g',s=[50,5,50])
    #Varying Na abundance
    # P.scatter([cvd13ind1[0][3],cvdvarind1[0][0]],
    #           [cvd13ind2[0][3],cvdvarind2[0][0]],
    #           marker='D',c='m',s=[50,5,50])
    #Plot lines for all
    P.plot(
           #IMFs
           cvd13ind1[0][1:], cvd13ind2[0][1:], 'k-',
           #Ages
        #    [cvd3ind1[0][1],cvd5ind1[0][1],cvd7ind1[0][1],cvd9ind1[0][1],cvd11ind1[0][1],cvd13ind1[0][1]],
        #    [cvd3ind2[0][1],cvd5ind2[0][1],cvd7ind2[0][1],cvd9ind2[0][1],cvd11ind2[0][1],cvd13ind2[0][1]],'r--',
           [cvd3ind1[0][3],cvd5ind1[0][3],cvd7ind1[0][3],cvd9ind1[0][3],cvd11ind1[0][3],cvd13ind1[0][3]],
           [cvd3ind2[0][3],cvd5ind2[0][3],cvd7ind2[0][3],cvd9ind2[0][3],cvd11ind2[0][3],cvd13ind2[0][3]],'r--',
           #[a/Fe]
           [cvd13ind1[0][3],cvdafe2ind1[0][3],cvdafe3ind1[0][3]],[cvd13ind2[0][3],cvdafe2ind2[0][3],cvdafe3ind2[0][3]],'b-',
           #Varying element abundances
           [cvdvarind1[0][vposm],cvd13ind1[0][3],cvdvarind1[0][vposp]],[cvdvarind2[0][vposm],cvd13ind2[0][3],cvdvarind2[0][vposp]],'g--',
        #    [cvdvarind1[0][vposm2],cvd13ind1[0][3],cvdvarind1[0][vposp2]],[cvdvarind2[0][vposm2],cvd13ind2[0][3],cvdvarind2[0][vposp2]],'g--'
           )
    P.xlabel(ind1+un1, fontsize=tfsize)
    P.ylabel(ind2+un2, fontsize=tfsize)
    if labels=='True':
        # P.annotate('13', xy=(cvd13ind1[0][1],cvd13ind2[0][1]),xycoords='data',xytext=(10,-20),textcoords='offset points', color='r', fontsize=gensize)
        # P.annotate('11', xy=(cvd11ind1[0][1],cvd11ind2[0][1]),xycoords='data',xytext=(-5,5),textcoords='offset points', color='r', fontsize=gensize)
        # P.annotate('9', xy=(cvd9ind1[0][1],cvd9ind2[0][1]),xycoords='data',xytext=(10,0),textcoords='offset points', color='r', fontsize=gensize)
        # P.annotate('7', xy=(cvd7ind1[0][1],cvd7ind2[0][1]),xycoords='data',xytext=(10,0),textcoords='offset points', color='r', fontsize=gensize)
        # P.annotate('5', xy=(cvd5ind1[0][1],cvd5ind2[0][1]),xycoords='data',xytext=(15,-5),textcoords='offset points', color='r', fontsize=gensize)
        # P.annotate('3 Gyr', xy=(cvd3ind1[0][1],cvd3ind2[0][1]),xycoords='data',xytext=(10,0),textcoords='offset points', color='r', fontsize=gensize)
        P.annotate('13', xy=(cvd13ind1[0][3],cvd13ind2[0][3]),xycoords='data',xytext=(10,-15),textcoords='offset points', color='r', fontsize=gensize)
        P.annotate('11', xy=(cvd11ind1[0][3],cvd11ind2[0][3]),xycoords='data',xytext=(-5,-20),textcoords='offset points', color='r', fontsize=gensize)
        P.annotate('9', xy=(cvd9ind1[0][3],cvd9ind2[0][3]),xycoords='data',xytext=(10,-5),textcoords='offset points', color='r', fontsize=gensize)
        P.annotate('7', xy=(cvd7ind1[0][3],cvd7ind2[0][3]),xycoords='data',xytext=(10,0),textcoords='offset points', color='r', fontsize=gensize)
        P.annotate('5', xy=(cvd5ind1[0][3],cvd5ind2[0][3]),xycoords='data',xytext=(15,-5),textcoords='offset points', color='r', fontsize=gensize)
        P.annotate('3 Gyr', xy=(cvd3ind1[0][3],cvd3ind2[0][3]),xycoords='data',xytext=(10,0),textcoords='offset points', color='r', fontsize=gensize)
        P.annotate('$x = 3$', xy=(cvd13ind1[0][1],cvd13ind2[0][1]),xycoords='data',xytext=(5,-5),textcoords='offset points', color='k', fontsize=gensize)
        P.annotate('$x = 2.35$', xy=(cvd13ind1[0][2],cvd13ind2[0][2]),xycoords='data',xytext=(-5,5),textcoords='offset points', color='k', fontsize=gensize)
        P.annotate('Chab', xy=(cvd13ind1[0][3],cvd13ind2[0][3]),xycoords='data',xytext=(-25,10),textcoords='offset points', color='k', fontsize=gensize)
        P.annotate('B-L', xy=(cvd13ind1[0][4],cvd13ind2[0][4]),xycoords='data',xytext=(-25,-20),textcoords='offset points', color='k', fontsize=gensize)
        # P.annotate('[$\\alpha$/Fe]=+0.2', xy=(cvdafe2ind1[0][3],cvdafe2ind2[0][3]),xycoords='data',xytext=(-125,-10),textcoords='offset points', color='b', fontsize=gensize)
        if ind1=='NaIsdss':
            P.annotate('[$\\alpha$/Fe]=+0.3', xy=(cvdafe3ind1[0][3],cvdafe3ind2[0][3]),xycoords='data',xytext=(-70,5),textcoords='offset points', color='b', fontsize=gensize)
        else:
            P.annotate('[$\\alpha$/Fe]=+0.3', xy=(cvdafe3ind1[0][3],cvdafe3ind2[0][3]),xycoords='data',xytext=(-110,0),textcoords='offset points', color='b', fontsize=gensize)
        # P.annotate('& [Z/H]=+0.28', xy=(cvdafe3ind1[0][3],cvdafe3ind2[0][3]),xycoords='data',xytext=(-120,-18),textcoords='offset points', color='b', fontsize=gensize)

        # P.annotate(labm, xy=(cvdvarind1[0][vposm],cvdvarind2[0][vposm2]),xycoords='data',xytext=(-40,-25),textcoords='offset points', color='g', fontsize=gensize)
        # P.annotate(labp, xy=(cvdvarind1[0][vposp],cvdvarind2[0][vposp2]),xycoords='data',xytext=(-80,-25),textcoords='offset points', color='g', fontsize=gensize)
        # P.annotate(labm2, xy=(cvdvarind1[0][vposm2],cvdvarind2[0][vposm2]),xycoords='data',xytext=(-60,-15),textcoords='offset points', color='g', fontsize=gensize)
        # P.annotate(labp2, xy=(cvdvarind1[0][vposp2],cvdvarind2[0][vposp2]),xycoords='data',xytext=(-60,15),textcoords='offset points', color='g', fontsize=gensize)
        P.annotate(labm, xy=(cvdvarind1[0][vposm],cvdvarind2[0][vposm]),xycoords='data',xytext=(-60,5),textcoords='offset points', color='g', fontsize=gensize)
        P.annotate(labp, xy=(cvdvarind1[0][vposp],cvdvarind2[0][vposp]),xycoords='data',xytext=(-60,-20),textcoords='offset points', color='g', fontsize=gensize)
        # if ind1 and ind2 not in ['NaI', 'NaIsdss', 'NaD']:
        #     P.annotate(labm2, xy=(cvdvarind1[0][vposm],cvdvarind2[0][vposm2]),xycoords='data',xytext=(-40,-45),textcoords='offset points', color='g', fontsize=gensize)
        #     P.annotate(labp2, xy=(cvdvarind1[0][vposp],cvdvarind2[0][vposp2]),xycoords='data',xytext=(-80,-45),textcoords='offset points', color='g', fontsize=gensize)

    ### Arrows for [Z/H] from V15 models
        if ind2 == 'TiO':
            if comdisp == 270.0:
                P.annotate("", xy=(0.15, 1.036), xytext=(0.15, 1.040),
                arrowprops=dict(width=2, frac=0.2, headwidth=10, color='k',))
                P.annotate('$\delta$[Z/H] = +0.2', xy=(0.1,1.040),xycoords='data',xytext=(5,0),textcoords='offset points', fontsize=gensize-4)
            else:
                P.annotate("", xy=(0.25, 1.038), xytext=(0.2, 1.042),
                arrowprops=dict(width=2, frac=0.2, headwidth=10, color='k',))
                P.annotate('$\delta$[Z/H] = +0.2', xy=(0.25,1.0375),xycoords='data',xytext=(-20,-12),textcoords='offset points', fontsize=gensize-4)
                P.annotate("", xy=(0.14, 1.0385), xytext=(0.2, 1.042),
                arrowprops=dict(width=2, frac=0.2, headwidth=10, color='k',))
                P.annotate('$\delta$[Na/Fe] = +0.3', xy=(0.14,1.038),xycoords='data',xytext=(-65,-12),textcoords='offset points', fontsize=gensize-4)
        if ind2 == 'MgI':
            # if comdisp == 185.0:
            #     ymin=0.64; ymax=0.685; xmin=0.25; xmax=.387
            if comdisp == 200.0:
                ymin=0.6; ymax=0.64; xmin=0.35; xmax=.487
            # elif comdisp == 270.0:
            #     ymin=0.5; ymax=0.535; xmin=0.25; xmax=.37
            # elif comdisp == 400.0:
            #     ymin=0.25; ymax=0.28; xmin=0.25; xmax=.35
            P.annotate("", xy=(xmax, ymax), xytext=(xmin, ymin),
            arrowprops=dict(width=2, frac=0.2, headwidth=10, color='k',))
            P.annotate('$\delta$[Z/H] = +0.2', xy=(xmin,ymin),xycoords='data',xytext=(30,-10),textcoords='offset points', fontsize=gensize-4)
        if ind1 == 'Mgb':
            if comdisp == 200.0:
                ymin=3.2; ymax=3.55; xmin=4.5; xmax=4.94
            P.annotate("", xy=(xmax, ymax), xytext=(xmin, ymin),
            arrowprops=dict(width=2, frac=0.1, headwidth=10, color='k',))
            P.annotate('$\delta$[Z/H] = +0.2', xy=(xmin,ymin),xycoords='data',xytext=(15,-10),textcoords='offset points', fontsize=gensize-4)
    ###

    #Code to plot GMP2921 NaIsdss v CaT
    # P.annotate(labm, xy=(cvdvarind1[0][vposm],cvdvarind2[0][vposm2]),xycoords='data',xytext=(-20,-25),textcoords='offset points', color='g', fontsize=gensize)
    # P.annotate(labm2, xy=(cvdvarind1[0][vposm],cvdvarind2[0][vposm2]),xycoords='data',xytext=(-20,-45),textcoords='offset points', color='g', fontsize=gensize)
    # P.annotate(labp, xy=(cvdvarind1[0][vposp],cvd13ind2[0][3]),xycoords='data',xytext=(-40,10),textcoords='offset points', color='g', fontsize=gensize)
    # P.annotate(labp2, xy=(cvdvarind1[0][vposp],cvd13ind2[0][3]),xycoords='data',xytext=(-40,-45),textcoords='offset points', color='g', fontsize=gensize)
    P.xticks(size=lfsize)
    P.yticks(size=lfsize)
    P.minorticks_on()
    P.tick_params(axis='both', direction='in', length=majtsize,
                  width=1.5, which='major', pad=10)
    P.tick_params(axis='both', direction='in', length=mintsize,
                  width=1.5, which='minor', pad=10)
    #Plot data from file if given
    gmpnames = ['GMP2921', 'GMP3329', 'GMP4928', 'GMP3367', 'M31']
    ngcnames = ['NGC4889', 'NGC4874', 'NGC4839', 'NGC4873', 'M31']
    gals = []
    plots = []
    ggals = []
    gplots = []
    indpos = {'NaI':1, 'NaIsdss':1, 'CaT':4, 'MgI':7, 'TiO':10, 'FeH':13}
    gindpos = {'NaI':0, 'NaIsdss':0, 'CaT':2, 'MgI':4, 'TiO':6, 'FeH':8}
    if data != 'None':
        gal = data.split('/')[-1].split('_')[3]
        sig = data.split('/')[-1].split('_')[4].split('kms')[0]
        dat = n.genfromtxt(data, delimiter=',')
        rs = dat[:,0]
        ind1vals = dat[:,indpos[ind1]]
        ind1vars = dat[:,indpos[ind1]+1]
        ind2vals = dat[:,indpos[ind2]]
        ind2vars = dat[:,indpos[ind2]+1]
        if gal == 'M31':
            ind1vars = ind1vars**2
            ind2vars = ind2vars**2
            rs = n.log10(rs)
        if (gal=='GMP3329' and ind1 == 'FeH') or (gal=='GMP3329' and ind2=='FeH'):
            ind1vals = n.concatenate((ind1vals[0:-2],n.array([ind1vals[-1]])))
            ind1vars = n.concatenate((ind1vars[0:-2],n.array([ind1vars[-1]])))
            ind2vals = n.concatenate((ind2vals[0:-2],n.array([ind2vals[-1]])))
            ind2vars = n.concatenate((ind2vars[0:-2],n.array([ind2vars[-1]])))
            rs = n.concatenate((rs[0:-2],n.array([rs[-1]])))
        if (gal=='GMP4928' and ind1 in ['NaIsdss', 'CaT']) or (gal=='GMP4928' and ind2 in ['NaIsdss', 'CaT']):
            print 'Ignoring outermost GMP4928 datapoint'
            ind1vals = ind1vals[0:-1]
            ind1vars = ind1vars[0:-1]
            ind2vals = ind2vals[0:-1]
            ind2vars = ind2vars[0:-1]
            rs = rs[0:-1]
        if (gal=='GMP3329' and ind1 in ['NaIsdss', 'CaT']) or (gal=='GMP3329' and ind2 in ['NaIsdss', 'CaT']):
            print 'Ignoring outermost GMP3329 datapoint'
            ind1vals = ind1vals[0:-1]
            ind1vars = ind1vars[0:-1]
            ind2vals = ind2vals[0:-1]
            ind2vars = ind2vars[0:-1]
            rs = rs[0:-1]
        print ind1vals
        print ind1vars
        print ind2vals
        print ind2vars
        print 'RS = ', rs
        markers = [220 -i*20. for i in xrange(len(ind1vals))]

        if data2 =='None' and data3 =='None' and data4=='None':
            plot = P.scatter(ind1vals, ind2vals, marker='o', s=[230,200,170,140,110,80,60,40,10,5], c=rs, linewidths=0, cmap=P.cm.copper, zorder=10)
            clb = P.colorbar(plot)
            clb.set_label('r [arcsec]', fontsize=tfsize)
            if gal=='M31':
                clb.set_label('log$_{10}$(r/arcsec)', fontsize=tfsize)
            #convert time to a color tuple using the colormap used for scatter
            r_colour = clb.to_rgba(rs)
            #pm0.05 RANDOM ERROR ADDED IN QUADRATURE - 28-05-15 - COMPUTE THIS FOR COMA GALAXIES!!! - 12-11-15
            ind1errfac = 0.00
            ind2errfac = 0.00
            ind1vars = ind1vars+ind1errfac**2
            ind2vars = ind2vars+ind2errfac**2
            a,b,c = P.errorbar(ind1vals, ind2vals, xerr=n.sqrt(ind1vars), yerr=n.sqrt(ind2vars),
                            marker='', ls='', capsize=0, barsabove=False, capthick=1.0,
                            elinewidth=1.4, zorder=10)
            lines = colorline(ind1vals, ind2vals, rs)
            #adjust the color of c[0], which is a LineCollection, to the colormap
            c[0].set_color(r_colour)
            c[1].set_color(r_colour)

            gals.append(ngcnames[gmpnames.index(gal)]+' ('+sig+' km s$^{-1}$)')
            plots.append(plot)
        else:
            plot = P.scatter(ind1vals, ind2vals, marker='', s=[0,0,0,0,0,0,0,0,0,0], c=rs, linewidths=0, cmap=P.cm.copper, zorder=0)
            clb = P.colorbar(plot)
            clb.set_label('r [arcsec]', fontsize=tfsize)
            if gal=='M31':
                clb.set_label('log$_{10}$(r/arcsec)', fontsize=tfsize)
            #convert time to a color tuple using the colormap used for scatter
            r_colour = clb.to_rgba(rs)

    #plot global index value if given
    if gdata != 'None':
        ggal = gdata.split('/')[-1].split('_')[0]
        sig = gdata.split('/')[-1].split('_')[5].split('kms')[0]
        ggals.append(ngcnames[gmpnames.index(ggal)]+' ('+sig+' km s$^{-1}$)')
        gdat = n.genfromtxt(gdata, delimiter=',')
        gind1vals = gdat[gindpos[ind1]]
        gind1errs = gdat[gindpos[ind1]+1]
        gind2vals = gdat[gindpos[ind2]]
        gind2errs = gdat[gindpos[ind2]+1]
        gpl = P.errorbar(gind1vals, gind2vals, xerr=gind1errs, yerr=gind2errs,
        fmt='o', ms=msize+4, mec='k', mew=1.5, mfc='w', ls='', ecolor='k',
        capsize=0, barsabove=False, capthick=1.0, elinewidth=1.6, zorder=10)
        gplots.append(gpl)

    #Plot data from second file if given
    if data2 != 'None':
        gal2 = data2.split('/')[-1].split('_')[3]
        sig2 = data2.split('/')[-1].split('_')[4].split('kms')[0]
        dat2 = n.genfromtxt(data2, delimiter=',')
        rs2 = dat2[:,0]
        ind1vals2 = dat2[:,indpos[ind1]]
        ind1vars2 = dat2[:,indpos[ind1]+1]
        ind2vals2 = dat2[:,indpos[ind2]]
        ind2vars2 = dat2[:,indpos[ind2]+1]
        if (gal2=='GMP3329' and ind1 == 'FeH') or (gal2=='GMP3329' and ind2=='FeH'):
            ind1vals2 = n.concatenate((ind1vals2[0:-2],n.array([ind1vals2[-1]])))
            ind1vars2 = n.concatenate((ind1vars2[0:-2],n.array([ind1vars2[-1]])))
            ind2vals2 = n.concatenate((ind2vals2[0:-2],n.array([ind2vals2[-1]])))
            ind2vars2 = n.concatenate((ind2vars2[0:-2],n.array([ind2vars2[-1]])))
            rs2 = n.concatenate((rs2[0:-2],n.array([rs2[-1]])))
        if (gal2=='GMP4928' and ind1 in ['NaIsdss', 'CaT']) or (gal2=='GMP4928' and ind2 in ['NaIsdss', 'CaT']):
            print 'Ignoring outermost GMP4928 datapoint'
            ind1vals2 = ind1vals2[:-1]
            ind1vars2 = ind1vars2[:-1]
            ind2vals2 = ind2vals2[:-1]
            ind2vars2 = ind2vars2[:-1]
            rs2 = rs2[:-1]
        if (gal2=='GMP3329' and ind1 in ['NaIsdss', 'CaT']) or (gal2=='GMP3329' and ind2 in ['NaIsdss', 'CaT']):
            print 'Ignoring outermost GMP3329 datapoint'
            ind1vals2 = ind1vals2[0:-1]
            ind1vars2 = ind1vars2[0:-1]
            ind2vals2 = ind2vals2[0:-1]
            ind2vars2 = ind2vars2[0:-1]
            rs2 = rs2[0:-1]
        print ind1vals2
        print ind1vars2
        print ind2vals2
        print ind2vars2
        print 'RS2 = ', rs2
        markers = [220 -i*20. for i in xrange(len(ind1vals2))]
        r_colour2 = clb.to_rgba(rs2)
        plot2 = P.scatter(ind1vals2, ind2vals2, marker='s', s=[230,200,170,140,110,80,60,40,10,5], c=r_colour2, linewidths=0, cmap=P.cm.copper_r, zorder=10)
        #convert time to a color tuple using the colormap used for scatter
        # r_colour2 = clb.to_rgba(rs2)
        #pm0.05 RANDOM ERROR ADDED IN QUADRATURE - 28-05-15 - COMPUTE THIS FOR COMA GALAXIES!!! - 12-11-15
        ind1errfac = 0.00
        ind2errfac = 0.00
        ind1vars2 = ind1vars2+ind1errfac**2
        ind2vars2 = ind2vars2+ind2errfac**2
        a2,b2,c2 = P.errorbar(ind1vals2, ind2vals2, xerr=n.sqrt(ind1vars2), yerr=n.sqrt(ind2vars2),
                           marker='', ls='', capsize=0, barsabove=False, capthick=1.0,
                           elinewidth=1.4, zorder=10)
        lines2 = colorline(ind1vals2, ind2vals2, rs2)

        #adjust the color of c[0], which is a LineCollection, to the colormap
        c2[0].set_color(r_colour2)
        c2[1].set_color(r_colour2)
        plots.append(plot2)
        gals.append(ngcnames[gmpnames.index(gal2)]+' ('+sig2+' km s$^{-1}$)')

    if gdata2 != 'None':
        ggal2 = gdata2.split('/')[-1].split('_')[0]
        sig = gdata2.split('/')[-1].split('_')[5].split('kms')[0]
        ggals.append(ngcnames[gmpnames.index(ggal2)]+' ('+sig+' km s$^{-1}$)')
        gdat2 = n.genfromtxt(gdata2, delimiter=',')
        gind1vals2 = gdat2[gindpos[ind1]]
        gind1errs2 = gdat2[gindpos[ind1]+1]
        gind2vals2 = gdat2[gindpos[ind2]]
        gind2errs2 = gdat2[gindpos[ind2]+1]
        gpl2 = P.errorbar(gind1vals2, gind2vals2, xerr=gind1errs2, yerr=gind2errs2,
        fmt='s', ms=msize+4, mec='k', mew=1.5, mfc='w', ls='', ecolor='k',
        capsize=0, barsabove=False, capthick=1.0, elinewidth=1.6, zorder=10)
        gplots.append(gpl2)

    if data3 != 'None':
        gal3 = data3.split('/')[-1].split('_')[3]
        sig3 = data3.split('/')[-1].split('_')[4].split('kms')[0]
        dat3 = n.genfromtxt(data3, delimiter=',')
        rs3 = dat3[:,0]
        ind1vals3 = dat3[:,indpos[ind1]]
        ind1vars3 = dat3[:,indpos[ind1]+1]
        ind2vals3 = dat3[:,indpos[ind2]]
        ind2vars3 = dat3[:,indpos[ind2]+1]
        if gal3 == 'M31':
            ind1vars3 = ind1vars3**2
            ind2vars3 = ind2vars3**2
            rs3 = n.log10(rs3)
        if (gal3=='GMP3329' and ind1 == 'FeH') or (gal3=='GMP3329' and ind2=='FeH'):
            ind1vals3 = n.concatenate((ind1vals3[0:-2],n.array([ind1vals3[-1]])))
            ind1vars3 = n.concatenate((ind1vars3[0:-2],n.array([ind1vars3[-1]])))
            ind2vals3 = n.concatenate((ind2vals3[0:-2],n.array([ind2vals3[-1]])))
            ind2vars3 = n.concatenate((ind2vars3[0:-2],n.array([ind2vars3[-1]])))
            rs3 = n.concatenate((rs3[0:-2],n.array([rs3[-1]])))
        if (gal3=='GMP4928' and ind1 in ['NaIsdss', 'CaT']) or (gal3=='GMP4928' and ind2 in ['NaIsdss', 'CaT']):
            print 'Ignoring outermost GMP4928 datapoint'
            ind1vals3 = ind1vals3[0:-1]
            ind1vars3 = ind1vars3[0:-1]
            ind2vals3 = ind2vals3[0:-1]
            ind2vars3 = ind2vars3[0:-1]
            rs3 = rs3[0:-1]
        if (gal3=='GMP3329' and ind1 in ['NaIsdss', 'CaT']) or (gal3=='GMP3329' and ind2 in ['NaIsdss', 'CaT']):
            print 'Ignoring outermost GMP3329 datapoint'
            ind1vals3 = ind1vals3[0:-1]
            ind1vars3 = ind1vars3[0:-1]
            ind2vals3 = ind2vals3[0:-1]
            ind2vars3 = ind2vars3[0:-1]
            rs3 = rs3[0:-1]
        print ind1vals3
        print ind1vars3
        print ind2vals3
        print ind2vars3
        print 'RS = ', rs3
        markers = [220 -i*20. for i in xrange(len(ind1vals3))]
        r_colour3 = clb.to_rgba(rs3)
        plot3 = P.scatter(ind1vals3, ind2vals3, marker='D', s=[230,200,170,140,110,80,60,40,10,5], c=rs3, linewidths=0, cmap=P.cm.copper, zorder=10)
        #pm0.05 RANDOM ERROR ADDED IN QUADRATURE - 28-05-15 - COMPUTE THIS FOR COMA GALAXIES!!! - 12-11-15
        ind1errfac = 0.00
        ind2errfac = 0.00
        ind1vars3 = ind1vars3+ind1errfac**2
        ind2vars3 = ind2vars3+ind2errfac**2
        a3,b3,c3 = P.errorbar(ind1vals3, ind2vals3, xerr=n.sqrt(ind1vars3), yerr=n.sqrt(ind2vars3),
                           marker='', ls='', capsize=0, barsabove=False, capthick=1.0,
                           elinewidth=1.4, zorder=10)
        lines3 = colorline(ind1vals3, ind2vals3, rs3)

        #adjust the color of c[0], which is a LineCollection, to the colormap
        c3[0].set_color(r_colour3)
        c3[1].set_color(r_colour3)
        gals.append(ngcnames[gmpnames.index(gal3)]+' ('+sig3+' km s$^{-1}$)')
        plots.append(plot3)

    #plot global index value if given
    if gdata3 != 'None':
        ggal3 = gdata3.split('/')[-1].split('_')[0]
        sig = gdata3.split('/')[-1].split('_')[5].split('kms')[0]
        ggals.append(ngcnames[gmpnames.index(ggal3)]+' ('+sig+' km s$^{-1}$)')
        gdat3 = n.genfromtxt(gdata3, delimiter=',')
        gind1vals3 = gdat3[gindpos[ind1]]
        gind1errs3 = gdat3[gindpos[ind1]+1]
        gind2vals3 = gdat3[gindpos[ind2]]
        gind2errs3 = gdat3[gindpos[ind2]+1]
        gpl3 = P.errorbar(gind1vals3, gind2vals3, xerr=gind1errs3, yerr=gind2errs3,
        fmt='D', ms=msize+4, mec='k', mew=1.5, mfc='w', ls='', ecolor='k',
        capsize=0, barsabove=False, capthick=1.0, elinewidth=1.6, zorder=10)
        gplots.append(gpl3)

    if data4 != 'None':
        gal4 = data4.split('/')[-1].split('_')[3]
        sig4 = data4.split('/')[-1].split('_')[4].split('kms')[0]
        dat4 = n.genfromtxt(data4, delimiter=',')
        rs4 = dat4[:,0]
        ind1vals4 = dat4[:,indpos[ind1]]
        ind1vars4 = dat4[:,indpos[ind1]+1]
        ind2vals4 = dat4[:,indpos[ind2]]
        ind2vars4 = dat4[:,indpos[ind2]+1]
        if gal4 == 'M31':
            ind1vars4 = ind1vars4**2
            ind2vars4 = ind2vars4**2
            rs4 = n.log10(rs4)
        if (gal4=='GMP3329' and ind1 == 'FeH') or (gal4=='GMP3329' and ind2=='FeH'):
            ind1vals4 = n.concatenate((ind1vals4[0:-2],n.array([ind1vals4[-1]])))
            ind1vars4 = n.concatenate((ind1vars4[0:-2],n.array([ind1vars4[-1]])))
            ind2vals4 = n.concatenate((ind2vals4[0:-2],n.array([ind2vals4[-1]])))
            ind2vars4 = n.concatenate((ind2vars4[0:-2],n.array([ind2vars4[-1]])))
            rs4 = n.concatenate((rs4[0:-2],n.array([rs4[-1]])))
        if (gal4=='GMP4928' and ind1 in ['NaIsdss', 'CaT']) or (gal4=='GMP4928' and ind2 in ['NaIsdss', 'CaT']):
            print 'Ignoring outermost GMP4928 datapoint'
            ind1vals4 = ind1vals4[0:-1]
            ind1vars4 = ind1vars4[0:-1]
            ind2vals4 = ind2vals4[0:-1]
            ind2vars4 = ind2vars4[0:-1]
            rs4 = rs4[0:-1]
        if (gal4=='GMP3329' and ind1 in ['NaIsdss', 'CaT']) or (gal4=='GMP3329' and ind2 in ['NaIsdss', 'CaT']):
            print 'Ignoring outermost GMP3329 datapoint'
            ind1vals4 = ind1vals4[0:-1]
            ind1vars4 = ind1vars4[0:-1]
            ind2vals4 = ind2vals4[0:-1]
            ind2vars4 = ind2vars4[0:-1]
            rs4 = rs4[0:-1]
        print ind1vals4
        print ind1vars4
        print ind2vals4
        print ind2vars4
        print 'RS = ', rs4
        markers = [220 -i*20. for i in xrange(len(ind1vals4))]
        r_colour4 = clb.to_rgba(rs4)
        plot4 = P.scatter(ind1vals4, ind2vals4, marker='*', s=[230,200,170,140,110,80,60,40,10,5], c=rs4, linewidths=0, cmap=P.cm.copper, zorder=10)
        #pm0.05 RANDOM ERROR ADDED IN QUADRATURE - 28-05-15 - COMPUTE THIS FOR COMA GALAXIES!!! - 12-11-15
        ind1errfac = 0.00
        ind2errfac = 0.00
        ind1vars4 = ind1vars4+ind1errfac**2
        ind2vars4 = ind2vars4+ind2errfac**2
        a4,b4,c4 = P.errorbar(ind1vals4, ind2vals4, xerr=n.sqrt(ind1vars4), yerr=n.sqrt(ind2vars4),
                           marker='', ls='', capsize=0, barsabove=False, capthick=1.0,
                           elinewidth=1.4, zorder=10)
        lines4 = colorline(ind1vals4, ind2vals4, rs4)

        #adjust the color of c[0], which is a LineCollection, to the colormap
        c4[0].set_color(r_colour4)
        c4[1].set_color(r_colour4)
        gals.append(ngcnames[gmpnames.index(gal4)]+' ('+sig4+' km s$^{-1}$)')
        plots.append(plot4)

    #plot global index value if given
    if gdata4 != 'None':
        ggal4 = gdata4.split('/')[-1].split('_')[0]
        sig = gdata4.split('/')[-1].split('_')[5].split('kms')[0]
        ggals.append(ngcnames[gmpnames.index(ggal4)]+' ('+sig+' km s$^{-1}$)')
        gdat4 = n.genfromtxt(gdata4, delimiter=',')
        gind1vals4 = gdat4[gindpos[ind1]]
        gind1errs4 = gdat4[gindpos[ind1]+1]
        gind2vals4 = gdat4[gindpos[ind2]]
        gind2errs4 = gdat4[gindpos[ind2]+1]
        gpl4 = P.errorbar(gind1vals4, gind2vals4, xerr=gind1errs4, yerr=gind2errs4,
        fmt='*', ms=msize+6, mec='k', mew=1.5, mfc='w', ls='', ecolor='k',
        capsize=0, barsabove=False, capthick=1.0, elinewidth=1.5, zorder=10)
        gplots.append(gpl4)

    if data != 'None' or data2 != 'None' or data3 != 'None' or data4 != 'None':
        if ind1=='NaIsdss':
            lpos = 'lower left'
        else:
            lpos = 'lower right'
        if gdata4 != 'None':
            leg = P.legend(n.concatenate((plots,gplots[-1])), n.concatenate((gals,[ggals[-1]])), loc=lpos, frameon=False,
            borderaxespad=0.2, labelspacing=0.1, numpoints=1, scatterpoints=1)
            for i in xrange(len(leg.legendHandles)):
                leg.legendHandles[i].set_color(r_colour[0])
        else:
            leg = P.legend(plots, gals, loc=lpos, frameon=False,
            borderaxespad=0.2, labelspacing=0.1, numpoints=1, scatterpoints=1)
            for i in xrange(len(leg.legendHandles)):
                leg.legendHandles[i].set_color(r_colour[0])

    if (gdata != 'None' and data == 'None') or (gdata2 != 'None' and data2 == 'None') or (gdata3 != 'None' and data3 == 'None') or (gdata3 != 'None' and data3 == 'None'):
        gleg = P.legend(gplots, ggals, loc='lower right', frameon=False,
        borderaxespad=0.3,labelspacing=0.2, numpoints=1, scatterpoints=1)

    if ind1 == 'FeH':
        P.xlim([0.0,1.0])
        P.ylim([1.03,1.09])
    if ind1 == 'NaIsdss':
        P.xlim([0.3,1.0])
        P.ylim([0.2,0.8])
    if ind2 == 'CaT':
        P.xlim([0.2,1.0])
        P.ylim([6.5,9.0])       

    #NaIsdss v NaD data: NaD measured from SDSS 3" spectra
    # gind2921nad = 5.178; gind2921naderr = 0.117
    # gind3329nad = 4.413; gind3329naderr = 0.136
    # gind4928nad = 4.663; gind4928naderr = 0.117
    # gind3367nad = 3.680; gind3367naderr = 0.132
    #NaIsdss from SWIFT data
    # gind2921nai = 0.712; gind2921naierr = 0.049
    # gind3329nai = 0.410; gind3329naierr = 0.096
    # gind4928nai = 0.640; gind4928naierr = 0.072
    # gind3367nai = 0.572; gind3367naierr = 0.068
    #NaIsdss from SDSS data
    # gind2921nai = 0.643; gind2921naierr = 0.135
    # gind3329nai = 0.742; gind3329naierr = 0.143
    # gind4928nai = 0.756; gind4928naierr = 0.125
    # gind3367nai = 0.413; gind3367naierr = 0.150

    # pl2921 = P.errorbar(gind2921nai, gind2921nad, xerr=gind2921naierr, yerr=gind2921naderr,
    # fmt='o', ms=msize+4, mec='k', mew=1.2, mfc='w', ls='', ecolor='k',
    # capsize=0, barsabove=False, capthick=1.0, elinewidth=1.4, zorder=10)
    # leg = P.legend([pl2921], ['GMP2921 (400 km s$^{-1}$)'], loc='lower right',
    # frameon=False, borderaxespad=0.6, labelspacing=0.2, numpoints=1, scatterpoints=1)

    # pl3329 = P.errorbar(gind3329nai, gind3329nad, xerr=gind3329naierr, yerr=gind3329naderr,
    # fmt='o', ms=msize+4, mec='k', mew=1.2, mfc='w', ls='', ecolor='k',
    # capsize=0, barsabove=False, capthick=1.0, elinewidth=1.4, zorder=10)
    # pl4928 = P.errorbar(gind4928nai, gind4928nad, xerr=gind4928naierr, yerr=gind4928naderr,
    # fmt='s', ms=msize+4, mec='k', mew=1.2, mfc='w', ls='', ecolor='k',
    # capsize=0, barsabove=False, capthick=1.0, elinewidth=1.4, zorder=10)
    # leg = P.legend([pl3329, pl4928], ['GMP3329 (270 km s$^{-1}$)', 'GMP4928 (270 km s$^{-1}$)'],
    # loc='lower right', frameon=False, borderaxespad=0.6, labelspacing=0.2, numpoints=1, scatterpoints=1)

    # pl3367 = P.errorbar(gind3367nai, gind3367nad, xerr=gind3367naierr, yerr=gind3367naderr,
    # fmt='o', ms=msize+4, mec='k', mew=1.2, mfc='w', ls='', ecolor='k',
    # capsize=0, barsabove=False, capthick=1.0, elinewidth=1.4, zorder=10)
    # leg = P.legend([pl3367], ['GMP3367 (200 km s$^{-1}$)'], loc='lower right',
    # frameon=False, borderaxespad=0.6, labelspacing=0.2, numpoints=1, scatterpoints=1)

##### Plot Coma BCG data from Loubser et al. 2009
    if comadata=='True':
        sys.path.append('/Users/zieleniewski/Documents/Oxford/D.Phil/My_Papers/Coma/')
        from Coma_plots import Fe, mgfe
        #First do Trager2008 data
        galname = []
        ref = []
        sig = []
        esig = []
        Hb = []
        eHb = []
        Mgb = []
        eMgb = []
        Fe52 = []
        eFe52 = []
        Fe53 = []
        eFe53 = []
        with open('/Users/zieleniewski/Data/swift/Coma_stellar_pops/Trager2008_indices.tsv') as f:
            next(f)
            for line in f:
                data = line.split(';')
                galname.append(data[0].split( )[0])
                ref.append(data[1])
                sig.append(float(data[2]))
                esig.append(float(data[3]))
                Hb.append(float(data[4]))
                eHb.append(float(data[5]))
                Mgb.append(float(data[6]))
                eMgb.append(float(data[7]))
                Fe52.append(float(data[8]))
                eFe52.append(float(data[9]))
                Fe53.append(float(data[10]))
                eFe53.append(float(data[11]))
        
        gmp2921s = []
        gmp3329s = []
        gmp4928s = []
        gmp3367s = []
        for i in xrange(len(galname)):
            if galname[i] == 'GMP2921':
                gmp2921s.append(i)
            if galname[i] == 'GMP3329':
                gmp3329s.append(i)
            if galname[i] == 'GMP4928':
                gmp4928s.append(i)
            if galname[i] == 'GMP3367':
                gmp3367s.append(i)
        
        print gmp2921s
        
        avFe, eavFe = Fe(n.array(Fe52), n.array(eFe52),n.array(Fe53), n.array(eFe53))
        MgFe, eMgFe = mgfe(n.array(Mgb), n.array(eMgb), n.array(Fe52), n.array(eFe52),n.array(Fe53), n.array(eFe53))
        Mgb = n.array(Mgb)
        eMgb = n.array(eMgb)
        
        gmp3367mgb = n.mean(Mgb[gmp3367s])
        gmp3367mgbv = n.sum(eMgb[gmp3367s]**2) / len(Mgb[gmp3367s])**2
        gmp3367fe = n.mean(avFe[gmp3367s])
        gmp3367fev = n.sum(eavFe[gmp3367s]**2) / len(Mgb[gmp3367s])**2
        
        #Loubser2009 data
        galname = []
        Hb = []
        eHb = []
        Mgb = []
        eMgb = []
        Fe52 = []
        eFe52 = []
        Fe53 = []
        eFe53 = []
        NaD = []
        eNaD = []
        with open('/Users/zieleniewski/Data/swift/Coma_stellar_pops/Loubser2009_indices.tsv') as f:
            next(f)
            for line in f:
                data = line.split(';')
                galname.append(data[0].split( )[0])
                Hb.append(float(data[1]))
                eHb.append(float(data[2]))
                Mgb.append(float(data[3]))
                eMgb.append(float(data[4]))
                Fe52.append(float(data[5]))
                eFe52.append(float(data[6]))
                Fe53.append(float(data[7]))
                eFe53.append(float(data[8]))
                NaD.append(float(data[9]))
                eNaD.append(float(data[10]))
        
        gmp2921s = []
        gmp3329s = []
        gmp4928s = []
        gmp3367s = []
        for i in xrange(len(galname)):
            if galname[i] == 'NGC4889':
                gmp2921s.append(i)
            if galname[i] == 'NGC4874':
                gmp3329s.append(i)
            if galname[i] == 'NGC4839':
                gmp4928s.append(i)
            if galname[i] == 'NGC4873':
                gmp3367s.append(i)
        
        print gmp2921s
        
        avFe, eavFe = Fe(n.array(Fe52), n.array(eFe52),n.array(Fe53), n.array(eFe53))
        MgFe, eMgFe = mgfe(n.array(Mgb), n.array(eMgb), n.array(Fe52), n.array(eFe52),n.array(Fe53), n.array(eFe53))
        Mgb = n.array(Mgb)
        eMgb = n.array(eMgb)
        
        pl2921 = P.errorbar(Mgb[gmp2921s], avFe[gmp2921s], xerr=eMgb[gmp2921s], yerr=eavFe[gmp2921s],
        fmt='s', ms=msize+3, mec='r', mew=1.2, mfc='r', ls='', ecolor='r',
        capsize=0, barsabove=False, capthick=1.0, elinewidth=1.4, zorder=15)
        
        pl3329 = P.errorbar(Mgb[gmp3329s], avFe[gmp3329s], xerr=eMgb[gmp3329s], yerr=eavFe[gmp3329s],
        fmt='*', ms=msize+6, mec='b', mew=1.2, mfc='b', ls='', ecolor='b',
        capsize=0, barsabove=False, capthick=1.0, elinewidth=1.4, zorder=15)
        
        pl4928 = P.errorbar(Mgb[gmp4928s], avFe[gmp4928s], xerr=eMgb[gmp4928s], yerr=eavFe[gmp4928s],
        fmt='^', ms=msize+6, mec='g', mew=1.2, mfc='g', ls='', ecolor='g',
        capsize=0, barsabove=False, capthick=1.0, elinewidth=1.4, zorder=15)
        
        pl3367_mean = P.errorbar(gmp3367mgb, gmp3367fe, xerr=n.sqrt(gmp3367mgbv), yerr=n.sqrt(gmp3367fev),
        fmt='D', ms=msize, mec='m', mew=1.2, mfc='m', ls='', ecolor='m',
        capsize=0, barsabove=False, capthick=1.0, elinewidth=1.4, zorder=15)
        
        P.ylabel('<Fe> [$\AA$]', fontsize=tfsize)
        P.annotate('[$\\alpha$/Fe]=+0.2', xy=(cvdafe2ind1[0][3],cvdafe2ind2[0][3]),xycoords='data',xytext=(-125,-20),textcoords='offset points', color='b', fontsize=gensize)
        P.annotate('[$\\alpha$/Fe]=+0.3', xy=(cvdafe3ind1[0][3],cvdafe3ind2[0][3]),xycoords='data',xytext=(-105,-20),textcoords='offset points', color='b', fontsize=gensize)
        
        leg = P.legend([pl2921, pl3329, pl4928, pl3367_mean], ['NGC4889', 'NGC4874', 'NGC4839', 'NGC4873'],
        loc='lower right', frameon=True, borderaxespad=0.8, labelspacing=0.5, numpoints=1, scatterpoints=1)
#####

#####
    if fehsys == 'True':
        #Add min and max systematic range as rectangle
        if ind1=='FeH':
            dat2921 = n.genfromtxt('/Users/zieleniewski/Data/swift/Coma_stellar_pops/GMP2921_FeH_min_max_dat.txt', delimiter=',')
            dat3329 = n.genfromtxt('/Users/zieleniewski/Data/swift/Coma_stellar_pops/GMP3329_FeH_min_max_dat_for_index-index_map.txt', delimiter=',')
            dat4928 = n.genfromtxt('/Users/zieleniewski/Data/swift/Coma_stellar_pops/GMP4928_FeH_min_max_dat.txt', delimiter=',')
            dats = [dat2921, dat3329, dat4928]
            cols = [r_colour, r_colour2, r_colour3]
            inds1 = [ind1vals, ind1vals2, ind1vals3]
            inds2 = [ind2vals, ind2vals2, ind2vals3]
            ginds1 = [gind1vals, gind1vals2, gind1vals3]
            ginds2 = [gind2vals, gind2vals2, gind2vals3]
            
            for i in xrange(len(dats)):
                print 'i = ', i
                print 'ind1vals = ', inds1[i]
                print 'gind1vals = ', ginds1[i]
                print 'ind2vals = ', inds2[i]
                print 'gind2vals = ', ginds2[i]
                alldatsx = n.concatenate((inds1[i],[ginds1[i]]))
                alldatsy = n.concatenate((inds2[i],[ginds2[i]]))
                print 'alldatsx = ', alldatsx
                print 'alldatsy = ', alldatsy
                print 'dats[i][:,1] = ', dats[i][:,1]
                print 'dats[i][:,2] = ', dats[i][:,2]
                print 'cols[i] = ', cols[i]
                makeErrorBoxes(alldatsx, alldatsy, xerror=n.row_stack((alldatsx-dats[i][:,1],dats[i][:,2]-alldatsx)),
                yerror=n.row_stack((n.array([0.0005]*len(dats[i][:,0])),n.array([0.0005]*len(dats[i][:,0])))), fc=n.row_stack((cols[i],cols[i])), ec=n.row_stack((cols[i], cols[i])))
#####

#####
    if fehsys == 'all':
        #Add min and max systematic range as rectangle
        if ind1=='FeH':
            dat2921 = n.genfromtxt('/Users/zieleniewski/Data/swift/Coma_stellar_pops/GMP2921_FeH_min_max_dat.txt', delimiter=',')
            dat3329 = n.genfromtxt('/Users/zieleniewski/Data/swift/Coma_stellar_pops/GMP3329_FeH_min_max_dat_for_index-index_map.txt', delimiter=',')
            dat4928 = n.genfromtxt('/Users/zieleniewski/Data/swift/Coma_stellar_pops/GMP4928_FeH_min_max_dat.txt', delimiter=',')
            dats = [dat2921, dat3329, dat4928]
            cols = [r_colour, r_colour2, r_colour3]
            inds1 = [ind1vals, ind1vals2, ind1vals3]
            inds2 = [ind2vals, ind2vals2, ind2vals3]
            ginds1 = [gind1vals, gind1vals2, gind1vals3]
            ginds2 = [gind2vals, gind2vals2, gind2vals3]
            
            for i in xrange(len(dats)):
                print 'i = ', i
                print 'ind1vals = ', inds1[i]
                print 'gind1vals = ', ginds1[i]
                print 'ind2vals = ', inds2[i]
                print 'gind2vals = ', ginds2[i]
                alldatsx = n.concatenate((inds1[i],[ginds1[i]]))
                alldatsy = n.concatenate((inds2[i],[ginds2[i]]))
                print 'alldatsx = ', alldatsx
                print 'alldatsy = ', alldatsy
                print 'dats[i][:,1] = ', dats[i][:,1]
                print 'dats[i][:,2] = ', dats[i][:,2]
                print 'cols[i] = ', cols[i]
                makeErrorBoxes(alldatsx, alldatsy, xerror=n.row_stack((alldatsx-dats[i][:,1],dats[i][:,2]-alldatsx)),
                yerror=n.row_stack((n.array([0.0005]*len(dats[i][:,0])),n.array([0.0005]*len(dats[i][:,0])))), fc=n.row_stack((cols[i],cols[i])), ec=n.row_stack((cols[i], cols[i])))
#####

#####
    if fehsys == 'global':
        #Add min and max systematic range as rectangle
        if ind1=='FeH':
            dat2921 = n.genfromtxt('/Users/zieleniewski/Data/swift/Coma_stellar_pops/GMP2921_FeH_min_max_dat.txt', delimiter=',')
            dat3329 = n.genfromtxt('/Users/zieleniewski/Data/swift/Coma_stellar_pops/GMP3329_FeH_min_max_dat_for_index-index_map.txt', delimiter=',')
            dat4928 = n.genfromtxt('/Users/zieleniewski/Data/swift/Coma_stellar_pops/GMP4928_FeH_min_max_dat.txt', delimiter=',')
            dats = [dat2921, dat3329, dat4928]
            cols = ['k', 'k', 'k']
            ginds1 = [gind1vals, gind1vals2, gind1vals3]
            ginds2 = [gind2vals, gind2vals2, gind2vals3]
            
            for i in xrange(len(dats)):
                print 'i = ', i
                print 'gind1vals = ', ginds1[i]
                print 'gind2vals = ', ginds2[i]
                alldatsx = n.array([ginds1[i]])
                alldatsy = n.array([ginds2[i]])
                print 'alldatsx = ', alldatsx
                print 'alldatsy = ', alldatsy
                print 'dats[i][-1,1] = ', dats[i][-1,1]
                print 'dats[i][-1,2] = ', dats[i][-1,2]
                print 'cols[i] = ', cols[i]
                makeErrorBoxes(alldatsx, alldatsy, xerror=n.row_stack((alldatsx-dats[i][-1,1],dats[i][-1,2]-alldatsx)),
                yerror=n.row_stack((n.array([0.0005]),n.array([0.0005]))))
#####

#####
    if fehsys == 'radial':
        #Add min and max systematic range as rectangle
        if ind1=='FeH':
            dats = []
            cols = []
            inds1 = []
            inds2 = []
            dat2921 = n.genfromtxt('/Users/zieleniewski/Data/swift/Coma_stellar_pops/GMP2921_FeH_min_max_dat.txt', delimiter=',')
            dat3329 = n.genfromtxt('/Users/zieleniewski/Data/swift/Coma_stellar_pops/GMP3329_FeH_min_max_dat_for_index-index_map.txt', delimiter=',')
            dat4928 = n.genfromtxt('/Users/zieleniewski/Data/swift/Coma_stellar_pops/GMP4928_FeH_min_max_dat.txt', delimiter=',')
            if data != 'None' and data2 =='None' and data3=='None' and data4=='None':
                dats.append(dat2921)
                cols.append(r_colour)
                inds1.append(ind1vals)
                inds2.append(ind2vals)
            if data2 != 'None':
                dats.append(dat3329)
                cols.append(r_colour2)
                inds1.append(ind1vals2)
                inds2.append(ind2vals2)
            if data3 != 'None':
                dats.append(dat4928)
                cols.append(r_colour3)
                inds1.append(ind1vals3)
                inds2.append(ind2vals3)
            # dats = [dat2921, dat3329, dat4928]
            # cols = [r_colour, r_colour2, r_colour3]
            # inds1 = [ind1vals, ind1vals2, ind1vals3]
            # inds2 = [ind2vals, ind2vals2, ind2vals3]
            
            for i in xrange(len(dats)):
                print 'i = ', i
                alldatsx = n.array(inds1[i])
                alldatsy = n.array(inds2[i])
                print 'alldatsx = ', alldatsx
                print 'alldatsy = ', alldatsy
                print 'cols[i] = ', cols[i]
                makeErrorBoxes(alldatsx, alldatsy, xerror=n.row_stack((alldatsx-dats[i][:-1,1],dats[i][:-1,2]-alldatsx)),
                yerror=n.row_stack((n.array([0.0005]*len(dats[i][:-1,0])),n.array([0.0005]*len(dats[i][:-1,0])))),
                fc=n.row_stack((cols[i],cols[i])), ec=n.row_stack((cols[i], cols[i])))
#####

    P.tight_layout()
    P.show()
#------------#


if __name__=="__main__":
    import sys
    print len(sys.argv)
    if len(sys.argv) not in [4,5,6,7,8,9,10,11,12,13,14,15]:
        print ''
        print '1. Index 1'
        print '2. Index 2'
        print '3. Velocity disp [km/s]'
        print '4. data [bool] - add Coma optical data'
        print '5. fehsys [all, global, radial] - add systematic errorbars for FeH'
        print '6. labels [bool] - add SPS labels on plot'
        print '[7. datafile]'
        print '[8. Global datafile]'
        print '[9. datafile2]'
        print '[10. Global datafile2]'
        print '[11. datafile3]'
        print '[12. Global datafile3]'
        print '[13. datafile4]'
        print '[14. Global datafile4]'
        print ''
    elif len(sys.argv) == 4:
        indexindexplot(sys.argv[1], sys.argv[2], sys.argv[3])
    elif len(sys.argv) == 5:
        indexindexplot(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
    elif len(sys.argv) == 6:
        indexindexplot(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
    elif len(sys.argv) == 7:
        indexindexplot(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6])
    elif len(sys.argv) == 8:
        indexindexplot(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7])
    elif len(sys.argv) == 9:
        indexindexplot(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7], sys.argv[8])
    elif len(sys.argv) == 10:
        indexindexplot(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7], sys.argv[8], sys.argv[9])
    elif len(sys.argv) == 11:
        indexindexplot(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7], sys.argv[8], sys.argv[9], sys.argv[10])
    elif len(sys.argv) == 12:
        indexindexplot(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7], sys.argv[8], sys.argv[9], sys.argv[10], sys.argv[11])
    elif len(sys.argv) == 13:
        indexindexplot(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7], sys.argv[8], sys.argv[9], sys.argv[10], sys.argv[11], sys.argv[12])
    elif len(sys.argv) == 14:
        indexindexplot(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7], sys.argv[8], sys.argv[9], sys.argv[10], sys.argv[11], sys.argv[12], sys.argv[13])
    elif len(sys.argv) == 15:
        indexindexplot(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7], sys.argv[8], sys.argv[9], sys.argv[10], sys.argv[11], sys.argv[12], sys.argv[13], sys.argv[14])

#------------#



#------------#
# if len(sys.argv) != 3:
#    print ' '
#    print '1. Index [NaD, NaI, CaT, CaTstar, MgI, TiO, FeH]'
#    print '2. Velocity dispersion [km/s]'
#    print ' '
#    sys.exit()
#
#
# sdss = s.loadSDSSFilters()
# comdisp = float(sys.argv[2])
#
# #Plot index v colour (or index) for models convolved to common veolocity dispersion
# cvd3 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/CvD12/t03.0_solar.ssp', varspecf=None, kins='CvD12', out_sig=comdisp, index='MgI')
# cvd5 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/CvD12/t05.0_solar.ssp', varspecf=None, kins='CvD12', out_sig=comdisp, index='MgI')
# cvd7 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/CvD12/t07.0_solar.ssp', varspecf=None, kins='CvD12', out_sig=comdisp, index='MgI')
# cvd9 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/CvD12/t09.0_solar.ssp', varspecf=None, kins='CvD12', out_sig=comdisp, index='MgI')
# cvd11 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/CvD12/t11.0_solar.ssp', varspecf=None, kins='CvD12', out_sig=comdisp, index='MgI')
# cvd13 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/CvD12/t13.5_solar.ssp', varspecf=None, kins='CvD12', out_sig=comdisp, index='MgI')
# cvdafe2 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/CvD12/t13.5_afe+0.2.ssp', varspecf=None, kins='CvD12', out_sig=comdisp, index='MgI')
# cvdafe3 = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/CvD12/t13.5_afe+0.3.ssp', varspecf=None, kins='CvD12', out_sig=comdisp, index='MgI')
# cvdvar = mri.measure_convolved_spec_index('/Users/zieleniewski/Data/stelpopmods/CvD12/t13.5_varelem.ssp', varspecf=None, kins='CvD12', out_sig=comdisp, index='MgI')
#
# cvd3mgi = cvd3.irindex(cvd3.lam[1]-cvd3.lam[0], index='MgI')
# cvd5mgi = cvd5.irindex(cvd5.lam[1]-cvd5.lam[0], index='MgI')
# cvd7mgi = cvd7.irindex(cvd7.lam[1]-cvd7.lam[0], index='MgI')
# cvd9mgi = cvd9.irindex(cvd9.lam[1]-cvd9.lam[0], index='MgI')
# cvd11mgi = cvd11.irindex(cvd11.lam[1]-cvd11.lam[0], index='MgI')
# cvd13mgi = cvd13.irindex(cvd13.lam[1]-cvd13.lam[0], index='MgI')
# cvdafe2mgi = cvdafe2.irindex(cvdafe2.lam[1]-cvdafe2.lam[0], index='MgI')
# cvdafe3mgi = cvdafe3.irindex(cvdafe3.lam[1]-cvdafe3.lam[0], index='MgI')
# cvdvarmgi = cvdvar.irindex(cvdvar.lam[1]-cvdvar.lam[0], index='MgI')
#
# cvd3tio = cvd3.TiOindex(cvd3.lam[1]-cvd3.lam[0])
# cvd5tio = cvd5.TiOindex(cvd5.lam[1]-cvd5.lam[0])
# cvd7tio = cvd7.TiOindex(cvd7.lam[1]-cvd7.lam[0])
# cvd9tio = cvd9.TiOindex(cvd9.lam[1]-cvd9.lam[0])
# cvd11tio = cvd11.TiOindex(cvd11.lam[1]-cvd11.lam[0])
# cvd13tio = cvd13.TiOindex(cvd13.lam[1]-cvd13.lam[0])
# cvdafe2tio = cvdafe2.TiOindex(cvdafe2.lam[1]-cvdafe2.lam[0])
# cvdafe3tio = cvdafe3.TiOindex(cvdafe3.lam[1]-cvdafe3.lam[0])
# cvdvartio = cvdvar.TiOindex(cvdvar.lam[1]-cvdvar.lam[0])
#
# cvd3nai = cvd3.irindex(cvd3.lam[1]-cvd3.lam[0], index='NaI')
# cvd5nai = cvd5.irindex(cvd5.lam[1]-cvd5.lam[0], index='NaI')
# cvd7nai = cvd7.irindex(cvd7.lam[1]-cvd7.lam[0], index='NaI')
# cvd9nai = cvd9.irindex(cvd9.lam[1]-cvd9.lam[0], index='NaI')
# cvd11nai = cvd11.irindex(cvd11.lam[1]-cvd11.lam[0], index='NaI')
# cvd13nai = cvd13.irindex(cvd13.lam[1]-cvd13.lam[0], index='NaI')
# cvdafe2nai = cvdafe2.irindex(cvdafe2.lam[1]-cvdafe2.lam[0], index='NaI')
# cvdafe3nai = cvdafe3.irindex(cvdafe3.lam[1]-cvdafe3.lam[0], index='NaI')
# cvdvarnai = cvdvar.irindex(cvdvar.lam[1]-cvdvar.lam[0], index='NaI')
#
# cvd3nais = cvd3.irindex(cvd3.lam[1]-cvd3.lam[0], index='NaIsdss')
# cvd5nais = cvd5.irindex(cvd5.lam[1]-cvd5.lam[0], index='NaIsdss')
# cvd7nais = cvd7.irindex(cvd7.lam[1]-cvd7.lam[0], index='NaIsdss')
# cvd9nais = cvd9.irindex(cvd9.lam[1]-cvd9.lam[0], index='NaIsdss')
# cvd11nais = cvd11.irindex(cvd11.lam[1]-cvd11.lam[0], index='NaIsdss')
# cvd13nais = cvd13.irindex(cvd13.lam[1]-cvd13.lam[0], index='NaIsdss')
# cvdafe2nais = cvdafe2.irindex(cvdafe2.lam[1]-cvdafe2.lam[0], index='NaIsdss')
# cvdafe3nais = cvdafe3.irindex(cvdafe3.lam[1]-cvdafe3.lam[0], index='NaIsdss')
# cvdvarnais = cvdvar.irindex(cvdvar.lam[1]-cvdvar.lam[0], index='NaIsdss')
#
# cvd3nad = cvd3.irindex(cvd3.lam[1]-cvd3.lam[0], index='NaD')
# cvd5nad = cvd5.irindex(cvd5.lam[1]-cvd5.lam[0], index='NaD')
# cvd7nad = cvd7.irindex(cvd7.lam[1]-cvd7.lam[0], index='NaD')
# cvd9nad = cvd9.irindex(cvd9.lam[1]-cvd9.lam[0], index='NaD')
# cvd11nad = cvd11.irindex(cvd11.lam[1]-cvd11.lam[0], index='NaD')
# cvd13nad = cvd13.irindex(cvd13.lam[1]-cvd13.lam[0], index='NaD')
# cvdafe2nad = cvdafe2.irindex(cvdafe2.lam[1]-cvdafe2.lam[0], index='NaD')
# cvdafe3nad = cvdafe3.irindex(cvdafe3.lam[1]-cvdafe3.lam[0], index='NaD')
# cvdvarnad = cvdvar.irindex(cvdvar.lam[1]-cvdvar.lam[0], index='NaD')
#
# cvd3ca = cvd3.irindex(cvd3.lam[1]-cvd3.lam[0], index='CaT')
# cvd5ca = cvd5.irindex(cvd5.lam[1]-cvd5.lam[0], index='CaT')
# cvd7ca = cvd7.irindex(cvd7.lam[1]-cvd7.lam[0], index='CaT')
# cvd9ca = cvd9.irindex(cvd9.lam[1]-cvd9.lam[0], index='CaT')
# cvd11ca = cvd11.irindex(cvd11.lam[1]-cvd11.lam[0], index='CaT')
# cvd13ca = cvd13.irindex(cvd13.lam[1]-cvd13.lam[0], index='CaT')
# cvdafe2ca = cvdafe2.irindex(cvdafe2.lam[1]-cvdafe2.lam[0], index='CaT')
# cvdafe3ca = cvdafe3.irindex(cvdafe3.lam[1]-cvdafe3.lam[0], index='CaT')
# cvdvarca = cvdvar.irindex(cvdvar.lam[1]-cvdvar.lam[0], index='CaT')
#
# cvd3cat = cvd3.CaT_star(cvd3.lam[1]-cvd3.lam[0])
# cvd5cat = cvd5.CaT_star(cvd5.lam[1]-cvd5.lam[0])
# cvd7cat = cvd7.CaT_star(cvd7.lam[1]-cvd7.lam[0])
# cvd9cat = cvd9.CaT_star(cvd9.lam[1]-cvd9.lam[0])
# cvd11cat = cvd11.CaT_star(cvd11.lam[1]-cvd11.lam[0])
# cvd13cat = cvd13.CaT_star(cvd13.lam[1]-cvd13.lam[0])
# cvdafe2cat = cvdafe2.CaT_star(cvdafe2.lam[1]-cvdafe2.lam[0])
# cvdafe3cat = cvdafe3.CaT_star(cvdafe3.lam[1]-cvdafe3.lam[0])
# cvdvarcat = cvdvar.CaT_star(cvdvar.lam[1]-cvdvar.lam[0])
#
# cvd3feh = cvd3.irindex(cvd3.lam[1]-cvd3.lam[0], index='FeH')
# cvd5feh = cvd5.irindex(cvd5.lam[1]-cvd5.lam[0], index='FeH')
# cvd7feh = cvd7.irindex(cvd7.lam[1]-cvd7.lam[0], index='FeH')
# cvd9feh = cvd9.irindex(cvd9.lam[1]-cvd9.lam[0], index='FeH')
# cvd11feh = cvd11.irindex(cvd11.lam[1]-cvd11.lam[0], index='FeH')
# cvd13feh = cvd13.irindex(cvd13.lam[1]-cvd13.lam[0], index='FeH')
# cvdafe2feh = cvdafe2.irindex(cvdafe2.lam[1]-cvdafe2.lam[0], index='FeH')
# cvdafe3feh = cvdafe3.irindex(cvdafe3.lam[1]-cvdafe3.lam[0], index='FeH')
# cvdvarfeh = cvdvar.irindex(cvdvar.lam[1]-cvdvar.lam[0], index='FeH')
#
# cvd3c = cvd3.calcABmag(sdss['g'])-cvd3.calcABmag(sdss['r'])
# cvd5c = cvd5.calcABmag(sdss['g'])-cvd5.calcABmag(sdss['r'])
# cvd7c = cvd7.calcABmag(sdss['g'])-cvd7.calcABmag(sdss['r'])
# cvd9c = cvd9.calcABmag(sdss['g'])-cvd9.calcABmag(sdss['r'])
# cvd11c = cvd11.calcABmag(sdss['g'])-cvd11.calcABmag(sdss['r'])
# cvd13c = cvd13.calcABmag(sdss['g'])-cvd13.calcABmag(sdss['r'])
# cvdafe2c = cvdafe2.calcABmag(sdss['g'])-cvdafe2.calcABmag(sdss['r'])
# cvdafe3c = cvdafe3.calcABmag(sdss['g'])-cvdafe3.calcABmag(sdss['r'])
# cvdvarc = cvdvar.calcABmag(sdss['g'])-cvdvar.calcABmag(sdss['r'])
#
# ####FONTS
# tfsize = 28
# lfsize = 24
# majtsize = 16
# mintsize = 10
# twidth = 2
# msize = 10
# ####
#
# #Plot
# P.rcParams['font.family'] = 'Times New Roman'
# P.rcParams.update({'font.size': 20})
# P.figure(figsize=(13,10), dpi=72)
# P.subplots_adjust(left=0.14, bottom=0.12)
#
# if sys.argv[1] == 'NaD':
#    ###NaI v NaD
#    #Changing IMF at 13.5 Gyr
#    P.scatter(cvd13nai[0][1:], cvd13nad[0][1:], marker='o',c='k',s=[90,70,50,30])
#    #Changing age for x=3 IMF
#    P.scatter([cvd3nai[0][1],cvd5nai[0][1],cvd7nai[0][1],cvd9nai[0][1],cvd11nai[0][1],cvd13nai[0][1]],
#              [cvd3nad[0][1],cvd5nad[0][1],cvd7nad[0][1],cvd9nad[0][1],cvd11nad[0][1],cvd13nad[0][1]],
#              marker='s',c='r', s=[110,90,70,50,30,10])
#    #Changing age for Chabrier IMF
#    P.scatter([cvd3nai[0][3],cvd5nai[0][3],cvd7nai[0][3],cvd9nai[0][3],cvd11nai[0][3],cvd13nai[0][3]],
#              [cvd3nad[0][3],cvd5nad[0][3],cvd7nad[0][3],cvd9nad[0][3],cvd11nad[0][3],cvd13nad[0][3]],
#              marker='s',c='r', s=[110,90,70,50,30,10])
#    #Changing [alpha/Fe] for 13.5 Gyr Chabrier IMF
#    P.scatter([cvd13nai[0][3],cvdafe2nai[0][3],cvdafe3nai[0][3]],
#              [cvd13nad[0][3],cvdafe2nad[0][3],cvdafe3nad[0][3]],
#              marker='D',c='b',s=[10,30,50])
#    #Varying [Na/Fe] for 13.5 Gyr Chabrier IMF
#    P.scatter([cvdvarnai[0][0],cvd13nai[0][3],cvdvarnai[0][1]],
#              [cvdvarnad[0][0],cvd13nad[0][3],cvdvarnad[0][1]],
#              marker='D',c='g',s=[50,5,50])
#    #Lines for each variation above
#    P.plot(cvd13nai[0][1:], cvd13nad[0][1:], 'k-',
#           [cvd3nai[0][1],cvd5nai[0][1],cvd7nai[0][1],cvd9nai[0][1],cvd11nai[0][1],cvd13nai[0][1]],
#           [cvd3nad[0][1],cvd5nad[0][1],cvd7nad[0][1],cvd9nad[0][1],cvd11nad[0][1],cvd13nad[0][1]],'r--',
#           [cvd3nai[0][3],cvd5nai[0][3],cvd7nai[0][3],cvd9nai[0][3],cvd11nai[0][3],cvd13nai[0][3]],
#           [cvd3nad[0][3],cvd5nad[0][3],cvd7nad[0][3],cvd9nad[0][3],cvd11nad[0][3],cvd13nad[0][3]],'r--',
#           [cvd13nai[0][3],cvdafe2nai[0][3],cvdafe3nai[0][3]],[cvd13nad[0][3],cvdafe2nad[0][3],cvdafe3nad[0][3]],'b-',
#           [cvdvarnai[0][0],cvd13nai[0][3],cvdvarnai[0][1]],[cvdvarnad[0][0],cvd13nad[0][3],cvdvarnad[0][1]],'g--')
#    P.xlabel('NaI0.82 [$\AA$]', fontsize=tfsize)
#    P.ylabel('NaD0.59 [$\AA$]', fontsize=tfsize)
#    P.annotate('13', xy=(cvd13nai[0][1],cvd13nad[0][1]),xycoords='data',xytext=(10,-10),textcoords='offset points', color='r', fontsize=lfsize)
#    P.annotate('11', xy=(cvd11nai[0][1],cvd11nad[0][1]),xycoords='data',xytext=(5,-10),textcoords='offset points', color='r', fontsize=lfsize)
#    P.annotate('9', xy=(cvd9nai[0][1],cvd9nad[0][1]),xycoords='data',xytext=(10,-15),textcoords='offset points', color='r', fontsize=lfsize)
#    P.annotate('7', xy=(cvd7nai[0][1],cvd7nad[0][1]),xycoords='data',xytext=(10,-15),textcoords='offset points', color='r', fontsize=lfsize)
#    P.annotate('5', xy=(cvd5nai[0][1],cvd5nad[0][1]),xycoords='data',xytext=(15,-15),textcoords='offset points', color='r', fontsize=lfsize)
#    P.annotate('3 Gyr', xy=(cvd3nai[0][1],cvd3nad[0][1]),xycoords='data',xytext=(10,-10),textcoords='offset points', color='r', fontsize=lfsize)
#    #P.annotate('x = 3.5', xy=(cvd13nai[0][0],cvd13nad[0][0]),xycoords='data',xytext=(-20,10),textcoords='offset points', color='k', fontsize=lfsize)
#    P.annotate('x = 3', xy=(cvd13nai[0][1],cvd13nad[0][1]),xycoords='data',xytext=(-10,10),textcoords='offset points', color='k', fontsize=lfsize)
#    P.annotate('x = 2.35', xy=(cvd13nai[0][2],cvd13nad[0][2]),xycoords='data',xytext=(-25,10),textcoords='offset points', color='k', fontsize=lfsize)
#    P.annotate('Chab', xy=(cvd13nai[0][3],cvd13nad[0][3]),xycoords='data',xytext=(-50,10),textcoords='offset points', color='k', fontsize=lfsize)
#    P.annotate('B-L', xy=(cvd13nai[0][4],cvd13nad[0][4]),xycoords='data',xytext=(-20,10),textcoords='offset points', color='k', fontsize=lfsize)
#    P.annotate('[$\\alpha$/Fe]=+0.2', xy=(cvdafe2nai[0][3],cvdafe2nad[0][3]),xycoords='data',xytext=(-130,-5),textcoords='offset points', color='b', fontsize=lfsize)
#    P.annotate('[$\\alpha$/Fe]=+0.3', xy=(cvdafe3nai[0][3],cvdafe3nad[0][3]),xycoords='data',xytext=(-65,-35),textcoords='offset points', color='b', fontsize=lfsize)
#    P.annotate('[Na/Fe]=+0.3', xy=(cvdvarnai[0][0],cvdvarnad[0][0]),xycoords='data',xytext=(-30,10),textcoords='offset points', color='g', fontsize=lfsize)
#    P.annotate('[Na/Fe]=-0.3', xy=(cvdvarnai[0][1],cvdvarnad[0][1]),xycoords='data',xytext=(-40,-25),textcoords='offset points', color='g', fontsize=lfsize)
#
# elif sys.argv[1] == 'NaI':
#    ###NaI v g-r
#    P.scatter(cvd13nai[0][1:], cvd13c[1:], marker='o',c='k',s=[90,70,50,30])
#    P.scatter([cvd3nai[0][1],cvd5nai[0][1],cvd7nai[0][1],cvd9nai[0][1],cvd11nai[0][1],cvd13nai[0][1]],
#              [cvd3c[1],cvd5c[1],cvd7c[1],cvd9c[1],cvd11c[1],cvd13c[1]],
#              marker='s',c='r', s=[110,90,70,50,30,10])
#    P.scatter([cvd3nai[0][3],cvd5nai[0][3],cvd7nai[0][3],cvd9nai[0][3],cvd11nai[0][3],cvd13nai[0][3]],
#              [cvd3c[3],cvd5c[3],cvd7c[3],cvd9c[3],cvd11c[3],cvd13c[3]],
#              marker='s',c='r', s=[110,90,70,50,30,10])
#    P.scatter([cvd13nai[0][3],cvdafe2nai[0][3],cvdafe3nai[0][3]],
#              [cvd13c[3],cvdafe2c[3],cvdafe3c[3]],
#              marker='D',c='b',s=[10,30,50])
#    #Varying [Na/Fe] abundance
#    P.scatter([cvdvarnai[0][0],cvd13nai[0][3],cvdvarnai[0][1]],
#              [cvdvarc[0],cvd13c[3],cvdvarc[1]],
#              marker='D',c='g',s=[50,5,50])
#    P.plot(cvd13nai[0][1:], cvd13c[1:], 'k-',
#           [cvd3nai[0][1],cvd5nai[0][1],cvd7nai[0][1],cvd9nai[0][1],cvd11nai[0][1],cvd13nai[0][1]],
#           [cvd3c[1],cvd5c[1],cvd7c[1],cvd9c[1],cvd11c[1],cvd13c[1]],'r--',
#           [cvd3nai[0][3],cvd5nai[0][3],cvd7nai[0][3],cvd9nai[0][3],cvd11nai[0][3],cvd13nai[0][3]],
#           [cvd3c[3],cvd5c[3],cvd7c[3],cvd9c[3],cvd11c[3],cvd13c[3]],'r--',
#           [cvd13nai[0][3],cvdafe2nai[0][3],cvdafe3nai[0][3]],[cvd13c[3],cvdafe2c[3],cvdafe3c[3]],'b-',
#           [cvdvarnai[0][0],cvd13nai[0][3],cvdvarnai[0][1]],[cvdvarc[0],cvd13c[3],cvdvarc[1]],'g--')
#    P.xlabel('NaI0.82 [$\AA$]', fontsize=tfsize)
#    P.ylabel('$g - r$', fontsize=tfsize)
#    P.annotate('13', xy=(cvd13nai[0][1],cvd13c[1]),xycoords='data',xytext=(5,-10),textcoords='offset points', color='r', fontsize=lfsize)
#    P.annotate('11', xy=(cvd11nai[0][1],cvd11c[1]),xycoords='data',xytext=(5,0),textcoords='offset points', color='r', fontsize=lfsize)
#    P.annotate('9', xy=(cvd9nai[0][1],cvd9c[1]),xycoords='data',xytext=(10,-10),textcoords='offset points', color='r', fontsize=lfsize)
#    P.annotate('7', xy=(cvd7nai[0][1],cvd7c[1]),xycoords='data',xytext=(10,-5),textcoords='offset points', color='r', fontsize=lfsize)
#    P.annotate('5', xy=(cvd5nai[0][1],cvd5c[1]),xycoords='data',xytext=(15,-10),textcoords='offset points', color='r', fontsize=lfsize)
#    P.annotate('3 Gyr', xy=(cvd3nai[0][1],cvd3c[1]),xycoords='data',xytext=(10,0),textcoords='offset points', color='r', fontsize=lfsize)
#    #P.annotate('x = 3.5', xy=(cvd13nai[0][0],cvd13c[0]),xycoords='data',xytext=(-5,10),textcoords='offset points', color='k', fontsize=lfsize)
#    P.annotate('x = 3', xy=(cvd13nai[0][1],cvd13c[1]),xycoords='data',xytext=(0,10),textcoords='offset points', color='k', fontsize=lfsize)
#    P.annotate('x = 2.35', xy=(cvd13nai[0][2],cvd13c[2]),xycoords='data',xytext=(-25,12),textcoords='offset points', color='k', fontsize=lfsize)
#    P.annotate('Chab', xy=(cvd13nai[0][3],cvd13c[3]),xycoords='data',xytext=(-20,10),textcoords='offset points', color='k', fontsize=lfsize)
#    P.annotate('B-L', xy=(cvd13nai[0][4],cvd13c[4]),xycoords='data',xytext=(-20,10),textcoords='offset points', color='k', fontsize=lfsize)
#    P.annotate('[$\\alpha$/Fe]=+0.2', xy=(cvdafe2nai[0][3],cvdafe2c[3]),xycoords='data',xytext=(-100,0),textcoords='offset points', color='b', fontsize=lfsize)
#    P.annotate('[$\\alpha$/Fe]=+0.3', xy=(cvdafe3nai[0][3],cvdafe3c[3]),xycoords='data',xytext=(-40,-35),textcoords='offset points', color='b', fontsize=lfsize)
#    P.annotate('[Na/Fe]=+0.3', xy=(cvdvarnai[0][0],cvdvarc[0]),xycoords='data',xytext=(-30,-20),textcoords='offset points', color='g', fontsize=lfsize)
#    P.annotate('[Na/Fe]=-0.3', xy=(cvdvarnai[0][1],cvdvarc[1]),xycoords='data',xytext=(-140,10),textcoords='offset points', color='g', fontsize=lfsize)
#
#    #####05-01-15: Plot M31 g-r and NaI datapoints!
#    m31coldat = n.genfromtxt('/Users/zieleniewski/Documents/Oxford/D.Phil/stelpops/Index_measurements/M31_sdss_colour_data/M31_SDSS_g-r_colour_data_Saglia2010_22-12-14.txt', delimiter=',')
#    #m31coldattempel = n.genfromtxt('/Users/zieleniewski/Documents/Oxford/D.Phil/stelpops/Index_measurements/M31_sdss_colour_data/M31_SDSS_colour_data_Tempel2011_22-12-14.txt', delimiter=',')
#
#    m31dat = n.genfromtxt('/Users/zieleniewski/Documents/Oxford/D.Phil/stelpops/Index_measurements/17-12-14/NaI_index_values_200kms_17-12-14.txt', delimiter=',')
#    m31moddat = n.genfromtxt('/Users/zieleniewski/Documents/Oxford/D.Phil/stelpops/Index_measurements/17-12-14/NaI_index_values_bestfitMODEL_200kms_17-12-14.txt', delimiter=',')
#    m31_1dat = n.genfromtxt('/Users/zieleniewski/Documents/Oxford/D.Phil/stelpops/Index_measurements/m31/NaI_index_values_mids_m31_1_21-12-14.txt', delimiter=',')
#    m31_1mdat = n.genfromtxt('/Users/zieleniewski/Documents/Oxford/D.Phil/stelpops/Index_measurements/m31/NaI_index_values_mids_bestfitMODEL_m31_1_21-12-14.txt', delimiter=',')
#
# #    m31dat = n.genfromtxt('/Users/zieleniewski/Documents/Oxford/D.Phil/stelpops/Index_measurements/20-04-15/NaI_index_values_200kms_20-04-15_newdefB.txt', delimiter=',')
# #    m31moddat = n.genfromtxt('/Users/zieleniewski/Documents/Oxford/D.Phil/stelpops/Index_measurements/20-04-15/NaI_index_values_bestfitMODEL_200kms_20-04-15_newdefB.txt', delimiter=',')
# #    m31_1dat = n.genfromtxt('/Users/zieleniewski/Documents/Oxford/D.Phil/stelpops/Index_measurements/20-04-15/NaI_index_values_mids_m31_1_20-04-15_newdefB.txt', delimiter=',')
# #    m31_1mdat = n.genfromtxt('/Users/zieleniewski/Documents/Oxford/D.Phil/stelpops/Index_measurements/20-04-15/NaI_index_values_mids_bestfitMODEL_m31_1_20-04-15_newdefB.txt', delimiter=',')
#
#
#    m31col = n.zeros(11, dtype=float)
#    m31col[:] = m31coldat[:,1]
#    #m31col[-1] = m31coldattempel[-1,1]
#    m31feh = n.zeros(11, dtype=float); m31feherrs = n.zeros(11, dtype=float)
#    m31feh[:6] = m31_1dat[:,1]; m31feherrs[:6] = m31_1dat[:,2]
#    m31feh[6:] = m31dat[:5,1]; m31feherrs[6:] = m31dat[:5,2]
#    m31feh[-1] = m31moddat[4,1]; m31feherrs[-1] = m31moddat[4,2]
#
#    rad_dat = n.log10(m31coldat[:,0])
#
#    plot = P.scatter(m31feh, m31col, marker='o', s=[220,200,180,160,140,120,100,80,60,40], c=rad_dat, linewidths=0, cmap=P.cm.copper_r, zorder=10)
#    clb = P.colorbar(plot)
#    clb.set_label('log$_{10}$(r/arcsec)', fontsize=tfsize)
#    #convert time to a color tuple using the colormap used for scatter
#    rad_color = clb.to_rgba(rad_dat)
#    #pm0.05 RANDOM ERROR ADDED IN QUADRATURE - 28-05-15
#    m31feherrs[-5:] = n.sqrt(m31feherrs[-5:]**2+0.05**2)
#    a,b,c = P.errorbar(m31feh, m31col, xerr=m31feherrs, #yerr=m31coldat[:,2],
#                       marker='', ls='', capsize=0, barsabove=False, capthick=1.0,
#                       elinewidth=1.8, zorder=10)
#
#    #adjust the color of c[0], which is a LineCollection, to the colormap
#    c[0].set_color(rad_color)
#    #c[1].set_color(rad_color)
#
#    P.xlim([0.0, 0.8])
#    P.ylim([0.65, 0.9])
#    #####End
#
# elif sys.argv[1] == 'NaIsdss':
#    ###NaI v g-r
#    P.scatter(cvd13nais[0][1:], cvd13c[1:], marker='o',c='k',s=[90,70,50,30])
#    P.scatter([cvd3nais[0][1],cvd5nais[0][1],cvd7nais[0][1],cvd9nais[0][1],cvd11nais[0][1],cvd13nais[0][1]],
#              [cvd3c[1],cvd5c[1],cvd7c[1],cvd9c[1],cvd11c[1],cvd13c[1]],
#              marker='s',c='r', s=[110,90,70,50,30,10])
#    P.scatter([cvd3nais[0][3],cvd5nais[0][3],cvd7nais[0][3],cvd9nais[0][3],cvd11nais[0][3],cvd13nais[0][3]],
#              [cvd3c[3],cvd5c[3],cvd7c[3],cvd9c[3],cvd11c[3],cvd13c[3]],
#              marker='s',c='r', s=[110,90,70,50,30,10])
#    P.scatter([cvd13nais[0][3],cvdafe2nais[0][3],cvdafe3nais[0][3]],
#              [cvd13c[3],cvdafe2c[3],cvdafe3c[3]],
#              marker='D',c='b',s=[10,30,50])
#    #Varying [Na/Fe] abundance
#    P.scatter([cvdvarnais[0][0],cvd13nais[0][3],cvdvarnais[0][1]],
#              [cvdvarc[0],cvd13c[3],cvdvarc[1]],
#              marker='D',c='g',s=[50,5,50])
#    P.plot(cvd13nais[0][1:], cvd13c[1:], 'k-',
#           [cvd3nais[0][1],cvd5nais[0][1],cvd7nais[0][1],cvd9nais[0][1],cvd11nais[0][1],cvd13nais[0][1]],
#           [cvd3c[1],cvd5c[1],cvd7c[1],cvd9c[1],cvd11c[1],cvd13c[1]],'r--',
#           [cvd3nais[0][3],cvd5nais[0][3],cvd7nais[0][3],cvd9nais[0][3],cvd11nais[0][3],cvd13nais[0][3]],
#           [cvd3c[3],cvd5c[3],cvd7c[3],cvd9c[3],cvd11c[3],cvd13c[3]],'r--',
#           [cvd13nais[0][3],cvdafe2nais[0][3],cvdafe3nais[0][3]],[cvd13c[3],cvdafe2c[3],cvdafe3c[3]],'b-',
#           [cvdvarnais[0][0],cvd13nais[0][3],cvdvarnais[0][1]],[cvdvarc[0],cvd13c[3],cvdvarc[1]],'g--')
#    P.xlabel('NaI$_{SDSS}$ [$\AA$]', fontsize=tfsize)
#    P.ylabel('$g - r$', fontsize=tfsize)
#    P.annotate('13', xy=(cvd13nais[0][1],cvd13c[1]),xycoords='data',xytext=(5,-10),textcoords='offset points', color='r', fontsize=lfsize)
#    P.annotate('11', xy=(cvd11nais[0][1],cvd11c[1]),xycoords='data',xytext=(5,0),textcoords='offset points', color='r', fontsize=lfsize)
#    P.annotate('9', xy=(cvd9nais[0][1],cvd9c[1]),xycoords='data',xytext=(10,-10),textcoords='offset points', color='r', fontsize=lfsize)
#    P.annotate('7', xy=(cvd7nais[0][1],cvd7c[1]),xycoords='data',xytext=(10,-5),textcoords='offset points', color='r', fontsize=lfsize)
#    P.annotate('5', xy=(cvd5nais[0][1],cvd5c[1]),xycoords='data',xytext=(15,-10),textcoords='offset points', color='r', fontsize=lfsize)
#    P.annotate('3 Gyr', xy=(cvd3nais[0][1],cvd3c[1]),xycoords='data',xytext=(10,0),textcoords='offset points', color='r', fontsize=lfsize)
#    #P.annotate('x = 3.5', xy=(cvd13nais[0][0],cvd13c[0]),xycoords='data',xytext=(-5,10),textcoords='offset points', color='k', fontsize=lfsize)
#    P.annotate('x = 3', xy=(cvd13nais[0][1],cvd13c[1]),xycoords='data',xytext=(0,10),textcoords='offset points', color='k', fontsize=lfsize)
#    P.annotate('x = 2.35', xy=(cvd13nais[0][2],cvd13c[2]),xycoords='data',xytext=(-25,12),textcoords='offset points', color='k', fontsize=lfsize)
#    P.annotate('Chab', xy=(cvd13nais[0][3],cvd13c[3]),xycoords='data',xytext=(-20,10),textcoords='offset points', color='k', fontsize=lfsize)
#    P.annotate('B-L', xy=(cvd13nais[0][4],cvd13c[4]),xycoords='data',xytext=(-20,10),textcoords='offset points', color='k', fontsize=lfsize)
#    P.annotate('[$\\alpha$/Fe]=+0.2', xy=(cvdafe2nais[0][3],cvdafe2c[3]),xycoords='data',xytext=(-100,0),textcoords='offset points', color='b', fontsize=lfsize)
#    P.annotate('[$\\alpha$/Fe]=+0.3', xy=(cvdafe3nais[0][3],cvdafe3c[3]),xycoords='data',xytext=(-40,-35),textcoords='offset points', color='b', fontsize=lfsize)
#    P.annotate('[Na/Fe]=+0.3', xy=(cvdvarnais[0][0],cvdvarc[0]),xycoords='data',xytext=(-30,-20),textcoords='offset points', color='g', fontsize=lfsize)
#    P.annotate('[Na/Fe]=-0.3', xy=(cvdvarnais[0][1],cvdvarc[1]),xycoords='data',xytext=(-140,10),textcoords='offset points', color='g', fontsize=lfsize)
#
# elif sys.argv[1] == 'CaT':
#    ###CaT v g-r
#    P.scatter(cvd13ca[0][1:], cvd13c[1:], marker='o',c='k',s=[90,70,50,30])
#    P.scatter([cvd3ca[0][1],cvd5ca[0][1],cvd7ca[0][1],cvd9ca[0][1],cvd11ca[0][1],cvd13ca[0][1]],
#              [cvd3c[1],cvd5c[1],cvd7c[1],cvd9c[1],cvd11c[1],cvd13c[1]],
#              marker='s',c='r', s=[110,90,70,50,30,10])
#    P.scatter([cvd3ca[0][3],cvd5ca[0][3],cvd7ca[0][3],cvd9ca[0][3],cvd11ca[0][3],cvd13ca[0][3]],
#              [cvd3c[3],cvd5c[3],cvd7c[3],cvd9c[3],cvd11c[3],cvd13c[3]],
#              marker='s',c='r', s=[110,90,70,50,30,10])
#    P.scatter([cvd13ca[0][3],cvdafe2ca[0][3],cvdafe3ca[0][3]],
#              [cvd13c[3],cvdafe2c[3],cvdafe3c[3]],
#              marker='D',c='b',s=[10,30,50])
#    #Varying [Ca/Fe] abundance
#    P.scatter([cvdvarca[0][2],cvd13ca[0][3],cvdvarca[0][3]],
#              [cvdvarc[2],cvd13c[3],cvdvarc[3]],
#              marker='D',c='g',s=[50,5,50])
#    P.plot(cvd13ca[0][1:], cvd13c[1:], 'k-',
#           [cvd3ca[0][1],cvd5ca[0][1],cvd7ca[0][1],cvd9ca[0][1],cvd11ca[0][1],cvd13ca[0][1]],
#           [cvd3c[1],cvd5c[1],cvd7c[1],cvd9c[1],cvd11c[1],cvd13c[1]],'r--',
#           [cvd3ca[0][3],cvd5ca[0][3],cvd7ca[0][3],cvd9ca[0][3],cvd11ca[0][3],cvd13ca[0][3]],
#           [cvd3c[3],cvd5c[3],cvd7c[3],cvd9c[3],cvd11c[3],cvd13c[3]],'r--',
#           [cvd13ca[0][3],cvdafe2ca[0][3],cvdafe3ca[0][3]],[cvd13c[3],cvdafe2c[3],cvdafe3c[3]],'b-',
#           [cvdvarca[0][2],cvd13ca[0][3],cvdvarca[0][3]],[cvdvarc[2],cvd13c[3],cvdvarc[3]],'g--')
#    P.xlabel('CaT [$\AA$]', fontsize=tfsize)
#    P.ylabel('$g - r$', fontsize=tfsize)
#    P.annotate('13', xy=(cvd13ca[0][1],cvd13c[1]),xycoords='data',xytext=(10,-25),textcoords='offset points', color='r', fontsize=lfsize)
#    P.annotate('11', xy=(cvd11ca[0][1],cvd11c[1]),xycoords='data',xytext=(5,5),textcoords='offset points', color='r', fontsize=lfsize)
#    P.annotate('9', xy=(cvd9ca[0][1],cvd9c[1]),xycoords='data',xytext=(10,0),textcoords='offset points', color='r', fontsize=lfsize)
#    P.annotate('7', xy=(cvd7ca[0][1],cvd7c[1]),xycoords='data',xytext=(10,5),textcoords='offset points', color='r', fontsize=lfsize)
#    P.annotate('5', xy=(cvd5ca[0][1],cvd5c[1]),xycoords='data',xytext=(15,0),textcoords='offset points', color='r', fontsize=lfsize)
#    P.annotate('3 Gyr', xy=(cvd3ca[0][1],cvd3c[1]),xycoords='data',xytext=(10,0),textcoords='offset points', color='r', fontsize=lfsize)
#    #P.annotate('x = 3.5', xy=(cvd13ca[0][0],cvd13c[0]),xycoords='data',xytext=(-5,10),textcoords='offset points', color='k', fontsize=lfsize)
#    P.annotate('x = 3', xy=(cvd13ca[0][1],cvd13c[1]),xycoords='data',xytext=(0,10),textcoords='offset points', color='k', fontsize=lfsize)
#    P.annotate('x = 2.35', xy=(cvd13ca[0][2],cvd13c[2]),xycoords='data',xytext=(-30,20),textcoords='offset points', color='k', fontsize=lfsize)
#    P.annotate('Chab', xy=(cvd13ca[0][3],cvd13c[3]),xycoords='data',xytext=(-15,10),textcoords='offset points', color='k', fontsize=lfsize)
#    P.annotate('B-L', xy=(cvd13ca[0][4],cvd13c[4]),xycoords='data',xytext=(15,5),textcoords='offset points', color='k', fontsize=lfsize)
#    P.annotate('[$\\alpha$/Fe]=+0.2', xy=(cvdafe2ca[0][3],cvdafe2c[3]),xycoords='data',xytext=(-50,5),textcoords='offset points', color='b', fontsize=lfsize)
#    P.annotate('[$\\alpha$/Fe]=+0.3', xy=(cvdafe3ca[0][3],cvdafe3c[3]),xycoords='data',xytext=(-40,2),textcoords='offset points', color='b', fontsize=lfsize)
#    P.annotate('[Ca/Fe]=+0.15', xy=(cvdvarca[0][2],cvdvarc[2]),xycoords='data',xytext=(-50,10),textcoords='offset points', color='g', fontsize=lfsize)
#    P.annotate('[Ca/Fe]=-0.15', xy=(cvdvarca[0][3],cvdvarc[3]),xycoords='data',xytext=(-100,10),textcoords='offset points', color='g', fontsize=lfsize)
#
#
# elif sys.argv[1] == 'CaTstar':
#    ###CaT v g-r
#    P.scatter(cvd13cat[0][1:], cvd13c[1:], marker='o',c='k',s=[90,70,50,30])
#    P.scatter([cvd3cat[0][1],cvd5cat[0][1],cvd7cat[0][1],cvd9cat[0][1],cvd11cat[0][1],cvd13cat[0][1]],
#              [cvd3c[1],cvd5c[1],cvd7c[1],cvd9c[1],cvd11c[1],cvd13c[1]],
#              marker='s',c='r', s=[110,90,70,50,30,10])
#    P.scatter([cvd3cat[0][3],cvd5cat[0][3],cvd7cat[0][3],cvd9cat[0][3],cvd11cat[0][3],cvd13cat[0][3]],
#              [cvd3c[3],cvd5c[3],cvd7c[3],cvd9c[3],cvd11c[3],cvd13c[3]],
#              marker='s',c='r', s=[110,90,70,50,30,10])
#    P.scatter([cvd13cat[0][3],cvdafe2cat[0][3],cvdafe3cat[0][3]],
#              [cvd13c[3],cvdafe2c[3],cvdafe3c[3]],
#              marker='D',c='b',s=[10,30,50])
#    #Varying [Ca/Fe] abundance
#    P.scatter([cvdvarcat[0][2],cvd13cat[0][3],cvdvarcat[0][3]],
#              [cvdvarc[2],cvd13c[3],cvdvarc[3]],
#              marker='D',c='g',s=[50,5,50])
#    P.plot(cvd13cat[0][1:], cvd13c[1:], 'k-',
#           [cvd3cat[0][1],cvd5cat[0][1],cvd7cat[0][1],cvd9cat[0][1],cvd11cat[0][1],cvd13cat[0][1]],
#           [cvd3c[1],cvd5c[1],cvd7c[1],cvd9c[1],cvd11c[1],cvd13c[1]],'r--',
#           [cvd3cat[0][3],cvd5cat[0][3],cvd7cat[0][3],cvd9cat[0][3],cvd11cat[0][3],cvd13cat[0][3]],
#           [cvd3c[3],cvd5c[3],cvd7c[3],cvd9c[3],cvd11c[3],cvd13c[3]],'r--',
#           [cvd13cat[0][3],cvdafe2cat[0][3],cvdafe3cat[0][3]],[cvd13c[3],cvdafe2c[3],cvdafe3c[3]],'b-',
#           [cvdvarcat[0][2],cvd13cat[0][3],cvdvarcat[0][3]],[cvdvarc[2],cvd13c[3],cvdvarc[3]],'g--')
#    P.xlabel('CaT* [$\AA$]', fontsize=tfsize)
#    P.ylabel('$g - r$', fontsize=tfsize)
#    P.annotate('13', xy=(cvd13cat[0][1],cvd13c[1]),xycoords='data',xytext=(10,-25),textcoords='offset points', color='r', fontsize=lfsize)
#    P.annotate('11', xy=(cvd11cat[0][1],cvd11c[1]),xycoords='data',xytext=(5,5),textcoords='offset points', color='r', fontsize=lfsize)
#    P.annotate('9', xy=(cvd9cat[0][1],cvd9c[1]),xycoords='data',xytext=(10,0),textcoords='offset points', color='r', fontsize=lfsize)
#    P.annotate('7', xy=(cvd7cat[0][1],cvd7c[1]),xycoords='data',xytext=(10,5),textcoords='offset points', color='r', fontsize=lfsize)
#    P.annotate('5', xy=(cvd5cat[0][1],cvd5c[1]),xycoords='data',xytext=(15,0),textcoords='offset points', color='r', fontsize=lfsize)
#    P.annotate('3 Gyr', xy=(cvd3cat[0][1],cvd3c[1]),xycoords='data',xytext=(10,0),textcoords='offset points', color='r', fontsize=lfsize)
#    #P.annotate('x = 3.5', xy=(cvd13cat[0][0],cvd13c[0]),xycoords='data',xytext=(-5,10),textcoords='offset points', color='k', fontsize=lfsize)
#    P.annotate('x = 3', xy=(cvd13cat[0][1],cvd13c[1]),xycoords='data',xytext=(0,10),textcoords='offset points', color='k', fontsize=lfsize)
#    P.annotate('x = 2.35', xy=(cvd13cat[0][2],cvd13c[2]),xycoords='data',xytext=(-30,20),textcoords='offset points', color='k', fontsize=lfsize)
#    P.annotate('Chab', xy=(cvd13cat[0][3],cvd13c[3]),xycoords='data',xytext=(-15,10),textcoords='offset points', color='k', fontsize=lfsize)
#    P.annotate('B-L', xy=(cvd13cat[0][4],cvd13c[4]),xycoords='data',xytext=(15,5),textcoords='offset points', color='k', fontsize=lfsize)
#    P.annotate('[$\\alpha$/Fe]=+0.2', xy=(cvdafe2cat[0][3],cvdafe2c[3]),xycoords='data',xytext=(-50,5),textcoords='offset points', color='b', fontsize=lfsize)
#    P.annotate('[$\\alpha$/Fe]=+0.3', xy=(cvdafe3cat[0][3],cvdafe3c[3]),xycoords='data',xytext=(-40,2),textcoords='offset points', color='b', fontsize=lfsize)
#    P.annotate('[Ca/Fe]=+0.15', xy=(cvdvarcat[0][2],cvdvarc[2]),xycoords='data',xytext=(-50,10),textcoords='offset points', color='g', fontsize=lfsize)
#    P.annotate('[Ca/Fe]=-0.15', xy=(cvdvarcat[0][3],cvdvarc[3]),xycoords='data',xytext=(-35,10),textcoords='offset points', color='g', fontsize=lfsize)
#
# elif sys.argv[1] == 'MgI':
#    ###MgI v g-r
#    P.scatter(cvd13mgi[0][1:], cvd13c[1:], marker='o',c='k',s=[90,70,50,30])
#    P.scatter([cvd3mgi[0][1],cvd5mgi[0][1],cvd7mgi[0][1],cvd9mgi[0][1],cvd11mgi[0][1],cvd13mgi[0][1]],
#              [cvd3c[1],cvd5c[1],cvd7c[1],cvd9c[1],cvd11c[1],cvd13c[1]],
#              marker='s',c='r', s=[110,90,70,50,30,10])
#    P.scatter([cvd3mgi[0][3],cvd5mgi[0][3],cvd7mgi[0][3],cvd9mgi[0][3],cvd11mgi[0][3],cvd13mgi[0][3]],
#              [cvd3c[3],cvd5c[3],cvd7c[3],cvd9c[3],cvd11c[3],cvd13c[3]],
#              marker='s',c='r', s=[110,90,70,50,30,10])
#    P.scatter([cvd13mgi[0][3],cvdafe2mgi[0][3],cvdafe3mgi[0][3]],
#              [cvd13c[3],cvdafe2c[3],cvdafe3c[3]],
#              marker='D',c='b',s=[10,30,50])
#    #Plot [Mg/Fe] abundance
#    P.scatter([cvdvarmgi[0][14],cvd13mgi[0][3],cvdvarmgi[0][15]],
#              [cvdvarc[14],cvd13c[3],cvdvarc[15]],
#              marker='D',c='g',s=[50,5,50])
#    P.plot(cvd13mgi[0][1:], cvd13c[1:], 'k-',
#           [cvd3mgi[0][1],cvd5mgi[0][1],cvd7mgi[0][1],cvd9mgi[0][1],cvd11mgi[0][1],cvd13mgi[0][1]],
#           [cvd3c[1],cvd5c[1],cvd7c[1],cvd9c[1],cvd11c[1],cvd13c[1]],'r--',
#           [cvd3mgi[0][3],cvd5mgi[0][3],cvd7mgi[0][3],cvd9mgi[0][3],cvd11mgi[0][3],cvd13mgi[0][3]],
#           [cvd3c[3],cvd5c[3],cvd7c[3],cvd9c[3],cvd11c[3],cvd13c[3]],'r--',
#           [cvd13mgi[0][3],cvdafe2mgi[0][3],cvdafe3mgi[0][3]],[cvd13c[3],cvdafe2c[3],cvdafe3c[3]],'b-',
#           [cvdvarmgi[0][14],cvd13mgi[0][3],cvdvarmgi[0][15]],[cvdvarc[14],cvd13c[3],cvdvarc[15]],'g--')
#    P.xlabel('MgI0.88 [$\AA$]', fontsize=tfsize)
#    P.ylabel('g - r', fontsize=tfsize)
#    #P.xlim([0.35,0.65])
#    P.annotate('13', xy=(cvd13mgi[0][1],cvd13c[1]),xycoords='data',xytext=(-30,-10),textcoords='offset points', color='r', fontsize=lfsize)
#    P.annotate('11', xy=(cvd11mgi[0][1],cvd11c[1]),xycoords='data',xytext=(5,5),textcoords='offset points', color='r', fontsize=lfsize)
#    P.annotate('9', xy=(cvd9mgi[0][1],cvd9c[1]),xycoords='data',xytext=(10,0),textcoords='offset points', color='r', fontsize=lfsize)
#    P.annotate('7', xy=(cvd7mgi[0][1],cvd7c[1]),xycoords='data',xytext=(10,5),textcoords='offset points', color='r', fontsize=lfsize)
#    P.annotate('5', xy=(cvd5mgi[0][1],cvd5c[1]),xycoords='data',xytext=(15,0),textcoords='offset points', color='r', fontsize=lfsize)
#    P.annotate('3', xy=(cvd3mgi[0][1],cvd3c[1]),xycoords='data',xytext=(10,0),textcoords='offset points', color='r', fontsize=lfsize)
#    #P.annotate('x = 3.5', xy=(cvd13mgi[0][0],cvd13c[0]),xycoords='data',xytext=(-5,10),textcoords='offset points', color='k', fontsize=lfsize)
#    P.annotate('x = 3', xy=(cvd13mgi[0][1],cvd13c[1]),xycoords='data',xytext=(0,10),textcoords='offset points', color='k', fontsize=lfsize)
#    P.annotate('x = 2.35', xy=(cvd13mgi[0][2],cvd13c[2]),xycoords='data',xytext=(-20,12),textcoords='offset points', color='k', fontsize=lfsize)
#    P.annotate('Chab', xy=(cvd13mgi[0][3],cvd13c[3]),xycoords='data',xytext=(-15,10),textcoords='offset points', color='k', fontsize=lfsize)
#    P.annotate('B-L', xy=(cvd13mgi[0][4],cvd13c[4]),xycoords='data',xytext=(20,-10),textcoords='offset points', color='k', fontsize=lfsize)
#    P.annotate('[$\\alpha$/Fe]=+0.2', xy=(cvdafe2mgi[0][3],cvdafe2c[3]),xycoords='data',xytext=(-40,10),textcoords='offset points', color='b', fontsize=lfsize)
#    P.annotate('[$\\alpha$/Fe]=+0.3', xy=(cvdafe3mgi[0][3],cvdafe3c[3]),xycoords='data',xytext=(-40,5),textcoords='offset points', color='b', fontsize=lfsize)
#    P.annotate('[Mg/Fe]=+0.3', xy=(cvdvarmgi[0][14],cvdvarc[14]),xycoords='data',xytext=(-90,10),textcoords='offset points', color='g', fontsize=lfsize)
#    P.annotate('[Mg/Fe]=-0.3', xy=(cvdvarmgi[0][15],cvdvarc[15]),xycoords='data',xytext=(10,10),textcoords='offset points', color='g', fontsize=lfsize)
#
# elif sys.argv[1] == 'TiO':
#    #####TiO v g-r
#    P.scatter(cvd13tio[0][1:], cvd13c[1:], marker='o',c='k',s=[90,70,50,30])
#    P.scatter([cvd3tio[0][1],cvd5tio[0][1],cvd7tio[0][1],cvd9tio[0][1],cvd11tio[0][1],cvd13tio[0][1]],
#              [cvd3c[1],cvd5c[1],cvd7c[1],cvd9c[1],cvd11c[1],cvd13c[1]],
#              marker='s',c='r', s=[110,90,70,50,30,10])
#    P.scatter([cvd3tio[0][3],cvd5tio[0][3],cvd7tio[0][3],cvd9tio[0][3],cvd11tio[0][3],cvd13tio[0][3]],
#              [cvd3c[3],cvd5c[3],cvd7c[3],cvd9c[3],cvd11c[3],cvd13c[3]],
#              marker='s',c='r', s=[110,90,70,50,30,10])
#    P.scatter([cvd13tio[0][3],cvdafe2tio[0][3],cvdafe3tio[0][3]],
#              [cvd13c[3],cvdafe2c[3],cvdafe3c[3]],
#              marker='D',c='b',s=[10,30,50])
#    #Varying [Ti/Fe] abundance
#    P.scatter([cvdvartio[0][12],cvd13tio[0][3],cvdvartio[0][13]],
#              [cvdvarc[12],cvd13c[3],cvdvarc[13]],
#              marker='D',c='g',s=[50,5,50])
#    P.plot(cvd13tio[0][1:], cvd13c[1:], 'k-',
#           [cvd3tio[0][1],cvd5tio[0][1],cvd7tio[0][1],cvd9tio[0][1],cvd11tio[0][1],cvd13tio[0][1]],
#           [cvd3c[1],cvd5c[1],cvd7c[1],cvd9c[1],cvd11c[1],cvd13c[1]],'r--',
#           [cvd3tio[0][3],cvd5tio[0][3],cvd7tio[0][3],cvd9tio[0][3],cvd11tio[0][3],cvd13tio[0][3]],
#           [cvd3c[3],cvd5c[3],cvd7c[3],cvd9c[3],cvd11c[3],cvd13c[3]],'r--',
#           [cvd13tio[0][3],cvdafe2tio[0][3],cvdafe3tio[0][3]],[cvd13c[3],cvdafe2c[3],cvdafe3c[3]],'b-',
#           [cvdvartio[0][12],cvd13tio[0][3],cvdvartio[0][13]],[cvdvarc[12],cvd13c[3],cvdvarc[13]],'g--')
#    P.xlabel('TiO0.89', fontsize=tfsize)
#    P.ylabel('$g - r$', fontsize=tfsize)
#    P.annotate('13', xy=(cvd13tio[0][1],cvd13c[1]),xycoords='data',xytext=(10,-30),textcoords='offset points', color='r', fontsize=lfsize)
#    P.annotate('11', xy=(cvd11tio[0][1],cvd11c[1]),xycoords='data',xytext=(5,5),textcoords='offset points', color='r', fontsize=lfsize)
#    P.annotate('9', xy=(cvd9tio[0][1],cvd9c[1]),xycoords='data',xytext=(10,0),textcoords='offset points', color='r', fontsize=lfsize)
#    P.annotate('7', xy=(cvd7tio[0][1],cvd7c[1]),xycoords='data',xytext=(10,5),textcoords='offset points', color='r', fontsize=lfsize)
#    P.annotate('5', xy=(cvd5tio[0][1],cvd5c[1]),xycoords='data',xytext=(20,0),textcoords='offset points', color='r', fontsize=lfsize)
#    P.annotate('3 Gyr', xy=(cvd3tio[0][1],cvd3c[1]),xycoords='data',xytext=(10,-5),textcoords='offset points', color='r', fontsize=lfsize)
#    #P.annotate('x = 3.5', xy=(cvd13tio[0][0],cvd13c[0]),xycoords='data',xytext=(-5,10),textcoords='offset points', color='k', fontsize=lfsize)
#    P.annotate('x = 3', xy=(cvd13tio[0][1],cvd13c[1]),xycoords='data',xytext=(0,10),textcoords='offset points', color='k', fontsize=lfsize)
#    P.annotate('x = 2.35', xy=(cvd13tio[0][2],cvd13c[2]),xycoords='data',xytext=(-30,20),textcoords='offset points', color='k', fontsize=lfsize)
#    P.annotate('Chab', xy=(cvd13tio[0][3],cvd13c[3]),xycoords='data',xytext=(-25,5),textcoords='offset points', color='k', fontsize=lfsize)
#    P.annotate('B-L', xy=(cvd13tio[0][4],cvd13c[4]),xycoords='data',xytext=(20,-10),textcoords='offset points', color='k', fontsize=lfsize)
#    P.annotate('[$\\alpha$/Fe]=+0.2', xy=(cvdafe2tio[0][3],cvdafe2c[3]),xycoords='data',xytext=(-40,8),textcoords='offset points', color='b', fontsize=lfsize)
#    P.annotate('[$\\alpha$/Fe]=+0.3', xy=(cvdafe3tio[0][3],cvdafe3c[3]),xycoords='data',xytext=(-40,5),textcoords='offset points', color='b', fontsize=lfsize)
#    P.annotate('[Ti/Fe]=+0.15', xy=(cvdvartio[0][12],cvdvarc[12]),xycoords='data',xytext=(-30,20),textcoords='offset points', color='g', fontsize=lfsize)
#    P.annotate('[Ti/Fe]=-0.15', xy=(cvdvartio[0][13],cvdvarc[13]),xycoords='data',xytext=(-90,10),textcoords='offset points', color='g', fontsize=lfsize)
#
# elif sys.argv[1] == 'FeH':
#    ###FeH v g-r
#    P.scatter(cvd13feh[0][1:], cvd13c[1:], marker='o',c='k',s=[90,70,50,30])
#    P.scatter([cvd3feh[0][1],cvd5feh[0][1],cvd7feh[0][1],cvd9feh[0][1],cvd11feh[0][1],cvd13feh[0][1]],
#              [cvd3c[1],cvd5c[1],cvd7c[1],cvd9c[1],cvd11c[1],cvd13c[1]],
#              marker='s',c='r', s=[110,90,70,50,30,10])
#    P.scatter([cvd3feh[0][3],cvd5feh[0][3],cvd7feh[0][3],cvd9feh[0][3],cvd11feh[0][3],cvd13feh[0][3]],
#              [cvd3c[3],cvd5c[3],cvd7c[3],cvd9c[3],cvd11c[3],cvd13c[3]],
#              marker='s',c='r', s=[110,90,70,50,30,10])
#    P.scatter([cvd13feh[0][3],cvdafe2feh[0][3],cvdafe3feh[0][3]],
#              [cvd13c[3],cvdafe2c[3],cvdafe3c[3]],
#              marker='D',c='b',s=[10,30,50])
#    #Varying [Fe/H] abundance
#    P.scatter([cvdvarfeh[0][4],cvd13feh[0][3],cvdvarfeh[0][5]],
#              [cvdvarc[4],cvd13c[3],cvdvarc[5]],
#              marker='D',c='g',s=[50,5,50])
#    P.plot(cvd13feh[0][1:], cvd13c[1:], 'k-',
#           [cvd3feh[0][1],cvd5feh[0][1],cvd7feh[0][1],cvd9feh[0][1],cvd11feh[0][1],cvd13feh[0][1]],
#           [cvd3c[1],cvd5c[1],cvd7c[1],cvd9c[1],cvd11c[1],cvd13c[1]],'r--',
#           [cvd3feh[0][3],cvd5feh[0][3],cvd7feh[0][3],cvd9feh[0][3],cvd11feh[0][3],cvd13feh[0][3]],
#           [cvd3c[3],cvd5c[3],cvd7c[3],cvd9c[3],cvd11c[3],cvd13c[3]],'r--',
#           [cvd13feh[0][3],cvdafe2feh[0][3],cvdafe3feh[0][3]],[cvd13c[3],cvdafe2c[3],cvdafe3c[3]],'b-',
#           [cvdvarfeh[0][4],cvd13feh[0][3],cvdvarfeh[0][5]],[cvdvarc[4],cvd13c[3],cvdvarc[5]],'g--')
#    P.xlabel('FeH0.99 [$\AA$]', fontsize=tfsize)
#    P.ylabel('$g - r$', fontsize=tfsize)
#    P.annotate('13', xy=(cvd13feh[0][1],cvd13c[1]),xycoords='data',xytext=(10,-20),textcoords='offset points', color='r', fontsize=lfsize)
#    P.annotate('11', xy=(cvd11feh[0][1],cvd11c[1]),xycoords='data',xytext=(10,0),textcoords='offset points', color='r', fontsize=lfsize)
#    P.annotate('9', xy=(cvd9feh[0][1],cvd9c[1]),xycoords='data',xytext=(10,0),textcoords='offset points', color='r', fontsize=lfsize)
#    P.annotate('7', xy=(cvd7feh[0][1],cvd7c[1]),xycoords='data',xytext=(10,0),textcoords='offset points', color='r', fontsize=lfsize)
#    P.annotate('5', xy=(cvd5feh[0][1],cvd5c[1]),xycoords='data',xytext=(15,0),textcoords='offset points', color='r', fontsize=lfsize)
#    P.annotate('3 Gyr', xy=(cvd3feh[0][1],cvd3c[1]),xycoords='data',xytext=(10,0),textcoords='offset points', color='r', fontsize=lfsize)
#    #P.annotate('x = 3.5', xy=(cvd13feh[0][0],cvd13c[0]),xycoords='data',xytext=(-5,10),textcoords='offset points', color='k', fontsize=lfsize)
#    P.annotate('x = 3', xy=(cvd13feh[0][1],cvd13c[1]),xycoords='data',xytext=(0,10),textcoords='offset points', color='k', fontsize=lfsize)
#    P.annotate('x = 2.35', xy=(cvd13feh[0][2],cvd13c[2]),xycoords='data',xytext=(-25,12),textcoords='offset points', color='k', fontsize=lfsize)
#    P.annotate('Chab', xy=(cvd13feh[0][3],cvd13c[3]),xycoords='data',xytext=(-15,10),textcoords='offset points', color='k', fontsize=lfsize)
#    P.annotate('B-L', xy=(cvd13feh[0][4],cvd13c[4]),xycoords='data',xytext=(-20,10),textcoords='offset points', color='k', fontsize=lfsize)
#    P.annotate('[$\\alpha$/Fe]=+0.2', xy=(cvdafe2feh[0][3],cvdafe2c[3]),xycoords='data',xytext=(-125,-10),textcoords='offset points', color='b', fontsize=lfsize)
#    P.annotate('[$\\alpha$/Fe]=+0.3', xy=(cvdafe3feh[0][3],cvdafe3c[3]),xycoords='data',xytext=(-120,-30),textcoords='offset points', color='b', fontsize=lfsize)
#    P.annotate('[Fe/H]=+0.3', xy=(cvdvarfeh[0][4],cvdvarc[4]),xycoords='data',xytext=(-40,-20),textcoords='offset points', color='g', fontsize=lfsize)
#    P.annotate('[Fe/H]=-0.3', xy=(cvdvarfeh[0][5],cvdvarc[5]),xycoords='data',xytext=(-10,12),textcoords='offset points', color='g', fontsize=lfsize)
#
# ##    #####07-01-15: Plot M31 g-r and FeH datapoints!
# ##    m31coldat = n.genfromtxt('/Users/zieleniewski/Documents/Oxford/D.Phil/stelpops/Index_measurements/M31_sdss_colour_data/M31_SDSS_g-r_colour_data_Saglia2010_22-12-14.txt', delimiter=',')
# ##    m31dat = n.genfromtxt('/Users/zieleniewski/Documents/Oxford/D.Phil/stelpops/Index_measurements/17-12-14/FeH_index_values_200kms_with_bestfit_models_17-12-14.txt', delimiter=',')
# ##    m31moddat = n.genfromtxt('/Users/zieleniewski/Documents/Oxford/D.Phil/stelpops/Index_measurements/17-12-14/FeH_index_values_bestfitMODEL_200kms_17-12-14.txt', delimiter=',')
# ##    m31_1dat = n.genfromtxt('/Users/zieleniewski/Documents/Oxford/D.Phil/stelpops/Index_measurements/m31/FeH_index_values_mids_m31_1_21-12-14.txt', delimiter=',')
# ##    m31feh = n.zeros(11, dtype=float)
# ##    m31feherrs = n.zeros(11, dtype=float)
# ##    m31feh[:6] = m31_1dat[:,1]; m31feherrs[:6] = m31_1dat[:,2]
# ##    m31feh[6:-2] = m31dat[:3,1]; m31feherrs[6:-2] = m31dat[:3,2]
# ##    m31feh[-2] = m31moddat[3,1]; m31feherrs[-2] = m31moddat[3,2]
# ##    m31feh[-1] = m31dat[4,1]; m31feherrs[-1] = m31dat[4,2]
# ##
# ##    rad_dat = n.log10(m31coldat[:,0])
# ##
# ##    plot = P.scatter(m31feh, m31coldat[:,1], marker='o', s=[220,200,180,160,140,120,100,80,60,40], c=rad_dat, linewidths=0, cmap=P.cm.copper_r, zorder=10)
# ##    clb = P.colorbar(plot)
# ##    clb.set_label('log$_{10}$(r/arcsec)', fontsize=tfsize)
# ##    #convert time to a color tuple using the colormap used for scatter
# ##    rad_color = clb.to_rgba(rad_dat)
# ##    #pm0.05 RANDOM ERROR ADDED IN QUADRATURE - 28-05-15
# ##    m31feherrs[-5:] = n.sqrt(m31feherrs[-5:]**2+0.05**2)
# ##    a,b,c = P.errorbar(m31feh, m31coldat[:,1], xerr=m31feherrs, #yerr=m31coldat[:,2],
# ##                       marker='', ls='', capsize=0, barsabove=False, capthick=1.0,
# ##                       elinewidth=1.8, zorder=10)
# ##
# ##    #adjust the color of c[0], which is a LineCollection, to the colormap
# ##    c[0].set_color(rad_color)
# ##    #c[1].set_color(rad_color)
# ##
# ##    P.xlim([0.18, 0.9])
# ##    P.ylim([0.65, 0.9])
# ##    #####End
#
# else:
#    print 'Please choose one of: NaD, NaI, CaT, MgI, TiO, FeH'
#    sys.exit()
#
# #P.ylim([0.65,0.95])
# P.xticks(size=lfsize)
# P.yticks(size=lfsize)
# P.tick_params(axis='both', direction='in', length=majtsize,
#              width=1.5, which='major')
# P.tick_params(axis='both', direction='in', length=mintsize,
#              width=1.5, which='minor')
# P.show()
###------------#



#--------------#
###Plot MgI0.88 v g-r
##P.figure(figsize=(11,9),dpi=72)
##
##print cvd13mgi[0].shape
##print cvd13c.shape
##P.scatter(cvd13mgi[0], cvd13c, marker='o',c='k',s=[90,70,50,30,10])
##
##P.scatter([cvd3mgi[0][1],cvd5mgi[0][1],cvd7mgi[0][1],cvd9mgi[0][1],cvd11mgi[0][1],cvd13mgi[0][1]],
##          [cvd3c[1],cvd5c[1],cvd7c[1],cvd9c[1],cvd11c[1],cvd13c[1]],
##          marker='s',c='r', s=[110,90,70,50,30,10])
##
##P.scatter([cvd3mgi[0][3],cvd5mgi[0][3],cvd7mgi[0][3],cvd9mgi[0][3],cvd11mgi[0][3],cvd13mgi[0][3]],
##          [cvd3c[3],cvd5c[3],cvd7c[3],cvd9c[3],cvd11c[3],cvd13c[3]],
##          marker='s',c='r', s=[110,90,70,50,30,10])
##
##P.scatter([cvd13mgi[0][3],cvdafe2mgi[0][3],cvdafe3mgi[0][3]],
##          [cvd13c[3],cvdafe2c[3],cvdafe3c[3]],
##          marker='D',c='b',s=[10,30,50])
###Plot [Mg/Fe] abundance
##P.scatter([cvdvarmgi[0][14],cvd13mgi[0][3],cvdvarmgi[0][15]],
##          [cvdvarc[14],cvd13c[3],cvdvarc[15]],
##          marker='D',c='g',s=[50,5,50])
##
##P.plot(cvd13mgi[0], cvd13c, 'k-', [cvd3mgi[0][1],cvd5mgi[0][1],cvd7mgi[0][1],
##                                  cvd9mgi[0][1],cvd11mgi[0][1],cvd13mgi[0][1]],
##       [cvd3c[1],cvd5c[1],cvd7c[1],cvd9c[1],cvd11c[1],cvd13c[1]],'r--',
##       [cvd3mgi[0][3],cvd5mgi[0][3],cvd7mgi[0][3],
##        cvd9mgi[0][3],cvd11mgi[0][3],cvd13mgi[0][3]],
##       [cvd3c[3],cvd5c[3],cvd7c[3],
##        cvd9c[3],cvd11c[3],cvd13c[3]],'r--',
##       [cvd13mgi[0][3],cvdafe2mgi[0][3],cvdafe3mgi[0][3]],
##       [cvd13c[3],cvdafe2c[3],cvdafe3c[3]],'b-',
##       [cvdvarmgi[0][14],cvd13mgi[0][3],cvdvarmgi[0][15]],
##          [cvdvarc[14],cvd13c[3],cvdvarc[15]],'g--')
##
##P.xlabel('MgI0.88 [$\AA$]', fontsize=22)
##
##P.ylabel('$g - r$', fontsize=22)
##
##P.xlim([0.35,0.65])
##
##P.xticks(size=18)
##
##P.yticks(size=18)
##
##P.tick_params(axis='both', which='major', length=12, width=1.5, direction='in')
##P.tick_params(axis='both', which='minor', length=8, width=1.5, direction='in')
##P.annotate('13 Gyr', xy=(cvd13mgi[0][1],cvd13c[1]),xycoords='data',xytext=(10,-30),textcoords='offset points',fontsize=14, color='r')
##P.annotate('11 Gyr', xy=(cvd11mgi[0][1],cvd11c[1]),xycoords='data',xytext=(5,5),textcoords='offset points',fontsize=14, color='r')
##P.annotate('9 Gyr', xy=(cvd9mgi[0][1],cvd9c[1]),xycoords='data',xytext=(10,0),textcoords='offset points',fontsize=14, color='r')
##P.annotate('7 Gyr', xy=(cvd7mgi[0][1],cvd7c[1]),xycoords='data',xytext=(10,5),textcoords='offset points',fontsize=14, color='r')
##P.annotate('5 Gyr', xy=(cvd5mgi[0][1],cvd5c[1]),xycoords='data',xytext=(15,0),textcoords='offset points',fontsize=14, color='r')
##P.annotate('3 Gyr', xy=(cvd3mgi[0][1],cvd3c[1]),xycoords='data',xytext=(10,0),textcoords='offset points',fontsize=14, color='r')
##P.annotate('x = 3.5', xy=(cvd13mgi[0][0],cvd13c[0]),xycoords='data',xytext=(-5,10),textcoords='offset points',fontsize=14, color='k')
##P.annotate('x = 3', xy=(cvd13mgi[0][1],cvd13c[1]),xycoords='data',xytext=(0,10),textcoords='offset points',fontsize=14, color='k')
##P.annotate('x = 2.35', xy=(cvd13mgi[0][2],cvd13c[2]),xycoords='data',xytext=(-25,12),textcoords='offset points',fontsize=14, color='k')
##P.annotate('x = Chabrier', xy=(cvd13mgi[0][3],cvd13c[3]),xycoords='data',xytext=(-15,10),textcoords='offset points',fontsize=14, color='k')
##P.annotate('x = B-L', xy=(cvd13mgi[0][4],cvd13c[4]),xycoords='data',xytext=(20,-10),textcoords='offset points',fontsize=14, color='k')
##P.annotate('[$\\alpha$/Fe]=+0.2', xy=(cvdafe2mgi[0][3],cvdafe2c[3]),xycoords='data',xytext=(-40,10),textcoords='offset points',fontsize=14, color='b')
##P.annotate('[$\\alpha$/Fe]=+0.3', xy=(cvdafe3mgi[0][3],cvdafe3c[3]),xycoords='data',xytext=(-40,10),textcoords='offset points',fontsize=14, color='b')
##P.annotate('[Mg/Fe]=+0.3', xy=(cvdvarmgi[0][14],cvdvarc[14]),xycoords='data',xytext=(-60,5),textcoords='offset points',fontsize=14, color='g')
##P.annotate('[Mg/Fe]=-0.3', xy=(cvdvarmgi[0][15],cvdvarc[15]),xycoords='data',xytext=(10,5),textcoords='offset points',fontsize=14, color='g')
##
##P.show()
##
###Plot NaI0.82 v g-r
##P.figure(figsize=(11,9),dpi=72)
##
##print cvd13mgi[0].shape
##print cvd13c.shape
##P.scatter(cvd13nai[0], cvd13c, marker='o',c='k',s=[90,70,50,30,10])
##
##P.scatter([cvd3nai[0][1],cvd5nai[0][1],cvd7nai[0][1],cvd9nai[0][1],cvd11nai[0][1],cvd13nai[0][1]],
##          [cvd3c[1],cvd5c[1],cvd7c[1],cvd9c[1],cvd11c[1],cvd13c[1]],
##          marker='s',c='r', s=[110,90,70,50,30,10])
##
##P.scatter([cvd3nai[0][3],cvd5nai[0][3],cvd7nai[0][3],cvd9nai[0][3],cvd11nai[0][3],cvd13nai[0][3]],
##          [cvd3c[3],cvd5c[3],cvd7c[3],cvd9c[3],cvd11c[3],cvd13c[3]],
##          marker='s',c='r', s=[110,90,70,50,30,10])
##
##P.scatter([cvd13nai[0][3],cvdafe2nai[0][3],cvdafe3nai[0][3]],
##          [cvd13c[3],cvdafe2c[3],cvdafe3c[3]],
##          marker='D',c='b',s=[10,30,50])
###Varying [Na/Fe] abundance
##print [cvdvarnai[0][0],cvd13nai[0][3],cvdvarnai[0][1]]
##print [cvdvarc[0],cvd13c[3],cvdvarc[1]]
##P.scatter([cvdvarnai[0][0],cvd13nai[0][3],cvdvarnai[0][1]],
##          [cvdvarc[0],cvd13c[3],cvdvarc[1]],
##          marker='D',c='g',s=[50,5,50])
##
##P.plot(cvd13nai[0], cvd13c, 'k-', [cvd3nai[0][1],cvd5nai[0][1],cvd7nai[0][1],
##                                  cvd9nai[0][1],cvd11nai[0][1],cvd13nai[0][1]],
##       [cvd3c[1],cvd5c[1],cvd7c[1],cvd9c[1],cvd11c[1],cvd13c[1]],'r--',
##       [cvd3nai[0][3],cvd5nai[0][3],cvd7nai[0][3],
##        cvd9nai[0][3],cvd11nai[0][3],cvd13nai[0][3]],
##       [cvd3c[3],cvd5c[3],cvd7c[3],
##        cvd9c[3],cvd11c[3],cvd13c[3]],'r--',
##       [cvd13nai[0][3],cvdafe2nai[0][3],cvdafe3nai[0][3]],
##       [cvd13c[3],cvdafe2c[3],cvdafe3c[3]],'b-',
##       [cvdvarnai[0][0],cvd13nai[0][3],cvdvarnai[0][1]],
##          [cvdvarc[0],cvd13c[3],cvdvarc[1]],'g--')
##
##
##P.xlabel('NaI0.82 [$\AA$]', fontsize=22)
##
##P.ylabel('$g - r$', fontsize=22)
##
##P.xticks(size=18)
##
##P.yticks(size=18)
##
##P.tick_params(axis='both', which='major', length=12, width=1.5, direction='in')
##P.tick_params(axis='both', which='minor', length=8, width=1.5, direction='in')
##P.annotate('13 Gyr', xy=(cvd13nai[0][1],cvd13c[1]),xycoords='data',xytext=(10,-30),textcoords='offset points',fontsize=14, color='r')
##P.annotate('11 Gyr', xy=(cvd11nai[0][1],cvd11c[1]),xycoords='data',xytext=(5,5),textcoords='offset points',fontsize=14, color='r')
##P.annotate('9 Gyr', xy=(cvd9nai[0][1],cvd9c[1]),xycoords='data',xytext=(10,0),textcoords='offset points',fontsize=14, color='r')
##P.annotate('7 Gyr', xy=(cvd7nai[0][1],cvd7c[1]),xycoords='data',xytext=(10,5),textcoords='offset points',fontsize=14, color='r')
##P.annotate('5 Gyr', xy=(cvd5nai[0][1],cvd5c[1]),xycoords='data',xytext=(15,0),textcoords='offset points',fontsize=14, color='r')
##P.annotate('3 Gyr', xy=(cvd3nai[0][1],cvd3c[1]),xycoords='data',xytext=(10,0),textcoords='offset points',fontsize=14, color='r')
##P.annotate('x = 3.5', xy=(cvd13nai[0][0],cvd13c[0]),xycoords='data',xytext=(-5,10),textcoords='offset points',fontsize=14, color='k')
##P.annotate('x = 3', xy=(cvd13nai[0][1],cvd13c[1]),xycoords='data',xytext=(0,10),textcoords='offset points',fontsize=14, color='k')
##P.annotate('x = 2.35', xy=(cvd13nai[0][2],cvd13c[2]),xycoords='data',xytext=(-25,12),textcoords='offset points',fontsize=14, color='k')
##P.annotate('x = Chabrier', xy=(cvd13nai[0][3],cvd13c[3]),xycoords='data',xytext=(-20,10),textcoords='offset points',fontsize=14, color='k')
##P.annotate('x = B-L', xy=(cvd13nai[0][4],cvd13c[4]),xycoords='data',xytext=(-20,10),textcoords='offset points',fontsize=14, color='k')
##P.annotate('[$\\alpha$/Fe]=+0.2', xy=(cvdafe2nai[0][3],cvdafe2c[3]),xycoords='data',xytext=(-60,0),textcoords='offset points',fontsize=14, color='b')
##P.annotate('[$\\alpha$/Fe]=+0.3', xy=(cvdafe3nai[0][3],cvdafe3c[3]),xycoords='data',xytext=(-40,-30),textcoords='offset points',fontsize=14, color='b')
##P.annotate('[Na/Fe]=+0.3', xy=(cvdvarnai[0][0],cvdvarc[0]),xycoords='data',xytext=(-30,-20),textcoords='offset points',fontsize=14, color='g')
##P.annotate('[Na/Fe]=-0.3', xy=(cvdvarnai[0][1],cvdvarc[1]),xycoords='data',xytext=(-90,0),textcoords='offset points',fontsize=14, color='g')
##
##
##P.show()
##
##
###Plot FeH0.99 v g-r
##P.figure(figsize=(11,9),dpi=72)
##
##print cvd13mgi[0].shape
##print cvd13c.shape
##P.scatter(cvd13feh[0], cvd13c, marker='o',c='k',s=[90,70,50,30,10])
##
##P.scatter([cvd3feh[0][1],cvd5feh[0][1],cvd7feh[0][1],cvd9feh[0][1],cvd11feh[0][1],cvd13feh[0][1]],
##          [cvd3c[1],cvd5c[1],cvd7c[1],cvd9c[1],cvd11c[1],cvd13c[1]],
##          marker='s',c='r', s=[110,90,70,50,30,10])
##
##P.scatter([cvd3feh[0][3],cvd5feh[0][3],cvd7feh[0][3],cvd9feh[0][3],cvd11feh[0][3],cvd13feh[0][3]],
##          [cvd3c[3],cvd5c[3],cvd7c[3],cvd9c[3],cvd11c[3],cvd13c[3]],
##          marker='s',c='r', s=[110,90,70,50,30,10])
##
##P.scatter([cvd13feh[0][3],cvdafe2feh[0][3],cvdafe3feh[0][3]],
##          [cvd13c[3],cvdafe2c[3],cvdafe3c[3]],
##          marker='D',c='b',s=[10,30,50])
###Varying [Fe/H] abundance
##P.scatter([cvdvarfeh[0][4],cvd13feh[0][3],cvdvarfeh[0][5]],
##          [cvdvarc[4],cvd13c[3],cvdvarc[5]],
##          marker='D',c='g',s=[50,5,50])
##
##P.plot(cvd13feh[0], cvd13c, 'k-', [cvd3feh[0][1],cvd5feh[0][1],cvd7feh[0][1],
##                                  cvd9feh[0][1],cvd11feh[0][1],cvd13feh[0][1]],
##       [cvd3c[1],cvd5c[1],cvd7c[1],cvd9c[1],cvd11c[1],cvd13c[1]],'r--',
##       [cvd3feh[0][3],cvd5feh[0][3],cvd7feh[0][3],
##        cvd9feh[0][3],cvd11feh[0][3],cvd13feh[0][3]],
##       [cvd3c[3],cvd5c[3],cvd7c[3],
##        cvd9c[3],cvd11c[3],cvd13c[3]],'r--',
##       [cvd13feh[0][3],cvdafe2feh[0][3],cvdafe3feh[0][3]],
##       [cvd13c[3],cvdafe2c[3],cvdafe3c[3]],'b-',
##       [cvdvarfeh[0][4],cvd13feh[0][3],cvdvarfeh[0][5]],
##          [cvdvarc[4],cvd13c[3],cvdvarc[5]],'g--')
##
##P.xlabel('FeH0.99 [$\AA$]', fontsize=22)
##
##P.ylabel('$g - r$', fontsize=22)
##
##P.xticks(size=18)
##
##P.yticks(size=18)
##
##P.tick_params(axis='both', which='major', length=12, width=1.5, direction='in')
##P.tick_params(axis='both', which='minor', length=8, width=1.5, direction='in')
##P.annotate('13 Gyr', xy=(cvd13feh[0][1],cvd13c[1]),xycoords='data',xytext=(10,-30),textcoords='offset points',fontsize=14, color='r')
##P.annotate('11 Gyr', xy=(cvd11feh[0][1],cvd11c[1]),xycoords='data',xytext=(5,5),textcoords='offset points',fontsize=14, color='r')
##P.annotate('9 Gyr', xy=(cvd9feh[0][1],cvd9c[1]),xycoords='data',xytext=(10,0),textcoords='offset points',fontsize=14, color='r')
##P.annotate('7 Gyr', xy=(cvd7feh[0][1],cvd7c[1]),xycoords='data',xytext=(10,5),textcoords='offset points',fontsize=14, color='r')
##P.annotate('5 Gyr', xy=(cvd5feh[0][1],cvd5c[1]),xycoords='data',xytext=(15,0),textcoords='offset points',fontsize=14, color='r')
##P.annotate('3 Gyr', xy=(cvd3feh[0][1],cvd3c[1]),xycoords='data',xytext=(10,0),textcoords='offset points',fontsize=14, color='r')
##P.annotate('x = 3.5', xy=(cvd13feh[0][0],cvd13c[0]),xycoords='data',xytext=(-5,10),textcoords='offset points',fontsize=14, color='k')
##P.annotate('x = 3', xy=(cvd13feh[0][1],cvd13c[1]),xycoords='data',xytext=(0,10),textcoords='offset points',fontsize=14, color='k')
##P.annotate('x = 2.35', xy=(cvd13feh[0][2],cvd13c[2]),xycoords='data',xytext=(-25,12),textcoords='offset points',fontsize=14, color='k')
##P.annotate('x = Chabrier', xy=(cvd13feh[0][3],cvd13c[3]),xycoords='data',xytext=(-15,10),textcoords='offset points',fontsize=14, color='k')
##P.annotate('x = B-L', xy=(cvd13feh[0][4],cvd13c[4]),xycoords='data',xytext=(20,-10),textcoords='offset points',fontsize=14, color='k')
##P.annotate('[$\\alpha$/Fe]=+0.2', xy=(cvdafe2feh[0][3],cvdafe2c[3]),xycoords='data',xytext=(-40,0),textcoords='offset points',fontsize=14, color='b')
##P.annotate('[$\\alpha$/Fe]=+0.3', xy=(cvdafe3feh[0][3],cvdafe3c[3]),xycoords='data',xytext=(-60,-10),textcoords='offset points',fontsize=14, color='b')
##P.annotate('[Fe/H]=+0.3', xy=(cvdvarfeh[0][4],cvdvarc[4]),xycoords='data',xytext=(-30,-20),textcoords='offset points',fontsize=14, color='g')
##P.annotate('[Fe/H]=-0.3', xy=(cvdvarfeh[0][5],cvdvarc[5]),xycoords='data',xytext=(-10,10),textcoords='offset points',fontsize=14, color='g')
##
##P.show()
##
###Plot CaT* v g-r
##P.figure(figsize=(11,9),dpi=72)
##
##P.scatter(cvd13cat[0], cvd13c, marker='o',c='k',s=[90,70,50,30,10])
##
##P.scatter([cvd3cat[0][1],cvd5cat[0][1],cvd7cat[0][1],cvd9cat[0][1],cvd11cat[0][1],cvd13cat[0][1]],
##          [cvd3c[1],cvd5c[1],cvd7c[1],cvd9c[1],cvd11c[1],cvd13c[1]],
##          marker='s',c='r', s=[110,90,70,50,30,10])
##
##P.scatter([cvd3cat[0][3],cvd5cat[0][3],cvd7cat[0][3],cvd9cat[0][3],cvd11cat[0][3],cvd13cat[0][3]],
##          [cvd3c[3],cvd5c[3],cvd7c[3],cvd9c[3],cvd11c[3],cvd13c[3]],
##          marker='s',c='r', s=[110,90,70,50,30,10])
##
##P.scatter([cvd13cat[0][3],cvdafe2cat[0][3],cvdafe3cat[0][3]],
##          [cvd13c[3],cvdafe2c[3],cvdafe3c[3]],
##          marker='D',c='b',s=[10,30,50])
###Varying [Ca/Fe] abundance
##P.scatter([cvdvarcat[0][2],cvd13cat[0][3],cvdvarcat[0][3]],
##          [cvdvarc[2],cvd13c[3],cvdvarc[3]],
##          marker='D',c='g',s=[50,5,50])
##
##P.plot(cvd13cat[0], cvd13c, 'k-', [cvd3cat[0][1],cvd5cat[0][1],cvd7cat[0][1],
##                                  cvd9cat[0][1],cvd11cat[0][1],cvd13cat[0][1]],
##       [cvd3c[1],cvd5c[1],cvd7c[1],cvd9c[1],cvd11c[1],cvd13c[1]],'r--',
##       [cvd3cat[0][3],cvd5cat[0][3],cvd7cat[0][3],
##        cvd9cat[0][3],cvd11cat[0][3],cvd13cat[0][3]],
##       [cvd3c[3],cvd5c[3],cvd7c[3],
##        cvd9c[3],cvd11c[3],cvd13c[3]],'r--',
##       [cvd13cat[0][3],cvdafe2cat[0][3],cvdafe3cat[0][3]],
##       [cvd13c[3],cvdafe2c[3],cvdafe3c[3]],'b-',
##       [cvdvarcat[0][2],cvd13cat[0][3],cvdvarcat[0][3]],
##          [cvdvarc[2],cvd13c[3],cvdvarc[3]],'g--')
##
##P.xlabel('CaT* [$\AA$]', fontsize=22)
##
##P.ylabel('$g - r$', fontsize=22)
##
##P.xticks(size=18)
##
##P.yticks(size=18)
##
##P.tick_params(axis='both', which='major', length=12, width=1.5, direction='in')
##P.tick_params(axis='both', which='minor', length=8, width=1.5, direction='in')
##P.annotate('13 Gyr', xy=(cvd13cat[0][1],cvd13c[1]),xycoords='data',xytext=(10,-30),textcoords='offset points',fontsize=14, color='r')
##P.annotate('11 Gyr', xy=(cvd11cat[0][1],cvd11c[1]),xycoords='data',xytext=(5,5),textcoords='offset points',fontsize=14, color='r')
##P.annotate('9 Gyr', xy=(cvd9cat[0][1],cvd9c[1]),xycoords='data',xytext=(10,0),textcoords='offset points',fontsize=14, color='r')
##P.annotate('7 Gyr', xy=(cvd7cat[0][1],cvd7c[1]),xycoords='data',xytext=(10,5),textcoords='offset points',fontsize=14, color='r')
##P.annotate('5 Gyr', xy=(cvd5cat[0][1],cvd5c[1]),xycoords='data',xytext=(15,0),textcoords='offset points',fontsize=14, color='r')
##P.annotate('3 Gyr', xy=(cvd3cat[0][1],cvd3c[1]),xycoords='data',xytext=(10,0),textcoords='offset points',fontsize=14, color='r')
##P.annotate('x = 3.5', xy=(cvd13cat[0][0],cvd13c[0]),xycoords='data',xytext=(-5,10),textcoords='offset points',fontsize=14, color='k')
##P.annotate('x = 3', xy=(cvd13cat[0][1],cvd13c[1]),xycoords='data',xytext=(0,10),textcoords='offset points',fontsize=14, color='k')
##P.annotate('x = 2.35', xy=(cvd13cat[0][2],cvd13c[2]),xycoords='data',xytext=(-25,12),textcoords='offset points',fontsize=14, color='k')
##P.annotate('x = Chabrier', xy=(cvd13cat[0][3],cvd13c[3]),xycoords='data',xytext=(-15,10),textcoords='offset points',fontsize=14, color='k')
##P.annotate('x = B-L', xy=(cvd13cat[0][4],cvd13c[4]),xycoords='data',xytext=(20,-10),textcoords='offset points',fontsize=14, color='k')
##P.annotate('[$\\alpha$/Fe]=+0.2', xy=(cvdafe2cat[0][3],cvdafe2c[3]),xycoords='data',xytext=(-40,10),textcoords='offset points',fontsize=14, color='b')
##P.annotate('[$\\alpha$/Fe]=+0.3', xy=(cvdafe3cat[0][3],cvdafe3c[3]),xycoords='data',xytext=(-40,10),textcoords='offset points',fontsize=14, color='b')
##P.annotate('[Ca/Fe]=+0.15', xy=(cvdvarcat[0][2],cvdvarc[2]),xycoords='data',xytext=(-30,10),textcoords='offset points',fontsize=14, color='g')
##P.annotate('[Ca/Fe]=-0.15', xy=(cvdvarcat[0][3],cvdvarc[3]),xycoords='data',xytext=(-15,10),textcoords='offset points',fontsize=14, color='g')
##
##P.show()
##
##
###Plot TiO0.89 v g-r
##P.figure(figsize=(11,9),dpi=72)
##
##P.scatter(cvd13tio[0], cvd13c, marker='o',c='k',s=[90,70,50,30,10])
##
##P.scatter([cvd3tio[0][1],cvd5tio[0][1],cvd7tio[0][1],cvd9tio[0][1],cvd11tio[0][1],cvd13tio[0][1]],
##          [cvd3c[1],cvd5c[1],cvd7c[1],cvd9c[1],cvd11c[1],cvd13c[1]],
##          marker='s',c='r', s=[110,90,70,50,30,10])
##
##P.scatter([cvd3tio[0][3],cvd5tio[0][3],cvd7tio[0][3],cvd9tio[0][3],cvd11tio[0][3],cvd13tio[0][3]],
##          [cvd3c[3],cvd5c[3],cvd7c[3],cvd9c[3],cvd11c[3],cvd13c[3]],
##          marker='s',c='r', s=[110,90,70,50,30,10])
##
##P.scatter([cvd13tio[0][3],cvdafe2tio[0][3],cvdafe3tio[0][3]],
##          [cvd13c[3],cvdafe2c[3],cvdafe3c[3]],
##          marker='D',c='b',s=[10,30,50])
###Varying [Ti/Fe] abundance
##P.scatter([cvdvartio[0][12],cvd13tio[0][3],cvdvartio[0][13]],
##          [cvdvarc[12],cvd13c[3],cvdvarc[13]],
##          marker='D',c='g',s=[50,5,50])
##
##P.plot(cvd13tio[0], cvd13c, 'k-', [cvd3tio[0][1],cvd5tio[0][1],cvd7tio[0][1],
##                                  cvd9tio[0][1],cvd11tio[0][1],cvd13tio[0][1]],
##       [cvd3c[1],cvd5c[1],cvd7c[1],cvd9c[1],cvd11c[1],cvd13c[1]],'r--',
##       [cvd3tio[0][3],cvd5tio[0][3],cvd7tio[0][3],
##        cvd9tio[0][3],cvd11tio[0][3],cvd13tio[0][3]],
##       [cvd3c[3],cvd5c[3],cvd7c[3],
##        cvd9c[3],cvd11c[3],cvd13c[3]],'r--',
##       [cvd13tio[0][3],cvdafe2tio[0][3],cvdafe3tio[0][3]],
##       [cvd13c[3],cvdafe2c[3],cvdafe3c[3]],'b-',
##       [cvdvartio[0][12],cvd13tio[0][3],cvdvartio[0][13]],
##          [cvdvarc[12],cvd13c[3],cvdvarc[13]],'g--')
##
##P.xlabel('TiO0.89', fontsize=22)
##
##P.ylabel('$g - r$', fontsize=22)
##
##P.xticks(size=18)
##
##P.yticks(size=18)
##
##P.tick_params(axis='both', which='major', length=12, width=1.5, direction='in')
##P.tick_params(axis='both', which='minor', length=8, width=1.5, direction='in')
##P.annotate('13', xy=(cvd13tio[0][1],cvd13c[1]),xycoords='data',xytext=(10,-30),textcoords='offset points',fontsize=14, color='r')
##P.annotate('11', xy=(cvd11tio[0][1],cvd11c[1]),xycoords='data',xytext=(5,5),textcoords='offset points',fontsize=14, color='r')
##P.annotate('9', xy=(cvd9tio[0][1],cvd9c[1]),xycoords='data',xytext=(10,0),textcoords='offset points',fontsize=14, color='r')
##P.annotate('7', xy=(cvd7tio[0][1],cvd7c[1]),xycoords='data',xytext=(10,5),textcoords='offset points',fontsize=14, color='r')
##P.annotate('5', xy=(cvd5tio[0][1],cvd5c[1]),xycoords='data',xytext=(15,0),textcoords='offset points',fontsize=14, color='r')
##P.annotate('3 Gyr', xy=(cvd3tio[0][1],cvd3c[1]),xycoords='data',xytext=(10,0),textcoords='offset points',fontsize=14, color='r')
##P.annotate('x = 3.5', xy=(cvd13tio[0][0],cvd13c[0]),xycoords='data',xytext=(-5,10),textcoords='offset points',fontsize=14, color='k')
##P.annotate('x = 3', xy=(cvd13tio[0][1],cvd13c[1]),xycoords='data',xytext=(0,10),textcoords='offset points',fontsize=14, color='k')
##P.annotate('x = 2.35', xy=(cvd13tio[0][2],cvd13c[2]),xycoords='data',xytext=(-25,12),textcoords='offset points',fontsize=14, color='k')
##P.annotate('Chab', xy=(cvd13tio[0][3],cvd13c[3]),xycoords='data',xytext=(-15,10),textcoords='offset points',fontsize=14, color='k')
##P.annotate('B-L', xy=(cvd13tio[0][4],cvd13c[4]),xycoords='data',xytext=(20,-10),textcoords='offset points',fontsize=14, color='k')
##P.annotate('[$\\alpha$/Fe]=+0.2', xy=(cvdafe2tio[0][3],cvdafe2c[3]),xycoords='data',xytext=(-40,10),textcoords='offset points',fontsize=14, color='b')
##P.annotate('[$\\alpha$/Fe]=+0.3', xy=(cvdafe3tio[0][3],cvdafe3c[3]),xycoords='data',xytext=(-40,10),textcoords='offset points',fontsize=14, color='b')
##P.annotate('[Ti/Fe]=+0.15', xy=(cvdvartio[0][12],cvdvarc[12]),xycoords='data',xytext=(-30,-20),textcoords='offset points',fontsize=14, color='g')
##P.annotate('[Ti/Fe]=-0.15', xy=(cvdvartio[0][13],cvdvarc[13]),xycoords='data',xytext=(-90,0),textcoords='offset points',fontsize=14, color='g')
##
##P.show()
##
##
##
##
###Plot NaI0.82 v MgI0.88
##P.figure(figsize=(11,9),dpi=72)
###Changing IMF at 13.5 Gyr
##P.scatter(cvd13nai[0], cvd13mgi[0], marker='o',c='k',s=[90,70,50,30,10])
###Changing age for x=3 IMF
##P.scatter([cvd3nai[0][1],cvd5nai[0][1],cvd7nai[0][1],cvd9nai[0][1],cvd11nai[0][1],cvd13nai[0][1]],
##          [cvd3mgi[0][1],cvd5mgi[0][1],cvd7mgi[0][1],cvd9mgi[0][1],cvd11mgi[0][1],cvd13mgi[0][1]],
##          marker='s',c='r', s=[110,90,70,50,30,10])
###Changing age for Chabrier IMF
##P.scatter([cvd3nai[0][3],cvd5nai[0][3],cvd7nai[0][3],cvd9nai[0][3],cvd11nai[0][3],cvd13nai[0][3]],
##          [cvd3mgi[0][3],cvd5mgi[0][3],cvd7mgi[0][3],cvd9mgi[0][3],cvd11mgi[0][3],cvd13mgi[0][3]],
##          marker='s',c='r', s=[110,90,70,50,30,10])
###Changing [alpha/Fe] for 13.5 Gyr Chabrier IMF
##P.scatter([cvd13nai[0][3],cvdafe2nai[0][3],cvdafe3nai[0][3]],
##          [cvd13mgi[0][3],cvdafe2mgi[0][3],cvdafe3mgi[0][3]],
##          marker='D',c='b',s=[10,30,50])
##
##P.plot(cvd13nai[0], cvd13mgi[0], 'k-', [cvd3nai[0][1],cvd5nai[0][1],cvd7nai[0][1],
##                                  cvd9nai[0][1],cvd11nai[0][1],cvd13nai[0][1]],
##       [cvd3mgi[0][1],cvd5mgi[0][1],cvd7mgi[0][1],cvd9mgi[0][1],cvd11mgi[0][1],cvd13mgi[0][1]],'r--',
##       [cvd3nai[0][3],cvd5nai[0][3],cvd7nai[0][3],
##        cvd9nai[0][3],cvd11nai[0][3],cvd13nai[0][3]],
##       [cvd3mgi[0][3],cvd5mgi[0][3],cvd7mgi[0][3],
##        cvd9mgi[0][3],cvd11mgi[0][3],cvd13mgi[0][3]],'r--',
##       [cvd13nai[0][3],cvdafe2nai[0][3],cvdafe3nai[0][3]],
##       [cvd13mgi[0][3],cvdafe2mgi[0][3],cvdafe3mgi[0][3]],'b-')
##
##P.xlabel('NaI0.82 [$\AA$]', fontsize=22)
##
##P.ylabel('MgI0.88 [$\AA$]', fontsize=22)
##
##P.xticks(size=18)
##
##P.yticks(size=18)
##
##P.tick_params(axis='both', which='major', length=12, width=1.5, direction='in')
##P.tick_params(axis='both', which='minor', length=8, width=1.5, direction='in')
##P.annotate('13 Gyr', xy=(cvd13nai[0][1],cvd13mgi[0][1]),xycoords='data',xytext=(10,-30),textcoords='offset points',fontsize=14, color='r')
##P.annotate('11 Gyr', xy=(cvd11nai[0][1],cvd11mgi[0][1]),xycoords='data',xytext=(5,5),textcoords='offset points',fontsize=14, color='r')
##P.annotate('9 Gyr', xy=(cvd9nai[0][1],cvd9mgi[0][1]),xycoords='data',xytext=(10,0),textcoords='offset points',fontsize=14, color='r')
##P.annotate('7 Gyr', xy=(cvd7nai[0][1],cvd7mgi[0][1]),xycoords='data',xytext=(10,5),textcoords='offset points',fontsize=14, color='r')
##P.annotate('5 Gyr', xy=(cvd5nai[0][1],cvd5mgi[0][1]),xycoords='data',xytext=(15,0),textcoords='offset points',fontsize=14, color='r')
##P.annotate('3 Gyr', xy=(cvd3nai[0][1],cvd3mgi[0][1]),xycoords='data',xytext=(10,0),textcoords='offset points',fontsize=14, color='r')
##P.annotate('x = 3.5', xy=(cvd13nai[0][0],cvd13mgi[0][0]),xycoords='data',xytext=(-5,10),textcoords='offset points',fontsize=14, color='k')
##P.annotate('x = 3', xy=(cvd13nai[0][1],cvd13mgi[0][1]),xycoords='data',xytext=(0,10),textcoords='offset points',fontsize=14, color='k')
##P.annotate('x = 2.35', xy=(cvd13nai[0][2],cvd13mgi[0][2]),xycoords='data',xytext=(-25,12),textcoords='offset points',fontsize=14, color='k')
##P.annotate('x = Chabrier', xy=(cvd13nai[0][3],cvd13mgi[0][3]),xycoords='data',xytext=(-15,10),textcoords='offset points',fontsize=14, color='k')
##P.annotate('x = B-L', xy=(cvd13nai[0][4],cvd13mgi[0][4]),xycoords='data',xytext=(20,-10),textcoords='offset points',fontsize=14, color='k')
##P.annotate('[$\\alpha$/Fe]=+0.2', xy=(cvdafe2nai[0][3],cvdafe2mgi[0][3]),xycoords='data',xytext=(-40,10),textcoords='offset points',fontsize=14, color='b')
##P.annotate('[$\\alpha$/Fe]=+0.3', xy=(cvdafe3nai[0][3],cvdafe3mgi[0][3]),xycoords='data',xytext=(-40,10),textcoords='offset points',fontsize=14, color='b')
##
##
##P.show()
##
##
##
###Plot TiO0.89 v MgI0.88
##P.figure(figsize=(11,9),dpi=72)
###Changing IMF at 13.5 Gyr
##P.scatter(cvd13tio[0], cvd13mgi[0], marker='o',c='k',s=[90,70,50,30,10])
###Changing age for x=3 IMF
##P.scatter([cvd3tio[0][1],cvd5tio[0][1],cvd7tio[0][1],cvd9tio[0][1],cvd11tio[0][1],cvd13tio[0][1]],
##          [cvd3mgi[0][1],cvd5mgi[0][1],cvd7mgi[0][1],cvd9mgi[0][1],cvd11mgi[0][1],cvd13mgi[0][1]],
##          marker='s',c='r', s=[110,90,70,50,30,10])
###Changing age for Chabrier IMF
##P.scatter([cvd3tio[0][3],cvd5tio[0][3],cvd7tio[0][3],cvd9tio[0][3],cvd11tio[0][3],cvd13tio[0][3]],
##          [cvd3mgi[0][3],cvd5mgi[0][3],cvd7mgi[0][3],cvd9mgi[0][3],cvd11mgi[0][3],cvd13mgi[0][3]],
##          marker='s',c='r', s=[110,90,70,50,30,10])
###Changing [alpha/Fe] for 13.5 Gyr Chabrier IMF
##P.scatter([cvd13tio[0][3],cvdafe2tio[0][3],cvdafe3tio[0][3]],
##          [cvd13mgi[0][3],cvdafe2mgi[0][3],cvdafe3mgi[0][3]],
##          marker='D',c='b',s=[10,30,50])
##
##P.plot(cvd13tio[0], cvd13mgi[0], 'k-', [cvd3tio[0][1],cvd5tio[0][1],cvd7tio[0][1],
##                                  cvd9tio[0][1],cvd11tio[0][1],cvd13tio[0][1]],
##       [cvd3mgi[0][1],cvd5mgi[0][1],cvd7mgi[0][1],cvd9mgi[0][1],cvd11mgi[0][1],cvd13mgi[0][1]],'r--',
##       [cvd3tio[0][3],cvd5tio[0][3],cvd7tio[0][3],
##        cvd9tio[0][3],cvd11tio[0][3],cvd13tio[0][3]],
##       [cvd3mgi[0][3],cvd5mgi[0][3],cvd7mgi[0][3],
##        cvd9mgi[0][3],cvd11mgi[0][3],cvd13mgi[0][3]],'r--',
##       [cvd13tio[0][3],cvdafe2tio[0][3],cvdafe3tio[0][3]],
##       [cvd13mgi[0][3],cvdafe2mgi[0][3],cvdafe3mgi[0][3]],'b-')
##
##P.xlabel('TiO0.89', fontsize=22)
##
##P.ylabel('MgI0.88 [$\AA$]', fontsize=22)
##
##P.xticks(size=18)
##
##P.yticks(size=18)
##
##P.tick_params(axis='both', which='major', length=12, width=1.5, direction='in')
##P.tick_params(axis='both', which='minor', length=8, width=1.5, direction='in')
##P.annotate('13 Gyr', xy=(cvd13tio[0][1],cvd13mgi[0][1]),xycoords='data',xytext=(10,-30),textcoords='offset points',fontsize=14, color='r')
##P.annotate('11 Gyr', xy=(cvd11tio[0][1],cvd11mgi[0][1]),xycoords='data',xytext=(5,5),textcoords='offset points',fontsize=14, color='r')
##P.annotate('9 Gyr', xy=(cvd9tio[0][1],cvd9mgi[0][1]),xycoords='data',xytext=(10,0),textcoords='offset points',fontsize=14, color='r')
##P.annotate('7 Gyr', xy=(cvd7tio[0][1],cvd7mgi[0][1]),xycoords='data',xytext=(10,5),textcoords='offset points',fontsize=14, color='r')
##P.annotate('5 Gyr', xy=(cvd5tio[0][1],cvd5mgi[0][1]),xycoords='data',xytext=(15,0),textcoords='offset points',fontsize=14, color='r')
##P.annotate('3 Gyr', xy=(cvd3tio[0][1],cvd3mgi[0][1]),xycoords='data',xytext=(10,0),textcoords='offset points',fontsize=14, color='r')
##P.annotate('x = 3.5', xy=(cvd13tio[0][0],cvd13mgi[0][0]),xycoords='data',xytext=(-5,10),textcoords='offset points',fontsize=14, color='k')
##P.annotate('x = 3', xy=(cvd13tio[0][1],cvd13mgi[0][1]),xycoords='data',xytext=(0,10),textcoords='offset points',fontsize=14, color='k')
##P.annotate('x = 2.35', xy=(cvd13tio[0][2],cvd13mgi[0][2]),xycoords='data',xytext=(-25,12),textcoords='offset points',fontsize=14, color='k')
##P.annotate('x = Chabrier', xy=(cvd13tio[0][3],cvd13mgi[0][3]),xycoords='data',xytext=(-15,10),textcoords='offset points',fontsize=14, color='k')
##P.annotate('x = B-L', xy=(cvd13tio[0][4],cvd13mgi[0][4]),xycoords='data',xytext=(20,-10),textcoords='offset points',fontsize=14, color='k')
##P.annotate('[$\\alpha$/Fe]=+0.2', xy=(cvdafe2tio[0][3],cvdafe2mgi[0][3]),xycoords='data',xytext=(-40,10),textcoords='offset points',fontsize=14, color='b')
##P.annotate('[$\\alpha$/Fe]=+0.3', xy=(cvdafe3tio[0][3],cvdafe3mgi[0][3]),xycoords='data',xytext=(-40,10),textcoords='offset points',fontsize=14, color='b')
##
##
##P.show()
##
##
##
###Plot TiO0.89 v CaT
##P.figure(figsize=(11,9),dpi=72)
###Changing IMF at 13.5 Gyr
##P.scatter(cvd13tio[0], cvd13cat[0], marker='o',c='k',s=[90,70,50,30,10])
###Changing age for x=3 IMF
##P.scatter([cvd3tio[0][1],cvd5tio[0][1],cvd7tio[0][1],cvd9tio[0][1],cvd11tio[0][1],cvd13tio[0][1]],
##          [cvd3cat[0][1],cvd5cat[0][1],cvd7cat[0][1],cvd9cat[0][1],cvd11cat[0][1],cvd13cat[0][1]],
##          marker='s',c='r', s=[110,90,70,50,30,10])
###Changing age for Chabrier IMF
##P.scatter([cvd3tio[0][3],cvd5tio[0][3],cvd7tio[0][3],cvd9tio[0][3],cvd11tio[0][3],cvd13tio[0][3]],
##          [cvd3cat[0][3],cvd5cat[0][3],cvd7cat[0][3],cvd9cat[0][3],cvd11cat[0][3],cvd13cat[0][3]],
##          marker='s',c='r', s=[110,90,70,50,30,10])
###Changing [alpha/Fe] for 13.5 Gyr Chabrier IMF
##P.scatter([cvd13tio[0][3],cvdafe2tio[0][3],cvdafe3tio[0][3]],
##          [cvd13cat[0][3],cvdafe2cat[0][3],cvdafe3cat[0][3]],
##          marker='D',c='b',s=[10,30,50])
##
##P.plot(cvd13tio[0], cvd13cat[0], 'k-', [cvd3tio[0][1],cvd5tio[0][1],cvd7tio[0][1],
##                                  cvd9tio[0][1],cvd11tio[0][1],cvd13tio[0][1]],
##       [cvd3cat[0][1],cvd5cat[0][1],cvd7cat[0][1],cvd9cat[0][1],cvd11cat[0][1],cvd13cat[0][1]],'r--',
##       [cvd3tio[0][3],cvd5tio[0][3],cvd7tio[0][3],
##        cvd9tio[0][3],cvd11tio[0][3],cvd13tio[0][3]],
##       [cvd3cat[0][3],cvd5cat[0][3],cvd7cat[0][3],
##        cvd9cat[0][3],cvd11cat[0][3],cvd13cat[0][3]],'r--',
##       [cvd13tio[0][3],cvdafe2tio[0][3],cvdafe3tio[0][3]],
##       [cvd13cat[0][3],cvdafe2cat[0][3],cvdafe3cat[0][3]],'b-')
##
##P.xlabel('TiO0.89', fontsize=22)
##
##P.ylabel('CaT* [$\AA$]', fontsize=22)
##
##P.xticks(size=18)
##
##P.yticks(size=18)
##
##P.tick_params(axis='both', which='major', length=12, width=1.5, direction='in')
##P.tick_params(axis='both', which='minor', length=8, width=1.5, direction='in')
##P.annotate('13 Gyr', xy=(cvd13tio[0][1],cvd13cat[0][1]),xycoords='data',xytext=(10,-30),textcoords='offset points',fontsize=14, color='r')
##P.annotate('11 Gyr', xy=(cvd11tio[0][1],cvd11cat[0][1]),xycoords='data',xytext=(5,5),textcoords='offset points',fontsize=14, color='r')
##P.annotate('9 Gyr', xy=(cvd9tio[0][1],cvd9cat[0][1]),xycoords='data',xytext=(10,0),textcoords='offset points',fontsize=14, color='r')
##P.annotate('7 Gyr', xy=(cvd7tio[0][1],cvd7cat[0][1]),xycoords='data',xytext=(10,5),textcoords='offset points',fontsize=14, color='r')
##P.annotate('5 Gyr', xy=(cvd5tio[0][1],cvd5cat[0][1]),xycoords='data',xytext=(15,0),textcoords='offset points',fontsize=14, color='r')
##P.annotate('3 Gyr', xy=(cvd3tio[0][1],cvd3cat[0][1]),xycoords='data',xytext=(10,0),textcoords='offset points',fontsize=14, color='r')
##P.annotate('x = 3.5', xy=(cvd13tio[0][0],cvd13cat[0][0]),xycoords='data',xytext=(-5,10),textcoords='offset points',fontsize=14, color='k')
##P.annotate('x = 3', xy=(cvd13tio[0][1],cvd13cat[0][1]),xycoords='data',xytext=(0,10),textcoords='offset points',fontsize=14, color='k')
##P.annotate('x = 2.35', xy=(cvd13tio[0][2],cvd13cat[0][2]),xycoords='data',xytext=(-25,12),textcoords='offset points',fontsize=14, color='k')
##P.annotate('x = Chabrier', xy=(cvd13tio[0][3],cvd13cat[0][3]),xycoords='data',xytext=(-15,10),textcoords='offset points',fontsize=14, color='k')
##P.annotate('x = B-L', xy=(cvd13tio[0][4],cvd13cat[0][4]),xycoords='data',xytext=(20,-10),textcoords='offset points',fontsize=14, color='k')
##P.annotate('[$\\alpha$/Fe]=+0.2', xy=(cvdafe2tio[0][3],cvdafe2cat[0][3]),xycoords='data',xytext=(-40,10),textcoords='offset points',fontsize=14, color='b')
##P.annotate('[$\\alpha$/Fe]=+0.3', xy=(cvdafe3tio[0][3],cvdafe3cat[0][3]),xycoords='data',xytext=(-40,10),textcoords='offset points',fontsize=14, color='b')
##
##
##P.show()
##
##
###Plot NaI0.82 v CaT
##P.figure(figsize=(11,9),dpi=72)
###Changing IMF at 13.5 Gyr
##P.scatter(cvd13nai[0], cvd13cat[0], marker='o',c='k',s=[90,70,50,30,10])
###Changing age for x=3 IMF
##P.scatter([cvd3nai[0][1],cvd5nai[0][1],cvd7nai[0][1],cvd9nai[0][1],cvd11nai[0][1],cvd13nai[0][1]],
##          [cvd3cat[0][1],cvd5cat[0][1],cvd7cat[0][1],cvd9cat[0][1],cvd11cat[0][1],cvd13cat[0][1]],
##          marker='s',c='r', s=[110,90,70,50,30,10])
###Changing age for Chabrier IMF
##P.scatter([cvd3nai[0][3],cvd5nai[0][3],cvd7nai[0][3],cvd9nai[0][3],cvd11nai[0][3],cvd13nai[0][3]],
##          [cvd3cat[0][3],cvd5cat[0][3],cvd7cat[0][3],cvd9cat[0][3],cvd11cat[0][3],cvd13cat[0][3]],
##          marker='s',c='r', s=[110,90,70,50,30,10])
###Changing [alpha/Fe] for 13.5 Gyr Chabrier IMF
##P.scatter([cvd13nai[0][3],cvdafe2nai[0][3],cvdafe3nai[0][3]],
##          [cvd13cat[0][3],cvdafe2cat[0][3],cvdafe3cat[0][3]],
##          marker='D',c='b',s=[10,30,50])
##
##P.plot(cvd13nai[0], cvd13cat[0], 'k-', [cvd3nai[0][1],cvd5nai[0][1],cvd7nai[0][1],
##                                  cvd9nai[0][1],cvd11nai[0][1],cvd13nai[0][1]],
##       [cvd3cat[0][1],cvd5cat[0][1],cvd7cat[0][1],cvd9cat[0][1],cvd11cat[0][1],cvd13cat[0][1]],'r--',
##       [cvd3nai[0][3],cvd5nai[0][3],cvd7nai[0][3],
##        cvd9nai[0][3],cvd11nai[0][3],cvd13nai[0][3]],
##       [cvd3cat[0][3],cvd5cat[0][3],cvd7cat[0][3],
##        cvd9cat[0][3],cvd11cat[0][3],cvd13cat[0][3]],'r--',
##       [cvd13nai[0][3],cvdafe2nai[0][3],cvdafe3nai[0][3]],
##       [cvd13cat[0][3],cvdafe2cat[0][3],cvdafe3cat[0][3]],'b-')
##
##P.xlabel('NaI0.82 [$\AA$]', fontsize=22)
##
##P.ylabel('CaT* [$\AA$]', fontsize=22)
##
##P.xticks(size=18)
##
##P.yticks(size=18)
##
##P.tick_params(axis='both', which='major', length=12, width=1.5, direction='in')
##P.tick_params(axis='both', which='minor', length=8, width=1.5, direction='in')
##P.annotate('13 Gyr', xy=(cvd13nai[0][1],cvd13cat[0][1]),xycoords='data',xytext=(10,-30),textcoords='offset points',fontsize=14, color='r')
##P.annotate('11 Gyr', xy=(cvd11nai[0][1],cvd11cat[0][1]),xycoords='data',xytext=(5,5),textcoords='offset points',fontsize=14, color='r')
##P.annotate('9 Gyr', xy=(cvd9nai[0][1],cvd9cat[0][1]),xycoords='data',xytext=(10,0),textcoords='offset points',fontsize=14, color='r')
##P.annotate('7 Gyr', xy=(cvd7nai[0][1],cvd7cat[0][1]),xycoords='data',xytext=(10,5),textcoords='offset points',fontsize=14, color='r')
##P.annotate('5 Gyr', xy=(cvd5nai[0][1],cvd5cat[0][1]),xycoords='data',xytext=(15,0),textcoords='offset points',fontsize=14, color='r')
##P.annotate('3 Gyr', xy=(cvd3nai[0][1],cvd3cat[0][1]),xycoords='data',xytext=(10,0),textcoords='offset points',fontsize=14, color='r')
##P.annotate('x = 3.5', xy=(cvd13nai[0][0],cvd13cat[0][0]),xycoords='data',xytext=(-5,10),textcoords='offset points',fontsize=14, color='k')
##P.annotate('x = 3', xy=(cvd13nai[0][1],cvd13cat[0][1]),xycoords='data',xytext=(0,10),textcoords='offset points',fontsize=14, color='k')
##P.annotate('x = 2.35', xy=(cvd13nai[0][2],cvd13cat[0][2]),xycoords='data',xytext=(-25,12),textcoords='offset points',fontsize=14, color='k')
##P.annotate('x = Chabrier', xy=(cvd13nai[0][3],cvd13cat[0][3]),xycoords='data',xytext=(-15,10),textcoords='offset points',fontsize=14, color='k')
##P.annotate('x = B-L', xy=(cvd13nai[0][4],cvd13cat[0][4]),xycoords='data',xytext=(20,-10),textcoords='offset points',fontsize=14, color='k')
##P.annotate('[$\\alpha$/Fe]=+0.2', xy=(cvdafe2nai[0][3],cvdafe2cat[0][3]),xycoords='data',xytext=(-40,10),textcoords='offset points',fontsize=14, color='b')
##P.annotate('[$\\alpha$/Fe]=+0.3', xy=(cvdafe3nai[0][3],cvdafe3cat[0][3]),xycoords='data',xytext=(-40,10),textcoords='offset points',fontsize=14, color='b')
##
##
##P.show()



###Plot NaI0.82 v NaD0.59
##P.rcParams['font.family'] = 'Times New Roman'
##P.figure(figsize=(11,9),dpi=72)
###Changing IMF at 13.5 Gyr
##P.scatter(cvd13nai[0], cvd13nad[0], marker='o',c='k',s=[90,70,50,30,10])
###Changing age for x=3 IMF
##P.scatter([cvd3nai[0][1],cvd5nai[0][1],cvd7nai[0][1],cvd9nai[0][1],cvd11nai[0][1],cvd13nai[0][1]],
##          [cvd3nad[0][1],cvd5nad[0][1],cvd7nad[0][1],cvd9nad[0][1],cvd11nad[0][1],cvd13nad[0][1]],
##          marker='s',c='r', s=[110,90,70,50,30,10])
###Changing age for Chabrier IMF
##P.scatter([cvd3nai[0][3],cvd5nai[0][3],cvd7nai[0][3],cvd9nai[0][3],cvd11nai[0][3],cvd13nai[0][3]],
##          [cvd3nad[0][3],cvd5nad[0][3],cvd7nad[0][3],cvd9nad[0][3],cvd11nad[0][3],cvd13nad[0][3]],
##          marker='s',c='r', s=[110,90,70,50,30,10])
###Changing [alpha/Fe] for 13.5 Gyr Chabrier IMF
##P.scatter([cvd13nai[0][3],cvdafe2nai[0][3],cvdafe3nai[0][3]],
##          [cvd13nad[0][3],cvdafe2nad[0][3],cvdafe3nad[0][3]],
##          marker='D',c='b',s=[10,30,50])
###Varying [Na/Fe] for 13.5 Gyr Chabrier IMF
##P.scatter([cvdvarnai[0][0],cvd13nai[0][3],cvdvarnai[0][1]],
##          [cvdvarnad[0][0],cvd13nad[0][3],cvdvarnad[0][1]],
##          marker='D',c='g',s=[50,5,50])
##
##P.plot(cvd13nai[0], cvd13nad[0], 'k-', [cvd3nai[0][1],cvd5nai[0][1],cvd7nai[0][1],
##                                  cvd9nai[0][1],cvd11nai[0][1],cvd13nai[0][1]],
##       [cvd3nad[0][1],cvd5nad[0][1],cvd7nad[0][1],cvd9nad[0][1],cvd11nad[0][1],cvd13nad[0][1]],'r--',
##       [cvd3nai[0][3],cvd5nai[0][3],cvd7nai[0][3],
##        cvd9nai[0][3],cvd11nai[0][3],cvd13nai[0][3]],
##       [cvd3nad[0][3],cvd5nad[0][3],cvd7nad[0][3],
##        cvd9nad[0][3],cvd11nad[0][3],cvd13nad[0][3]],'r--',
##       [cvd13nai[0][3],cvdafe2nai[0][3],cvdafe3nai[0][3]],
##       [cvd13nad[0][3],cvdafe2nad[0][3],cvdafe3nad[0][3]],'b-',
##       [cvdvarnai[0][0],cvd13nai[0][3],cvdvarnai[0][1]],
##       [cvdvarnad[0][0],cvd13nad[0][3],cvdvarnad[0][1]],'g--')
##
##P.xlabel('NaI0.82 [$\AA$]', fontsize=24)
##P.ylabel('NaD0.59 [$\AA$]', fontsize=24)
##P.xticks(size=22)
##P.yticks(size=22)
##P.tick_params(axis='both', which='major', length=12, width=1.5, direction='in')
##P.tick_params(axis='both', which='minor', length=8, width=1.5, direction='in')
##P.xlim([0.0, 1.2])
##
##P.annotate('13', xy=(cvd13nai[0][1],cvd13nad[0][1]),xycoords='data',xytext=(10,-10),textcoords='offset points',fontsize=14, color='r')
##P.annotate('11', xy=(cvd11nai[0][1],cvd11nad[0][1]),xycoords='data',xytext=(5,0),textcoords='offset points',fontsize=14, color='r')
##P.annotate('9', xy=(cvd9nai[0][1],cvd9nad[0][1]),xycoords='data',xytext=(10,-5),textcoords='offset points',fontsize=14, color='r')
##P.annotate('7', xy=(cvd7nai[0][1],cvd7nad[0][1]),xycoords='data',xytext=(10,0),textcoords='offset points',fontsize=14, color='r')
##P.annotate('5', xy=(cvd5nai[0][1],cvd5nad[0][1]),xycoords='data',xytext=(15,-5),textcoords='offset points',fontsize=14, color='r')
##P.annotate('3 Gyr', xy=(cvd3nai[0][1],cvd3nad[0][1]),xycoords='data',xytext=(10,0),textcoords='offset points',fontsize=14, color='r')
##P.annotate('x = 3.5', xy=(cvd13nai[0][0],cvd13nad[0][0]),xycoords='data',xytext=(-20,10),textcoords='offset points',fontsize=14, color='k')
##P.annotate('x = 3', xy=(cvd13nai[0][1],cvd13nad[0][1]),xycoords='data',xytext=(-10,10),textcoords='offset points',fontsize=14, color='k')
##P.annotate('x = 2.35', xy=(cvd13nai[0][2],cvd13nad[0][2]),xycoords='data',xytext=(-25,10),textcoords='offset points',fontsize=14, color='k')
##P.annotate('Chab', xy=(cvd13nai[0][3],cvd13nad[0][3]),xycoords='data',xytext=(-30,5),textcoords='offset points',fontsize=14, color='k')
##P.annotate('B-L', xy=(cvd13nai[0][4],cvd13nad[0][4]),xycoords='data',xytext=(-20,10),textcoords='offset points',fontsize=14, color='k')
##P.annotate('[$\\alpha$/Fe]=+0.2', xy=(cvdafe2nai[0][3],cvdafe2nad[0][3]),xycoords='data',xytext=(-75,-5),textcoords='offset points',fontsize=14, color='b')
##P.annotate('[$\\alpha$/Fe]=+0.3', xy=(cvdafe3nai[0][3],cvdafe3nad[0][3]),xycoords='data',xytext=(-45,-42),textcoords='offset points',fontsize=14, color='b')
##P.annotate('[Na/Fe]=+0.3', xy=(cvdvarnai[0][0],cvdvarnad[0][0]),xycoords='data',xytext=(-30,10),textcoords='offset points',fontsize=14, color='g')
##P.annotate('[Na/Fe]=-0.3', xy=(cvdvarnai[0][1],cvdvarnad[0][1]),xycoords='data',xytext=(-40,-20),textcoords='offset points',fontsize=14, color='g')
##P.show()
