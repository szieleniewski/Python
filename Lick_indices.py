'''Code to handle Lick index measurements from Sanchez-Blazquez et al. 2006b and Loubser et al. 2009.

Author: Simon Zieleniewski

Last updated: 18-12-15

'''

import numpy as n
import pylab as P
import astropy.io.fits as p

def mehlert(ind, gal, sigma):
    '''Code to return value for indices
    from Mehlert et al. 2000 for given
    velocity disperison value

    Inputs:

        - gal: Galaxy name: GMP2921, GMP3329, GMP4928
        - ind: index [Mgb, Fe5270, Fe5335, <Fe>, NaD]
        - sigma: vel. disp. [km/s]

    Outputs:

        - index: index vlue [A]
        - err: error on index [A]

    '''
    params = {'Hb': [1.0, 7.444E-5, 2.667E-7, 1.136E-9],
              'Mgb': [1.0, -9.333E-5, 2.800E-6, -1.481E-9],
              'Fe5270': [1.0, 4.000E-5, 2.667E-6, -1.481E-9],
              'Fe5335': [1.0, -5.667E-5, 5.444E-6, 2.963E-10],
              'NaD': [1.0, 5.222E-5, 2.000E-7, 1.975E-9]}
    indvals = n.genfromtxt('./Loubser2009_line_indices.tsv', delimiter='; ')
    coeffs = params[ind]
    C = coeffs[0] + coeffs[1]*sigma + coeffs[2]*sigma**2 + coeffs[3]*sigma**3




def calc_index(author, ind, gal, sigma):
    '''Code to return value for <Fe> Index
    from Loubser et al. 2009 for given
    velocity disperison value

    Inputs:

        - author: [loubser, sanchez]
        - gal: Galaxy name: GMP2921, GMP3329, GMP4928
        - ind: index [Mgb, Fe5270, Fe5335, <Fe>, NaD]
        - sigma: vel. disp. [km/s]

    Outputs:

        - index: index vlue [A]
        - err: error on index [A]

    '''

    sigL = 200.0
    sigma = n.float(sigma)

    if author == 'sanchez':
        params = {'Hb': [0.9907, 0.0056, 0.0037],
                  'Mgb': [0.9645, -0.0749, 0.1104],
                  'Fe5270': [0.8253, 0.1228, 0.0518],
                  'Fe5335': [0.8432, -0.0814, 0.2382]}
        indvals = n.genfromtxt('./Sanchez-Blazquez2006b_line_indices.tsv', delimiter='; ')
        coeffs = params[ind]
        C = coeffs[0] + coeffs[1]*(sigL/sigma) + coeffs[2]*(sigL/sigma)**2

    elif author == 'loubser':
        params = {'Hb': [1.0, 4.553E-5, 2.4319E-9, 1.096E-9],
                  'Mgb': [1.0, -4.1861E-5, 2.02E-6, 4.1426E-11],
                  'Fe5270': [1.0, -8.9369E-6, 2.8442E-6, -2.4675E-9],
                  'Fe5335': [1.0, 2.0324E-5, 4.199E-6, 1.1364E-9],
                  'NaD': [1.0, 3.003E-5, 3.4807E-7, 9.2784E-10]}
        indvals = n.genfromtxt('./Loubser2009_line_indices.tsv', delimiter='; ')
        coeffs = params[ind]
        C = coeffs[0] + coeffs[1]*sigma + coeffs[2]*sigma**2 + coeffs[3]*sigma**3

    #load Lick index values [Hbeta;e_Hbeta;Mgb;e_Mgb;Fe5270;e_Fe5270;Fe5335;e_Fe5335;NaD;e_NaD;<Fe>;e_<Fe>;[MgFe]';e_[MgFe]']
    galpos = {'GMP2921':2, 'GMP3329':1, 'GMP4928':0}
    pos = {'Hb':0, 'Mgb':2, 'Fe5270':4, 'Fe5335':6, 'NaD':8}

    index = indvals[galpos[gal],pos[ind]]
    inderr = indvals[galpos[gal],pos[ind]+1]

    newind = index/C
    newerr = inderr/C

    print author, ' index = ', n.round(newind,3), u'\xb1', n.round(newerr,3), ' A'
    return newind, newerr


def Fe(fe52, fe52_e, fe53, fe53_e):
    ind = 0.5*(fe52+fe53)
    indsig = n.sqrt((0.5*fe52_e)**2+(0.5*fe53_e)**2)
    return ind, indsig


def mgfe(mgb, mgb_e, fe52, fe52_e, fe53, fe53_e):
    ind = n.sqrt(mgb*(0.72*fe52 + 0.28*fe53))
    indsig = ind*(0.5*(n.sqrt((mgb_e/mgb)**2 + (n.sqrt((0.72*fe52_e)**2 + (0.28*fe53_e)**2)/(0.72*fe52+0.28*fe53))**2)))
    return ind, indsig


def Fe_index(author, gal, sigma):

    fe5270 = calc_index(author, 'Fe5270', gal, sigma)
    fe5335 = calc_index(author, 'Fe5335', gal, sigma)

    fe_index = Fe(fe5270[0], fe5270[1], fe5335[0], fe5335[1])

    print '<Fe> = ', n.round(fe_index[0],3), u'\xb1', n.round(fe_index[1],3), ' A'
    return fe_index[0], fe_index[1]


def MgFe_index(author, gal, sigma):

    Mgb = calc_index(author, 'Mgb', gal, sigma)
    fe5270 = calc_index(author, 'Fe5270', gal, sigma)
    fe5335 = calc_index(author, 'Fe5335', gal, sigma)

    mgfe_index = mgfe(Mgb[0], Mgb[1], fe5270[0], fe5270[1], fe5335[0], fe5335[1])

    print "[MgFe]' = ", n.round(mgfe_index[0],3), u'\xb1', n.round(mgfe_index[1],3), ' A'
    return mgfe_index[0], mgfe_index[1]



if __name__=="__main__":
    import sys

    names = ['Fe', 'MgFe']
    funcs = [Fe_index, MgFe_index]
    print len(sys.argv)
    if len(sys.argv) <= 1 or sys.argv[1] not in names:
        print ''
        print '1.  index [Fe, MgFe], author [sanchez, loubser], gal [GMP2921, GMP3329, GMP4928], sigma [km/s]'
        print ''
        sys.exit()

    if sys.argv[1] == names[0]:
        funcs[0](sys.argv[2], sys.argv[3], sys.argv[4])
    if sys.argv[1] == names[1]:
        funcs[1](sys.argv[2], sys.argv[3], sys.argv[4])
