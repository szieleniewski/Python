''' Code to experiment with different IMF forms

Author: Simon Zieleniewski

Last Updated: 10-03-16

'''

import numpy as n
import scipy.interpolate as si
import pylab as P


####FONTS
tfsize = 28
lfsize = 24
gensize = 20
majtsize = 14
mintsize = 6
twidth = 2
msize = 10
####


def log_interp1d(xx, yy, kind='linear'):
    logx = n.log10(xx)
    logy = n.log10(yy)
    lin_interp = si.interp1d(logx, logy)
    log_interp = lambda zz: n.power(10.0, lin_interp(n.log10(zz)))
    return log_interp


#IMF defined as number of stars per unit mass per unit volume
#formed at one time.


def unimodal(masses, mu):
    imf = masses**(-(mu+1))
    #Normalise
    # beta = 1./((masses[-1]**(-x+1)/(-x+1))-(masses[0]**(-x+1)/(-x+1)))
    beta = n.trapz(imf, x=masses)

    logimf = imf*masses
    logbeta = n.trapz(logimf, x=masses)

    return n.column_stack((masses, imf/beta, logimf/logbeta))


def bimodal(masses, mu):
    m0 = 0.2
    m1 = 0.4
    m2 = 0.6
    cut1 = n.where(masses > m0)[0][0]
    cut2 = n.where(masses > m2)[0][0]

    finalimf = n.zeros(len(masses),dtype=float)

    #Less than 0.2
    imf1 = m1**(-mu)/masses[:cut1]
    finalimf[:cut1] = imf1

    #above 0.6
    imf3 = masses[cut2:]**(-(mu+1))
    finalimf[cut2:] = imf3

    arr = n.concatenate((masses[:cut1], masses[cut2:]))
    imfs = n.concatenate((finalimf[:cut1], finalimf[cut2:]))
    # interp = si.interp1d(arr, imfs, 1)
    interp = log_interp1d(arr, imfs, 6)
    imf2 = interp(masses[cut1:cut2])

    finalimf[cut1:cut2] = imf2

    #Normalise
    # beta = 1/(0.4**(-x+1)*(n.log(masses[cut1])-n.log(masses[0]))
    #           +((masses[-1]**(-x+1)/(-x+1))-(masses[cut2]**(-x+1)/(-x+1)))
    #           +n.trapz(finalimf[cut1:cut2], x=masses[cut1:cut2]))
    beta = n.trapz(finalimf, x=masses)

    logimf = finalimf*masses
    logbeta = n.trapz(logimf, x=masses)

    return n.column_stack((masses, finalimf/beta, logimf/logbeta))


def kroupa(masses):
    m0 = 0.08
    m1 = 0.5
    cut1 = n.where(masses > m0)[0][0]
    cut2 = n.where(masses > m1)[0][0]
    finalimf = n.zeros(len(masses),dtype=float)

    #Less than 0.08
    imf1 = (masses[:cut1]/m0)**(-0.3)
    finalimf[:cut1] = imf1

    #Between 0.08 and 0.5
    imf2 = (masses[cut1:cut2]/m0)**(-1.3)
    finalimf[cut1:cut2] = imf2

    #Above 0.5
    imf3 = (m1/m0)**(-1.3)*(masses[cut2:]/m1)**(-2.3)
    finalimf[cut2:] = imf3

    #Normalise
    # beta = 1/((0.6696*(masses[cut1]**0.7-masses[0]**0.7))
    #           +(-0.125*(masses[cut2]**-0.3-masses[cut1]**-0.3))
    #           +(-0.01442*(masses[-1]**-1.3 - masses[cut2]**-1.3)))
    beta = n.trapz(finalimf, x=masses)

    logimf = finalimf*masses
    logbeta = n.trapz(logimf, x=masses)

    return n.column_stack((masses, finalimf/beta, logimf/logbeta))


# def krouparev(masses):
#     cut1 = n.where(masses > 0.08)[0][0]
#     cut2 = n.where(masses > 0.5)[0][0]
#     cut3 = n.where(masses > 1.0)[0][0]
#     finalimf = n.zeros(len(masses),dtype=float)
#
#     #Less than 0.08
#     imf1 = (masses[:cut1]/0.08)**(-0.3)
#     finalimf[:cut1] = imf1
#
#     #Between 0.08 and 0.5
#     imf2 = (masses[cut1:cut2]/0.08)**(-1.8)
#     finalimf[cut1:cut2] = imf2
#
#     #Between 0.5 and 1.0
#     imf3 = (0.5/0.08)**(-1.8)*(masses[cut2:cut3]/0.5)**(-2.7)
#     finalimf[cut2:cut3] = imf3
#
#     #Above 1.0
#     imf4 = (0.5/0.08)**(-1.8)*(1.0/0.5)**(-2.7)*(masses[cut3:]/1.0)**(-2.3)
#     finalimf[cut3:] = imf4
#
#     #Normalise
#     beta = 1/((0.6696*(masses[cut1]**0.7-masses[0]**0.7))
#               +(-0.01326*(masses[cut2]**-0.8-masses[cut1]**-0.8))
#               +(-0.00334*(masses[cut3]**-1.7-masses[cut2]**-1.7))
#               +(-0.00437*(masses[-1]**-1.3 - masses[cut2]**-1.3)))
#
#     return n.column_stack((masses, beta*finalimf))


def chabrier(masses):
    cut1 = n.where(masses > 1.0)[0][0]
    finalimf = n.zeros(len(masses),dtype=float)
    #Below 1.0
    imf1 = 0.086*(1/masses[:cut1])*n.exp(-(n.log10(masses[:cut1])-n.log10(0.22))**2/(2.*.57**2))
    finalimf[:cut1] = imf1
    #Above 1.0
    imf2 = 4.43E-2*masses[cut1:]**-2.3
    finalimf[cut1:] = imf2

    #Normalise
    beta = n.trapz(finalimf, x=masses)

    logimf = finalimf*masses
    logbeta = n.trapz(logimf, x=masses)

    return n.column_stack((masses, finalimf/beta, logimf/logbeta))


if __name__=="__main__":

    lcut = 0.01
    hcut = 120.0
    #create array of stellar masses
    masses = n.logspace(n.log10(lcut), n.log10(hcut), 100000)

    salp = 1.3
    bh = 2.0
    un2p35 = unimodal(masses, salp)
    un3 = unimodal(masses, bh)
    bm2p35 = bimodal(masses, salp)
    bm3 = bimodal(masses, bh)
    kp = kroupa(masses)
    ch = chabrier(masses)

    print n.sum(un2p35[:,0]*un2p35[:,1])
    print n.sum(un3[:,0]*un3[:,1])
    print n.sum(bm2p35[:,1])
    print n.sum(bm3[:,1])
    print n.sum(kp[:,1])
    print n.sum(ch[:,1])

    # un2p35[:,1] /=n.sum(un2p35[:,1])
    # un3[:,1] /=n.sum(un3[:,1])
    # bm2p35[:,1] /=n.sum(bm2p35[:,1])
    # bm3[:,1] /=n.sum(bm3[:,1])
    # kp[:,1] /=n.sum(kp[:,1])
    # ch[:,1] /=n.sum(ch[:,1])

    #Plot
    P.rcParams['font.family'] = 'Times New Roman'
    P.rcParams.update({'font.size': gensize})
    P.figure(figsize=(11,9), dpi=72)
    P.subplots_adjust(left=0.14, bottom=0.12)

    p1, = P.plot(n.log10(un2p35[:,0]), n.log10(un2p35[:,1]), 'k-', lw=2.0)
    p2, = P.plot(n.log10(un3[:,0]), n.log10(un3[:,1]), 'k--', lw=2.0)
    p3, = P.plot(n.log10(bm2p35[:,0]), n.log10(bm2p35[:,1]), 'r-', lw=2.0)
    p4, = P.plot(n.log10(bm3[:,0]), n.log10(bm3[:,1]), 'r--', lw=2.0)
    p5, = P.plot(n.log10(kp[:,0]), n.log10(kp[:,1]), 'b-', lw=2.0, alpha=0.6)
    p6, = P.plot(n.log10(ch[:,0]), n.log10(ch[:,1]), 'm-', lw=2.0, alpha=0.6)

    P.xlabel('log$_{10}(m)$', fontsize=tfsize)
    P.ylabel('log$_{10}$[$\\xi(m)$]', fontsize=tfsize)
    P.xticks(size=lfsize)
    P.yticks(size=lfsize)
    P.minorticks_on()
    P.tick_params(axis='both', direction='in', length=majtsize,
                  width=1.5, which='major')
    P.tick_params(axis='both', direction='in', length=mintsize,
                  width=1.5, which='minor')

    leg = P.legend([p1,p2,p3,p4,p5,p6],
                     ['Salp (x='+str(salp+1)+')', 'Un (x='+str(bh+1)+')', 'Bi (x='+str(salp+1)+')', 'Bi (x='+str(bh+1)+')', 'Kroupa01', 'Chabrier03'],
                     loc='upper right', frameon=False, borderaxespad=0.5, labelspacing=0.5, handlelength=2.2)
    # leg = P.legend([p1,p2,p3,p4,p5],
    #                  ['Salp (x='+str(salp)+')', 'Un (x='+str(bh)+')', 'Bi (x='+str(salp)+')', 'Bi (x='+str(bh)+')', 'Kroupa01'],
    #                  loc='upper right', frameon=False, borderaxespad=0.5, labelspacing=0.5, handlelength=2.2)
    ltext = leg.get_texts()

    P.setp(ltext[0], fontsize=20, color='k')
    P.setp(ltext[1], fontsize=20, color='k')
    P.setp(ltext[2], fontsize=20, color='r')
    P.setp(ltext[3], fontsize=20, color='r')
    P.setp(ltext[4], fontsize=20, color='b')
    P.setp(ltext[5], fontsize=20, color='m')

    P.xlim([-2, 2])
    P.ylim([-6,1.5])

    P.tight_layout()
    P.show()

    dat = n.genfromtxt('/Users/zieleniewski/Desktop/Scalo1986_data_26-02-16.txt', delimiter=',')

    #Plot
    P.rcParams['font.family'] = 'Times New Roman'
    P.rcParams.update({'font.size': gensize})
    P.figure(figsize=(11,9), dpi=72)
    P.subplots_adjust(left=0.14, bottom=0.12)

    p1, = P.plot(n.log10(un2p35[:,0]), n.log10(un2p35[:,2])+0.4, 'k-', lw=2.0)
    # p2, = P.plot(n.log10(un3[:,0]), n.log10(un3[:,2]), 'k--', lw=2.0)
    # p3, = P.plot(n.log10(bm2p35[:,0]), n.log10(bm2p35[:,2]), 'r-', lw=2.0)
    # p4, = P.plot(n.log10(bm3[:,0]), n.log10(bm3[:,2]), 'r--', lw=0.0)
    p5, = P.plot(n.log10(kp[:,0]), n.log10(kp[:,2]), 'r-', lw=2.0, alpha=0.6)
    p6, = P.plot(n.log10(ch[:,0]), n.log10(ch[:,2]), 'b-', lw=2.0, alpha=0.6)

    p7 = P.errorbar(dat[:,0], dat[:,1]-2, yerr=dat[:,2], ls='', fmt='ko', mec='k',
    capthick=2, elinewidth=1.5)

    P.xlabel('log$_{10}$(m)', fontsize=tfsize)
    P.ylabel('log$_{10}$[$\\xi$(log$_{10}$(m))]', fontsize=tfsize)
    P.xticks(size=lfsize)
    P.yticks(size=lfsize)
    P.minorticks_on()
    P.tick_params(axis='both', direction='in', length=majtsize,
                  width=1.5, which='major')
    P.tick_params(axis='both', direction='in', length=mintsize,
                  width=1.5, which='minor')

    leg = P.legend([p1,p5,p6,p7],
                     ['Salpeter (x='+str(salp+1)+')', 'Kroupa01', 'Chabrier03', 'Scalo1986'],
                     loc='lower center', frameon=False, borderaxespad=0.6, labelspacing=0.5, handlelength=2.2, numpoints=1)
    # leg = P.legend([p1,p2,p3,p4,p5],
    #                  ['Salp (x='+str(salp)+')', 'Un (x='+str(bh)+')', 'Bi (x='+str(salp)+')', 'Bi (x='+str(bh)+')', 'Kroupa01'],
    #                  loc='upper right', frameon=False, borderaxespad=0.5, labelspacing=0.5, handlelength=2.2)
    ltext = leg.get_texts()

    P.setp(ltext[0], fontsize=20, color='k')
    P.setp(ltext[1], fontsize=20, color='r')
    P.setp(ltext[2], fontsize=20, color='b')
    P.setp(ltext[3], fontsize=20, color='k')

    P.xlim([-1.1, 0.6])
    P.ylim([-3,1.1])

    P.tight_layout()
    P.show()
