'''Program to compute line indices from stellar
population model SEDs.
First attempt: Read in SEDs from V03, plot spectra
and calculate CaT* index.
Written by Simon Zieleniewski

Last updated: 20-04-16
'''


from scipy import interpolate, optimize
import scipy as s
import astropy.io.fits as p
import numpy as n
import math as m
import pylab as P
import os
import sys
sys.path.append('~/Documents/Oxford/D.Phil/Python/Ryan_code/')
import spectools as t

L_sun = 3.826E33 # the L_sun in erg/s defined by Vazdekis et al. (2012)

class cattoolset:

    def __init__(self):
#VACUUM WAVELENGTHS
        #Bandpasses for NaD (angstroms)
        self.nad = n.array([5878.5, 5911.0])
        self.nadcont = n.array([5862.2, 5877.2, 5923.7, 5949.7])
        #Band passes for CaT* (angstroms)
        self.ca = n.array((8484., 8513., 8522., 8562., 8642., 8682.))
        self.pa = n.array((8461., 8474., 8577., 8619., 8730., 8772.))
        self.cacont = n.array((8474.,8484.,8563.,8577.,8619.,
                        8642.,8700.,8725.,8776.,8792.))
        #Bandpasses for MgI, NaI, FeH (angstroms)
        self.na = n.array((8177., 8205.))
        self.nacont = n.array((8170., 8177., 8205., 8215.))
        self.mg = n.array((8801.9, 8816.9))
        self.mgcont = n.array((8777.4, 8789.4, 8847.4, 8857.4))
        self.feh = n.array((9905., 9935.))
        self.fehcont = n.array((9855., 9880., 9940., 9970.))
        self.tio = n.array((0., 0.))
        self.tiocont = n.array((8835., 8855., 8870., 8890.))

        #Bandpasses for H-beta, Mgb, Fe5270, Fe5335 (angstroms)
#        self.hb = n.array((4847.875, 4876.625))
#        self.hbcont = n.array((4827.875, 4847.875, 4876.625, 4891.625))
        self.mgb = n.array((5160.125, 5192.625))
        self.mgbcont = n.array((5142.625, 5161.375, 5191.375, 5206.375))
        self.fe52 = n.array((5245.650, 5285.650))
        self.fe52cont = n.array((5233.150, 5248.150, 5285.650, 5318.150))
        self.fe53 = n.array((5312.125, 5352.125))
        self.fe53cont = n.array((5304.625, 5315.875, 5353.375, 5363.375))

        #Sausage H-alpha index band passes (feature band pass from Groves et al. 2011, continuum=sausage)
        self.ha = n.array([6553.,6573.])
#        self.hacont = n.array([6500.,6520.,6600.,6620.])
        self.hacont = n.array([6546.,6551.,6575.,6580.])
        #Sausage H-beta index band pass (same as H-alpha)
        self.hb = n.array([4851.,4871.])
#        self.hbcont = n.array([4835.,4840.,4900.,4920.])
        self.hbcont = n.array([4844.,4849.,4873.,4878.])

#AIR WAVELENGTHS
        #Bandpasses for NaD (angstroms)
        self.anad = n.array([5876.87110607978, 5909.362402982217])
        self.anadcont = n.array([5860.575470321278, 5875.571454165305, 5922.059001582479, 5948.052037227286])
        #Band passes for CaT* (angstroms)
        self.aca = n.array((8481.669445901154, 8510.661605850613, 8519.659172639998, 8559.648357851342, 8639.626725770459, 8679.615908502386))
        self.apa = n.array((8458.67566354994, 8471.672149261918, 8574.644302088736, 8616.632945331574, 8727.60292672214, 8769.591566731815))
        self.acacont = n.array((8471.672149261918, 8481.669445901154,8560.64808747083,8574.644302088736,8616.632945331574,
                        8639.626725770459, 8697.611040469177, 8722.60427904426, 8773.590484783086, 8789.58615691068))
        #Bandpasses for MgI, NaI, FeH (angstroms)
        self.ana = n.array((8174.752413483058, 8202.744848675571))
        self.anacont = n.array((8167.75430461091, 8174.752413483058, 8202.744848675571, 8212.742146844472))
        self.amg = n.array((8799.483478977758, 8814.479421413944))
        self.amgcont = n.array((8774.990106099198, 8786.986860198365, 8844.971170703131, 8854.968465456062))
        self.afeh = n.array((9902.284831158933, 9932.27670274))
        self.afehcont = n.array((9852.298377857542, 9877.291604612945, 9437.410781562552, 9967.26721921015))
        self.atio = n.array((0., 0.))
        self.atiocont = n.array((8832.574525143811, 8852.56911471966, 8867.565056777688, 8887.559646024309))
        #Initialise dispersion [A/pixel]
        self.disp = None

    def getspec(self, filename):
        '''Function to read SED file and return array of wavelength and SED.
Returns array of form (wavelength, data).'''

        dat, head = p.getdata(filename, 0, header=True)
        wavel = n.linspace(head['CRVAL1'], head['CRVAL1']+(len(dat)-1)\
                           *head['CDELT1'],len(dat))
        #Dispersion value (wavelength/pixel)
        self.disp = head['CDELT1']
        #Wavelength correction for SWIFT data (which is in microns)
        try:
            if head['INSTRUME'] == 'SWIFT':
                wavel *= 10000. #Puts lambda in angstroms
                self.disp *= 10000. #Puts dispersion in angstroms/pixel
        except:
            pass
        SED = n.column_stack((wavel,dat))

        return SED


    def specplot(self, SED):
        '''Function to plot SED'''
        P.plot(SED[:,0],SED[:,1],'k-')
        #For units, 0 = Lsun, 1 = erg/s/cm^2/Ang
        P.xlabel('Wavelength [$\AA$]')
        P.ylabel('SED [Lsun]')
        P.show()



    def irindex(self, SED, disp=1.0, index='MgI', var_SED=None, vac_or_air='vac'):
        '''Function to calculate several near-IR line indices.
Function takes as argument an SED of form (wavlength col, data col)
and a choice of line index (NaI, MgI, FeH).

        Inputs:

            SED: spectral energy distribution of the form (lambda column, flux column)
            disp: float - spectral dispersion [A/pixel]
            index: Choice of 'MgI', 'NaI', 'FeH', 'CaT', 'PaT', 'Ha', 'Hb', NaD'.
                   Sets the index measurement.
            var_SED: variance spectrum corresponding to SED of the same form as SED.
            vac_or_air: vacuum or air wavelengths for index definitions.

        Outputs:

            ind: value of chosen index
            ind_var_tot: variance on index measurement.
        '''

        if vac_or_air == 'vac':
            if index=='MgI':
                cont = self.mgcont; band = self.mg
            elif index=='NaI':
                cont = self.nacont; band = self.na
            elif index=='FeH':
                cont = self.fehcont; band = self.feh
            elif index=='CaT':
                cont = self.cacont; band = self.ca
            elif index=='PaT':
                cont = self.cacont; band = self.pa
            elif index=='Ha':
                cont = self.hacont; band = self.ha
##            elif index=='Hb':
##                cont = self.hbcont; band = self.hb
            elif index=='NaD':
                cont = self.nadcont; band = self.nad
            elif index=='Mgb':
                cont = self.mgbcont; band = self.mgb
            elif index=='Hb':
                cont = self.hbcont; band = self.hb
            elif index=='Fe52':
                cont = self.fe52cont; band = self.fe52
            elif index=='Fe53':
                cont = self.fe53cont; band = self.fe53
            else:
                raise ValueError('Please choose one of the following indices:'+
                                 'MgI, NaI, FeH, CaT, PaT, NaD')
        elif vac_or_air == 'air':
            if index=='MgI':
                cont = self.amgcont; band = self.amg
            elif index=='NaI':
                cont = self.anacont; band = self.ana
            elif index=='FeH':
                cont = self.afehcont; band = self.afeh
            elif index=='CaT':
                cont = self.acacont; band = self.aca
            elif index=='PaT':
                cont = self.acacont; band = self.apa
            elif index=='NaD':
                cont = self.anadcont; band = self.anad
            else:
                raise ValueError('Please choose one of the following indices:'+
                                 'MgI, NaI, FeH, CaT, PaT, NaD')

        #Interpolation of SED
        itp = interpolate.interp1d(SED[:,0], SED[:,1], kind='linear')

        eps1, eps2, eps3, eps4, eps5 = 0., 0., 0., 0., 0.

        #Method for calculating pseudo-continuum using equations from
        #Cenarro et al (2001) Appendix A.
        if var_SED==None:
            for j in xrange(0, len(cont)-1, 2):
                #Find first and last data indices within bandpass
                a = n.where(SED[:,0] > cont[j])[0][0]
                b = n.where(SED[:,0] < cont[j+1])[0][-1]
                #Determine which pixel is closest to start and end of bandpass
                if (cont[j] - SED[a-1,0]) < (SED[a,0] - cont[j]):
                    a -= 1
                if (cont[j+1] - SED[b,0]) > (SED[b+1,0] - cont[j+1]):
                    b += 1

                vals = n.copy(SED[a:b+1,:])

                e1 = n.ones(len(vals[:,0]))
                eps1 += n.sum(e1)
                e2 = vals[:,0]
                eps2 += n.sum(e2)
                e3 = vals[:,0]**2
                eps3 += n.sum(e3)
                e4 = vals[:,1]
                eps4 += n.sum(e4)
                e5 = vals[:,0]*vals[:,1]
                eps5 += n.sum(e5)

            delta = eps1*eps3 - eps2**2
            alpha1 = (1/float(delta))*(eps3*eps4 - eps2*eps5)
            alpha2 = (1/float(delta))*(eps1*eps5 - eps2*eps4)

        #If real spectrum with variance spectrum:
        elif var_SED != None:
            for j in xrange(0, len(cont)-1, 2):
                #Find first and last data indices within bandpass
                a = n.where(SED[:,0] > cont[j])[0][0]
                b = n.where(SED[:,0] < cont[j+1])[0][-1]
                #Determine which pixel is closest to start and end of bandpass
                if (cont[j] - SED[a-1,0]) < (SED[a,0] - cont[j]):
                    a -= 1
                if (cont[j+1] - SED[b,0]) > (SED[b+1,0] - cont[j+1]):
                    b += 1

                vals = n.copy(SED[a:b+1,:])
                err_vals = n.copy(var_SED[a:b+1,:])

                e1 = n.ones(len(vals[:,0]))/err_vals[:,1]
                eps1 += n.sum(e1)
                e2 = vals[:,0]/err_vals[:,1]
                eps2 += n.sum(e2)
                e3 = vals[:,0]**2/err_vals[:,1]
                eps3 += n.sum(e3)
                e4 = vals[:,1]/err_vals[:,1]
                eps4 += n.sum(e4)
                e5 = vals[:,0]*vals[:,1]/err_vals[:,1]
                eps5 += n.sum(e5)

            delta = eps1*eps3 - eps2**2
            alpha1 = (1/float(delta))*(eps3*eps4 - eps2*eps5)
            alpha2 = (1/float(delta))*(eps1*eps5 - eps2*eps4)

        #Test to check that continua works
        xs = n.arange(cont[0], cont[-1], 1)
        ys = alpha2*xs + alpha1
        P.plot(SED[:,0], SED[:,1], 'k-', xs, ys, 'r--')
        P.show()

        #Next stage, calculate index
        ind = 0.; ind_var_tot = 0.
        ind_var_1 = 0.; ind_var_2 = 0.
        for j in xrange(0, len(band)-1, 2):
            #Find first and last data indices within bandpass
            a = n.where(SED[:,0] > band[j])[0][0]
            b = n.where(SED[:,0] < band[j+1])[0][-1]
            #Determine which pixel is closest to start and end of bandpass
            if (band[j] - SED[a-1,0]) < (SED[a,0] - band[j]):
                a -= 1
            if (band[j+1] - SED[b,0]) > (SED[b+1,0] - band[j+1]):
                b += 1
            #Multiplicative factors for start and end pixels
            Cstart_c = (SED[a,0] - band[j] + 0.5*disp)/disp
            Cend_c = (band[j+1] - SED[b,0] + 0.5*disp)/disp

            Cvals = SED[a:b+1,:]

            contys = alpha2*Cvals[:,0] + alpha1
            Ccontvals = n.column_stack((Cvals[:,0],contys))

            array = (1-Cvals[:,1]/Ccontvals[:,1])
            array[0] *= Cstart_c
            array[-1] *= Cend_c
            value = disp*n.sum(array)
            ind += value

            ###ERRORS:
            if var_SED != None:
                #Index error:
                Cerrvals = var_SED[a:b+1,:]

                ind_cont_var = n.zeros(len(Ccontvals),dtype=float)
                #Loop over continuum feature pixels
                for x in xrange(len(ind_cont_var)):
                    #Calculate continuum variance:
                    Cvar_cont = 0.
                    for i in xrange(0, len(cont)-1, 2):
                        #Find first and last data indices within bandpass
                        ia = n.where(SED[:,0] > cont[i])[0][0]
                        ib = n.where(SED[:,0] < cont[i+1])[0][-1]
                        #Determine which pixel is closest to start and end of bandpass
                        if (cont[i] - SED[ia-1,0]) < (SED[ia,0] - cont[i]):
                            ia -= 1
                        if (cont[i+1] - SED[ib,0]) > (SED[ib+1,0] - cont[i+1]):
                            ib += 1
                        #Extract relevent region of SED and variance SED
                        contvals = n.copy(SED[ia:ib+1,:])
                        var_contvals = n.copy(var_SED[ia:ib+1,:])
                        #Apply equation
                        contsum = ((1/delta)*((1/var_contvals[:,1])*eps3 - (contvals[:,0]/var_contvals[:,1])*eps2)\
                                  + (Ccontvals[x,0]/delta)*((contvals[:,0]/var_contvals[:,1])*eps1\
                                                            - (1/var_contvals[:,1])*eps2))**2.*var_contvals[:,1]
                        Cvar_cont += n.sum(contsum)
                    ind_cont_var[x] = Cvar_cont

                Carray_var_1 = (Ccontvals[:,1]**2*Cerrvals[:,1] + Cvals[:,1]**2*ind_cont_var)/\
                             Ccontvals[:,1]**4.
                Carray_var_1[0] *= Cstart_c
                Carray_var_1 *= Cend_c
                ind_var_1 += n.sum(Carray_var_1)

                #For part 2 of equation which includes covariance matrix:
                #Loop over pixels in spectral feature
                for y in xrange(len(Cvals[:,0])):
                    #Loop over first spectral feature
                    for z in xrange(0, len(band)-1, 2):
                        #Find first and last data indices within bandpass
                        iia = n.where(SED[:,0] > band[z])[0][0]
                        iib = n.where(SED[:,0] < band[z+1])[0][-1]
                        #Determine which pixel is closest to start and end of bandpass
                        if (band[z] - SED[iia-1,0]) < (SED[iia,0] - band[z]):
                            iia -= 1
                        if (band[z+1] - SED[iib,0]) > (SED[iib+1,0] - band[z+1]):
                            iib += 1
                        #Multiplicative factors for start and end pixels
                        Fstart_c = (SED[iia,0] - band[z] + 0.5*disp)/float(disp)
                        Fend_c = (band[z+1] - SED[iib,0] + 0.5*disp)/float(disp)

                        Fvals = SED[iia:iib+1,:]
                        Fcontys = alpha2*Fvals[:,0] + alpha1
                        Fcontvals = n.column_stack((Fvals[:,0],Fcontys))
                        #Loop over pixels in second spectral feature
                        for zz in xrange(len(Fvals[:,0])):
                            cov_part = ((1/delta**2)*(eps1*eps3*eps3-eps2*eps2*eps3))+\
                                       ((1/delta**2)*(eps2*eps2*eps2-eps1*eps2*eps3))*(Cvals[y,0] + Fvals[zz,0])+\
                                       ((1/delta**2)*(eps1*eps1*eps3-eps1*eps2*eps2))*(Cvals[y,0]*Fvals[zz,0])

                            part_2_num = (Cvals[y,1]*Fvals[zz,1]*cov_part)/(Ccontvals[y,1]**2*Fcontvals[zz,1]**2)
                            #edge pixel factors
                            if y == 0.:
                                part_2_num *= Cstart_c
                            if zz == 0.:
                                part_2_num *= Fstart_c
                            if y == (len(Cvals[:,0])-1):
                                part_2_num *= Cend_c
                            if zz == (len(Fvals[:,0])-1):
                                part_2_num *= Fend_c

                            ind_var_2 += part_2_num

                ind_var_tot = disp**2*(ind_var_1 + ind_var_2)

        print  index+' = %.3f pm %.3f Angstroms' % (ind, n.sqrt(ind_var_tot))

        return ind, ind_var_tot



##################################################
    def CaT(self, SED, disp=1.0, var_SED=None, vac_or_air='vac'):
        '''Function to calculate the CaT* line index strength from an SED.

        Inputs:

            SED: spectral energy distribution of the form (lambda column, flux column)
            disp: float - spectral dispersion [A/pixel]
            var_SED: variance spectrum corresponding to SED of the same form as SED.
            vac_or_air: vacuum or air wavelengths for index definitions.

        Outputs:

            CaT_star: value of CaT* index
            err_CaT_star: uncertainty on CaT* measurement.
        '''
        if vac_or_air == 'vac':
            cont = self.cacont; ca = self.ca; pa = self.pa
        elif vac_or_air == 'air':
            cont = self.acacont; ca = self.aca; pa = self.apa
        eps1, eps2, eps3, eps4, eps5 = 0., 0., 0., 0., 0.
        #Method for calculating pseudo-continuum using equations from
        #Cenarro et al (2001) Appendix A.

        #If model spectrum with no error:
        if var_SED == None:
            for j in xrange(0, len(cont)-1, 2):
                #Find first and last data indices within bandpass
                a = n.where(SED[:,0] > cont[j])[0][0]
                b = n.where(SED[:,0] < cont[j+1])[0][-1]
                #Determine which pixel is closest to start and end of bandpass
                if (cont[j] - SED[a-1,0]) < (SED[a,0] - cont[j]):
                    a -= 1
                if (cont[j+1] - SED[b,0]) > (SED[b+1,0] - cont[j+1]):
                    b += 1

                vals = n.copy(SED[a:b+1,:])

                e1 = n.ones(len(vals[:,0]))
                eps1 += n.sum(e1)
                e2 = vals[:,0]
                eps2 += n.sum(e2)
                e3 = vals[:,0]**2
                eps3 += n.sum(e3)
                e4 = vals[:,1]
                eps4 += n.sum(e4)
                e5 = vals[:,0]*vals[:,1]
                eps5 += n.sum(e5)

            delta = eps1*eps3 - eps2**2
            alpha1 = (1/float(delta))*(eps3*eps4 - eps2*eps5)
            alpha2 = (1/float(delta))*(eps1*eps5 - eps2*eps4)

        #If real spectrum with variance spectrum:
        elif var_SED != None:
            for j in xrange(0, len(cont)-1, 2):
                #Find first and last data indices within bandpass
                a = n.where(SED[:,0] > cont[j])[0][0]
                b = n.where(SED[:,0] < cont[j+1])[0][-1]
                #Determine which pixel is closest to start and end of bandpass
                if (cont[j] - SED[a-1,0]) < (SED[a,0] - cont[j]):
                    a -= 1
                if (cont[j+1] - SED[b,0]) > (SED[b+1,0] - cont[j+1]):
                    b += 1

                vals = n.copy(SED[a:b+1,:])
                err_vals = n.copy(var_SED[a:b+1,:])

                e1 = n.ones(len(vals[:,0]))/err_vals[:,1]
                eps1 += n.sum(e1)
                e2 = vals[:,0]/err_vals[:,1]
                eps2 += n.sum(e2)
                e3 = vals[:,0]**2/err_vals[:,1]
                eps3 += n.sum(e3)
                e4 = vals[:,1]/err_vals[:,1]
                eps4 += n.sum(e4)
                e5 = vals[:,0]*vals[:,1]/err_vals[:,1]
                eps5 += n.sum(e5)

            delta = eps1*eps3 - eps2**2
            alpha1 = (1/float(delta))*(eps3*eps4 - eps2*eps5)
            alpha2 = (1/float(delta))*(eps1*eps5 - eps2*eps4)

##        #Test to check that continua works
##        xs = n.arange(cont[0], cont[-1], 1)
##        ys = alpha2*xs + alpha1
##        P.figure(figsize=(11,9), dpi=72)
##        P.plot(SED[:,0], SED[:,1], 'k-', xs, ys, 'r--', lw=1.5)
##        P.show()

        #Next stage, calculate CaT and PaT indices!
        Cat = 0; Pat = 0
        CaT_var_1 = 0.; PaT_var_1 = 0.; CaT_var_2 = 0.; PaT_var_2 = 0.
        CaT_var_tot = 0.; PaT_var_tot = 0.
        for j in xrange(0, len(ca)-1, 2):
            #For CaT
            #Find first and last data indices within bandpass
            Ca = n.where(SED[:,0] > ca[j])[0][0]
            Cb = n.where(SED[:,0] < ca[j+1])[0][-1]
            #Determine which pixel is closest to start and end of bandpass
            if (ca[j] - SED[Ca-1,0]) < (SED[Ca,0] - ca[j]):
                Ca -= 1
            if (ca[j+1] - SED[Cb,0]) > (SED[Cb+1,0] - ca[j+1]):
                Cb += 1
            #Multiplicative factors for start and end pixels
            Cstart_c = (SED[Ca,0] - ca[j] + 0.5*disp)/float(disp)
            Cend_c = (ca[j+1] - SED[Cb,0] + 0.5*disp)/float(disp)

            Cvals = SED[Ca:Cb+1,:]

            Ccontys = alpha2*Cvals[:,0] + alpha1
            Ccontvals = n.column_stack((Cvals[:,0],Ccontys))

            Carray = (1-Cvals[:,1]/Ccontvals[:,1])
            Carray[0] *= Cstart_c
            Carray[-1] *= Cend_c
            Cvalue = disp*n.sum(Carray)
            Cat += Cvalue

            #For PaT
            #Find first and last data indices within bandpass
            Pa = n.where(SED[:,0] > pa[j])[0][0]
            Pb = n.where(SED[:,0] < pa[j+1])[0][-1]
            #Determine which pixel is closest to start and end of bandpass
            if (pa[j] - SED[Pa-1,0]) < (SED[Pa,0] - pa[j]):
                Pa -= 1
            if (pa[j+1] - SED[Pb,0]) > (SED[Pb+1,0] - pa[j+1]):
                Pb += 1
            #Multiplicative factors for start and end pixels
            Pstart_c = (SED[Pa,0] - pa[j] + 0.5*disp)/float(disp)
            Pend_c = (pa[j+1] - SED[Pb,0] + 0.5*disp)/float(disp)

            Pvals = SED[Pa:Pb+1,:]

            Pcontys = alpha2*Pvals[:,0] + alpha1
            Pcontvals = n.column_stack((Pvals[:,0],Pcontys))

            Parray = (1-Pvals[:,1]/Pcontvals[:,1])
            Parray[0] *= Pstart_c
            Parray[-1] *= Pend_c
            Pvalue = disp*n.sum(Parray)
            Pat += Pvalue

            ###ERRORS:
            if var_SED != None:
                #CaT error:
                Cerrvals = var_SED[Ca:Cb+1,:]

                CaT_cont_var = n.zeros(len(Ccontvals),dtype=float)
                #Loop over continuum feature pixels
                for x in xrange(len(CaT_cont_var)):
                    #Calculate continuum variance:
                    Cvar_cont = 0.
                    for i in xrange(0, len(cont)-1, 2):
                        #Find first and last data indices within bandpass
                        a = n.where(SED[:,0] > cont[i])[0][0]
                        b = n.where(SED[:,0] < cont[i+1])[0][-1]
                        #Determine which pixel is closest to start and end of bandpass
                        if (cont[i] - SED[a-1,0]) < (SED[a,0] - cont[i]):
                            a -= 1
                        if (cont[i+1] - SED[b,0]) > (SED[b+1,0] - cont[i+1]):
                            b += 1
                        #Extract relevent region of SED and variance SED
                        contvals = n.copy(SED[a:b+1,:])
                        var_contvals = n.copy(var_SED[a:b+1,:])
                        #Apply equation
                        contsum = ((1/delta)*((1/var_contvals[:,1])*eps3 - (contvals[:,0]/var_contvals[:,1])*eps2)\
                                  + (Ccontvals[x,0]/delta)*((contvals[:,0]/var_contvals[:,1])*eps1\
                                                            - (1/var_contvals[:,1])*eps2))**2.*var_contvals[:,1]
                        Cvar_cont += n.sum(contsum)
                    CaT_cont_var[x] = Cvar_cont

                Carray_var_1 = (Ccontvals[:,1]**2*Cerrvals[:,1] + Cvals[:,1]**2*CaT_cont_var)/\
                             Ccontvals[:,1]**4.
                Carray_var_1[0] *= Cstart_c
                Carray_var_1 *= Cend_c
                CaT_var_1 += n.sum(Carray_var_1)

                #For part 2 of equation which includes covariance matrix:
                #Loop over pixels in spectral feature
                for y in xrange(len(Cvals[:,0])):
                    #Loop over first spectral feature
                    for z in xrange(0, len(ca)-1, 2):
                        #Find first and last data indices within bandpass
                        Fa = n.where(SED[:,0] > ca[z])[0][0]
                        Fb = n.where(SED[:,0] < ca[z+1])[0][-1]
                        #Determine which pixel is closest to start and end of bandpass
                        if (ca[z] - SED[Fa-1,0]) < (SED[Fa,0] - ca[z]):
                            Fa -= 1
                        if (ca[z+1] - SED[Fb,0]) > (SED[Fb+1,0] - ca[z+1]):
                            Fb += 1
                        #Multiplicative factors for start and end pixels
                        Fstart_c = (SED[Fa,0] - ca[z] + 0.5*disp)/float(disp)
                        Fend_c = (ca[z+1] - SED[Fb,0] + 0.5*disp)/float(disp)

                        Fvals = SED[Fa:Fb+1,:]
                        Fcontys = alpha2*Fvals[:,0] + alpha1
                        Fcontvals = n.column_stack((Fvals[:,0],Fcontys))
                        #Loop over pixels in second spectral feature
                        for zz in xrange(len(Fvals[:,0])):
                            cov_part = ((1/delta**2)*(eps1*eps3*eps3-eps2*eps2*eps3))+\
                                       ((1/delta**2)*(eps2*eps2*eps2-eps1*eps2*eps3))*(Cvals[y,0] + Fvals[zz,0])+\
                                       ((1/delta**2)*(eps1*eps1*eps3-eps1*eps2*eps2))*(Cvals[y,0]*Fvals[zz,0])

                            part_2_num = (Cvals[y,1]*Fvals[zz,1]*cov_part)/(Ccontvals[y,1]**2*Fcontvals[zz,1]**2)
                            #edge pixel factors
                            if y == 0.:
                                part_2_num *= Cstart_c
                            if zz == 0.:
                                part_2_num *= Fstart_c
                            if y == (len(Cvals[:,0])-1):
                                part_2_num *= Cend_c
                            if zz == (len(Fvals[:,0])-1):
                                part_2_num *= Fend_c

                            CaT_var_2 += part_2_num

                CaT_var_tot = disp**2*(CaT_var_1 + CaT_var_2)

                #PaT error:
                Perrvals = var_SED[Pa:Pb+1,:]

                PaT_cont_var = n.zeros(len(Pcontvals),dtype=float)
                for x in xrange(len(PaT_cont_var)):
                    #Calculate continuum variance:
                    Pvar_cont = 0.
                    for i in xrange(len(cont)-1, 2):
                        #Find first and last data indices within bandpass
                        a = n.where(SED[:,0] > cont[i])[0][0]
                        b = n.where(SED[:,0] < cont[i+1])[0][-1]
                        #Determine which pixel is closest to start and end of bandpass
                        if (cont[i] - SED[a-1,0]) < (SED[a,0] - cont[i]):
                            a -= 1
                        if (cont[i+1] - SED[b,0]) > (SED[b+1,0] - cont[i+1]):
                            b += 1
                        #Extract relevent region of SED and variance SED
                        contvals = n.copy(SED[a:b+1,:])
                        var_contvals = n.copy(var_SED[a:b+1,:])
                        #Apply equation
                        contsum = ((1/delta)*((1/var_contvals[:,1])*eps3 - (contvals[:,0]/var_contvals[:,1])*eps2)\
                                  + (Pcontvals[x,0]/delta)*((contvals[:,0]/var_contvals[:,1])*eps1\
                                                            - (1/var_contvals[:,1])*eps2))**2.*var_contvals[:,1]
                        Pvar_cont += contsum
                    PaT_cont_var[x] = Pvar_cont

                Parray_var_1 = 0.93**2*(Pcontvals[:,1]**2*Perrvals[:,1] + Pvals[:,1]**2*PaT_cont_var)/\
                             Pcontvals[:,1]**4.
                Parray_var_1[0] *= Pstart_c
                Parray_var_1 *= Pend_c
                PaT_var_1 += n.sum(Parray_var_1)

                #For part 2 of equation which includes covariance matrix:
                #Loop over pixels in spectral feature
                for y in xrange(len(Pvals[:,0])):
                    #Loop over first spectral feature
                    for z in xrange(0, len(pa)-1, 2):
                        #Find first and last data indices within bandpass
                        Fa = n.where(SED[:,0] > pa[z])[0][0]
                        Fb = n.where(SED[:,0] < pa[z+1])[0][-1]
                        #Determine which pixel is closest to start and end of bandpass
                        if (pa[z] - SED[Fa-1,0]) < (SED[Fa,0] - pa[z]):
                            Fa -= 1
                        if (pa[z+1] - SED[Fb,0]) > (SED[Fb+1,0] - pa[z+1]):
                            Fb += 1
                        #Multiplicative factors for start and end pixels
                        Fstart_c = (SED[Fa,0] - pa[z] + 0.5*disp)/float(disp)
                        Fend_c = (pa[z+1] - SED[Fb,0] + 0.5*disp)/float(disp)

                        Fvals = SED[Fa:Fb+1,:]
                        Fcontys = alpha2*Fvals[:,0] + alpha1
                        Fcontvals = n.column_stack((Fvals[:,0],Fcontys))
                        #Loop over pixels in second spectral feature
                        for zz in xrange(len(Fvals[:,0])):
                            cov_part = ((1/delta**2)*(eps1*eps3*eps3-eps2*eps2*eps3))+\
                                       ((1/delta**2)*(eps2*eps2*eps2-eps1*eps2*eps3))*(Pvals[y,0] + Fvals[zz,0])+\
                                       ((1/delta**2)*(eps1*eps1*eps3-eps1*eps2*eps2))*(Pvals[y,0]*Fvals[zz,0])

                            part_2_num = 0.93**2*(Pvals[y,1]*Fvals[zz,1]*cov_part)/(Pcontvals[y,1]**2*Fcontvals[zz,1]**2)
                            #edge pixel factors
                            if y == 0.:
                                part_2_num *= Pstart_c
                            if zz == 0.:
                                part_2_num *= Fstart_c
                            if y == (len(Cvals[:,0])-1):
                                part_2_num *= Pend_c
                            if zz == (len(Fvals[:,0])-1):
                                part_2_num *= Fend_c

                            PaT_var_2 += part_2_num

                PaT_var_tot = disp**2*(PaT_var_1 + PaT_var_2)



        CaT_star = Cat - 0.93*Pat
        CaT_star_var = CaT_var_tot + 0.93**2*PaT_var_tot

        print 'CaT = %.3f pm %.3f Angstroms' % (Cat, n.sqrt(CaT_var_tot))
        print 'PaT = %.3f pm %.3f Angstroms' % (Pat, n.sqrt(PaT_var_tot))
        print 'CaT* = %.3f pm %.3f Angstroms' % (CaT_star, n.sqrt(CaT_star_var))

        return CaT_star, CaT_star_var

################################################


def loadV03spec(filepath):

    ''' Function to read in V03 SSP SEDs and return as a
        spectrum class (created by RCWH).

        Inputs:

            filepath: Path and filename string of FITS file for V03 SED

        Outputs:

            spectrum: A spectrum class for the given SSP SED.
                      Initialised with units (lambda=A, flux=erg/s/cm2/A) @ D=10pc

    '''

    dat, head = p.getdata(filepath, 0, header=True)
    lambs = n.linspace(head['CRVAL1'], head['CRVAL1']+(len(dat)-1)\
                       *head['CDELT1'],len(dat))

    #Metallicity, age and IMF
    fname = filepath.split('/')[-1]
    z = fname.split('Z')[1]
    if 'm' in z:
        met = -1*float(z[1:5])
    elif 'm' not in z:
        met = float(z[1:4])

    Age = [float(fname.split('T')[1].split('.fits')[0])] #in Gyrs

    imf = fname[1:7]

    # Convert: 1/(L_sun*A) => / (L_sun[erg/s] /(4.*pi*D[cm]**2) => erg/s/cm**2/A (@10pc)
    # Note that D = 10pc!
    factor= L_sun/(4.*n.pi*(10.0*t.pc*100.0)**2.0)
    dat *= factor

    spectra = t.spectrum(lamspec=dat, lam=lambs, age=Age,
                         Z=met, model="V03")
    return spectra





def Gauss(sigma, delta_v):
    '''Function that creates an array of a normalised
    Gaussian distribution for use in convolutions.

    Inputs:

        sigma: Gaussian dispersion or velocity dispersion [km/s]
        delta_v: resolution [km/s/pixel]

    Outputs:

        column array of (x, y)
    '''
    #Create Gaussian
    num_vs = n.round(sigma*10/float(delta_v),decimals=0)
    if n.mod(num_vs,2) == 0.:
        num_vs += 1.
    xs = n.linspace(-5*sigma, 5*sigma, num=num_vs)
    ys = (1/(n.sqrt(2.*n.pi)*sigma))*n.exp(-0.5*(xs/float(sigma))**2)
    #normalise
    #integral = m.erf(5/n.sqrt(2.))
    #ys = ys * delta_v / (2.*float(integral))
    ys /= n.sum(ys)

    return n.column_stack((xs, ys))


def convolvespec(SED, sig, z=0):
    '''Function that convolves a spectrum with a Gaussian
    of given sigma (dispersion). Delta_v is calculated from
    the log spaced wavelength values.

    Inputs:

        SED: Spectrum in format (wavelength [A], flux)
        sig: velocity dispersion [km/s]
        z: redshift

    Outputs:

        convolved spectrum in format (wavelength [A], flux)

    '''

    #Logspaced wavelength array of same length as original wavelength array
    lambs = 10**n.linspace(n.log10(SED[0,0]), n.log10(SED[-1,0]), len(SED[:,0]))
    #Interpolate SED
    interspec = s.interpolate.interp1d(n.linspace(lambs[0], lambs[-1], len(lambs)),
                                       SED[:,1], kind='linear')

    #find SED values for log(lambda) bins
    logspec = interspec(lambs)

    #Find velocity spacing
    delv = (10**(n.log10(lambs[1])-n.log10(lambs[0]))-1)*2.998E5
    print 'delv = %.3f (km/s)/pixel' % delv


    #Create normalised Gaussian using sigma and delta_v
    if sig != 0.:
        gs = Gauss(sig, delv)
    else:
        gs = n.array([[1.,1.]])
    #Convolve flux array
    conspec = n.convolve(logspec, gs[:,1], mode='same')


#    #Redshift correction
#    logambs = loglambs - n.log10(1 + z)


    #Reinterpolate wavelengths back onto linear spaced wavelength array
    linlambs = n.linspace(lambs[0], lambs[-1], len(SED[:,0]))
    interspec2 = s.interpolate.interp1d(lambs, conspec, kind='linear')
    finalspec = interspec2(linlambs)

    return n.column_stack((linlambs, finalspec))


def cat_v_sig_interp():
    '''Function that returns interpolated function for CaT(sigma=0)/CaT(sigma)
    as a function of sigma. This fits data of a Z=0, 10Gyr, 1.30 Unimodal IMF
    star from Vazdekis model with an 8th order polynomial and return function.

    This will be used to calculate intrinsic CaT strength for observational data
    of Coma galaxies and M31/M32, to then populate CaT v Colour plots.
    '''

    dat = n.loadtxt('./sigma_v_cat_imfCun1.30z+0.00t10.00.txt', delimiter=',')
    coefs8 = n.polyfit(dat[:,0], dat[:,1], 8)
    dat8 = coefs8[0]*dat[:,0]**8 + coefs8[1]*dat[:,0]**7 + coefs8[2]*dat[:,0]**6 +\
    coefs8[3]*dat[:,0]**5 + coefs8[4]*dat[:,0]**4 + coefs8[5]*dat[:,0]**3 +\
    coefs8[6]*dat[:,0]**2 + coefs8[7]*dat[:,0] + coefs8[8]

    interp = interpolate.interp1d(dat[:,0], dat8, kind=3)

    return interp



def z_correct(SED, z):
    '''Function that takes an SED and corrects for redshift.

    Inputs:
        SED: Spectral energy distribution of form (Wavelength [A], Flux)
        z: Redshift - use negative value for blueshift

    Outputs:
        corrected_SED: redshift corrected SED
    '''

    interspec = s.interpolate.interp1d(SED[:,0]/float((1+z)), SED[:,1], kind='linear')
    lambs = n.linspace(SED[0,0]/float((1+z)), SED[-1,0]/float((1+z)), len(SED[:,0]))

    newflux = interspec(lambs)

    return n.column_stack((lambs, newflux))



def apply_redshift(SED, z, mode='correct'):
    '''Function that takes an SED and applies a redshift.
       Can choose to either correct for a shift or apply one.

    Inputs:
        SED: Spectral energy distribution of form (Wavelength [A], Flux)
        z: Redshift - use negative value for blueshift
        mode: Choose 'correct' or 'apply'. Either corrects for or applies a redshift.

    Outputs:
        corrected_SED: redshift corrected SED
    '''
    if mode == 'correct':
	factor = 1/float(1+z)
    elif mode == 'apply':
	factor = float(1+z)
    else:
	raise ValueError("mode must either be 'correct' or 'apply'.")


    interspec = s.interpolate.interp1d(SED[:,0]*factor, SED[:,1], kind='linear')
    lambs = n.linspace(SED[0,0]*factor, SED[-1,0]*factor, len(SED[:,0]))

    newflux = interspec(lambs)

    return n.column_stack((lambs, newflux))


def air_to_vac(air_wave):
    '''Function to convert air wavelengths to vacuum wavelengths

    Inputs:

        air_wave: air wavelength in Angstroms.

    Outputs:

        vac_wave: vacuum wavelength in Angstroms

    '''

    sigma2 = (10000./air_wave)**2.
    fact = 1.0 + 0.05792105/(238.0185 - sigma2) + 0.00167917/(57.362 - sigma2)
    vac_wave = air_wave*fact
    return vac_wave


def vac_to_air(vac_wave):
    '''Function to convert vacuum wavelengths to air wavelengths

    Inputs:

        vac_wave: air wavelength in Angstroms.

    Outputs:

        air_wave: vacuum wavelength in Angstroms

    '''

    sigma2 = (10000./vac_wave)**2.
    fact = 1.0 + 0.05792105/(238.0185 - sigma2) + 0.00167917/(57.362 - sigma2)
    air_wave = vac_wave/fact
    return air_wave


if __name__=="__main__":

    from scipy.io.idl import readsav

    data = readsav('/Users/SimonZ/Data/swift/M31-M32-Coma_stelpops/ppxf_M32.idl')
    spec = n.column_stack((n.round(data["wave"]*10000., 3), data["altbinflux"][:,0]))
    var_spec = n.column_stack((n.round(data["wave"]*10000., 3), data["binvar"][:,0]))
    v03 = cattoolset()
    v03.CaT(spec, var_SED=var_spec)
