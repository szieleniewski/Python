import numpy as n
from scipy import interpolate
import spectools as t


#last updated 21-04-16


L_sun = 3.826E33 # the L_sun in erg/s defined by Vazdekis et al. (2012)


def loadCvD12spec(filepath):
    ''' Function to read in CvD12 SSP files and return spectra as a
        spectrum class (created by RCWH).

        Inputs:

            filepath: Path and filename string of file for CvD spectra

        Outputs:

            spectrum: A spectrum class for the given SSP SED.
                      Initialised with units (lambda=A, flux=erg/s/cm2/A) @ D=10pc

    '''

    dat = n.genfromtxt(filepath)

    #Get filename
    fname = filepath.split('/')[-1]

    #Wavelenghts in A
    lambs = dat[:,0].copy()

    #Set flux units to erg/s/cm**2/A at D = 10 pc. CvD flux in units of L_sun/um
    flux = n.transpose(dat[:,1:].copy())
    factor = L_sun/(10000.*(4.*n.pi*(10.0*t.pc*100.0)**2.0))
    #factor.shape = (1, len(lambs))
    flux *= factor

    #Age of spectra in Gyrs
    Age = [float(fname.split('_')[0].split('t')[1])]*flux.shape[0]

    #Interpolate to get linear dispersion
    newlambs = n.linspace(lambs[0], lambs[-1], len(lambs))
    finterp = interpolate.interp1d(lambs, flux, kind='linear', axis=-1)
    newflux = finterp(newlambs)

    #Get mass file
    masspath = filepath.split(fname)[0]
    masses = n.genfromtxt(masspath+'mass_ssp.dat')
    masspos = {13.5:6, 11.0:5, 9.0:4, 7.0:3, 5.0:2, 3.0:1}
    mass = masses[:,masspos[Age[0]]]

    #Depending on filename, spectra correspond to different IMFs, ages etc
    if 'solar' in fname:
        #IMFs = x=3.5, 3.0, 2.35, Chabrier, bottom-light
        IMFs = 'x = 3.5, x = 3.0, x = 2.35, Chabrier, bottom-light'
        return t.spectrum(lamspec=newflux, lam=newlambs, age=Age,
                          Z=0.2, IMF=IMFs, model='CvD12', mass=mass)

    if 'afe' in fname:
        met = 0.0
        IMFs = 'x = 3.5, x = 3.0, x = 2.35, Chabrier, bottom-light'
        afes = {'afe': float(fname.split('+')[1][0:3])}
        return t.spectrum(lamspec=newflux, lam=newlambs, age=Age,
                          Z=met, IMF=IMFs, model='CvD12', userdict=afes,
                          mass=mass)

    if 'varelem' in fname:
        IMFs = 'Chabrier'
        #Need to finish! - 08-12-13
        met = 0.0
        abunlist = {'abundances': ['[Na/Fe] = +0.3', '[Na/Fe] = -0.3','[Ca/Fe] = +0.15', '[Ca/Fe] = -0.15',
                '[Fe/H] = +0.3', '[Fe/H] = -0.3', '[C/Fe] = +0.15', '[C/Fe] = -0.15',
                '[a/Fe] = +0.2', '[as/Fe] = +0.2', '[N/Fe] = +0.3', '[N/Fe] = -0.3',
                '[Ti/Fe] = +0.3', '[Ti/Fe] = -0.3', '[Mg/Fe] = +0.3', '[Mg/Fe] = -0.3',
                '[Si/Fe] = +0.3', '[Si/Fe] = -0.3']}

##        afelist = {'[Na/Fe]':+0.3, '[Na/Fe]':-0.3, '[Ca/Fe]':+0.15, '[Ca/Fe]':-0.15,
##                '[Fe/H]':+0.3, '[Fe/H]':-0.3, '[C/Fe]':+0.15, '[C/Fe]':-0.15,
##                '[a/Fe]':+0.2, '[as/Fe]':+0.2, '[N/Fe]':+0.3, '[N/Fe]':-0.3,
##                '[Ti/Fe]':+0.3, '[Ti/Fe]':-0.3, '[Mg/Fe]':+0.3, '[Mg/Fe]':-0.3,
##                '[Si/Fe]':+0.3, '[Si/Fe]':-0.3}
        return t.spectrum(lamspec=newflux, lam=newlambs, age=Age,
                          Z=met, IMF=IMFs, model='CvD12', userdict=abunlist,
                          mass=mass)

    else:
        raise ValueError('Did not input correct CvD12 file [as of 03-04-14]')
