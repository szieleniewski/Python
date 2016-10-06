import numpy as n
import astropy.io.fits as p
import spectools as t


#last updated 14-10-15

L_sun = 3.826E33 # the L_sun in erg/s defined by Vazdekis et al. (2012)

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



def loadV15spec(filepath):

    ''' Function to read in V15 SSP SEDs and return as a
        spectrum class (created by RCWH).

        Inputs:

            filepath: Path and filename string of FITS file for V15 SED

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

    Age = [float(fname.split('T')[1].split('_')[0])] #in Gyrs

    imf = fname[1:7]

    # Convert: 1/(L_sun*A) => / (L_sun[erg/s] /(4.*pi*D[cm]**2) => erg/s/cm**2/A (@10pc)
    # Note that D = 10pc!
    factor= L_sun/(4.*n.pi*(10.0*t.pc*100.0)**2.0)
    dat *= factor

    spectra = t.spectrum(lamspec=dat, lam=lambs, age=Age,
                         Z=met, model="V03")
    return spectra
