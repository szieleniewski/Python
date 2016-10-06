'''Functions to convert between air and vacuum wavelengths

Author: Simon Zieleniewski

Last updated: 06-04-16

'''

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
