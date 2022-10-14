"""
Get lower hybrid frequency based on various calculations (Currently only for H+, O+ and mixed H+,O+ plasmas):

    def flhr_IonMassFractions(ne, fce, nH_ne, nO_ne)
        Lower hybrid frequency (not in high density limit) for fractional percentages of H+ and O+
    def flhr_singleion(ni, Bo, species)
        Lower Hybrid freq (not high density limit) for single ion species
    def flhr_HighDensityLimitTest
        Test to see if high density limit is applicable (fpi^2 >> fce*fci) 
    def flhr_HighDensityLimit(fce, fci)
        lower hybrid frequency in high density limit for single ion species

All frequencies in Hz, densities in cm-3

NOTE: all input quantities must be on same cadence!!!!
***Tested against the code plasma_params_get_density_from flhr_freq.py
and the inversion works perfectly (ne->flh vs flh->ne)

"""


import numpy as np
import plasmapy
#from astropy.constants import e
from astropy import units as u  


Oplus = plasmapy.particles.Particle("O 1+")
Hplus = plasmapy.particles.Particle("H 1+")
elec = plasmapy.particles.Particle("electron")

me = elec.mass
mH = Hplus.mass
mO = Oplus.mass



"""
Test to see if high density limit is applicable (fpi^2 >> fce*fci) 
"""
def flhr_HighDensityLimitTest(ne, Bo):

    #...use H+ to calculate ratio (same value for every ion species)
    fpi = [plasmapy.formulary.plasma_frequency(i, particle='H+', to_hz=True) for i in ne]
    fci = [plasmapy.formulary.gyrofrequency(Bo[i], particle='H+', to_hz=True) for i in range(len(ne))]
    fpi_fcefci = [fpi[i]**2/(fce[i]*fci[i]) for i in range(len(fce))]

    return {"fpi_fcefci":fpi_fcefci, "info":"High density limit only if fpi^2 >> fce*fci. Same value for every ion species"}


"""
lower hybrid frequency in high density limit for single ion species
"""
def flhr_HighDensityLimit(fce, fci):
    return [np.sqrt(fce[i]*fci[i]) for i in range(len(fce))]



"""
Lower hybrid frequency (not in high density limit) for fractional percentages of H+ and O+
"""
def flhr_IonMassFractions(ne, fce, nH_ne, nO_ne):

    fpe = [plasmapy.formulary.plasma_frequency(i, particle='electron', to_hz=True) for i in ne]

    nH = [nH_ne[i]*ne[i] for i in range(len(ne))]
    nO = [nO_ne[i]*ne[i] for i in range(len(ne))]
    
    Meff = [1./((me/ne[i])* ((nH[i]/mH) + (nO[i]/mO))) for i in range(len(ne))]

    goo = [np.sqrt((1./Meff[i])*((fce[i]**2.*fpe[i]**2.)/(fpe[i]**2. + fce[i]**2.))) for i in range(len(fce))] 
    return goo
    #return [np.sqrt((1./Meff[i])*((fce[i]**2.*fpe[i]**2.)/(fpe[i]**2. + fce[i]**2.))) for i in range(len(fce))] 



"""
Lower Hybrid freq (not high density limit) for single ion species
"""

def flhr_singleion(ni, Bo, species):

    #Kludgy, but if I don't do take the "value" and then reassign the units 
    #the below part doesn't work properly. 
    goo = [plasmapy.formulary.lower_hybrid_frequency(Bo[i].value * u.nT, ni[i].value * u.cm**-3, ion=species,to_hz = True) for i in range(len(ni))]
    return goo



 
if __name__ == '__main__': 
    print("Running as script")

    ne = [17632. * u.cm**-3, 18892. * u.cm**-3]
    Bo = [51908. * u.nT, 47280 * u.nT]

    flhrH = flhr_singleion(ne, Bo, 'H+')
    flhrO = flhr_singleion(ne, Bo, 'O+')



    fce = [plasmapy.formulary.gyrofrequency(Bo, 'e-', to_hz=True)]
    fcH = [plasmapy.formulary.gyrofrequency(Bo, 'H+', to_hz=True)]
    fcO = [plasmapy.formulary.gyrofrequency(Bo, 'O+', to_hz=True)]

    flhr = flhr_IonMassFractions(ne, fce, [0.0], [1.0])

    print(flhr, flhrH, flhrO)

    #Test high density limit 
#    flhrHD_tst = flhr_HighDensityLimitTest(ne, Bo)
#    flhrHD = flhr_HighDensityLimit(fce, fcH)
#    print(flhrHD_tst, flhrHD)


