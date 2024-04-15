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


#Oplus = plasmapy.particles.Particle("O 1+")
#Hplus = plasmapy.particles.Particle("H 1+")
#elec = plasmapy.particles.Particle("electron")

#me = elec.mass
#mH = Hplus.mass
#mO = Oplus.mass

me = 9.1093837e-31  #kg 
mH = 1.67291244e-27
mO = 2.65660536e-26

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
    fce = [np.asarray(i) for i in fce]
    fci = [np.asarray(i) for i in fci]

    
    return [np.sqrt(fce[i]*fci[i]) for i in range(len(fce))] 



"""
Lower hybrid frequency (not in high density limit) for fractional percentages of H+ and O+
"""
def flhr_IonMassFractions(ne, fce, nH_ne, nO_ne):

    ne = [np.asarray(i.value) for i in ne]
    fce = [np.asarray(i.value) for i in fce]
    fpe = 8980 * np.sqrt(ne)
    #nH_ne = [np.asarray(i.value) for i in nH_ne]
    #nO_ne = [np.asarray(i.value) for i in nO_ne]

    #fce = np.asarray(fce)
    #fpe = 8980 * np.sqrt(ne)
    #nH_ne = np.asarray(nH_ne)
    #nO_ne = np.asarray(nO_ne)
    nH = [nH_ne[i]*ne[i] for i in range(len(ne))]
    nO = [nO_ne[i]*ne[i] for i in range(len(ne))]

    #nH = [nH_ne[i]*ne[i] for i in range(len(ne))]
    #nO = [nO_ne[i]*ne[i] for i in range(len(ne))]
    #Meff = [1./((me/ne[i])* ((nH[i]/mH) + (nO[i]/mO))) for i in range(len(ne))]
    #goo = [np.sqrt((1/Meff[i])*((fce[i]**2*fpe[i]**2.)/(fpe[i]**2 + fce[i]**2))) for i in range(len(fce))] 
    
    Meff = [1./((me/ne[i]) * ((nH[i]/mH) + (nO[i]/mO))) for i in range(len(ne))]
    goo = [np.sqrt((1/Meff[i])*((fce[i]**2 * fpe[i]**2)/(fpe[i]**2 + fce[i]**2))) for i in range(len(ne))]

    return goo




"""
Lower Hybrid freq (not high density limit) for single ion species
"""

def flhr_singleion(ni, Bo, species):

    #Kludgy, but if I don't do take the "value" and then reassign the units 
    #the below part doesn't work properly. 
    goo = [plasmapy.formulary.lower_hybrid_frequency(Bo[i].value * u.nT, ni[i].value * u.cm**-3, ion=species,to_hz = True) for i in range(len(ni))]
    #goo = [plasmapy.formulary.lower_hybrid_frequency(Bo[i], ni[i], ion=species,to_hz = True) for i in range(len(ni))]
    return goo



 
if __name__ == '__main__': 
    print("Running as script")

    ne = [220058 * u.cm**-3, 246920 * u.cm**-3]
    Bo = [49375 * u.nT, 47669 * u.nT]
    #ne = [220058, 246920]
    #Bo = [49375, 47669]

    flhrH = flhr_singleion(ne, Bo, 'H+')
    flhrO = flhr_singleion(ne, Bo, 'O+')



    fce = [plasmapy.formulary.gyrofrequency(i, 'e-', to_hz=True) for i in Bo]
    fcH = [plasmapy.formulary.gyrofrequency(i, 'H+', to_hz=True) for i in Bo]
    fcO = [plasmapy.formulary.gyrofrequency(i, 'O+', to_hz=True) for i in Bo]

    flhr = flhr_IonMassFractions(ne, fce, [1.0,1.0], [0.0,0.0])

    print(flhr, flhrH, flhrO)

    #Test high density limit 
    flhrHD_tst = flhr_HighDensityLimitTest(ne, Bo)
    flhrHD = flhr_HighDensityLimit(fce, fcH)
    print(flhrHD_tst, flhrHD)


    print('h')
