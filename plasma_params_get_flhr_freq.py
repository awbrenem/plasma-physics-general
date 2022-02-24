"""
Get lower hybrid frequency based on various calculations:

    flhr_HighDensityLimitTest -> Test to see if high density limit is applicable (fpi^2 >> fce*fci) 

    flhr_HighDensityLimit(fce, fci) -> lower hybrid frequency in high density limit for single ion species

    flhr_IonMassFractions(ne, fce, nH_ne, nO_ne) -> Lower hybrid frequency (not in high density limit) for fractional percentages of H+ and O+

    flhr_H(ne, fce, fcH) -> Lower hybrid frequency (not in high density limit) for full H+
    flhr_O(ne, fce, fcO) -> Lower hybrid frequency (not in high density limit) for full O+

All frequencies in Hz, densities in cm-3
"""


import numpy as np

me     = 9.1093897e-31
mp     = 1.6726231e-27


#Test to see if high density limit is applicable (fpi^2 >> fce*fci) 
def flhr_HighDensityLimitTest(ne, fce, fcH):

    fpe = 8980*np.sqrt(ne) #Hz
    fpH = fpe/(np.sqrt(1.)*43.)

    #...use H+ to calculate ratio (same value for every ion species)
    fpH_fcefcH = [fpH**2/(fce[i]*fcH[i]) for i in range(len(fce))]

    return {"fpi_fcefci":fpH_fcefcH, "info":"High density limit only if fpi^2 >> fce*fci. Same value for every ion species"}


#lower hybrid frequency in high density limit for single ion species
def flhr_HighDensityLimit(fce, fci):
    return [np.sqrt(fce[i]*fci[i]) for i in range(len(fce))]



#Lower hybrid frequency (not in high density limit) for fractional percentages of H+ and O+
def flhr_IonMassFractions(ne, fce, nH_ne, nO_ne):

    fpe = 8980.*np.sqrt(ne)

    nH = nH_ne*ne #cm-3 
    nO = nO_ne*ne #cm-3 
    mH = mp 
    mO = 15.*mp
    
    Meff = 1./((me/ne)* ((nH/mH) + (nO/mO)))
    ftmp = [np.sqrt((1./Meff)*((fce[i]**2.*fpe**2.)/(fpe**2. + fce[i]**2.))) for i in range(len(fce))] 
    return ftmp 



#Lower hybrid frequency (not in high density limit) for full H+
def flhr_H(ne, fce, fcH):

    fpe = 8980*np.sqrt(ne) #Hz
    fpH = fpe/(np.sqrt(1.)*43.)
    return [np.sqrt(fpH**2*fce[i]*fcH[i])/np.sqrt(fce[i]*fcH[i]+fpH**2) for i in range(len(fce))]


#Lower hybrid frequency (not in high density limit) for full O+
def flhr_O(ne, fce, fcO):
    fpe = 8980*np.sqrt(ne) #Hz
    fpO = fpe/(np.sqrt(15.)*43.)
    return [np.sqrt(fpO**2*fce[i]*fcO[i])/np.sqrt(fce[i]*fcO[i]+fpO**2) for i in range(len(fce))]




if __name__ == '__main__': 
    print("Running as script")

