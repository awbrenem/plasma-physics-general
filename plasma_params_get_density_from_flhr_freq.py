"""
Get plasma density from identification of lower hybrid frequency (Currently only for H+, O+ and mixed H+,O+ plasmas):

    def flhr_IonMassFractions(ne, fce, nH_ne, nO_ne)
        Lower hybrid frequency (not in high density limit) for fractional percentages of H+ and O+
    def flhr_H(ne, fce, fcH)
        Lower hybrid frequency (not in high density limit) for full H+
    def flhr_O(ne, fce, fcO)
        Lower hybrid frequency (not in high density limit) for full O+

All frequencies in Hz, densities in cm-3

NOTE: all input quantities must be on same cadence!!!!
"""


import numpy as np

me     = 9.1093897e-31
mp     = 1.6726231e-27
mH = mp         #H+ mass
mO = 15.*mp     #O+ mass




"""
Density based on lower hybrid frequency composed of fractional percentages of H+ and O+
"""
def dens_IonMassFractions(flh, fce, nH_ne, nO_ne):

    rH = nH_ne 
    rO = nO_ne

    num = (fce/8980)**2
    den = ((fce**2 * me) * (rH/mH + rO/mO)) / flh**2  
    ne = num/(den - 1)


    """
    fpe = [8980.*np.sqrt(ne[i]) for i in range(len(ne))]

    nH = [nH_ne[i]*ne[i] for i in range(len(ne))] #cm-3 
    nO = [nO_ne[i]*ne[i] for i in range(len(ne))] #cm-3 
    
    Meff = [1./((me/ne[i])* ((nH[i]/mH) + (nO[i]/mO))) for i in range(len(ne))]
    return [np.sqrt((1./Meff[i])*((fce[i]**2.*fpe[i]**2.)/(fpe[i]**2. + fce[i]**2.))) for i in range(len(fce))] 
    """
    return ne

"""
Lower hybrid frequency (not in high density limit) for full H+
"""
def dens_H(flh, fce, fcH):

    num = (fce/8980)**2
    den = (me/mH) * (fce**2/flh**2)
    ne = num/(den-1)

    #fpe = [8980.*np.sqrt(ne[i]) for i in range(len(ne))] #Hz
    #fpH = [fpe[i]/(np.sqrt(1.)*43.) for i in range(len(ne))]
    #return [np.sqrt(fpH[i]**2*fce[i]*fcH[i])/np.sqrt(fce[i]*fcH[i]+fpH[i]**2) for i in range(len(fce))]

    return ne 


"""
Lower hybrid frequency (not in high density limit) for full O+
"""
def dens_O(flh, fce, fcO):

    num = (fce/8980)**2
    den = (me/mO) * (fce**2/flh**2)
    ne = num/(den-1)

    #fpe = [8980*np.sqrt(ne[i]) for i in range(len(ne))] #Hz
    #fpO = [fpe[i]/(np.sqrt(15.)*43.) for i in range(len(ne))]
    #return [np.sqrt(fpO[i]**2*fce[i]*fcO[i])/np.sqrt(fce[i]*fcO[i]+fpO[i]**2) for i in range(len(fce))]

    return ne 



if __name__ == '__main__': 
    print("Running as script")

