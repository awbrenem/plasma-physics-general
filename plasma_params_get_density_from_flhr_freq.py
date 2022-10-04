"""
Get plasma density from identification of lower hybrid frequency (Currently only for H+, O+ and mixed H+,O+ plasmas):

    def flhr_IonMassFractions(ne, fce, nH_ne, nO_ne)
        Lower hybrid frequency (not in high density limit) for fractional percentages of H+ and O+
    def flhr_H(ne, fce, fcH)
        Lower hybrid frequency (not in high density limit) for 100% H+
    def flhr_O(ne, fce, fcO)
        Lower hybrid frequency (not in high density limit) for 100% O+

All frequencies in Hz, densities in cm-3

NOTE: all input quantities must be on same cadence!!!!
"""


import numpy as np

me     = 9.1093897e-31
mp     = 1.6726231e-27
mH = mp         #H+ mass
mO = 15.*mp     #O+ mass
Hplus2e_mass = mp/me



"""
Density based on lower hybrid frequency composed of fractional percentages of H+ and O+
"""
def dens_IonMassFractions(flh, fce, nH_ne, nO_ne):

    rH = nH_ne 
    rO = nO_ne

    num = (fce/8980)**2
    den = ((fce**2 * me) * (rH/mH + rO/mO)) / flh**2  
    ne = num/(den - 1)

    return ne


"""
Lower hybrid frequency (not in high density limit) for 100% H+
"""
def dens_H(flh, fce, fcH):

    #Version 1 based on reduced form of nO and nH full version
    num = (fce/8980)**2
    den = (me/mH) * (fce**2/flh**2)
    ne = num/(den-1)

    ##Version 2 (**Checked - gives same answer as Version 1)
    #gama = np.sqrt(1)*43
    #num = fce*fcH*(gama/8980)**2
    #den = fce*fcH/(flh**2)
    #ne = num/(den - 1)

    return ne 


"""
Lower hybrid frequency (not in high density limit) for 100% O+
"""
def dens_O(flh, fce, fcO):

    num = (fce/8980)**2
    den = (me/mO) * (fce**2/flh**2)
    ne = num/(den-1)

    #Version 2 (**Checked - gives same answer as Version 1)
    gama = np.sqrt(15)*43
    num = fce*fcO*(gama/8980)**2
    den = fce*fcO/(flh**2)
    ne = num/(den - 1)


    return ne 



if __name__ == '__main__': 
    print("Running as script")
    flh = 1000.
    fce = 1e6 
    nH_ne = 0. 
    nO_ne = 1.
    ne1 = dens_IonMassFractions(flh, fce, nH_ne, nO_ne)

    fcH = fce/(Hplus2e_mass*1)
    ne2 = dens_H(flh, fce, fcH)


    fcO = fce/(Hplus2e_mass*15)
    ne3 = dens_O(flh, fce, fcO)


