"""
Get plasma density from identification of lower hybrid frequency (Currently only for H+, O+ and mixed H+,O+ plasmas):

    def flhr_IonMassFractions(ne, fce, nH_ne, nO_ne)
        Lower hybrid frequency (not in high density limit) for fractional percentages of H+ and O+
    def dens_singleion(flh, Bo, species)
        Density based on single ion species


All frequencies in Hz, densities in cm-3

NOTE: all input quantities must be on same cadence!!!!
***Tested against the code plasma_params_get_flhr_freq.py
and the inversion works perfectly (ne->flh vs flh->ne)


If densities are negative this means that there is no density that can support 
the input lower hybrid frequency. NaN values are returned. 
"""


import numpy as np
#import plasmapy

#The following are the species in the IRI 2016 model
mH = 1.6729e-27  #kg
mO = 16*mH
mN = 14*mH
mHe = 4*mH
mO2 = 32*mH
mNO2 = (14+32)*mH

me = 9.1094e-31 

"""
Density based on lower hybrid frequency composed of fractional percentages of ion species 
NOTE: ions are those from the IRI 2016 model.

flh = lower hybrid freq (Hz)
fce = electron cyclotron freq (Hz)

nH_ne = fraction of ions that are H+
nO_ne = fraction of ions that are O+
nN_ne = fraction of ions that are N+

The following are only important at low altitudes (e.g. < 200 km)
nHe_ne = fraction of ions that are He+
nO2_ne = fraction of ions that are O2+
nNO2_ne = fraction of ions that are NO2+

(NOTE: that since O+ and N+ have nearly the same weight (16*mH vs 14*mH), modifying their factional percentages relative to each other
doesn't change the density much)

"""
def dens_IonMassFractions(flh, fce, nH_ne=1, nO_ne=0, nN_ne=0, nHe_ne=0, nO2_ne=0, nNO2_ne=0):

    flh = np.asarray(flh)
    fce = np.asarray(fce)
    nH_ne = np.asarray(nH_ne)
    nO_ne = np.asarray(nO_ne)
    nN_ne = np.asarray(nN_ne)
    nHe_ne = np.asarray(nHe_ne)
    nO2_ne = np.asarray(nO2_ne)
    nNO2_ne = np.asarray(nNO2_ne)

    #Turn single values into arrays with single value
    if not np.shape(flh):
        flh = [flh]
    if not np.shape(fce):
        fce = [fce]
    if not np.shape(nH_ne):
        nH_ne = [nH_ne]
    if not np.shape(nO_ne):
        nO_ne = [nO_ne]
    if not np.shape(nN_ne):
        nN_ne = [nN_ne]
    if not np.shape(nHe_ne):
        nHe_ne = [nHe_ne]
    if not np.shape(nO2_ne):
        nO2_ne = [nO2_ne]
    if not np.shape(nNO2_ne):
        nNO2_ne = [nNO2_ne]


    #Make sure all input arrays are of same length
    maxlen = len(max([[nH_ne],[nO_ne],[nN_ne],[nHe_ne],[nO2_ne],[nNO2_ne]], key=len)[0])
    if len(nH_ne) != maxlen:
        nH_ne = np.zeros(maxlen)
    if len(nO_ne) != maxlen:
        nO_ne = np.zeros(maxlen)
    if len(nN_ne) != maxlen:
        nN_ne = np.zeros(maxlen)
    if len(nHe_ne) != maxlen:
        nHe_ne = np.zeros(maxlen)
    if len(nO2_ne) != maxlen:
        nO2_ne = np.zeros(maxlen)
    if len(nNO2_ne) != maxlen:
        nNO2_ne = np.zeros(maxlen)



    M = [1./((me/mH)*nH_ne[i] + (me/mO)*nO_ne[i] + (me/mN)*nN_ne[i] + (me/mHe)*nHe_ne[i] + (me/mO2)*nO2_ne[i] + (me/mNO2)*nNO2_ne[i]) for i in range(len(fce))]
    C = 8980.**2

    n1 = [-M[i]*fce[i]**2 * flh[i]**2 for i in range(len(fce))]
    d1 = [C * (M[i]*flh[i]**2 - fce[i]**2) for i in range(len(fce))]

    ne = [n1[i]/d1[i] for i in range(len(fce))]

    ne = np.asarray(ne)
    bad = np.where(ne < 0)
    ne[bad] = np.nan


    """ Alternative calculation - gives same result 
    num = [(i/8980.)**2 for i in fce]
    den = [((fce[i]**2 * me) * (nH_ne[i]/mH + nO_ne[i]/mO)) / flh[i]**2 for i in range(len(fce))]
    ne2 = [num[i]/(den[i]-1) for i in range(len(fce))]
    """


    return ne


"""
Density based on single ion species
"""
def dens_singleion(flh, Bo, species):

    if not np.shape(Bo):
        Bo = [Bo]
    if not np.shape(flh):
        flh = [flh]


    fce = [28*i for i in Bo]

    Z = 1    #charge state (number of unmatched e-; qs/e)

    #atomic mass unit (mi/mp)
    if species == 'H+':
        muu = 1 
    elif species == 'O+':
        muu = 16
    elif species == 'He+':
        muu = 4
    elif species == 'N+':
        muu = 14

    Hplus2e_mass = mH/me       #H+ to electron mass ratio

    fci = [i*Z/(Hplus2e_mass*muu) for i in fce]


    gama = np.sqrt(muu*mH/me)

    num = [fce[i]*fci[i]*(gama/8980.)**2 for i in range(len(fce))]
    den = [fce[i]*fci[i]/(flh[i]**2) for i in range(len(fce))]
    ne = [num[i]/(den[i] - 1) for i in range(len(fce))]

    ne = np.asarray(ne)
    bad = np.where(ne < 0)
    ne[bad] = np.nan

    return ne




    """
    #Invert the equation for flh for fractional mass to solve for ne. 
    #(from plasma_params_flhr_freq.py)
    from sympy.solvers import solve
    from sympy import Symbol
    n = Symbol('n')
    flh = Symbol('flh')
    fce = Symbol('fce')
    C = Symbol('C')  
    M = Symbol('M')
    x = solve(C*n*(flh**2 - fce**2/M) + flh**2*fce**2, n)

    ne = -M*fce**2*flh**2/(C*(M*flh**2 - fce**2))
    """
