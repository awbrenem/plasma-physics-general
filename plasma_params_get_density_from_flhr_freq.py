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
import plasmapy

mH = 1.6729e-27  #kg
mO = 2.6566e-26  
me = 9.1094e-31 

"""
Density based on lower hybrid frequency composed of fractional percentages of H+ and O+
"""
def dens_IonMassFractions(flh, fce, nH_ne, nO_ne):


    M = [1./((me/mH)*nH_ne[i] + (me/mO)*nO_ne[i]) for i in range(len(fce))]
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
