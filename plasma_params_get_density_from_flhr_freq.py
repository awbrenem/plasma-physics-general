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

"""


import numpy as np
import plasmapy
from astropy.constants import m_e, m_p
from astropy import units as u  


Oplus = plasmapy.particles.Particle("O 1+")
Hplus = plasmapy.particles.Particle("H 1+")
elec = plasmapy.particles.Particle("electron")


mH = Hplus.mass
mO = Oplus.mass
me = elec.mass



"""
Density based on lower hybrid frequency composed of fractional percentages of H+ and O+
"""
def dens_IonMassFractions(flh, fce, nH_ne, nO_ne):

    num = [(i/8980)**2 for i in fce]
    den = [((fce[i]**2 * me) * (nH_ne[i]/mH + nO_ne[i]/mO)) / flh[i]**2 for i in range(len(fce))]
    ne = [num[i]/(den[i]-1) for i in range(len(fce))]

    return ne


"""
Density based on single ion species
"""
def dens_singleion(flh, Bo, species):

    fce = [plasmapy.formulary.gyrofrequency(i, particle='electron', to_hz=True) for i in Bo]
    fci = [plasmapy.formulary.gyrofrequency(i, particle=species, to_hz=True) for i in Bo]

    ion = plasmapy.particles.Particle(species)
    gama = np.sqrt(ion.mass/me)


    num = [fce[i]*fci[i]*(gama/8980)**2 for i in range(len(fce))]
    den = [fce[i]*fci[i]/(flh[i]**2) for i in range(len(fce))]
    ne = [num[i]/(den[i] - 1) for i in range(len(fce))]

    return ne 



if __name__ == '__main__': 
    print("Running as script")
    flh = [4000.] * u.Hz
    Bo = [43000.] * u.nT
    nH_ne = [1.] * u.dimensionless_unscaled 
    nO_ne = [0.] * u.dimensionless_unscaled



    fce = [plasmapy.formulary.gyrofrequency(i, particle='electron', to_hz=True) for i in Bo]
    fcH = [plasmapy.formulary.gyrofrequency(i, particle='H+', to_hz=True) for i in Bo]
    fcO = [plasmapy.formulary.gyrofrequency(i, particle='O+', to_hz=True) for i in Bo]

    ne1 = dens_IonMassFractions(flh, fce, nH_ne, nO_ne)
    ne2 = dens_singleion(flh, Bo, 'H+')
    ne3 = dens_singleion(flh, Bo, 'O+')

    print(ne1, ne2, ne3)




    print('done')

