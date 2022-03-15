"""
Crib sheet with examples to determine whether a wave/structure is electrostatic or electromagnetic 
1) Calculate E/B ratio to determine wave group velocity
2) Compare to shear Alfven wave group velocity

**Note that the E/B ratios are invariant along a field line


GROUP VELOCITY AND RESONANCE ENERGY FOR ALFVEN WAVES FROM EW AND BW
See Chaston 2021 (Frontiers); Stasiewicz+00; Goertz and Boswell+79 for Alfven wave impedance (dispersion) relation
Eperp/Bperp = VA * sqrt((1 + kperp**2 * lambda_e**2)(1 + kperp**2 * rho_i**2))
However, the sqrt() terms are only relevant for inertial Alfven waves, so the usual expression is just 
-->Eperp/Bperp = VA

E in V/m, B in Tesla
Tesla = kg/s^2/A
Volt = kg*m^2/s^3/A
--> E/B = m/s

"""

#Rough rule of thumb in ionosphere: 500 mV/m = 10 km/s
Bo = 50000./1e9 #Tesla
Vsc = 10*1000. #m/s
Emot = Vsc * Bo #V/m
Emot = Emot * 1000. # = 500 mV/m


# Example 1 - RBSP at perigee (what motional E-field magnitude do we expect?)
Bo = 30000./1e9 #Tesla
Vsc = 10*1000. #m/s
Emot = Vsc * Bo #V/m
Emot = Emot * 1000. # = 300 mV/m (about what we'd expect)


# Example 2 - Dovner+94 - Test to see if solitary structures are electrostatic
# by comparing their E/B wave velocity to Alfven velocity of medium
Bw = 0.2/1e9 #T (max value)
Ew = 4./1000. #V/m
V = Ew/Bw # = 2d7  m/s

n = 4500. #cm-3
Bo = 20000. #nT
BnT2BG = 1e-5
#H+ Alfven velocity (m/s)
BnT2BG = 1e-5               #nT to Gauss (multiply by this number)
VA = 1000.*2.18e11*Bo*BnT2BG/sqrt(n)/100./1000. #=6.5d6  m/s
#Note that this is an upper limit b/c we only consider pure H+. In reality, probably
#lots of O+ at Freja location. 

#--> Thus we conclude that E/B velocity of structure is >> VA, and they are likely electrostatic



