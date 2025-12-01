#Read IRI outfile files from https://kauai.ccmc.gsfc.nasa.gov/instantrun/iri/
#Returns a numpy dictionary with values:    km  Ne/cm-3 Ne/NmF2 Tn/K  Ti/K  Te/K  O+  N+  H+ He+ O2+ NO+ Clust TEC t/%

#NOTE: to convert IRI data into data at rocket mission times use, for example, gir_cal_iri_vals_to_mission_times.py


#Below is example IRI output: 
"""

yyyy/mmdd(or -ddd)/hh.h):2025/ -33/ 8.6UT  geog Lat/Long/Alt= 67.0/ 213.0/ 300.0

IRIcor2 is used for topside Ne profile
OTSB2012 plasmasphere model w/o  plasmapause
URSI maps are used for the F2 peak density (NmF2)
foF2 STORM model is turned on 
Shubin2015 model is used for F2 peak height (hmF2)
ABT-2009 option is used for the bottomside thickness parameter B0
Scotto-97 no L   option is used for the F1 occurrence probability
foE auroral storm model is turned off
IRI-1990   option is used for D-region
TBT-2012 with solar dependence is used for the electron temperature
TBKS2021 option is used for the ion temperature
RBV10+TBT15 option is used for ion composition
CGM coordinate computation is turned off

Peak Densities/cm-3: NmF2= 154939.1   NmF1=      0.0   NmE=   3726.8
Peak Heights/km:     hmF2=   308.71   hmF1=     0.00   hmE=   110.00

Solar Zenith Angle/degree                             127.5
Dip (Magnetic Inclination)/degree                      78.48
Modip (Modified Dip)/degree                            65.47
Solar Sunspot Number (12-months running mean) Rz12     99.7           
Ionospheric-Effective Solar Index IG12                119.1           
Solar radio flux F10.7 (daily)                        209.5           
Solar radio flux F10.7 (81-day average)               183.5           

TEC [1.E16 m-2] is obtained by numerical integration in 1km steps
 from    65.0 to   600.0 km. t is the percentage of TEC above the F peak.

-
   H   ELECTRON DENSITY   TEMPERATURES      ION PERCENTAGES[%]*10    1E16m-2
   km  Ne/cm-3 Ne/NmF2 Tn/K  Ti/K  Te/K  O+  N+  H+ He+ O2+ NO+ Clust TEC t/%
    0.0     -1 -1.000    -1    -1    -1  -1  -1  -1  -1  -1  -1  -1   5.0  81
   10.0     -1 -1.000    -1    -1    -1  -1  -1  -1  -1  -1  -1  -1   5.0  81
   20.0     -1 -1.000    -1    -1    -1  -1  -1  -1  -1  -1  -1  -1   5.0  81

"""


import pandas as pd
import matplotlib.pyplot as plt
import numpy as np 


def read_iri_output(filename):
    # Find the start of the data table
    with open(filename, 'r') as f:
        lines = f.readlines()
    # Find the line index where the table starts (look for header line)
    for i, line in enumerate(lines):
        if line.strip().startswith('H   ELECTRON DENSITY'):
            header_idx = i
            break
    # The next line is the column names, skip header line
    data_start = header_idx + 2
    # Read the data into a DataFrame
    df = pd.read_csv(filename, delim_whitespace=True, skiprows=data_start,
                        names=['H_km', 'Ne_cm-3', 'Ne_NmF2', 'Tn_K', 'Ti_K', 'Te_K',
                            'O+', 'N+', 'H+', 'He+', 'O2+', 'NO+', 'Clust', 'TEC_1E16m-2', 't_%'],
                        comment='-')


    # Convert DataFrame to a numpy dictionary
    np_dict = {col: df[col].to_numpy() for col in df.columns}

    # Turn the 10*percentage values into fractions
    np_dict['O+'] = np_dict['O+']/1000.
    np_dict['N+'] = np_dict['N+']/1000.
    np_dict['H+'] = np_dict['H+']/1000.
    np_dict['He+'] = np_dict['He+']/1000.
    np_dict['O2+'] = np_dict['O2+']/1000.
    np_dict['NO+'] = np_dict['NO+']/1000.

    return np_dict




#filename = '/Users/abrenema/Desktop/Research/Rocket_missions/GIRAFF/data/IRI/iri_output_36381.txt'
#df = read_iri_output(filename)

#print(np_dict)
#plt.plot(np.asarray(df['H_km']), np.asarray(df['Ne_cm-3']))


#print('h')