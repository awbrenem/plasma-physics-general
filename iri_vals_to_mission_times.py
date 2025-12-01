#Convert IRI model values (a function of altitude) to mission times

#iridat - data from the IRI model read in from iri_read_output.py
#mission_times - array of mission times (s)
#mission_altitudes - array of mission altitudes (km)


def iri_vals_to_mission_times(iridat, mission_times, mission_altitudes, test=0):
        
    import numpy as np 

    #For each time of ephem altitude data find the value of fpe from IRI at that altitude
    for i in range(len(mission_times)):
        alt = mission_altitudes[i]
        #Find nearest altitude in IRI data
        idx = (np.abs(iridat['H_km'] - alt)).argmin()
 
        if i == 0:
            H_km_interp = [iridat['H_km'][idx]]
            Ne_cm3_interp = [iridat['Ne_cm-3'][idx]]
            Ne_NmF2_interp = [iridat['Ne_NmF2'][idx]]
            Tn_K_interp = [iridat['Tn_K'][idx]]
            Ti_K_interp = [iridat['Ti_K'][idx]]
            Te_K_interp = [iridat['Te_K'][idx]]
            Oplus_interp = [iridat['O+'][idx]]
            Nplus_interp = [iridat['N+'][idx]]
            Hplus_interp = [iridat['H+'][idx]]
            Heplus_interp = [iridat['He+'][idx]]
            O2plus_interp = [iridat['O2+'][idx]]
            NOplus_interp = [iridat['NO+'][idx]]
        else:
            H_km_interp = np.append(H_km_interp, iridat['H_km'][idx])
            Ne_cm3_interp = np.append(Ne_cm3_interp, iridat['Ne_cm-3'][idx])
            Ne_NmF2_interp = np.append(Ne_NmF2_interp, iridat['Ne_NmF2'][idx])
            Tn_K_interp = np.append(Tn_K_interp, iridat['Tn_K'][idx])
            Ti_K_interp = np.append(Ti_K_interp, iridat['Ti_K'][idx])
            Te_K_interp = np.append(Te_K_interp, iridat['Te_K'][idx])
            Oplus_interp = np.append(Oplus_interp, iridat['O+'][idx])
            Nplus_interp = np.append(Nplus_interp, iridat['N+'][idx])
            Hplus_interp = np.append(Hplus_interp, iridat['H+'][idx])
            Heplus_interp = np.append(Heplus_interp, iridat['He+'][idx])
            O2plus_interp = np.append(O2plus_interp, iridat['O2+'][idx])
            NOplus_interp = np.append(NOplus_interp, iridat['NO+'][idx])



    if test:
        import matplotlib.pyplot as plt
        fig, axs = plt.subplots(5, 1, figsize=(10, 8))

        axs[0].plot(mission_times, H_km_interp, label='alt IRI', color='blue', linestyle='--')
        axs[0].plot(mission_times, mission_altitudes, label='alt mission', color='red', linestyle='-')
 
        axs[1].plot(mission_times, Oplus_interp, label='O+ IRI', color='blue', linestyle='--')
        axs[1].plot(mission_times, Nplus_interp, label='N+ IRI', color='magenta', linestyle='--')
        axs[1].plot(mission_times, Hplus_interp, label='H+ IRI', color='red', linestyle='--')
        axs[1].plot(mission_times, Heplus_interp, label='He+ IRI', color='green', linestyle='--')
        axs[1].plot(mission_times, O2plus_interp, label='O2+ IRI', color='orange', linestyle='--')
        axs[1].plot(mission_times, NOplus_interp, label='NO+ IRI', color='brown', linestyle='--')

        axs[2].plot(mission_times, Ne_cm3_interp, label='Ne IRI', color='black', linestyle='--')
        axs[3].plot(mission_times, Ne_NmF2_interp, label='Ne_NmF2 IRI', color='black', linestyle='--')

        axs[4].plot(mission_times, Tn_K_interp, label='Tn IRI', color='black', linestyle='--')
        axs[4].plot(mission_times, Ti_K_interp, label='Ti IRI', color='blue', linestyle='-')
        axs[4].plot(mission_times, Te_K_interp, label='Te IRI', color='red', linestyle='--')
        for ax in axs:
            ax.legend()
        axs[0].set_ylabel('Altitude (km)')
        axs[1].set_ylabel('Ion Fractions')
        axs[2].set_ylabel('Ne (cm-3)')
        axs[3].set_ylabel('Ne_NmF2 (cm-3)')
        axs[4].set_ylabel('Tn, Ti, Te (K)')
        plt.xlabel('Time (s)')
        plt.ylabel('Density (cm-3)')
        plt.title('IRI Model Values at Mission Times')
        plt.legend()
        plt.show()  

    return {
        'times': np.array(mission_times),
        'H_km': np.array(H_km_interp),
        'Ne_cm-3': np.array(Ne_cm3_interp),
        'Ne_NmF2': np.array(Ne_NmF2_interp),
        'Tn_K': np.array(Tn_K_interp),
        'Ti_K': np.array(Ti_K_interp),
        'Te_K': np.array(Te_K_interp),
        'O+': np.array(Oplus_interp),
        'N+': np.array(Nplus_interp),
        'H+': np.array(Hplus_interp),
        'He+': np.array(Heplus_interp),
        'O2+': np.array(O2plus_interp),
        'NO+': np.array(NOplus_interp)
    }



