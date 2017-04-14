;Crib sheet for calling cycl_energies_parameter_space.pro



;--------------------------------------------------
;TEST OF NEW CYCLOTRON RES CALCULATION (compare to Lorentzen01, plate 7)
;--------------------------------------------------


ps = 0      ;save to postscript?

pa=5.        ;e- pitch angle
theta_k=0.   ;wave normal angle
density_range = [1,10]   ;cm-3
fce_range = [800,25000]  ;Hz
freq_range = [0,8000]    ;Hz
minval = 1.    ;minimum energy plotted (keV)
maxval = 1000.   ;maximum energy plotted (keV)

ndens = 20  &  nfce = 20  &  nfreq = 20  ;number of contours

;scheme 0:  fce vs density (at constant wave freq)
cycl_energies_parameter_space,pa,theta_k,scheme=0,ps=ps,minval=minval,maxval=maxval,$
                              freqv=1500.,fce_range=fce_range,density_range=density_range,$
                              ndens=ndens,nfce=nfce,nfreq=nfreq,harmonic=1


;scheme 1:  fce vs freq (at constant density)
cycl_energies_parameter_space,pa,theta_k,scheme=1,ps=ps,minval=minval,maxval=maxval,$
                              densv=5.,fce_range=fce_range,freq_range=freq_range,$
                              ndens=ndens,nfce=nfce,nfreq=nfreq;,type='landau'


;scheme 2:  density vs freq (at constant fce)
cycl_energies_parameter_space,pa,theta_k,scheme=2,ps=ps,minval=minval,maxval=maxval,$
                              fcev=7000.,density_range=density_range,freq_range=freq_range,$
                              ndens=ndens,nfce=nfce,nfreq=nfreq;,type='landau'



;--------------------------------------------------
;For Cindy's droopy whistlers
;--------------------------------------------------


ps = 0      ;save to postscript?

pa=20.        ;e- pitch angle
theta_k=10.   ;wave normal angle
density_range = [1,60]   ;cm-3
fce_range = [1000,7000]  ;Hz
freq_range = [80,300]    ;Hz
minval = 0.001    ;minimum energy plotted (keV)
maxval = 100.   ;maximum energy plotted (keV)

ndens = 20  &  nfce = 20  &  nfreq = 20  ;number of contours

;scheme 0:  fce vs density (at constant wave freq)
cycl_energies_parameter_space,pa,theta_k,scheme=0,ps=ps,minval=minval,maxval=maxval,$
                              freqv=100.,fce_range=fce_range,density_range=density_range,$
                              ndens=ndens,nfce=nfce,nfreq=nfreq


;scheme 1:  fce vs freq (at constant density)
cycl_energies_parameter_space,pa,theta_k,scheme=1,ps=ps,minval=minval,maxval=maxval,$
                              densv=10.,fce_range=fce_range,freq_range=freq_range,$
                              ndens=ndens,nfce=nfce,nfreq=nfreq


;scheme 2:  density vs freq (at constant fce)
cycl_energies_parameter_space,pa,theta_k,scheme=2,ps=ps,minval=minval,maxval=maxval,$
                              fcev=5000.,density_range=density_range,freq_range=freq_range,$
                              ndens=ndens,nfce=nfce,nfreq=nfreq


;--------------------------------------------------
;For Lois Keller's L=2 whistlers
;--------------------------------------------------



ps = 1     ;save to postscript?

pa=0.        ;e- pitch angle
theta_k=2.   ;wave normal angle
density_range = [2000,15000]   ;cm-3
fce_range = [22000,26400]  ;Hz
freq_range = [16000.,22000.]    ;Hz
minval = 0.0001
maxval = 10.   ;maximum energy plotted (keV)

ndens = 20  &  nfce = 20  &  nfreq = 20  ;number of contours


;--------------------------------------------------
;Test values to be sure of validity of energy calculations
;--------------------------------------------------

x = cold_dispersion(epol=1,freq=22361.,dens=10000.,Bo=800.)
function cold_dispersion,epol=epol,freq=freq,dens=dens,Bo=Bo,H_plus=pH,He_plus=pHe,O_plus=pO

res_angle_max = acos(max(freq_range)/min(fce_range))/!dtor  ;resonance cone angle
res_angle_min = acos(min(freq_range)/max(fce_range))/!dtor  ;resonance cone angle
;--------------------------------------------------

;scheme 0:  fce vs density (at constant wave freq)
cycl_energies_parameter_space,pa,theta_k,scheme=0,ps=ps,minval=minval,maxval=maxval,$
                              freqv=22361.,fce_range=fce_range,density_range=density_range,$
                              ndens=ndens,nfce=nfce,nfreq=nfreq,type='landau'


;scheme 1:  fce vs freq (at constant density)
cycl_energies_parameter_space,pa,theta_k,scheme=1,ps=ps,minval=minval,maxval=maxval,$
                              densv=2000.,fce_range=fce_range,freq_range=freq_range,$
                              ndens=ndens,nfce=nfce,nfreq=nfreq;,type='landau'


;scheme 2:  density vs freq (at constant fce)
cycl_energies_parameter_space,pa,theta_k,scheme=2,ps=ps,minval=minval,maxval=maxval,$
                              fcev=22400.,density_range=density_range,freq_range=freq_range,$
                              ndens=ndens,nfce=nfce,nfreq=nfreq;,type='landau'
