;Crib sheet for calling cycl_energies_parameter_space.pro

;TIP: **** if you don't know what energy range to request, then don't request one. 
;Program will automatically figure it out. 


;-------------------------------
;Scheme 3 testing
;-------------------------------

dens = 10. 
nres = 1.
pa = 0.
fce = 28.*100.

;scheme 3:  freq vs theta_kb (at constant fce, dens)
cycl_energies_parameter_space,pa,scheme=3,ps=ps,$
  fcev=fce,density_range=density_range,freq_range=freq_range,$
  ndens=ndens,nfce=nfce,nfreq=nfreq,harmonic=1.,type='counterstream';,$
;  minval=1.,maxval=2.,minzval=minzval,maxzval=500.





;--------------------------------------------------
;For Cindy's PSP whistlers
;--------------------------------------------------


ps = 0      ;save to postscript?

pa=15.        ;e- pitch angle
theta_k=60.   ;wave normal angle
density_range = [200,500]   ;cm-3
fce_range = [1000,2000]  ;Hz
freq_range = [50,300]    ;Hz
minval = 0.00    ;minimum energy plotted (keV)
maxval = 10.   ;maximum energy plotted (keV)
minzval = 0.0    ;minimum energy plotted (keV) (FA energy plot)
maxzval = 10.   ;maximum energy plotted (keV)

ndens = 20  &  nfce = 20  &  nfreq = 20  ;number of contours

;                scheme 0:  fce vs density (at constant wave freq)
;                scheme 1:  fce vs freq (at constant density)
;                scheme 2:  density vs freq (at constant fce)
;                scheme 3:  freq vs theta_kb (at constant fce, dens)


cycl_energies_parameter_space,pa,theta_k,scheme=0,ps=ps,$
  density_range=density_range,$
  ndens=ndens,nfce=nfce,nfreq=nfreq,harmonic=1.,type='costream',$
  minval=0,maxval=30,minzval=0,maxzval=30

cycl_energies_parameter_space,pa,theta_k,scheme=1,ps=ps,$
  freq_range=freq_range,$
  ndens=ndens,nfce=nfce,nfreq=nfreq,harmonic=1.,type='costream',$
  minval=1.,maxval=2000.,minzval=minzval,maxzval=2000.

cycl_energies_parameter_space,pa,theta_k,scheme=2,ps=ps,$
  fcev=1500.,density_range=density_range,freq_range=freq_range,$
  ndens=ndens,nfce=nfce,nfreq=nfreq,harmonic=1.,type='costream',$
  minval=0.,maxval=1.,minzval=0.,maxzval=1.

;scheme 3:  freq vs theta_kb (at constant fce, dens)
;Can you run the angle one for n=250/cc and fce of 1100? and n=325 and fce=2000?
cycl_energies_parameter_space,pa,scheme=3,ps=ps,$
  fcev=2000.,densv=325.,freq_range=freq_range,$
  ndens=ndens,nfce=nfce,nfreq=nfreq,harmonic=1.,type='costream',$
  minval=0.,maxval=2.5,minzval=0.,maxzval=2.5







;--------------------------------------------------
;Aaron BARREL Paper 3 calculation for waves near flh
;--------------------------------------------------

;(n~1-5, B~57-80).
;FREQS = 0.2fce and 0.3fce

ps = 0     ;save to postscript?

pa=5.        ;e- pitch angle
theta_k=0.   ;wave normal angle
density_range = [1.,10.]   ;cm-3
fce_range = 28d*[74.,76.]  ;Hz
freq_range = [40,80]    ;Hz (or f/fce if keyword is set)
;freq_range = [0.02,0.04]    ;Hz (or f/fce if keyword is set)
minval = 0.0      ;minimum energy plotted (keV) (total energy plot)
maxval = 100.   ;maximum energy plotted (keV)
minzval = 0.0    ;minimum energy plotted (keV) (FA energy plot)
maxzval = 100.   ;maximum energy plotted (keV)

ndens = 20  &  nfce = 20  &  nfreq = 20  ;number of contours



;scheme 2:  density vs freq (at constant fce)
cycl_energies_parameter_space,pa,theta_k,scheme=2,ps=ps,$
  fcev=fce_range[1],density_range=density_range,freq_range=freq_range,$
  ndens=ndens,nfce=nfce,nfreq=nfreq,harmonic=1.,type='counterstream';,$
;  minval=minval,maxval=maxval,minzval=minzval,maxzval=maxzval,$

  minval = 0.0      ;minimum energy plotted (keV) (total energy plot)
  maxval = 2.   ;maximum energy plotted (keV)
  minzval = 0.0    ;minimum energy plotted (keV) (FA energy plot)
  maxzval = 2.   ;maximum energy plotted (keV)

  cycl_energies_parameter_space,pa,theta_k,scheme=2,ps=ps,$
    fcev=fce_range[1],density_range=density_range,freq_range=freq_range,$
    ndens=ndens,nfce=nfce,nfreq=nfreq,harmonic=1.,type='landau';,$
;    minval=minval,maxval=maxval,minzval=minzval,maxzval=maxzval,$




;--------------------------------------------------
;CINDY AND LINDSAY LWS PROPOSAL (SW ENERGIES)
;--------------------------------------------------

;(n~5-20, B~5-15).
;FREQS = 0.2fce and 0.3fce

ps = 1     ;save to postscript?

pa=45.        ;e- pitch angle
theta_k=45.   ;wave normal angle
density_range = [5.,20.]   ;cm-3
fce_range = 28d*[5.,15.]  ;Hz
freq_range = [0.2*fce_range[0],0.3*fce_range[1]]    ;Hz (or f/fce if keyword is set)
minval = 0.0      ;minimum energy plotted (keV) (total energy plot)
maxval = 100.   ;maximum energy plotted (keV)
minzval = 0.0    ;minimum energy plotted (keV) (FA energy plot)
maxzval = 100.   ;maximum energy plotted (keV)

ndens = 20  &  nfce = 20  &  nfreq = 20  ;number of contours

;Test using direct call
;fpe = 8980.*sqrt(density_range)
;c=3d5
;kvec = sqrt(4*!pi^2*fpe[0]^2*28./(c^2*(420.*cos(0.*!dtor)-28.)))
;evals = cycl_energies(28.,0.,0.,420.,kvec,30.,1.)


;scheme 0:  fce vs density (at constant wave freq)
cycl_energies_parameter_space,pa,theta_k,scheme=0,ps=ps,$
  minval=minval,maxval=maxval,minzval=minzval,maxzval=maxzval,$
  freqv=freq_range[0],fce_range=fce_range,density_range=density_range,$
  ndens=ndens,nfce=nfce,nfreq=nfreq,harmonic=5.


;scheme 1:  fce vs freq (at constant density)
cycl_energies_parameter_space,pa,theta_k,scheme=1,ps=ps,$
  minval=minval,maxval=maxval,minzval=minzval,maxzval=maxzval,$
  densv=5.,fce_range=fce_range,freq_range=freq_range,$
  ndens=ndens,nfce=nfce,nfreq=nfreq,harmonic=5.,/f_fce;,type='landau'


;scheme 2:  density vs freq (at constant fce)
cycl_energies_parameter_space,pa,theta_k,scheme=2,ps=ps,$
  minval=minval,maxval=maxval,minzval=minzval,maxzval=maxzval,$
  fcev=fce_range[1],density_range=density_range,freq_range=freq_range,$
  ndens=ndens,nfce=nfce,nfreq=nfreq,harmonic=5.,/f_fce;,type='landau'






;--------------------------------------------------
;JOHN FOSTER'S (2017) VLF BUBBLE CALCULATION
;--------------------------------------------------


ps = 0      ;save to postscript?

pa=60.        ;e- pitch angle
theta_k=60.   ;wave normal angle
density_range = [30,35]   ;cm-3
fce_range = [38000.,50000.]  ;Hz
freq_range = [21000.,22000.]/4.    ;Hz (or f/fce if keyword is set)
minval = 0.01      ;minimum energy plotted (keV) (total energy plot)
maxval = 5000.   ;maximum energy plotted (keV)
minzval = 0.01    ;minimum energy plotted (keV) (FA energy plot)
maxzval = 5000.   ;maximum energy plotted (keV)

ndens = 20  &  nfce = 20  &  nfreq = 20  ;number of contours

;scheme 0:  fce vs density (at constant wave freq)
cycl_energies_parameter_space,pa,theta_k,scheme=0,ps=ps,$
  minval=minval,maxval=maxval,minzval=minzval,maxzval=maxzval,$
  freqv=21400.,fce_range=fce_range,density_range=density_range,$
  ndens=ndens,nfce=nfce,nfreq=nfreq,harmonic=1


;scheme 1:  fce vs freq (at constant density)
cycl_energies_parameter_space,pa,theta_k,scheme=1,ps=ps,$
  minval=minval,maxval=maxval,minzval=minzval,maxzval=maxzval,$
  densv=30.,fce_range=fce_range,freq_range=freq_range,$
  ndens=ndens,nfce=nfce,nfreq=nfreq;,type='landau'


;scheme 2:  density vs freq (at constant fce)
cycl_energies_parameter_space,pa,theta_k,scheme=2,ps=ps,$
  minval=minval,maxval=maxval,minzval=minzval,maxzval=maxzval,$
  fcev=5000.,density_range=density_range,freq_range=freq_range,$
  ndens=ndens,nfce=nfce,nfreq=nfreq,/f_fce;,type='landau'





;--------------------------------------------------
;TEST OF NEW CYCLOTRON RES CALCULATION (compare to Lorentzen01, plate 7)
;--------------------------------------------------


ps = 0      ;save to postscript?

pa=5.        ;e- pitch angle
theta_k=0.   ;wave normal angle
density_range = [1,10]   ;cm-3
fce_range = [800,25000]  ;Hz
;freq_range = [0,8000]    ;Hz
;freq_range = [0,1]    ;f/fce
freq_range = [0,25000.]    ;f/fce
minval = 1.    ;minimum energy plotted (keV)
maxval = 1000.   ;maximum energy plotted (keV)
minzval = 0.1    ;minimum energy plotted (keV)
maxzval = 100.   ;maximum energy plotted (keV)

ndens = 20  &  nfce = 20  &  nfreq = 20  ;number of values for axes

;scheme 0:  fce vs density (at constant wave freq)
cycl_energies_parameter_space,pa,theta_k,scheme=0,ps=ps,$
  minval=minval,maxval=maxval,minzval=minzval,maxzval=maxzval,$
  freqv=1500.,fce_range=fce_range,density_range=density_range,$
  ndens=ndens,nfce=nfce,nfreq=nfreq,harmonic=1


;scheme 1:  fce vs freq (at constant density)
cycl_energies_parameter_space,pa,theta_k,scheme=1,ps=ps,$
  minval=minval,maxval=maxval,minzval=minzval,maxzval=maxzval,$
  densv=5.,fce_range=fce_range,freq_range=freq_range,$
  ndens=ndens,nfce=nfce,nfreq=nfreq;,type='landau'


;scheme 2:  density vs freq (at constant fce)
cycl_energies_parameter_space,pa,theta_k,scheme=2,ps=ps,$
  minval=minval,maxval=maxval,minzval=minzval,maxzval=maxzval,$
  fcev=5000.,density_range=density_range,freq_range=freq_range,$
  ndens=ndens,nfce=nfce,nfreq=nfreq,/f_fce;,type='landau'



  ;--------------------------------------------------
  ;ALEXA'S ENERGY CALCULATIONS
  ;--------------------------------------------------
  ;Bw = 2.5d-2   ;nT
  ;Bo = 167      ;nT
  ;dens = 12
  ;fpe = 31107.6
  ;fce = 4676.00
  ;f = 0.56*fce

  ps = 1      ;save to postscript?

  pa=5.        ;e- pitch angle
  theta_k=50.   ;wave normal angle
  density_range = [5,15]   ;cm-3
  fce_range = [4650,4700]  ;Hz
  freq_range = [2100,3100]    ;Hz
  ;freq_range = [0.5,0.6]    ;f/fce
  minval = 0.001    ;minimum energy plotted (keV)
  maxval = 10.   ;maximum energy plotted (keV)
  minzval = 0.001    ;minimum energy plotted (keV)
  maxzval = 10.   ;maximum energy plotted (keV)

  ndens = 20  &  nfce = 20  &  nfreq = 20  ;number of contours

  ;scheme 0:  fce vs density (at constant wave freq)
  cycl_energies_parameter_space,pa,theta_k,scheme=0,ps=ps,$
    minval=minval,maxval=maxval,minzval=minzval,maxzval=maxzval,$
    freqv=2618.,fce_range=fce_range,density_range=density_range,$
    ndens=ndens,nfce=nfce,nfreq=nfreq,harmonic=1


  ;scheme 1:  fce vs freq (at constant density)
  cycl_energies_parameter_space,pa,theta_k,scheme=1,ps=ps,$
    minval=minval,maxval=maxval,minzval=minzval,maxzval=maxzval,$
    densv=12.,fce_range=fce_range,freq_range=freq_range,$
    ndens=ndens,nfce=nfce,nfreq=nfreq,harmonic=1;,/f_fce;,type='landau'


  ;scheme 2:  density vs freq (at constant fce)
  cycl_energies_parameter_space,pa,theta_k,scheme=2,ps=ps,$
    minval=minval,maxval=maxval,minzval=minzval,maxzval=maxzval,$
    fcev=4676.,density_range=density_range,freq_range=freq_range,$
    ndens=ndens,nfce=nfce,nfreq=nfreq,/f_fce;,type='landau'



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
;function cold_dispersion,epol=epol,freq=freq,dens=dens,Bo=Bo,H_plus=pH,He_plus=pHe,O_plus=pO

res_angle_max = acos(max(freq_range)/min(fce_range))/!dtor  ;resonance cone angle
res_angle_min = acos(min(freq_range)/max(fce_range))/!dtor  ;resonance cone angle
;--------------------------------------------------

;scheme 0:  fce vs density (at constant wave freq)
theta_k = 0.
cycl_energies_parameter_space,pa,theta_k,scheme=0,ps=ps,$
                              freqv=22361.,fce_range=fce_range,density_range=density_range,$
                              ndens=ndens,nfce=nfce,nfreq=nfreq,type='landau'
;                              minval=minval,maxval=maxval


;scheme 1:  fce vs freq (at constant density)
cycl_energies_parameter_space,pa,theta_k,scheme=1,ps=ps,$
                              densv=2000.,fce_range=fce_range,freq_range=freq_range,$
                              ndens=ndens,nfce=nfce,nfreq=nfreq;,$
;                              minval=minval,maxval=maxval

;scheme 2:  density vs freq (at constant fce)
cycl_energies_parameter_space,pa,theta_k,scheme=2,ps=ps,$
                              fcev=22400.,density_range=density_range,freq_range=freq_range,$
                              ndens=ndens,nfce=nfce,nfreq=nfreq;,$
;                              minval=minval,maxval=maxval
