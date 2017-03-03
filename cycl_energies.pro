;+
; NAME: cycl_energies
; SYNTAX:
; PURPOSE: Calculates relativistic cyclotron (and Landau) resonance energies (keV) for input arrays
;          of wave freq, fce, |k|, nres, etc... Solves the quadratic eqn for parallel
;          e- velocity that comes from Doppler-shifted cyclotron resonance condition.
;          Note that user needs to choose their own dispersion relation to determine |k|
; INPUT: (see cycl_energies_test.pro)
; OUTPUT: Total electron energy (keV) and field-aligned energy for both +/- resonance, as
;         well as respective velocities in km/s
;
; NOTE: This has been tested against Lorentzen01 Plate 7 showing close agreement (see cycl_energies_test.pro)
;
; KEYWORDS:
;          pitchangle -> e- pitch angle (deg)  [array]
;          theta_kb -> wave normal angle (deg) [array]
;          freq -> wave frequency (Hz)         [array]
;          fce -> electron cyclotron freq (Hz) [array]
;          nres -> absolute value of resonance harmonic (integer) [array]
;          kmag -> |k| (1/km) [array]
;
;
; REQUIRES:
; HISTORY: Written by AWB 2016-11-29
; VERSION:
;-



function cycl_energies,freq,theta_kb,pitchangle,fce,kmag,dens,nres

  n = n_elements(theta_kb[*,0])
  nrays = n_elements(freq[0,*])

  if nrays eq 0 then nrays = 1
  ;returned values
  vz_cycl = fltarr(n,nrays)
  vtots_cycl = fltarr(n,nrays)
  Ez_cycl = fltarr(n,nrays)
  Etots_cycl = fltarr(n,nrays)
  vz_landau = fltarr(n,nrays)
  E_landau = fltarr(n,nrays)
  vzp = fltarr(n,nrays)
  vzm = fltarr(n,nrays)
  vzl = fltarr(n,nrays)
  vtots_p = fltarr(n,nrays)
  vtots_m = fltarr(n,nrays)
  vtots_l = fltarr(n,nrays)
  etots_cyclp = fltarr(n,nrays)
  etots_cyclm = fltarr(n,nrays)
  etots_landau = fltarr(n,nrays)
  ez_cyclp = fltarr(n,nrays)
  ez_cyclm = fltarr(n,nrays)
  ez_landau = fltarr(n,nrays)



  me_ev  = 0.51d6       ; -Electron mass in eV/c^2
  c2t = 1  ;speed of light squared. Cancels out from me_ev
  fpe = 8980.*sqrt(dens)
  w = 2d0*!dpi*freq
  wce = 2d0*!dpi*fce
  cos_t = cos(theta_kb*!dtor)
  cos_pa = cos(pitchangle*!dtor)

  c_ms = 2.99792458d8      ; -Speed of light in vacuum (m/s)
  c_kms = c_ms/1000d0


  for qq=0,nrays-1 do begin

  ;set up quadratic eqn
  a1 = nres*wce[*,qq]/c_kms/cos_pa[qq]
  a1 = a1^2
  a2 = kmag[*,qq]*cos_t[*,qq]
  a2 = a2^2
  a = a1 + a2


  b = 2.*w[qq]*kmag[*,qq]*cos_t[*,qq]

  c1 = w[qq]^2
  c1 = replicate(c1,n)
  c2 = nres*wce[*,qq]
  c2 = c2^2
  c = c1 - c2


  ;plus solution to quadratic eqn
  vzp[*,qq] = (-b + sqrt(b^2 - 4.*a*c))/(2.*a)
  ;minus solution to quadratic eqn
  vzm[*,qq] = (-b - sqrt(b^2 - 4.*a*c))/(2.*a)

  vzp[*,qq] = 1000.*abs(vzp[*,qq])   ;m/s
  vzm[*,qq] = 1000.*abs(vzm[*,qq])   ;m/s

  vtots_p[*,qq] = vzp[*,qq]/cos_pa[qq]
  vtots_m[*,qq] = vzm[*,qq]/cos_pa[qq]

  ;; Relativistic energy in keV (e.g. p37 in "Modern Physics, 2nd edition")
  Etots_cyclp[*,qq] = (0.511d6/sqrt(1-(vtots_p[*,qq]^2/c_ms^2)) - 0.511d6)/1000.
  Ez_cyclp[*,qq] = (0.511d6/sqrt(1-(vzp[*,qq]^2/c_ms^2)) - 0.511d6)/1000.

  Etots_cyclm[*,qq] = (0.511d6/sqrt(1-(vtots_m[*,qq]^2/c_ms^2)) - 0.511d6)/1000.
  Ez_cyclm[*,qq] = (0.511d6/sqrt(1-(vzm[*,qq]^2/c_ms^2)) - 0.511d6)/1000.


  ;-----------------------------------------
  ;Same calculation but for Landau resonance (note: the +/- solutions are the same,
  ;so I'll just use the +)
  ;-----------------------------------------


  Ez_landau[*,qq] = 0.5*me_ev*c2t*((fce[*,qq]^2)/(fpe[*,qq]^2))*(freq[*,qq]/fce[*,qq])*(cos(theta_kb[*,qq]*!dtor)-(freq[*,qq]/fce[*,qq]))*(1/(cos(theta_kb[*,qq]*!dtor)^2))/1000.

  fac1 = 0.511d6/(0.511d6 + 1000.*Ez_landau[*,qq])
  fac1 = 1. - fac1^2
  vz2 = c_ms^2*fac1
  vz_landau[*,qq] = sqrt(vz2)

;  ;set up quadratic eqn
;  a1 = 0d
;  a2 = double(kmag[*,qq]*cos_t[*,qq])
;  a2 = a2^2
;  a = a1 + a2
;
;  b = 2d*w[qq]*kmag[*,qq]*cos_t[*,qq]
;
;  c1 = w[qq]^2
;  c2 = replicate(0d,n)
;  c = c1 - c2
;
;
;  ;plus solution to quadratic eqn
;  vzl[*,qq] = (-1d*b + sqrt(b^2 - 4d*a*c))/(2d*a)
;  vzl[*,qq] = 1000d*abs(vzl[*,qq])   ;m/s
;  vtots_l[*,qq] = vzl[*,qq]/cos_pa[qq]
;
;  ;; Relativistic energy in keV (e.g. p37 in "Modern Physics, 2nd edition")
;  Etots_landau[*,qq] = (0.511d6/sqrt(1-(vtots_l[*,qq]^2/c_ms^2)) - 0.511d6)/1000d
;  Ez_landau[*,qq] = (0.511d6/sqrt(1-(vzl[*,qq]^2/c_ms^2)) - 0.511d6)/1000d


endfor


  return,{$
  E_cycl_counterstream:Etots_cyclp,$           ;keV
  E_cycl_costream:Etots_cyclm,$
;  E_landau:Etots_landau,$
  Ez_cycl_counterstream:Ez_cyclp,$
  Ez_cycl_costream:Ez_cyclm,$
  Ez_landau:Ez_landau,$
  vtotal_cycl_counterstream:vtots_p/1000.,$    ;km/s
  vtotal_cycl_costream:vtots_m/1000.,$
  vtotal_landau:vtots_l/1000.,$
  vz_cycl_counterstream:vzp/1000.,$
  vz_cycl_costream:vzm/1000.,$
  vz_landau:vz_landau/1000.}

end
