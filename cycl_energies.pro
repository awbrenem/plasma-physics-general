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



function cycl_energies,freq,theta_kb,pitchangle,fce,kmag,nres

  n = n_elements(freq)

  ;returned values
  vz_cycl = fltarr(n)
  vtots_cycl = fltarr(n)
  Ez_cycl = fltarr(n)
  Etots_cycl = fltarr(n)
  v_landau = fltarr(n)
  E_landau = fltarr(n)


  w = 2d0*!dpi*freq
  wce = 2d0*!dpi*fce
  cos_t = cos(theta_kb*!dtor)
  cos_pa = cos(pitchangle*!dtor)

  c_ms = 2.99792458d8      ; -Speed of light in vacuum (m/s)
  c_kms = c_ms/1000d0


  ;set up quadratic eqn
  a1 = nres*wce/c_kms/cos_pa
  a1 = a1^2
  a2 = kmag*cos_t
  a2 = a2^2
  a = a1 + a2


  b = 2.*w*kmag*cos_t

  c1 = w^2
  c2 = nres*wce
  c2 = c2^2
  c = c1 - c2


  ;plus solution to quadratic eqn
  vzp = (-b + sqrt(b^2 - 4.*a*c))/(2.*a)
  ;minus solution to quadratic eqn
  vzm = (-b - sqrt(b^2 - 4.*a*c))/(2.*a)

  vzp = 1000.*abs(vzp)   ;m/s
  vzm = 1000.*abs(vzm)   ;m/s

  vtots_p = vzp/cos_pa
  vtots_m = vzm/cos_pa

  ;; Relativistic energy in keV (e.g. p37 in "Modern Physics, 2nd edition")
  Etots_cyclp = (0.511d6/sqrt(1-(vtots_p^2/c_ms^2)) - 0.511d6)/1000.
  Ez_cyclp = (0.511d6/sqrt(1-(vzp^2/c_ms^2)) - 0.511d6)/1000.

  Etots_cyclm = (0.511d6/sqrt(1-(vtots_m^2/c_ms^2)) - 0.511d6)/1000.
  Ez_cyclm = (0.511d6/sqrt(1-(vzm^2/c_ms^2)) - 0.511d6)/1000.


  ;-----------------------------------------
  ;Same calculation but for Landau resonance (note: the +/- solutions are the same,
  ;so I'll just use the +)
  ;-----------------------------------------


  ;set up quadratic eqn
  a1 = 0.
  a2 = kmag*cos_t
  a2 = a2^2
  a = a1 + a2

  b = 2.*w*kmag*cos_t

  c1 = w^2
  c2 = 0.
  c = c1 - c2


  ;plus solution to quadratic eqn
  vzl = (-b + sqrt(b^2 - 4.*a*c))/(2.*a)
  vzl = 1000.*abs(vzl)   ;m/s
  vtots_l = vzl/cos_pa

  ;; Relativistic energy in keV (e.g. p37 in "Modern Physics, 2nd edition")
  Etots_landau = (0.511d6/sqrt(1-(vtots_l^2/c_ms^2)) - 0.511d6)/1000.
  Ez_landau = (0.511d6/sqrt(1-(vzl^2/c_ms^2)) - 0.511d6)/1000.


  return,{$
  E_cycl_normal:Etots_cyclp,$           ;keV
  E_cycl_anom:Etots_cyclm,$
  E_landau:Etots_landau,$
  Ez_cycl_normal:Ez_cyclp,$
  Ez_cycl_anom:Ez_cyclm,$
  Ez_landau:Ez_landau,$
  vtotal_cycl_normal:vtots_p/1000.,$    ;km/s
  vtotal_cycl_anom:vtots_m/1000.,$
  vtotal_landau:vtots_l/1000.,$
  vz_cycl_normal:vzp/1000.,$
  vz_cycl_anom:vzm/1000.,$
  vz_landau:vzl/1000.}

end
