;+
;*****************************************************************************************
;
;  FUNCTION : cold_plasma_dispersion.pro
;
;****NOTE: Calculated Stix parameters nearly identical to Lynn Wilson's cold_plasma_params.pro 
;	
;  PURPOSE  : Input wave freq, Emax/Eint ratio, dens, Bo.
;			  Returns a structure with the following:
;				 --Useful cold plasma parameters from full cold plasma dispersion
;				   (R,L,D,S,P)
;				 --plasma, cyclotron and related frequencies
;				 --Dispersion relation (index of refr, wavelength, kvect,
;				   polarization, etc...)
;				 --Resonance energies (cyclotron, landau) for range of
;				   pitch angles
;				 --Wave mode from CMA diagram (NOT FINISHED....)
;
;  NOTE: identify wave mode (2 component plasma) by using Stix Figs. 2-2 and 2-1 along with the signs of the seven 
;  quantities R,L,D,S,P,P/S,(RL-PS) 
;
;
;
;  CALLED BY: N/A
;
;  CALLS:
;
;  REQUIRES:
;
;
;  INPUT: 	epol -> Ex/Ey ratio (plane parallel to z-hat (Bo) in Stix coord.
;		freq -> array of freqs (Hz)
;		dens -> density (cm-3)
;		Bo   -> background magnetic field (nT).
;		mr -> mass ratio of proton to electron. Default is 1836. Useful to change if
;				 you'd like to compare results to CMA diagrams with artificial mr.
;		H_plus, He_plus, O_plus -> fractions of H+, He+ and O+. If not set
;				program assumes H+ plasma.
;
;  NOTES:       -Valid for any cold plasma with no density or magnetic field gradients on the scale
;		  of wavelength.
;		  -Whistlers, ion cyclotron, X-mode, O-mode, hydromagnetic modes....
;		  -On the nightside the plasma should be completely
;			H+ by about 1400 km.
;	 	  -Hot plasma effects not important for beta<<1 and for freqs above but not too close
;			to the ion cyclotron freq.
;
;
;
;		  All equation references are from Stix 1992.
;
;		  If we know a-priori the wave normal angle then we'd just use eqn 34 to calculate
;		  the index of refraction. This would be the way to proceed when you have Bw
;		  measurements and can determine theta_kb directly. However, with Ew measurements
;		  we have to find the index of refraction via eqn 42. Then we can figure out
;		  theta_kb, etc. Both methods are exact for a cold plasma.
;
;  ASSUMPTIONS: Assumes plane waves. i.e. that each value of w has a single k and that
;			    there is no spatial dependence of density or magnetic field.
;
;  OUTPUT:  Structure that contains arrays of |k| (1/km), theta_kb (deg), n, wavelength (km),
;			resonance cone angle, phase velocity (km/s), cyclotron res energies (eV)
;
;   CHANGED:  1)  NA [MM/DD/YYYY   v1.0.0]
;
;   CREATED:  11/02/2009
;   CREATED BY:  Aaron W. Breneman
;   LAST MODIFIED:    v1.1.0
;   MODIFIED BY: Hacked from J. Wygant's cold plasma dispersion solver.
;				07/14/2011 - AWB - the speed of light variable is now defined as "c_ms". This
;				 was being overwritten by the cold plasma param "C" causing
;				 the resonance energies to be incorrect. The cyclotron
;				 resonance energies correspond well to Fig2 in Abel and Thorne98
;
; VERSION:
;   $LastChangedBy: $
;   $LastChangedDate: $
;   $LastChangedRevision: $
;   $URL: $
;-
;*****************************************************************************************
;-
;******************************************************************************************

function cold_plasma_dispersion,epol=epol,$
	freq=freq,$
	dens=dens,$
	Bo=Bo,$
	H_plus=pH,$
	He_plus=pHe,$
	O_plus=pO



	;-----------------------------------------------------------------------
	;SETUP AND OTHER STUFF
	;-----------------------------------------------------------------------

	notes = ['freq in Hz',$
	'kmag in 1/km',$
	'theta_kb in deg',$
	'n=index of refraction',$
	'wavelength in km',$
	'resangle in deg',$
	'phasevel in km/sec',$
	'bpol is magnetic polarization ratio',$
	'energy in eV for nonrelativistic energies',$
	'resonance and cutoff values in Hz, calculated without ion contributions',$
	'cBw_to_Ew is the ratio cB/E']

	nfreqs = n_elements(freq)


	;Constants
	c_ms = 2.99792458d8           ; -Speed of light in vacuum (m/s)
	me = 9.1093897d-31            ; -Electron mass (kg)
	mp = 1.6726231d-27            ; -Proton mass (kg)
	e1eV = 1.6e-19
	erg2joule = 1e-7              ; -ergs to joules
	e1eV = 1.6e-19                ;joules in 1 eV
	BnT2BG = 1e-5                 ;nT to Gauss (multiply by this number)
	mr1 = 1836.                   ;mass ratio of H+ (proton) to electron
	mr2 = 1836.*4.                ;              He+ to electron
	mr3 = 1836.*16.               ;              O+  to electron
	if ~keyword_set(pHe) and n_elements(pHe) eq 0. then pHe = 0.
	if ~keyword_set(pO) and n_elements(pO) eq 0. then pO = 0.
	if ~keyword_set(pH) and n_elements(pH) eq 0. then pH = 1.


	;density ratios
	nr1 = pH
	nr2 = pHe
	nr3 = pO

	dtmp = ''
	Botmp = ''
	epoltmp = ''
	freqtmp = ''

	if not keyword_set(dens) then begin
		read,dtmp,prompt='Enter density in cm-3: '
		dens=dtmp
	endif
	if not keyword_set(Bo) then begin
		read,Botmp,prompt='Enter Magnitude of B-field (nT)'
		Bo=Botmp
	endif
	if not keyword_set(epol) then begin
		read,epoltmp,prompt='Enter Ex/Ey ratio (Stix coord)'
		epol=epoltmp
	endif
	if not keyword_set(freq) then begin
		read,freqtmp,prompt='Enter wave freq in plasma frame (Hz)'
		freq=freqtmp
	endif

	w = 2.*!pi*freq

	;consider just H+ (which doesn't have a neutron -- deuterium does)
	w_ce = -28.*Bo*2.*!pi
	w_cH = -1*w_ce/mr1
	w_cHe = -1*w_ce/mr2
	w_cO = -1*w_ce/mr3

	w_pe = 8980.*sqrt(dens)*2.*!pi
	w_pH = w_pe/sqrt(mr1/nr1)
	w_pHe = w_pe/sqrt(mr2/nr2)
	w_pO = w_pe/sqrt(mr3/nr3)


	fce = abs(w_ce/2./!pi)
	fpe = w_pe/2./!pi

	;------------------------------------------------------------------------------------------------
	; We now have all the angular  frequencies in place
	; (electron gyrofreq is negative)
	; -----------------------------------------------------------------------------------------------
	; The next step is to calculate R, L, and P as defined in Stix page 7 (1992 version) equations
	; 20, 21, 22 . There is an ion and electron contribution for each. Remember that in
	; Stix the gyrofrequency has a negative sign for the electrons and a positive
	; sign for the ions (eqn 12). By using this convention, the contributions to
	; R and L from the ions and electrons have the same algebraic form. We use
	; this fact when we sum over the ion and electron contributions to R
	; the ion term is R_i_tmp and the electron term is R_e_tmp
	; the total expression is R= 1-(R_i_tmp + R_e_tmp)
	;------------------------------------------------------------------------------------------------

	R_e_tmp= (w_pe^2)/(w*(w+w_ce))
	R_H_tmp= (w_pH^2)/(w*(w+w_cH))
	R_He_tmp= (w_pHe^2)/(w*(w+w_cHe))
	R_O_tmp= (w_pO^2)/(w*(w+w_cO))
	R = 1 - (R_e_tmp + R_H_tmp + R_He_tmp + R_O_tmp)

	;we now do the same thing for L where L_e_tmp and L_i_tmp are the electron and
	;ion contributions
	L_e_tmp = (w_pe^2)/(w*(w-w_ce))
	L_H_tmp = (w_pH^2)/(w*(w-w_cH))
	L_He_tmp = (w_pHe^2)/(w*(w-w_cHe))
	L_O_tmp = (w_pO^2)/(w*(w-w_cO))
	L = 1 - (L_e_tmp + L_H_tmp + L_He_tmp + L_O_tmp)


	; Now for P which is the along B part
	P_e_tmp = (w_pe^2)/(w^2)
	P_H_tmp = (w_pH^2)/(w^2)
	P_He_tmp = (w_pHe^2)/(w^2)
	P_O_tmp = (w_pO^2)/(w^2)
	P = 1 - (P_e_tmp + P_H_tmp + P_He_tmp + P_O_tmp)


	; we now have R,L, and P
	;Now following Stix Eqn 19, pg 7 we define S and D
	S=0.5*(R+L)
	D=0.5*(R-L)


	;------------------------------------------------------------------------------------------------
	; Now in this first use of these relation we know the polarization of E.
	; For Stix, the Bfield is along z. k is in the xz plane, and there is a
	; relation for the polarization i(Ex/Ey)=(n**2-S)/D (eqn 42 page 10).
	; Now we know what (Ex/Ey)
	; is and so we can figure out what n = kc/w is. Then we can go back and
	; calculate other stuff like the k vector etc.
	;------------------------------------------------------------------------------------------------



	;  if keyword_set(bpol) then begin
	;	;testing
	;	theta2 = acos(bpol/epol)/!dtor
	;
	;	theta2 = acos(1/bpol)
	;	thetadeg2 = theta2/!dtor
	;	epol = bpol/cos(theta2)
	;print,thetadeg2
	;print,epol
	; endif




	n2 = epol*D+S                  ;eqn 42 Stix
	n=sqrt(n2)
	kvect=n*w/(3.0d5)             ;1/km
	wavelength=2.*!pi/kvect       ;km


	;Calculate angle for resonance cone Stix page 12 eq 1-45
	;Note that NaN values mean no resonance cone on dispersion surface
	tantheta2=-P/S
	tantheta=sqrt(tantheta2)
	resangle= (360./6.2830)* atan(tantheta)

	tmp = where(finite(resangle) eq 0.)
	if tmp[0] ne -1 then resangle[tmp] = 90.

	phasevel=3d5/n                ;km/sec


	;---------------------------------------------------------------------------
	; Now that we know n2=n**2 (the index of refraction) we can use equation 36
	; in Stix to get the tan(theta) were theta is the angle between k and B.
	;---------------------------------------------------------------------------

	tan2=-1*P*(n2-R)*(n2-L)/((S*n2-R*L)*(n2-P))
	tanx=sqrt(tan2)
	theta=atan(tanx)
	thetadeg = theta/!dtor
	kvalue = 1000.*sqrt(4.*!pi^2*fpe^2*freq/(c_ms^2*(fce*cos(thetadeg*!dtor)-freq)))


	A = S*sin(thetadeg*!dtor)^2 + P*cos(thetadeg*!dtor)^2
	B = R*L*sin(thetadeg*!dtor)^2 + P*S*(1 + cos(thetadeg*!dtor)^2)
	C = P*R*L
	F2 = (R*L - P*S)^2 * sin(thetadeg*!dtor)^4 + 4*P^2*D^2*cos(thetadeg*!dtor)^2
	F = sqrt(F2)


	;----------------------------------------------------------------------------
	;Frequencies of cutoffs and 1s
	;----------------------------------------------------------------------------

	;R=0 cutoff --> no ion term
	acoeff = 1.
	bcoeff = w_ce
	ccoeff = -1*w_pe^2
	w_ro_e = (-bcoeff + sqrt(bcoeff^2 - 4*acoeff*ccoeff))/(2*acoeff)
	f_Rcutoff_electrons_only = w_ro_e/2./!pi
	;L=0 cutoff --> no ion term, but otherwise no approximations
	acoeff = 1.
	bcoeff = -1*w_ce
	ccoeff = -1*w_pe^2
	w_lo_e = (-bcoeff + sqrt(bcoeff^2 - 4*acoeff*ccoeff))/(2*acoeff)
	f_Lcutoff_electrons_only = w_lo_e/2./!pi
  ;P=0 cutoff 
  f_Pcutoff_electrons_only = w_pe/2./!pi
  f_Pcutoff = sqrt(w_pe^2 + w_pH^2 + w_pHe^2 + w_pO^2)/(2*!pi)

	;Gurnett's version....gives the same value as mine
	;w_lo_e2 = -1*abs(w_ce)/2. +  sqrt((w_ce/2.)^2 + w_pe^2)
	;f_Lcutoff2 = w_lo_e2/2./!pi


  ;----------------------------------------------------
  ;Solve for f_Rcutoff and f_Lcutoff including ion terms 
  ;----------------------------------------------------
  nfreqs = 1000

  ;f_Rcutoff
  ;Use the value of f_Lcutoff (e- term only) to narrow the range of frequencies to search for solution
  deltav = 0.1*f_Rcutoff_electrons_only
  maxfreq = 2*!pi * (f_Rcutoff_electrons_only + deltav)
  minfreq = 2*!pi * (f_Rcutoff_electrons_only - deltav)
  omegas = indgen(nfreqs)*maxfreq / (nfreqs-1) + minfreq ;2*!pi*f_Rcutoff/2

  lhs = omegas
  rhs = w_pe^2/(omegas+w_ce) + w_pH^2/(omegas+w_cH) + w_pHe^2/(omegas+w_cHe) + w_pO^2/(omegas+w_cO)
  diff = rhs - lhs
  goo = min(abs(diff),wh)
  f_Rcutoff = omegas[wh]/(2.*!pi)

  ;!p.multi = [0,0,2]
  ;plot,omegas,lhs
  ;oplot,omegas,rhs,color=240
  ;plot,omegas,diff


  ;f_Lcutoff
  ;Use the value of f_Lcutoff (e- term only) to narrow the range of frequencies to search for solution
  deltav = 0.1*f_Lcutoff_electrons_only
  maxfreq = 2*!pi * (f_Lcutoff_electrons_only + deltav)
  minfreq = 2*!pi * (f_Lcutoff_electrons_only - deltav)
  omegas = indgen(nfreqs)*maxfreq / (nfreqs-1) + minfreq

  lhs = omegas 
  rhs = w_pe^2/(omegas-w_ce) + w_pH^2/(omegas-w_cH) + w_pHe^2/(omegas-w_cHe) + w_pO^2/(omegas-w_cO)
  diff = rhs - lhs
  goo = min(abs(diff),wh)
  f_Lcutoff = omegas[wh]/(2.*!pi)

  ;!p.multi = [0,0,2]
  ;plot,omegas,lhs
  ;oplot,omegas,rhs,color=240
  ;plot,omegas,diff



	s_o = sqrt(w_ce^2 + w_pe^2)/2./!pi



	;; --------------------------------------------------------------------------
	;; Find relativistic cyclotron resonance energies for a range of pitch angles
	;; --------------------------------------------------------------------------


	;; test pitch angles
	pa = [0.,10.,20.,30.,40.,50.,60.,70.,80.,89.]
	freqtmp = replicate(freq,10)
	kvecttmp = replicate(kvect,10)
	fcetmp = replicate(fce,10)
	thetadegtmp = replicate(thetadeg,10)
	nres = replicate(1,10)



	;cyclotron and costream cyclotron
	evals = cycl_energies(freqtmp,thetadegtmp,pa,fcetmp,kvecttmp,dens,nres)
	vz_cycl = evals.vz_cycl_counterstream
	vz_costream = evals.vz_cycl_costream
	vz_landau = evals.vz_landau
	vtots_cycl = evals.vtotal_cycl_counterstream
	vtots_costream = evals.vtotal_cycl_costream
	vtots_landau = evals.vtotal_landau
	Ez_cycl = evals.Ez_cycl_counterstream
	Ez_costream = evals.Ez_cycl_costream
	Ez_landau = evals.Ez_landau
	Etots_cycl = evals.E_cycl_counterstream
	Etots_costream = evals.E_cycl_costream
	Etots_landau = evals.Ez_landau
	vc_ratio_cycl = 1000.*vtots_cycl/c_ms
	vc_ratio_costream = 1000.*vtots_costream/c_ms
	vc_ratio_landau = 1000.*vtots_landau/c_ms


	notes2 = ['values are for each pitch angle','resonance energies in keV','velocities in km/s','vc_ratio_xxx is the ratio v/c']
	cyclo_counterstream_res = {vz:vz_cycl,vtots:vtots_cycl,vc_ratio:vc_ratio_cycl,Ez:Ez_cycl,Etots:Etots_cycl,pitch_angles:pa,notes:notes2}
	cyclo_costream_res = {vz:vz_costream,vtots:vtots_costream,vc_ratio:vc_ratio_costream,Ez:Ez_costream,Etots:Etots_costream,pitch_angles:pa,notes:notes2}
	landau_res = {vz:vz_landau,vtots:vtots_landau,vc_ratio:vc_ratio_landau,Ez:Ez_landau,Etots:Etots_landau,pitch_angles:pa,notes:notes2}



	;; ------------------------------------------------------------------------------------------------
	;; VARIOUS PLASMA QUANTITIES
	;; ------------------------------------------------------------------------------------------------


	fpe = w_pe/2./!pi
	fpH = w_pH/2./!pi
	fpHe = w_pHe/2./!pi
	fpO= w_pO/2./!pi

  ;fce = abs(w_ce/2./!pi)
  fce = w_ce/2./!pi
	fcH = w_cH/2./!pi
	fcHe = w_cHe/2./!pi
	fcO = w_cO/2./!pi

  fuh = sqrt(fce^2 + fpe^2)

	meff = 1/(pH/1. + pHe/4. + pO/16.)
	flhr = sqrt(fpe^2*fce^2/(fpe^2 + fce^2)*(1/1836./meff))


	e_inertial = 5.31e5/sqrt(dens)/100./1000.        ;e- inertial length (skin depth) (km)
	ion_inertial = 2.28e7/sqrt(dens)/100./1000.      ;H+ inertial length (km)
	VA = 2.18e11*Bo*BnT2BG/sqrt(dens)/100./1000.     ;H+ Alfven vel (km/s)
	VA2 = fcH/fpH*3e5                                ;Alternate H+ Alfven vel (km/s)
	VAe = VA * sqrt(1836.)                           ;electron Alfven vel (km/sec)
	mepp = (Bo*BnT2BG)^2/8./!pi/dens * erg2joule /e1ev ;characteristic magnetic energy per particle (eV)



	;X=0 cutoff
	goo1 = 4*(fpe/fce)^2
	fx = (fce/2.)*(1 + sqrt(1 + goo1))
  ;Z=0 cutoff 
  fz = fx - fce
  ;Z=infinity resonance (upper oblique resonance)
  goo1 = sqrt(fuh^4 - 4*(fce*fpe*cos(theta))^2)
  fzi = (1/sqrt(2.))*sqrt(fuh^2 + goo1)


  ;Return values of the cold plasma parameters. Note that these can be used 
  ;along with Stix Fig. 2-2 (and 2-1) to determine which region (1-13) a wave falls within. 
  cp_params = {R:R,L:L,D:D,S:S,P:P,$
    P_S:P/S,$
    RLminusPS:(R*L - P*S),$
    A:A,B:B,C:C,F:F,$
    RL_S:R*L/S,$
    R_cutoff:f_Rcutoff,$
    R_cutoff_electrons_only:f_Rcutoff_electrons_only,$
    L_cutoff:f_Lcutoff,$
    L_cutoff_electrons_only:f_Lcutoff_electrons_only,$
    P_cutoff:f_Pcutoff,$
    P_cutoff_electrons_only:f_Pcutoff_electrons_only,$
    X_cutoff:fx,$
    Z_cutoff:fz,$   
    S_resonance:s_o,$
    R_resonance:abs(w_ce)/2./!pi,$
    Z_oblique_resonance:fzi,$
    fuh_resonance:fuh}



	;;------------------------------------------------------------------------------------------------
	;; now we calculate the magnetic field polarization in the plane perpendicular
	;; to the k vector. We use the fact that the plane of the electric
	;; field variations is perpendicular to B and we rotate to a plane perpendicular
	;; to K (by theta) and use the phase velocity measurement in the relation between; E and B.
	;; This give a simple relation between electric and magnetic field polarization:
	;;------------------------------------------------------------------------------------------------

	bpol=epol*cos(theta)


	;;Calculate Ex,Ey,Ez
	;;Since I only know epol=Ex/Ey I'll assume that Ey=1 and Ex=epol. Ez will then be related to these.
	Ex2Ey = 1/epol   ;use the inverse value b/c Ex should be semimajor axis and Ey semiminor
	Ez2Ex = (n^2*cos(theta)*sin(theta))/(n^2*sin(theta)^2-P)


	;------------------------------------------------------------------------------
	;Find cB/E for cold plasma waves
	;To convert to B in nT use 
	;BnT = 1d9*(cB/E)*EmV/m/3d8
	;------------------------------------------------------------------------------

	;Oblique waves Eqn 2 and 3 Agapitov 2013 ("Statistics of whistler mode waves...")
	num = sin(theta)*(n^2 - P)
	den1 = (n^2*sin(theta)^2 - P)^2
	den2 = (D^2*(n^2*sin(theta)^2 - P)^2) / ((n^2 - S^2)^2)
	den3 = n^4 * sin(theta)^2 * cos(theta)^2
	cosbeta = num/sqrt(den1 + den2 + den3)
	betatmp = acos(cosbeta)  ;angle (radians) b/t Ew and k
	cBw_to_Ew_agapitov = n*sin(betatmp)


	;From Hartley16 Eqn1 (parallel whistler mode only)
	;(Bw in Tesla and E in V/m)
	t1 = (fpe^2)/(freq*(freq-fce))
	cBw_to_Ew_hartley = sqrt(1. - t1)


	;Oblique waves from Hartley16 Eqn2 (doi:10.1002/2016JA022501)
	;(Bw in Tesla and E in V/m)
	t1 = (D/(S-n^2))^2
	t2 = (P - n^2*sin(theta)^2)^2
	t3 = P^2*cos(theta)^2
	t4 = t2 
	t5 = t1 + 1.
	t6 = (n^2*cos(theta)*sin(theta))^2
	num = t1*t2 + t3 
	den = t4*t5 + t6
	cBw_to_Ew_hartley2 = n*sqrt(num/den) 

	;------------------------------------------------------------------------------

	struct =   {freq:freq, $
	kmag:kvect, $
	theta_kb:thetadeg, $
	n:n, $
	dens:dens,$
	Bo:Bo,$
	wavelength:wavelength, $
	resangle:resangle, $
	phasevel:phasevel, $
	bpol:bpol, $
	cyclo_counterstream_res:cyclo_counterstream_res,$
	cyclo_costream_res:cyclo_costream_res,$
	landau_res:landau_res,$
	energy_notes:notes2,$
	fpe:fpe,$
	fpH:fpH,$
	fpHe:fpHe,$
	fpO:fpO,$
	fce:fce,$
	fcH:fcH,$
	fcHe:fcHe,$
	fcO:fcO,$
	flhr:flhr,$
	e_inertial:e_inertial,$
	ion_inertial:ion_inertial,$
	VA_ion:VA,$
	VA_electron:VAe,$
	mepp:mepp,$
	Ex2Ey:epol,$
	Ez2Ex:Ez2Ex,$
	cB2E_obliquepropagation_v1:cBw_to_Ew_agapitov,$
	cB2E_obliquepropagation_v2:cBw_to_Ew_hartley2,$
	cB2E_parallelpropagation:cBw_to_Ew_hartley,$
	betaangle_E_k:betatmp,$
	cp_params:cp_params,$
	notes:notes}

	return,struct
end
