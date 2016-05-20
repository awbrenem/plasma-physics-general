;Various plasma parameters, etc...


;--------------------------
;PLASMA PARAMETERS
;SHOCK QUANTITIES
;Plasma drift values
;Instability calculations
;Dipole field stuff (see dipole.pro)
;EARTH'S DENSITY PROFILE
;DISPERSION RELATIONS
	;WHISTLERS (e- only)
	;quasi-electrostatic whistlers (Bell and Ngo, 1990)
	;EMIC (parallel only)
	;EIC
	;Alfven waves (ideal, parallel only)
	;Ion acoustic waves
;Doppler shift plot
;RELATIVISTIC CYCLOTRON RESONANCE ENERGIES (NEED TO CHOOSE PARTICLE PITCH ANGLE)
;TRANSVERSE RESONANCE ENERGY
;NONRELATIVISTIC CYCLOTRON RES ENERGIES (NOTE: FOR ANYTHING IN THE 
;VICINITY OF 100 KEV YOU DEFINITELY NEED TO USE RELATIVISTIC EQNS ABOVE
;POYNTING FLUX CALCULATIONS (MONOCHROMATIC, PLANE WAVE)
;Radiation calculations
;E and B amplitudes
;ATMOSPHERIC STUFF.....
		;LIGHTNING CALCULATIONS
		;COLLISION CALCULATIONS
		;MINIMUM PENETRATION FREQS OF IONOSPHERE



;-------------------------------
;THINGS TO ADD: note that the units should be the same as those in plastic files: n=cm-3, Tp=K, Bt=nT
 ;entropy=log((tp**1.5)/np) ) - 6*log(10);
    ;beta=np*10**(-5)*(469895.8 + 4.048*tp) / bt**2;
;-------------------------------



;--------------------------
;INPUT VALUES
;--------------------------
B = 3700. ;nT
n = 44000. ;cm-3
Z = 1    ;charge state (number of unmatched e-)
muu = 1  ;mi/mp


;Define polytropic index

;Isothermal process (gamma = 1)
;Adiabatic (isentropic) process
	;gamma = (alpha + 1)/alpha, where alpha is the number of degrees of freedom divided by 2
	;So, for full 3D degrees of freedom for monotonic gas gamma=5/3
	;For diatomic gas gamma=7/5
	;For monotonic gas undergoing 1D compression (like sound wave) alpha = 1/2, gamma = 3

gamai = 3.
gamae = 1.


Te = 0.1 ;eV
Ti = 1. ;eV
f=99. ;Hz  wave freq in sc frame
Vsw = 400  ;solar wind vel (km/sec)
theta_k = 0.  ;degrees
thetas = indgen(90)

;-------------------------
;CONSTANTS (SI units)
;-------------------------
epo    = 8.854187817d-12   ; -Permittivity of free space (F/m)
muo    = 4d0*!DPI*1d-7     ; -Permeability of free space (N/A^2 or H/m)
me     = 9.1093897d-31     ; -Electron mass (kg)
mp     = 1.6726231d-27
ma     = 4d0*(mp + me)     ; -Alpha-Particle mass (kg)
qq     = 1.60217733d-19    ; -Fundamental charge (C)
kB     = 1.380658d-23      ; -Boltzmann Constant (J/K)
c      = 2.99792458d8      ; -Speed of light in vacuum (m/s)
h	   = 4.135e-15		   ; -Planck constant (eV s)
;-------------------------
;CONSTANTS (Gaussian (cgs) units)
;-------------------------
;epo    = 1.
;muo    = 1.
;me     = 9.1093897d-31     ; -Electron mass (kg)
;mp     = 1.6726231d-27     ; -Proton mass (kg)
;ma     = 4d0*(mp + me)     ; -Alpha-Particle mass (kg)
;qq     = 4.8032d-10        ; -statcoulombs
;kB     = 1.3807d-16        ; -Boltzmann Constant (erg/K)
;-------------------------
;CONVERSION FACTORS
;-------------------------
K_eV   = 1.160474d4        ; -Conversion = degree Kelvin/eV
mm     = [-1d0,0d0,1d0]    ; -[Normal Cyclotron, Landau, Anomalous Cyclotron]
me_ev  = 0.51d6       ; -Electron mass in eV/c^2
me_3dp = me_ev/(c*1d-3)^2  ; -Electron mass in eV/(km/s)^2
mp_ev  = 938.27231d6       ; -Proton mass in eV/c^2
erg2joule = 1e-7           ; -ergs to joules
e1eV = 1.6e-19              ;joules in 1 eV
BnT2BG = 1e-5               ;nT to Gauss (multiply by this number)
Hplus2e_mass = mp/me       ;H+ to electron mass ratio


;-----------------------
;YARDSTICKS
;-------------------------

;Transform from E, B to velocity
E=VxB, where E(mV/m) = V(km/s) * B(nT)/1000
;Test this. We see at RBSP perigee ~200-300 mV/m Vsc x B fields
;Vsc ~ 10 km/s,   B = 3d4 nT
;E = 10.*3d4/1000. = 300 nT




;-------------------------
;USEFUL UNITS
;-------------------------
;-----MKS--------------------------------------------CGS
Tesla (T) = Wb m-2 = kg s-2 A-1.....................Gauss (G)
Volt (V) = W/A = kg m2 s-3 A-1......................statvolt
Watt (W) = J s-1 = kg m-2 s-3.......................erg
Farad (F) = A s V-1 = s4 A2 kg-1 m-2



;--------------------------
;CALCULATED VALUES
	;note that I follow the convention of the NRL formulary:	
	;Input density in cm-3, B in nT (I convert to gauss in the equations) and temp in eV
;--------------------------

fpe = 8980.*sqrt(n)                                  ;electron plasma freq (Hz)
fpi = fpe*Z/(sqrt(muu)*43.)                          ;Ion plasma freq (Hz)
fce = 28.*B                                          ;electron cyclotron freq (Hz)
fci = fce*Z/(1836.*muu)                              ;ion cyclotron freq (Hz)
debye = 7.43e2*sqrt(Te)/sqrt(n)/100./1000.           ;Debye length (km) (for e- and ions)
flhr = sqrt(fpi*fpi*fce*fci)/sqrt(fce*fci+(fpi^2))   ;Lower Hybrid Res freq (Hz)
flhr2 = sqrt(abs(fce)*fci)							 ;Lower Hybrid Res freq (Hz) --> high dens limit - no plasma freq terms (Gurnett 4.4.51)
;;flhr3 = fpi									     ;Lower Hybrid Res freq (Hz) --> low dens limit 
fuh = sqrt(fpe^2 + fce^2)							 ;Upper hybrid res freq (Hz)
fx = fce/2. + 0.5*sqrt(fce^2 + 4*fpe^2)              ;x-mode frequency (Hz)
pe = 2.38*sqrt(Te)/(B*BnT2BG)/100./1000.             ;e- thermal gyroradius (km)
pi = 102*sqrt(muu)*sqrt(Ti)/(B*BnT2BG)/100./1000./Z  ;ion thermal gyroradius (km)     
VA = 2.18e11*B*BnT2BG/sqrt(n)/100./1000.             ;H+ Alfven vel (km/s)
VA2 = fci/fpi*3e5                                    ;Alternate H+ Alfven vel (km/s)
VAe = VA * sqrt(1836.)                               ;electron Alfven vel (km/sec)
e_inertial = 5.31/sqrt(n)                            ;e- inertial length (skin depth) (km) - VERIFIED
e_inertial2 = Vae/2./!pi/fce					     ;alternate formulation (km)
ion_inertial = 2.28e7*sqrt(muu)/(Z*sqrt(n))/100./1000. ;Ion inertial length (km)
ion_inertial2 = Va/2./!pi/fci                        ;alternate formulation (km/s) (ONLY CORRECT FOR H+)
Ve = 4.19e7*sqrt(Te)/100./1000.                      ;e- thermal vel (km/s)
Vi = 9.79e5*sqrt(Ti/muu)/100./1000.                  ;Ion thermal vel (km/s)
Cs_i = 9.79e5*sqrt(gamai*Z*Te/muu)/100./1000.        ;Ion sound vel (km/s)
Cs_i2 = sqrt((Te+3*Ti)/mp_ev)
;vel_ion = 1.4d*sqrt(ion_eV)                          ;velocity of ion in terms of ion energy in eV (km/s)
;vel_elec = 5.9d2*sqrt(elec_eV)                       ; ''         electron ''        ''            (km/s)
mepp = (B*BnT2BG)^2/8./!pi/n * erg2joule /e1ev       ;characteristic magnetic energy per particle (eV) - which particles will cyclotron interact
beta_p = 4.03e-11*n*Ti/((B*BnT2BG)^2)                ;Proton beta
beta_e = 4.03e-11*n*Te/((B*BnT2BG)^2)                ;electron beta
beta_t = 4.03e-11*n*(Ti+Te)/((B*BnT2BG)^2)     	     ;total beta
pmag = 0.3979*B^2                                    ;magnetic pressure (pPa) (B must be in nT)

;I DON'T THINK THIS IS RIGHT. THE VALUE IS MANY ORDERS OF MAG LESS THAN PMAG
pkin = (100.^3) *n*1.38e-23*Ti*10^12                 ;kinetic pressure in pPa (B in nT, n in cm-3) 
;-------------
;THIS VALUE TENDS TO BE CLOSE TO PMAG, BUT IT INCLUDES 4% ALPHA PARTICLE CORRECTION
pkin = nU*1e-3*(1869.6979 + 0.0161082*TiU) ;kinetic pressure


gendrin_angle = acos(2.*f/fce)/!dtor                 ;whistler mode Gendrin angle
res_angle = acos(f/fce)/!dtor					     ;resonance cone angle


;pkin = n*1e-3*(1869.6979 + 0.0161082*Ti) ;kinetic pressure

print,'fpe = ' + strtrim(fpe,2) + ' Hz'
print,'fpi = ' + strtrim(fpi,2) + ' Hz'
print,'fce = ' + strtrim(fce,2) + ' Hz'
print,'fci = ' + strtrim(fci,2) + ' Hz'
print,'flhr = ' + strtrim(flhr,2) + ' Hz'
print,'fuh = ' + strtrim(fuh,2) + ' Hz'
print,'pe = ' + strtrim(pe,2) + ' km'
print,'pi= ' + strtrim(pi,2) + ' km'
print,'e- inertial length = ' + strtrim(e_inertial,2) + ' km'
print,'ion inertial length = ' + strtrim(ion_inertial,2) + ' km'
print,'Debye length = ' + strtrim(Debye,2) + ' km'
print,'Ve = ' + strtrim(Ve,2) + ' km/sec'
print,'Vi = ' + strtrim(Vi,2) + ' km/sec'
print,'VA = ' + strtrim(VA,2) + ' km/sec'
print,'VA2 = ' + strtrim(VA2,2) + ' km/sec'
print,'VAe = ' + strtrim(vae,2) + ' km/sec'
print,'Cs_i = ' + strtrim(Cs_i,2) + ' km/sec'
print,'beta_e = ' + strtrim(beta_e)
print,'beta_p = ' + strtrim(beta_p)
print,'beta_t = ' + strtrim(beta_t)




;---------------------------------------------------------------------------------
;SHOCK QUANTITIES
;---------------------------------------------------------------------------------

	
	;density (cm-3)
	;U = upstream
	;D = downstream
	nU = 1
	nD = 2
	;temp (eV)
	TiU = 1
	TiD = 2
	TeU = 3
	TeD = 4
	;velocity and shock normal unit vector (use the same coord for both)
	vswU = [1,2,3]
	vswD = [1,2,3]
	shnorm = [1,1,1]   
	
	;Magnetic field (nT)
	BmagU = 20
	BmagD = 30
	
	;project into shock normal direction
	VprojU = vswU*shnorm
	VprojD = vswD*shnorm
	
	VprojU_mag = reform(sqrt(VprojU[0]^2 + VprojU[1]^2 + VprojU[2]^2))
	VprojD_mag = reform(sqrt(VprojD[0]^2 + VprojD[1]^2 + VprojD[2]^2))
	
	gamae = 1.  ;isothermal
	gamai = 3.  ;adiabatic
	
	
	CsU = 9.79d5*sqrt(gamai*TeU)/100./1000.   ;km/s
	VAU = 2.18d11*BmagU*bnt2bg/sqrt(nU)/100./1000.    ;km/s
	CmsU = sqrt(VAU^2 + CsU^2)
	
	CsD = 9.79d5*sqrt(gamai*TeD)/100./1000.   ;km/s
	VAD = 2.18d11*BmagD*bnt2bg/sqrt(nD)/100./1000.    ;km/s
	CmsD = sqrt(VAD^2 + CsD^2)
	
	
	;------------------------------------------
	print,'Upstream quantities:'
	print,'Magnetosonic mach # = ',VprojU_mag/CmsU
	print,'Alfven mach # = ',VprojU_mag/VAU
	
	print,'Downstream quantities:'
	print,'Magnetosonic mach # = ',VprojD_mag/CmsD
	print,'Alfven mach # = ',VprojD_mag/VAD
	
	print,'For the polytropic index I used gamma_e = ',gamae
	print,'...and gamma_i = ',gamai
	
	print,'Upstream ion beta',4.03e-11*nU*TiU/((BmagU*BnT2BG)^2)          
	print,'Downstream ion beta',4.03e-11*nD*TiD/((BmagD*BnT2BG)^2)          
	
	print,'Upstream e- beta',4.03e-11*nU*TeU/((BmagU*BnT2BG)^2)          
	print,'Downstream e- beta',4.03e-11*nD*TeD/((BmagD*BnT2BG)^2)          
	
	print,'Upstream total beta',4.03e-11*nU*(TiU+TeU)/((BmagU*BnT2BG)^2) 
	print,'Downstream total beta',4.03e-11*nD*(TiD+TeD)/((BmagD*BnT2BG)^2) 
	;-------------------------------------------
	
	
	
	;relativistic or non-relativistic shock?
		;cond1
		val = 1.4d4*sqrt(n)/Bmag_nt     ;Eqn 1 in Treumann09 shock bible
		;if val >> MA then shock is non-relativistic
		;cond2
		;T << 0.511 MeV
	
	
	Mc_max = 2.76   ;eqn29 Treumann09 (largest critical mach number. See Fig2 in this paper
					;for Mc as a function of theta_bn and beta)
	

	;-------------
	;Thickness of quasi-perp laminar bow shock (Russell82, eqn1; Galeev76)
	;Value should typically be ~1 ion inertial length for this type of shock.
	;Note that this may not apply to supercritical shocks. 
	
	Mms = VprojU_mag/CmsU
	betaU_e = 4.03e-11*nU*TeU/((BmagU*BnT2BG)^2)    
	Lsh = ion_inertial*(me/mi)^0.25 * (Mms^2 - 1)/sqrt(betaU_e)
	;-------------


;----------------------------------------------------------------------------------
;Plasma drift values
;----------------------------------------------------------------------------------


	;--------------------------------------------
	;Drift period of particle around Earth (Kivelson eqn 10.6)
	;Dipole field, grad/curvature drift
	r = 2.6   ;Radius period in RE
	B = 1700.  ;Bo in nT
	keV = 100. ;electron energy
	mlt_extent = 1    ;hours
	
	Td = 60. * 56.*(r/5)^2 * (B/100.) * (1/keV) * (mlt_extent/24.)  ;minutes


	;-------------------
	;ExB drift    (km/s)
	Emag = 5.    ;(mV/m)
	Bmag = 150.    ;(nT)
	angle = 90.
	
	Ve = sin(angle*!dtor)*1000*Emag/Bmag

	;-------------------
	;Grad-B drift  (km/s)   [Wperp * (B x gradB)/(qB^3)]
	m = me       ;mass (kg)
	q = e        ;charge (C)
	Bmag = 10.   ;(nT)
	
	;now enter either perp velocity or perp energy of gyrating particle
		vperp = xxx  ;Perp velocity of particle (km/s)
		Wperp = 0.5*m*vperp^2  
	
	angle = 0.   ;angle b/t Bo and grad-B (deg)
	L = 10  ;gradient scale size (km)

	Vgb = Wperp*sin(angle*!dtor)/(q*Bmag*L)
	

	;-------------------
	;Curvature drift (km/s)  [2*Wpar * (rc x B)/(q*Rc*B^2)]
	m = me       ;mass (kg)
	q = e        ;charge (C)
	Bmag = 10.   ;(nT)
	angle = 10.  ;angle b/t Bo and grad-B (or b/t Bo and Rc)

	;now enter either par velocity or energy of gyrating particle
		vpar = xxxx  ;guiding-center velocity of gyrating particle (km/s)
		Wpar = 0.5*m*vpar^2

	;Method 1 (B*grad-B)
	L = xx       ;curvature scale size (km)
	Vc1 = -2*Wpar/(q*L*B)*cos(angle*!dtor)

	;Method 2 (Rc x B)
	Rc = xx      ;Radius of curvature (optional - for Vc2 calculation)
	Vc2 = 2*Wpar/(q*B*Rc)*sin(angle*!dtor)
	
	Vc = Vc1
 
	;-------------------
	;Polarization drift (km/s)
	m = me         ;mass (kg)
	q = e          ;charge (C)
	Bmag =         ;(nT)
	t_EB =         ;timescale for change of ExB drift
	t_gb_B =       ;timescale for change of grad-B drift
	t_cB =         ;timescale for change of curvature drift
	angle_EB =     ;angle b/t ExB drift and Bo
	angle_gb_B     ;angle b/t grad-B drift and Bo
	angle_cB =     ;angle b/t curvature drift and Bo
	
	
	Vp = -m/(q*Bmag^2)*((Ve*Bmag*sin(angle_EB)/t_EB) + (Vgb*Bmag*sin(angle_gb_B)/t_gb_B) + (Vc*Bmag*sin(angle_cB)/t_cB))


;-----------------------------------------------------------------------------------
;Instability calculations
;-----------------------------------------------------------------------------------

	wpe = 2*!pi*12300.
	wpi = 2*!pi*280.
	wce = 2*!pi*87.
	wci = 2*!pi*0.05
	wlh = 2*!pi*2.
	Debye = 0.015
	pe = 2.2
	pi = 267.
	mi = mp
	Z = 1
	Te = 8.4 ;eV
	Ti = 66.3 ;eV
	Cs_i = 9.79e5*sqrt(gamai*Z*Te/muu)/100./1000.        ;Ion sound vel (km/s)
	Cs = Cs_i

	k = 0.05*1000.  ;1/km

	;Buneman
	
		;-  important for wpe>>wce (unmagnetized e- and ions generally, or magnetized
		;       electrons (Hall current) and unmagnetized ions)
	
		wr = 0.5*(me/(2*mi))^0.333 * wpe
		km = wpe/Vd
	
	;(electron current driven) Ion acoustic
	
		;-  kpar mostly
		;-  Vd/Ve >= Ti/Te
	
		;for k*Debye << 1
			wr = k*Cs   
		
		;max growth occurs for:
			wr = wpi/sqrt(3)
			km = (1/sqrt(2))/Debye
			wi = (1/3.)*sqrt((!pi/6.)*(me/mi))*Vd*wpi/Cs
			coll_freq = 10d-2*(Te/Ti)*(Vd/Ve)*wpe   ;effective resistivity
	
	
	;Beam cyclotron (e- current driven ion cyclotron instability)
	
		;-  Threshold given by  Vd > max(Cs, wce/wpe * Ve)
		;-  Tends to stabilize by even small Bo inhomogeneities. 
		;-  Provides little resistivity unless Vd -> Ve. May see waves but unlikely to
		;		be important.
	
		;-  For cold ions (Vi << Vd) max growth occurs for:
			km = (1/sqrt(2))/Debye
			wr = km*(Vd + Cs)
			wi = sqrt((!pi/18.)*(me/mi))*Vd*wpe/Cs
			coll_freq = 10d-2*(Vd/Ve)^3*wce   ;effective resistivity
	
	;Modified two stream instability
	
		;-  Almost identical to LHDI except this one doesn't need grad-n
		;-  Can exist where Ti >= Te
		;-  Insensitive to Te/Ti
		;-  kpar/kperp = sqrt(me/mi)
		;-  should be important in space shocks
	
			wr = wlh
			k = 1/sqrt(pe*pi)   
			coll_freq = 0.1*wlh


	;ECDI
	
		;Forslund72 instability condition (Eqn15)
		;Vd/Ve ge val
		;theta = angle b/t k and Vd
		;n=cyclotron harmonic
		theta = 0.
		n=1
		val = n*wce/(cos(theta*!dtor)*sqrt(2)*wpe)

;-----------------------------------------------------------------------------------
;Dipole field stuff (see dipole.pro)
;-----------------------------------------------------------------------------------


;Basic equations:


ilat = acos(sqrt(1/L))/!dtor  ;invariant latitude
Lshell = rad/(cos(!dtor*mlat)^2)  ;L-shell in centered dipole
x_sm = rad*cos(mlat*!dtor)
z_sm = rad*sin(mlat*!dtor)


;Mlat calculated with SM coord (NOT GSM!!!)
dr2 = sqrt(SMx^2 + SMy^2)
dz2 = SMz[*,2]
mlat = atan(dz2,dr2)/!dtor


;MLT (I think this is technically calculated from GSM but GSE, GSM and SM give nearly 
;identical values...within 0.1 deg)
angle_tmp = atan(pos_gse_a.y[*,1],pos_gse_a.y[*,0])/!dtor
goo = where(angle_tmp lt 0.)
if goo[0] ne -1 then angle_tmp[goo] = 360. - abs(angle_tmp[goo])
angle_rad_a = angle_tmp * 12/180. + 12.
goo = where(angle_rad_a ge 24.)
if goo[0] ne -1 then angle_rad_a[goo] = angle_rad_a[goo] - 24




	;Create a structure with info along the entire length of a field line 
	;(From Kris's dipole.pro) returns a structure with the colatitude, field 
	;line length, and magnetic field magnitude along a field line defined by input L value

	L=2
	RE=6378.d ; earth radius [km]
	bsurf=30000.d ; equatorial surface field [nT]
	
	lambda=dindgen(9000)/100. ; colatitude [degrees]
	lambdarad=lambda*!dpi/180.d
	
	; arcsinh(x) = ln(x+sqrt(1+x^2))
	x=sqrt(3)*sin(lambdarad)
	s=RE*L/(2.*sqrt(3.))*(alog(x+sqrt(1.+x^2.)) + x*sqrt(x^2+1))
	r=RE*L*cos(lambdarad)^2
	goodindex=where(r gt RE)

	r=r[goodindex]
	s=s[goodindex]   ;arc length (last value is field line half-length)
	lambda=lambda[goodindex]
	lambdarad=lambdarad[goodindex]
	
	B=bsurf/(L^3.)*sqrt(1+x^2)/(cos(lambdarad)^6.)

	;another eq for Bo dipole
	Beq = 3.11d4  ;nT
	Bo = Beq*sqrt(1+3*cos(theta)^2)/r^3
	
	fce=.028d*b ; fce [kHz] for B [nT]
	
	struct={colat:lambda,r:r,s:s,B:B,fce:fce}

	plot,struct.R/6370.,struct.colat*!dtor,/polar
	plot,struct.R/6370.,-1*struct.colat*!dtor,/polar

;---------------------------



;SM longitude to MLT
longs = [25.3,32.5,46.7,129.0,278.8,300.7]
mlts = longs/15. - 12.
;now adjust so that the negative values are correct
foo = where(mlts lt 0.)
if foo[0] ne -1 then mlts[foo] = 24 + mlts[foo]


;---L,ilat,Lshell for specific spacecraft meridian
L = [3.18,2.49]
mlat = [30.74,29.67]
rad = [14970.,11952.]/6370.


;Mlat calculated with SM coord (NOT GSM!!!)
dr2 = sqrt(SMx^2 + SMy^2)
dz2 = SMz[*,2]
mlat = atan(dz2,dr2)/!dtor


ilat = acos(sqrt(1/L))/!dtor  ;invariant latitude
Lshell = rad/(cos(!dtor*mlat)^2)  ;L-shell in centered dipole
x_sm = rad*cos(mlat*!dtor)
z_sm = rad*sin(mlat*!dtor)

;----values for Earth/Sun meridian
L = [3.18,2.49]
mlat = [30.74,29.67]
colat = 90. - mlat
smLong = [74.39,93.94]
rad = [14970.,11952.]/6370.

ilat = acos(sqrt(1/L))/!dtor  ;invariant latitude
Lshell = rad/(cos(!dtor*mlat)^2)  ;L-shell in centered dipole
x_sm = rad*sin(colat*!dtor)*cos(smLong*!dtor)
y_sm = rad*sin(colat*!dtor)*sin(smLong*!dtor)
z_sm = rad*cos(colat*!dtor)


;---field strength at magnetic equator Bo=0.34 Gauss 
Bo=0.34 ;constant
Re=1.   ;Earth's radius 
req=3.1  ;L-shell at equator
Br = Bo*Re^3/req^3
;--equatorial radius at particular field strength
Req = (2000*30.4/Bnt2Bg)^0.3333
;determine field line of depending on value of fce
;Nose whistler (fastest vg freq) is at 0.25fce for homogeneous medium. 




;---Radius of curvature (RE) of dipole field line. Curvature drift becomes important when this is comparable
;---to cyclotron radius (Hamlin61)
Rc = L/3. * sin((90-mlat)*!dtor)*(1+3*cos((90-mlat)*!dtor)^2)^(3/2.)/(1+cos((90-mlat)*!dtor)^2)

;---Ratio of perp to parallel e- velocities in dipole field (Hamlin61)
pa = 85d ;equatorial particle pitch angle
mlat = 20d
mu = sin(pa*!dtor)
mu2 = mu^2  
B2Bo = (1+3*cos((90-mlat)*!dtor)^0.5)/(sin((90-mlat)*!dtor)^6)
vper2vpar = mu2 * B2Bo/(1-B2Bo*mu2)  ;NOTE THAT NEGATIVE VALUES MEAN THAT THE PARTICLE HAS REFLECTED

;---magnetic latitude of turning point (approximation from Hamlin61. Note that approximation (Fig2) is quite accurate)
mlat_turn = 90. - asin(mu^0.25)/!dtor


;-----------------------------------------------------------------------------------------------
;---bounce period for dipole field (http://farside.ph.utexas.edu/teaching/plasma/lectures/node22.html)
ang = 65.   ;PA at mag eq
L = 5.5
E = 20.      ;Energy in keV

;electrons
Tbe = 5.62d-2 * L * (1-0.43*sin(ang*!dtor))/sqrt(E/1000.)    ;seconds

;protons
Tbi = 2.41 * L * (1-0.43*sin(ang*!dtor))/sqrt(E/1000.)       ;seconds
;---------------------------------------------------------------------------------------------------









;r = [20679.449,17326.589,17326.589,17326.589,17326.589,17326.589,15182.513,15182.513,15182.513,9801.0049,8574.1678,8780.2046,9662.5222,15719.414,9958.0506,8011.7603,11796.83,18312.416]/6370.
;mlat = [18.580712,20.715154,20.715154,20.715154,20.715154,20.715154,22.004295,22.004295,22.004295,21.511855,17.034984,-9.2381702,-14.294564,21.698302,21.790091,0.40410432,-20.19017,-24.346044]


M = 30.4  ;magnetic moment in microTeslas*Re^3 = 30.4 microteslas for Re=1
R = (7878.)/6370.
mlat = 3.
theta = 90. - mlat  ;colat

;---dipole field in SM spherical coord (Kivelson and Russell pg 165)
Br = (1d3)*2*M*cos(!dtor*theta)/R^3   ;in nT
Btheta = (1d3)*M*sin(!dtor*theta)/R^3 
Btots = (1d3)*(M/R^3)*(1 + 3*cos(!dtor*theta)^2)^0.5  ;nT

;fce(mlat,R) for dipole field (Helliwell65 pg 181)
fc0=886000.  ;kHz  cycl freq at mag eq on Earth's surface
RE = 1.
fce = fc0*(1+3*sin(mlat*!dtor)^2)^0.5 /(R/RE)^3  ;Hz

;--equatorial cyclotron frequency value (worked out by Aaron from Btots equation)
fce_eq = fce*cos(2*mlat*!dtor)^3/sqrt(1+3*sin(mlat*!dtor)^2)
;--equatorial cyclotron frequency from footpoint (Helliwell65 pg 182)
fce_eq = fc0*cos(lat_foot*!dtor)^6
;;--from cold_dispersion_loop.pro (I DON'T THINK THIS IS CORRECT)
;Beq = Bo*(cos(latsc*!dtor)^6)/sqrt(1+3*sin(latsc*!dtor)^2)
;fce_eq = 28*Beq/1000.




;--------------------------------------------------
;LOSS CONE
;--------------------------------------------------

;sin(alpha)^2 = Bo/Bmax    ;Gurnett eqn 3.4.23
alpha = asin(sqrt(Bo/Bmax))/!dtor    

bdip = dipole(5.)

alpha_5 = asin(sqrt(0.004))/!dtor
alpha_3 = asin(sqrt(0.02))/!dtor

;Loss cone increases for decreasing L
;; RBSP_EFW> print,alpha_5
;;       3.62612
;; RBSP_EFW> print,alpha_3
;;       8.13010






;-----------------------
;EARTH'S DENSITY PROFILE
;-----------------------

	;Here's a very simple equatorial PS density model assuming that density 
	;falls off as 1/r^3 (Helliwell65 p 191,196)
	No = 10000. 
	N = No/(L/1.)^3  ;1/cm^3    ;THESE VALUES SEEM WAY TOO LOW
		
	;inner ps equatorial profile 2.24<L<Lppi from Carpenter and Anderson
	t1 = -0.3145*L + 3.9043
	t2 = -0.35*exp(-0.5/1.5)
	N = 10^(t1+t2)
	
	;e- density profile up to 90 km (Wait and Spies64)  -> also see Rodger and Nunn99
	;For higher altitude models use the IRI.
	n = 1.4265d13*exp(-0.15*h)*exp((bet-0.15)*(z-h_prime))
	



;----------------------
;COLD PLASMA PARAMETERS (see cold_dispersion.pro)
;----------------------

	P=1-(fpe^2/f^2)-(fpi^2/f^2)
	S=1-(fpe^2/(f^2-fce^2))-(fpi^2/(f^2-fci^2))
	
	theta_res = atan(sqrt(-1*P/S)) * 180/!pi

;--------------------
;DISPERSION RELATIONS
;--------------------

	;----------------------------
	;WHISTLERS (electrons only)
	;----------------------------

		f=1000.  ;Hz 
		fpe = 8980.*sqrt(30) 
		fce = 28.*100.
		thet = 0.
		c=3e5
		Ew = 14.37  ;mV/m
		;whistler dispersion for oblique, low freq whistlers
		k = sqrt(4*!pi^2*fpe^2*f/(c^2*(fce*cos(thet*!dtor)-f)))
		;Bw and Ew relation (nT and mV/m), c in km/sec for a FA, monochromatic plane wave whistler
		ckm = 3d5
		Bw = (1000.*Ew/f/ckm)*sqrt(f*fce^2/(fce-f))
		;in terms of vphase
		Bw = 1000.*Ew/vphase ;vphase in km/sec
		;in terms of vectors (not sure about units here)
		Bw_vec = c/(2*!pi*f)*(crossp(kvec,Ew_vec))
		
		
		Ew = Bw*f*ckm/1000./sqrt(f*fce^2/(fce-f))
		
		
		;------------------
		;QL WAVE DIFFUSION STUFF
		;the energy and pitch angle diffusion rates in QL diffusion scale as SB=Bw^2 
		;(magnetic power spectral density)
		;See Meredith04,07, Kennel and Engelmann66, Bell84
		;Sb=Bw^2 in nT^2
		;Se=Ew^2 in (mV/m)^2
		;For FA whistlers only!!!!!
		Sb = Se*1d6*fce^2/f/(ckm^2)/(fce-f) ;proportional to energy/pitch-angle diffusion rates
	
		;---Helliwell's angles for determining whistler trapping ("Theory of whistlers" pg. 49)
	
	;	frat = 0.04     ;f/fce ratio
	;	frat = 0.14
		frat = 10/347.
		theta_0 = 73.3.   ;initial wave normal angle (assume in center of crest/duct)
		theta_1 = acos((sqrt(frat/(1-frat)) + (frat/(1-frat) - 2*frat*sqrt(frat/(1-frat)))^0.5))/!dtor
		theta_2 = acos(2*frat)/!dtor ;Gendrin angle
		theta_3 = acos(frat/(1-frat))/!dtor
		theta_4 = acos(frat)/!dtor  ;resonance cone angle
	
		;crest/trough enhancement eq 3.39 (for crest n2>n1; for trough n1>n2)
		n1_2_n2=(cos(!dtor*theta_0)-frat)/(cos(!dtor*theta_0)^2 * (1-frat))
		;tough enhancement eq 3.44
		nmax2nmin=(cos(!dtor*theta_0))/(4*frat*(cos(!dtor*theta_0)-frat))
	
		;---find nose freq based on sc L-value (req)
		req = 1.4
		fnose = (0.34/req^3)*(28/(4.*bnt2bg))
	
	
		;whistler group velocity from vgroup.pro (km/s)
		;Equation for group velocity from Helliwell65 p30,31.
		;"Equation describes the propagation of a beat between two infinite plane waves of 
		;slightly different frequency whose wave-normal directions are the same"
		
		;Note that for f<0.5fce the group velocity is faster than phase velocity.
		Vg = 2*c*sqrt(f)*(fce*cos(theta_kb*!dtor)-f)^(3/2.)/(fpe*fce*cos(theta_kb*!dtor))
	
		;parallel component of Ew (to Bo) - Omura09 eqn42
		epsilon = sqrt(f*(fce - f)/fpe^2)	
		delta = 1/(1 + epsilon^2)             ;in terms of freq
		Ew_par = f*sin(theta_kb*!dtor)*Ew/(delta^2 * fce - f) ;Ew is the component in transverse plane (Ex-stix)
	
	
		;Whistler electric field polarizations (Stix coord)
		Ez = Ex * n^2*cos(th)*sin(th)/(n^2*sin(th)^2 - P)
		Ey = Ex * D/(n^2 - S)
	
		;another formulation - from eqn6 "Properties of obliquely propagating chorus", 2010
					;Olga P. Verkhoglyadova,1,2 Bruce T. Tsurutani,1 and Gurbax S. Lakhina3 
		
		num = 2*!pi*freq*sin(theta)*cos(theta)/(2*!pi*fc*cos(theta)-2*!pi*freq)
		den	= = 1 + 2*!pi*freq*sin(theta)^2/(2*!pi*fc*cos(theta)-2*!pi*freq)
		Ez2Ex = num/den
	
	
	
		;Whistler magnetic field polarizations (Stix coord) - see Verkhoglyadova and Tsurutani 2009
		a1_num = w^2*wpe^2/c^2/(w^2-wce^2) + k^2
		a1_den = w*wpe^2*wce/c^2/(w^2-wce^2)
		A1 = -1*a1_num/a1_den
		
		a2_num = wpe^2/c^2 + kperp^2
		a2_den = kperp*kpar
		A2 = a2_num/a2_den
		
		Bx = -kpar/kperp * Bz
		Bx = kpar*By/A1/(kpar - kperp/A2)


	;----------------------------
	;quasi-electrostatic whistlers (Bell and Ngo, 1990)
	;----------------------------


		;theta is the wavenormal angle of the incident EM whistler
	
	
		f = 21400. 
		fce = 731000.
		fpe = 1889000. 
		flhr = 11580.
		theta = indgen(90)
	
		nx = (fpe*sqrt(fce*cos(!dtor*theta)))/(sqrt(f*(f^2-flhr^2)))
		lambda = c/(nx*f)
	

	;----------------------------
	;EMIC (parallel only)
	;----------------------------

		c = 3e8   ;light speed
		n = 50000. ;cm^-3
		B = 10000. ;nT
		BG = B*BnT2BG  ;Gauss
		VA = 2.18e11*BG/sqrt(n)/100.  ;Alfven speed m/s
	
		f = 200d  ;Hz
		fpe = 8980*sqrt(n)
		fpi = fpe/43.
		fce = 28*B
		fci = fce/1836.
	
		c1 = 4*!pi*!pi/c^2
		k = sqrt(c1*(1-(fpe^2 + fpi^2)/((f+fce)*(f-fci))))  ;1/m
		lambda = 2*!pi/k   ;m


	;---------------------
	;EIC  
	;---------------------

	;	w^2 = wce^2 + k^2*vs^2
		
		c = 3e8   ;light speed
		n = 50000. ;cm^-3
		B = 10000. ;nT
		BG = B*BnT2BG  ;Gauss
	
		f = 200d  ;Hz
		fpe = 8980*sqrt(n)
		fpi = fpe/43.
		fce = 28*B
		fci = fce/1836.
		
		w = 2*!pi*f
		wci = 2*!pi*fci
		gama_e = 1  ;isothermal e-
		gama_i = 3  ;1D motion
	;	TeeV = 1.
	;	TieV = 0.1
		Te = 2000.  ;K
		Ti = 2000.   ;K
		M = 1.6d-27
		Z = 1
	
		vs = sqrt(gama_e*kb*Te + gama_i*Z*kb*Ti)/sqrt(M)   ;m/s
		
		k = sqrt((w^2 - wci^2)/(vs^2))     ;1/m
		lambda = 2*!pi/k                 ;m
		
		print,lambda
	


	;------------------------------
	;Alfven waves (ideal, parallel only)
	;--------------------------------
	;ka = sqrt(4*!pi*!pi*f^2/VA^2)
	
		ka = f*2*!pi*sqrt(1+(VA^2/c^2))/VA
		lambda_a = 2*!pi/ka
		theta = 90*!dtor
		L = 1.8


	;------------------------------
	;Ion acoustic waves
	;--------------------------------

		gamma_e = 1
		gamma_i = 3
		kB     = 1.380658d-23      ; -Boltzmann Constant (J/K)
		Te = 1000  ;K
		Ti = 800   ;K
		Z = 1      ;charge state
		M = 1.67d-27  ;ion mass
		;sound velocity
		Cs = sqrt((gamma_e*kb*Te + Z*gamma_i*kb*Ti)/M)
	
		w = 2*!pi*21400.	
		;IAWs are dispersionless so w = vs*k 
		k = w/Cs    ;1/m
		lambda = 2*!pi/k


	;--------------------------------
	;Bernstein waves (Chen p278)
	;--------------------------------

	kperp = 14.   ;1/km
	fce = 560 	  ;Hz
	fpe = 28000.	;Hz
	Te = 20       ;eV
	Ve = 4.19e7*sqrt(Te)/100./1000.                      ;e- thermal vel (km/s)
	pe = 0.5     ;km  (Larmour radius)

	wce = 2*!pi*fce
	wpe = 2*!pi*fpe

	kD2 = 2*wpe^2/(Ve^2)
	b = kperp^2 * Ve^2/(2.*wce^2)

	alpha = kperp^2/kD2   ;both alpha and kperp are functions of frequency


	pe*kperp = 7    ;From Chen Figure 7-34 this means that the harmonics should
					;be almost exact multiples of n*fce. Better yet, Fig 1 of
					;Crawford 65 shows the dispersion curves for wpe/wce=infinity, which
					;is a good approximation for the bow shock (assumes a Maxwellian
					;transverse e- velocity distribution)
					
					;Every harmonic must have a doubled value of k. Here's why:
					;The phase velocity of each harmonic will be the same. Because
					;the freq doubles so must k










;------------------
;Doppler shift plot
;------------------

	t0 = (2*!pi*f/(3e8))^2
	t1 = fpe^2/(f*(f+fce))
	t2 = fpi^2/(f*(f+fci))
	k = sqrt(abs(t0*(1-(t1+t2))))
	DS = k*cos(theta_k*!dtor)*Vsw*1000./2./!pi
	DSS = k*cos(thetas*!dtor)*Vsw*1000./2./!pi
	print,'Doppler shift = ' + strtrim(DS,2) + ' Hz'
	fwave = f - DS
	print,'k = ' + strtrim(k,2) + ' 1/m' + ' --> wave vector mag for whistler wave in sc frame'
	print,'fwave = ' + strtrim(fwave,2) + ' Hz --> freq in wave frame'
	plot,DSS,xtitle = 'Theta (degrees) -- theta res (cold plasma) = ' + strtrim(theta_res,2),ytitle='Doppler shift freq (Hz)'
	
	device,decomposed=0
	loadct,39
	thetas_10Hz = acos(10/k/indgen(1000)/1000.)*180/!pi
	thetas_30Hz = acos(30/k/indgen(1000)/1000.)*180/!pi
	thetas_50Hz = acos(50/k/indgen(1000)/1000.)*180/!pi
	thetas_70Hz = acos(70/k/indgen(1000)/1000.)*180/!pi
	thetas_90Hz = acos(90/k/indgen(1000)/1000.)*180/!pi
	thetas_110Hz = acos(110/k/indgen(1000)/1000.)*180/!pi
	title = 'Vsw and theta_k when Doppler shift = XX Hz (k= ' + strtrim(k,2) + ' 1/m)'
	plot,thetas_10Hz,xtitle = 'Solar Wind Velocity (km/sec)',ytitle='wave normal angle (deg)',title=title,color=255
	oplot,thetas_30Hz,color=225
	oplot,thetas_50Hz,color=195
	oplot,thetas_70Hz,color=165
	oplot,thetas_90Hz,color=135
	oplot,thetas_110Hz,color=105
	
	;--------Range of possible k-values-----------
	;This is equation 4 from Coroniti82. When the conditions w'<<wce*cos(theta_k) and k'=kc/wp<<1 are satisfied
	;then this equation can be used. The limiting values are: k1, the max red-shifted value and k2, the blue-shifted value
	;theta_k can be estimated with cold_dispersion_loop.pro which plots the wave normal angles for all reasonable
	;wave freqs in SW frame. How well this works depends on how small the range of possible theta_k values is 
	;because, among other things, I use this to calculate delta_kv. 
	
	delta_kv = 45.   ;angle b/t Vsw and k
	
	t1 = vsw*abs(cos(delta_kv*!dtor))/2./vae/abs(cos(theta_k*!dtor))
	t2 = vsw^2*cos(delta_kv*!dtor)^2/4./vae^2/cos(theta_k*!dtor)^2
	t3 = f/fce/abs(cos(theta_k*!dtor))
	k1p = t1 + sqrt(t2 + t3)
	k1 = k1p*fpe*2*!pi/c/1e-3   ;red-shifted
	k2p = -1*t1 + sqrt(t2 + t3)
	k2 = k2p*fpe*2*!pi/c/1e-3   ;blue-shifted

;-------------------------------------------------------------------------------
;RELATIVISTIC CYCLOTRON RESONANCE ENERGIES (NEED TO CHOOSE PARTICLE PITCH ANGLE)
;-------------------------------------------------------------------------------

	;freq in Hz
	;all distances in meters
	
	
	pa = 0.  ;test pitch angle in deg
	theta_k = 20.0  ;wave normal angle
	freq = 8000.
	fce = 68000.
	fpe = 17000.
	c=3d8
	kvec = sqrt(4*!pi^2*fpe^2*freq/(c^2*(fce*cos(theta_k*!dtor)-freq)))
	;kvec = 1./1000.   ;1/m
	
	kz = abs((kvec)*cos(theta_k*!dtor))
	
	;test velocities
	vz = indgen(10000)*c/9999. * cos(pa*!dtor)
	
	gama = 1/sqrt(1-(vz/c/cos(pa*!dtor))^2)  ;relativistic gamma factor
	
	f1cycl = vz
	f2cycl = (2*!pi/kz)*(1*fce/gama - freq)
	
	
	diff = abs(f1cycl-f2cycl)
	tmp = min(diff,val)
	
	vz_cycl = vz[val]   ;m/s
	vtots_cycl = vz_cycl/cos(pa*!dtor)  ;electron velocity
	vc_ratio_cycl = vtots_cycl/c
	if vc_ratio_cycl ge 0.1 then relativistic_cycl = 'yes' else relativistic_cycl = 'no'
	;------------------------------------------------
	;anomalous resonance energy
	f1anom = -1*vz
	f2anom = (2*!pi/kz)*(-1*fce/gama - freq)
	diff = abs(f1anom-f2anom)
	tmp = min(diff,val)
	
	vz_anom = vz[val]   ;m/s
	vtots_anom = vz_anom/cos(pa*!dtor)  ;ion velocity
	vc_ratio_anom = vtots_anom/c
	if vc_ratio_anom ge 0.1 then relativistic_anom = 'yes' else relativistic_anom = 'no'
	;-----------------------------------------------
	;landau resonance energy
	vz_landau = 2*!pi*freq/kz
	vtots_landau = vz_landau/cos(pa*!dtor)
	vc_ratio_landau = vtots_landau/c
	if vc_ratio_landau ge 0.1 then relativistic_landau = 'yes' else relativistic_landau = 'no'
	 
	
	;Relativistic energy in eV (e.g. p37 in "Modern Physics, 2nd edition")
	Etots_cycl = 0.511d6/sqrt(1-(vtots_cycl^2/c^2)) - 0.511d6
	Etots_anom = 0.511d6/sqrt(1-(vtots_anom^2/c^2)) - 0.511d6
	Etots_landau = 0.511d6/sqrt(1-(vtots_landau^2/c^2)) - 0.511d6
	
	Ez_cycl = 0.511d6/sqrt(1-(vz_cycl^2/c^2)) - 0.511d6
	Ez_anom = 0.511d6/sqrt(1-(vz_anom^2/c^2)) - 0.511d6
	Ez_landau = 0.511d6/sqrt(1-(vz_landau^2/c^2)) - 0.511d6




;----------------------------
;TRANSVERSE RESONANCE ENERGY
;----------------------------

	;This occurs for waves whose perp wavelength is less than or equal to the gyroradius. 
	;The condition to be satisfied is  w~k_perp*v_perp
	;Lower hybrid waves can resonate with high energy ions this way. (See Bell et al., 2004 - CLUSTER....)
	
        f = 5000.
        kperp = 2.*!pi/400.


	vperp = 2*!pi*f/kperp   ;perp velocity. k in 1/m, f in Hz
	
	;relativistic energy in eV
	E_trans_electron = me_ev/sqrt(1-(vperp^2/c^2)) - me_ev
	E_trans_proton = mp_ev/sqrt(1-(vperp^2/c^2)) - mp_ev

        print,E_trans_proton/1000.  ;keV

;----------------------------------------------------------------------------------------------------------------------------------------
;NONRELATIVISTIC CYCLOTRON RES ENERGIES (NOTE: FOR ANYTHING IN THE 
;VICINITY OF 100 KEV YOU DEFINITELY NEED TO USE RELATIVISTIC EQNS ABOVE
;----------------------------------------------------------------------------------------------------------------------------------------


	;1) IN TERMS OF Bo, k
	;kvec = 0.34  ;1/km
	;k_unitless = kvec*c*1e-3/2./!pi/fpe
	;Elandau = mepp*k_unitless^2/(1+k_unitless^2)^2/e1ev  ;eV  (Kennel66)
	;Vlandau = sqrt(2*qq*Elandau/me)/1000.                ;km/s
	;Ecycl = mepp/k_unitless^2/(cos(theta_k*!dtor)^2)/e1ev  ;eV (only for w<<wce) (KennelandPetschek66)
	;Vcycl = sqrt(2.*qq*Ecycl/me)/1000.                    ;km/sec
	
	
	
	
	;2) IN TERMS OF FREQ, k
	;NOTE THAT THE CYCL RES ENERGY EQUATION FROM KENNEL66 IS WRONG (K SHOULD BE OUTSIDE OF THE PARENTHESIS AND IN THE DENOMINATOR)
	kvec = kvec/1000.
	;The below equations come from solving for v_par in the Doppler shifted resonance condition. 
	;See Gurnett p 9-19
	Ecycl = 0.5*me*4*!pi^2 * (f - fce)^2/(kvec*cos(theta_k*!dtor))^2/e1ev  ;eV
	Eanom = 0.5*me*4*!pi^2 * (f + fce)^2/(kvec*cos(theta_k*!dtor))^2/e1ev  ;eV
	Elandau = 0.5*me*4*!pi^2*(f)^2/(kvec*cos(theta_k*!dtor))^2/e1ev        ;eV
	
	vcycl = 2*!pi*(f-fce)/((kvect)*cos(theta_k*!dtor))   ;m/s
	
	;Check to make sure that the velocities are not relativistic
	if vcycl/c ge 0.9 then Ecycl = !values.f_nan & Eanom = !values.f_nan



;##################################################
;plot theta_k as a function of Ez/Ey
;##################################################

plot_theta_k = 'no'
if plot_theta_k eq 'yes' then begin
	plot,[0,0],/nodata,xrange=[0,100],yrange=[0,100]
	epol_tmp = indgen(100)
	thetadeg_tmp = fltarr(100)
	for i=0,99 do begin 								$
		n2= epol_tmp[i]*vals.cp_params.D+vals.cp_params.S								& $
		n=sqrt(n2)										& $
		kvect=n*2*!pi*vals.freq/(3*100000.)							& $
		tan2=-1*vals.cp_params.P*(n2-vals.cp_params.R)*(n2-vals.cp_params.L)/((vals.cp_params.S*n2-vals.cp_params.R*vals.cp_params.L)*(n2-vals.cp_params.P))		& $
		tan=sqrt(tan2)									& $
		theta=atan(tan)									& $
		thetadeg_tmp[i]=(360./6.283)*theta				
;	endfor 
	
	oplot,epol_tmp,thetadeg_tmp
endif




;--------------------------
;POYNTING FLUX CALCULATIONS (MONOCHROMATIC, PLANE WAVE)
;--------------------------

	;--in terms of Bw
	c=c
	n=15. ;index of refraction
	muo=muo
	Bnt = 2  ;Bw in nanoteslas (perp component only I think)
	S_b = (1d-9)*(1d-9)*(c/n/muo)*Bnt^2  ;Poynting flux from magnetic field in Watts/m^2 (for a single freq)
	;--in terms of Ew (useful for electrostatic waves)
	epo=epo
	Ew = 100 ;mV/m
	S_e = (c*epo/2./1000./1000.)*Ew^2 ;Poynting flux from electric field in Watts/m^2 (for a single freq)



;---------------------------------
;Radiation calculations
;---------------------------------

	b=2.9e-3  ;Wien's displacement constant (m*K)
	;Tbb=5800.   ;Blackbody temp (K)
	Tbb=1e6
	
	lambda_max = 1e9*b/Tbb 		;wavelength of max radiation (nm)
	fmax = 1e9*c/lambda_max		;freq of max radiation (Hz)
	Emax=h*fmax					;energy (eV)
	
	
	lambda_max = 1e-9  ;x-rays
	;lambda_max = 1e-12 ;gamma rays
	;lambda_max = 1e-13 ;cosmic rays
	;lambda_max = 1e1   ;radio waves
	Tbb = b/lambda_max

;------------------------------
;E and B amplitudes
;------------------------------

	;Assuming a plane wave, k=kz, Bz=0 we get the simple relation:
	;Ey=w/k * Bx  in SI units. E is in V/m and B is in T. 
	
	
	freq=100.  ;Hz
	k=0.0125  ;1/km --> appropriate for DS=5 Hz and vsw = 400 km/s
	k=k/1000. ;1/m
	BnT = 0.1  ;nT
	
	Eamp = 1000.*2.*!pi*freq/k  * BnT/1e9   ; in mV/m
	
	;More generally the equation is E=(w*Bd)/(k*sin(d)), where d is the angle b/t k and E, and Bd is the 
	;projection of B in the kxE direction. k is no longer necessarily along Bo. 
	
	d = 5.
	
	Eamp2 = 1000.*2.*!pi*freq/k  * BnT/1e9/sin(!dtor*d)   ; in mV/m
	
	;Calculate phase velocity from this
	Vp = (Eamp2/1000.)/(BnT/1e9)/1000.   ;km/sec
	
	;Plot E vs d assuming minimal Doppler-shift
	;f~fobs
	alpha = 1000.*2.*!pi*freq * BnT/1e9   ;constant
	vsw = 400.*1000.  ;b/c k is in 1/m in the Eamp2 equation
	DS_max = 5.  ;max DS (Hz)
	k_max = 2*!pi*DS_max/vsw
	
	;The inequality is...
	;E(d) >= vsw*alpha/DS_max/sin(d) 
	
	d = indgen(90)
	
	Emin = vsw*alpha/DS_max/sin(!dtor*d)
	

;-----------------------
;ATMOSPHERIC STUFF.....
;-----------------------


	;------------------------
	;LIGHTNING CALCULATIONS
	;------------------------
	
		;Peak of lightning amplitude angle for radiation field for a given altitude (Rowland98)
		theta_max = 0.5*acos(b^2/(2-b^2))
		
		;Lightning current as a function of time (Cho and Rycroft98)
		xxx
		
		;Electric field from dipole radiator (Kelley90)
		E = THETA*p/(4*!pi*e_o*r^3) * (2*cos(theta) + sin(theta))  ;warning, don't use this without looking it up. It is a 
																	;function of both r and theta (colat).
																   ;THETA = joule extinction effects.
		
		;e- density perturbation due to EMP (Moore03) -> also see Inan91, Cheng and Cummer05
		dn = A1*exp(-(h-h1)/sigma1)^2 - A2*exp(-(h-h2)/sigma2)^2
		
		;Nondeviative absorption coefficient (dB/km) (Davies90)
		K1 = 4.6d4 * (n_z*v_z)/((w+w_L)^2 + v_z^2)
		K2 = 4.6d4 * (n_z*v_z)/((w-w_L)^2 + v_z^2)   ;z = altitude (km)
													 ;n_z = density as function of alt (cm-3)
													 ;v_z = e-/neutral collision freq as function of alt (s-1)

	;-----------------------
	;COLLISION CALCULATIONS
	;-----------------------
	
		;e-/e- collision freq (Schunk and Nagy00)
		vee = 1.5d-4 * E^1.5 * N_e    ;E=energy in eV
		
		;e-/neutral collision freq as a function of altitude ("standard approx" -> Morfitt and Shellman76)
									;Also see: Rishbeth and Garriott69, pg. 133
									;Schunk and Nagy00, pg 131
		v_z = 1.82d11*exp(-0.15*z)
		
		;e-/neutral collision freq (elastic)
		ven = sigma_en * V * N     ;sigma_en = cross section (Banks and Kocharts73)
								   ;N = neutral density
								   
		;e-/neutral collision freq (inelastic)
		ven = sigma_kj * V * N_k   ;sigma_kj = cross section of jth state for kth species
								   ;N_k - density of kth species
		
		;perturbation of e-/neutral collision freq from e- heating (Rodger98)
		ven_p = ven*(Te/To)^0.5   ;See Inan91 for Te profile (Fig2b)
	

	;--------------------------------
	;MINIMUM PENETRATION FREQS OF IONOSPHERE
	;--------------------------------
	
		;First a definition: The penetration freq is that which makes n=0 at the max plasma density
		;of the layer (Budden '61 section 13.2)
		
		;plasma and cyclotron freqs at max plasma freq
		fceM = xx
		fpeM = xx
		
		;O-mode
		fminO = fpeM    
		;X-mode
		fminX1 = 0.5*(sqrt(fceM^2 + 4*fpeM^2) + fceM)  ;first root (normally observed when f<fceM)
		fminX2 = 0.5*(sqrt(fceM^2 + 4*fpeM^2) - fceM)  ;second root


