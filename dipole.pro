;+
;*****************************************************************************************
;
;  FUNCTION : dipole_field
;  PURPOSE  : Returns a structure with various centered-dipole field values
;
;  CALLED BY:
;
;
;  CALLS:
;
;
;  REQUIRES:
;
;
;  INPUT:       L -> L-shell
;
;
;  EXAMPLES:    To plot use:
;					plot,struct.R/6370.,struct.lat*!dtor,/polar
;					plot,struct.R/6370.,-1*struct.lat*!dtor,/polar
;
;
;
; RETURNS:      lat (deg)
;				R -> radius (km)
;				S -> arc length
;				B -> Bfield strength (nT)
;				fce -> (kHz)
;				Br -> radial B-component  (dipole coord)
;				Btheta -> latitude component (dipole coord)
;
;
;  KEYWORDS:
;
;
;   CHANGED:  1)  NA [MM/DD/YYYY   v1.0.0]
;
;   NOTES:      See Kivelson and Russell pg 165 for definition
;
;
;   CREATED:  07/26/2010
;   CREATED BY:  Aaron W. Breneman
;    LAST MODIFIED:  MM/DD/YYYY   v1.0.0
;    MODIFIED BY: Aaron W. Breneman
;
;*****************************************************************************************
;-

function dipole,L

; returns a structure with the latitude, field line length, and
; magnetic field magnitude along a field line defined by input L value


	RE=6370.d ; earth radius [km]
	bsurf=30000.d ; equatorial surface field [nT]

	theta=dindgen(9000)/100. ; magnetic latitude [degrees]
	thetarad=theta*!dpi/180.d

	; arcsinh(x) = ln(x+sqrt(1+x^2))
	x=sqrt(3)*sin(thetarad)
	s=RE*L/(2.*sqrt(3.))*(alog(x+sqrt(1.+x^2.)) + x*sqrt(x^2+1))
	r=RE*L*cos(thetarad)^2
	goodindex=where(r ge RE)

	if goodindex[0] ne -1 then begin
		r=r[goodindex]
		s=s[goodindex]
		theta=theta[goodindex]
		thetarad=thetarad[goodindex]
	endif else begin
		r = !values.f_nan
		s = !values.f_nan
		theta = !values.f_nan
		thetarad = !values.f_nan
	endelse

	B=bsurf/(L^3.)*sqrt(1+x^2)/(cos(thetarad)^6.)
	fce=.028d*b ; fce [kHz] for B [nT]


	M = 30.4  ;magnetic moment in microTeslas*Re^3 = 30.4 microteslas for Re=1
	R2 = (r)/6370.
;	theta = 90. - theta

	;---dipole field in SM spherical coord (Kivelson and Russell pg 165)
	Br = (1d3)*2*M*cos(!dtor*theta)/R2^3   ;in nT
	Btheta = (1d3)*M*sin(!dtor*theta)/R2^3

	notes = ['latitude in deg','r=radius in km','s=arc length in km','B in nT','fce in kHz',$
			'Br=radial Bfield component','Btheta=theta Bfield component']
	struct={lat:theta,r:r,s:s,B:B,fce:fce,Br:Br,Btheta:Btheta,notes:notes}
	
	return,struct
end
