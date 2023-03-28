;Calculate the plasma frame frequency via Doppler-shift equation
;for whistler mode waves. Note that single-point measurements require the 
;assumption of a dispersion relation
;
;For example, see the Appendix of 
    ;Wilson III, L. B., et al. (2013), Electromagnetic waves and electron anisotropies 
    ;downstream of supercritical interplanetary shocks, 
    ;J. Geophys. Res. Space Physics, 118, 5–16, doi:10.1029/2012JA018167.
;and
    ;Agapitov et al, 2020: Sunward-propagating Whistler Waves Collocated 
    ;with Localized Magnetic Field Holes in the Solar Wind: Parker Solar Probe Observations at 35.7 Re Radi
;and
    ;Kretzschmar+21; Whistler waves observed by Solar Orbiter/RPW between 0.5 AU and 1 AU
;and 
    ;Columban+23: submitted to JGR --> they develop a methodology for whistlers after Encounter 1 where the Bwu SCM channel breaks. 


;To use, first calculate the wave normal vector using 3x Bw and min variance analysis 
;(a method using only 2x Bw from PSP encounters after Encounter 1, where the Bw_u SCM channel malfunctions, are derived in Colomban+23)
;Then input:
;fsc --> observed (spacecraft frame) frequency (Hz)
;tkv --> angle b/t the wave normal vector and plasma flow velocity (deg)
;tkb --> angle b/t the wave normal vector and Bo (deg)
    ;***NOTE: It is possible to use the SC frame tkv and tkb because the Vsw << c Lorentz transformations for Bw are
    ;Bw_par = Bw_par'
    ;Bw_perp = Bw_perp' 
;dens --> density in 1/cm-3
;Bo --> DC magnetic field in nT 
;Vsw --> solar wind velocity in km/s 

;The result will be one of two returned roots (e.g. sunward and anti-sunward flow). To determine which is the correct 
;plasma frame frequency you have to know the sign of k-hat. One way to do this is to 
;calculate the Poynting flux. The correct k-hat will have a component in this direction. 



;Solve the equation a*x^3 + b*x^2 + c*x + d
function newtfunc, v  
    common coeff, VN, ctkb, wscN
    vr = VN*v^3 + (ctkb - wscN)*v^2 + VN*v - wscN
    return, vr
end



function doppler_shift_calculate_whistlermode,fsc,tkv,tkb,dens,Bo,Vsw

    common coeff, VN, ctkb, wscN

    ;fsc = 50.  ;Hz
    ;tkv = 0.    ;deg
    ;tkb = 45.   ;deg
    ;dens = 300. ;cm-3
    ;Bo = 80.   ;nT
    ;Vsw = 500. ;km/s


    ;Some constants
    c = 3d5 ;km/s 
    BnT2BG = 1e-5               ;nT to Gauss (multiply by this number)
    me     = 9.1093897d-31     ; -Electron mass (kg)
    muo    = 4d0*!DPI*1d-7     ; -Permeability of free space (N/A^2 or H/m)


    ;Calculate some basic quantities
    ctkb = cos(tkb*!dtor)
    ctkv = cos(tkv*!dtor)
    wsc = 2.*!pi*fsc
    wpe = 2.*!pi*8980.*sqrt(dens)
    wce = 2*!pi*28.*Bo
    VA = 2.18e11*Bo*BnT2BG/sqrt(dens)/100./1000.             ;H+ Alfven vel (km/s)
    VAe = VA * sqrt(1836.)                               ;electron Alfven vel (km/sec)



    ;normalized quantities
    wscN = wsc/wce
    VN = Vsw*ctkv/VAe


    ;---------------------------------------------------
    ;This is not the main calculation, but just for kicks we'll also run the calculation
    ;in the low frequency approximation
    ;kc/wpe << 1 and w << wce (Coroniti82)

    T1n = Vsw*abs(ctkv)
    T1d = 2*VAe*abs(ctkb)
    T1 = T1n/T1d

    T2n = Vsw*ctkv 
    T2d = 2*VAe*ctkb
    T2 = (T2n/T2d)^2 

    T3n = wsc 
    T3d = wce*abs(ctkb)
    T3 = T3n/T3d

    ;Positive/negative solutions
    kNp = T1 + sqrt(T2 + T3)
    kNm = -1*T1 + sqrt(T2 + T3)

    kmagp = kNp * wpe/c  ;1/km
    kmagm = kNm * wpe/c  ;1/km


    ;Now find frequency in plasma frame 
    wplasmap = wsc - kmagp * Vsw * ctkv
    wplasmam = wsc - kmagm * Vsw * ctkv


    vpp = wplasmap/kmagp
    vpm = wplasmam/kmagm

    print,'Coroniti + fplasma (Hz) = ',wplasmap/2./!pi
    print,'Coroniti - fplasma (Hz) = ',wplasmam/2./!pi


    ;-------------------------------------------------------
    ;Full equation with no approximations

    kguess = [-0.3,0.3]  ;guess at wave normal number  (1/km)
                        ;Allow for both +/- values
    kN = NEWTON(kguess,'newtfunc')


    k = kN*wpe/c
    DSval = k * Vsw * ctkv
    wplasma = wsc - DSval

    vp = wplasma/k

    print,'fplasma (Hz, no approximations) = ',wplasma/2./!pi

    vals = {fplasma_Hz:wplasma/2./!pi, vp_km_s:vp, kmag_1_km:k, DSamount_deltaHz:-1*DSval/2./!pi}

    return,vals


end