;Load and plot OMNI quantities, including:
;   clock and cone angles, as well as the "IMF orientation". The third 
;   is from Petrinic et al., 2013, who shows that coherence tends to exist under
;   specific IMF conditions. 

;   Dynamic pressure
;   Dynamic pressure assuming constant density
;   Dynamic pressure assuming constant Vsw
;   Dynamic pressure: comparison of the three above quantities to show whether density or Vsw
;       is driving the dynamic pressure
;   Density 
;   Velocity
;   etc...

;Calculate OMNI dynamic pressure as n*v^2
;From OMNIWeb:
;Flow pressure = (2*10**-6)*Np*Vp**2 nPa (Np in cm**-3,
;Vp in km/s, subscript "p" for "proton")
;(NOTE THAT THIS CAN DIFFER STRONGLY FROM OMNI DATA, PRESUMABLY DUE TO SW EVOLUTION)

;KEYWORDS:
;   t0avg, t1avg --> start and stop times (strings like 2018-01-01/01:00:00) for averaging
;       solar wind density and velocity for the dynamic pressure comparisons (see above).
;   smoothtime --> if set, smooth the clock, cone, IMF orientation, Bz_rat values over this amount of time (hrs)
;
;
;Returns tplot variables of useful quantities. 

pro plot_omni_quantities,noplot=noplot,t0avg=t0tmp,t1avg=t1tmp,smoothtime=smoothtime

    rbsp_efw_init
    omni_hro_load


    ;If t0avg and t1avg aren't set, then set them to the first and last time of currently
    ;selected timerange 
    if not keyword_set(t0tmp) or keyword_set(t1tmp) then begin 
        trtmp = timerange()
        t0tmp = trtmp[0] & t1tmp = trtmp[1]
    endif
    t0tmp = time_double(t0tmp)
    t1tmp = time_double(t1tmp)

    if ~keyword_set(smoothtime) then smoothstr = 'no smoothing' else smoothstr = 'Smoothing='+string(smoothtime*60.,format='(f8.2)')+ ' min'

    ;Create |B| variable
    get_data,'OMNI_HRO_1min_BX_GSE',ttmp,bx
    get_data,'OMNI_HRO_1min_BY_GSE',ttmp,by
    get_data,'OMNI_HRO_1min_BZ_GSE',ttmp,bz
    bmag = sqrt(bx^2 + by^2 + bz^2)
    store_data,'OMNI_HRO_1min_Bmag',ttmp,bmag


    ;Average values of solar wind velocity (x-dir) and density for pressure comparison.
    vtmp = tsample('OMNI_HRO_1min_Vx',[t0tmp,t1tmp],times=ttt)
    ntmp = tsample('OMNI_HRO_1min_proton_density',[t0tmp,t1tmp],times=ttt)

    get_data,'OMNI_HRO_1min_Vx',data=vv
    get_data,'OMNI_HRO_1min_proton_density',data=dd


    vsw = vv.y ;change velocity to m/s
    dens = dd.y ;change number density to 1/m^3


    ;Pressure in nPa (rho*v^2)
    press_proxy = 2d-6 * dens * vsw^2
    store_data,'omni_press_dyn',data={x:vv.x,y:press_proxy}
    ;calculate pressure using averaged Vsw value
    vsw_mean = mean(vtmp,/nan)
    press_proxy = 2d-6 * dens * vsw_mean^2
    store_data,'omni_press_dyn_constant_vsw',data={x:vv.x,y:press_proxy}
    ;calculate pressure using averaged density value
    dens_mean = mean(ntmp,/nan)
    press_proxy = 2d-6 * dens_mean * vsw^2
    store_data,'omni_press_dyn_constant_dens',data={x:vv.x,y:press_proxy}

    store_data,'omni_pressure_dyn_compare',data=['omni_press_dyn','omni_press_dyn_constant_dens','omni_press_dyn_constant_vsw']
    ylim,'omni_pressure_dyn_compare',0,0
    options,'omni_pressure_dyn_compare','colors',[0,50,250]



    ;Detrend
    rbsp_detrend,'omni_press_dyn',60.*2.
    rbsp_detrend,'omni_press_dyn_smoothed',80.*60.
    get_data,'omni_press_dyn_smoothed_detrend',ttmp,dtmp
    store_data,'omni_press_dyn_smoothed_detrend',ttmp,dtmp






    ;-----------------------------------------------------------------

    ;^^Check the SW clock angle (GSE coord)
    ;Clockangle: zero deg is along zGSE, 90 deg is along yGSE 
    ;Coneangle: zero deg is along xGSE, 90 along r=sqrt(yGSE^2+zGSE^2)

    store_data,'clockangle',ttmp,atan(by,bz)/!dtor
    store_data,'coneangle',ttmp,acos(bx/bmag)/!dtor
    store_data,'IMF_orientation',ttmp,atan(by,bx)/!dtor
    store_data,'Bz_rat',ttmp,abs(bz)/(sqrt(bx^2 + by^2)/2.)


    if keyword_set(smoothtime) then begin 
        rbsp_detrend,'clockangle',60.*60.*smoothtime & copy_data,'clockangle_smoothed','clockangle'
        rbsp_detrend,'coneangle',60.*60.*smoothtime & copy_data,'coneangle_smoothed','coneangle'
        rbsp_detrend,'IMF_orientation',60.*60.*smoothtime & copy_data,'IMF_orientation_smoothed','IMF_orientation'
        rbsp_detrend,'Bz_rat',60.*60.*smoothtime & copy_data,'Bz_rat_smoothed','Bz_rat'

    endif

    ;Various lines for overplotting on "IMF orientation", described in Petrinic13
    store_data,'unityline',ttmp,replicate(1.,n_elements(bx))
    store_data,'nearperplineL',ttmp,replicate(70.,n_elements(bx))
    store_data,'nearperplineH',ttmp,replicate(110.,n_elements(bx))
    store_data,'nearperplinemL',ttmp,replicate(-70.,n_elements(bx))
    store_data,'nearperplinemH',ttmp,replicate(-110.,n_elements(bx))
    options,['nearperpline*'],'colors',250
    store_data,'nearPSline1',ttmp,replicate(125.,n_elements(bx))
    store_data,'nearPSline2',ttmp,replicate(145.,n_elements(bx))
    store_data,'nearPSline3',ttmp,replicate(-35.,n_elements(bx))
    store_data,'nearPSline4',ttmp,replicate(-55.,n_elements(bx))

    store_data,'betweenlinem1',ttmp,replicate(-10.,n_elements(bx))
    store_data,'betweenlinem2',ttmp,replicate(-30.,n_elements(bx))
    store_data,'betweenline1',ttmp,replicate(150.,n_elements(bx))
    store_data,'betweenline2',ttmp,replicate(170.,n_elements(bx))
    options,'betweenline*','colors',50

    store_data,'conepar1',ttmp,replicate(0.,n_elements(bx))    ;parallel 1
    store_data,'conepar2',ttmp,replicate(24.,n_elements(bx))   ;parallel 2
    store_data,'coneOPS1',ttmp,replicate(25.,n_elements(bx))  ;ortho Parker spiral1
    store_data,'coneOPS2',ttmp,replicate(65.,n_elements(bx))  ;ortho Parker spiral1
    store_data,'coneperp1',ttmp,replicate(66.,n_elements(bx))   ;Perp line 1
    store_data,'coneperp2',ttmp,replicate(114.,n_elements(bx))  ;Perp line 2
    store_data,'conePS1',ttmp,replicate(115.,n_elements(bx))    ;Parker spiral 1
    store_data,'conePS2',ttmp,replicate(155.,n_elements(bx))    ;Parker spiral 2
    store_data,'conepar3',ttmp,replicate(156.,n_elements(bx))    ;parallel 3
    store_data,'conepar4',ttmp,replicate(180.,n_elements(bx))   ;parallel 4




    ;Various lines for overplotting on clock, cone angles

    store_data,'PSline',ttmp,replicate(135.,n_elements(bx))  ;Parker spiral
    store_data,'OrthoPSline',ttmp,replicate(45.,n_elements(bx)) ;Ortho Parker spiral
    store_data,'90line',ttmp,replicate(90.,n_elements(bx))
    store_data,'0line',ttmp,replicate(0.,n_elements(bx))
    store_data,'m90line',ttmp,replicate(-90.,n_elements(bx))
    store_data,'180line',ttmp,replicate(180.,n_elements(bx))
    store_data,'clockangle_comb',data=['clockangle','90line','0line','m90line','180line']
;    store_data,'coneangle_comb',data=['coneangle','90line','0line','m90line','180line','PSline','OrthoPSline']
    store_data,'coneangle_comb',data=['coneangle','conepar1','conepar2','coneOPS1','coneOPS2','coneperp1','coneperp2','conePS1','conePS2','conepar3','conepar4']
    options,['conepar1','conepar2','conepar3','conepar4'],'colors',250
    options,['coneOPS1','coneOPS2'],'colors',200
    options,['coneperp1','coneperp2'],'colors',50
    options,['conePS1','conePS2'],'colors',0

    options,['90line','m90line','0line','180line'],'color',250
    options,'OrthoPSline','color',250 & options,'PSline','color',85
    options,'OrthoPSline','linestyle',3 & options,'PSline','linestyle',4
    options,'clockangle_comb','ytitle','Clockangle!Catan(by/bz)!C0 faces N!C180 faces S!C'+smoothstr
    options,'coneangle_comb','ytitle','Coneangle!Cacos(bx/|B|)!C0 faces sun!C 180 faces Earth!C135=ParkerSpiral!C45=OrthoParkerSpiral!C'+smoothstr


    store_data,'Bz_rat_comb',data=['Bz_rat','unityline']
    options,'Bz_rat_comb','ytitle','|Bz|/sqrt(bx2+by2)/2!C'+smoothstr
    ylim,'Bz_rat_comb',0.001,10,1


    store_data,'IMF_orientation_comb',data=['IMF_orientation','nearperplineL','nearperplineH','nearperplinemL','nearperplinemH',$
    'nearPSline1','nearPSline2','nearPSline3','nearPSline4','betweenline1','betweenline2','betweenlinem1','betweenlinem2']
    options,'IMF_orientation_comb','ytitle','IMForientation!Catan(by/bx)!CBlack=nearParkerSpiral!CRed=nearPerp!CBlue=nearParallel!C'+smoothstr


    if ~keyword_set(noplot) then begin
        tplot,['IMF_orientation_comb','Bz_rat_comb','clockangle_comb','coneangle_comb','omni_press_dyn_smoothed'];,'omni_press_dyn_smoothed_detrend']
    endif

end 