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

;Returns tplot variables of useful quantities. 

pro plot_omni_quantities,noplot=noplot,t0avg=t0tmp,t1avg=t1tmp

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



    ;Various lines for overplotting on "IMF orientation", described in Petrinic13
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

    ;Various lines for overplotting on clock, cone angles
    store_data,'90line',ttmp,replicate(90.,n_elements(bx))
    store_data,'0line',ttmp,replicate(0.,n_elements(bx))
    store_data,'m90line',ttmp,replicate(-90.,n_elements(bx))
    store_data,'180line',ttmp,replicate(180.,n_elements(bx))
    store_data,'clockangle_comb',data=['clockangle','90line','0line','m90line','180line']
    store_data,'coneangle_comb',data=['coneangle','90line','0line','m90line','180line']
    options,['90line','m90line','0line','180line'],'color',250


    store_data,'IMF_orientation_comb',data=['IMF_orientation','nearperplineL','nearperplineH','nearperplinemL','nearperplinemH',$
    'nearPSline1','nearPSline2','nearPSline3','nearPSline4','betweenline1','betweenline2','betweenlinem1','betweenlinem2']


    if ~keyword_set(noplot) then begin
        tplot,['IMF_orientation_comb','Bz_rat'] & stop
        tplot,['IMF_orientation_comb','coneangle_comb'] & stop
        tplot,['clockangle_comb','coneangle_comb','kyoto_ae','kyoto_dst','OMNI_HRO_1min_flow_speed'] & stop
        tplot,['omni_press_dyn_smoothed','omni_press_dyn_smoothed_detrend'] & stop
    endif

end 