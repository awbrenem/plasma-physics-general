;General routine to calculate ExB velocity for THEMIS data.
;Best to grab data from CDAWeb. The online THEMIS files from the Berkeley site
;can't be read by cdf2tplot.

;See specific crib sheets for various sc, like exb_velocity_rbsp_crib.pro
;to call this routine.


;pro exb_velocity_themis_crib

    rbsp_efw_init

    date = '2014-01-11'
    sc = 'a'

    smootime = 8.  ;smooth over this many minutes

    timespan,date

    path = '/Users/aaronbreneman/Desktop/Research/RBSP_hiss_precip2_coherence_survey/Analysis_major_events_campaign2/Jan11/'

    ;load ephemeris data
    cdf2tplot,path+'tha_l1s_state_20140111000000_20140113000000.cdf'
    ;load spinfit, calibrated E, B data
    cdf2tplot,path+'tha_l2s_fit_20140111000001_20140112235959.cdf'

    tplot,['tha_efs_dot0_gse','tha_fgs_gse','tha_pos_gse']


    intedt = total(Ey,/nan,/cumulative)*find_datarate(etimes)
    store_data,'Time-Int-Ey!Cgse!CmV/m-s',data={x:etimes,y:intedt},dlim={constant:[0]}



;interp B to E times
    get_data,'tha_efs_dot0_gse',etimes,d
    tinterpol_mxn,'tha_fgs_gse',etimes
    get_data,'tha_fgs_gse_interp',data=B_gse

    B_mag = SQRT(B_gse.y[*,0]^2+B_gse.y[*,1]^2+B_gse.y[*,2]^2)


    rate = find_datarate(etimes[0:1000])
    smoo = round((smootime*60.)/rate)
    smoothed_int = string(format='(i0.1)',round(smoo*rate/60.))
    print,'Smoothing (averaging) the data '+smoothed_int+ ' min'

    rbsp_detrend,'tha_fgs_gse_interp',60.*smootime
    rbsp_detrend,'tha_efs_dot0_gse',60.*smootime


    get_data,'tha_fgs_gse_interp_smoothed',data=B_bkgrd
    get_data,'tha_efs_dot0_gse_smoothed',data=E_bkgrd

    Bx_bkgrd = B_bkgrd.y[*,0]
    By_bkgrd = B_bkgrd.y[*,1]
    Bz_bkgrd = B_bkgrd.y[*,2]
    Ex_smooth = E_bkgrd.y[*,0]
    Ey_smooth = E_bkgrd.y[*,1]
    Ez_smooth = E_bkgrd.y[*,2]

    B_field_mag = sqrt(Bx_bkgrd^2 + By_bkgrd^2 + Bz_bkgrd^2)


    vx = 1000.*(Ey_smooth*Bz_bkgrd - Ez_smooth*By_bkgrd)/B_field_mag^2
    vy = 1000.*(Ez_smooth*Bx_bkgrd - Ex_smooth*Bz_bkgrd)/B_field_mag^2
    vz = 1000.*(Ex_smooth*By_bkgrd - Ey_smooth*Bx_bkgrd)/B_field_mag^2

    edb_stat='E-dot-B-where-OK'

    store_data,'Vx!CExB-drift!Ckm/s',data={x:etimes,y:vx}
    store_data,'Vy!CExB-drift!C'+edb_stat+'!Ckm/s',data={x:etimes,y:vy}
    store_data,'Vz!CExB-drift!C'+edb_stat+'!Ckm/s',data={x:etimes,y:vz}

    options,['E-dot-B!Cflag','Vx!CExB-drift!Ckm/s','Vy!CExB-drift!CE-dot-B-where-OK!Ckm/s','Vz!CExB-drift!CE-dot-B-where-OK!Ckm/s'],thick=2,constant=0
    tplot,['E-dot-B!Cflag',$
    'Vx!CExB-drift!Ckm/s',$
    'Vy!CExB-drift!CE-dot-B-where-OK!Ckm/s',$
    'Vz!CExB-drift!CE-dot-B-where-OK!Ckm/s']




    ;put E field and ExB velocities in FAC coordinates
    nb= n_elements(Bx_bkgrd)
    bg_field = fltarr(nb,3)
    bg_field[*,0]  = Bx_bkgrd/B_field_mag & bg_field[*,1]  = By_bkgrd/B_field_mag & bg_field[*,2]  = Bz_bkgrd/B_field_mag

    nef = n_elements(Ex_smooth)
    dE_field = fltarr(nef,3)
    dE_field[*,0] = Ex_smooth & dE_field[*,1] = Ey_smooth & dE_field[*,2] = Ez_smooth

    nbf = n_elements(vx)
    exb_flow = fltarr(nbf,3)
    exb_flow[*,0] = vx & exb_flow[*,1] = vy & exb_flow[*,2] = vz



    get_data,'tha_pos_gse',data=gse_pos
    mptimes = gse_pos.x
    xgse = gse_pos.y[*,0] & ygse = gse_pos.y[*,1] & zgse = gse_pos.y[*,2]

    radial_pos = SQRT(xgse^2+ygse^2+zgse^2)

    xgse = interp(xgse,mptimes,etimes,/no_extrap)
    ygse = interp(ygse,mptimes,etimes,/no_extrap)
    zgse = interp(zgse,mptimes,etimes,/no_extrap)

    radial_pos = interp(radial_pos,mptimes,etimes,/no_extrap)

    r_dir_vec = fltarr(nbf,3) ;the vectors along the spin axis
    r_dir_vec[*,0] = xgse/radial_pos   ;REPLACE WITH RADIAL VECTOR gse
    r_dir_vec[*,1] = ygse/radial_pos
    r_dir_vec[*,2] = zgse/radial_pos


    ;define orthogonal perpendicular unit vectors
    perp1_dir = fltarr(nef,3)
    for xx=0L,nef-1 do perp1_dir[xx,*] = crossp(bg_field[xx,*],r_dir_vec[xx,*])  ;azimuthal, east
    perp2_dir = fltarr(nef,3)
    for xx=0L,nef-1 do perp2_dir[xx,*] = crossp(perp1_dir[xx,*],bg_field[xx,*]) ;radial, outward


    ;need to normalize perp 1 and perp2 direction
    bdotr = fltarr(nef)
    for xx=0L,nef-1 do bdotr[xx] = bg_field[xx,0]*r_dir_vec[xx,0]+bg_field[xx,1]*r_dir_vec[xx,1]+bg_field[xx,2]*r_dir_vec[xx,2]

    one_array = fltarr(nef)
    one_array[*] = 1.0
    perp_norm_fac1 = SQRT(one_array - (bdotr*bdotr))
    perp_norm_fac = fltarr(nef,3)
    perp_norm_fac[*,0] = perp_norm_fac1 & perp_norm_fac[*,1] = perp_norm_fac1 & perp_norm_fac[*,2] = perp_norm_fac1

    perp1_dir = perp1_dir/(perp_norm_fac) & perp2_dir = perp2_dir/(perp_norm_fac)


    store_data,'perp1_dir!Cazimuthal!Ceastward!Cunit-vec',data={x:etimes,y:perp1_dir},dlim={colors:[2,4,6],labels:['x','y','z']}
    store_data,'perp2_dir!Cradial!Coutward!Cunit-vec',data={x:etimes,y:perp2_dir},dlim={colors:[2,4,6],labels:['x','y','z']}



    ;take dot product of E perp into the two perp unit vecs to find perp E in FAC
    E_perp_1  = fltarr(nef)
    for xx=0L,nef-1 do E_perp_1[xx] = dE_field[xx,0]*perp1_dir[xx,0] +  dE_field[xx,1]*perp1_dir[xx,1] +  dE_field[xx,2]*perp1_dir[xx,2]
    E_perp_2  = fltarr(nef)
    for xx=0L,nef-1 do E_perp_2[xx] = dE_field[xx,0]*perp2_dir[xx,0] +  dE_field[xx,1]*perp2_dir[xx,1] +  dE_field[xx,2]*perp2_dir[xx,2]

    edb_stat='E-dot-B-where-OK'
    store_data,'E-field!CAzimuthal(East)!C'+edb_stat+'!C'+smoothed_int+'-min-ave!CmV/m',data = {x:etimes,y:E_perp_1},dlim={constant:[0],colors:[0],labels:['']}
    store_data,'E-field!CRadial(Outward)!C'+edb_stat+'!C'+smoothed_int+'-min-ave!CmV/m',data = {x:etimes,y:E_perp_2},dlim={constant:[0],colors:[0],labels:['']}


    inteadt = total(E_perp_1,/nan,/cumulative)*find_datarate(etimes)

    store_data,'Time-Int!CE-azimuthal!CmV/m-s',data={x:etimes,y:inteadt},dlim={constant:[0]}



    ; Put ExB flow using E dot B = 0 to find Ex into FAC
    dV_perp_1  = fltarr(nef)
    for xx=0L,nef-1 do dV_perp_1[xx] = exb_flow[xx,0]*perp1_dir[xx,0] +  exb_flow[xx,1]*perp1_dir[xx,1] +  exb_flow[xx,2]*perp1_dir[xx,2]
    dV_perp_2  = fltarr(nef)
    for xx=0L,nef-1 do dV_perp_2[xx] = exb_flow[xx,0]*perp2_dir[xx,0] +  exb_flow[xx,1]*perp2_dir[xx,1] +  exb_flow[xx,2]*perp2_dir[xx,2]

    edb_stat='E-dot-B-where-OK'
    store_data,'TH'+sc+'ExB-flow!CAzimuthal(East)!C'+smoothed_int+'-min-ave!C'+edb_stat+'!Ckm/s',data = {x:etimes,y:dV_perp_1},dlim={constant:[0],colors:[0],labels:['']}
    store_data,'TH'+sc+'ExB-flow!CRadial(Outward)!C'+smoothed_int+'-min-ave!C'+edb_stat+'!Ckm/s',data = {x:etimes,y:dV_perp_2},dlim={constant:[0],colors:[0],labels:['']}



    title='FAC Perp1 is in the B cross X R direction (azimuthal) and!CPerp2 is in the perp1 cross B direction (radial)'
    print,title


    tplot,['TH'+sc+'ExB-flow!CAzimuthal(East)!C'+smoothed_int+'-min-ave!C'+edb_stat+'!Ckm/s',$
    'TH'+sc+'ExB-flow!CRadial(Outward)!C'+smoothed_int+'-min-ave!C'+edb_stat+'!Ckm/s']


    get_data,'tha_fgs_gse',data = B_gse
    get_data,'tha_vel_gse',data = v_gse

    vx = v_gse.y[*,0] & vy = v_gse.y[*,1] & vz = v_gse.y[*,2]

    vx=interp(vx,v_gse.x,B_gse.x,/no_extrap)
    vy=interp(vy,v_gse.x,B_gse.x,/no_extrap)
    vz=interp(vz,v_gse.x,B_gse.x,/no_extrap)

    bx = B_gse.y[*,0] & by = B_gse.y[*,1] & bz = B_gse.y[*,2]

    bmag = SQRT(bx^2+by^2+bz^2)
    vmag = SQRT(vx^2+vy^2+vz^2)

    bdotv  = (vx*bx+vy*by+vz*bz)/(bmag)

    store_data,'B-dot-V!isc!n!Cnormalized!Ckm/s',data={x:B_gse.x,y:bdotv}

;    tplot,['E-dot-B!Cflag','Vx!CExB-drift!Ckm/s','Vy!CExB-drift!CE-dot-B-where-OK!Ckm/s','Vz!CExB-drift!CE-dot-B-where-OK!Ckm/s',$
;    'ExB-flow!CAzimuthal(East)!CE-dot-B-where-OK!Ckm/s',$
;    'ExB-flow!CRadial(Outward)!CE-dot-B-where-OK!Ckm/s','E-dot-B!Cflag']


stop


    copy_data,'tha_fgs_gse_interp_smoothed','Bfield_for_ExB'
    copy_data,'tha_efs_dot0_gse_smoothed','Efield_for_ExB'

    savevars = ['Bfield_for_ExB','Efield_for_ExB','Vx!CExB-drift!Ckm/s']

    tplot_save,savevars,filename='~/Desktop/exb_20140111_tha'

end
