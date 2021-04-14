  
;input for doppler_shift_whistlermode_calculate2.pro for the STEREO burst captures



;pro doppler_shift_whistlermode_calculate2,Bo,vflow,fsc,dens,epol,plot_str=plot_str,coordsys=coordsys,plotdir=plotdir,nlevels=nlevels,plot_ps=plot_ps,titleroot=titleroot




;###############################
;For April 08, 2007 STA at 21:04:50.617
vflow = [-364.,35.,43.9]
Bo = [-0.21,2.6,-14.0]
;field_struct = read_tds(filenames='~/Desktop/swidl-0.1.3/tds_data/STEREO_A-TDS_b-08Apr2007_210450.617_SC.txt')
field_struct = read_tds(filenames='~/Desktop/swidl-0.1.3/tds_data/STEREO_A-TDS_b-08Apr2007_210524.023_SC.txt')
field= [[field_struct.ex_cor],[field_struct.ey_cor],[field_struct.ez_cor]]

fsc = 85.  ;freq in spacecraft frame (Hz)
dens = 10.9 ;cm-3
epol = 7.
EorB = 'E'



plot_struct = {plot_str:'',$
			   coordsys:inputcoord,$
			   plotdir:'~/Desktop/dstest/',$
			   titleroot:'-sc' + sc + '-' + datetime,$
			   datetime:'2007-04-08_210450.617',$
			   sc:'a'}




doppler_shift_whistlermode_calculate2,field,Bo,vflow,fsc,dens,EorB,epol=epol,/plot_ps,plot_struct=plot_struct








;###############################
;For April 08, 2007 STA at 22:05:54.957  (MOST CIRCULARLY POLARIZED WAVE IN EX-EY PLANE)
;inputcoord = 'sc'   
;sc = 'a'
;datetime = '2007-04-08_220554.957'
;Vsw = [-403.5,8.3,51.5]
;B = [-10.7,-2.4,-6.3]
;field_struct = read_tds(filenames='~/Desktop/swidl-0.1.3/tds_data/STEREO_A-TDS_b-08Apr2007_220554.957_SC.txt')
;fsc = 99.9  ;freq in spacecraft frame (Hz)
;dens = 6.9 ;cm-3
;epol = 1.2
;###################################
;For April 22, 2007 STA at 04:07:44.648
;inputcoord = 'sc'   
;sc = 'a'
;datetime = '2007-04-22_040744.648'
;Vsw = [-325.0,-1.7,7.6]
;B = [-2.5,0.8,-5.8]
;field_struct = read_tds(filenames='~/Desktop/swidl-0.1.3/tds_data/STEREO_A-TDS_b-22Apr2007_040744.648_SC.txt')
;fsc = 32.9  ;freq in spacecraft frame (Hz)
;dens = 16.4 ;cm-3
;epol = 1 ;temporary
;###################################
;For April 22, 2007 STA at 11:10:49.547
;inputcoord = 'sc'   
;sc = 'a'
;datetime = '2007-04-22_111049.547'
;Vsw = [-356.0,-2.5,30.9]
;B = [7.8,-2.9,-0.9]
;field_struct = read_tds(filenames='~/Desktop/swidl-0.1.3/tds_data/STEREO_A-TDS_b-22Apr2007_111049.547_SC.txt')
;fsc = 35.3  ;freq in spacecraft frame (Hz)
;dens = 30.9 ;cm-3
;epol = 1 ;temporary
;###############################
;For April 27, 2007 STA at 23:29:47.762
;inputcoord = 'sc'   
;sc = 'a'
;datetime = '2007-04-27_232947.762'
;B = [4.43,-1.59,7.34]
;Vsw = [-566.7,-49.6,-30.0]
;field_struct = read_tds(filenames='~/Desktop/swidl-0.1.3/tds_data/STEREO_A-TDS_b-27Apr2007_232947.762_SC.txt')
;fsc = 38.5  ;freq in spacecraft frame (Hz)
;dens = 11.6 ;cm-3
;epol= 5.9
;##################################
;For May 07, 2007 STA at 12:15:44.730
;inputcoord = 'sc'   
;sc = 'a'
;datetime = '2007-05-07_121544.730'
;Vsw = [-332.1,15.9,46.6]
;B = [-8.0,-7.3,-6.6]
;field_struct = read_tds(filenames='~/Desktop/swidl-0.1.3/tds_data/STEREO_A-TDS_b-07May2007_121544.730_SC.txt')
;fsc = 34.6  ;freq in spacecraft frame (Hz)
;dens = 41.9 ;cm-3
;epol = 3.1
;##################################
;For May 18, 2007 STA at 16:18:45.711
;inputcoord = 'sc'   
;sc = 'a'
;datetime = '2007-05-18_161845.711'
;Vsw = [-469.7,-36.1,-31.8]
;B = [-9.2,6.4,-8.5]
;field_struct = read_tds(filenames='~/Desktop/swidl-0.1.3/tds_data/STEREO_A-TDS_b-18May2007_161845.711_SC.txt')
;fsc = 65.9  ;freq in spacecraft frame (Hz)
;dens = 16.8 ;cm-3
;epol = 1 ;temporary
;###############################
;For Jun 03, 2007 STA at 06:22:40.621 - NO DATA
;inputcoord = 'sc'   
;sc = 'a'
;datetime = '2007-06-10_062240.621'
;Vsw = []
;B = []
;field_struct = read_tds(filenames='~/Desktop/swidl-0.1.3/tds_data/STEREO_A-TDS_b-03Jun2007_062240.621_SC.txt')
;fsc =   ;freq in spacecraft frame (Hz)
;dens =  ;cm-3
;epol=  ;temporary
;###############################
;For Jun 04, 2007 STA at 14:27:59.176 - NO DATA
;inputcoord = 'sc'   
;sc = 'a'
;datetime = '2007-06-04_142759.176'
;Vsw = []
;B = []
;field_struct = read_tds(filenames='~/Desktop/swidl-0.1.3/tds_data/STEREO_A-TDS_b-04Jun2007_142759.176_SC.txt')
;fsc =   ;freq in spacecraft frame (Hz)
;dens =  ;cm-3
;epol=  ;temporary
;###############################
;For Jun 09, 2007 STA at 15:31:51.852
;inputcoord = 'sc'   
;sc = 'a'
;datetime = '2007-06-09_153151.852'
;Vsw = [-364.8,-11.,43.7]
;B = [4.8,4.0,3.3]
;field_struct = read_tds(filenames='~/Desktop/swidl-0.1.3/tds_data/STEREO_A-TDS_b-09Jun2007_153151.852_SC.txt')
;fsc = 13.8 ;freq in spacecraft frame (Hz)
;dens = 13.3 ;cm-3
;epol= 1. ;temporary
;###############################
;For Jun 10, 2007 STA at 12:53:46.160
;inputcoord = 'sc'   
;sc = 'a'
;datetime = '2007-06-10_125346.160'
;Vsw = [-439.9,-1.1,7.1]
;B = [1.0,-0.4,4.5]
;field_struct = read_tds(filenames='~/Desktop/swidl-0.1.3/tds_data/STEREO_A-TDS_b-10Jun2007_125346.160_SC.txt')
;fsc = 46.6  ;freq in spacecraft frame (Hz)
;dens = 3.1 ;cm-3
;epol = 1 ;temporary
;###############################
;For Jun 14, 2007 STA at 02:58:18.125
;inputcoord = 'sc'   
;sc = 'a'
;datetime = '2007-06-14_025818.125'
;Vsw = [-325.8,11.4,16.9]
;B = [-4.2,-1.1,-6.8]
;field_struct = read_tds(filenames='~/Desktop/swidl-0.1.3/tds_data/STEREO_A-TDS_b-14Jun2007_025818.125_SC.txt')
;fsc = 49.9  ;freq in spacecraft frame (Hz)
;dens = 7.8 ;cm-3
;epol=3.9  ;temporary
;###############################
;For Jun 22, 2007 STA at 18:45:42.793
;inputcoord = 'sc'   
;sc = 'a'
;datetime = '2007-06-22_184542.793'
;Vsw = [-526.,-17.7,-33.4]
;B = [-2.9,-7.6,-2.8]
;field_struct = read_tds(filenames='~/Desktop/swidl-0.1.3/tds_data/STEREO_A-TDS_b-22Jun2007_184542.793_SC.txt')
;fsc = 38.5  ;freq in spacecraft frame (Hz)
;dens = 8.7 ;cm-3
;epol=5.0
;###############################
;For Jun 29, 2007 STA at 19:22:01.363
;inputcoord = 'sc'   
;sc = 'a'
;datetime = '2007-06-29_192201.363'
;Vsw = [-332.6,-7.7,28.9]
;B = [1.6,-5.3,1.8]
;field_struct = read_tds(filenames='~/Desktop/swidl-0.1.3/tds_data/STEREO_A-TDS_b-29Jun2007_192201.363_SC.txt')
;fsc = 34.6  ;freq in spacecraft frame (Hz)
;dens = 24.5 ;cm-3
;epol=5.4
;#################################
;For Jul 11, 2007 STA at 11:43:05.980
;inputcoord = 'sc'   
;sc = 'a'
;datetime = '2007-07-11_114305.980'
;Vsw = [-417.8,28.5,3.3]
;B = [-7.0,-11.5,-8.5]
;field_struct = read_tds(filenames='~/Desktop/swidl-0.1.3/tds_data/STEREO_A-TDS_b-11Jul2007_114305.980_SC.txt')
;fsc = 76.9  ;freq in spacecraft frame (Hz)
;dens = 28.8 ;cm-3
;epol=5.4
;###############################

;###############################
;For Jun 30, 2007 STB at 04:07:07.992
;NO VELOCITY DATA FOR ENTIRE DAY!!!
;inputcoord = 'sc'   
;sc = 'b'
;datetime = '2007-06-30_040707.992'
;Vsw = []
;B = []
;field_struct = read_tds(filenames='~/Desktop/swidl-0.1.3/tds_data/STEREO_B-TDS_b-30Jun2007_040707.992_SC.txt')
;fsc = 46.2  ;freq in spacecraft frame (Hz)
;dens =  ;cm-3
;###############################


;#################################
;SETUP FOR WHISTLER EVENT ON WIND
;For Mar 11, 2007 WIND at 11:08:18.164
;restore,'~/Desktop/code/Aaron/datafiles/wind_whistler_20070311_110818x164.sav'
;inputcoord = 'gse'   
;sc = 'w'
;datetime = '2007-03-11_110818.164'
;Vsw = vsw_gse
;B = bkg_field
;field_struct = bs
;fsc = 50.  ;freq in spacecraft frame (Hz)
;dens = 20. ;cm-3   ;TEMPORARY!!!!!!!!!!!
;epol=1.           ;TEMPORARY!!!!!!!!!!!!!
;###############################


;#################################
;SETUP FOR 2006-11-06 NPM event
;For Nov 06, 2006 STEREO at 09:07:01.898
;inputcoord = 'sc'   
;sc = 'a'
;datetime = '2006-11-06_090701.898'
;Vsc = []
;Vsw = Vsc ;determine the DS due to sc motion
;B = 26147.
;field_struct = read_tds(filenames=~/Desktop/Research/inner_ps_density_cavity/ascii_files/SCTDS_coord/sta_nov06/STEREO_A-TDS_b-06Nov2006_090701.898_SC.txt')
;fsc = 24100.  ;freq in spacecraft frame (Hz)
;dens = 2000. ;cm-3   ;TEMPORARY!!!!!!!!!!!
;epol=1.           ;TEMPORARY!!!!!!!!!!!!!
;###############################










;--------------------------
;PS PLOTS
;--------------------------

if plot_ps eq 'yes' then begin

;------ plot plasma frame freq

set_plot,'ps'
device,filename = plot_struct.plotdir + 'pframe_freq' + plot_struct.titleroot + '-Dopplershift' + '.ps',/color,_extra=winspecs
	!p.multi = [0,0,1,0,0]

	cs = bytscl(f_pframe2,min=min(f_pframe2,/nan),max=max(f_pframe2,/nan),/nan)
	toob = where(cs eq 255.)
	if toob[0] ne -1 then cs[toob] = 254.
	polar_contour,cs*multiplier,theta_kb_gen*!dtor,kmag_gen,levels=levels,/fill,position=[0.10,0.15,0.80,0.90], $
		title='z=Freq_pframe (Hz) (from DS relation), theta=theta_kb, r=|k|',$  ;max_value=(fce-1.),min_value = 0., $
		xrange=[-1*max(kmag_gen),max(kmag_gen)],_extra=xwinc1,/isotropic
	for i=0,3 do oplot,quadrants[i]*dispvals2_x[*,i],quadrants[i]*dispvals2_y[*,i],color=220,thick=4
	oplot,xvals3,yvals3
	oplot,xvals4,yvals4
	oplot,[0,vs2[0]],[0,vs2[1]],thick=3,linestyle=2
	oplot,[0,bs2[0]*fce],[0,bs2[1]*fce],thick=3,linestyle=2
	oplot,[0,fmax2_1[0]],[0,fmax2_1[1]],thick=3,color=50,linestyle=2
	oplot,[0,fmax2_2[0]],[0,fmax2_2[1]],thick=3,color=50,linestyle=2
	colorbar,position=[0.85,0.15,0.90,0.90],/vertical,/right ,range=[min(f_pframe2,/nan),max(f_pframe2,/nan)],charsize=1.2
	plotspecs,struct,fmax_stix1=struct.fmax_stix1,fmax_stix2=struct.fmax_stix2,fint_stix=fint_stix
device,/close
set_plot,'x'
;------ plot phase velocity 

set_plot,'ps'
device,filename = plot_struct.plotdir + 'vphase' + plot_struct.titleroot + '-Dispersion' + '.ps',/color,_extra=winspecs
	!p.multi = [0,0,1,0,0]
	polar_contour,vphase,theta_kb_gen*!dtor,freqs_gen,c_colors=c_colors,$
		levels=vphase_lvls,position=[0.10,0.15,0.80,0.90],/cell_fill,title = 'z=vphase (km/s), theta=theta_kb, r=f_pframe',$
		xrange=[-fce,fce],xstyle=4,ystyle=4,_extra=psc1,/isotropic
	for i=0,3 do oplot,quadrants[i]*dispvals_x[*,i],quadrants[i]*dispvals_y[*,i],color=220,thick=4
	oplot,xvals1,yvals1
	oplot,xvals2,yvals2
	oplot,ftmp*cos(tr1*!dtor),ftmp*sin(tr1*!dtor),thick=2
	oplot,ftmp*cos(tr2*!dtor),ftmp*sin(tr2*!dtor),thick=2
	oplot,ftmp*cos(tr3*!dtor),ftmp*sin(tr3*!dtor),thick=2
	oplot,ftmp*cos(tr4*!dtor),ftmp*sin(tr4*!dtor),thick=2
	oplot,[0,vs2[0]*fce],[0,vs2[1]*fce],thick=3,linestyle=2
	oplot,[0,bs2[0]*fce],[0,bs2[1]*fce],thick=3,linestyle=4
	oplot,[0,fmax2_1[0]*fce],[0,fmax2_1[1]*fce],thick=3,color=50,linestyle=2
	oplot,[0,fmax2_2[0]*fce],[0,fmax2_2[1]*fce],thick=3,color=50,linestyle=2
	colorbar,ncolors=255,position=[0.85,0.15,0.90,0.90],/vertical,/right,range=[min(vphase_lvls),max(vphase_lvls)] 
	plotspecs,struct,fmax_stix1=struct.fmax_stix1,fmax_stix2=struct.fmax_stix2,fint_stix=fint_stix
device,/close
set_plot,'x'
;------ plot |k| for spectrum of freqs in plasma frame and wave normal angles
set_plot,'ps'
device,filename = plot_struct.plotdir + 'kmag_gen' + plot_struct.titleroot + '-Dispersion' + '.ps',/color,_extra=winspecs
	!p.multi = [0,0,1,0,0]
	polar_contour,kvalues,theta_kb_gen*!dtor,freqs_gen,levels=klvls,/cell_fill,position=[0.10,0.15,0.80,0.90], $
		xtitle=tkb,ytitle='freq (plasma frame) (Hz)',title='z=|k| (1/km), theta=theta_kb, r=f_pframe',charsize=1.1,xrange=[-fce,fce],xstyle=4,ystyle=4,_extra=psc1,/isotropic
	for i=0,3 do oplot,quadrants[i]*dispvals_x[*,i],quadrants[i]*dispvals_y[*,i],color=220,thick=4
	oplot,xvals1,yvals1
	oplot,xvals2,yvals2
	oplot,ftmp*cos(tr1*!dtor),ftmp*sin(tr1*!dtor),thick=2
	oplot,ftmp*cos(tr2*!dtor),ftmp*sin(tr2*!dtor),thick=2
	oplot,ftmp*cos(tr3*!dtor),ftmp*sin(tr3*!dtor),thick=2
	oplot,ftmp*cos(tr4*!dtor),ftmp*sin(tr4*!dtor),thick=2
	oplot,[0,vs2[0]*fce],[0,vs2[1]*fce],thick=3,linestyle=2
	oplot,[0,bs2[0]*fce],[0,bs2[1]*fce],thick=3,linestyle=4
	oplot,[0,fmax2_1[0]*fce],[0,fmax2_1[1]*fce],thick=3,color=50,linestyle=2
	oplot,[0,fmax2_2[0]*fce],[0,fmax2_2[1]*fce],thick=3,color=50,linestyle=2
	colorbar,ncolors=255,/vertical,range=[min(klvls),max(klvls)],/right,position=[0.85,0.15,0.90,0.90],format='(F5.2)'
	plotspecs,struct,fmax_stix1=struct.fmax_stix1,fmax_stix2=struct.fmax_stix2,fint_stix=fint_stix
device,/close
set_plot,'x'
;------ plot cyclotron res energy
set_plot,'ps'
device,filename = plot_struct.plotdir + 'cycl_res' + plot_struct.titleroot + '-Dispersion' + '.ps',/color,_extra=winspecs
	!p.multi = [0,0,1,0,0]
	polar_contour,Ecycl,theta_kb_gen*!dtor,freqs_gen,c_colors=c_colors,levels=ecycl_lvls,position=[0.10,0.15,0.80,0.90], $
		xtitle=tkb,ytitle='Freq (pframe) '+'(Hz)',title='z=Ecycl (eV), theta=theta_kb, r=f_pframe',/cell_fill,charsize=1.1,xrange=[-fce,fce],xstyle=4,ystyle=4,_extra=psc1,/isotropic
	for i=0,3 do oplot,quadrants[i]*dispvals_x[*,i],quadrants[i]*dispvals_y[*,i],color=220,thick=4
	oplot,xvals1,yvals1
	oplot,xvals2,yvals2
	oplot,ftmp*cos(tr1*!dtor),ftmp*sin(tr1*!dtor),thick=2
	oplot,ftmp*cos(tr2*!dtor),ftmp*sin(tr2*!dtor),thick=2
	oplot,ftmp*cos(tr3*!dtor),ftmp*sin(tr3*!dtor),thick=2
	oplot,ftmp*cos(tr4*!dtor),ftmp*sin(tr4*!dtor),thick=2
	oplot,[0,vs2[0]*fce],[0,vs2[1]*fce],thick=3,linestyle=2
	oplot,[0,bs2[0]*fce],[0,bs2[1]*fce],thick=3,linestyle=4
	oplot,[0,fmax2_1[0]*fce],[0,fmax2_1[1]*fce],thick=3,color=50,linestyle=2
	oplot,[0,fmax2_2[0]*fce],[0,fmax2_2[1]*fce],thick=3,color=50,linestyle=2
	colorbar,ncolors=255,/vertical,/right,range=[min(ecycl_lvls),max(ecycl_lvls)],position=[0.85,0.15,0.90,0.90] 
device,/close
set_plot,'x'
;------ plot Landau res energy
set_plot,'ps'
device,filename = plot_struct.plotdir + 'landau_res' + plot_struct.titleroot + '-Dispersion' + '.ps',/color,_extra=winspecs
	!p.multi = [0,0,1,0,0]
	polar_contour,Elandau,theta_kb_gen*!dtor,freqs_gen,c_colors=c_colors,levels=elandau_lvls,position=[0.10,0.15,0.80,0.90], $
		xtitle=tkb,ytitle='Freq (pframe) '+'(Hz)',title='z=Elandau (eV), theta=theta_kb, r=f_pframe',/cell_fill,charsize=1.1,xrange=[-fce,fce],xstyle=4,ystyle=4,_extra=psc1,/isotropic
	for i=0,3 do oplot,quadrants[i]*dispvals_x[*,i],quadrants[i]*dispvals_y[*,i],color=220,thick=4
	oplot,xvals1,yvals1
	oplot,xvals2,yvals2
	oplot,ftmp*cos(tr1*!dtor),ftmp*sin(tr1*!dtor),thick=2
	oplot,ftmp*cos(tr2*!dtor),ftmp*sin(tr2*!dtor),thick=2
	oplot,ftmp*cos(tr3*!dtor),ftmp*sin(tr3*!dtor),thick=2
	oplot,ftmp*cos(tr4*!dtor),ftmp*sin(tr4*!dtor),thick=2
	oplot,[0,vs2[0]*fce],[0,vs2[1]*fce],thick=3,linestyle=2
	oplot,[0,bs2[0]*fce],[0,bs2[1]*fce],thick=3,linestyle=4
	oplot,[0,fmax2_1[0]*fce],[0,fmax2_1[1]*fce],thick=3,color=50,linestyle=2
	oplot,[0,fmax2_2[0]*fce],[0,fmax2_2[1]*fce],thick=3,color=50,linestyle=2
	colorbar,ncolors=255,/vertical,/right,range=[min(elandau_lvls),max(elandau_lvls)],position=[0.85,0.15,0.90,0.90]
device,/close
set_plot,'x'
;----- plot anomalous res energy 
set_plot,'ps'
device,filename = plot_struct.plotdir + 'anom_res' + plot_struct.titleroot + '-Dispersion' + '.ps',/color,_extra=winspecs
	!p.multi = [0,0,1,0,0]
	polar_contour,Eanom,theta_kb_gen*!dtor,freqs_gen,c_colors=c_colors,levels=eanom_lvls,position=[0.10,0.15,0.80,0.90], $
		xtitle=tkb,ytitle='Freq (pframe) '+'(Hz)',title='z=Eanom res energy (eV), theta=theta_kb, r=f_pframe',/cell_fill,charsize=1.1,ystyle=4,xrange=[-fce,fce],xstyle=4,_extra=psc1,/isotropic
	for i=0,3 do oplot,quadrants[i]*dispvals_x[*,i],quadrants[i]*dispvals_y[*,i],color=220,thick=4
	oplot,xvals1,yvals1
	oplot,xvals2,yvals2
	oplot,ftmp*cos(tr1*!dtor),ftmp*sin(tr1*!dtor),thick=2
	oplot,ftmp*cos(tr2*!dtor),ftmp*sin(tr2*!dtor),thick=2
	oplot,ftmp*cos(tr3*!dtor),ftmp*sin(tr3*!dtor),thick=2
	oplot,ftmp*cos(tr4*!dtor),ftmp*sin(tr4*!dtor),thick=2
	oplot,[0,vs2[0]*fce],[0,vs2[1]*fce],thick=3,linestyle=2
	oplot,[0,bs2[0]*fce],[0,bs2[1]*fce],thick=3,linestyle=4
	oplot,[0,fmax2_1[0]*fce],[0,fmax2_1[1]*fce],thick=3,color=50,linestyle=2
	oplot,[0,fmax2_2[0]*fce],[0,fmax2_2[1]*fce],thick=3,color=50,linestyle=2
	colorbar,ncolors=255,/vertical,/right,range=[min(eanom_lvls),max(eanom_lvls)],position=[0.85,0.15,0.90,0.90]
device,/close
set_plot,'x'

endif








;--------------------------
;PS PLOTS - TO THE Z-BUFFER
;--------------------------

plot_z = 'no'
if plot_z eq 'yes' then begin

zres = [1000,1000]
thisDevice = !d.name


;------ plot plasma frame freq
set_plot,'z',/copy
device,set_resolution=zres,z_buffer=0,_extra=winspecs
device,_extra=winspecs
erase
!p.position=[0,0,1,1]                                  

!p.multi = [0,0,1,0,0]
polar_contour,f_pframe2*multiplier,theta_kb_gen*!dtor,kmag_gen,c_colors=c_colors,levels=pframe_lvls,/cell_fill,position=[0.10,0.15,0.80,0.90], $
	title='z=Freq_pframe (Hz) (from DS relation), theta=theta_kb, r=|k|',max_value=(fce-1.),min_value = 0., $
	xrange=[-1*max(kmag_gen)/2,max(kmag_gen)/2],xstyle=4,ystyle=4,/isotropic,_extra=psc1
for i=0,3 do oplot,quadrants[i]*dispvals2_x[*,i],quadrants[i]*dispvals2_y[*,i],color=220,thick=4
oplot,xvals3,yvals3
oplot,xvals4,yvals4
oplot,[0,vs2[0]],[0,vs2[1]],thick=3,linestyle=2
oplot,[0,bs2[0]*fce],[0,bs2[1]*fce],thick=3,linestyle=2
oplot,[0,fmax2_1[0]],[0,fmax2_1[1]],thick=3,color=50,linestyle=2
oplot,[0,fmax2_2[0]],[0,fmax2_2[1]],thick=3,color=50,linestyle=2
colorbar,ncolors=255,position=[0.85,0.15,0.90,0.90],/vertical,range=[min(pframe_lvls),max(pframe_lvls)],/right

a=TVRD()
device,/close
set_plot,thisDevice ;return device to previous

plotspecs,struct,fmax_stix1=struct.fmax_stix1,fmax_stix2=struct.fmax_stix2,fint_stix=fint_stix

set_plot,'ps'
device,filename = plot_struct.plotdir + 'pframe_freq' + plot_struct.titleroot + '-Dopplershift' + '.ps',/color,_extra=winspecs
!p.multi = [0,0,1,0,0]
tvimage,a,Position=Aspect(1.2/1.0)
plotspecs,struct,fmax_stix1=struct.fmax_stix1,fmax_stix2=struct.fmax_stix2,fint_stix=fint_stix
device,/close
set_plot,'x'
;------ plot phase velocity 
set_plot,'z',/copy
device,set_resolution=zres,z_buffer=0
erase
!p.position=[0,0,1,1]                                  
	!p.multi = [0,0,1,0,0]
	polar_contour,vphase,theta_kb_gen*!dtor,freqs_gen,c_colors=c_colors,$
		levels=vphase_lvls,position=[0.10,0.15,0.80,0.90],/cell_fill,title = 'z=vphase (km/s), theta=theta_kb, r=f_pframe',$
		xrange=[-fce,fce],xstyle=4,ystyle=4,_extra=psc1,/isotropic
	for i=0,3 do oplot,quadrants[i]*dispvals_x[*,i],quadrants[i]*dispvals_y[*,i],color=220,thick=4
	oplot,xvals1,yvals1
	oplot,xvals2,yvals2
	oplot,ftmp*cos(tr1*!dtor),ftmp*sin(tr1*!dtor),thick=2
	oplot,ftmp*cos(tr2*!dtor),ftmp*sin(tr2*!dtor),thick=2
	oplot,ftmp*cos(tr3*!dtor),ftmp*sin(tr3*!dtor),thick=2
	oplot,ftmp*cos(tr4*!dtor),ftmp*sin(tr4*!dtor),thick=2
	oplot,[0,vs2[0]*fce],[0,vs2[1]*fce],thick=3,linestyle=2
	oplot,[0,bs2[0]*fce],[0,bs2[1]*fce],thick=3,linestyle=4
	oplot,[0,fmax2_1[0]*fce],[0,fmax2_1[1]*fce],thick=3,color=50,linestyle=2
	oplot,[0,fmax2_2[0]*fce],[0,fmax2_2[1]*fce],thick=3,color=50,linestyle=2
	colorbar,ncolors=255,position=[0.85,0.15,0.90,0.90],/vertical,/right,range=[min(vphase_lvls),max(vphase_lvls)] 

a=TVRD() 
device,/close
set_plot,thisDevice ;return device to previous


set_plot,'ps'
device,filename = plot_struct.plotdir + 'vphase' + plot_struct.titleroot + '-Dispersion' + '.ps',/color
	tvimage,a
	plotspecs,struct,fmax_stix1=struct.fmax_stix1,fmax_stix2=struct.fmax_stix2,fint_stix=fint_stix
device,/close
set_plot,'x'
;------ plot |k| for spectrum of freqs in plasma frame and wave normal angles
set_plot,'z',/copy
device,set_resolution=zres,z_buffer=0
erase
!p.position=[0,0,1,1]                                  

	!p.multi = [0,0,1,0,0]
	polar_contour,kvalues,theta_kb_gen*!dtor,freqs_gen,levels=klvls,/cell_fill,position=[0.10,0.15,0.80,0.90], $
		xtitle=tkb,ytitle='freq (plasma frame) (Hz)',title='z=|k| (1/km), theta=theta_kb, r=f_pframe',charsize=1.1,xrange=[-fce,fce],xstyle=4,ystyle=4,_extra=psc1,/isotropic
	for i=0,3 do oplot,quadrants[i]*dispvals_x[*,i],quadrants[i]*dispvals_y[*,i],color=220,thick=4
	oplot,xvals1,yvals1
	oplot,xvals2,yvals2
	oplot,ftmp*cos(tr1*!dtor),ftmp*sin(tr1*!dtor),thick=2
	oplot,ftmp*cos(tr2*!dtor),ftmp*sin(tr2*!dtor),thick=2
	oplot,ftmp*cos(tr3*!dtor),ftmp*sin(tr3*!dtor),thick=2
	oplot,ftmp*cos(tr4*!dtor),ftmp*sin(tr4*!dtor),thick=2
	oplot,[0,vs2[0]*fce],[0,vs2[1]*fce],thick=3,linestyle=2
	oplot,[0,bs2[0]*fce],[0,bs2[1]*fce],thick=3,linestyle=4
	oplot,[0,fmax2_1[0]*fce],[0,fmax2_1[1]*fce],thick=3,color=50,linestyle=2
	oplot,[0,fmax2_2[0]*fce],[0,fmax2_2[1]*fce],thick=3,color=50,linestyle=2
	colorbar,ncolors=255,/vertical,range=[min(klvls),max(klvls)],/right,position=[0.85,0.15,0.90,0.90],format='(F5.2)'

a=TVRD() 
device,/close
set_plot,thisDevice ;return device to previous

set_plot,'ps'
device,filename = plot_struct.plotdir + 'kmag_gen' + plot_struct.titleroot + '-Dispersion' + '.ps',/color
	tvimage,a
	plotspecs,struct,fmax_stix1=struct.fmax_stix1,fmax_stix2=struct.fmax_stix2,fint_stix=fint_stix
device,/close
set_plot,'x'
;------ plot cyclotron res energy
set_plot,'z',/copy
device,set_resolution=zres,z_buffer=0
erase
!p.position=[0,0,1,1]                                  
	!p.multi = [0,0,1,0,0]
	polar_contour,Ecycl,theta_kb_gen*!dtor,freqs_gen,c_colors=c_colors,levels=ecycl_lvls,position=[0.10,0.15,0.80,0.90], $
		xtitle=tkb,ytitle='Freq (pframe) '+'(Hz)',title='z=Ecycl (eV), theta=theta_kb, r=f_pframe',/cell_fill,charsize=1.1,xrange=[-fce,fce],xstyle=4,ystyle=4,_extra=psc1,/isotropic
	for i=0,3 do oplot,quadrants[i]*dispvals_x[*,i],quadrants[i]*dispvals_y[*,i],color=220,thick=4
	oplot,xvals1,yvals1
	oplot,xvals2,yvals2
	oplot,ftmp*cos(tr1*!dtor),ftmp*sin(tr1*!dtor),thick=2
	oplot,ftmp*cos(tr2*!dtor),ftmp*sin(tr2*!dtor),thick=2
	oplot,ftmp*cos(tr3*!dtor),ftmp*sin(tr3*!dtor),thick=2
	oplot,ftmp*cos(tr4*!dtor),ftmp*sin(tr4*!dtor),thick=2
	oplot,[0,vs2[0]*fce],[0,vs2[1]*fce],thick=3,linestyle=2
	oplot,[0,bs2[0]*fce],[0,bs2[1]*fce],thick=3,linestyle=4
	oplot,[0,fmax2_1[0]*fce],[0,fmax2_1[1]*fce],thick=3,color=50,linestyle=2
	oplot,[0,fmax2_2[0]*fce],[0,fmax2_2[1]*fce],thick=3,color=50,linestyle=2
	colorbar,ncolors=255,/vertical,/right,range=[min(ecycl_lvls),max(ecycl_lvls)],position=[0.85,0.15,0.90,0.90] 

a=TVRD() 
device,/close
set_plot,thisDevice ;return device to previous

set_plot,'ps'
device,filename = plot_struct.plotdir + 'cycl_res' + plot_struct.titleroot + '-Dispersion' + '.ps',/color
	tvimage,a
	plotspecs,struct,fmax_stix1=struct.fmax_stix1,fmax_stix2=struct.fmax_stix2,fint_stix=fint_stix
device,/close
set_plot,'x'
;------ plot Landau res energy
set_plot,'z',/copy
device,set_resolution=zres,z_buffer=0
erase
!p.position=[0,0,1,1]                                  
	!p.multi = [0,0,1,0,0]
	polar_contour,Elandau,theta_kb_gen*!dtor,freqs_gen,c_colors=c_colors,levels=elandau_lvls,position=[0.10,0.15,0.80,0.90], $
		xtitle=tkb,ytitle='Freq (pframe) '+'(Hz)',title='z=Elandau (eV), theta=theta_kb, r=f_pframe',/cell_fill,charsize=1.1,xrange=[-fce,fce],xstyle=4,ystyle=4,_extra=psc1,/isotropic
	for i=0,3 do oplot,quadrants[i]*dispvals_x[*,i],quadrants[i]*dispvals_y[*,i],color=220,thick=4
	oplot,xvals1,yvals1
	oplot,xvals2,yvals2
	oplot,ftmp*cos(tr1*!dtor),ftmp*sin(tr1*!dtor),thick=2
	oplot,ftmp*cos(tr2*!dtor),ftmp*sin(tr2*!dtor),thick=2
	oplot,ftmp*cos(tr3*!dtor),ftmp*sin(tr3*!dtor),thick=2
	oplot,ftmp*cos(tr4*!dtor),ftmp*sin(tr4*!dtor),thick=2
	oplot,[0,vs2[0]*fce],[0,vs2[1]*fce],thick=3,linestyle=2
	oplot,[0,bs2[0]*fce],[0,bs2[1]*fce],thick=3,linestyle=4
	oplot,[0,fmax2_1[0]*fce],[0,fmax2_1[1]*fce],thick=3,color=50,linestyle=2
	oplot,[0,fmax2_2[0]*fce],[0,fmax2_2[1]*fce],thick=3,color=50,linestyle=2
	colorbar,ncolors=255,/vertical,/right,range=[min(elandau_lvls),max(elandau_lvls)],position=[0.85,0.15,0.90,0.90]

a=TVRD() 
device,/close
set_plot,thisDevice ;return device to previous

set_plot,'ps'
device,filename = plot_struct.plotdir + 'landau_res' + plot_struct.titleroot + '-Dispersion' + '.ps',/color
	tvimage,a
	plotspecs,struct,fmax_stix1=struct.fmax_stix1,fmax_stix2=struct.fmax_stix2,fint_stix=fint_stix
device,/close
set_plot,'x'
;----- plot anomalous res energy 
set_plot,'z',/copy
device,set_resolution=zres,z_buffer=0
erase
!p.position=[0,0,1,1]                                  
	!p.multi = [0,0,1,0,0]
	polar_contour,Eanom,theta_kb_gen*!dtor,freqs_gen,c_colors=c_colors,levels=eanom_lvls,position=[0.10,0.15,0.80,0.90], $
		xtitle=tkb,ytitle='Freq (pframe) '+'(Hz)',title='z=Eanom res energy (eV), theta=theta_kb, r=f_pframe',/cell_fill,charsize=1.1,ystyle=4,xrange=[-fce,fce],xstyle=4,_extra=psc1,/isotropic
	for i=0,3 do oplot,quadrants[i]*dispvals_x[*,i],quadrants[i]*dispvals_y[*,i],color=220,thick=4
	oplot,xvals1,yvals1
	oplot,xvals2,yvals2
	oplot,ftmp*cos(tr1*!dtor),ftmp*sin(tr1*!dtor),thick=2
	oplot,ftmp*cos(tr2*!dtor),ftmp*sin(tr2*!dtor),thick=2
	oplot,ftmp*cos(tr3*!dtor),ftmp*sin(tr3*!dtor),thick=2
	oplot,ftmp*cos(tr4*!dtor),ftmp*sin(tr4*!dtor),thick=2
	oplot,[0,vs2[0]*fce],[0,vs2[1]*fce],thick=3,linestyle=2
	oplot,[0,bs2[0]*fce],[0,bs2[1]*fce],thick=3,linestyle=4
	oplot,[0,fmax2_1[0]*fce],[0,fmax2_1[1]*fce],thick=3,color=50,linestyle=2
	oplot,[0,fmax2_2[0]*fce],[0,fmax2_2[1]*fce],thick=3,color=50,linestyle=2
	colorbar,ncolors=255,/vertical,/right,range=[min(eanom_lvls),max(eanom_lvls)],position=[0.85,0.15,0.90,0.90]

a=TVRD() 
device,/close
set_plot,thisDevice ;return device to previous

set_plot,'ps'
device,filename = plot_struct.plotdir + 'anom_res' + plot_struct.titleroot + '-Dispersion' + '.ps',/color
	tvimage,a
	plotspecs,struct,fmax_stix1=struct.fmax_stix1,fmax_stix2=struct.fmax_stix2,fint_stix=fint_stix
device,/close
set_plot,'x'






end