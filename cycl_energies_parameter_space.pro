;+
; NAME: cycl_energies_parameter_space
; SYNTAX:
; PURPOSE: create a parameter space plot of relativistic cyclotron resonance energies as a function
; of wave freq, fce, density
; INPUT: (see cycl_energies_parameter_space_call.pro.
; OUTPUT:
; NOTES: Tested against Lorentzen01, plate7 results. See cycl_energies_test.pro
; KEYWORDS:
;          pa -> e- pitch angle (deg)
;          theta_k -> wave normal angle (deg)
;          scheme -> 0, 1 or 2
;                scheme 0:  fce vs density (at constant wave freq)
;                scheme 1:  fce vs freq (at constant density)
;                scheme 2:  density vs freq (at constant fce)
;          ps -> plot to postscript (outputs to desktop)
;          maxval -> max energy. Sets colorbar resolution (keV)
;          minval -> min energy. Sets colorbar resolution (keV)
;          maxzval -> same but for Ez energy
;          minzval ->
;          density_range -> range of axis densities (cm-3)
;          fce_range     -> range of axis fce (Hz)
;          freq_range    -> range of axis frequencies (Hz or unitless
;                           if f_fce keyword is set)
;          ndens; nfce; nfreq -> number of densities, fces and freqs
;          (color resolution)
;          densv         -> value of input density
;          fcev          -> value of input fce (Hz)
;          freqv         -> value of input freq (Hz)
;          type          -> defaults to counter-streaming ('counterstream') cyclotron resonance. Other
;                           options are 'landau' and 'costream'
;          harmonic      -> |harmonic| of the resonance cyclotron resonance
;          f_fce -> set for scheme 2 to plot f/fce instead of f on x-axis
;
;
; REQUIRES: updated version of dfanning's colorbar.pro
; HISTORY: Written by AWB 2015-01-30
; VERSION:
;   $LastChangedBy: $
;   $LastChangedDate: $
;   $LastChangedRevision: $
;   $URL: $
;-



pro cycl_energies_parameter_space,pa,theta_k,$
  scheme=scheme,ps=ps,$
  minval=minval,maxval=maxval,$
  minzval=minzval,maxzval=maxzval,$
  density_range=density_range,$
  fce_range=fce_range,freq_range=freq_range,$
  ndens=ndens,nfce=nfce,nfreq=nfreq,$
  densv=dens,fcev=fce,freqv=freq,$
  type=type,harmonic=nres,f_fce=f_fce;,zlog=zlog


  rbsp_efw_init

  if ~KEYWORD_SET(dens) then dens = 1  ;means density not used
  if ~KEYWORD_SET(nres) then nres = 1
  if ~keyword_set(type) then type = 'counterstream'
  if keyword_set(density_range) then density_range = float(density_range)
  if keyword_set(fce_range) then fce_range = float(fce_range)
  if keyword_set(freq_range) then freq_range = float(freq_range)

  ;Number of density and fce array elements (resolution)
  if ~keyword_set(scheme) then scheme = 0


  ;Choose upper colorbar range (keV) (don't waste the colors plotting energies above this)
  if ~keyword_set(maxval) then maxval = 500.                ;keV
  if ~keyword_set(minval) and KEYWORD_SET(zlog) then minval = 0.1                ;keV
  if ~keyword_set(maxzval) then maxzval = 500.                ;keV
  if ~keyword_set(minzval) and KEYWORD_SET(zlog) then minzval = 0.1                ;keV
  if ~keyword_set(minzval) and ~KEYWORD_SET(zlog) then minzval = 0.                ;keV


;  if ~KEYWORD_SET(zlog) then zlog = 0. else zlog = 1.


  ;Choose range of parameters
  if scheme eq 0 then begin
    if ~keyword_set(density_range) then density_range=[1.,50.] ;cm-3
    if ~keyword_set(fce_range) then fce_range = [1000.,6000.]  ;Hz

    ;Number of density and fce array elements (resolution)
    if ~keyword_set(ndens) then ndens = 20.
    if ~keyword_set(nfce) then nfce = 20.
    if ~keyword_set(freq) then freq = 100.  ;Hz

    dens = (density_range[1]-density_range[0])*indgen(ndens)/(ndens-1) + density_range[0]
    fce = (fce_range[1]-fce_range[0])*indgen(nfce)/(nfce-1) + fce_range[0]
  endif

  if scheme eq 1 then begin

    if ~keyword_set(freq_range) then freq_range = [80.,300.]       ;Hz
    if ~keyword_set(fce_range) then fce_range = [1000.,6000.]  ;Hz

    ;Number of freq and fce array elements (resolution)
    if ~keyword_set(nfreq) then nfreq = 20.
    if ~keyword_set(nfce) then nfce = 20.
    if ~keyword_set(dens) then dens = 10.

    freq = (freq_range[1]-freq_range[0])*indgen(nfreq)/(nfreq-1) + freq_range[0]
    fce = (fce_range[1]-fce_range[0])*indgen(nfce)/(nfce-1) + fce_range[0]
  endif

  if scheme eq 2 then begin

    if ~keyword_set(freq_range) then begin
      if ~keyword_set(f_fce) then freq_range = [80.,300.]       ;Hz
      if keyword_set(f_fce) then  freq_range = [0.,2.]  ;unitless f/fce
    endif
    if ~keyword_set(density_range) then density_range = [1.,50.] ;cm-3

    ;Number of freq and fce array elements (resolution)
    if ~keyword_set(nfreq) then nfreq = 20.
    if ~keyword_set(ndens) then ndens = 20.
    if ~keyword_set(fce) then fce = 2000.

    freq = (freq_range[1]-freq_range[0])*indgen(nfreq)/(nfreq-1) + freq_range[0]
    ;if keyword_set(f_fce) then freq *= fce
    dens = (density_range[1]-density_range[0])*indgen(ndens)/(ndens-1) + density_range[0]

  endif


  Etots = fltarr(nfreq,ndens)
  Ez = Etots

  fpe = 8980d*sqrt(dens)
  c=3d5


  if scheme eq 0 then begin
    for i=0L,n_elements(dens)-1 do begin
      for j=0L,n_elements(fce)-1 do begin

        kvec = sqrt(4*!pi^2*fpe[i]^2*freq/(c^2*(fce[j]*cos(theta_k*!dtor)-freq)))
        ;kvec = 60.
        evals = cycl_energies(freq,theta_k,pa,fce[j],kvec,dens,nres)

        if type eq 'counterstream' then begin
          Ez[i,j] = evals.ez_cycl_counterstream
          Etots[i,j] = evals.e_cycl_counterstream
        endif
        if type eq 'costream' then begin
          Ez[i,j] = evals.ez_cycl_costream
          Etots[i,j] = evals.e_cycl_costream
        endif
        if type eq 'landau' then begin
          Ez[i,j] = evals.ez_landau
          Etots[i,j] = evals.ez_landau
        endif
      endfor
    endfor
  endif


  if scheme eq 1 then begin
    for i=0L,n_elements(freq)-1 do begin
      for j=0L,n_elements(fce)-1 do begin

        kvec = sqrt(4.*!pi^2*fpe^2*freq[i]/(c^2*(fce[j]*cos(theta_k*!dtor)-freq[i])))
        ;kvec = 60.
        evals = cycl_energies(freq[i],theta_k,pa,fce[j],kvec,dens,nres)
        if type eq 'counterstream' then begin
          Ez[i,j] = evals.ez_cycl_counterstream
          Etots[i,j] = evals.e_cycl_counterstream
        endif
        if type eq 'costream' then begin
          Ez[i,j] = evals.ez_cycl_costream
          Etots[i,j] = evals.e_cycl_costream
        endif
        if type eq 'landau' then begin
          Ez[i,j] = evals.ez_landau
          Etots[i,j] = evals.ez_landau
        endif

      endfor
    endfor
  endif

  if scheme eq 2 then begin
    for i=0L,n_elements(freq)-1 do begin
      for j=0L,n_elements(dens)-1 do begin

        kvec = sqrt(4.*!pi^2*fpe[j]^2*freq[i]/(c^2*(fce*cos(theta_k*!dtor)-freq[i])))
        ;kvec = 60.
        evals = cycl_energies(freq[i],theta_k,pa,fce,kvec,dens,nres)
        if type eq 'counterstream' then begin
          Ez[i,j] = evals.ez_cycl_counterstream
          Etots[i,j] = evals.e_cycl_counterstream
        endif
        if type eq 'costream' then begin
          Ez[i,j] = evals.ez_cycl_costream
          Etots[i,j] = evals.e_cycl_costream
        endif
        if type eq 'landau' then begin
          Ez[i,j] = evals.ez_landau
          Etots[i,j] = evals.ez_landau
        endif
      endfor
    endfor
  endif


  ;---------------------------------------------------------
  ;plot a spectra of the cyclotron
  ;energies in parameter space
  ;---------------------------------------------------------


  pastr = strtrim(string(pa,format='(F5.1)'),2)
  theta_kstr = strtrim(string(theta_k,format='(F5.1)'),2)
  if scheme eq 0 then freqstr = strtrim(string(freq/1000.,format='(F9.3)'),2)
  if scheme eq 1 then densstr = strtrim(string(dens,format='(F7.1)'),2)
  if scheme eq 2 then fcestr = strtrim(string(fce/1000.,format='(F7.1)'),2)



  if scheme eq 0 then titlestr = 'Etotal '+type+' energy (keV)!Cfor e- w/ PA='+pastr+$
  ' deg!CFor theta_kb='+theta_kstr+$
  ' deg!Cfor whistler mode wave of '+freqstr+' kHz'+$
  '!Cfrom cycl_energies_parameter_space.pro'
  if scheme eq 1 then titlestr = 'Etotal '+type+' energy (keV)!Cfor e- w/ PA='+pastr+$
  ' deg!CFor theta_kb='+theta_kstr+$
  ' deg!Cfor density of '+densstr+' cm-3'+$
  '!Cfrom cycl_energies_parameter_space.pro'
  if scheme eq 2 then titlestr = 'Etotal '+type+' energy (keV)!Cfor e- w/ PA='+pastr+$
  ' deg!CFor theta_kb='+theta_kstr+$
  ' deg!Cfor fce of '+fcestr+' kHz'+$
  '!Cfrom cycl_energies_parameter_space.pro'


  titlestr2 = 'Ez '+type+' energy (keV)'

  if scheme eq 0 then typestr = 'fce_vs_dens--freq='+freqstr+'kHz--PA='+pastr+'deg--TBK='+theta_kstr+'deg'
  if scheme eq 1 then typestr = 'fce_vs_freq--dens='+densstr+'cm3--PA='+pastr+'deg--TBK='+theta_kstr+'deg'
  if scheme eq 2 and ~KEYWORD_SET(f_fce) then typestr = 'dens_vs_freq--fce='+fcestr+'kHz--PA='+pastr+'deg--TBK='+theta_kstr+'deg'
  if scheme eq 2 and  KEYWORD_SET(f_fce) then typestr = 'dens_vs_fdivfce--fce='+fcestr+'kHz--PA='+pastr+'deg--TBK='+theta_kstr+'deg'



  plottitle = type+'energy_spec--'+typestr
  if keyword_set(ps) then popen,'~/Desktop/'+plottitle+'.ps',/landscape
  if keyword_set(ps) then !p.charsize = 1.


  maxval_print = max(Etots)
  maxzval_print = max(Ez)

  ;Make values above and below max/min values equal to the max/min values.
  goo = where(Etots gt maxval)
  if goo[0] ne -1 then Etots[goo] = maxval
  goo = 0.
  goo = where(Etots lt minval)
  if goo[0] ne -1 then Etots[goo] = minval


  ;Define levels for the colors (from dfanning website:
  ;page 144: http://www.idlcoyote.com/books/tg/samples/tg_chap5.pdf
  ;Letting IDL manually define the colors based on the nlevels keyword
  ;often leads to bad color scales.
  nlevels = 12
  LoadCT, 33, NColors=nlevels, Bottom=1
  step = (maxval-minval) / nlevels
  levels = IndGen(nlevels) * step + minval   ;in keV

  stepz = (maxzval-minzval)/nlevels
  levelsz = IndGen(nlevels) * stepz + minzval

  SetDecomposedState, 0, CurrentState=currentState


  !p.multi = [0,0,2]

  if scheme eq 0 then begin
    SetDecomposedState, currentState
    contour,Etots,dens,fce/1000.,$
    xtitle='dens (cm-3)',ytitle='fce (kHz)',$
    /fill,$
    C_Colors=IndGen(nlevels)+1,$
    levels=levels,$
    title=titlestr,$
    ymargin=[4,8],xmargin=[10,20],$
    xrange=[min(dens),max(dens)],$
    yrange=[min(fce/1000.),max(fce/1000.)],$
    xstyle=1,ystyle=1

    SetDecomposedState, currentState
    contour,Ez,dens,fce/1000.,$
    xtitle='dens (cm-3)',ytitle='fce (kHz)',$
    /fill,$
    C_Colors=IndGen(nlevels)+1,$
    levels=levelsz,$
    title=titlestr2,$
    ymargin=[4,8],xmargin=[10,20],$
    xrange=[min(dens),max(dens)],$
    yrange=[min(fce/1000.),max(fce/1000.)],$
    xstyle=1,ystyle=1
  endif


  if scheme eq 1 then begin
    if ~KEYWORD_SET(f_fce) then xtitle = 'freq (kHz)' else xtitle = 'f/fce'
    if ~KEYWORD_SET(f_fce) then tmpy=freq/1000. else tmpy=freq/fce
    if ~KEYWORD_SET(f_fce) then xr_tmp=[min(freq/1000.),max(freq/1000.)] else $
      xr_tmp=[min(freq/fce),max(freq/fce)]

    SetDecomposedState, currentState
    contour,Etots,tmpy,fce/1000.,$
    xtitle=xtitle,ytitle='fce (kHz)',$
    /fill,$
    C_Color=IndGen(nlevels)+1,$
    levels=levels,$
    title=titlestr,$
    ymargin=[4,8],xmargin=[10,20],$
    xrange=xr_tmp,yrange=[min(fce/1000.),max(fce/1000.)],$
    xstyle=1,ystyle=1

    SetDecomposedState, currentState
    contour,Ez,tmpy,fce/1000.,$
    xtitle=xtitle,ytitle='fce (kHz)',$
    /fill,$
    C_Colors=IndGen(nlevels)+1,$
    levels=levelsz,$
    title=titlestr2,$
    ymargin=[4,8],xmargin=[10,20],$
    xrange=xr_tmp,yrange=[min(fce/1000.),max(fce/1000.)],$
    xstyle=1,ystyle=1
  endif

  if scheme eq 2 then begin

    if ~KEYWORD_SET(f_fce) then xtitle = 'freq (kHz)' else xtitle = 'f/fce'
    if ~KEYWORD_SET(f_fce) then tmpy=freq/1000. else tmpy=freq/fce
    if ~KEYWORD_SET(f_fce) then xr_tmp=[min(freq/1000.),max(freq/1000.)] else $
      xr_tmp=[min(freq/fce),max(freq/fce)]
    SetDecomposedState, currentState
    contour,Etots,tmpy,dens,$
    xtitle=xtitle,ytitle='dens (cm-3)',$
    /fill,$
    C_Colors=IndGen(nlevels)+1,$
    levels=levels,$
    title=titlestr,$
    ymargin=[4,8],xmargin=[10,20],$
    xrange=xr_tmp,yrange=[min(dens),max(dens)],$
    xstyle=1,ystyle=1

    SetDecomposedState, currentState
    contour,Ez,tmpy,dens,$
    xtitle=xtitle,ytitle='dens (cm-3)',$
    /fill,$
    C_Colors=IndGen(nlevels)+1,$
    levels=levelsz,$
    title=titlestr2,$
    ymargin=[4,8],xmargin=[10,20],$
    xrange=xr_tmp,yrange=[min(dens),max(dens)],$
    xstyle=1,ystyle=1
  endif


  SetDecomposedState, currentState

  cgColorbar, Range=[minval,maxval], $
    Divisions=nlevels, XTicklen=1, XMinor=0, $
    AnnotateColor='black', NColors=nlevels, Bottom=1, $
    Position=[0.94, 0.55, 0.97, 0.90],/vertical,$
    Charsize=0.75;,_extra={yminor:9,ytickformat:'(f9.3)'}


  SetDecomposedState, currentState
  cgColorbar, Range=[minzval,maxzval], $
    Divisions=nlevels, XTicklen=1, XMinor=0, $
    AnnotateColor='black', NColors=nlevels, Bottom=1, $
    Position=[0.94, 0.05, 0.97, 0.4], /vertical,$
    Charsize=0.75;,_extra={yminor:9,ytickformat:'(f9.3)'}



  print,etots
  print,'Etots_max = ', max(maxval_print)
  print,'Ez_max = ', max(maxzval_print)
  if keyword_set(ps) then pclose

end
