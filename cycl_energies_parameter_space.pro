;+
; NAME: cycl_energies_parameter_space
; SYNTAX:
; PURPOSE: create a parameter space plot of relativistic cyclotron resonance energies as a function
; of wave freq, fce, density
; INPUT: (see cycl_energies_parameter_space_call.pro. Tested against Lorentzen01, plate7 results)
; OUTPUT:
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
;          fcev          -> value of input fce
;          freqv         -> value of input freq
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
  type=type,harmonic=nres,f_fce=f_fce


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
  if ~keyword_set(minval) then minval = 0.1                ;keV
  if ~keyword_set(maxzval) then maxzval = 500.                ;keV
  if ~keyword_set(minzval) then minzval = 0.1                ;keV


  zlog = 1


  ;Choose range of parameters
  if scheme eq 0 then begin
    if ~keyword_set(density_range) then density_range=[1.,50.]      ;cm-3
    if ~keyword_set(fce_range) then fce_range = [1000.,6000.]  ;Hz

    ;Number of density and fce array elements (resolution)
    if ~keyword_set(ndens) then ndens = 20.
    if ~keyword_set(nfce) then nfce = 20.
    if ~keyword_set(freq) then freq = 100.

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
      if keyword_set(f_fce) then  freq_range = [0,2]  ;unitless f/fce
    endif
    if ~keyword_set(density_range) then density_range = [1.,50.] ;cm-3

    ;Number of freq and fce array elements (resolution)
    if ~keyword_set(nfreq) then nfreq = 20.
    if ~keyword_set(ndens) then ndens = 20.
    if ~keyword_set(fce) then fce = 2000.

    freq = (freq_range[1]-freq_range[0])*indgen(nfreq)/(nfreq-1) + freq_range[0]
    if keyword_set(f_fce) then freq *= fce
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

        kvec = sqrt(4*!pi^2*fpe^2*freq[i]/(c^2*(fce[j]*cos(theta_k*!dtor)-freq[i])))
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

        kvec = sqrt(4*!pi^2*fpe[j]^2*freq[i]/(c^2*(fce*cos(theta_k*!dtor)-freq[i])))
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

  fpe /= 1000.    ;put into kHz
  if ~keyword_set(f_fce) then freq /= 1000.  ;kHz
  if KEYWORD_SET(f_fce) then freq /= fce  ;f/fce
  fce /= 1000.

  ;---------------------------------------------------------
  ;plot a spectra of the cyclotron
  ;energies in parameter space
  ;---------------------------------------------------------


  if zlog then EtotsBS = bytscl(alog10(Etots),min=alog10(minval),max=alog10(maxval)) $
  else         EtotsBS = bytscl(Etots,min=minval,max=maxval)
  if zlog then EzBS = bytscl(alog10(Ez),min=alog10(minzval),max=alog10(maxzval)) $
  else         EzBS = bytscl(Ez,min=minzval,max=maxzval)


  loadct,39

  pastr = strtrim(string(pa,format='(F5.1)'),2)
  theta_kstr = strtrim(string(theta_k,format='(F5.1)'),2)
  if scheme eq 0 then freqstr = strtrim(string(freq,format='(F7.1)'),2)
  if scheme eq 1 then densstr = strtrim(string(dens,format='(F7.1)'),2)
  if scheme eq 2 then fcestr = strtrim(string(fce,format='(F7.1)'),2)



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
  if scheme eq 2 and  KEYWORD_SET(f_fce) then typestr = 'dens_vs_f/fce--fce='+fcestr+'kHz--PA='+pastr+'deg--TBK='+theta_kstr+'deg'



  plottitle = type+'energy_spec--'+typestr

  if keyword_set(ps) then popen,'~/Desktop/'+plottitle+'.ps',/landscape
  if keyword_set(ps) then !p.charsize = 1.

  !p.multi = [0,0,2]


  if scheme eq 0 then begin
    contour,EtotsBS,dens,fce,$
    xtitle='dens (cm-3)',ytitle='fce (kHz)',$
    /fill,nlevels=40,$
    title=titlestr,$
    ymargin=[4,8],xmargin=[10,20],$
    xrange=[min(dens),max(dens)],yrange=[min(fce),max(fce)],$
    xstyle=1,ystyle=1

    contour,EzBS,dens,fce,$
    xtitle='dens (cm-3)',ytitle='fce (kHz)',$
    /fill,nlevels=40,$
    title=titlestr2,$
    ymargin=[4,8],xmargin=[10,20],$
    xrange=[min(dens),max(dens)],yrange=[min(fce),max(fce)],$
    xstyle=1,ystyle=1
  endif


  if scheme eq 1 then begin
    contour,EtotsBS,freq,fce,$
    xtitle='freq (kHz)',ytitle='fce (kHz)',$
    /fill,nlevels=40,$
    title=titlestr,$
    ymargin=[4,8],xmargin=[10,20],$
    xrange=[min(freq),max(freq)],yrange=[min(fce),max(fce)],$
    xstyle=1,ystyle=1

    contour,EzBS,freq,fce,$
    xtitle='freq (kHz)',ytitle='fce (kHz)',$
    /fill,nlevels=40,$
    title=titlestr2,$
    ymargin=[4,8],xmargin=[10,20],$
    xrange=[min(freq),max(freq)],yrange=[min(fce),max(fce)],$
    xstyle=1,ystyle=1
  endif

  if scheme eq 2 then begin
    if ~KEYWORD_SET(f_fce) then xtitle = 'freq (kHz)' else xtitle = 'f/fce'

    contour,EtotsBS,freq,dens,$
    xtitle=xtitle,ytitle='dens (cm-3)',$
    /fill,nlevels=40,$
    title=titlestr,$
    ymargin=[4,8],xmargin=[10,20],$
    xrange=[min(freq),max(freq)],yrange=[min(dens),max(dens)],$
    xstyle=1,ystyle=1

    contour,EzBS,freq,dens,$
    xtitle=xtitle,ytitle='dens (cm-3)',$
    /fill,nlevels=40,$
    title=titlestr2,$
    ymargin=[4,8],xmargin=[10,20],$
    xrange=[min(freq),max(freq)],yrange=[min(dens),max(dens)],$
    xstyle=1,ystyle=1
  endif

  ;.compile ~/Desktop/code/Aaron/RBSP/coyote/colorbar.pro
  colorbar,range=[minval,maxval],position=[0.94, 0.55, 0.97, 0.90],/vertical,/ylog,$
  _extra={yminor:9,ytickformat:'(f9.3)'}

  colorbar,range=[minzval,maxzval],position=[0.94, 0.05, 0.97, 0.4],/vertical,/ylog,$
  _extra={yminor:9,ytickformat:'(f9.3)'}



  if keyword_set(ps) then pclose

end
