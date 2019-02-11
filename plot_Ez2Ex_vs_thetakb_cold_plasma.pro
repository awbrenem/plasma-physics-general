;Plot variation of Ez/Ex as a function of theta_kb for various
;plasma conditions




freq = 100  ;Hz
dens = 7.  ;cm-3
bo = 8.   ;nT



freq = 1200  ;Hz
dens = 10.  ;cm-3
bo = 200.   ;nT


epol = 1+indgen(50.)/5.
thetakb = fltarr(50.)
erat = fltarr(50)

for i=0,n_elements(erat)-1 do begin $
  vals = cold_plasma_dispersion(epol=epol[i], freq=freq,dens=dens,Bo=bo) & $
  print,i  & $
  thetakb[i] = vals.theta_kb & $
  erat[i] = vals.Ez2Ex



print,thetakb
print,erat

plot,thetakb,erat,title='Ez/Ex vs wave normal angle (deg)',yrange=[0,1]
oplot,[vals.resangle,vals.resangle],[0,10],linestyle=2

oplot,thetakb,erat,color=250
oplot,[vals.resangle,vals.resangle],[0,10],linestyle=2,color=250
