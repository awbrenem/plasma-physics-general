;Test relativistic cyclotron resonance energy calculations

rbsp_efw_init

;--------------------------------------------------------
;Create Plate 7a plot from Lorentzen01 that shows Ecycl as a function
;of freq for three latitudes (0, 15, 30)
;--------------------------------------------------------


c_ms = 2.99792458d8      ; -Speed of light in vacuum (m/s)
c_kms = c_ms/1000.



;get fce values from dipole model
.compile /Users/aaronbreneman/Desktop/code/Aaron/IDL/analysis/dipole.pro
L = dipole(4.9)

plot,l.colat,l.b

fce0 = 28*l.b[0]
goo = where(l.colat ge 15.)
fce15 = 28.*l.b[goo[0]]
goo = where(l.colat ge 30.)
fce30 = 28.*l.b[goo[0]]

freq = indgen(1000)*8000./999.
pitchangle = replicate(5.,1000.)
theta_kb = replicate(30.,1000.)
nres = replicate(1.,1000.)


;determine k from whistler dispersion relation for low freq whistlers
fce = replicate(fce0,1000.)
dens = replicate(5.,1000.)
fpe = 8980.*sqrt(dens)
wpe = 2.*!pi*fpe
w = 2.*!pi*freq
wce = 2.*!pi*fce
index_ref2 = 1 + wpe^2/(w*(wce*cos(theta_kb*!dtor)-w))
kmag = sqrt(index_ref2*w^2/c_kms^2)
vals_0lat = cycl_energies(freq,theta_kb,pitchangle,fce,kmag,nres)

fce = replicate(fce15,1000.)
dens = replicate(5.,1000.)
wce = 2.*!pi*fce
index_ref2 = 1 + wpe^2/(w*(wce*cos(theta_kb*!dtor)-w))
kmag = sqrt(index_ref2*w^2/c_kms^2)
vals_15lat = cycl_energies(freq,theta_kb,pitchangle,fce,kmag,nres)

fce = replicate(fce30,1000.)
dens = replicate(5.,1000.)
wce = 2.*!pi*fce
index_ref2 = 1 + wpe^2/(w*(wce*cos(theta_kb*!dtor)-w))
kmag = sqrt(index_ref2*w^2/c_kms^2)
vals_30lat = cycl_energies(freq,theta_kb,pitchangle,fce,kmag,nres)


;Plot the values to recreate Lorentzen01 Plate 7a
mevline = replicate(1d6,1000)

!p.charsize = 2
plot,freq,vals_0lat.E_cycl_normal*1000.,/ylog,yrange=[1d3,1d7],xtitle='freq(Hz)',ytitle='Energy(eV)'
oplot,freq,vals_15lat.E_cycl_normal*1000
oplot,freq,vals_30lat.E_cycl_normal*1000
oplot,freq,mevline,linestyle=2



;--------------------------------------------------------
;Create Plate 7b plot from Lorentzen01 that shows Ecycl as a function
;of freq for three harmonics at lat = 0
;--------------------------------------------------------

freq = indgen(1000)*8000./999.
pitchangle = replicate(5.,1000.)
theta_kb = replicate(30.,1000.)

fce = replicate(fce0,1000.)
dens = replicate(5.,1000.)
fpe = 8980.*sqrt(dens)
wpe = 2.*!pi*fpe
w = 2.*!pi*freq
wce = 2.*!pi*fce
index_ref2 = 1 + wpe^2/(w*(wce*cos(theta_kb*!dtor)-w))
kmag = sqrt(index_ref2*w^2/c_kms^2)
nres = replicate(1.,1000.)
vals_s1 = cycl_energies(freq,theta_kb,pitchangle,fce,kmag,nres)
nres = replicate(2.,1000.)
vals_s2 = cycl_energies(freq,theta_kb,pitchangle,fce,kmag,nres)
nres = replicate(3.,1000.)
vals_s3 = cycl_energies(freq,theta_kb,pitchangle,fce,kmag,nres)



!p.charsize = 2
plot,freq,vals_s1.E_cycl_normal*1000.,/ylog,yrange=[1d3,1d7],xtitle='freq(Hz)',ytitle='Energy(eV)'
oplot,freq,vals_s2.E_cycl_normal*1000
oplot,freq,vals_s3.E_cycl_normal*1000
oplot,freq,mevline,linestyle=2















;Attempt to recreate Abel and Thorne98 Fig2

L = dipole(3.6)
fce0 = 28*l.b[0]

freq = indgen(10000)*8000./9999.
pitchangle = replicate(5.,10000.)
theta_kb = replicate(0.,10000.)
fce = replicate(fce0,10000.)


;determine k from whistler dispersion relation for low freq whistlers
dens = replicate(200.,10000.)
fpe = 8980.*sqrt(dens)
wpe = 2.*!pi*fpe
w = 2.*!pi*freq
wce = 2.*!pi*fce
index_ref2 = 1 + wpe^2/(w*(wce*cos(theta_kb*!dtor)-w))
kmag = sqrt(index_ref2*w^2/c_kms^2)

nres = replicate(1.,10000.)
vals_1 = cycl_energies(freq,theta_kb,pitchangle,fce,kmag,nres)
nres = replicate(2.,10000.)
vals_2 = cycl_energies(freq,theta_kb,pitchangle,fce,kmag,nres)
nres = replicate(4.,10000.)
vals_4 = cycl_energies(freq,theta_kb,pitchangle,fce,kmag,nres)
nres = replicate(10.,10000.)
vals_10 = cycl_energies(freq,theta_kb,pitchangle,fce,kmag,nres)



!p.charsize = 2
plot,freq,vals_1.E_cycl,/ylog,yrange=[1d0,1d4],xrange=[4000,5000]
oplot,freq,vals_2.E_cycl
oplot,freq,vals_4.E_cycl
oplot,freq,vals_10.E_cycl
