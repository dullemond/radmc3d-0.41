;
; Some natural constants
;
AU  = 1.49598d13     ; Astronomical Unit       [cm]
pc  = 3.08572d18     ; Parsec                  [cm]
MS  = 1.98892d33     ; Solar mass              [g]
TS  = 5.78d3         ; Solar temperature       [K]
LS  = 3.8525d33      ; Solar luminosity        [erg/s]
RS  = 6.96d10        ; Solar radius            [cm]
;
; Monte Carlo parameters
;
nphot    = 1000000
;
; Grid parameters
;
nx       = 32
ny       = 32
nz       = 32
sizex    = 10*AU
sizey    = 10*AU
sizez    = 10*AU
;
; Model parameters
;
radius   = 5*AU
rho0     = 1d-16
temp0    = 20.d0                ; Put this to 0 if you dont want temp to be set
tnoise   = 0.d0                 ; Noise in temperature = 0
;tnoise   = 10.d0               ; To clearly see the AMR refinement, make noise
;
; Star parameters
;
mstar    = MS
rstar    = RS
tstar    = TS
pstar    = [0.,0.,0.]
;
; Write the wavelength_micron.inp file
;
lambda1 = 0.1d0
lambda2 = 7.0d0
lambda3 = 25.d0
lambda4 = 1.0d4
n12     = 20
n23     = 100
n34     = 30
lam12   = lambda1 * (lambda2/lambda1)^(dindgen(n12)/(1.d0*n12))
lam23   = lambda2 * (lambda3/lambda2)^(dindgen(n23)/(1.d0*n23))
lam34   = lambda3 * (lambda4/lambda3)^(dindgen(n34)/(1.d0*(n34-1.d0)))
lambda  = [lam12,lam23,lam34]
nlam    = n_elements(lambda)
;
; Write the wavelength file
;
openw,1,'wavelength_micron.inp'
printf,1,nlam
for ilam=0,nlam-1 do printf,1,lambda[ilam]
close,1
;
; Write the stars.inp file
;
openw,1,'stars.inp'
printf,1,2
printf,1,1,nlam
printf,1,' '
printf,1,rstar,mstar,pstar[0],pstar[1],pstar[2]
printf,1,' '
for ilam=0,nlam-1 do printf,1,lambda[ilam]
printf,1,' '
printf,1,-tstar
close,1
;
; Dust opacity control file
;
openw,1,'dustopac.inp'
printf,1,'2               Format number of this file'
printf,1,'1               Nr of dust species'
printf,1,'============================================================================'
printf,1,'1               Way in which this dust species is read'
printf,1,'0               0=Thermal grain'
printf,1,'silicate        Extension of name of dustkappa_***.inp file'
printf,1,'----------------------------------------------------------------------------'
close,1
;
; Write the radmc3d.inp control file
;
openw,1,'radmc3d.inp'
printf,1,'nphot = ',nphot
printf,1,'scattering_mode_max = 0'
printf,1,'iranfreqmode = 1'
printf,1,'userdef_nx = ',nx
printf,1,'userdef_ny = ',ny
printf,1,'userdef_nz = ',nz
printf,1,'userdef_sizex = ',sizex
printf,1,'userdef_sizey = ',sizey
printf,1,'userdef_sizez = ',sizez
printf,1,'userdef_radius =',radius
printf,1,'userdef_rho0 = ',rho0
printf,1,'userdef_levelmax = ',10
printf,1,'userdef_nrefinefact = ',4.0
printf,1,'userdef_amr_relcellsize = ',0.05
printf,1,'userdef_amr_refregion = ',0.5
printf,1,'userdef_temp0 = ',temp0
printf,1,'userdef_tempnoise = ',tnoise
close,1
;
end
