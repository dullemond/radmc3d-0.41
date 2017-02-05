@flash2radmc.pro 
;;
;; Some constants
;;
GG  = 6.672d-8                  ; Gravitational constant
mp  = 1.6726d-24                ; Mass of proton          [g]
kk  = 1.3807d-16                ; Bolzmann's constant     [erg/K]
cc  = 2.9979245800000d10        ; Light speed             [cm/s]
ss  = 5.6703d-5                 ; Stefan-Boltzmann const  [erg/cm^2/K^4/s]
LS  = 3.8525d33                 ; Solar luminosity        [erg/s]
RS  = 6.96d10                   ; Solar radius            [cm]
MS  = 1.99d33                   ; Solar mass              [g]
TS  = 5.78d3                    ; Solar temperature       [K]
AU  = 1.496d13                  ; Astronomical Unit       [cm]
pc  = 3.08572d18                ; Parsec                  [cm]
;;
;; Put the FLASH model onto the RADMC-3D grid
;;
filename='BE_collapse_sinks_hdf5_plt_cnt_0050'
flash2radmc,filename,/writetemp,/center
;;
;; Put a star somewhere
;;
mstar    = 20*MS
rstar    = 5*RS
tstar    = 20000.
;;
;; Make a wavelength grid which has some refinement in the dust feature region
;;
lambda1 = 0.05d0
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
;;
;; Write the wavelength file
;;
print,'Writing wavelength_micron.inp...'
openw,1,'wavelength_micron.inp'
printf,1,nlam
for ilam=0,nlam-1 do printf,1,lambda[ilam]
close,1
;;
;; Write the stars.inp file
;;
nstars=1
posstar=[0.0859,0.6640,-0.1328]*1d16
;;
print,'Writing stars.inp...'
openw,1,'stars.inp'
printf,1,2
printf,1,nstars,nlam
printf,1,' '
printf,1,rstar,mstar,posstar[0],posstar[1],posstar[2]
printf,1,' '
for ilam=0,nlam-1 do printf,1,lambda[ilam]
printf,1,' '
printf,1,-tstar
close,1
;;
print,'Writing dustopac.inp'
openw,1,'dustopac.inp'
printf,1,'2     Format number'
printf,1,'1     Number of dust species'
printf,1,'----------------------------'
printf,1,'1     Way in which dust species is read'
printf,1,'0     Means thermal grain'
printf,1,'1     Extension of dustkappa file'
printf,1,'----------------------------'
close,1

end

