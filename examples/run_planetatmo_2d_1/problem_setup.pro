;============================================================================
;            EXTREMELY SIMPLISTIC DUSTY PLANETARY ATMOSPHERE MODEL
;                            VERSION 13.01.10
;============================================================================
;
; Some natural constants
;
GG  = 6.673d-8       ; Gravitational constant
kk  = 1.3807d-16     ; Bolzmann's constant     [erg/K]
mp  = 1.6726d-24     ; Mass of proton          [g]
rearth   = 6366.d0 * 1d5 
tearth   = 255.            ;; http://classes.geology.uiuc.edu/055prgClass/geo116/8-1.pdf
mearth   = 5.97d27         ;; http://nssdc.gsfc.nasa.gov/planetary/planetfact.html
LS  = 3.8525d33      ; Solar luminosity        [erg/s]
RS  = 6.96d10        ; Solar radius            [cm]
MS  = 1.98892d33     ; Solar mass              [g]
TS  = 5.78d3         ; Solar temperature       [K]
AU  = 1.49598d13       ; Astronomical Unit       [cm]
pc  = 3.08572d18     ; Parsec                  [cm]
;
; Monte Carlo parameters
;
nphot    = 1000000L
;
; Grid parameters
;
nx       = 200L
ny       = 60L
nz       = 1L
;
; Model parameters
;
rin      = rearth
rout     = rearth+300*1d5
p0       = 1d6               ;; == 1 bar
t0       = 293.d0            ;; == 20 C
mugas    = 28.97d0           ;; http://www.engineeringtoolbox.com/molecular-mass-air-d_679.html
n0       = p0 / (kk*t0)      ;; p = n k T
rho0     = n0 * mugas * mp  
hp       = 20*1d5
;
; Planet parameters
;
mplan    = mearth
rplan    = rearth-1d3        ;; Sphere of Earth a fraction inside of grid
tplan    = t0
pplan    = [0.,0.,0.]
;
; Star parameters
;
mstar    = 0.1*MS
rstar    = RS
tstar    = 4000.
pstar    = [0.,0.,0.2*AU]    ;; Put star at 0.2 AU distance from planet
;
; Make the coordinates
;
xi       = rin * (rout/rin)^(dindgen(nx+1)/(nx*1.d0))
yi       = !dpi*dindgen(ny+1)/(ny*1.d0)
zi       = [0.,0.]
xc       = 0.5d0 * ( xi[0:nx-1] + xi[1:nx] )
yc       = 0.5d0 * ( yi[0:ny-1] + yi[1:ny] )
;
; Make the dust density model
;
rr       = xc
rhog     = rho0 * exp(-(rr-rin)/hp)
rhog     = rebin(rhog,nx,ny)
rhod     = 1d-5 * rhog
;
; Make the wavelength grid
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
printf,1,2,nlam
printf,1,' '
printf,1,rplan,mplan,pplan[0],pplan[1],pplan[2]
printf,1,rstar,mstar,pstar[0],pstar[1],pstar[2]
printf,1,' '
for ilam=0,nlam-1 do printf,1,lambda[ilam]
printf,1,' '
printf,1,-tplan
printf,1,-tstar
close,1
;
; Write the grid file
;
openw,1,'amr_grid.inp'
printf,1,1                      ; iformat
printf,1,0                      ; AMR grid style  (0=regular grid, no AMR)
printf,1,100                    ; Coordinate system
printf,1,0                      ; gridinfo
printf,1,1,1,0                  ; Include x,y,z coordinate
printf,1,nx,ny,1                 ; Size of grid
for i=0,nx do printf,1,xi[i]    ; X coordinates (cell walls)
for i=0,ny do printf,1,yi[i]    ; Y coordinates (cell walls)
for i=0,nz do printf,1,zi[i]    ; Z coordinates (cell walls)
close,1
;
; Write the density file
;
openw,1,'dust_density.inp'
printf,1,1                      ; Format number
printf,1,nx*ny*nz               ; Nr of cells
printf,1,1                      ; Nr of dust species
for iz=0,nz-1 do begin
   for iy=0,ny-1 do begin
      for ix=0,nx-1 do begin
         printf,1,rhod[ix,iy,iz]
      endfor
   endfor
endfor
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
; Write the radmc.inp control file
;
openw,1,'radmc.inp'
printf,1,'nphot = ',nphot          ;; Nr of photon packages for thermal MC
printf,1,'scattering_mode_max = 0' ;; For simplicity switch off scattering for now
printf,1,'istar_sphere = 1'        ;; Signal that "stars" to be treated as spheres
close,1
;
end
