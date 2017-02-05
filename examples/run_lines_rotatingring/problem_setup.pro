@natconst
;
; Grid parameters
;
nx       = 4L      ;; 1 cell in radius
ny       = 4L      ;; 2 cells in theta
;nz       = 1L      ;; Uncomment to try with 1 azimuthal cell
nz       = 360L    ;; Uncomment to try with 360 azimuthal cells
ri       = AU      ;; Inner edge of cell
ro       = 1.1*AU  ;; Outer edge of cell
hr       = 0.1     ;; Vertical thickness (H/R) of cell
;
; Model parameters
;
rhogas0  = 3d-19
temp0g   = 300.d0
temp0d   = 300.d0 
dusttogas= 0.01
vturb0   = 0.1*1d5
mstar    = MS
rc       = 0.5*(ri+ro)
;
; Star parameters
; NOTE: Star is commented out in this setup (if you test it out, and afterward
; want to make sure you switch the star off again, make sure to delete the
; stars.inp file).
;
;mstar    = MS
;rstar    = RS
;tstar    = TS
;pstar    = [0.,0.,0.]
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
; Make grid with size interpreted as half-width size here
;
xi      = ri + (ro-ri)*dindgen(nx+1)/(1.d0*nx)
yi      = !dpi/2.d0 -hr + 2*hr*dindgen(ny+1)/(1.d0*ny)
zi      = 0.d0 + 2*!dpi*dindgen(nz+1)/(1.d0*nz)
xc      = 0.5 * ( xi[0:nx-1] + xi[1:nx] )
yc      = 0.5 * ( yi[0:ny-1] + yi[1:ny] )
zc      = 0.5 * ( zi[0:nz-1] + zi[1:nz] )
;
; Kepler
;
vk       = sqrt(GG*mstar/xc)
;
; Model
;
rhogas  = dblarr(nx,ny,nz) + rhogas0
tgas    = dblarr(nx,ny,nz) + temp0g
tdust   = dblarr(nx,ny,nz) + temp0d
vx      = dblarr(nx,ny,nz)
vy      = dblarr(nx,ny,nz)
vz      = dblarr(nx,ny,nz) + rebin(vk,nx,ny,nz)
vturb   = dblarr(nx,ny,nz) + vturb0
;
; Write the stars.inp file
;
;openw,1,'stars.inp'
;printf,1,2
;printf,1,1,nlam
;printf,1,' '
;printf,1,rstar,mstar,pstar[0],pstar[1],pstar[2]
;printf,1,' '
;for ilam=0,nlam-1 do printf,1,lambda[ilam]
;printf,1,' '
;printf,1,-tstar
;close,1
;
; Write the grid file
;
openw,1,'amr_grid.inp'
printf,1,1                      ; iformat
printf,1,0                      ; AMR grid style  (0=regular grid, no AMR)
printf,1,100                    ; Coordinate system
printf,1,0                      ; gridinfo
if nz eq 1 then begin
   printf,1,1,1,0                              ; Include x,y coordinate
endif else begin
   printf,1,1,1,1                              ; Include x,y,z coordinate
endelse
printf,1,nx,ny,nz               ; Size of grid
for i=0,nx do printf,1,xi[i],format='(E19.12)'    ; X coordinates (cell walls)
for i=0,ny do printf,1,yi[i],format='(E24.17)'    ; Y coordinates (cell walls)
for i=0,nz do printf,1,zi[i],format='(E19.12)'    ; Z coordinates (cell walls)
close,1
;
; Write the dust density file. Here we use a dust-to-gas ratio of 0.01
;
openw,1,'dust_density.inp'
printf,1,1                      ; Format number
printf,1,nx*ny*nz               ; Nr of cells
printf,1,1                      ; Nr of species
for iz=0,nz-1 do begin
   for iy=0,ny-1 do begin
      for ix=0,nx-1 do begin
         printf,1,rhogas[ix,iy,iz]*dusttogas
      endfor
   endfor
endfor
close,1
;
; Write the molecule number density file. 
;
abunco = 1d-4
factco = abunco/(2.3*mp)
openw,1,'numberdens_co.inp'
printf,1,1                      ; Format number
printf,1,nx*ny*nz               ; Nr of cells
for iz=0,nz-1 do begin
   for iy=0,ny-1 do begin
      for ix=0,nx-1 do begin
         printf,1,rhogas[ix,iy,iz]*factco
      endfor
   endfor
endfor
close,1
;
; Write the gas velocity field
;
openw,1,'gas_velocity.inp'
printf,1,1                      ; Format number
printf,1,nx*ny*nz               ; Nr of cells
for iz=0,nz-1 do begin
   for iy=0,ny-1 do begin
      for ix=0,nx-1 do begin
         printf,1,vx[ix,iy,iz],vy[ix,iy,iz],vz[ix,iy,iz] 
      endfor
   endfor
endfor
close,1
;
; Write the microturbulence file
;
openw,1,'microturbulence.inp'
printf,1,1                      ; Format number
printf,1,nx*ny*nz               ; Nr of cells
for iz=0,nz-1 do begin
   for iy=0,ny-1 do begin
      for ix=0,nx-1 do begin
         printf,1,vturb[ix,iy,iz]
      endfor
   endfor
endfor
close,1
;
; Write the gas temperature file
;
openw,1,'gas_temperature.inp'
printf,1,1                      ; Format number
printf,1,nx*ny*nz               ; Nr of cells
for iz=0,nz-1 do begin
   for iy=0,ny-1 do begin
      for ix=0,nx-1 do begin
         printf,1,tgas[ix,iy,iz]
      endfor
   endfor
endfor
close,1
;
; Write the dust temperature file
;
openw,1,'dust_temperature.dat'
printf,1,1                      ; Format number
printf,1,nx*ny*nz               ; Nr of cells
printf,1,1                      ; Nr of dust species
for iz=0,nz-1 do begin
   for iy=0,ny-1 do begin
      for ix=0,nx-1 do begin
         printf,1,tdust[ix,iy,iz]
      endfor
   endfor
endfor
close,1
;
; Write the wavelength file
;
openw,1,'wavelength_micron.inp'
printf,1,nlam
for ilam=0,nlam-1 do printf,1,lambda[ilam]
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
; Write the lines.inp control file
;
openw,1,'lines.inp'
printf,1,'2'
printf,1,'1'
printf,1,'co    leiden    0    0    0'
close,1
;
; Write the radmc3d.inp control file
;
openw,1,'radmc3d.inp'
;printf,1,'nphot = 1000000'
;printf,1,'scattering_mode_max = 0'
;printf,1,'tgas_eq_tdust   = 1'
printf,1,'lines_mode = 1'   ;; LTE
close,1
;
end
