@natconst
;
; Grid parameters
;
nx       = 10L
ny       = 1L
nz       = 1L
sizex    = 1000*AU
sizey    = 1000*AU/nx
sizez    = 1000*AU/nx
;
; Model parameters
;
abun_h2  = 0.5
abun_he  = 0.1
mgas     = mp*(2*abun_h2+4*abun_he)/abun_h2                ; Mass of gas per H2-molecule
mugas    = mp*(2*abun_h2+4*abun_he)/(abun_h2+abun_he)      ; Mass of gas per particle (H2,He)
nh2      = 1d5
rhogas0  = nh2*mgas
temp0    = 30.
tdust0   = 30.
dusttogas= 1d-2
vturb0   = 1.*1d5
dvdau    = 1d-2*1d5              ; Velocity gradient in (km/s) per AU
;dvdau    = 1d-4*1d5              ; Velocity gradient in (km/s) per AU
;dvdau    = 1d-1*1d5              ; Velocity gradient in (km/s) per AU
;dvdau    = 0.d0              ; Velocity gradient in (km/s) per AU
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
xi      = -sizex + 2*sizex*dindgen(nx+1)/(1.d0*nx)
yi      = -sizey + 2*sizey*dindgen(ny+1)/(1.d0*ny)
zi      = -sizez + 2*sizez*dindgen(nz+1)/(1.d0*nz)
xc      = 0.5 * ( xi[0:nx-1] + xi[1:nx] )
yc      = 0.5 * ( yi[0:ny-1] + yi[1:ny] )
zc      = 0.5 * ( zi[0:nz-1] + zi[1:nz] )
;
; Model
;
rhogas  = dblarr(nx,ny,nz) + rhogas0
tgas    = dblarr(nx,ny,nz) + temp0
tdust   = dblarr(nx,ny,nz) + tdust0
vx      = rebin(dvdau*xc/au,nx,ny,nz)
vy      = dblarr(nx,ny,nz)
vz      = dblarr(nx,ny,nz)
vturb   = dblarr(nx,ny,nz) + vturb0
;
; Write the grid file
;
openw,1,'amr_grid.inp'
printf,1,1                      ; iformat
printf,1,0                      ; AMR grid style  (0=regular grid, no AMR)
printf,1,0                      ; Coordinate system
printf,1,0                      ; gridinfo
printf,1,1,1,1                  ; Include x,y,z coordinate
printf,1,nx,ny,nz               ; Size of grid
for i=0,nx do printf,1,xi[i]    ; X coordinates (cell walls)
for i=0,ny do printf,1,yi[i]    ; Y coordinates (cell walls)
for i=0,nz do printf,1,zi[i]    ; Z coordinates (cell walls)
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
abun = 1d-4
fact = abun/mgas
openw,1,'numberdens_co.inp'
printf,1,1                      ; Format number
printf,1,nx*ny*nz               ; Nr of cells
for iz=0,nz-1 do begin
   for iy=0,ny-1 do begin
      for ix=0,nx-1 do begin
         printf,1,rhogas[ix,iy,iz]*fact
      endfor
   endfor
endfor
close,1
;
; Write the number density file for the collision partner 
;
fact = 1.d0/mgas                ; Assume all H2 is para
openw,1,'numberdens_h2.inp'
printf,1,1                      ; Format number
printf,1,nx*ny*nz               ; Nr of cells
for iz=0,nz-1 do begin
   for iy=0,ny-1 do begin
      for ix=0,nx-1 do begin
         printf,1,rhogas[ix,iy,iz]*fact
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
printf,1,'co    leiden    0    0    1'
printf,1,'h2'
close,1
;
; Write the radmc3d.inp control file
;
openw,1,'radmc3d.inp'
;printf,1,'nphot = 1000000'
;printf,1,'scattering_mode_max = 0'
;printf,1,'tgas_eq_tdust   = 1'
printf,1,'lines_mode = 3'   ;; LVG mode
close,1
;
end
