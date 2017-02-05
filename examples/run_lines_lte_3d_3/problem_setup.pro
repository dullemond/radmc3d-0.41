@natconst
;
; Grid parameters
;
nx       = 16L
ny       = 16L
nz       = 16L
sizex    = 100*AU
sizey    = 100*AU
sizez    = 100*AU
;
; Model parameters
;
dvdy     = 100*1d5/(100*AU)
rhogas0  = 1d-16
;rhogas0  = 1d-66
temp0    = 50.d0
dusttogas= 0.01
vturb0   = 00.1*1d5
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
; Make grid with size interpreted as half-width size here
;
xi      = -sizex + 2*sizex*dindgen(nx+1)/(1.d0*nx)
yi      = -sizey + 2*sizey*dindgen(ny+1)/(1.d0*ny)
zi      = -sizez + 2*sizez*dindgen(nz+1)/(1.d0*nz)
xc      = 0.5 * ( xi[0:nx-1] + xi[1:nx] )
yc      = 0.5 * ( yi[0:ny-1] + yi[1:ny] )
zc      = 0.5 * ( zi[0:nz-1] + zi[1:nz] )
xx      = rebin(xc,nx,ny,nz)
yy      = transpose(rebin(yc,ny,nx,nz),[1,0,2])
zz      = transpose(rebin(zc,nz,ny,nx),[2,1,0])
rrcyl   = sqrt(xx^2+yy^2)
;
; Make a simple solid-body rotating gas flow
;
rhogas  = dblarr(nx,ny,nz) + rhogas0
tgas    = dblarr(nx,ny,nz) + temp0
vx      = dblarr(nx,ny,nz)
vy      = dvdy*yy 
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
         printf,1,tgas[ix,iy,iz]
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
; Write the lines.inp control file
;
openw,1,'lines.inp'
printf,1,'1'
printf,1,'1'
printf,1,'co    leiden    0    0'
close,1
;
; Write the radmc3d.inp control file
;
openw,1,'radmc3d.inp'
printf,1,'nphot = 1000000'
printf,1,'scattering_mode_max = 0'
;printf,1,'tgas_eq_tdust   = 1'
printf,1,'subbox_nx = ',64
printf,1,'subbox_ny = ',64
printf,1,'subbox_nz = ',64
printf,1,'subbox_x0 = ',-sizex
printf,1,'subbox_x1 = ',sizex
printf,1,'subbox_y0 = ',-sizey
printf,1,'subbox_y1 = ',sizey
printf,1,'subbox_z0 = ',-sizez
printf,1,'subbox_z1 = ',sizez
close,1
;
end
