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
nx       = 32L
ny       = 32L
nz       = 32L
sizex    = 10*AU
sizey    = 10*AU
sizez    = 10*AU
;
; The layer patch
;
laynlev  = 1
laynr    = 1
lix      = nx/4+1
liy      = ny/4+1
liz      = nz/4+1
lnx      = nx/2
lny      = ny/2
lnz      = nz/2
lnnx     = lnx*2
lnny     = lny*2
lnnz     = lnz*2
;
; Model parameters
;
radius   = 5*AU
rho0     = 1d-16
;
; Star parameters
;
mstar    = MS
rstar    = RS
tstar    = TS
pstar    = [0.,0.,0.]
;
; Make the coordinates
;
xi       = -sizex + 2*sizex*dindgen(nx+1)/(1.d0*nx)
yi       = -sizey + 2*sizey*dindgen(ny+1)/(1.d0*ny)
zi       = -sizez + 2*sizez*dindgen(nz+1)/(1.d0*nz)
xc       = 0.5d0 * ( xi[0:nx-1] + xi[1:nx] )
yc       = 0.5d0 * ( yi[0:ny-1] + yi[1:ny] )
zc       = 0.5d0 * ( zi[0:nz-1] + zi[1:nz] )
;
; Make the coordinates of the layer patch of refinement
;
xi_l1    = -sizex/2 + 2*sizex*dindgen(lnnx+1)/(1.d0*lnnx)/2
yi_l1    = -sizey/2 + 2*sizey*dindgen(lnny+1)/(1.d0*lnny)/2
zi_l1    = -sizez/2 + 2*sizez*dindgen(lnnz+1)/(1.d0*lnnz)/2
xc_l1    = 0.5d0 * ( xi_l1[0:lnnx-1] + xi_l1[1:lnnx] )
yc_l1    = 0.5d0 * ( yi_l1[0:lnny-1] + yi_l1[1:lnny] )
zc_l1    = 0.5d0 * ( zi_l1[0:lnnz-1] + zi_l1[1:lnnz] )
;
; Make the dust density model
;
xx       = rebin(xc,nx,ny,nz)
yy       = transpose(rebin(yc,ny,nx,nz),[1,0,2])
zz       = transpose(rebin(zc,nz,ny,nx),[2,1,0])
rr       = sqrt(xx^2+yy^2+zz^2)
rhod     = rho0 * exp(-(rr^2/radius^2)/2.d0)
;
; Make the dust density model in the layer
;
xx_l1    = rebin(xc_l1,lnnx,lnny,lnnz)
yy_l1    = transpose(rebin(yc_l1,lnny,lnnx,lnnz),[1,0,2])
zz_l1    = transpose(rebin(zc_l1,lnnz,lnny,lnnx),[2,1,0])
rr_l1    = sqrt(xx_l1^2+yy_l1^2+zz_l1^2)
rhod_l1  = rho0 * exp(-(rr_l1^2/radius^2)/2.d0)
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
; Write the grid file
;
openw,1,'amr_grid.inp'
printf,1,1                      ; iformat
printf,1,10                     ; AMR grid style  (10=layer-style AMR)
printf,1,0                      ; Coordinate system
printf,1,0                      ; gridinfo
printf,1,1,1,1                  ; Include x,y,z coordinate
printf,1,nx,ny,nz               ; Size of grid
printf,1,laynlev,laynr          ; Layers: nr of levels, nr of layers
for i=0,nx do printf,1,xi[i]    ; X coordinates (cell walls)
for i=0,ny do printf,1,yi[i]    ; Y coordinates (cell walls)
for i=0,nz do printf,1,zi[i]    ; Z coordinates (cell walls)
printf,1,0,lix,liy,liz,lnx,lny,lnz  ; Info for layer
close,1
;
; Write the density file
;
openw,1,'dust_density.inp'
printf,1,1                      ; Format number
printf,1,nx*ny*nz+lnnx*lnny*lnnz  ; Nr of spatial data points (incl redundant ones)
printf,1,1                      ; Nr of dust species
for iz=0,nz-1 do begin
   for iy=0,ny-1 do begin
      for ix=0,nx-1 do begin
         printf,1,rhod[ix,iy,iz]
      endfor
   endfor
endfor
for iz=0,lnnz-1 do begin
   for iy=0,lnny-1 do begin
      for ix=0,lnnx-1 do begin
         printf,1,rhod_l1[ix,iy,iz]
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
; Write the radmc3d.inp control file
;
openw,1,'radmc3d.inp'
printf,1,'nphot = ',nphot
printf,1,'scattering_mode_max = 0'
close,1
;
end
