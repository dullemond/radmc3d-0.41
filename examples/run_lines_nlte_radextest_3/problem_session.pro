@readradmc
@natconst
;
; Grid parameters
;
nx       = 1L
ny       = 1L
nz       = 1L
sizex    = 1000*AU   ; HALF-width
sizey    = 1000*AU
sizez    = 1000*AU
;
; Model parameters that have to be fed into RADEX by van der Tak
;
mol      = 'co'
tbg      = 2.73
temp     = 30.
nh2      = 1d-1
ncol     = 1d20
vturb    = 1.*1d5    ; Velocity in cm/s, hence the *1d5 (=1km/s)
;
; Compute some stuff
;
length   = 2*sizex          ; because sizex is half width
nmol     = ncol/(2*sizex)   ; because sizex is half width
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
; Write the grid file
;
openw,1,'amr_grid.inp'
printf,1,1                      ; iformat
printf,1,0                      ; AMR grid style  (0=regular grid, no AMR)
printf,1,0                      ; Coordinate system
printf,1,0                      ; gridinfo
printf,1,1,1,1                  ; Include x,y,z coordinate
printf,1,1,1,1                  ; Size of grid
printf,1,xi[0],xi[1]            ; X coordinates (cell walls)
printf,1,yi[0],yi[1]            ; Y coordinates (cell walls)
printf,1,zi[0],zi[1]            ; Z coordinates (cell walls)
close,1
;
; Write the molecule number density file. 
;
openw,1,'numberdens_'+mol+'.inp'
printf,1,1                      ; Format number
printf,1,1                      ; Nr of cells
printf,1,nmol
close,1
;
; Write the number density file for the collision partner 
;
openw,1,'numberdens_h2.inp'
printf,1,1                      ; Format number
printf,1,1                      ; Nr of cells
printf,1,nh2
close,1
;
; Write the microturbulence file
;
openw,1,'microturbulence.inp'
printf,1,1                      ; Format number
printf,1,1                      ; Nr of cells
printf,1,vturb
close,1
;
; Write the gas temperature file
;
openw,1,'gas_temperature.inp'
printf,1,1                      ; Format number
printf,1,1                      ; Nr of cells
printf,1,temp
close,1
;
; Write the length scale file
;
openw,1,'escprob_lengthscale.inp'
printf,1,1                      ; Format number
printf,1,1                      ; Nr of cells
printf,1,length
close,1
;
; Write the wavelength file
;
openw,1,'wavelength_micron.inp'
printf,1,nlam
for ilam=0,nlam-1 do printf,1,lambda[ilam]
close,1
;
; Write the lines.inp control file
;
openw,1,'lines.inp'
printf,1,'2'
printf,1,'1'
printf,1,mol+'    leiden    0    0    1'
printf,1,'h2'
close,1
;
; Write the radmc3d.inp control file
;
openw,1,'radmc3d.inp'
printf,1,'lines_mode = 3'               ;; LVG mode
printf,1,'lines_nonlte_maxiter = 100'    ;; Max nr of iterations
close,1
;
; Call RADMC3D
;
spawn,'radmc3d calcpop'
;
; Read the results
;
a=read_levelpop('co',/ext)
print,a.ext
;
end
