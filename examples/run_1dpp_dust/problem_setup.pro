;
; Some natural constants
;
AU  = 1.49598d13     ; Astronomical Unit       [cm]
pc  = 3.08572d18     ; Parsec                  [cm]
MS  = 1.98892d33     ; Solar mass              [g]
TS  = 5.78d3         ; Solar temperature       [K]
LS  = 3.8525d33      ; Solar luminosity        [erg/s]
RS  = 6.96d10        ; Solar radius            [cm]
cc  = 2.9979245800000d10      ; Light speed             [cm/s]
pi  = 3.1415926535897932385d0 
;
; Monte Carlo parameters
;
nphot    = 100000
;
; Grid parameters
;
nx       = 1L
ny       = 1L
nz       = 200L
sizex    = 1d90
sizey    = 1d90
sizez    = 10*AU
;
; Model parameters
;
rho0     = 1d-16
;
; Illuminating parallel flux beams
;
;   *** Uncomment one of the options below ***
;
;  Test 1: No illuminating beams from the top.
;
inclillum = [0.]
tempillum = [0.]
dilillum  = [0.]
;
;  Test 2: Vertical incident beam of T=100 K. We should see that
;          the temperature drops a bit at the edge, because even
;          though the flux is consistent with a blackbody of 100 K,
;          it is not the same as being in a thermal bath of 100 K.
;
;          You can see this as follows: A grey dust grain of radius
;          a in a parallel beam of flux of F_nu = pi * B_nu(100K) 
;          will absorb int_0^infty pi*B_nu(100K)*pi*a^2 dnu erg/s.
;          It will emit 4*pi*a^2*int_0^infty pi*B_nu(T) dnu erg/s.
;          This leads to T=(1/4)^(1/4)*100K=70K. In our example
;          the flux from this beam only needs to account for half
;          of the radiation (the other half comes from below). But
;          even if you have two of these beams, you still get just
;          84 K. This means that if you put the sigma_SB * T^4 
;          flux in a parallel beam, you will not get the same
;          effect as a thermal boundary. 
; 
;inclillum = [0.]
;tempillum = [100.]
;dilillum  = [1.d0]
;
;  Test 3: According to the Rybicki & Lightman book, the two-stream
;          approximation is a decent approximation of radiative 
;          transfer in plane-parallel geometries. The two rays have
;          angle cos(theta) = 1/sqrt(3) or cos(theta) = -1/sqrt(3).
;          According to this logic, if we make the incident beam
;          under such an angle (which is 54.735 degrees), and if
;          we boost the unprojected flux such that the projected
;          flux is still equal to the thermal flux of 100 K, then
;          this should be as good as a thermal boundary. 
;
;inclillum = [54.735d0]
;tempillum = [100.]
;dilillum  = [1.d0*sqrt(3.)]
;
;  Test 4: Take a much stronger flux that the 100 K flux in the
;          previous test cases. Put the inclination angle 0, i.e.
;          the radiation goes vertically downward into the atmosphere.
;          You see that the temperature on the top of the atmosphere
;          is clearly higher than the base, but at the very top the
;          temperature drops a bit again. That drop is because due to
;          the incidence being vertical, the energy is injected deep
;          into the atmosphere, where, due to random direction
;          re-emission, it gets a bit "stuck", and thus heats the
;          interior more than the surface.
;
;inclillum = [0.]
;tempillum = [400.]
;dilillum  = [1.d0]
;
;  Test 5: Same as Test 4, but now with grazing incidence angle. Here
;          you see the opposite happening: the very surface layer has
;          a temperature increase. This is because of the fact that 
;          dust grains at the surface see the full blast of the flux
;          while dust grains in the interior only "feel" the projected
;          flux.
;
;inclillum = [80.]
;tempillum = [400.]
;dilillum  = [1.d0]
;
;----------------------------------------------------------------------
;
; Count the illumination beams
;
nillum    = n_elements(inclillum)
if tempillum[0] eq 0.d0 then nillum=0
;
; Make the coordinates
;
zi       = sizez*dindgen(nz+1)/(1.d0*nz)
zc       = 0.5d0 * ( zi[0:nz-1] + zi[1:nz] )
;
; Make the dust density model
;
rhod     = rho0 
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
; Write the grid file
;
openw,1,'amr_grid.inp'
printf,1,1                      ; iformat
printf,1,0                      ; AMR grid style  (0=regular grid, no AMR)
printf,1,10                     ; Coordinate system
printf,1,0                      ; gridinfo
printf,1,0,0,1                  ; Include only z coordinate
printf,1,1,1,nz                 ; Size of grid
printf,1,-1d90                  ; X coordinates (cell walls)
printf,1,1d90                   ; X coordinates (cell walls)
printf,1,-1d90                  ; Y coordinates (cell walls)
printf,1,1d90                   ; Y coordinates (cell walls)
for i=0,nz do printf,1,zi[i]    ; Z coordinates (cell walls)
close,1
;
; Write the density file
;
openw,1,'dust_density.inp'
printf,1,1                      ; Format number
printf,1,nz                     ; Nr of cells
printf,1,1                      ; Nr of dust species
for iz=0,nz-1 do begin
   printf,1,rhod
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
; Add an illuminating beam
;
if nillum gt 0 then begin
   openw,1,'illum.inp'
   printf,1,2
   printf,1,nillum,nlam
   for illum=0,nillum-1 do printf,1,inclillum[illum],0.d0
   printf,1,' '
   for ilam=0,nlam-1 do printf,1,lambda[ilam]
   for illum=0,nillum-1 do begin
      printf,1,' '
      for ilam=0,nlam-1 do begin
         nu  = 1d4*cc/lambda[ilam]
         bpl = bplanck(nu,tempillum[illum]) > 1d-90
         printf,1,dilillum[illum]*pi*bpl[0]
      endfor
   endfor
   close,1
endif else begin
   spawn,'rm illum.inp'
endelse
;
; Write the radmc3d.inp control file
;
openw,1,'radmc3d.inp'
printf,1,'nphot = ',nphot
printf,1,'scattering_mode_max = 0'
printf,1,'iranfreqmode = 1'
printf,1,'thermal_boundary_zl = 100.'
;printf,1,'thermal_boundary_zr = 100.'    ; Uncomment this to do the T=100 K test
close,1
;
end
