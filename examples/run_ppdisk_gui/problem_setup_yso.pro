;pro problem_setup_yso
@problem_params.pro
;
; Put the star into its proper position in the model grid
;

pstar    = [px,py,pz]     ; Position of the star in the grid (x,y,z)

;
; Create the coordinate system
;
ri       = grid_rin * (grid_rout/grid_rin)^(dindgen(grid_nr+1)/(grid_nr-1))
ti       = rotate((!dpi/2.0 - !pi/2.*grid_tmax*dindgen(grid_nt+1)/(grid_nt*1.0)),2)
phi      = [0.,0.]
rc       = 0.5d0 * ( ri[0:grid_nr-1] + ri[1:grid_nr] )
tc       = 0.5d0 * ( ti[0:grid_nt-1] + ti[1:grid_nt] )
;
; Make the dust density model (ppdisk)
;
rho_ppdisk = dblarr(grid_nr, grid_nt)
if ppdisk_enable eq 1 then begin
   sigma = ppdisk_sig0 * (rc / ppdisk_rout)^ppdisk_plsig1
   hr    = ppdisk_hrdisk * (rc / grid_rout)^ppdisk_plh * rc
   for ir=0, grid_nr-1 do begin
      if rc[ir] ge ppdisk_rin and rc[ir] le ppdisk_rout then begin
         rho0 = sigma[ir] / (hr[ir] * sqrt(2d0*!dpi))
         rho1 = exp(-0.5*(tan(!dpi/2. - tc(*))*rc[ir]/hr[ir])^2)
         rho_ppdisk[ir,*] = rho0 * rho1
      endif
   endfor
endif
;
; Make the dust density model (envelope)
;
rho_envelope = dblarr(grid_nr, grid_nt)
if envelope_enable eq 1 then begin
   for ir=0, grid_nr-1 do begin
      if rc[ir] ge envelope_rin and rc[ir] le envelope_rout then $
         rho_envelope[ir,*] = envelope_rho0 * (rc[ir] / envelope_rout)^envelope_plrho
   endfor
   if envelope_cavrad gt 0. and envelope_cavrfact ge 0 then begin
      cavrad_rad = envelope_cavrad/180.*!dpi
      ic         = where(tc le cavrad_rad)
      for ir=0, grid_nr-1 do begin
         rho_envelope[ir,ic] = rho_envelope[ir,ic] * envelope_cavrfact
      endfor
   endif
endif
;
; Make background density
;
rho_bg = dblarr(grid_nr, grid_nt) + 1d-22

;
; Make the global density structure
;
rho = dblarr(grid_nr, grid_nt)

for ir=0, grid_nr-1 do begin
   for it=0, grid_nt-1 do begin
      dum = [rho_ppdisk[ir,it], rho_envelope[ir,it], rho_bg[ir,it]]
      ii  = where(dum eq max(dum))
      
      if ii(0) eq 0 then rho[ir,it] = rho_ppdisk[ir,it]
      if ii(0) eq 1 then rho[ir,it] = rho_envelope[ir,it]
      if ii(0) eq 2 then rho[ir,it] = rho_bg[ir,it]
   endfor
endfor

;
; Plot the density structure (debugging purpose)
;

;r = rc
;theta = !dpi/2.-tc
;dlev = reverse(findgen(18)*(-1.)-12.)

;window,0,retain=2, xsize=700, ysize=600
;loadct, 0
;contour, alog10(rho), r/au, theta, /xl, lev=dlev
;polar_contour, alog10(transpose(rho)), theta, r/au, nlev=80, /dither,  $
;               xr=[0,20], yr=[0,20]
;stop

;
; Write the wavelength_micron.inp file
;
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
printf,1,0                      ; AMR grid style  (0=regular grid, no AMR)
printf,1,100                    ; Coordinate system
printf,1,0                      ; gridinfo
printf,1,1,1,0                  ; Include x,y,z coordinate
printf,1,grid_nr,grid_nt,grid_np ; Size of grid
for i=0, grid_nr do printf,1,ri[i]    ; X coordinates (cell walls)
for i=0, grid_nt do printf,1,ti[i]    ; Y coordinates (cell walls)
for i=0, grid_np do printf,1,phi[i]    ; Z coordinates (cell walls)
close,1
;
; Write the density file
;
openw,1,'dust_density.inp'
printf,1,1                                ; Format number
print, grid_nr*grid_nt*grid_np
printf,1,grid_nr*grid_nt*grid_np          ; Nr of cells
printf,1,1                                ; Nr of dust species
for ip=0,grid_np-1 do begin
   for it=0,grid_nt-1 do begin
      for ir=0,grid_nr-1 do begin
         ;printf,1,rho[ip,it,ir]
;
; TODO : I need to change the order of indices in rho
;      to be [phi, theta, r] as originally used by Kees
;
         printf,1,rho[ir,it]
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

;dens =read_data(/ddens)
;window,1,retain=2, xsize=700, ysize=600
;loadct, 0
;contour, alog10(dens.rho), r/au, theta, /xl, lev=dlev

end
