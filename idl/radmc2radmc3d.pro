;==========================================================================
;                    CONVERT FROM RADMC TO RADMC-3D
;==========================================================================
;
; Some constants
;
cc  = 2.9979245800000d10      ; Light speed             [cm/s]
;
; Read the radial grid of RADMC
;
openr,1,'radius.inp'
nr=0LL
readf,1,nr
r=dblarr(nr)
readf,1,r
close,1
;
; Read the theta-grid of RADMC
;
openr,1,'theta.inp'
nt=0LL
idum=0
readf,1,nt,idum
theta=dblarr(nt)
readf,1,theta
close,1
;
; Create the ri and thetai
;
ri=dblarr(nr+1)
ri[1:nr-1] = sqrt(r[0:nr-2]*r[1:nr-1])
ri[0]=r[0]
ri[nr]=r[nr-1]
thetai=dblarr(nt+1)
thetai[1:nt-1] = 0.5d0 * (theta[0:nt-2]+theta[1:nt-1])
thetai[0] = theta[0]
thetai[nt] = !dpi/2.d0
;
; Create the amr_grid.inp file for RADMC-3D
;
openw,1,'amr_grid.inp'
printf,1,1                      ; iformat
printf,1,0                      ; AMR grid style  (0=regular grid, no AMR)
printf,1,100                    ; Coordinate system
printf,1,0                      ; gridinfo
printf,1,1,1,0                  ; Include x,y,z coordinate
printf,1,nr,nt,1                ; Size of grid
for i=0,nr do printf,1,ri[i],format='(E23.15)'      ; R coordinates (cell walls)
for i=0,nt do printf,1,thetai[i],format='(E23.15)'  ; Theta coordinates (cell walls)
printf,1,0.d0,format='(E23.15)'                     ; Phi coordinates (cell walls)
printf,1,2*!dpi,format='(E23.15)'                    ; Phi coordinates (cell walls)
close,1
;
; Read the starinfo.inp of RADMC
;
openr,1,'starinfo.inp'
idum=0
readf,1,idum
rstar=0.d0
mstar=0.d0
tstar=0.d0
readf,1,rstar
readf,1,mstar
readf,1,tstar
close,1
;
; Read the star spectrum of RADMC
;
openr,1,'starspectrum.inp'
nf=0
readf,1,nf
data=dblarr(2,nf)
readf,1,data
close,1
lambda=1d4*cc/transpose(data[0,*])
starspec=transpose(data[1,*])
;
; Create the stars.inp file for RADMC-3D
;
openw,1,'stars.inp'
printf,1,2
printf,1,1,nf
printf,1,' '
printf,1,rstar,mstar,0.d0,0.d0,0.d0
printf,1,' '
for ilam=0,nf-1 do printf,1,lambda[ilam]
printf,1,' '
for inu=0,n_elements(lambda)-1 do begin
   printf,1,starspec[inu]
endfor
close,1
;
; Read the dust density of RADMC
;
openr,1,'dustdens.inp'
nspec=0
nrr=0
ntt=0
idum=0
readf,1,nspec,nrr,ntt,idum
if nrr ne nr or ntt ne nt then stop
rhodust=dblarr(nt,nr,nspec)
readf,1,rhodust
close,1
;
; Transpose this file, because RADMC-3D has ir as inner loop
;
if nspec gt 1 then begin
   rhodust=transpose(rhodust,[1,0,2])
endif else begin
   rhodust=transpose(rhodust,[1,0])
endelse
;
; Write the dust density file of RADMC-3D
;
openw,1,'dust_density.inp'
printf,1,1                      ; Format number
printf,1,nr*nt                  ; Nr of cells
printf,1,nspec                  ; Nr of dust species
for ispec=0,nspec-1 do begin
   for it=0,nt-1 do begin
      for ir=0,nr-1 do begin
         printf,1,rhodust[ir,it,ispec]
      endfor
   endfor
endfor
close,1
;
; Read the dustopac.inp file
; NOTE: RADMC-3D also uses this file
;
openr,1,'dustopac.inp'
iformat=0
nspecc=0
readf,1,iformat
readf,1,nspecc
if nspecc ne nspec then stop
str=' '
iext=intarr(nspec)
for ispec=0,nspec-1 do begin
   readf,1,str
   idum=0
   readf,1,idum
   readf,1,idum
   if idum ne 0 then begin
      printf,1,'Sorry, for now I can convert only thermal grain opacities'
      stop
   endif
   idum=0
   readf,1,idum
   iext[ispec] = idum
endfor
close,1
;
; Open the frequency.inp file of RADMC
;
openr,1,'frequency.inp'
nff=0
readf,1,nff
if nff ne nf then stop
freq=dblarr(nf)
readf,1,freq
close,1
;
; Do a test
;
dum = abs(1d4*cc/freq/lambda-1.d0)
if max(dum) gt 1d-4 then begin
   print,'ERROR: frequency.inp incompatible with starspectrum.inp'
   stop
endif
;
; Do a slight adaption of the frequency grid, just in case
;
if freq[0] lt freq[nf-1] then begin
   freq[0] = freq[0]       * 0.999999d0
   freq[nf-1] = freq[nf-1] * 1.000001d0
endif else begin
   freq[0] = freq[0]       * 1.000001d0
   freq[nf-1] = freq[nf-1] * 0.999999d0
endelse
;;
;; Now do a loop over the opacities and read the opacities and write them
;; to RADMC-3D style
;;
;for ispec=0,nspec-1 do begin
;   ;;
;   ;; Make the opacity file name
;   ;;
;   filein='dustopac_'+strcompress(string(ispec+1),/remove_all)+'.inp'
;   fileout='dustkappa_'+strcompress(string(ispec+1),/remove_all)+'.inp'
;   ;;
;   ;; Read the opacity from old RADMC style
;   ;;
;   openr,1,filein
;   nff=0
;   idum=0
;   readf,1,nff,idum
;   if nff ne nf then stop
;   kap_abs = dblarr(nf)
;   kap_sca = dblarr(nf)
;   readf,1,kap_abs
;   readf,1,kap_sca
;   close,1
;   ;;
;   ;; Write the opacity
;   ;;
;   openw,1,fileout
;   printf,1,2
;   printf,1,nf
;   for inu=0,nf-1 do begin
;      printf,1,1d4*cc/freq[inu],kap_abs[inu],kap_sca[inu]
;   endfor
;   close,1
;endfor
;
; Write radmc3d.inp (this overrides radmc.inp in RADMC-3D)
;
print,'Warning: Creating entirely new radmc3d.inp file, independent of radmc.inp'
openw,1,'radmc3d.inp'
printf,1,'nphot = 1000000'
printf,1,'istar_sphere = 1'
printf,1,'scattering_mode_max = 1'
close,1
;
; Done.
;
end







