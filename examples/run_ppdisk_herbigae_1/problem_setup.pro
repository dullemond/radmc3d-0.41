;==========================================================================
;               MAIN SETUP ROUTINE FOR THE 2-D DISK MODELS
;          EXECUTE THIS SCRIPT IN ORDER TO SET UP A SIMULATION
;                 (PARAMETERS ARE IN PROBLEM_PARAMS.PRO)
;==========================================================================
@problem_subroutines.pro
@problem_models.pro
@problem_useopac.pro
@problem_natconst.pro
;
close,/all
;
; Defaults
;
mugas  = 2.3
close,/all
;
; Read the parameters
;
@problem_params.pro
;
; Checks
;
if rdisk gt rout then begin
   print,'PROBLEM: Outer disk radius is larger than the outer grid radius!'
   stop
endif
;
; Mix dust opacities
;
if keyword_set(mixspecs) then begin
   sz = size(mixspecs)
   if sz[0] eq 1 then nmix=1 else nmix=sz[2]
   for imix=0,nmix-1 do mixopacities,mixspecs[*,imix],$
                   mixnames[imix],mixabun[*,imix]
endif
;
; Make the wavelength grid
;
lambda0  = 1d4
lambda1  = 40.0
lambda2  = 5.5
lambda3  = 0.1  
nfrnew01 = 20
nfrnew12 = 80
nfrnew23 = 30 
cc       = 2.9979e10
freq0    = 1e4*cc/lambda0
freq1    = 1e4*cc/lambda1
freq2    = 1e4*cc/lambda2
freq3    = 1e4*cc/lambda3
dlfreq01 = alog(freq1)-alog(freq0)
dlfreq12 = alog(freq2)-alog(freq1)
dlfreq23 = alog(freq3)-alog(freq2)
freqnw01 = exp((findgen(nfrnew01)/(nfrnew01))*dlfreq01+alog(freq0))
freqnw12 = exp((findgen(nfrnew12)/(nfrnew12))*dlfreq12+alog(freq1))
freqnw23 = exp((findgen(nfrnew23)/(nfrnew23-1.e0))*dlfreq23+alog(freq2))
freqnew  = [freqnw01,freqnw12,freqnw23]
freqnew  = rotate(freqnew,2)
nf       = n_elements(freqnew)
nlam     = nf
lambda   = 1d4*cc/freqnew
;
; Write the wavelength data file
;
openw,1,'wavelength_micron.inp'
printf,1,nf
printf,1,' '
for i=0,nf-1 do printf,1,lambda[i]
close,1
;
; Make the dustopac.inp file
;
nspec = n_elements(opacs)
spawn,'rm -f dustopac.inp'
openw,1,'dustopac.inp'
printf,1,'2               Format number of this file'
printf,1,strcompress(string(nspec,format='(I4)'),$
         /remove_all)+'               Nr of dust species'
printf,1,'============================================================================'
for iopac = 1,nspec do begin
    printf,1,'1               Way in which this dust species is read (1=kappa file)'
    printf,1,'0               0=Thermal grain'
    printf,1,opacs[iopac-1]+'               Extension of name of dustkappa_***.inp file'
    printf,1,'----------------------------------------------------------------------------'
endfor
close,1
;
; Write the star spectrum and info
;
if keyword_set(kurucz) then begin
   lstar = LS * (rstar/RS)^2 * (tstar/TS)^4
   kurucz,tstar,lstar,mstar,/worig,kuruczdir=kurdir
   openr,1,'starspectrum.inp'
   nf=0
   data=dblarr(2,nf)
   readf,1,data
   close,1
   starspec = transpose(data[1,*])
endif else begin
   starspec = dblarr(nlam)
   for inu=0,nlam-1 do begin
      starspec[inu] = 3.14159265359 * (rstar^2/pc^2) $
                      * bplanck(freqnew[inu],tstar)
   endfor
endelse
openw,1,'stars.inp'
printf,1,2
printf,1,1,nlam
printf,1,' '
printf,1,rstar,mstar,0.d0,0.d0,0.d0
printf,1,' '
for ilam=0,nlam-1 do printf,1,lambda[ilam]
printf,1,' '
for inu=0,n_elements(lambda)-1 do begin
   printf,1,starspec[inu]
endfor
close,1
;
; Find the opacity at the peak of the stellar spectrum
;
lampk = 0.55d0
kappa = dblarr(nspec)
for i=0,nspec-1 do begin
    o = readopacmaster('dustkappa_'+opacs[i]+'.inp')
    kappa[i] = findopac(o,lambda=lampk)
endfor
print,"lambda_peak = ",lampk," micron,  kappa = ",kappa
;
; Estimate rin from tin using a simple blackbody wall thing
;
; NOTE: In contrast to the original D&D models we now do not iterate
;       on the self-irradiation.
;
; NOTE: We now have a new option: pre-compute the optically thin
;       temp, and from there try to find the rin
;
if not keyword_set(thintin) then begin
    ;;
    ;; Use the D&D04 formula, but this time without self-irradiation
    ;;
    if keyword_set(tin) then begin
        rin  = rstar*(tstar/tin)^2
    endif
endif else begin
    ;;
    ;; Use the optically thin dust temperature to estimate the 
    ;; radius rin belonging to tin.
    ;;
    print,'PROBLEM: The optically thin real temperature mode for tin'
    print,'         is still under construction....'
    stop
endelse
;;
;; Compute the puffing-up radius
;;
if keyword_set(rpfrin) then begin
    if(rpfrin lt 1.d0) then begin
        print,'ERROR: rpfrin should not be smaller than 1 '
        print,'       (unless you put it to 0, i.e. deactivate puffing)'
        stop
    endif
    rpuff = rpfrin * rin
endif else begin
    rpuff = 0.d0
endelse
;;
;; Make the R-grid
;;
r      = make_rgrid(rin,rout,nr,rrefine=rrefine)
;;
;; Make the Theta-grid
;;
theta  = make_tgrid(hrgrid,nt,hrgmax=hrgmax,ntex=ntex,$
                    hrlg=hrlg,zrefine=zrefine)
;;
;; Just for safety, get nnr and nnt
;;
nnr    = n_elements(r)
nnt    = n_elements(theta)
;;
;; Now make a disk model
;;
if keyword_set(sigdust0) then begin
   ;;
   ;; Directly from sigdust0
   ;;
   rhodusttot = disk_model_1(r,theta,rdisk,sigdust0,plsig1,plsig2,$
                          hrdisk,plh,hrmin=hrmin,hrstore=hrstore,$
                          hrpuff=hrpuff,rpuff=rpuff,sigdust=sigdust)
endif else begin
   ;;
   ;; Compute sigdust0 from mdisk
   ;;
   sigdust00  = 1.d0
   rhodusttot = disk_model_1(r,theta,rdisk,sigdust00,plsig1,plsig2,$
                            hrdisk,plh,hrmin=hrmin,hrstore=hrstore,$
                            hrpuff=hrpuff,rpuff=rpuff,sigdust=sigdusttot)
   mddum      = integrate(r,2*pi*r*sigdusttot)*gastodust
   sigdust00  = mdisk / mddum
   rhodusttot = disk_model_1(r,theta,rdisk,sigdust00,plsig1,plsig2,$
                         hrdisk,plh,hrmin=hrmin,hrstore=hrstore,$
                         hrpuff=hrpuff,rpuff=rpuff,sigdust=sigdusttot)
endelse
;;
;; Compute the mass
;;
md = mass(rhodusttot,r,theta)
print,'Mdisk = ',md/MS,' Msun (using gas-to-dust=100)'
;;
;; Now split this dust density up in the various abundances
;;
if nspec gt 1 then begin
   abun = dblarr(nnr,nnt,nspec) + 1.d0
   for ispec=1,nspec-1 do begin
      q  = dblarr(nnr) + ab_ab0[ispec-1]
      ii = where(r gt ab_r0[ispec-1])
      q[ii] = ab_ab0[ispec-1]*(r[ii]/ab_r0[ispec-1])^ab_pl[ispec-1]
      q = q > ab_min[ispec-1]
      abun[*,*,ispec] = rebin(q,nnr,nnt)
      abun[*,*,0] = abun[*,*,0] - abun[*,*,ispec]
   endfor
   if min(abun[*,*,0]) lt 0.d0 then begin
      print,'ERROR: Total abundances dont add up!'
      stop
   endif
endif
;;
;; Now convert this into dustdens
;;
if keyword_set(abun) then begin
   rhodust = dblarr(nr,nt,nspec)
   for ispec=0,nspec-1 do begin
      rhodust[*,*,ispec] = rhodusttot * abun[*,*,ispec]
   endfor
endif else begin
   rhodust = rhodusttot
endelse
;;
;; Now smooth the inner rim if necessary
;;
if keyword_set(drsm) then begin
   smooth_rim,r,theta,rhodust,kappa,sigdust=sigdust,$
              drsm=drsm,tautol=tautol
endif
;
; Make the grid into RADMC-3D style
;
nx = nnr
ny = nnt
nz = 1
xi = dblarr(nx+1)
xi[1:nx-1] = sqrt(r[0:nx-2]*r[1:nx-1])
xi[0] = xi[1]^2/xi[2]
xi[nx] = xi[nx-1]^2/xi[nx-2]
yi = dblarr(ny+1)
yi[1:ny-1] = 0.5d0 * ( theta[0:ny-2] + theta[1:ny-1] )
yi[0] = 0.d0
yi[ny] = !dpi/2
zi = [0.d0,2*!dpi]
;
; Write the grid file
;
openw,1,'amr_grid.inp'
printf,1,1                      ; iformat
printf,1,0                      ; AMR grid style  (0=regular grid, no AMR)
printf,1,100                    ; Coordinate system
printf,1,0                      ; gridinfo
printf,1,1,1,0                  ; Include x,y,z coordinate
printf,1,nx,ny,nz                ; Size of grid
for i=0,nx do printf,1,xi[i],format='(E23.15)'    ; X coordinates (cell walls)
for i=0,ny do printf,1,yi[i],format='(E23.15)'    ; Y coordinates (cell walls)
for i=0,nz do printf,1,zi[i],format='(E23.15)'    ; Z coordinates (cell walls)
close,1
;
; Write the density file
;
openw,1,'dust_density.inp'
printf,1,1                      ; Format number
printf,1,nx*ny*nz               ; Nr of cells
printf,1,nspec                  ; Nr of dust species
for ispec=0,nspec-1 do begin
;for iz=0,nz-1 do begin
   for iy=0,ny-1 do begin
      for ix=0,nx-1 do begin
         printf,1,rhodust[ix,iy,ispec]
      endfor
   endfor
;endfor
endfor
close,1
;
; Write the radmc3d.inp control file
;
openw,1,'radmc3d.inp'
printf,1,'nphot = ',nphot
printf,1,'scattering_mode_max = ',scat
close,1
;;
;; Make optical depth array
;;
b  = {r:r,theta:theta,rho:rhodust}
tt =maketau(b,kappa)
;;
;; Write a message
;;
ddr=r(nr-1)/r(nr-2)-1.d0
print,'TAUR  = ',tt.taur(nr-1,nt-1)
print,'DR/R  = ',ddr
;;
;; Check the optical depth at the inner edge
;;
tr0 = tt.taur(1,nt-1)
if tr0 gt 1.d0 then begin
    print,'Optical depth error! Inner grid cell at equator optically thick'
    print,' Is this okay? (type 1)'
    read,i
    if i ne 1 then stop
endif
;
;
; Make plots and write info
;
if show eq 1 then begin
   pos = [0.1,0.55,0.95,0.95]
   xtitle = 'R [AU]'
   ytitle = '!4p!X/2-!4H!X = H/R'
   nnt = n_elements(theta)
   window,1,xsize=800,ysize=700
   !p.multi(2)=2
;   surface,tt.taur gt 1,r,theta,/xl,ax=50,/xs,/ys
   surface,tt.taur gt 1,r/au,!pi/2-theta,/xl,ax=90,az=0,/xs,/ys,position=pos,xtitle=xtitle,ytitle=ytitle
   contour,tt.taur gt 1,r/au,!pi/2-theta,/xl,/xs,/ys,position=pos,/noerase,/fill,$
          title='Grid and the !4s!X!DV!N>1 region'
;   surface,tt.taur gt 1,r/au,!pi/2-theta,xtitle=xtitle,ytitle=ytitle,/xl,ax=90,az=0,/xs,/ys,position=pos,/noerase
   plot,tt.taur(0:20,nnt-1),yrange=[0,6],psym=6,yminor=1,yticks=6,$
            xtitle='Radial grid points from inner edge',$ 
            ytitle='!4s!X!DV!N at midplane'
   oplot,tt.taur(0:20,nnt-1)
   !p.multi(2)=1
endif
;
end
