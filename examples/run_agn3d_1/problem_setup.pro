;==========================================================================
;               MAIN SETUP ROUTINE FOR THE AGN TORUS MODELS
;                 (PARAMETERS ARE IN PROBLEM_PARAMS.PRO)
;==========================================================================
@problem_subroutines.pro
@problem_build.pro
@readradmc.pro
;;@readopac.pro
;
; Natural constants
;
AU     = 1.496d13
cc     = 2.9979d10
RS     = 6.96d10
MS     = 1.99d33
LS     = 3.8525d33     
pc     = 3.08572d18
ss     = 5.6703d-5  
kk     = 1.3807d-16   
mp     = 1.6726d-24  
GG     = 6.672d-8 
TS     = 5.78d3
;
; Other constants
;
mugas  = 2.3
;
; Read the parameters
;
@problem_params.pro
;
; Self-consistency checks on parameters
;
if rin ne 0.d0 and tin ne 0.d0 then begin
   print,'PROBLEM: Cannot specify rin AND tin. Must specify'
   print,'         rin OR tin...'
   stop
endif
if mdisk gt 0.d0 and sig0 gt 0.d0 then begin
   print,'PROBLEM: In problem_params.pro one must specify EITHER'
   print,'    the disk mass mdisk, OR the Sigma(r0) sig0. Not both!'
   stop
endif
;
; Write the wavelength_micron.inp file
;
lambda1 = 0.001d0
lambda2 = 7.0d0
lambda3 = 25.d0
lambda4 = 1.0d4
n12     = 40
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
; Read the opacity
;
;o         = readopacmaster('dustkappa_laordraine_0.1mic.inp')
o         = readopac(spec='laordraine_0.1mic')
;
; Create AGN spectrum of Granato & Danese (1994)
;
nu        = 1d4*cc/lambda
lnu       = nu*0.d0
;
nu1       = 10^(13.5)
nu2       = 10^(14.5)
pl        = 2.0
a         = 1.d0
ii        = where(nu gt nu1 and nu le nu2)
lnu[ii]   = a * (nu[ii]/nu1)^pl
lnu2      = a * (nu2/nu1)^pl
;
nu1       = 10^(14.5)
nu2       = 10^(15.4)
pl        = -0.5
a         = lnu2
ii        = where(nu gt nu1 and nu le nu2)
lnu[ii]   = a * (nu[ii]/nu1)^pl
lnu2      = a * (nu2/nu1)^pl
;
nu1       = 10^(15.4)
nu2       = 10^(16.0)
pl        = -1.0
a         = lnu2
ii        = where(nu gt nu1 and nu le nu2)
lnu[ii]   = a * (nu[ii]/nu1)^pl
lnu2      = a * (nu2/nu1)^pl
;
nu1       = 10^(16.0)
nu2       = 10^(20.0)
pl        = -2.2
a         = lnu2
ii        = where(nu gt nu1 and nu le nu2)
lnu[ii]   = a * (nu[ii]/nu1)^pl
lnu2      = a * (nu2/nu1)^pl
;
lnutot    = integrate(nu,lnu)
lnu       = lnu*(lagn/lnutot)
;
; Write the stars.inp file
;
openw,1,'stars.inp'
printf,1,2
printf,1,1,nlam
printf,1,' '
printf,1,rstar,mstar,0.d0,0.d0,0.d0
printf,1,' '
for ilam=0,nlam-1 do printf,1,lambda[ilam]
printf,1,' '
for inu=0,n_elements(lambda)-1 do begin
   printf,1,lnu[inu]/1.1965280793d38
endfor
close,1
;
; Write the dustopac.inp
;
openw,1,'dustopac.inp'
printf,1,'2               Format number of this file'
printf,1,'1               Nr of dust species'
printf,1,'============================================================================'
printf,1,'1               Way in which this dust species is read'
printf,1,'0               0=Thermal grain'
printf,1,'laordraine_0.1mic     Extension of name of dustkappa_***.inp file'
printf,1,'----------------------------------------------------------------------------'
close,1
;
; Find the opacity at this frequency (wavelength)
;
inupk  = find_freq(o.freq,1d4*cc/0.55d0)
kappa  = o.cabs(inupk) + o.csca(inupk)
print,"lambda = 0.55 micron -->  kappa = ",kappa
;
; Do a big loop for the setup. 
;
if mdisk gt 0.d0 then sig0=1.d0
a = problem(nr,nt,np,rin,rout,r0,hrgrid,hrdisk,plh,sig0,plsig1,plsig2,$
            nlevr=nlevr,nspanr=nspanr,nstepr=nstepr,hrmin=hrmin,$
            ntextra=ntex,hrgmax=hrgmax,nblob=nblob,sblob=sblob,$
            iseed=iseed,istyle=istyle,tq=tq,cdens=cdens,$
            m_random=m_random,r_random=r_random,t_random=t_random,$
            nomirror=nomirror)
nnr   = n_elements(a.r)
nnt   = n_elements(a.theta)
rr    = rebin(a.r,nnr,nnt)
tt    = transpose(rebin(a.theta,nnt,nnr))
;;xx    = rr*sin(tt)
;;yy    = rr*cos(tt)
if mdisk gt 0.d0 then begin 
   md    = mass(a.rho,a.r,a.theta,a.phi)
   sig0  = mdisk/md
   a.sigma=a.sigma*mdisk/md
   a.rho = a.rho*mdisk/md
   if nblob gt 0 and sblob gt 0.d0 then begin
      m_random = m_random*mdisk/md
      tq       = tq*mdisk/md
      for i=0,n_elements(tq)-1 do print,'r = ',r_random[i]/pc,$
        ' m = ',m_random[i]/MS,$
        ' tau_V = ',tq[i]*kappa
      openw,1,'blobinfo.dat'
      printf,1,n_elements(tq)
      for i=0,n_elements(tq)-1 do printf,1,r_random[i],$
              m_random[i],tq[i]*kappa
      close,1
   endif
endif
;
; A-posteriori smoothing of density
;
if drsm ne 0.d0 then begin
   ;;
   ;; Make a smooth curve from 0 to 1, saturating at 1 around
   ;; R=R_0*(1+drsm), i.e. the point where the reduction should stop.
   ;; Using this smooth function we can prevent the reduction from
   ;; being abrupt.
   ;;
   r1 = a.r[0]*(1.d0+drsm)
   x  = (a.r-a.r[0])/(r1-a.r[0])
   q1 = 2.d1
   q2 = 1.d0
   z  = -alog10(exp(-x*q1)+10^(-q2))/q2   ; Make smooth curve
   z  = 1-(1-z)/(1-z[0])                  ; Make sure z[0]==0
   ;;
   ;; Now reduce the density.
   ;;
   if not keyword_set(rdxsm) then rdxsm=1d-2
   for ir=0,nnr-1 do begin
      x = 1.d0-z[ir]
      redux=rdxsm^x
      a.rho(ir,*)=a.rho(ir,*)*redux
      a.sigma(ir)=a.sigma(ir)*redux
   endfor
endif
;
; Compute the mass
;
md = mass(a.rho,a.r,a.theta,a.phi)
print,'Mdisk = ',md/MS
;;
;; Make optical depth array
;;
;t=maketau(a,kappa)
;;
;; Write a message
;;
ddr=a.r(nnr-1)/a.r(nnr-2)-1.d0
;print,'TAUR  = ',t.taur(nnr-1,nnt-1)
;print,'DR/R  = ',ddr
;
; Do tests
;
;if ddr lt 0.05 then begin
;   print,'Radial grid too finely spaced... This will take too much'
;   print,'computational time... Are you sure this is okay? (Type 1)'
;   read,i
;   if i ne 1 then stop
;endif
if ddr gt 0.15 then begin
   print,'Radial grid too coarse! Are you sure this is okay? (Type 1)'
   read,i
   if i ne 1 then stop
endif
;;
;; Check the optical depth at the inner edge
;;
;tr0 = t.taur(1,nnt-1)
;if tr0 gt 1.d0 then begin
;   print,'Optical depth error! Inner grid cell at equator optically thick'
;   print,' Is this okay? (type 1)'
;   read,i
;   if i ne 1 then stop
;endif
;
; Make the grid into RADMC-3D style
;
nx = nnr
ny = nnt
nz = np
xi = dblarr(nx+1)
xi[1:nx-1] = sqrt(a.r[0:nx-2]*a.r[1:nx-1])
xi[0] = xi[1]^2/xi[2]
xi[nx] = xi[nx-1]^2/xi[nx-2]
yi = dblarr(ny+1)
yi[1:ny-1] = 0.5d0 * ( a.theta[0:ny-2] + a.theta[1:ny-1] )
yi[0] = 0.d0
if a.theta[ny-1] gt !dpi/2.d0 then begin
   yi[ny] = !dpi
endif else begin
   yi[ny] = !dpi/2
endelse
zi = dblarr(nz+1)
zi[0] = 0.d0
zi[np] = 2*!dpi
zi[1:nz-1] = 0.5d0 * ( a.phi[1:nz-1] + a.phi[0:nz-2] ) 
;
; Write the grid file
;
print,'Writing amr_grid.inp...'
openw,1,'amr_grid.inp'
printf,1,1                      ; iformat
printf,1,0                      ; AMR grid style  (0=regular grid, no AMR)
printf,1,100                    ; Coordinate system
printf,1,0                      ; gridinfo
printf,1,1,1,1                  ; Include x,y,z coordinate
printf,1,nx,ny,nz                ; Size of grid
for i=0,nx do printf,1,xi[i],format='(E23.15)'    ; X coordinates (cell walls)
for i=0,ny do printf,1,yi[i],format='(E23.15)'    ; Y coordinates (cell walls)
for i=0,nz do printf,1,zi[i],format='(E23.15)'    ; Z coordinates (cell walls)
close,1
;
; Write the density file
;
print,'Writing dust_density.inp...'
openw,1,'dust_density.inp'
printf,1,1                      ; Format number
printf,1,nx*ny*nz               ; Nr of cells
printf,1,1                      ; Nr of dust species
for iz=0,nz-1 do begin
   for iy=0,ny-1 do begin
      for ix=0,nx-1 do begin
         printf,1,a.rho[ix,iy,iz]
      endfor
   endfor
endfor
close,1
;
; Write the radmc3d.inp control file
;
openw,1,'radmc3d.inp'
printf,1,'nphot = ',nphot
printf,1,'scattering_mode_max = ',scat
close,1
;
; Make plots and write info
;
;if show eq 1 then begin
;   window,1,xsize=500,ysize=700
;   !p.multi(2)=2
;   surface,t.taur gt 1,a.r,a.theta,/xl,ax=50,/xs,/ys
;   plot,t.taur(0:20,nnt-1),yrange=[0,6],psym=6,yminor=1,yticks=6
;   oplot,t.taur(0:20,nnt-1)
;   !p.multi(2)=1
;;   if nblob ne 0 then begin
;;      window,2,xsize=500,ysize=500
;;;      tvscl,rebin(sqrt(a.rho*rr^2),4*nr,4*(nt+ntex))
;;      xmax=1.2*r0
;;      theta=a.theta
;;      theta[0]=0.
;;      theta[n_elements(theta)-1]=!pi/2
;;      img=polar_surface(-sqrt(a.rho*a.rr^2),a.r,!pi/2-theta,$
;;           /grid,bounds=[0,0,xmax,xmax],$
;;           spacing=[1,1]*xmax/500)
;;      tvscl,img
;;      ;;
;;      ;; Now to postscript
;;      ;;
;;      thick=4.0                 ; Postscript line thickness
;;      size=1.4                  ; Postscript char size
;;      !p.thick=thick
;;      !x.thick=thick
;;      !y.thick=thick
;;      !z.thick=thick
;;      !p.charthick=thick
;;      !p.charsize=size
;;      set_plot,'ps'
;;      device,xsize=15,ysize=15,file='densdist.ps'
;;      plot,[0,1]*xmax/pc,[0,1]*xmax/pc,position=[0.15,0.15,0.9,0.9],$
;;             /nodata,/xs,/ys,xminor=1,yminor=1,xtitle='R!Dc!N [pc]',ytitle='Z [pc]'
;;      tv,bytscl(img),0.15,0.15,xsize=0.75,ysize=0.75,/norm
;;      phi=(dindgen(30)/29.)*!pi/2
;;      xc=[0,cos(phi),0]*1.01*min(a.r)/pc
;;      yc=[0,sin(phi),0]*1.01*min(a.r)/pc
;;      polyfill,xc,yc,color=255
;;      plot,[0,1]*xmax/pc,[0,1]*xmax/pc,position=[0.15,0.15,0.9,0.9],$
;;             /nodata,/xs,/ys,xminor=1,yminor=1,/noerase
;;      device,/close
;;      set_plot,'x'
;;      !p.thick=1
;;      !x.thick=1
;;      !y.thick=1
;;      !z.thick=1
;;      !p.charthick=1
;;      !p.charsize=1
;;      ;;
;;   endif
;;endif
;
end
