@problem_subroutines.pro

;-------------------------------------------------------------------------
;                      SET UP THE KURUCZ SPECTRUM
;
; This routine sets up the Kurucz spectrum for a given stellar effective
; temperature [K], stellar luminosity [erg/s] and mass [g]. The output
; is starspectrum.inp. The original Kurucz model (interpolated between
; the discrete parameters, but still on the high-resolution frequency 
; grid) can be written to a file by '/worig', or extracted with 
; 'model=model'. If the user wishes that automatic plots are made of the
; original model, the coarse-freq-grid version and the corresponding 
; blackbody curve, this can be done with 'silent=0'.
;-------------------------------------------------------------------------
pro kurucz,tstar,lstar,mstar,model=model,silent=silent,worig=worig,$
               kuruczdir=kuruczdir
;
; Defaults
;
if n_elements(kuruczdir) eq 0 then begin
   kuruczdir = '/home/dullemon/science/kurucz/'
endif
if n_elements(silent) eq 0 then silent = 1
;
; Natural constants
;
LS       = 3.86d33
ss       = 5.671d-5
GG       = 6.673d-8
cc       = 2.998d10
pc       = 3.08572d18
;
; For safety: close all
;
close,/all
;
; Read the frequency table
;
openr,1,'wavelength_micron.inp'
nf=0
readf,1,nf
freq=dblarr(nf)
readf,1,freq
close,1
freq = 1d4*cc/freq
;
; Compute certain quantities for the star
; 
Rstar = sqrt(Lstar/4./!DPI/ss/Tstar^4)
logg = alog10(GG*Mstar/Rstar^2)
;
; First clear possible old files
;
spawn,'rm -f flux.dat >& /dev/null'
;
; Copy the interpolator here
;
spawn,'cp -f '+kuruczdir+'interpolate_kurucz ./ >& /dev/null' 
;
; Command line for the kurucz interpolator
;
cmd = "./interpolate_kurucz -q "+string(format='(F6.0,1x,F5.1,1x,F5.1,A10)', $
                 Tstar,logg,0.0, $
                 "> flux.dat")
spawn,cmd
;
; Read the flux from the Kurucz interpolator
;
openr,1,'flux.dat'
a=dblarr(2,1221)
readf,1,a
close,1
loglam = a[0,*]
logflux = a[1,*]
flux = 10^logflux
lam = 10^loglam*1e-4
nu = cc/lam
s1 = nu
s2 = flux
LL = int_tabulated(s1,s2,/sort,/double)
lnu = flux*Lstar/LL
s1 = nu
s2 = lnu
LL1 = int_tabulated(s1,s2,/double,/sort)
;print,LL,LL1,LL1/Lstar
lnusmooth = lnu
imin = 0
imax=1200
;
; Print info, and plot the Kurucz model
;
if silent ne 1 then begin
   ;;print,nu[0],nu[imax]
   ;;print,lam[0],lam[imax]
   dumm=nu*lnusmooth/LS
   plot,3d14/nu,dumm,/xlog,/ylog,yrange=[1d-14,1]*max(dumm)
endif
;
; If the user wishes, the `real' Kurucz model is returned via the
; keywords
; 
lam_kur=transpose(1d4*cc/nu)
nulnu_kur=transpose(nu*lnusmooth)
model={lambda:lam_kur,nulnu:nulnu_kur}
;
; Write this real Kurucz model
;
if keyword_set(worig) then begin
   factor=1.d0/(4*!dpi*pc^2)
   data=dblarr(2,n_elements(nu))
   data[0,*]=nu[0,*]
   data[1,*]=lnusmooth[0,*]*factor
   openw,1,'kuruczmodel.dat'
   printf,1,n_elements(nu)
   printf,1,data,format='(2(E13.6,1X))'
   close,1
endif else begin
   spawn,'rm -f kuruczmodel.dat >& /dev/null'
endelse
;
; Interpolate on the coarse grid for the simulations
;
lint = interpol(lnusmooth,nu,freq)
ii=where(lint lt 0.d0)
if ii[0] ge 0 then lint[ii]=0.d0
;
; For wavelengths longer than 100 micron we simply extrapolate
;
ii=where( freq lt 1d4*cc/100d0 )
if ii[0] ge 0 then begin
   if ii[0] ne 0 then begin
      print,'ERROR: I think that the frequency grid is in the wrong'
      print,'  order... It should be in increasing frequency.'
      stop
   endif
   l0=lint[max(ii)+1]
   n0=freq[max(ii)+1]
   lint[ii]=l0*(freq[ii]/n0)^2
endif
;
; Overplot this interpolated spectrum
;
if silent ne 1 then begin
   oplot,3d14/freq,freq*lint/LS,psy=2
endif
;
; Overplot a blackbody
;
if silent ne 1 then begin
   bnu = freq
   for k=0,n_elements(freq)-1 do bnu(k) =bplanck(freq[k],Tstar)
   oplot,3d14/freq,freq*bnu*4*!DPI^2*Rstar^2/LS,line=1
endif
;
; Write the spectrum ot starspectrum.inp
;
filename='starspectrum.inp'
openw,2,filename
printf,2,n_elements(freq)
factor=1.d0/(4*!dpi*pc^2)
for j = 0,n_elements(lint)-1 do begin
   printf,2,freq[j],lint[j]*factor
endfor
close,2
;
; Check the error
;
s1 = freq
s2 = lint
LL1 = int_tabulated(s1,s2,/double,/sort)
print,' '
print, "   (Kurucz model relative error = ",abs(LL1/Lstar-1.d0),")"
print,' '
; 
end
  
