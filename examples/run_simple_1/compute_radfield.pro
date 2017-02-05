@readradmc
@bplanck
@natconst
;;
;; This is an example of how to compute the mean intensity radiation field
;; at any given location using the "radmc3d mcmono" command.
;;
;; First we have to decide at which wavelengths we wish to compute the
;; radiation field. Let us take to wavelengths: 3 micron and 10 micron.
;; 
lambda = [3.d0,10.d0]
;;
;; Now write these to a file called mcmono_wavelength_micron.inp
;;
nlam = n_elements(lambda)
openw,1,'mcmono_wavelength_micron.inp'
printf,1,nlam
for i=0,nlam-1 do begin
   printf,1,lambda[i]
endfor
close,1
;;
;; Now call RADMC-3D to compute the mean intensities
;;
;;          1  /
;;  J_nu = --- 0 I_nu(Omega) dOmega
;;         4pi /
;;
;; at these wavelengths.
;;
spawn,'radmc3d mcmono nphot_mono 1000000'
;;
;; Read the mean intensity file
;;
a=read_data(/meanint)
;;
;; Plot a slice
;;
surface,a.jnu[*,*,16,0],/zl
;;
;; If the model is optically thin (NOTE: The original run_simple_1
;; is not!), then we can compare J_nu to what you can compute analytically
;; from the stellar flux.
;;
print,'RADMC-3D gives: ',a.jnu[0,16,16,0],a.jnu[0,16,16,1]
r   = abs(a.grid.x[0])
nu  = 1d4*cc/lambda
bnu = bplanck(nu,TS)
lnu = !dpi*bnu*4*!dpi*RS^2
fnu = lnu/(4*!dpi*r^2)
jnu = fnu/(4*!dpi)
print,'The star alone: ',jnu
;;
;; For the case rho0 = 1d-30 (in problem_setup.pro) we get as output:
;;
;;  RADMC-3D gives:   6.8192669e-13   7.9571554e-14
;;  The star alone:   6.5635271e-13   8.1052783e-14
;;
;; which I think is an ok agreement.
;;
;; For the original case of rho0 = 1d-16 we get:
;;
;;  RADMC-3D gives:    1.3994859e-12   9.2372747e-14
;;  The star alone:    6.5635271e-13   8.1052783e-14
;;
;; which shows that the emissivity by the multiple scattering and
;; thermal emissivity dominates at 3 micron.
;;
;; 13.03.2012
;;
end
