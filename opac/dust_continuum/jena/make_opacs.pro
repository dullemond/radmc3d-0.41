@readkappa.pro
@natconst.pro
;
; Make a wavelength grid. Here we take it regular, but with many wavelength
; points, because for the dustkappa_*.inp files the wavelength grid will
; anyway be reduced to the model wavelength grid.
;
nf      = 2000
lam0    = 0.1
lam1    = 10000.0
lambda  = lam0*(lam1/lam0)^(dindgen(nf)/(nf-1.d0))
;
; Here we make a size distribution for the particles
;
na      = 5
amin    = 0.1
amax    = 1000.
amicron = amin*(amax/amin)^(dindgen(na)/(na-1.d0))
;
; We choose the optical constants and the specific weight of the material,
; which you can download from:
; http://www.astro.uni-jena.de/Laboratory/OCDB/newsilicates.html
; or more specifically from:
; http://www.astro.uni-jena.de/Laboratory/OCDB/data/silicate/pyrmg70.lnk
;
optname = 'pyrmg70.lnk'
swgt    = 3.01
;
; We write the wavelength grid for the dust opacity code
;
openw,1,'frequency.inp'
printf,1,nf
printf,1,' '
for inu=0,nf-1 do printf,1,1d4*cc/lambda[inu]
close,1
;
; We tell the dust opacity code which dust opacities to make
;
openw,1,'dust.inp'
for ia=0,na-1 do begin
   printf,1,optname
   printf,1,'"MIE"'
   q = alog10(amicron[ia])
   printf,1,'1  0.0 ',q,q,' 1 -3.5 ',swgt,' -2.0'
endfor
close,1
;
; We call the dust opacity code
;
spawn,'./makedust'
;
; In principle we are now done. The opacities dustkappa_*.inp have been
; written. But we may want to generate an average opacity dustkappa_av.inp
; for a given size distribution powerlaw. Here plaw=-3.5 is the standard
; MRN 1977 distribution, which is equivalent to Kolmogorov. 
;
plaw = -3.5
;
; Now make the weights
;
pwgt = (plaw-2.d0)/3.d0 + 2.d0
mass = (4.*!dpi/3.)*(amicron/1d4)^3*swgt
dum  = (mass/mass[0])^pwgt
mwgt = dum / total(dum)
;
; Read all the opacities and add them all up, using the mwgt as weight
;
lamb = dblarr(nf)
kabs = dblarr(nf)
ksca = dblarr(nf)
for ia=0,na-1 do begin
   o=readkappa(spec=ia+1)
   if ia eq 0 then begin
      plotkappa,o
      plotkappa,o,/abs,line=2,/oplot
      plotkappa,o,/sca,line=1,/oplot
   endif else begin
      plotkappa,o,/oplot
      plotkappa,o,/abs,line=2,/oplot
      plotkappa,o,/sca,line=1,/oplot
   endelse
   lamb = o.wave
   kabs = kabs + mwgt[ia] * o.cabs
   ksca = ksca + mwgt[ia] * o.csca
endfor
;
; Now write the average opacity to file
;
openw,1,'dustkappa_av.inp'
printf,1,2
printf,1,nf
for i=0,nf-1 do begin
   printf,1,lamb[i],kabs[i],ksca[i]
endfor
close,1
;
end
