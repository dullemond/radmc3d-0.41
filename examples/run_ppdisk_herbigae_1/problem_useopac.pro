@interpol.pro
;@problem_pah.pro

;-------------------------------------------------------------------------
;                      SET UP THE DUST OPACITY
;
; This routine sets up the dust opacity files dustopac.inp and 
; dustopac_XXX.inp for the RADMC and RADICAL codes. It takes opacity
; tables from files given by the input variable 'inputfiles', and remaps
; them onto the grid specified here. Each of the input files can have 
; its own frequency grid. 
;
; ARGUMENTS:
;   inputfiles               Array of files of opacity tables to be used
;   pllongs                  Array of power law indices for the opacities
;                            at long wavelengths.
;
; KEYWORDS:
;   fresmd                   Frequency resolution mode index
;   scat                     Include scattering opacity?
;   pnc                      (For inclusion of PAHs only)
;                               Nr of C-atoms of PAH
;   pnh                      (For inclusion of PAHs only)
;                               Nr of H-atoms of PAH
;   pz                       (For inclusion of PAHs only)
;                               Charge of PAH
;
; RETURNS:
;   nf                       Nr of frequencies used
;   nspec                    Nr of dust species
;   ntherm                   Nr of thermal dust species (without
;                              PAHs this equals nspec)
;
; NOTE: When onlynf is set, then only return nf! No opacity is then
;       generated. This is used for the problem_compilecodes.pro
;
;-------------------------------------------------------------------------
pro useopac,inputfiles,pllongs,fresmd=fresmd,scat=scat,$
            onlynf=onlynf,nf=nf,pnc=pnc,pnh=pnh,pz=pz,$
            nspec=nspec,npah=npah,ntherm=ntherm
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
;
; Defaults
;
incscat  = 0.d0
;
; Include scattering?
;
if keyword_set(scat) then incscat=1.d0
;
; Make the frequency array
;   (There is a discrete set of possible frequency grid resolution modes)
; 
if n_elements(fresmd) eq 0 then fresmd=0
case fresmd of
   0: begin
      lambda0  = 1d4
      lambda1  = 50.0
      lambda2  = 6.0
      lambda3  = 0.1
      nfrnew01 = 24
      nfrnew12 = 22
      nfrnew23 = 20
   end
   1: begin
      lambda0  = 1d4
      lambda1  = 50.0
      lambda2  = 6.0
      lambda3  = 0.1
      nfrnew01 = 20
      nfrnew12 = 80
      nfrnew23 = 20
   end
   2: begin            ; Freq-grid used for the group I/II simuls (BB-star)
      lambda0  = 1d4
      lambda1  = 30.0
      lambda2  = 6.0
      lambda3  = 0.1
      nfrnew01 = 20
      nfrnew12 = 26
      nfrnew23 = 19
   end
   3: begin            ; Freq-grid used for the group I/II simuls (Kurucz)
      lambda0  = 1d4
      lambda1  = 30.0
      lambda2  = 5.5
      lambda3  = 0.1  
      nfrnew01 = 20
      nfrnew12 = 30
      nfrnew23 = 30 
   end
   4: begin            ; Freq-grid used for the group I/II simuls (Kurucz)
      lambda0  = 1d4
      lambda1  = 30.0
      lambda2  = 7.0
      lambda3  = 0.1  
      nfrnew01 = 20
      nfrnew12 = 30
      nfrnew23 = 30 
      freq_hardcode = [4.28E+13,4.54E+13,4.84E+13,5.0E+13,5.17E+13,5.450723E+13,6.258446E+13,7.185862E+13,7.89E+13,8.250724E+13,8.8e+13,9.473368E+13,9.677E+13,1.087719E+14,1.248904E+14,1.38888E+14,1.646470E+14,1.8181818E+14,2.170594E+14,2.400000E+14,2.861575E+14,3.285621E+14,3.772505E+14,4.331539E+14,4.973413E+14,5.710405E+14,6.556609E+14,7.528209E+14,8.643787E+14,9.924678E+14,1.139538E+15,1.308407E+15,1.502295E+15,1.724914E+15,1.980523E+15,2.274010E+15,2.610987E+15,2.997899E+15]
   end
   30: begin            ; For the crystalline features
      lambda0  = 1d4
      lambda1  = 40.0
      lambda2  = 5.5
      lambda3  = 0.1  
      nfrnew01 = 20
      nfrnew12 = 80
      nfrnew23 = 30 
   end
   31: begin            ; For the crystalline features
      lambda0  = 1d4
      lambda1  = 75.0
      lambda2  = 8.0
      lambda3  = 0.1  
      nfrnew01 = 20
      nfrnew12 = 200
      nfrnew23 = 30 
   end
   else: begin
      print,"Do not know this frequency resolution mode"
      stop
   endelse
endcase
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
IF fresmd EQ 4 then freqnw23=freq_hardcode  ;Forcing specific frequency points to get exact jhk photometry.
freqnew  = [freqnw01,freqnw12,freqnw23]
nf       = n_elements(freqnew)
;
; If /onlynf then simply return nf
;
if keyword_set(onlynf) then begin
    if n_elements(npah) eq 0 then npah = 0
    nspec = n_elements(inputfiles) + npah
    return
endif
;
; Write the wavelength data file
;
openw,1,'wavelength_micron.inp'
printf,1,nf
printf,1,' '
for i=0,nf-1 do printf,1,1d4*cc/freqnew[i]
close,1
;
; Read the input files one by one
;
nfiles = n_elements(inputfiles)
if n_elements(pllongs) ne nfiles then begin
    print,'ERROR: Nr of elements of inputfiles not equal to '
    print,'       nr of elements of pllongs.'
    stop
endif
for ifile=0,nfiles-1 do begin
    ;;
    ;; Read the external opacity input file, which is a list of
    ;; wavelengths and opacities
    ;;
    q=findfile(inputfiles[ifile],count=count)
    if count eq 0 then begin
        print,'ERROR: Could not find file ',inputfiles[ifile]
        stop
    endif
    openr,1,inputfiles[ifile]
    iformat=0
    nlam=0
    readf,1,iformat
    readf,1,nlam
    case iformat of 
        1:begin   ;; Only absorption cross section
            data = dblarr(2,nlam)
        end
        2:begin   ;; Absorption and scattering cross section
            data = dblarr(3,nlam)
        end
        3:begin   ;; Absorption and scattering cross section
            data = dblarr(4,nlam)
        end
        else:begin
            print,'ERROR: Could not interpret iformat = ',iformat
            stop
        end
    endcase
    readf,1,data
    close,1
    case iformat of 
        1:begin   ;; Only absorption cross section
            lam  = transpose(data[0,*])
            kabs = transpose(data[1,*])
            ksca = kabs*0.d0
        end
        2:begin   ;; Absorption and scattering cross section
            lam  = transpose(data[0,*])
            kabs = transpose(data[1,*])
            ksca = transpose(data[2,*])
        end
        3:begin   ;; Absorption and scattering cross section
            lam  = transpose(data[0,*])
            kabs = transpose(data[1,*])
            ksca = transpose(data[2,*])
        end
        else:begin
            print,'ERROR: Could not interpret iformat = ',iformat
            stop
        end
    endcase
    ;;
    ;; Check monotonicity of lam
    ;;
    for ilam=1,nlam-1 do begin
        if lam[ilam] le lam[ilam-1] then begin
            print,'ERROR: Input files must have monotonically'
            print,'       increasing wavelength grid'
            stop
        endif
    endfor
    ;;
    ;; Now map this opacity onto the current frequency grid
    ;;
    kap_a = interpol(kabs,3d14/lam,freqnew)
    kap_s = interpol(ksca,3d14/lam,freqnew)
    ;;
    ;; At long wavellength, make powerlaw
    ;;
    if max(lam) lt max(3d14/freqnew) then begin
        ii  = where(freqnew le 3d14/max(lam))
        dum = (3d14/freqnew[ii])^pllongs[ifile]
        dum = dum / dum[n_elements(dum)-1]
        dum = dum * kabs[nlam-1]
        kap_a[ii] = dum
        kap_s[ii] = 0.d0
    endif
    ;;
    ;; At short wavelength, make constant opacity
    ;;
    if min(lam) gt min(3d14/freqnew) then begin
        ii  = where(freqnew ge 3d14/min(lam))
        kap_a[ii] = kabs[0]
        kap_s[ii] = ksca[0]
    endif
    ;;
    ;; Make sure it does not get zero or less
    ;;
    kap_a = kap_a > 1d-10
    kap_s = kap_s > 1d-10
    ;;
    ;; Write this to opacity file
    ;;
    file = 'dustkappa_'+strcompress(string(ifile+1,format='(I4)'),$
           /remove_all)+'.inp'
    openw,1,file
    printf,1,3
    printf,1,nf
    for inu=0,nf-1 do begin
        printf,1,lam[inu],kap_a[inu],kap_s[inu],0.d0
    endfor
    close,1
endfor
ntherm = nfiles
;
; Total number of opacities
;
nopac = ntherm 
;
; Write the dustopac.inp file
;
spawn,'rm -f dustopac.inp'
openw,1,'dustopac.inp'
printf,1,'2               Format number of this file'
printf,1,strcompress(string(nopac,format='(I4)'),$
         /remove_all)+'               Nr of dust species'
printf,1,'============================================================================'
for iopac = 1,ntherm do begin
    printf,1,'1               Way in which this dust species is read (1=kappa file)'
    printf,1,'0               0=Thermal grain'
    printf,1,strcompress(string(iopac,format='(I4)'),$
         /remove_all)+'               Extension of name of dustopac_***.inp file'
    printf,1,'----------------------------------------------------------------------------'
endfor
close,1
;
; One more tiny thing
;
nspec = nopac
;
; End...
;
end
