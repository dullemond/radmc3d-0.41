@problem_subroutines.pro

;------------------------------------------------------------------------
;            A SIMPLE DISK MODEL WITH PARAMETERIZED SIGMA(R)
;
; ARGUMENTS:
;   r          The radial grid of the 2-D density profile for RADMC
;   theta      The theta grid of the 2-D density profile for RADMC
;   rdisk      The outer radius of the disk
;   sigmadust0 The surface density of the DUST at r=rdisk
;   p1         The powerlaw of Sigmadust(R) for r<rdisk
;   p2         The powerlaw of Sigmadust(R) for r>rdisk (to avoid extreme
;                    cut-off of the disk, so choose p2=12 or so)
;   hr0        The height of the disk at r=rdisk
;                  NOTE: When iterated on structure this info gets lost
;   plh        The powerlaw index for the relative height Hp(R)/R
;                  NOTE: When iterated on structure this info gets lost
;
; KEYWORDS:
;   hrmin      The floor value for Hp(R)/R
;   hrpuff     The value of Hp(R)/R at R=R[0] for mimicing the
;                  puffing-up of the inner rim
;                  NOTE: Only useful when not iterating on structure      
;   rpuff      The radius at which the puffed-up rim structure goes
;                  over into the flaring disk stucture
;                  NOTE: Only useful when not iterating on structure      
;
; RETURNS:
;   <value>    The rho_dust(r,theta) 
;   hrstore    The Hp(R)/R values 
;   sigdust    The Sigma_dust(R)
;
;------------------------------------------------------------------------
function disk_model_1,r,theta,rdisk,sigmadust0,p1,p2,hr0,plh,hrmin=hrmin,$
                  hrstore=hrstore,hrpuff=hrpuff,rpuff=rpuff,sigdust=sigdust
;
; Some preparations
;
thmax = !pi/2.d0
nr      = n_elements(r)
nt      = n_elements(theta)
rhodust = dblarr(nr,nt)
sigdust = dblarr(nr)
hrstore = dblarr(nr)
if not keyword_set(hrmin) then hrmin=0.d0
;
; Now some stuff for the artificially puffed-up rim
; NOTE: This is not necessary when iterating on the vertical structure
;
if not keyword_set(rpuff) then begin
   rpuff=0.d0
   hrpuff0=0.d0
endif else begin
   hrpuff0=hr0 * (rpuff/rdisk)^plh
   if not keyword_set(hrpuff) then begin
      print,'ERROR: If you specify rpuff, you must also specify hrpuff'
      stop
   endif
endelse
;
; Produce the vertical pressure scale height
;
iwarn_negpuff=0
for ir=0,nr-1 do begin
   hr = hr0 * (r[ir]/rdisk)^plh
   hr = sqrt(hr^2+hrmin^2)
   if r[ir] lt rpuff then begin
      eps  = (alog(r[ir])-alog(r[0]))/(alog(rpuff)-alog(r[0]))
      hr00 = (1-eps)*hrpuff + eps*hrpuff0 
      if hr00 lt hr then iwarn_negpuff=1
      hr   = hr00
   endif
   hrstore[ir] = hr
endfor
if iwarn_negpuff ne 0 then begin
   print,'Warning: The parameterized puffed-up inner rim has thickness'
   print,'         smaller than the normal flaring part. Are you sure'
   print,'         that you want this?'
endif
;
; Produce the 2-D density distribution
;
for ir=0,nr-1 do begin
   if r[ir] le rdisk then pl = p1 else pl = p2
   hr = hrstore[ir]
   sigdust[ir]   = sigmadust0 * (r[ir]/rdisk)^pl
   rhodust[ir,*] = 3.98942280401d-1 * sigdust[ir] * $
                 exp(-0.5*((thmax-theta[*])/hr)^2 ) / (hr*r[ir])
endfor
return,rhodust
end



;------------------------------------------------------------------------
;              A SIMPLE DISK MODEL FOR GIVEN SIGMA_i(R)
;
; This routine produces a simple 2-D disk density structure (i.e.
; in gram/cm^3 in dust!!) for a given surface density profile as a
; function of radius. This profile can be on its own grid: interpolation
; in R will be done to get it on the 2-D computational grid. Also this
; routine allows multiple dust species with each its own Sigma. One
; can, with this routine and an appropriate Sigma(R), reproduce the
; same result as with disk_model_1(), but disk_model_2() is much more
; versatile.
; 
; ARGUMENTS:
;   r          The radial grid of the 2-D density profile for RADMC
;   theta      The theta grid of the 2-D density profile for RADMC
;   sigmadustdisk  The surface density of the dust species on their own rc-grid
;   rcyl       The cylindric radius-grid of the disk, i.e. sigmadustdisk(rcyl)
;   hpdisk     [optional; see below] The pressure scale height array of
;                  the disk:  hpdisk(rcyl). If not set, then the hr0, plh
;                  and hrmin are being used to internally produce hr.
;
; KEYWORDS:
;   r0         The radius where hr0 is fixed
;                  NOTE: Only if hpdisk is not set
;   hr0        The height of the disk at r=r0
;                  NOTE: Only if hpdisk is not set
;                  NOTE: When iterated on structure this info gets lost
;   plh        The powerlaw index for the relative height Hp(R)/R
;                  NOTE: Only if hpdisk is not set
;                  NOTE: When iterated on structure this info gets lost
;   hrmin      The floor value for Hp(R)/R
;                  NOTE: Only if hpdisk is not set
;   hrpuff     The value of Hp(R)/R at R=R[0] for mimicing the
;                  puffing-up of the inner rim
;                  NOTE: Only useful when not iterating on structure      
;   rpuff      The radius at which the puffed-up rim structure goes
;                  over into the flaring disk stucture
;                  NOTE: Only useful when not iterating on structure      
;
; RETURNS:
;   <value>    The rho_dust(r,theta) 
;   hrstore    The Hp(R)/R values 
;   sigmadust  The sigma of the dust on the new (RADMC) grid
;
;------------------------------------------------------------------------
function disk_model_2,r,theta,rcyl,sigmadustdisk,hpdisk,$
                  r0=r0,hr0=hr0,plh=plh,hrmin=hrmin,$
                  hrstore=hrstore,hrpuff=hrpuff,rpuff=rpuff,$
                  sigmadust=sigmadust
;
; Some preparations
;
thmax   = !pi/2.d0
nr      = n_elements(r)
nt      = n_elements(theta)
nspec   = n_elements(sigmadustdisk[0,0,*])
hrstore = dblarr(nr)
if not keyword_set(hrmin) then hrmin=0.d0
;
; Now some stuff for the artificially puffed-up rim
; NOTE: This is not necessary when iterating on the vertical structure
;
if not keyword_set(rpuff) then begin
   rpuff=0.d0
   hrpuff0=0.d0
   if not keyword_set(hrpuff) then begin
      print,'ERROR: If you specify hrpuff, you must also specify rpuff'
      stop
   endif
endif else begin
   hrpuff0=hr0 * (rpuff/r0)^plh
   if not keyword_set(hrpuff) then begin
      print,'ERROR: If you specify rpuff, you must also specify hrpuff'
      stop
   endif
endelse
;
; Defaults
;
irdisk = nr-1
ir00   = 0
;
; Interpolate the sigmadustdisk of the disk to the new grid
;
sigmadust  = dblarr(nr,nspec)
for ispec=0,nspec-1 do begin
    sigmadust  = interpol(sigmadustdisk,rcyl,r)
    ii       = where(r gt max(rcyl))
    if ii[0] ge 0 then begin
        sigmadust[ii,ispec] = 0.d0
        irdisk = min(ii)-1
    endif else begin
        irdisk = nr-1
    endelse
    ii       = where(r lt min(rcyl*0.99)) ; The 0.99 is to avoid num problems
    if ii[0] ge 0 then begin
        sigmadust[ii,ispec] = 0.d0
        ir00 = max(ii)+1
    endif else begin
        ir00 = 0
    endelse
endfor
;
; If hp is defined, then also put that onto the new grid, else use
; the usual parameters to do that.
;
if keyword_set(hpdisk) then begin
    if n_elements(hpdisk) ne n_elements(rcyl) then stop
    hp       = interpol(hpdisk,rcyl,r)
    ii       = where(r gt max(rcyl))
    if ii[0] ge 0 then begin
        hp[ii] = 0.d0
    endif
    ii       = where(r lt min(rcyl*0.99))
    if ii[0] ge 0 then begin
        hp[ii] = 0.d0
    endif
endif else begin
    for ir=0,nr-1 do begin
        hr = hr0 * (r[ir]/r0)^plh
        hr = sqrt(hr^2+hrmin^2)
        if r[ir] lt rpuff then begin
            eps = (alog(r[ir])-alog(r[0]))/(alog(rpuff)-alog(r[0]))
            hr  = (1-eps)*hrpuff + eps*hrpuff0 
        endif
        hp[ir] = hr * r[ir]
    endfor
endelse
hrstore = hp
;
; Now create the density structure
;
rhodust = dblarr(nr,nt,nspec)
for ispec=0,nspec-1 do begin
    for ir=ir00,irdisk do begin
        rhodust[ir,*,ispec] = 3.98942280401d-1 * sigmadust[ir,ispec] * $
                     exp(-0.5*(r[ir]*(!pi/2-theta[*])/hp[ir])^2 ) / hp[ir]
    endfor
endfor
;
; Done
;
return,rhodust
end


