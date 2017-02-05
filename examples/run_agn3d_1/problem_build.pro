;========================================================================
;                   PROBLEM GENERATOR MODULE
;
; This module is meant for creating the setup for the 2-D radiative 
; transfer problem of the test cases.
;========================================================================


;------------------------------------------------------------------------
;          HELPER ROUTINE FOR MAIN PROBLEM GENERATOR
; New version: a broken power law possible.
;------------------------------------------------------------------------
function makedens,r,theta,r0,sigma0,hr0,p1,p2,plh,hrmin=hrmin,$
                  hrstore=hrstore,cdens=cdens
thmax = !pi/2.d0
nr=n_elements(r)
nt=n_elements(theta)
rho=fltarr(nr,nt)
hrstore=fltarr(nr)
if not keyword_set(hrmin) then hrmin=0.d0
for ir=0,nr-1 do begin
   if r[ir] le r0 then pl = p1 else pl = p2
   hr = hr0 * (r[ir]/r0)^plh
   hr = sqrt(hr^2+hrmin^2)
   if not keyword_set(cdens) then begin
      rho(ir,*) = 3.98942280401d-1 * sigma0 * (r[ir]/r0)^pl * $
                  exp(-0.5*((thmax-theta(*))/hr)^2 ) / (hr*r[ir])
   endif else begin
      rho(ir,*) = ((!pi/2-theta[*]) lt hr) * sigma0 * $
                  (r[ir]/r0)^pl / (2*hr*r[ir])
   endelse
   hrstore[ir] = hr
endfor
return,rho
end



;------------------------------------------------------------------------
;              PUT MASS ON RANDOM POSITIONS
;------------------------------------------------------------------------
function random_blobs_3d,r,theta,phi,r0,sigma0,hr0,p1,p2,plh,$
             nblob,relsizeblob,istyle,iseed=iseed,hrmin=hrmin,$
             masstot=masstot,sblob=sblob,tq=tq,cdens=cdens,$
             m_random=m_random,r_random=r_random,t_random=t_random,$
             p_random=p_random
;
; First make a high-resolution 2-D density distribution as if no 
; blobs were used (we use this as the probability function)
;
rin = min(r)
rout= max(r)
nra = 200
nta = 100
ra  = rin * (rout/rin)^(dindgen(nra)/(nra-1.d0))
;ta  = 0.5*!pi*dindgen(nta)/(nta-1.d0)
ta  = !pi*dindgen(nta)/(nta-1.d0)
qq  = makedens(ra,ta,r0,sigma0,hr0,p1,p2,plh,hrmin=hrmin,cdens=cdens)
;
; Compute from this the Sigmar2(R)
;
sigma=dblarr(nra)
sigmar2=dblarr(nra)
for ir=0,nra-1 do begin
   sigma[ir]   = 2*integrate(ra[ir]*ta,transpose(qq[ir,*]))
   sigmar2[ir] = sigma[ir]*ra[ir]^2
endfor
;
; Compute the total dust mass of the torus
;
masstorus = diskmass(ra,ta,qq,gastodust=1.)
;
; Make the random positions in r, theta and phi, and the corresponding masses
; The random positions are equally distributed in log(R) and cos(theta). 
; The masses are proportional to rho*r^3 (note: they will be normalized later)
;
case istyle of
   1: begin
      r_random = dblarr(nblob)
      m_random = dblarr(nblob)
      t_random = dblarr(nblob)
      p_random = dblarr(nblob)
      for iblob=0,nblob-1 do begin
         r_random[iblob] = rin * (rout/rin)^(randomu(iseed,/uniform))
         t_random[iblob] = acos(randomu(iseed,/uniform))
         p_random[iblob] = 2*!dpi*randomu(iseed,/uniform)
         hunt,ra,r_random[iblob],ir,eps_r
         hunt,ta,t_random[iblob],it,eps_t
         dum0 = (1-eps_t)*qq[ir,it]+eps_t*qq[ir,it+1]
         dum1 = (1-eps_t)*qq[ir+1,it]+eps_t*qq[ir+1,it+1]
         dumm = (1-eps_r)*dum0+eps_r*dum1
         m_random[iblob] = (r_random[iblob]^3)*dumm
         print,'mblob = ',m_random[iblob]/1.99d33,' dumm = ',dumm
      endfor
   end
   2: begin
      r_random = dblarr(nblob)
      m_random = dblarr(nblob)
      t_random = dblarr(nblob)
      p_random = dblarr(nblob)
      for iblob=0,nblob-1 do begin
         r_random[iblob] = exp(randomdistr(alog(ra),sigmar2,iseed))
         hunt,ra,r_random[iblob],ir,eps
         if ir lt 0 or ir ge nra then stop
         rh = (1-eps)*qq[ir,*]+eps*qq[ir+1,*]
         t_random[iblob] = acos(randomdistr(cos(ta),rh,iseed))
         p_random[iblob] = 2*!dpi*randomu(iseed,/uniform)
         m_random[iblob] = 1.d0
      endfor
   end
   10: begin
      r_random = dblarr(nblob)
      m_random = dblarr(nblob)
      t_random = dblarr(nblob)
      p_random = dblarr(nblob)
      for iblob=0,nblob-1 do begin
         r_random[iblob] = exp(randomdistr(alog(ra),sigma,iseed))
         hunt,ra,r_random[iblob],ir,eps
         if ir lt 0 or ir ge nra then stop
         rh = (1-eps)*qq[ir,*]+eps*qq[ir+1,*]
         t_random[iblob] = acos(randomdistr(cos(ta),rh,iseed))
         p_random[iblob] = 2*!dpi*randomu(iseed,/uniform)
         m_random[iblob] = (r_random[iblob]/ra[0])^2*sin(t_random[iblob])
      endfor
   end
endcase
;
; Now normalize
;
masstot=0.d0
for iblob=0,nblob-1 do begin
   masstot = masstot + 2*m_random[iblob]
endfor
m_random = m_random * (masstorus/masstot)
masstot  = masstorus
;
; Send a message to user
;
tq = m_random * 0.d0
if istyle eq 10 then begin
   for iblob=0,nblob-1 do begin
;;;      rhoavblob = m_random[iblob]/$
;;;           (2*!dpi^2*sblob^2*r_random[iblob]^3*sin(t_random[iblob]))
;;;      tq[iblob] = rhoavblob * sblob * r_random[iblob]
      tq[iblob] = m_random[iblob] / $
              ( (2*!dpi)^(1.5) * r_random[iblob]^2 * sblob * $
                sin(t_random[iblob]) )
   endfor
endif
;
; Now put the blobs on the grid
;
if n_elements(iseed) eq 0 then iseed=324234
nr=n_elements(r)
nt=n_elements(theta)
np=n_elements(phi)
rho=dblarr(nr,nt,np)
rr = rebin(r,nr,nt,np)
tt = transpose(rebin(theta,nt,nr,np),[1,0,2])
pp = transpose(rebin(phi,np,nt,nr),[2,1,0])
for iblob=0,nblob-1 do begin
   ;;
   ;; Get the random log(R) and cos(theta)
   ;;
   rblob = r_random[iblob]
   tblob = t_random[iblob]
   pblob = p_random[iblob]
   print,rblob,(!dpi/2.d0-tblob)/(!dpi/2),pblob/(!dpi*2.d0)
   ;;
   ;; Put the blob on a grid
   ;;
   xx      = (rr-rblob)/rblob
   yy      = (tt-tblob)
   zz      = (pp-pblob)
   rhoblob = exp(-0.5*(xx^2+yy^2+zz^2)/relsizeblob^2)
   ;;
   ;; If the blob is close to phi=0 then add it to the phi=2*pi region, too,
   ;; to ensure correct cyclic boundaries
   ;;
   if pblob lt 0.7d0 then begin
      xx      = (rr-rblob)/rblob
      yy      = (tt-tblob)
      zz      = (pp-(pblob+2*!dpi))
      rhoblob = rhoblob + exp(-0.5*(xx^2+yy^2+zz^2)/relsizeblob^2)
   endif
   ;;
   ;; If the blob is close to phi=2*pi then add it to the phi=0 region, too,
   ;; to ensure correct cyclic boundaries
   ;; 
   if pblob gt 2*!dpi-0.7d0 then begin
      xx      = (rr-rblob)/rblob
      yy      = (tt-tblob)
      zz      = (pp-(pblob-2*!dpi))
      rhoblob = rhoblob + exp(-0.5*(xx^2+yy^2+zz^2)/relsizeblob^2)
   endif
   ;;
   ;; Normalize the mass of this blob
   ;;
   m = mass(rhoblob,r,theta,phi,vol=vol,gastodust=1.)
   rhoblob=rhoblob*(m_random[iblob]/m)
   ;;
   ;; Add this blob to the total rho
   ;; 
   rho = rho + rhoblob
endfor
return,rho
end


;------------------------------------------------------------------------
;                     MAKE OPTICAL DEPTHS
;------------------------------------------------------------------------
function maketau,a,kappa
sz=size(a.rho)
nr=sz(1)
nt=sz(2)
taur=fltarr(nr,nt)
taut=fltarr(nr,nt)
for it=0,nt-1 do begin
    taur(0,it) = 0.d0
    for ir=1,nr-1 do begin
        taur(ir,it) = taur(ir-1,it) + $
            0.5*(a.rho(ir,it)+a.rho(ir-1,it))*kappa * $
            (a.r(ir)-a.r(ir-1))
    endfor
endfor
for ir=0,nr-1 do begin
    taut(ir,0) = 0.d0
    for it=1,nt-1 do begin
        taut(ir,it) = taut(ir,it-1) + $
            0.5*(a.rho(ir,it)+a.rho(ir,it-1))*kappa * $
            (a.theta(it)-a.theta(it-1))*a.r(ir)
    endfor
endfor
return,{taur:taur,taut:taut}
end



;------------------------------------------------------------------------
;                 PROBLEM MAIN FUNCTIONS
;------------------------------------------------------------------------
function problem,nr,nt,np,rin,rout,r0,hrgrid,hrdisk,plh,sigma0,p1,p2,$
    nlevr=nlevr,nspanr=nspanr,nstepr=nstepr,hrlg=hrlg,hrmin=hrmin,$
    ntextra=ntextra,hrgmax=hrgmax,nblob=nblob,sblob=sblob,istyle=istyle,$
    iseed=iseed,tq=tq,m_random=m_random,r_random=r_random,t_random=t_random,$
    cdens=cdens,nomirror=nomirror
;
; R-grid:
;
if not keyword_set(nlevr) then begin
    r = rin*(rout/rin)^(findgen(nr)/(nr-1.d0))
endif else begin
    if not keyword_set(nspanr) then nspanr=1
    if not keyword_set(nstepr) then nstepr=1
    nextra = nspanr*nlevr*(2^nstepr-1)
    if nextra gt nr-5 then begin
       print,"Sorry not enough nr for this refinement level"
       stop
    endif
    r = rin*(rout/rin)^(findgen(nr-nextra)/(nr-1.d0-nextra))
    for ilev=1,nlevr do begin
        for ispan=nspanr,1,-1 do begin
            rins=fltarr(2^nstepr-1)
            fr=(r(ispan)/r(ispan-1))^(0.5^nstepr)
            for i=1,2^nstepr-1 do rins(i-1)=r(ispan-1)*fr^i
            r=[r(0:ispan-1),rins,r(ispan:*)]
        endfor
    endfor
endelse
;
; Theta-grid:
;
if not keyword_set(hrlg) then begin
   thmax = !pi/2.d0
   thmin = !pi/2.d0 - hrgrid
   theta = (thmax-thmin)*findgen(nt)/(nt-1.d0+0.5d0)+thmin
   if keyword_set(hrgmax) then begin
      if not keyword_set(ntextra) then begin
         print,'PROBLEM: If you define hrgmax, must also define nextra'
         stop
      endif 
      if hrgmax le hrgrid then begin
         print,'PROBLEM: hrgmax must be larger than hrgrid'
         stop
      endif
      thmax = !pi/2.d0 - hrgrid
      thmin = !pi/2.d0 - hrgmax
      thex=(thmax-thmin)*findgen(ntextra)/(ntextra*1.d0)+thmin
      theta = [thex,theta]
   endif
endif else begin
   B     = alog(hrgrid)
   A     = (B-alog(hrlg))/(nt-1.d0)
   theta = !pi/2.d0-exp(-A*findgen(nt)+B)
endelse
if keyword_set(nomirror) then begin
   theta = [theta,rotate(!dpi-theta,2)]
endif
nnt = n_elements(theta)
;
; Phi-grid:
;
phi = 2.d0*!dpi*(dindgen(np)+0.5d0)/(1.d0*np)
;
; The density distribution:
;
if nblob eq 0 or sblob eq 0 then begin
   ;;
   ;; Classical (non-blobby) manner
   ;;
   rho=makedens(r,theta,r0,sigma0,hrdisk,p1,p2,plh,hrmin=hrmin,$
                hrs=hrs,cdens=cdens)
   rho=rebin(rho,nr,nnt,np)
endif else begin
   ;;
   ;; Blobby manner
   ;;
   if nblob eq 0 or sblob eq 0 then begin
      print,'ERROR: If blobby: then define nblob AND sblob'
      stop
   endif
   if n_elements(istyle) eq 0 then istyle=1
   rho=random_blobs_3d(r,theta,phi,r0,sigma0,hrdisk,p1,p2,plh,nblob,sblob,istyle,$
        iseed=iseed,sblob=sblob,tq=tq,m_random=m_random,r_random=r_random,$
        t_random=t_random,cdens=cdens)
endelse
;
; Make the sigma array
;
nr    = n_elements(r)
sigma = dblarr(nr)
for ir=0,nr-1 do begin
   if r[ir] le r0 then plsig = p1 else plsig = p2
   sigma[ir] = sigma0 * (r[ir]/r0)^plsig
endfor
;
;rr=rebin(r,nr,nnt)
;tt=transpose(rebin(theta,nnt,nr))
;
;return,{r:r,theta:theta,phi:phi,rr:rr,tt:tt,rho:rho,sigma:sigma}
return,{r:r,theta:theta,phi:phi,rho:rho,sigma:sigma}
end


;-----------------------------------------------------------------
;                         PLOT THE GRID
;-----------------------------------------------------------------
pro plotgrid,a,data=data,rmax=rmax,ymax=ymax,ax=ax,az=az,$
           contour=contour,nlevels=nlevels,_extra=_extra
if not keyword_set(rmax) then rmax=max(a.r)
if not keyword_set(data) then data=alog10(a.rho)
if not keyword_set(az) then az=0
if not keyword_set(ax) then ax=90
if not keyword_set(ymax) then ymax=rmax
if not keyword_set(nlevels) then nlevels=100
if not keyword_set(contour) then begin
   surface,data,a.rr*sin(a.tt),a.rr*cos(a.tt),ax=ax,az=az,$
     xrange=[0,rmax],yrange=[0,ymax],_extra=_extra
endif else begin 
   contour,data,a.rr*sin(a.tt),a.rr*cos(a.tt),$
     xrange=[0,rmax],yrange=[0,ymax],nlevels=nlevels,_extra=_extra
endelse
end


;-----------------------------------------------------------------
;                 FIND A FREQUENCY IN A TABLE
;-----------------------------------------------------------------
function find_freq,freq,nu
nf=n_elements(freq)
idx=-100
for inu=0,nf-2 do begin
   if (nu ge freq(inu) and nu le freq(inu+1)) or $
     (nu ge freq(inu+1) and nu le freq(inu)) then begin
      idx = inu
   endif
endfor
if idx lt 0 then begin
   print,"ERROR: nu out of range"
   stop
endif
return,idx
end

;-----------------------------------------------------------------
;                COMPUTE THE MASS OF THE DISK
;-----------------------------------------------------------------
function mass,rhodust,r,t,p,vol=vol,gastodust=gastodust
if n_elements(gastodust) eq 0 then gastodust=100.
if max(t) le !pi/2 then begin
   imirt=1 
endif else begin 
   imirt=0
;   print,'Warning: mirror theta not active'
endelse
nr=n_elements(r)
nt=n_elements(t)
np=n_elements(p)
ri=dblarr(nr+1)
ti=dblarr(nt+1)
pi=dblarr(np+1)
ri[0]=r[0]
ri[nr]=r[nr-1]
ti[0]=0.
pi[0]=0.d0
pi[np]=2*!dpi
if imirt eq 1 then ti[nt]=!dpi/2 else ti[nt]=!dpi
for ir=1,nr-1 do ri[ir]=sqrt(r[ir]*r[ir-1])
for it=1,nt-1 do ti[it]=0.5*(t[it]+t[it-1])
for ip=1,np-1 do pi[ip]=0.5*(p[ip]+p[ip-1])
if n_elements(vol) le 1 then begin
   dri = ri[1:nr]^3 - ri[0:nr-1]^3
   dti = abs(cos(ti[0:nt-1])-cos(ti[1:nt]))
   dpi = pi[1:np] - pi[0:np-1]
   vol = (1.d0/3.d0)*rebin(dri,nr,nt,np)*$
         transpose(rebin(dti,nt,nr,np),[1,0,2])*$
         transpose(rebin(dpi,np,nt,nr),[2,1,0])
endif
m = total(rhodust*vol)
m = m*gastodust
if imirt eq 1 then m=m*2
return,m
end


