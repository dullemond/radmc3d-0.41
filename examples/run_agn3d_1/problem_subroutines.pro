;------------------------------------------------------------------------
;                   FIND PEAK OF STELLAR RADIATION 
;
; In order to determine the relevant kappa for the determination of the
; dust evaporation radius we need to find the peak of the stellar 
; spectrum. 
;------------------------------------------------------------------------
function find_peak_starspec,nu,nulnu
nf=n_elements(nu)
inu=-1
val=0.d0
for i=0,nf-1 do begin
   if nulnu[i] gt val then begin
      inu=i
      val=nulnu[i]
   endif
endfor
return,inu
end


;------------------------------------------------------------------------
;              FUNCTION: FIND THE SURFACE HEIGHT
;
; Solve the equation:
;
;                            2 * alpha
;   1 - erf(chi/sqrt(2)) = -------------
;                          Sigma * kappa
;
; where Sigma is the total surface density from z=-infty..infty.
; The solver goes down to chi=1.0, below which it will simply 
; return 0.d0.
;------------------------------------------------------------------------
function findchi,kappa,alpha,sigma
pi    = !dpi
rhs   = 2.d0 * alpha / ( kappa * sigma )
chi   = 4.0
iloop = 1
if 0.5*sigma*kappa/alpha lt 3.0 then begin
   return,0.d0
endif
lim = 0.1
repeat begin
   chiold = chi
   dum   = -2.d0*alog(rhs*exp(-chi^2/2.d0)/(1.d0-errorf(chi*sqrt(0.5d0))))
   chi   = sqrt(abs(dum))
   iloop = iloop+1      
   if iloop gt 20 then begin
      print,'No convergence in Hs solver'
      stop
   endif
endrep until abs(chi-chiold) lt 1.d-4*abs(chi+chiold) 
return,chi
end


;------------------------------------------------------------------------
;                       FUNCTION: PLANCK
;------------------------------------------------------------------------
function bplanck,nnu,TT
nu=nnu*1.d0
T=TT*1.d0
cc=2.9979d10
hh=6.6262d-27
kk=1.3807e-16
n=n_elements(nu)
bpl=dblarr(n)
if TT eq 0.d0 then return,dblarr(n)
for i=0,n-1 do begin
  x=hh*nu[i]/(kk*T)
  ;print,x
  if x gt 100.0 then b=0.d0
  if x gt 1.d-3 then b=(2.d0*hh*nu[i]^3/cc^2)/(exp(x)-1.d0) $
    else b=2.0*nu[i]^2*kk*T/cc^2
  bpl[i]=b
endfor
return,bpl
end


;------------------------------------------------------------------------
;              FUNCTION: FIND THE PLANCK MEAN KAPPA 
;------------------------------------------------------------------------
function kappaplanck,opac,temp,ispec=ispec
;
; Constants
;
pi    = !dpi
ss    = 5.6703d-5
AU    = 1.496d13
GG    = 6.672d-8
Msun  = 1.99d33
Lsun  = 3.86d33
mugas = 2.3d0
mp    = 1.6726d-24
kk    = 1.3807d-16
cc    = 2.9979d10
hh    = 6.6262d-27
;
; First make a planck spectrum
;
nu    = opac.freq
bnu   = bplanck(nu,temp)
;
; Get the kappa table for this temperature
;
kappa = opac.cabs
;
; Now make the average
;
nf=n_elements(nu)
aa=0.d0
bb=0.d0
for i=1,nf-2 do begin
   dnu = 0.5 * nu[i] * (alog(nu[i+1])-alog(nu[i-1]))
   aa  = aa + bnu(i)*kappa[i]*dnu
   bb  = bb + bnu(i)*dnu
endfor
kappa = aa/bb
;
; Find the corresponding kappa
;
return,kappa
;
end


;------------------------------------------------------------------------
;               SOLVE RIN INCLUDING THE SELF-IRRADIATION
;
; This is a way to include the self-irradiation of the inner wall into
; the determination of the radius of the inner wall for a given wall 
; temperature (use solvehpinselfirr() instead if you want to give Rin
; instead of Tin). So we give the Hp0 which is calculated without these 
; self irradiation corrections, and then make the corrections. This
; involves an iteration, because as Rin moves outwards, the Sigma(Rin)
; changes (simply because Sigma(R)=Sigma0*(R/Rsig)^plsig, where Rsig
; is NOT the Rin, but a fixed radius). NOTE: The argument 'kap'
; contains already the factor of 8 due to the structure of the inner
; rim (see paper). 
;------------------------------------------------------------------------
function solverinselfirr,rin0,Sigma0,rsig,plsigma,kap,Hp0,$
        fixhs=fixhs,chi=chi,hsback=hsback,chback=chback
rin = rin0
for i=0,10 do begin
   if not keyword_set(chi) then begin  
       ;; 
       ;; Modify Sigma, according to Sigma(R) powerlaw
       ;;
       Sig = Sigma0*(rin/rsig)^plsigma
       chi = findchi(kap,1.d0,Sig)
   endif
   if not keyword_set(fixhs) then begin
       ;;
       ;; Since Tin is fixed, we can use the equation of vertical
       ;; pressure balance to compute Hpin.
       ;;
       Hpin = Hp0 * (rin/rin0)^(1.5)
       Hsin = chi * Hpin
   endif else begin
       Hsin = fixhs
       Hpin = Hsin / chi
   endelse
   ;;
   ;; Now modify Rin by the effect of self-irradiation
   ;;
   rinold = rin
   rin = rin0 * sqrt(1.d0+Hsin/Rin)
   if abs(rin-rinold) lt 1.d-3*(rin+rinold) then begin
      hsback=hsin
      chback=chi
      return,rin
   endif   
endfor
print,"MAJOR PROBLEM: solve Rin with self-irradiation!"
stop
end


;------------------------------------------------------------------------
;              STANDARD INTEGRATION ROUTINE
;------------------------------------------------------------------------
function integrate,x,f,cons=cons,prim=prim
n=n_elements(x)
if keyword_set(prim) then ff=dblarr(n)
if(n_elements(f) ne n) then begin
   print,"Error in function integrate: x and f not equally long"
   stop 
endif
sign=2.d0*(x[n-1] gt x[0])-1.d0
int=0.d0
iact=0
if keyword_set(cons) then begin
   for i=1,n-1 do begin
      int=int+f[i-1]*(x[i]-x[i-1])
      if keyword_set(prim) then ff[i]=int
   endfor
endif else begin
   for i=1,n-1 do begin
      int=int+0.5d0*(f[i]+f[i-1])*(x[i]-x[i-1])
      if keyword_set(prim) then ff[i]=int
   endfor
endelse
if not keyword_set(prim) then return,int*sign else begin
   if x[0] lt x[n-1] then return,ff else return,ff-ff[n-1]
endelse
end

;-----------------------------------------------------------------
;                COMPUTE THE MASS OF THE 2-D DISK
;-----------------------------------------------------------------
function diskmass,r,t,rhodust,gastodust=gastodust
if n_elements(gastodust) eq 0 then gastodust=100.
if max(t) le !pi/2 then begin
   imirt=1 
endif else begin 
   imirt=0
;   print,'Warning: mirror theta not active'
endelse
nr=n_elements(r)
nt=n_elements(t)
ri=dblarr(nr+1)
ti=dblarr(nt+1)
ri[0]=r[0]
ri[nr]=r[nr-1]
ti[0]=0.
if imirt eq 1 then ti[nt]=!dpi/2 else ti[nt]=!dpi
for ir=1,nr-1 do ri[ir]=sqrt(r[ir]*r[ir-1])
for it=1,nt-1 do ti[it]=0.5*(t[it]+t[it-1])
sz=size(rhodust)
if sz[0] eq 3 then begin
   na=sz[3]
endif else begin
   na=1
endelse
m=0.d0
for ia=0,na-1 do begin
   for ir=0,nr-1 do begin
      for it=0,nt-1 do begin
         vol = (2.*!dpi/3.)*(ri[ir+1]^3-ri[ir]^3)*abs(cos(ti[it])-cos(ti[it+1]))
         m   = m + vol * rhodust[ir,it,ia]
      endfor
   endfor
endfor
m=m*gastodust
if imirt eq 1 then m=m*2
return,m
end

;-------------------------------------------------------
;     SIMILAR TO NRECIP HUNT, BUT MUCH SIMPLER
;-------------------------------------------------------
pro hunt,xx,x,jlo,eps
n=n_elements(xx)
if xx[n-1] lt xx[0] then rev=-1 else rev=1
if rev eq 1 then begin
   q1=xx[1:n-1] gt x
   q2=xx[0:n-2] le x
endif else begin
   q1=xx[1:n-1] lt x
   q2=xx[0:n-2] ge x
endelse
qq=q1*q2
jlo=where(qq eq 1)
if n_elements(jlo) gt 1 then stop 
jlo=jlo[0]
if jlo lt 0 then begin
   if rev eq 1 then begin
      if x lt xx[0] then jlo=-1
      if x ge xx[n-1] then jlo=n
   endif else begin
      if x gt xx[0] then jlo=-1
      if x le xx[n-1] then jlo=n
   endelse
   eps = 0.d0
   return
endif
eps = (x-xx[jlo]) / (xx[jlo+1]-xx[jlo])
eps = eps[0]
return
end


;-------------------------------------------------------
;      RANDOM NUMBER FROM A GIVEN DISTRIBUTION
;-------------------------------------------------------
function randomdistr,x,p,iseed
nx=n_elements(x)
pp=dblarr(nx+1)
pp[0]=0.d0
pp[1]=0.5*(x[1]-x[0])*p[0]
for i=2,nx-1 do pp[i]=pp[i-1]+0.5d0*(x[i]-x[i-2])*p[i-1]
pp[nx]=pp[nx-1]+0.5*(x[nx-1]-x[nx-2])*p[nx-1]
pp=pp/pp[nx]
xx=randomu(iseed,/uniform)
hunt,pp,xx,ip,eps
if ip eq 0 then begin
   xr=(1-eps)*x[0]+eps*0.5*(x[1]+x[0])
endif else begin
   if ip eq nx-1 then begin
      xr=(1-eps)*0.5*(x[nx-2]+x[nx-1])+eps*x[nx-1]
   endif else begin
      xr=(1-eps)*0.5*(x[ip]+x[ip-1])+eps*0.5*(x[ip+1]+x[ip])
   endelse
endelse
return,xr
end
