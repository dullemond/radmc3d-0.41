;--------------------------------------------------------------------
;                         MAKE THE R-GRID
;--------------------------------------------------------------------
function make_rgrid,rin,rout,nr,rrefine=rrefine
if not keyword_set(rrefine) then begin
    r = rin*(rout/rin)^(findgen(nr)/(nr-1.d0))
endif else begin
    nextra = rrefine.nspanr*rrefine.nlevr*(2^rrefine.nstepr-1)
    if nextra gt nr-5 then begin
       print,"Sorry not enough nr for this refinement level"
       stop
    endif
    r = rin*(rout/rin)^(findgen(nr-nextra)/(nr-1.d0-nextra))
    for ilev=1,rrefine.nlevr do begin
        for ispan=rrefine.nspanr,1,-1 do begin
            rins=fltarr(2^rrefine.nstepr-1)
            fr=(r(ispan)/r(ispan-1))^(0.5^rrefine.nstepr)
            for i=1,2^rrefine.nstepr-1 do rins(i-1)=r(ispan-1)*fr^i
            r=[r(0:ispan-1),rins,r(ispan:*)]
        endfor
    endfor
 endelse
;
; Check the DeltaR/R
;
nnr  = n_elements(r)
if nnr ne nr then stop
drr  = (r[nr-1]-r[nr-2])/r[nr-2]
if drr gt 0.15 then begin
   print,'ERROR: Radial grid too coarse...'
   stop
endif
if drr lt 0.05 then begin
    print,'Radial grid too finely spaced... This will take too much'
    print,'computational time... Are you sure this is okay? (Type 1)'
    read,i
    if i ne 1 then stop
endif
return,r
end


;--------------------------------------------------------------------
;                       MAKE THE THETA-GRID
;
; Produce the theta grid for the RADMC simulation
;
; NOTE: Contrary to an earlier version (not used here) this routine
;       assumes that nt is the total numberr of theta points above the
;       equator. That is: n_elements(theta).eq.nt
; EXCEPTION: This is no longer true if zrefine is used. That is 
;       because zrefine is really necessary for a very fine midplane
;       layer of dust or so. You then don't want to lose other 
;       theta points.
;--------------------------------------------------------------------
function make_tgrid,hrgrid,nt,hrgmax=hrgmax,ntex=ntex,$
                    hrlg=hrlg,zrefine=zrefine
if keyword_set(ntex) then ntt=nt-ntex else ntt=nt
if ntt lt 4 then stop
if not keyword_set(hrlg) then begin
    thmax = !pi/2.d0
    thmin = !pi/2.d0 - hrgrid
    theta = (thmax-thmin)*findgen(ntt)/(ntt-1.d0+0.5d0)+thmin
    if keyword_set(hrgmax) then begin
        if not keyword_set(ntex) then begin
            print,'PROBLEM: If you define hrgmax, must also define nextra'
            stop
        endif 
        if hrgmax le hrgrid then begin
            print,'PROBLEM: hrgmax must be larger than hrgrid'
            stop
        endif
        thmax = !pi/2.d0 - hrgrid
        thmin = !pi/2.d0 - hrgmax
        thex=(thmax-thmin)*findgen(ntex)/(ntex*1.d0)+thmin
        theta = [thex,theta]
    endif
endif else begin
    B     = alog(hrgrid)
    A     = (B-alog(hrlg))/(nt-1.d0)
    theta = !pi/2.d0-exp(-A*findgen(nt)+B)
endelse
if keyword_set(zrefine) then begin
    hunt,theta[*],!pi/2-zrefine.zrref,iz
    theta = [theta[0:iz-1],!pi/2-rotate(dindgen(zrefine.nzref)+0.5,2)*$
             (!pi/2-theta[iz-1])/(1.d0*zrefine.nzref+1.)]
endif
return,theta
end



;--------------------------------------------------------------------
;                  ROUTINES FOR MIXING OPACITIES
;--------------------------------------------------------------------

pro mixopacities,files,fileo,abun
;
nfiles0 = n_elements(files)
nfiles  = nfiles0
for ifile=nfiles0-1,0,-1 do begin
   if files[ifile] eq '' then nfiles=ifiles
endfor
if nfiles lt 2 then return
if n_elements(abun) lt nfiles then stop
ab     = abun
dum    = total(ab)
ab     = ab / dum
;
; Read the first file (which sets the lambda grid)
;
openr,1,files[0]
iformat=0
nf=0
readf,1,iformat
case iformat of
    1: idum=2
    2: idum=3
    3: idum=4
    else: stop
endcase
readf,1,nf
data1=dblarr(idum,nf)
readf,1,data1
close,1
;
; Set the wavelength grid and the main kappa file
;
lambda     = transpose(data1[0,*])
kappa_abs  = ab[0] * transpose(data1[1,*])
if iformat gt 1 then begin
    kappa_sca  = ab[0] * transpose(data1[2,*])
endif else begin
    kappa_sca  = dblarr(nf)
endelse
;
; Read the other files
;
for ifile=1,nfiles-1 do begin
    ;;
    ;; Read opacity file
    ;; 
    openr,1,files[ifile]
    iformat=0
    nf2=0
    readf,1,iformat
    case iformat of
        1: idum=2
        2: idum=3
        3: idum=4
        else: stop
    endcase
    readf,1,nf2
    data2=dblarr(idum,nf2)
    readf,1,data2
    close,1
    ;;
    ;; Interpolate onto the first
    ;;
    kappa_abs = kappa_abs + ab[ifile] * exp(interpol(alog(data2[1,*]),$
                 alog(data2[0,*]),alog(data1[0,*])))
    if iformat gt 1 then begin
        kappa_sca  = kappa_sca + ab[ifile] * exp(interpol(alog(data2[2,*]),$
                 alog(data2[0,*]),alog(data1[0,*])))
    endif
    ;;
endfor
;
; Write the opacity mixture
;
openw,1,fileo
printf,1,3
printf,1,nf
printf,1,' '
for i=0,nf-1 do printf,1,lambda[i],kappa_abs[i],$
     kappa_sca[i],0.d0,format='(4(E13.6,1X))'
close,1
;
end


function readkappa,file
openr,1,file
iformat=0
readf,1,iformat
case iformat of
    1: begin
        nf=0
        readf,1,nf
        data = dblarr(2,nf)
        readf,1,data
        lambda    = transpose(data[0,*])
        kappa_abs = transpose(data[1,*])
        kappa_sca = dblarr(nf)
        g_sca     = dblarr(nf)
    end
    2: begin
        nf=0
        readf,1,nf
        data = dblarr(3,nf)
        readf,1,data
        lambda    = transpose(data[0,*])
        kappa_abs = transpose(data[1,*])
        kappa_sca = transpose(data[2,*])
        g_sca     = dblarr(nf)
    end
    3: begin
        nf=0
        readf,1,nf
        data = dblarr(4,nf)
        readf,1,data
        lambda    = transpose(data[0,*])
        kappa_abs = transpose(data[1,*])
        kappa_sca = transpose(data[2,*])
        g_sca     = transpose(data[3,*])
    end
endcase
close,1	
freq = 2.9979d14/lambda
return,{freq:freq,lambda:lambda,kappa_abs:kappa_abs,kappa_sca:kappa_sca,$
           g_sca:g_sca}
end




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
;                COMPUTE THE MASS OF THE DISK
;-----------------------------------------------------------------
function mass,rhodust,r,t,gastodust=gastodust
if n_elements(gastodust) eq 0 then gastodust=100.
if max(t) le !pi/2 then begin
   imirt=1 
endif else begin 
   imirt=0
   print,'Warning: mirror theta not active'
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
m=0.d0
for ir=0,nr-1 do begin
   for it=0,nt-1 do begin
      vol = (2.*!dpi/3.)*(ri[ir+1]^3-ri[ir]^3)*abs(cos(ti[it])-cos(ti[it+1]))
      m   = m + vol * rhodust[ir,it]
   endfor
endfor
m=m*gastodust
if imirt eq 1 then m=m*2
return,m
end


;-----------------------------------------------------------------
;                 READ THE DUST OPACITY FILES 
;-----------------------------------------------------------------
function readopac,nr=nr
close,1
if not keyword_set(nr) then nr=1
filename = 'dustopac_'+strcompress(string(nr),/remove_all)+'.inp'
print,"Reading ",filename
openr,1,filename
nf=0
ns=0
readf,1,nf,ns
if nf ge 1 then begin
    cabs = fltarr(nf,ns)
    csca = fltarr(nf,ns)
    dum  = 0.e0
    for k=1,nf do begin
        for is=1,ns do begin
            readf,1,dum
            cabs(k-1,is-1) = dum
        endfor
    endfor
    for k=1,nf do begin
        for is=1,ns do begin
            readf,1,dum
            csca(k-1,is-1) = dum
        endfor
    endfor
    nrtrange=1
    trange=[0.d0,0.d0]
endif else begin
    nf=ns
    ismooth=0
    nrtrange=0
    readf,1,ns,ismooth,nrtrange
    if(ismooth ne 0) then stop   ; No smoothing yet allowed.
    cabs = fltarr(nf,ns,nrtrange)
    csca = fltarr(nf,ns,nrtrange)
    dum  = 0.e0
    trange=fltarr(nrtrange+1)
    for ir=1,nrtrange do begin
        readf,1,a,b
        trange(ir)=b
        for k=1,nf do begin
            for is=1,ns do begin
                readf,1,dum
                cabs(k-1,is-1,ir-1) = dum
            endfor
        endfor
        for k=1,nf do begin
            for is=1,ns do begin
                readf,1,dum
                csca(k-1,is-1,ir-1) = dum
            endfor
        endfor
    endfor
endelse
close,1
file='frequency.inp'
dd=findfile(file,count=count)
if(count le 0) then begin
    print,"Could not find frequency.inp. Taking frequency.dat"
    file='frequency.dat'
    dd=findfile(file,count=count)
    if(count le 0) then begin
       print,"Could not find frequency.dat either"
       stop
    endif
endif
openr,1,file
nnf=0
readf,1,nnf
if nnf ne nf then begin
   print,"ERROR: frequency file has different nr of points as dustopac file"
   stop
endif
freq = fltarr(nf)
wave = fltarr(nf)
for k=1,nf do begin
    readf,1,dum
    freq(k-1) = dum
    wave(k-1) = 2.9979d14 / dum
endfor
close,1
return,{nf:nf,ns:ns,freq:freq,wave:wave,cabs:cabs,csca:csca,nrt:nrtrange,trange:trange}
end


;------------------------------------------------------------------------
;                     MAKE OPTICAL DEPTHS
;------------------------------------------------------------------------
function maketau,a,kappa
sz=size(a.rho)
nr=sz(1)
nt=sz(2)
nspec=n_elements(kappa)
if sz[0] le 2 then nspecc=1 else nspecc=sz[3]
if nspec ne nspecc then begin
    print,'ERROR: Nr of opacities must equal nr of dust components'
    stop
endif
taur=fltarr(nr,nt)
taut=fltarr(nr,nt)
for it=0,nt-1 do begin
    taur(0,it) = 0.d0
    for ir=1,nr-1 do begin
        dummy = 0.5*(a.rho(ir,it,*)+a.rho(ir-1,it,*))*kappa[*] * $
                (a.r(ir)-a.r(ir-1))
        taur(ir,it) = taur(ir-1,it) + total(dummy)
    endfor
endfor
for ir=0,nr-1 do begin
    taut(ir,0) = 0.d0
    for it=1,nt-1 do begin
        dummy = 0.5*(a.rho(ir,it,*)+a.rho(ir,it-1,*))*kappa[*] * $
                (a.theta(it)-a.theta(it-1))*a.r(ir)
        taut(ir,it) = taut(ir,it-1) + total(dummy)
    endfor
endfor
return,{taur:taur,taut:taut}
end


;------------------------------------------------------------------------
;                             SMOOTH INNER RIM 
;
; Sometimes the inner rim region is so optically thick that even grid
; refinement of any reasonable level is not able to get to an optically
; thin inner rim. In that case a certain level of smoothing is 
; required. This will be done automatically to the right level.
;------------------------------------------------------------------------
pro smooth_rim,r,theta,rhodust,kappa,sigdust=sigdust,$
                    drsm=drsm,tautol=tautol
nr    = n_elements(r)
nt    = n_elements(theta)
nspec = n_elements(rhodust[0,0,*])
if not keyword_set(drsm) then drsm=0.1
if not keyword_set(tautol) then tautol=0.5
if n_elements(sigdust) gt 1 and n_elements(sigdust) ne nr then begin
    print,'ERROR: smooth_rim: Sigdust has wrong number of elements'
    stop
endif
;
; Find the optical depth of the first cell at the midplane
;
tau   = 0.d0
for ispec=0,nspec-1 do begin
    tau = tau + 0.5 * ( rhodust[0,nt-1,ispec] + rhodust[1,nt-1,ispec] ) * $
          kappa[ispec] * ( r[1] - r[0] )
endfor
;
; If this is too high, then smooth the inner rim a bit
;
if tau gt tautol then begin
    ;;
    ;; Calculate the required reduction factor
    ;;
    rdxsm = tautol / tau
    ;;
    ;; Message
    ;;
    print,'Smooting inner rim with a factor ',rdxsm
    ;;
    ;; Make a smooth curve from 0 to 1, saturating at 1 around
    ;; R=R_0*(1+drsm), i.e. the point where the reduction should stop.
    ;; Using this smooth function we can prevent the reduction from
    ;; being abrupt.
    ;;
    r1 = r[0]*(1.d0+drsm)
    x  = (r-r[0])/(r1-r[0])
    q1 = 2.d1
    q2 = 1.d0
    z  = -alog10(exp(-x*q1)+10^(-q2))/q2 ; Make smooth curve
    z  = 1-(1-z)/(1-z[0])       ; Make sure z[0]==0
    ;;
    ;; Now reduce the density.
    ;; BUGFIX 04.04.07: Forgot that rhodust has 3 indices!
    ;;
    for ir=0,nr-1 do begin
        x = 1.d0-z[ir]
        redux=rdxsm^x
        rhodust[ir,*,*] = rhodust[ir,*,*]*redux
        if n_elements(sigdust) gt 1 then begin
            sigdust[ir,*] = sigdust[ir,*]*redux
        endif
    endfor
endif
;;
end


;------------------------------------------------------------------------
;                A simple integration routine
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


;--------------------------------------------------------------------
;                    READ THE GRID INFORMATION
;--------------------------------------------------------------------
function read_grid,dat=dat
nr=0
ntheta=0
nthsmall=0
nfr=0
imirt=0
close,1
;
; Read the radius data
;
file='radius.inp'
dd=findfile(file,count=count)
if keyword_set(dat) then count=0
if(count le 0) then begin
    file='radius.dat'
    dd=findfile(file,count=count)
    if(count le 0) then stop
endif
openr,1,file
readf,1,nr
r=dblarr(nr)
readf,1,r
close,1
;
; Read the theta data
;
file='theta.inp'
dd=findfile(file,count=count)
if keyword_set(dat) then count=0
if(count le 0) then begin
    file='theta.dat'
    dd=findfile(file,count=count)
    if(count le 0) then stop
endif
openr,1,file
readf,1,ntheta,imirt
nthsmall=ntheta
theta=dblarr(ntheta)
readf,1,theta
close,1
if imirt eq 1 then begin
   thbk=dblarr(2*ntheta)
   thbk(0:ntheta-1) = theta(*)
   thbk(ntheta:2*ntheta-1) = !pi - rotate(theta(*),2)
   theta = thbk
   ntheta = ntheta*2
endif
;
; Now read the frequencies
;
file='frequency.inp'
dd=findfile(file,count=count)
if keyword_set(dat) then count=0
if(count le 0) then begin
    file='frequency.dat'
    dd=findfile(file,count=count)
    if(count le 0) then begin
       print,'Could not find frequency.inp/.dat!!'
       nu=[1.0,1.0]
       nnu=0
    endif
endif
filename=file
str=findfile(filename)
if(str(0) eq filename) then begin
    openr,1,filename
    nnu=1
    readf,1,nnu
    nu=dblarr(nnu+1)
    dum=dblarr(nnu)
    readf,1,dum
    nu(1:*) = dum(*)
    close,1
endif else begin
    nu=[0.d0,0.d0]
    nnu=0
endelse
;
return,{nr:nr,ntheta:ntheta,nthsmall:nthsmall,nnu:nnu,r:r,theta:theta,nu:nu}
end


;-----------------------------------------------------------------
;                   READ DUST DENSITIES
;-----------------------------------------------------------------
function read_dustdens,ext=ext,dat=dat,col=col
mass = 1

grid = read_grid()
dum=0

nrspec=0
sizer=0
sizet=0
imirt=0

filebase='dustdens'
if n_elements(ext) gt 0 then begin
   filebase=filebase+'_'+strcompress(string(ext),/remove_all)
endif
if not keyword_set(dat) then begin
   filename=filebase+'.inp' 
endif else begin
   filename=filebase+'.dat' 
endelse
dd=findfile(filename,count=count)
if(count eq 0) then begin
   filename=filebase+'.dat'
   dd=findfile(filename,count=count)
   if(count eq 0) then stop
endif
openr,1,filename
readf,1,nrspec,sizer,sizet,imirt
sizetbig=sizet*(1+imirt)
dens=dblarr(sizer,sizetbig,nrspec)
idum=0
for ispec=1,nrspec do begin
    for ir=1,sizer do begin
        for it=1,sizet do begin
            dum = 0.d0
            readf,1,dum
            dens(ir-1,it-1,ispec-1) = dum
            if(imirt eq 1) then begin
                itt = 2*sizet + 1 - it
                dens(ir-1,itt-1,ispec-1) = dum
            endif
        endfor
    endfor
endfor
close,1
rr=rebin(grid.r,grid.nr,grid.ntheta)
if grid.ntheta gt 1 then $
   tt=transpose(rebin(grid.theta,grid.ntheta,grid.nr)) $ 
   else tt=0
close,1
if keyword_set(mass) then begin
   sigmadust=dblarr(grid.nr,nrspec)
   massdust=dblarr(nrspec)
   for ispec=0,nrspec-1 do begin
      qq=sigmamass(dens[*,*,ispec],grid.r,grid.theta,0)
      sigmadust[*,ispec]=qq.sigma
      massdust[ispec]=qq.m
   endfor
endif else begin
   sigmadust=-1
   massdust=-1
endelse
if keyword_set(col) then begin
   qq=columns(dens,grid.r,grid.theta)
   colr=qq.colr
   colt=qq.colt
endif else begin
   colr=-1
   colt=-1
endelse
return,{nr:grid.nr,ntheta:grid.ntheta,nnu:grid.nnu,$
        r:grid.r,theta:grid.theta,nu:grid.nu(1:*),$
        nspec:nrspec,rr:rr,tt:tt,rho:dens,$
        sigmadust:sigmadust,massdust:massdust,colr:colr,colt:colt}
end


;-----------------------------------------------------------------
;                   READ DUST TEMPERATURES 
;-----------------------------------------------------------------
function read_dusttemp,ext=ext,all=all
grid = read_grid()
niter=0
dum=0
nrsizemax=0
openr,1,'dusttemp.info'
readf,1,niter
readf,1,dum
readf,1,nrsizemax
readf,1,lastsave
close,1

iter=0
nrspec=0
sizer=0
sizet=0
imirt=0

if not keyword_set(all) then begin
   if(niter eq -2) then niter='final'
   if n_elements(ext) gt 0 then niter=ext
   filename='dusttemp_'+strcompress(string(niter),/remove_all)+'.dat'
endif else begin
   filename=['']
   for i=1,lastsave do begin
      filename=[filename,'dusttemp_'+strcompress(string(i),/remove_all)+'.dat']
   endfor
   if niter eq -2 then begin
      filename=[filename,'dusttemp_final.dat']
   endif
   filename=filename(1:*)
endelse

nn=n_elements(filename)

for in=1,nn do begin
   openr,1,filename[in-1]
   readf,1,nrspec,sizer,sizet,imirt
   sizetbig=sizet*(1+imirt)
   if in eq 1 then begin
      temp=dblarr(sizer,sizetbig,nrsizemax,nrspec,nn)
      nrgrainsize=intarr(nrspec)
   endif
   idum=0
   for ispec=1,nrspec do begin
      readf,1,idum
      nrgrainsize(ispec-1)=idum
      for isize=1,nrgrainsize(ispec-1) do begin
         for ir=1,sizer do begin
            for it=1,sizet do begin
               dum = 0.d0
               readf,1,dum
               temp(ir-1,it-1,isize-1,ispec-1,in-1) = dum
               if(imirt eq 1) then begin
                  itt = 2*sizet + 1 - it
                  temp(ir-1,itt-1,isize-1,ispec-1,in-1) = dum
               endif
            endfor
         endfor
      endfor
   endfor
   close,1
endfor
rr=rebin(grid.r,grid.nr,grid.ntheta)
if grid.ntheta gt 1 then $
   tt=transpose(rebin(grid.theta,grid.ntheta,grid.nr)) $ 
   else tt=0
close,1
;
; If the starinfo file is there, then read this
;
dd=findfile('starinfo.inp',count=count)
if(count gt 0) then begin
   openr,1,'starinfo.inp'
   readf,1,iformat
   readf,1,rstar
   readf,1,mstar
   readf,1,tstar
   close,1
   tthin = sqrt(0.5d0*rstar/rr)*tstar   
endif else begin
   tthin=0.d0
endelse
;
return,{nr:grid.nr,ntheta:grid.ntheta,nnu:grid.nnu,$
        r:grid.r,theta:grid.theta,nu:grid.nu(1:*),$
        rr:rr,tt:tt,nspec:nrspec,nsize:nrgrainsize,$
        temp:temp,tthin:tthin}
end



;-----------------------------------------------------------------
;                        READ THE SPECTRUM
;-----------------------------------------------------------------
function read_spectrum,ext=ext,file=file,dpc=dpc
fr_units=0
;
; Read possible dust information
;
z=findfile('dustopac.inp',count=count)
dust=count
if dust gt 0 then begin
    print,"Found dustopac.inp, so assuming dust spectrum"
    fr_units = -1
endif

;
; Read the total spectrum
;
if not keyword_set(file) then begin
   if n_elements(ext) eq 0 then begin
      filename='spectrum.dat'
   endif else begin
      filename='spectrum_'+strcompress(string(ext),/remove_all)+'.dat'
   endelse
endif else begin
   filename=file
endelse
openr,1,filename
nrfr=1
readf,1,nrfr
dum=dblarr(2,nrfr)
freq=dblarr(nrfr+1)
spectrum=dblarr(nrfr+1)
readf,1,dum
freq(1:nrfr) = dum(0,0:nrfr-1)
spectrum(1:nrfr) = dum(1,0:nrfr-1)
close,1
;
; If dpc is specified, then the spectrum is not the standard d=1pc
; spectrum, but in fact is a spectrum at a given distance (for instance,
; if you want to use the new option of RADICAL to include an aperture,
; then RADICAL returns the actual spectrum at d=dpc distance). In order
; not to confuse the plot_spectrum() routine one can rescale it back
; to a spectrum at 1pc by specifying dpc here in this read routine. 
;
if keyword_set(dpc) then spectrum=spectrum*dpc^2
;
; Return all...
;
a={spectrum:spectrum,nfr:nrfr,freq:freq,fr_units:fr_units,obs:0}
return,a
end


;-----------------------------------------------------------------
;                       PLOT SPECTRUM 
;
; ARGUMENTS:
;   ev              = 1 --> frequency in electronvolt (default=Hz)
;   kev             = 1 --> frequency in kiloelectronvolt (default=Hz)
;   lnu             = 1 --> L_nu (default nu*L_nu)
;   rpc             = Distance of observer in units of parsec
;                   = 0 (or undef) --> Lum = Total Lum in erg/s
;                     instead of erg / s cm^2
;   lumtot          = 1 --> total luminosity (default=flux at 1 parsec)
;
;-----------------------------------------------------------------
pro plot_spectrum,a,ev=ev,kev=kev,hz=hz,micron=micron,$
        lnu=lnu,fnu=fnu,nulnu=nulnu,nufnu=nufnu,rpc=rpc,$
        dpc=dpc,xlg=xlg,oplot=oplot,jy=jy,lsun=lsun,$
        itheta=itheta,ldat=ldat,ylin=ylin,lum=lum,$
        thick=thick,xthick=xthick,ythick=ythick,$
        charthick=charthick,charsize=charsize,_extra=_extra
if not keyword_set(itheta) then itheta=0
fr_units = a.fr_units
if keyword_set(micron) then begin
    fr_units = -1
endif
if keyword_set(hz) then begin
    fr_units = 0
endif
if keyword_set(ev) then begin
    fr_units = 1
endif
if keyword_set(kev) then begin
    fr_units = 2
endif
case fr_units of
    -1: begin
        xcoord = 2.9979d14 / a.freq(1:*)
        xtitle = '!4k!X [!4l!Xm]'
    end
    0: begin
        xcoord = a.freq(1:*)
        xtitle = '!4m!X [Hz]'
    end
    1: begin
        xcoord = 4.13568842841d-15 * a.freq(1:*)
        xtitle = '!4m!X [eV]'
    end
    2: begin
        xcoord = 4.13568842841d-18 * a.freq(1:*)
        xtitle = '!4m!X [KeV]'
    end
endcase

if keyword_set(dpc) then rpc=dpc
;
; Plot nuFnu or Fnu (same with Lnu)? And what about Fnu vs Lnu?
;
sed=1
ylum=0
if keyword_set(jy) then begin
   sed=0
endif
if keyword_set(fnu) then begin
   ylum=0
   sed=0
endif
if keyword_set(lnu) then begin
   ylum=1
   sed=0
endif
if keyword_set(nufnu) then begin
   ylum=0
   sed=1
endif
if keyword_set(nulnu) then begin
   ylum=1
   sed=1
endif
if keyword_set(lum) then ylum=1 
if keyword_set(jy) then begin
   ylum=0
   ;sed=0
endif
if keyword_set(lsun) then begin
   ylum=1
   sed=1
endif
;
; If the distance is not given, then for an observed spectrum
; we can only plot the flux, while for a computed spectrum we
; can only plot the luminosity
;
if a.obs eq 0 then begin
   if not keyword_set(rpc) and keyword_set(jy) then begin
      print,"Cannot use Jy units of theoretical spectrum if dpc not known."
      return
   endif
   if not keyword_set(rpc) and ylum eq 0 then begin
      ylum=1
      print,"Do not know distance to source. Plotting Luminosity instead."
   endif
endif else begin
   if not keyword_set(rpc) and ylum eq 1 then begin
      ylum=0
      print,"Do not know distance to source. Plotting Flux instead."
   endif
endelse
;
; Which plot to make? Lum or flux?
;
if ylum eq 0 then begin
   ;
   ; Plot spectrum as flux at a certain distance
   ;
   if a.obs eq 0 then begin
      rpc=1.d0*rpc
      distfact = 1.d0/ (rpc^2)
   endif else begin
      distfact = 1.d0
   endelse
   if not keyword_set(jy) then begin
      if sed eq 0 then begin
         lumfact=1.d0
         ytitle='F!I!4m!X!N  [erg cm!E-2!N Hz!E-1!N s!E-1!N]'
      endif else begin
         lumfact=a.freq(1:*)
         ytitle='!4m!XF!I!4m!X!N [erg cm!E-2!N s!E-1!N]'
      endelse
   endif else begin
      if sed eq 0 then begin
         lumfact=1d+23
         ytitle='F!I!4m!X!N [Jy]'
      endif else begin
         lumfact=1d+23*a.freq(1:*)
         ytitle='!4m!XF!I!4m!X!N [JyHz]'
      endelse
   endelse
endif else begin
   ;
   ; Plot spectrum as luminosity
   ;
   if a.obs eq 0 then begin
      distfact = 1.1965280793d38 ; = (1 parsec)^2 = 1.19d38 cm^2
   endif else begin
      rpc=1.d0*rpc
      distfact = rpc^2 * 1.1965280793d38
   endelse
   if sed eq 0 then begin
      lumfact=1.d0
      ytitle='L!I!4m!X!N  [erg Hz!E-1!N s!E-1!N]'
   endif else begin
      if not keyword_set(lsun) then begin
         lumfact=a.freq(1:*)
         ytitle='!4m!XL!I!4m!X!N [erg s!E-1!N]'
      endif else begin
         lumfact=a.freq(1:*)*2.5956986d-34
;         ytitle='!4m!XL!I!4m!X!N [Lsun]'
         ytitle='!4m!XL!I!4m!X!N [L!D!9n!X!N]'
      endelse
   endelse
endelse
if not keyword_set(xlg) then begin
    xlog=1
endif else begin
    xcoord=alog10(xcoord)
    xtitle = 'log '+xtitle
    xlog=0
endelse
if not keyword_set(ylin) then begin
    ylog=1
endif else begin
    ylog=0
endelse
if not keyword_set(oplot) then begin
    plot,xcoord,distfact*lumfact*a.spectrum(1:*,itheta),$
         xtitle=xtitle,ytitle=ytitle,$
         xlog=xlog,ylog=ylog,thick=thick,xthick=xthick,ythick=ythick,$
        charthick=charthick,charsize=charsize,_extra=_extra
endif else begin
    oplot,xcoord,distfact*lumfact*a.spectrum(1:*,itheta),$
         thick=thick,_extra=_extra
endelse
end




;------------------------------------------------------;
; READ THE DUST OPACITY FILES: THE MASTER FILES        ;
;------------------------------------------------------;
function readopacmaster,filename
close,1
print,"Reading ",filename
openr,1,filename
iformat=0
readf,1,iformat
nf=0
readf,1,nf
case iformat of
   1: ncol=2
   2: ncol=3
   3: ncol=4
endcase
data=dblarr(ncol,nf)
readf,1,data
close,1
lambda    = dblarr(nf)
cabs      = dblarr(nf)
csca      = dblarr(nf)
heny      = dblarr(nf)
lambda[*] = data[0,*]
cabs[*]   = data[1,*]
if iformat ge 2 then csca[*]   = data[2,*]
if iformat ge 3 then heny[*]   = data[3,*]
freq      = 2.9979d14/lambda
;
return,{nf:nf,ns:1,freq:freq,wave:lambda,cabs:cabs,csca:csca,nrt:1,trange:0.d0}
end

;-----------------------------------------------------------------
;                 FIND A FREQUENCY IN A TABLE
;-----------------------------------------------------------------
function find_freq,freq,nu,eps=eps
nf=n_elements(freq)
idx=-100
for inu=0,nf-2 do begin
   if (nu ge freq(inu) and nu le freq(inu+1)) or $
     (nu ge freq(inu+1) and nu le freq(inu)) then begin
      idx = inu
      eps = (nu-freq[inu]) / (freq[inu+1]-freq[inu])
   endif
endfor
if idx lt 0 then begin
   print,"ERROR: nu out of range"
   stop
endif
return,idx
end

;-----------------------------------------------------------------
;                 FIND OPACITY AT CERTAIN WAVELENGTH
;-----------------------------------------------------------------
function findopac,o,lambda=lambda,freq=freq,abs=abs,sca=sca
if n_elements(abs) eq 0 and n_elements(sca) eq 0 then begin
   abs=1
   sca=1
endif else begin
   if n_elements(abs) eq 0 then abs=0
   if n_elements(sca) eq 0 then sca=0
endelse
if keyword_set(lambda) then freq=2.9979e14/lambda
if not keyword_set(freq) then stop
if freq ge max(o.freq) then return,0.d0
if freq le min(o.freq) then return,0.d0
idx=find_freq(o.freq,freq)
eps=(freq-o.freq[idx])/(o.freq[idx+1]-o.freq[idx])
kappa=0.d0
if abs gt 0 then kappa = kappa + $
   (1.d0-eps)*o.cabs(idx)+eps*o.cabs(idx+1)
if sca gt 0 then kappa = kappa + $
   (1.d0-eps)*o.csca(idx)+eps*o.csca(idx+1)
return,kappa
end

