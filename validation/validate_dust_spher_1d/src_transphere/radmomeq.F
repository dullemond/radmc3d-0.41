c     ==============================================================
c             MODULE FOR INTEGRATING THE MOMENT EQUATIONS 
c                     IN SPHERICAL COORDINATES
c
c     This module integrates the moment equations of radiative 
c     transfer in spherical coordinates, using flux-mean opacities.
c     In spherical transfer problems we are not faced with the
c     zero-flux problems encountered in passive disks (see paper
c     Dullemond, Van Zadelhoff & Natta 2002), so we can use the
c     flux-averaged opacity without problems.
c
c     ==============================================================

#include "main.h"


c     --------------------------------------------------------------
c                 DO INTEGRATION OF MOMENT EQUATIONS
c
c     This routine has the stellar luminosity as boundary condition.
c     For passive clouds there is no energy production within the
c     cloud, so that the luminosity is constant. But it is possible
c     to include energy sources or sinks. 
c
c     ARGUMENTS:
c       r          The R-grid
c       rho        The density
c       qfpi       The source term per unit gram per sterrad
c       l0         Boundary condition ir=1: the luminosity there
c       hj1        Boundary condition ir=nr: the ratio H/J there
c       fj         The eddington factor (f=0 means isotropic)
c       kaph       The H-mean opacity
c
c     OUTPUT:
c       jj         The mean intensity
c       hh         The Eddington flux
c     
c     NOTE: The f is chosen such that K = ( 1/3 + f ) J. 
c     --------------------------------------------------------------
      subroutine do_integr_momeq(nr,r,rho,qfpi,l0,hj1,fj,kaph,
     %                           jj,hh,hhrt)
      implicit none
      integer nr
      doubleprecision r(FRSIZE_R),qfpi(FRSIZE_R),l0,hj1
      doubleprecision rho(FRSIZE_R),fj(FRSIZE_R)
      doubleprecision kaph(FRSIZE_R)
      doubleprecision jj(FRSIZE_R),hh(FRSIZE_R),hhrt(FRSIZE_R)
c
      doubleprecision ll(FRSIZE_R),rk(FRSIZE_R),dum(FRSIZE_R),fpisq
      doubleprecision pi,rkorig
      parameter(pi  = 3.1415926535897932385d0)
      integer ir,iter
c
c     Constant
c
      fpisq = (4.*pi)**2
c
c     The inner boundary condition:  L(R_0) = L
c
      ll(1) = l0
      hh(1) = ll(1) / ( fpisq*r(1)**2 )
c
c     Integrate L(z) from R_0 to R_out:
c
c      dL(R)
c      ----- = (4*pi)^2 R^2 q(R)
c       dR
c
      do ir=2,nr
          ll(ir) = ll(ir-1) + (4.*pi)**2 * 0.5d0 *
     %                 ( rho(ir) * qfpi(ir) * r(ir) + 
     %                   rho(ir-1) * qfpi(ir-1) * r(ir-1) ) *
     %                 ( r(ir) - r(ir-1) ) 
          hh(ir) = ll(ir) / ( fpisq*r(ir)**2 )
          if(hh(ir).lt.0.d0) stop 633
      enddo
c
c     The outer boundary condition: J = H / hj1
c     
      jj(nr)  = abs(hh(nr)) / hj1
c
c     Now compute R^2 K = ( 1/3 + f ) R^2 J here
c
      rk(nr)  = ( (1.d0/3.d0) + fj(nr) ) * jj(nr) * r(nr)**2
c
c     The first (outermost) value of the integrand
c
      dum(nr) = rk(nr) * (2.-3.*fj(nr)) / 
     %               ( (1.+3.*fj(nr)) * r(nr) ) -
     %              rho(nr) * kaph(nr) * 
     %              hh(nr) * r(nr)**2
c     
c     Now integrate K R^2 back to the inner boundary
c
c       d (R^2 K)   (2-3f) R^2 K                 L(R)
c       --------- - ------ ----- = - rho kappa -------
c          dR       (1+3f)   R                 (4pi)^2
c
c     (Sorry for the stupid "implicit integration". I did it since it
c     is very simple, and this part of the program is not time critical.)
c
      do ir=nr-1,1,-1
          rkorig  = 0.d0
          rk(ir)  = rk(ir+1) - dum(ir+1) * ( r(ir+1) - r(ir) )
          do iter=1,20 
              dum(ir) = rk(ir) * (2.-3.*fj(ir)) / 
     %                  ( (1.+3.*fj(ir)) * r(ir) ) -
     %                  rho(ir) * kaph(ir) * 
     %                  hh(ir) * r(ir)**2
              rk(ir)  = rk(ir+1) - 0.5d0 * 
     %                  ( dum(ir) + dum(ir+1) ) *
     %                  ( r(ir+1) - r(ir) )
              if(abs(rk(ir)/rkorig-1.d0).lt.1.d-4) goto 30
              rkorig = rk(ir)
          enddo
          write(*,*) 'WEIRD ERROR: No convergence in integr iteration.'
          write(*,*) ir,rkorig,rk(ir),rk(ir+1),dum(ir),dum(ir+1)
          write(*,*) rho(ir),kaph(ir),ll(ir)
          stop 13
 30       continue
          jj(ir)  = rk(ir) / 
     %        ( ( (1.d0/3.d0) + fj(ir) ) * r(ir)**2 )
      enddo
c     
      end







c     --------------------------------------------------------------
c           COMPUTE THE DUST TEMPERATURE FROM MEAN INTENSITY
c
c     This routine will solve for all the temperatures, given the
c     mean intensity, the J-mean opacity, and an external source
c     qfpi.
c
c     ARGUMENTS:
c        qfpi         Source term in erg/s/g/ster
c        jj           Mean intensity in erg/s/cm^3/ster
c        kapj         J-mean opacity
c
c     OUTPUT:
c        temp         Dust temperature
c
c     --------------------------------------------------------------
      subroutine compute_dust_temp_from_jj(nr,qfpi,jj,kapj,temp)
      implicit none
      integer nr
      doubleprecision qfpi(FRSIZE_R),jj(FRSIZE_R),kapj(FRSIZE_R)
      doubleprecision temp(FRSIZE_R),kappa,ti,tiold,find_planckopac
c
      integer ir,i
c
c     Loop over space
c
      do ir=1,nr
c
c         Make a first guess for the temperature
c
          ti = 15.3421d0 * ( jj(ir) + qfpi(ir)/kapj(ir) )**0.25
c
c         Now do a convergence loop for the temperature
c
          do i=1,60
              kappa = find_planckopac(ti)
              tiold = ti
              ti    = 15.3421d0 * ( ( jj(ir)*kapj(ir) + qfpi(ir) ) 
     %                  / kappa )**0.25
              if(abs(ti-tiold)/(ti+tiold).lt.TEMPACCUR) goto 10
          enddo
          write(*,*) 'No convergence in temperature'
          write(*,*) jj(ir),kapj(ir)
          stop 733
 10       continue
          temp(ir) = ti
c
      enddo
c
      end


c     --------------------------------------------------------------
c                     COMPUTE THE PLANCK MEAN OPACITY
c     --------------------------------------------------------------
      subroutine compute_planckopac(nf,nt,nu,t0,t1,kappa)
      implicit none
      integer nt,nf
      doubleprecision t0,t1,nu(FRSIZE_FREQ)
      doubleprecision kappa(FRSIZE_FREQ)
c
#     include "common_transfer.h"
c
      doubleprecision enum,denom,bnucur,bnuold,bplanck
      integer it,inu
c
c     Do a loop over the discretized temperature
c
      dopac_nt = nt
      if(dopac_nt.gt.NR_TEMPERATURE) then
          write(*,*) 'Too many temperatures in requested for planck'
          write(*,*) ' Increase NR_TEMPERATURE and recompile'
          stop 13
      endif
      do it=1,dopac_nt
          dopac_temp(it)  = t0 * (t1/t0)**((it-1.d0)/(nt-1.d0))
          bnucur          = bplanck(dopac_temp(it),nu(1))
          enum            = 0.d0
          denom           = 0.d0
          do inu=2,nf
              bnuold  = bnucur
              bnucur  = bplanck(dopac_temp(it),nu(inu))
              enum    = enum + 0.5 * ( kappa(inu) * bnucur +
     %                    kappa(inu-1) * bnuold ) *
     %                    abs( nu(inu) - nu(inu-1) )
              denom   = denom + 0.5 * ( bnucur + bnuold ) *
     %                    abs( nu(inu) - nu(inu-1) )
          enddo
          dopac_kappa(it) = enum / denom
      enddo
c
      end


c     --------------------------------------------------------------
c                  INTERPOLATE THE PLANCK MEAN OPACITY
c     --------------------------------------------------------------
      function find_planckopac(temp)
      implicit none
      doubleprecision find_planckopac,temp
c
#     include "common_transfer.h"
c
      doubleprecision dum1,dum2,dt
      integer it
c
      if(dopac_nt.le.1) then 
          write(*,*) 'Must prepare planck opacity table first'
          write(*,*) dopac_nt
          stop 13
      endif
c
      if(temp.lt.dopac_temp(1)) then 
          find_planckopac = dopac_kappa(1)
      elseif(temp.gt.dopac_temp(dopac_nt)) then
c          find_planckopac = dopac_kappa(dopac_nt)
          write(*,*) 'Temperature exceeds ',dopac_temp(dopac_nt)
          write(*,*) 'If you want this to work, then increase'
          write(*,*) 'the upper temperature of planck opac table'
          stop 13
      else
          call hunt(dopac_temp,dopac_nt,temp,it)
          dt = (temp-dopac_temp(it))/
     %            (dopac_temp(it+1)-dopac_temp(it))
          if((dt.lt.0.d0).or.(dt.gt.1.d0)) stop 552
          find_planckopac = (1.d0-dt)*dopac_kappa(it)+
     %       dt*dopac_kappa(it+1)
      endif
c
      return
      end





