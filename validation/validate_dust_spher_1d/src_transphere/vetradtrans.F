c     ==============================================================
c        MODULE FOR VARIABLE EDDINGTON FACTOR RADIATIVE TRANSFER
c                       IN SPHERICAL COORDINATES
c
c     This module performs the radiative transfer using a variable
c     Eddington factor and mean opacity approach. It uses the 
c     transfer.F module for the full freq-mu-dependent transfer 
c     integral, and the radmomeq.F module for solving the moment
c     equation. In this module the transfer.F and radmomeq.F are
c     alternately called until convergence in the temperatures is
c     reached.
c
c     Note: This module assumes no scattering!
c     ==============================================================

#include "main.h"

#define HFLUXMIN 1.d-8

c     --------------------------------------------------------------
c            MAIN ROUTINE FOR THE VARIABLE EDD TENSOR TRANSFER
c
c     This routine performs and iterates the full frequency and
c     angle dependent radiative transfer equation until a convergent
c     temperature profile is reached. The algorithm uses the full 
c     transfer equation to compute all the average quantities used
c     for a moment-equation-integrator. By alternatively calling
c     the moment equation subroutine and the formal transfer 
c     subroutine, one has the convergence speed of the moment eqs,
c     but the accuracy of full RT.
c     
c     ARGUMENTS:
c        nr       = Nr of radial grid points
c        nf       = Nr of frequencies
c        ncst     = Nr of rays hitting the star
c        ncex     = Nr of rays between star and Rin
c        ncnr     = Nr of rays per radial grid point
c        r        = R-grid
c        nu       = frequency-grid
c        rho      = Density profile
c        rstar    = Stellar radius
c        starint  = Intensity at stellar surface
c        itermax  = Maximum nr of iterations
c        convcrit = Convergence criterion (dimless nr, say 1.d-3)
c        kappa    = The full (temperature INdependent) opacity
c        ing      = 0=no Ng acceleration, 1=Ng acceleration
c        itypemw  = Type of mu-weighting for mu-integrals
c        idump    = 0=no dump, >=1 = save conv hist every idump steps
c        idint    = 0=no dump, >=1 = dump full intensity array at the end
c
c     OUTPUT:
c        ierror   = Error signal (0=no error, 1=no convergence)
c        temp     = The gas/dust temperature 
c        jjme     = The freq-int mean intensity from moment equation
c        hhme     = The freq-int eddington flux from moment equation
c        jj       = The freq-integrated mean intensity
c        hh       = The freq-integrated eddington flux H = F / 4 pi
c        jnu      = Mean intensity J_nu 
c        Hnu      = Eddington flux H_nu = F_nu / 4 pi
c        fj       = The freq-integrated eddington factor f = K / J
c        fnu      = The eddington factor f_nu = K_nu / J_nu
c        kapt     = The Planck opacity for temperature temp
c        kapj     = The J-mean opacity
c        kaph     = The H-mean opacity
c     --------------------------------------------------------------
      subroutine do_vet_transfer(nr,nf,ncst,ncex,ncnr,r,nu,rho,qfpi,
     %              rstar,starint,itermax,convcrit,ing,intorder,
     %              itypemw,idump,
     %              idint,init,ierror,kappa,temp,jjme,hhme,
     %              jj,hh,jnu,hnu,fj,fnu,kapt,kapj,kaph)
      implicit none
c
#include "common_transfer.h"      
c
      integer nr,nf,ncst,ncex,ncnr,idump,ing,intorder,init,idint
      doubleprecision r(FRSIZE_R),nu(FRSIZE_FREQ),rstar
      doubleprecision starint(FRSIZE_FREQ)
      doubleprecision kappa(FRSIZE_FREQ),temp(FRSIZE_R)
      doubleprecision jj(FRSIZE_R),hh(FRSIZE_R),qfpi(FRSIZE_R)
      doubleprecision jjme(FRSIZE_R),hhme(FRSIZE_R)
      doubleprecision jnu(FRSIZE_FREQ,FRSIZE_R)
      doubleprecision hnu(FRSIZE_FREQ,FRSIZE_R)
      doubleprecision fj(FRSIZE_R),rho(FRSIZE_R)
      doubleprecision fnu(FRSIZE_FREQ,FRSIZE_R)
      doubleprecision kapt(FRSIZE_R),kapj(FRSIZE_R),kaph(FRSIZE_R)
      doubleprecision convcrit,ierror
      doubleprecision error,tempold(FRSIZE_R)
      integer itermax
c
      doubleprecision l0,hj1,find_planckopac,bplanck
      doubleprecision alp(0:FRSIZE_R)
      doubleprecision src(0:FRSIZE_R)
      doubleprecision bcout,str
      doubleprecision t1(FRSIZE_R)
      doubleprecision t2(FRSIZE_R)
      doubleprecision t3(FRSIZE_R)
      integer iter,ir,inu,imu,itdmp,savenr,icnt_ng,itypemw
      doubleprecision pi
      parameter(pi=3.1415926535897932385d0)
c
c     Reset error signal
c
      ierror = 0
c
c     Reset dump counter
c
      itdmp  = 0
      savenr = 1
c     
c     Reset the Ng counter. The Ng counter is meant to count to 4
c     and then reset, triggering the Ng accelleration at each 4-th
c     iteration step.
c
      icnt_ng = 0
c
c     Set some parameters
c
      iter_qdr = intorder
c
c     Put the Ng stuff to zero for safety
c
      do ir=1,nr
          t3(ir) = 0.d0
          t2(ir) = 0.d0
          t1(ir) = 0.d0
      enddo
c
c     Set inner boundary condition
c
      l0  = 0.d0
      do inu=2,nf
c >>> Bugfix 9-oct-2002 <<<
          l0 = l0 + 0.5 * ( starint(inu) + starint(inu-1) ) * 
     %               abs( nu(inu) - nu(inu-1) )
      enddo
      l0  = l0 * pi * ( 4.*pi*rstar**2 )
c
c     Set the outer boundary condition (initial guess)
c
      hj1 = 1.d0 
c
c     If dump is active, then clear the convergence history dump file
c
      if(idump.gt.0) then
          open(unit=1,file='convhist.info',status='unknown')
          write(1,*) 0,1
          close(1)
          open(unit=1,file='convhist.dat',status='unknown')
          write(1,*) nr
          close(1)
      endif
c
c     Put the variable eddington factor to 2/3 as a first guess, since
c     2/3 is free streaming (note: K=(1/3+f)J).
c
      do ir=1,nr
          fj(ir) = 2.d0/3.d0
      enddo
c
c     Compute the planck mean opacity, and put them as initial guesses 
c     for kapj and kaph
c
      do ir=1,nr
          kapt(ir) = find_planckopac(max(temp(ir),1.d0))
          kapj(ir) = kapt(ir)  ! Initial guess
          kaph(ir) = kapt(ir)  ! Initial guess
      enddo
c
c     Put the 'old' temperature to zero, so that the iteration will
c     take at least two steps
c
      do ir=1,nr
          tempold(ir) = 0.d0
      enddo
c
c     Reset the jnu and hnu
c
      do ir=1,nr
          hh(ir) = 0.d0
          do inu=1,nf
              jnu(inu,ir) = 0.d0
              hnu(inu,ir) = 0.d0
          enddo
      enddo
c
c     Set up the rays for the formal transfer integral
c
      call create_rays(ncst,ncex,ncnr,nr,r,rstar,itypemw)
c
c     Set up the image dimensions (for image and spectrum output)
c
      call setup_image_and_spectrum(rstar)
c
c     Do the main iteration
c
      do iter=1,itermax
c
c         Message
c     
          write(*,*) 'Iteration ',iter
c
c         If initial guess, then jump right to the transfer at iter.eq.1
c
          if((init.ne.0).and.(iter.eq.1)) goto 300
c
c         Increase the Ng counter
c
          icnt_ng = icnt_ng + 1
c
c         Integrate the moment equations
c          
          call do_integr_momeq(nr,r,rho,qfpi,l0,hj1,fj,
     %                         kaph,jjme,hhme,hh)
c
c         Compute the temperatures from these equations
c
          call compute_dust_temp_from_jj(nr,qfpi,jjme,kapj,temp)
c
c         Provide for the Ng acceleration
c
          if(ing.ne.0) then
              if(icnt_ng.eq.4) then
                  write(*,*) '   Do the Ng accelleration'
                  call do_ng_accell(nr,temp,t1,t2,t3)
                  icnt_ng = 0
              endif
c    
c             Now memorize the temperature for Ng-accelleration 
c             later. This is because Ng accelleration uses the memory
c             of the previous steps to make its prediction
c             
              do ir=1,nr
                  t3(ir) = t2(ir)
                  t2(ir) = t1(ir)
                  t1(ir) = temp(ir)
              enddo
          endif
c
 300      continue
c
c         Do transfer for each frequency and compute the moments
c
          do inu=1,nf
c
c             Compute the source terms and extinction terms for the
c             full transfer module
c
              do ir=1,nr
                  alp(ir) = rho(ir) * kappa(inu)
                  src(ir) = alp(ir) * bplanck(temp(ir),nu(inu))
              enddo
c
c             Set the boundary conditions for the full radiative
c             transfer: no incoming radiation
c
              bcout = 0.d0
c
c             Now perform the full frequency-mu-dependent transfer
c
              str = starint(inu)
              call do_transfer(nr,r,str,bcout,src,alp)
c
c             Compute all the moments
c
              call compute_moments(nr,inu,r,jnu,hnu,fnu)
c
c             Make image and spectrum
c
              call do_image_and_spectrum(inu)
c
          enddo
c
c         Now compute all the frequency averaged VET factors
c
          call compute_jh_meanopac_and_vet(nr,nf,r,nu,kappa,temp,
     %              jnu,hnu,fnu,1,jj,hh,fj,kapt,kapj,kaph)
c
c         Compute the H/J ratio at z=zmax
c          
          hj1 = abs(hh(nr)/jj(nr))
c
c         Safety check
c
          if(hj1.lt.1.d-4) stop 555
c
c         If its time, then dump the conv hist
c
          itdmp = itdmp + 1
          if((itdmp.ge.idump).and.(idump.gt.0)) then
              call append_vet_convhist(nr,temp,jjme,hhme,jj,hh,
     %            kapt,kapj,kaph,fj,savenr)
              savenr = savenr+1
              itdmp = 0
          endif
c
c         Now check the error
c
          error = 0.d0
          do ir=1,nr
              error=max(error,abs(tempold(ir)-temp(ir))/
     %                           (tempold(ir)+temp(ir)))
          enddo
          write(*,*) '     Error = ',error
          if(error.lt.convcrit) goto 10
          do ir=1,nr
              tempold(ir)=temp(ir)
          enddo
      enddo
c
c     Not converged...
c
      ierror = 1
      return
 10   continue
c
c     Converged!!!
c
      return
c
      end


c     --------------------------------------------------------------
c                   SETUP OF THE IMAGES AND SPECTRA
c     --------------------------------------------------------------
      subroutine setup_image_and_spectrum(rstar)
      implicit none
      doubleprecision rstar
c
#include "common_transfer.h"
c
      integer ib,inu
c
      if(rays_nr.eq.0) then
          write(*,*) 'INTERNAL ERROR: NEED TO SETUP RAYS FIRST'
          stop 13
      endif
c
c     First the centers of the image pixels
c
      do ib=1,rays_nr
          spim_b(ib) = rays_b(ib)
      enddo
c
c     Then the pixel boundaries
c
      spim_bi(0) = 0.d0
      do ib=1,rays_ncst-1
          spim_bi(ib) = 0.5 * ( spim_b(ib) + spim_b(ib+1) )
      enddo
      spim_bi(rays_ncst) = rstar
      do ib=rays_ncst+1,rays_nr-1
          spim_bi(ib) = sqrt(spim_b(ib)*spim_b(ib+1))
      enddo
      spim_bi(rays_nr) = spim_b(rays_nr)**2 / spim_b(rays_nr-1)
c
c     Check the pixels
c
      do ib=1,rays_nr
          if(spim_bi(ib).le.spim_bi(ib-1)) then
              stop 23473
          endif
      enddo
c
c     Finally, copy the array size
c
      spim_nb = rays_nr
c
c     Put spectrum to -1 to see later if no errors occured
c
      do inu=1,FRSIZE_FREQ
          spim_spectrum(inu) = -1.d0
      enddo
c     
      end



c     --------------------------------------------------------------
c                     COMPUTE IMAGE AND SPECTRUM
c
c     This subroutine takes the intensities found at the end-points
c     of the rays and makes these into spherical images. These are
c     then integrated to obtain the flux for the spectrum.
c     --------------------------------------------------------------
      subroutine do_image_and_spectrum(inu)
      implicit none
      integer inu
c
#include "common_transfer.h"
c
      doubleprecision flux
      integer ib
c
c     Check...
c
      if(spim_b(2).eq.0.d0) then
          write(*,*) 'INTERNAL ERROR: SPECTRUM GRID HAS TO BE SET'
          stop 13
      endif
c
c     Make image at this frequency
c
      do ib=1,spim_nb
          spim_intens(ib,inu) = rays_intensity(rays_ns(ib),ib)
      enddo
c
c     Integrate this to obtain the flux as seen at 1 parsec
c
      flux = 0.d0
      do ib=1,spim_nb
          flux = flux + ( spim_bi(ib)**2 - spim_bi(ib-1)**2 ) *
     %                   spim_intens(ib,inu)
      enddo
c
c     The norm:  pi / pc^2 = 3.29941422082d-37
c
      flux = flux * 3.29941422082d-37
c
c     Put into spectrum
c
      spim_spectrum(inu) = flux
c
      end



c     --------------------------------------------------------------
c            COMPUTE THE MOMENTS FROM THE INTENSITIES
c
c     This subroutine computes from the direct radiation field
c     I(mu,R) the moments J_nu, H_nu and f_nu. One has to call
c     this subroutine after each transfer call for each frequency.
c
c     ARGUMENTS:
c        r        = R-grid
c        inu      = frequency-grid index in which to put the values
c
c     OUTPUT:
c        jnu      = Mean intensity J_nu 
c        Hnu      = Eddington flux H_nu = F_nu / 4 pi
c        fnu      = The eddington factor f_nu = K_nu / J_nu
c 
c     --------------------------------------------------------------
      subroutine compute_moments(nr,inu,r,jnu,hnu,fnu)
      implicit none
      integer nr,inu
      doubleprecision r(FRSIZE_R)
      doubleprecision jnu(FRSIZE_FREQ,FRSIZE_R)
      doubleprecision hnu(FRSIZE_FREQ,FRSIZE_R)
      doubleprecision fnu(FRSIZE_FREQ,FRSIZE_R)
c
#include "common_transfer.h"
c
      integer imu,nmu,ir
c
c     First compute the J_nu, H_nu and their freq-integrated 
c     counterparts J and H
c
      do ir=1,nr
          nmu = node_nmu(ir)
          jnu(inu,ir) = 0.d0
          hnu(inu,ir) = 0.d0
          do imu=1,nmu
              jnu(inu,ir) = jnu(inu,ir) + 
     %             node_intensity(-imu,ir) * node_dmu(-imu,ir) + 
     %             node_intensity(imu,ir) * node_dmu(imu,ir)
              hnu(inu,ir) = hnu(inu,ir) +
     %             node_intensity(-imu,ir) * 
     %                     node_mu(-imu,ir) * node_dmu(-imu,ir) + 
     %             node_intensity(imu,ir) * 
     %                     node_mu(imu,ir) * node_dmu(imu,ir) 
          enddo
          if(node_iray(0,ir).gt.0) then
              jnu(inu,ir) = jnu(inu,ir) + 
     %             node_intensity(0,ir) * node_dmu(0,ir)
          endif
          jnu(inu,ir) = 0.5 * jnu(inu,ir) 
          hnu(inu,ir) = 0.5 * hnu(inu,ir) 
      enddo
c
c     Compute the frequency-dependent Eddington factor (the traceless 
c     one, i.e. isotropic --> fnu=0), and the freq-integrated
c     counterpart f
c
      do ir=1,nr
          nmu = node_nmu(ir)
          fnu(inu,ir) = 0.d0
          do imu=1,nmu
              fnu(inu,ir) = fnu(inu,ir) + 0.5 * 
     %           ( node_intensity(-imu,ir) * 
     %             node_mu(-imu,ir)**2 * node_dmu(-imu,ir) + 
     %             node_intensity(imu,ir) * 
     %             node_mu(imu,ir)**2 * node_dmu(imu,ir) ) 
          enddo
          fnu(inu,ir) = fnu(inu,ir) / jnu(inu,ir) - (1.d0/3.d0)
      enddo
c
      end


c     --------------------------------------------------------------
c       COMPUTE J,H, THE MEAN OPACITIES AND THE VARIABLE EDD FACT           
c
c     This subroutine computes all that is necessary for the mean
c     variable eddington radiative transfer. One has to give it the
c     freq-dependent moments jnu,hnu,fnu, the opacity kappa(nu) and
c     the gas temperature T(z), and it will compute the frequency
c     integrated moments J(z) and H(z), the variable eddington factor 
c     f(z), and the mean opacities kappa(T), kappa(J) and kappa(H).
c
c     ARGUMENTS:
c        r        = R-grid
c        nu       = frequency-grid
c        kappa    = The full (temperature INdependent) opacity
c        temp     = The gas/dust temperature 
c        idokaph  = 0= do not compute kaph, 1= yes, do kaph as well
c        jnu      = Mean intensity J_nu 
c        Hnu      = Eddington flux H_nu = F_nu / 4 pi
c        fnu      = The eddington factor f_nu = K_nu / J_nu
c
c     OUTPUT:
c        jj       = The freq-integrated mean intensity
c        hh       = The freq-integrated eddington flux H = F / 4 pi
c        fj       = The freq-integrated eddington factor f = K / J
c        kapt     = The Planck opacity for temperature temp
c        kapj     = The J-mean opacity
c        kaph     = The H-mean opacity
c 
c     --------------------------------------------------------------
      subroutine compute_jh_meanopac_and_vet(nr,nf,
     %              r,nu,kappa,temp,jnu,hnu,fnu,idokaph,
     %              jj,hh,fj,kapt,kapj,kaph)
      implicit none
c
      integer nr,nf,idokaph
      doubleprecision r(FRSIZE_R),nu(FRSIZE_FREQ)
      doubleprecision kappa(FRSIZE_FREQ),temp(FRSIZE_R)
      doubleprecision jj(FRSIZE_R),hh(FRSIZE_R)
      doubleprecision jnu(FRSIZE_FREQ,FRSIZE_R)
      doubleprecision hnu(FRSIZE_FREQ,FRSIZE_R)
      doubleprecision fj(FRSIZE_R)
      doubleprecision fnu(FRSIZE_FREQ,FRSIZE_R)
      doubleprecision kapt(FRSIZE_R),kapj(FRSIZE_R),kaph(FRSIZE_R)
c
      integer hflag(FRSIZE_R)
      doubleprecision enum,denom,bnuold,bnucur,bplanck
      integer ir,inu
c
c     Compute the frequency-averaged moments
c
      do ir=1,nr
          jj(ir) = 0.d0
          hh(ir) = 0.d0
          do inu=2,nf
              jj(ir) = jj(ir) + 0.5 * ( jnu(inu,ir) + jnu(inu-1,ir) ) *
     %                abs( nu(inu) - nu(inu-1) ) 
              hh(ir) = hh(ir) + 0.5 * ( hnu(inu,ir) + hnu(inu-1,ir) ) *
     %                abs( nu(inu) - nu(inu-1) ) 
          enddo
          enum = 0.d0
          do inu=2,nf
              enum = enum + 0.5 * ( fnu(inu,ir) * jnu(inu,ir) + 
     %              fnu(inu-1,ir) * jnu(inu-1,ir) ) *
     %              abs( nu(inu) - nu(inu-1) ) 
          enddo
          fj(ir) = enum / jj(ir)
      enddo
c
c     Check if |H| <<< J. If so, we need to use the Rosseland 
c     mean opacity for the kaph, otherwise things might blow up
c     
      if(idokaph.ne.0) then
          do ir=1,nr
              hflag(ir) = 0
              if(abs(hh(ir)).lt.HFLUXMIN*jj(ir)) then
                  hflag(ir) = 1
                  write(*,*) 'H=0 at ',ir
                  stop 13       ! For the moment we simply stop!
              endif
          enddo
      endif
c     
c     Compute the Planck mean opacity
c
      do ir=1,nr
          bnucur = bplanck(temp(ir),nu(1))
          enum   = 0.d0
          denom  = 0.d0
          do inu=2,nf
              bnuold  = bnucur
              bnucur  = bplanck(temp(ir),nu(inu))
              enum    = enum + 0.5 * ( kappa(inu) * bnucur +
     %                    kappa(inu-1) * bnuold ) *
     %                    abs( nu(inu) - nu(inu-1) )
              denom   = denom + 0.5 * ( bnucur + bnuold ) *
     %                    abs( nu(inu) - nu(inu-1) )
          enddo
          kapt(ir) = enum / denom
      enddo
c
c     Compute the J-mean opacity
c
      do ir=1,nr
          enum = 0.d0
          do inu=2,nf
              enum    = enum + 0.5 * ( kappa(inu) * jnu(inu,ir) +
     %                    kappa(inu-1) * jnu(inu-1,ir) ) *
     %                    abs( nu(inu) - nu(inu-1) )
          enddo
          kapj(ir) = enum / jj(ir)
      enddo
c
c     Compute the H-mean opacity
c
      if(idokaph.ne.0) then
          do ir=1,nr
              enum = 0.d0
              do inu=2,nf
                  enum = enum + 0.5 * ( kappa(inu) * hnu(inu,ir) +
     %                 kappa(inu-1) * hnu(inu-1,ir) ) *
     %                 abs( nu(inu) - nu(inu-1) )
              enddo
              kaph(ir) = enum / hh(ir)
          enddo
      endif
c
      end


c     --------------------------------------------------------------
c                           NG ACCELLERATION 
c
c     This routine extrapolates the temperature, given the current
c     temperature and the temperature at three previous iteration
c     steps. This trick should only be applied every 4 iteration
c     steps. And the calling routine should keep track of the 
c     previous temperatures (t1,t2,t3) itself.
c     --------------------------------------------------------------
      subroutine do_ng_accell(nz,temp,t1,t2,t3)
      implicit none
c
      integer nz
      doubleprecision temp(FRSIZE_R)
      doubleprecision t1(FRSIZE_R)
      doubleprecision t2(FRSIZE_R)
      doubleprecision t3(FRSIZE_R)
c
      doubleprecision a1,a2,b1,b2,c1,c2,a,b,c,x0,x1,x2,x3,w
      doubleprecision dum1,dum2,dum3,dummy
      integer iz
c
      a1 = 0.d0
      a2 = 0.d0
      b1 = 0.d0
      b2 = 0.d0
      c1 = 0.d0
      c2 = 0.d0
      do iz=1,nz
          x0   = temp(iz)**4
          x1   = t1(iz)**4
          x2   = t2(iz)**4
          x3   = t3(iz)**4
          w    = 1.d0 / temp(iz)**4
          dum1 = x0 - 2.d0*x1 + x2
          dum2 = x0 - x1 - x2 + x3
          dum3 = x0 - x1
          a1   = a1 + dum1 * dum1 * w
          b1   = b1 + dum1 * dum2 * w
          c1   = c1 + dum1 * dum3 * w
          a2   = b1
          b2   = b2 + dum2 * dum2 * w
          c2   = c2 + dum2 * dum3 * w
      enddo
      a = ( c1*b2 - c2*b1 ) / ( a1*b2 - a2*b1 ) 
      b = ( c2*a1 - c1*a2 ) / ( a1*b2 - a2*b1 ) 
      c = 1.d0 - a - b
c
c     Write a message
c
      write(*,*) '       Ng : a = ',a,'  b = ',b
c
c     Now perform the accelleration step
c
      do iz=1,nz
          dummy = c * temp(iz)**4 + a * t1(iz)**4 + b * t2(iz)**4
          if(dummy.lt.0.d0) then
              write(*,*) 'Ng ERROR: negative T^4 at iz=',iz
              write(*,*) '   a,b,c = ',a,b,c
              write(*,*) 'Skipping this Ng step'
          else
              temp(iz) = sqrt(sqrt(dummy))
          endif
      enddo
c
      end



c     --------------------------------------------------------------
c                APPEND THE CONVERGENCE HISTORY FILE
c
c     This subroutine is only for convenience: it saves the stuff
c     as it converges towards a solution. This is handy for debug
c     and demonstration purposes.
c     --------------------------------------------------------------
      subroutine append_vet_convhist(nr,temp,jjme,hhme,jj,hh,
     %                    kapt,kapj,kaph,fj,savenr)
      implicit none
      integer nr,savenr
      doubleprecision temp(FRSIZE_R),jj(FRSIZE_R),hh(FRSIZE_R)
      doubleprecision kapt(FRSIZE_R),kapj(FRSIZE_R),kaph(FRSIZE_R)
      doubleprecision fj(FRSIZE_R),jjme(FRSIZE_R),hhme(FRSIZE_R)
c
      integer ir
c
      open(unit=1,file='convhist.dat',status='old',access='append')
      do ir=1,nr 
          write(1,10) temp(ir),jjme(ir),hhme(ir),jj(ir),hh(ir),
     %              kapt(ir),kapj(ir),kaph(ir),fj(ir)
 10       format(9(E13.6,1X))
      enddo
      write(1,*)
      close(1)
      open(unit=1,file='convhist.info',status='unknown')
      write(1,*) savenr,1
      close(1)
c
      end



c     --------------------------------------------------------------
c                         WRITE SPECTRUM
c     --------------------------------------------------------------
      subroutine write_spectrum
      implicit none
c
#include "common_grid.h"
#include "common_transfer.h"
c
      integer inu
c
      open(unit=1,file='spectrum.dat',status='unknown')
      write(1,*) freq_nr
      write(1,*) 
      do inu=1,freq_nr
          write(1,10) freq_nu(inu),spim_spectrum(inu)+1d-90
 10       format(2(E13.6,1X))
      enddo
c
      end
