c     ==============================================================
c
c                         >>>> TRANSPHERE <<<<
c
c     A DUST CONTINUUM RADIATIVE TRANSFER CODE FOR SPHERICAL SYMMETRY
c
c              written by C.P. Dullemond (C) March 2002
c
c     This code solves the radiative transfer problems in spherical
c     symmetry of a dust envelope surrounding a star with a certain
c     temperature and luminosity.  The dust grains are assumed to be in
c     radiative equilibrium, and the total flux is conserved. The code
c     solves the dust temperature at every radius and returns the
c     spectral energy distribution of the system. The code can solve
c     problems of any optical depth within about 10 iterations
c     only. This powerful capability is made possible through the use of
c     variable eddington factors and variable mean opacities.
c     
c     In the present (downscaled) form, the code can handle only a
c     single dust species and size, and dust scattering is not
c     included. Also, the code only reads in the stellar temperature and
c     size, and assumes the star to be a blackbody. The full version of
c     the code does not have these limitations and can input a realistic
c     spectrum, can handle multiple species and can include scattering
c     (in the isotropic limit).
c     
c     If you like the present downscaled version of the code, and wish
c     to use it for scientific use, you can obtain a full version of the
c     code by contacting the author via email:
c     dullemon@mpa-garching.mpg.de. Conditions for use can then be
c     agreed upon through email exchange.
c     
c                  ----------------------------------
c                  !!!!!!!  Copyright notice: !!!!!!!
c                  ----------------------------------
c     
c     The present code is public domain, but remains under copyright.
c     The code or any part of the code) may NOT be used for (any part
c     of) a publication in a refereed scientific journal without written
c     permission from the author. You can reach the author via email:
c     dullemon@mpa-garching.mpg.de.
c     
c     This package may be distributed to others, under the condition
c     that it it kept intact and complete, including this README file.
c     
c     Note: 
c     
c       The package contains a few subroutines from Numerical Recipes. 
c       This is copyrighted software, so note the copyright restrictions 
c       of that software. For more information on Numerical Recipes, see 
c       their webpage: 
c           www.ulib.org/webRoot/Books/Numerical_Recipes/bookf.html
c     
c     Note:
c     
c       The program TRANSPHERE is written in fortran 77. An example setup 
c       is created by a program written in the IDL programming language 
c       (see www.rsinc.com). Also the example read and plot routines for 
c       the output are written in IDL. IDL is licenced software, and can 
c       be purchased from Research Systems Inc. 
c     
c     Website:
c     
c        More information about various astrophysical software packages 
c        written by C.P. Dullemond can be found on the following homepage:
c     
c            http://www.mpa-garching.mpg.de/~dullemon
c     
c            (NOTE: This page remains here until 2004; After leaving MPA,
c             this page will be moved. Unfortunately, MPA will not support
c             a forward link to my new address.)
c     
c     ==============================================================

#include "main.h"

#define PLANCKTAB_NT 100
#define PLANCKTAB_T0 (1.D-1)
#define PLANCKTAB_T1 (1.D+4)

c     --------------------------------------------------------------
c                             MAIN ROUTINE
c     --------------------------------------------------------------
      program main
      implicit none
c
      call read_all_data()
      call do_all_radtrans()
      call write_temperatures()
      call write_aux_data()
      call write_spectrum()
c
      end



c     --------------------------------------------------------------
c                         READ ALL THE DATA IN
c     --------------------------------------------------------------
      subroutine read_all_data()
      implicit none
c
#include "common_grid.h"
#include "common_dust.h"
#include "common_radtrans.h"
c
      doubleprecision rstar,a,b,c,d,e
      integer nr,ir,nrp,nf,inu
      integer iformat
      logical f_ex,fexf,fexo
c
c     Read the fullvrt input file
c
      global_init=0
      global_dump_intens=0
      open(unit=1,file='transphere.inp',status='old')
      read(1,*) iformat
      read(1,*) vet_itermax
      read(1,*) vet_convcrit
      read(1,*) global_ncst
      read(1,*) global_ncex
      read(1,*) global_ncnr
      if(iformat.ge.2) then
          read(1,*) global_itypemw
      else
          global_itypemw = 1
      endif
      read(1,*) global_idump
      read(1,*,end=100,err=100) global_init
      read(1,*,end=100,err=100) global_dump_intens
      close(1)
 100  continue
c
c     Read the star information
c
      open(unit=1,file='starinfo.inp',status='old')
      read(1,*) iformat
      read(1,*) star_r
      read(1,*) star_m
      read(1,*) star_t
      close(1)
c
c     Read the primary temperature and the density
c
      open(unit=1,file='envstruct.inp',status='old')
      read(1,*) nr
      if(nr.gt.FRSIZE_R) then
          write(*,*) 'ERROR:'
          write(*,*) 'file envstruct.inp: R grid too large'
          stop 13
      endif
      grid_nr = nr
      do ir=1,nr
          read(1,*) a,b,c
          grid_r(ir)=a
          srt_rho(ir)=b
          srt_temp(ir)=c
      enddo
      close(1)
c
c     Read the frequency array and the opacities
c
      inquire(file='frequency_hz.inp',exist=fexf)
      inquire(file='dustopac_1.inp',exist=fexo)
      if((.not.fexf).or.(.not.fexo)) then
          write(*,*) 'Need both frequency_hz.inp and dustopac_1.inp'
          stop 13
      endif
      open(unit=1,file='frequency_hz.inp',status='old')
      read(1,*) nf
      if(nf.gt.FRSIZE_FREQ) then
          write(*,*) 'ERROR in main:'
          write(*,*) 'Frequency grid too large'
          stop 13
      endif
      do inu=1,nf
          read(1,*) a
          freq_nu(inu) = a
      enddo
      dust_setup_nrspecies=1
      dust_setup_nrsizes(1)=1
      freq_nr = nf
      call read_dustdata()
c
      end



c     --------------------------------------------------------------
c                         PERFORM ALL RT
c     --------------------------------------------------------------
      subroutine do_all_radtrans()
      implicit none
c
#include "common_grid.h"
#include "common_dust.h"
#include "common_radtrans.h"
#include "common_transfer.h"
c
      doubleprecision dumtemp,find_dust_kappa,t0,t1
      integer ir,inu,nt
      parameter(nt=PLANCKTAB_NT)
      parameter(t0=PLANCKTAB_T0)
      parameter(t1=PLANCKTAB_T1)
c
c     First make a table of opacities, both the freq-dependent
c     ones, and the Planck-mean ones.
c
      if(dust_opacity_tempdep.ne.0) then
          write(*,*) 'Sorry, but fullrt cannot use temperature'
          write(*,*) 'dependent opacities for the moment...'
          stop 13
      endif
      dumtemp = 10.d0
      do inu=1,freq_nr
          optab_kappa(inu) = find_dust_kappa(inu,1,1,t0,1,0)
      enddo
      call compute_planckopac(freq_nr,nt,freq_nu,t0,t1,optab_kappa)
c
c     Do RT and check energy conservation
c
      call do_radtrans()
c      call check_encons()
c
      end




c     --------------------------------------------------------------
c                DO THE VERTICAL TRANSFER FOR ONE SHELL
c     --------------------------------------------------------------
      subroutine do_radtrans()
      implicit none
c
#include "common_grid.h"
#include "common_transfer.h"
#include "common_radtrans.h"
c
      doubleprecision bplanck,starint(FRSIZE_FREQ)
      integer ierror,ing,iali,imethod,intdump,ir,inu,intorder
c
c     Compute the source term for the vertical transfer
c
      do ir=1,grid_nr
          srt_qfpi(ir) = 0.d0
      enddo
c
c     Set up the stellar spectrum as surface intensity
c
      do inu=1,freq_nr
          starint(inu) = bplanck(star_t,freq_nu(inu))
      enddo
c
c     Now do the full transfer 
c
      ing = 0                   ! Yes Ng acceleration
      intdump=0
      intorder=1
      call do_vet_transfer(grid_nr,freq_nr,global_ncst,global_ncex,
     %              global_ncnr,grid_r,freq_nu,srt_rho,srt_qfpi,
     %              star_r,starint,vet_itermax,vet_convcrit,ing,
     %              intorder,global_itypemw,global_idump,intdump,
     %              global_init,ierror,optab_kappa,srt_temp,
     %              vet_jjme,vet_hhme,vet_jj,vet_hh,
     %              vet_jnu,vet_hnu,vet_fj,vet_fnu,vet_kapt,
     %              vet_kapj,vet_kaph)
c
c     Check error message
c
      if(ierror.ne.0) then
          write(*,*) '===== PROBLEM IN FULLRT AT IR=',ir,' ====='
      endif
c
      end




c     --------------------------------------------------------------
c                THE BLACKBODY PLANCK FUNCTION B_nu(T)
c
c     This function computes the Blackbody function 
c
c                    2 h nu^3 / c^2
c        B_nu(T)  = -------------------    [ erg / cm^2 s ster Hz ]
c                   exp(-h nu / kT) - 1
c
c     ARGUMENTS:
c        nu    [Hz]            = Frequency
c        temp  [K]             = Electron temperature
c     --------------------------------------------------------------
      function bplanck(temp,nu)
      implicit none
      doubleprecision temp
      doubleprecision nu
      doubleprecision bplanck
c
      if(temp.eq.0.d0) then 
          bplanck = 0.d0
          return
      endif
c
      bplanck = 1.47455d-47 * nu * nu * nu / 
     %      (exp(4.7989d-11 * nu / temp)-1.d0) + 1.d-290
c
      return
      end





c     --------------------------------------------------------------
c                    WRITE THE MAIN OUTPUT FILE
c     --------------------------------------------------------------
      subroutine write_temperatures()
      implicit none
c
#include "common_grid.h"
#include "common_radtrans.h"
c
      integer ir
c
      open(unit=1,file='envstruct.dat',status='unknown')
      write(1,*) grid_nr
      write(1,*)
      do ir=1,grid_nr
          write(1,10) grid_r(ir),srt_rho(ir),srt_temp(ir)
 10       format(3(E13.6,1X))
      enddo
      close(1)
c
      end


c     --------------------------------------------------------------
c                  WRITE SOME OF THE AUXILIARY DATA
c     --------------------------------------------------------------
      subroutine write_aux_data()
      implicit none
c
#include "common_grid.h"
#include "common_radtrans.h"
c
      integer ir
c
      open(unit=1,file='transphere.dat',status='unknown')
      write(1,*) grid_nr
      write(1,*)
      do ir=1,grid_nr
          write(1,10) vet_hh(ir),vet_jj(ir),
     %                vet_fj(ir),srt_qfpi(ir)
 10       format(4(E13.6,1X))
      enddo
      close(1)
c
      end


