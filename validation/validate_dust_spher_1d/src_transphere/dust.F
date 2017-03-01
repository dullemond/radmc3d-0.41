c     =============================================================
c                2-D SHORT CHARACTERISTIC TRANSFER ALGORITHM
c                 FOR POLAR COORDINATES AND AXIAL SYMMETRY
c
c                    ( MODULE: DUST OPACITY PACKAGE ) 
c
c                           Leiden, March 1999
c                    C.P. Dullemond / Alex de Koter
c
c     =============================================================
c
c     The general idea of this module:
c     
c     Rather than computing the emission and opacities from the
c     medium_rho and medium_temp, here we have a collection of
c     different dust species, each consisting of grain size 
c     distribution which is modeled by sampling the distribution.
c     The density of these species must be specified in the array
c     dust_rho(), specifying the density of each of the species
c     of dust particles. Note that the density is specified for
c     each species separately, but NOT for each grainsize separately.
c     The grainsizes are integrated to yield an average grain 
c     weight, so that the average grain density can be computed,
c     which is supposed to be done by the routines in the module
c     called "dustopacities.F", which at the time of writing this
c     is still a separate program that produces an input file
c     called "dustopac.inp", which the dust.F module can read in.
c     The dust temperatures, on the other hand, _are_ specified for
c     each species AND grain size. 
c
c     This module can compute the scattering and absorbtion for 
c     setups involving varying dust species abundances over space. 
c     Also it can solve the temperatures of the dust particles.
c
c     NOTE: The solving of the temperatures is not yet compatible
c           with ALI.
c
c       Leiden, 11 april 1999
c     =============================================================

#include "main.h"


c     =============================================================
c
c                 THE DUST DATA AND DUST OPACITY PART
c
c     This part of the dust module one can use also in other 
c     programs. It requires externally:
c
c       - A frequency array freq_nu(1..freq_nr) in which the 
c         frequencies of the dust opacity files are given. This
c         must be located in the common_grid.h file.
c       - A subroutine that reads these frequencies (if you dont
c         have such a routine, uncomment the line with the 
c         #define DUST_PRIVATE_EXTERNALS, and the routine
c         read_frequencies() will be provided to you, see below). 
c       - Define the following constants:
c           DUST_SPECIES_MAX   Maximum nr of dust species
c           DUST_SIZE_MAX      Maximum nr of grain sizes per species
c           DUST_TRANGE_MAX    Maximum nr of temperature ranges 
c                              for temperature-dependent opacities
c           INCLUDE_DUST       Just define this to activate module
c
c     Then, to use the module, do these things:
c
c       1) Read the frequency grid [call read_frequencies()]
c       2) Read the dust opacity data [call read_dustdata()]
c       3) To get the dust kappa for a certain species ispec,
c          a grain size isize, a temperature temp, call the
c          function find_dust_kappa().
c
c     =============================================================



c     --------------------------------------------------------------
c                       READ THE DUST DATA
c
c     This is the main routine for the input of all relevant dust
c     data, and the preprocessing of this data, if necessary.
c
c     REMARK: dust_done_read is a flag signalling to the program
c             that we have (fortunately) not forgotten to read in
c             the dust data. This is only for internal consistency
c             checking.
c     --------------------------------------------------------------
      subroutine read_dustdata()
      implicit none
c
#include "common_grid.h"
#include "common_dust.h"
c
      integer ispec,inu,isize,iformat,idum,idustfile
      doubleprecision dummy,temp0,temp1,dtemp0,dtemp1
      character*80 comstring,filename,base,ext
      logical fex_fwgt
c
c     Now open the master file for the dust
c
      open(unit=3,file='dustopac.inp',status='old',err=701)
      read(3,*) iformat
      read(3,*) dust_nr_species
      read(3,*) comstring      
c
c     Now a loop over the dust species
c
      do ispec=1,dust_nr_species
          read(3,*) idum
          read(3,*) idustfile
          base='dustopac_'
          ext ='.inp'
          call make_indexed_filename(base,idustfile,ext,filename)
c
          if(idum.eq.-1) then
c
c             Simple temperature-independent opacities
c
c              write(*,*) filename
              call read_dustopac_file(ispec,filename)
          elseif(idum.eq.-2) then
c
c             Include the minimal and maximal temperature for
c             this species.
c

              read(3,*) temp0
              read(3,*) temp1
              dust_tmin(ispec)  = temp0
              dust_tmax(ispec)  = temp1
              dust_dtmin(ispec) = 0.d0
              dust_dtmax(ispec) = 0.d0
c              if((temp0.ge.0.d0).or.(temp1.gt.0.d0)) then
c                  iradproc_dust_tminmax = 1
c              endif
              call read_dustopac_file(ispec,filename)
          elseif(idum.eq.-3) then
c
c             Not only tmin/tmax, but also smooth switch-off
c

              read(3,*) temp0
              read(3,*) dtemp0
              read(3,*) temp1
              read(3,*) dtemp1
              dust_tmin(ispec)  = temp0
              dust_tmax(ispec)  = temp1
              dust_dtmin(ispec) = abs(dtemp0)
              dust_dtmax(ispec) = abs(dtemp1)
c              if((temp0.ge.0.d0).or.(temp1.gt.0.d0)) then
c                  iradproc_dust_tminmax = 1
c              endif
              call read_dustopac_file(ispec,filename)
          elseif(idum.eq.-4) then
c
c             This is a near-future mode in which the dust opacities
c             can be really temperature-dependent
c
              write(*,*) 'Sorry, dust input type ',idum,' not yet',
     %             'finished.'
              stop 55
          else
              write(*,*) 'While reading dustopac.inp: other ',
     %         'input modes than -1 not yet implemented'
              stop 13
          endif
          read(3,*) comstring                
      enddo
c
      close(3)
c
c     Check if the nr of dust species found here agrees with
c     the number that the transfer routines expect. If the
c     dust_setup_nrspecies=-1 then this is a signal of the
c     transfer routine that the number of dust species is in
c     fact determined by the dustopac.inp file itself.
c
      if((dust_setup_nrspecies.gt.0).and.
     %    (dust_setup_nrspecies.ne.dust_nr_species)) then
          write(*,*) 'INCONSISTENCY, read_dustdata():'
          write(*,*) '  The dust file dustopac.inp specifies ',
     %                 dust_nr_species,' dust species.'
          write(*,*) '  The radiative transfer setup expects ',
     %                 dust_setup_nrspecies,' dust species.'
          write(*,*) '  (=dust_setup_nrspecies).'
          write(*,*) '  So: simulation and dust files do not '
          write(*,*) '  agree on the number of dust species.'
          write(*,*) '  modify either of the two. Now I quit.'
          stop 99
      endif
      if(dust_setup_nrspecies.lt.0) then
          dust_setup_nrspecies = dust_nr_species
      endif
c
      goto 710
  701 continue
      write(*,*) 'Could not open file dustopac.inp'
      stop 13
  710 continue
c
c     Now read the dust opacity frequency weights, if they have 
c     not yet been specified (which can be seed from dust_frwgt_read)
c
      inquire(file='dustopac_fwgt.inp',EXIST=fex_fwgt)
      if(fex_fwgt.and.(dust_frwgt_read.eq.0)) then
          open(unit=1,file='dustopac_fwgt.inp',status='old',err=701)
          read(1,*) idum
          if(idum.ne.freq_nr) then
              write(*,*) 'ERROR: file dustopac_fwgt.inp ',
     %            'is incompatible with frequency.dat'
              stop 13
          endif
          do inu=1,idum
              read(1,*) dummy
              dust_freq_wgt(inu) = dummy
          enddo
          close(1)
          dust_frwgt_read = 1
      endif
c
c     If the frequency weights have still not been specified (for 
c     whatever reason), then we generate them here
c
#     ifndef EXCLUDE_DUST_ITERATION
      if(dust_frwgt_read.eq.0) then
          call make_dust_freq_weights()
          dust_frwgt_read = 1
      endif
#     endif
c
      dust_done_read = 1
      return
      end



c     --------------------------------------------------------------
c                     READ THE DUST OPACITY FILES 
c
c     This routines read the ready-to-use data from the file
c     dustopac.inp. This file contains the dust opacities and
c     weights for the frequencies which are used in this simulation.
c     There are two possible ways these tables can be given. One
c     way (the old way) is to give a table that is INdependent off
c     the temperature. The other is to specify the opacities for
c     a series of temperature ranges (in each range the opacity 
c     remains independent off temperature, but may differ from range
c     to range). The latter option is automatically selected by
c     starting the file with the number -1. See below how the format
c     is then.
c
c     NOTE:   The dust opacities are assumed to have been computed 
c             elsewhere and written readily into this file. 
c
c     REMARK: This routine will be called by read_dustdata(), which
c             will keep unit=3 open.
c     --------------------------------------------------------------
      subroutine read_dustopac_file(ispec,filename)
      implicit none
      character*80 filename
      integer ispec
c
#include "common_dust.h"
#include "common_grid.h"
c
      integer ifr,isize,nsize,itemp,ntemp,idum
      doubleprecision dummy,a,b,bold
      logical readtrange
c
      ntemp = 1
      readtrange = .false.
c
#     ifdef DUST_RADICAL
      if(freq_grid_type.ne.-1) then
          write(*,*) 'ERROR: In radical.inp the frequency grid type'
          write(*,*) '       is specified as ',freq_grid_type
          write(*,*) '       But when dust opacity files are used,'
          write(*,*) '       the frequency grid MUST be read from '
          write(*,*) '       the file "frequency.inp", since the '
          write(*,*) '       opacity files dustopac_*.inp contain '
          write(*,*) '       the opacities at those frequencies. '
          write(*,*) '          In the future perhaps we will allow'
          write(*,*) '       other frequency grids, and use some kind'
          write(*,*) '       of interpolation to map the opactities '
          write(*,*) '       from the frequency.inp grid onto another '
          write(*,*) '       frequency grid. But for now: specify '
          write(*,*) '       frequency grid type -1 in radical.inp'
          stop 13
      endif
#     endif
c
      open(unit=1,file=filename,status='old',err=701)
      read(1,*) ifr,nsize
      if(ifr.eq.-1) then
c
c         ifr=-1 is a signal to say that the dust opacity file consists
c         of opacities at a series of temperature points. 
c
          readtrange = .true.
          ifr=nsize
          read(1,*) nsize,idum,ntemp
c
c         Check also if the DUST_OPAC_TEMPDEP is set.
c
#         ifndef DUST_OPAC_TEMPDEP
          write(*,*) 'PROBLEM: The code is compiled without the'
          write(*,*) '         DUST_OPAC_TEMPDEP option. But the '
          write(*,*) '         dust opacity file',filename,
     %          ' in fact has'
          write(*,*) '         multiple temperature ranges.'
          write(*,*) '         Recompile with DUST_OPAC_TEMPDEP in'
          write(*,*) '         the configure.h file!'
          stop 13
#         endif
      else
c
c         Else, there is only one (fixed) opacity table for the dust,
c         independent on temperature (with the exception that there is
c         a Tmin and Tmax within which this species can survive). 
c
          ntemp = 1
      endif
      dust_nr_size(ispec) = nsize
      dust_nr_temp(ispec) = ntemp
      if(ifr.ne.freq_nr) then
          write(*,*) 'Number of frequencies in ',filename
          write(*,*) 'not equal to the number in frequency.inp'
          write(*,*) ifr,freq_nr
          stop 13
      endif
      if(nsize.gt.DUST_SIZE_MAX) then
          write(*,*) 'Number of dust sizes for species ',ispec
          write(*,*) 'exceeds maximum as specified in Makefile'
          stop 13
      endif
      if(ispec.gt.DUST_SPECIES_MAX) then
          write(*,*) 'Dust species identifier is larger than '
          write(*,*) 'maximum as specified in Makefile'
          stop 13
      endif
      if(ntemp.gt.DUST_TRANGE_MAX) then
          write(*,*) 'Sorry, the nr of dust opacity temperature'
          write(*,*) 'points is larger than DUST_TRANGE_MAX'
          stop 13
      endif
c
c     Now read the ready-to-use kappa*wgt values. These values
c     include already the relative abundances of the various
c     dust sizes.
c
c     NOTE: the temperature ranges should (for the moment) be 
c           exactly matching: the upper temp of the previous range
c           should equal the lower temp of the next. In the future
c           this will no longer be necessary. Then there may be a
c           gap in between: linear interpolation will then be used
c           to bridge the gap. By taking the lower and upper temp
c           of the same temp range to be equal, one naturally 
c           reproduces the other interesting possibility: gradually
c           varying opacities as a function of temperature, with
c           linear interpolation. 
c
c     Make sure that, if the dust opacity table is simply a constant
c     opacity (as a function of temperature), that the first temperature
c     range is from 0 to infinity (=1d33). 
c
      dust_temprange_low(1,ispec)      = 0.d0
      dust_temprange_high(1,ispec)     = 1.d33
c
c     Now read tables in
c
      bold = 0.d0
      do itemp=1,ntemp
          if(readtrange) then
              read(1,*) a,b
              dust_temprange_low(itemp,ispec)  = a
              dust_temprange_high(itemp,ispec) = b
              dust_opacity_tempdep = 1
              if(a.lt.bold) then
                  write(*,*) 'ERROR while reading ',filename
                  write(*,*) '    for dust species ',ispec
                  write(*,*) '    The minimum temperature for '
                  write(*,*) '    opacity temp range ',itemp,' is '
                  write(*,*) '    smaller than the maximum for '
                  write(*,*) '    opacity temp range ',itemp-1
                  stop 13
              endif
              if(abs((a-bold)/(a+bold)).gt.1.d-3) then
                  write(*,*) 'ERROR in dust.F'
                  write(*,*) '    Sorry: I can not yet handle'
                  write(*,*) '    temperature-dependent dust opacity'
                  write(*,*) '    tables with gaps in temperature.'
                  write(*,*) '    So Tmax of one temp-range must equal'
                  write(*,*) '    Tmin of the next temp-range, in the'
                  write(*,*) '    file ',filename
                  stop 13
              endif
              if(a.lt.0.d0) then
                  write(*,*) 'PROBLEM: while reading ',filename
                  write(*,*) '    found negative minimum temperature'
                  stop 13
              endif
              if(b.lt.0.d0) then
                  write(*,*) 'PROBLEM: while reading ',filename
                  write(*,*) '    found negative maximum temperature'
                  stop 13
              endif
              if(b.lt.a) then
                  write(*,*) 'PROBLEM: while reading ',filename
                  write(*,*) '    minimum temp is larger than max'
                  stop 13
              endif
              bold = b
          else
              if(itemp.gt.1) then
                  write(*,*) 'BUG: Inconsistency in dust.F/',
     %                            'read_dustopac_file()'
                  stop 13
              endif
          endif
          do ifr=1,freq_nr
              do isize=1,nsize
                  read(1,*) dummy
                  dust_kappawgt_abs(ifr,itemp,isize,ispec) = dummy
              enddo
          enddo
          do ifr=1,freq_nr
              do isize=1,nsize
                  read(1,*) dummy
                  dust_kappawgt_scat(ifr,itemp,isize,ispec) = dummy
              enddo
          enddo
      enddo
c
c     Now, if there was a minimum or maximum temperature of the
c     dust selected, put this into the t-range database
c
      if(dust_tmax(ispec).gt.1.d-2) then
          dust_temprange_high(ntemp,ispec) = min(dust_tmax(ispec),
     %                          dust_temprange_high(ntemp,ispec))
          dust_opacity_tempdep = 1
      endif
      if(dust_tmin(ispec).gt.0.d0) then
c       >>> bugfix 06-08-00: .ge. --> .gt. <<<
          dust_temprange_low(1,ispec) = max(dust_tmin(ispec),
     %                          dust_temprange_low(1,ispec))
          dust_opacity_tempdep = 1
      endif
c
c     Done...
c     
      close(1)
c
      goto 710
  701 continue
      write(*,*) 'Could not open file ',filename
      stop 13
  710 continue
      return
      end


      
c     --------------------------------------------------------------
c                         FIND DUST OPACITY
c     
c     This function returns the dust opacity of a given species,
c     size, and at a given frequency-index. It takes consistently
c     into account the temperature range in which 
c
c     --------------------------------------------------------------
      function find_dust_kappa(inu,isize,ispec,temp,iabs,iscat)
      implicit none
      doubleprecision temp,find_dust_kappa
      integer inu,isize,ispec,iabs,iscat
c
#include "common_grid.h"
#include "common_dust.h"
c
      doubleprecision tmin,tmax,dtmin,dtmax,plindex
      doubleprecision condense,fact,omfact
      integer ntemp,itlo,ithi
c
      condense = 1.d0
      ntemp = dust_nr_temp(ispec)
      tmin  = dust_temprange_low(1,ispec)
      tmax  = dust_temprange_high(ntemp,ispec)
      dtmin = dust_dtmin(ispec)
      dtmax = dust_dtmax(ispec)
      itlo  = 1
      ithi  = 1
      fact  = 1.d0
c
c     Now take care of temperature regimes
c
      if(temp.lt.tmin) then
          if(dtmin.eq.0.d0) then
c             
c             Hard switch-off at low temperature
c         
              condense = 0.d0
          else
c         
c             Smooth switch-off at low temperature
c         
              plindex  = log(1.d0-(dtmin/tmin))
              condense = (tmin/temp)**plindex
              if(condense.gt.1.d0) then
                  stop 9932
              endif
          endif
      elseif(temp.ge.tmax) then
          if((dtmax.eq.0.d0).and.(tmax.gt.0.d0)) then
c         
c             Hard switch-off at high temperature
c         
              condense = 0.d0
          else
c         
c             Smooth switch-off at high temperature
c         
              plindex  = log(1.d0+(dtmax/tmax))
              condense = (tmin/temp)**plindex
              if(condense.gt.1.d0) then
                  stop 9932
              endif              
          endif
      else
          if(ntemp.gt.1) then
              call hunt(dust_temprange_low(1,ispec),ntemp,temp,itlo)
              if(temp.gt.dust_temprange_high(itlo,ispec)) then
c
c                 Consistency check
c
                  if(itlo.ge.ntemp) then
                      write(*,*) 'BUG: Inconsistency in dust opacity'
                      write(*,*) '     temperature range....'
                      stop 13
                  endif
c
c                 We are in the linear interpolation regime
c
                  ithi = itlo + 1
                  fact = ( temp - dust_temprange_high(itlo,ispec) ) 
     %                   / ( dust_temprange_low(ithi,ispec) -
     %                       dust_temprange_high(itlo,ispec) )
              else
c
c                 We are within a range of constant opacity
c
                  ithi = itlo
                  fact = 1.d0
              endif
          endif
      endif
      omfact = 1.d0 - fact
c
c     Now get the opacities
c
      find_dust_kappa = 0.d0
      if(iabs.ne.0) then
          find_dust_kappa = find_dust_kappa +
     %           omfact * dust_kappawgt_abs(inu,itlo,isize,ispec) +
     %           fact * dust_kappawgt_abs(inu,ithi,isize,ispec)
      endif
      if(iscat.ne.0) then
          find_dust_kappa = find_dust_kappa +
     %           omfact * dust_kappawgt_scat(inu,itlo,isize,ispec) +
     %           fact * dust_kappawgt_scat(inu,ithi,isize,ispec)
      endif
c
c     If either iscat or iabs .lt.0 then do not incorporate the
c     evaporation of the grain. This can be handy to ensure that
c     at every temperature (even for the unphysically large ones)
c     a dust temperature can be computed. Simply for consistency.
c
      if((iabs.ge.0).and.(iscat.ge.0)) then
          find_dust_kappa = find_dust_kappa * condense
      endif
c
      return
      end




c     --------------------------------------------------------------
c                COMPUTE THE DUST FREQUENCY WEIGHTS
c
c     This routine computes the integration weights for the integrals
c     over the frequency. In principle this is a very general thing
c     and should not be here in this specialized module. But for now
c     this is the way we do it. Be sure to have generated the frequency
c     grid before calling this routine.
c
c     The frequency weight is such that
c                        __
c         /              \
c         | f(nu) dnu  =  >  f(nu_i) wgt_i
c         /              /_
c
c     We use the simple trapezium rule for the integral. 
c
c     --------------------------------------------------------------
      subroutine make_dust_freq_weights()
      implicit none
c
#include "common_grid.h"
#include "common_dust.h"
c
      integer ifreq
c
      if(freq_nr.lt.1) then
          write(*,*) 'Zero or negative Number of freqs impossible'
          stop 99
      elseif(freq_nr.eq.1) then
cc          write(*,*) 'PROBLEM in dust module:'
cc          write(*,*) '  Attempt to construct frequency integration'
cc          write(*,*) '  weights for case of 1-frequency computation.'
cc          write(*,*) '  This is inconsistent. Are you sure you have '
cc          write(*,*) '  not accidently forgotten to switch off the '
cc          write(*,*) '  dust-temperature solver? This can be done '
cc          write(*,*) '  by setting iradproc_dust to -1 instead of 1.'
cc          stop 13
          return
      elseif((freq_nr.ge.2).and.(freq_nr.le.5)) then
          write(*,*) 'WARNING for dust temperature solver: '
          write(*,*) '  Too few frequency bins for consistently solving'
          write(*,*) '  dust temperatures! Continuuing while hoping '
          write(*,*) '  for the best....'
          dust_warn_few_freqs = 1
      else
          dust_freq_wgt(1) = 0.5d0 * abs( freq_nu(2) - freq_nu(1) )
          dust_freq_wgt(freq_nr) = 0.5d0 * abs( freq_nu(freq_nr) 
     %                              - freq_nu(freq_nr-1) )
          do ifreq=2,freq_nr-1
              dust_freq_wgt(ifreq) = 0.5d0 * 
     %             abs( freq_nu(ifreq+1) - freq_nu(ifreq-1) )
          enddo
          dust_frwgt_read = 1
      endif
c     
      end





c     --------------------------------------------------------------
c         READ THE FREQUENCY ARRAY (FOR USE OUTSIDE OF RADICAL)
c
c     If this module is used outside of RADICAL, then this module
c     itself will provide a routine for reading the frequency grid.
c     --------------------------------------------------------------
      subroutine read_frequencies()
      implicit none
c
#include "common_grid.h"
c
      integer inu
      doubleprecision frnu,dfrnu
c
      open(unit=1,file='frequency.inp',status='old',err=701)
      read(1,*) freq_nr
      if(freq_nr.gt.FRSIZE_FREQ) then
          write(*,*) 'ERROR: frequency.inp has more gridpoints'
          write(*,*) '       than FRSIZE_FREQ. '
          write(*,*) '       Recompile with larger FRSIZE_FREQ'
          stop 13
      endif          
      do inu=1,freq_nr
          read(1,*) freq_nu(inu)
      enddo
      close(1)
      goto 710
  701 continue
      write(*,*) 'Could not open file frequency.inp'
      stop 13
  710 continue
      return
      end

