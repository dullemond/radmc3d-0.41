c     ==============================================================
c              SUBPROGRAM FOR COMPUTING THE DUST OPACITIES
c                        Author: Alex de Koter
c               Modified with permission: C.P. Dullemond   
c                             19-04-99
c
c     ==============================================================
c
c     Example input file (strip off first 6 character columns):
c 
c     *
c     * SELECT  ABUNDANCE  AMIN  AMAX  NA  MA  RHO  ALPHA
c     *
c     ** family: pyroxene
c     ** latice structure: C
c     ** chem. comp.: mg0.9fe0.1
c     ** ID: 001.lnk
c     *
c     /afs/mpa/home/dullemon/science/database/dust/dbase2/pyroxene/C/mg0.9fe0.1_02.80_0300_001.lnk
c     "MIE"
c     1     0.00    -2.00     0.0    10    -3.50    02.80    -1.0
c
c     ==============================================================

#include "configure.h"

#define PICONST 3.1415926535898D0

c     --------------------------------------------------------------
c                           MAIN ROUTINE
c
c     This program computes the opacities of a collection of dust
c     species, each with a certain size distribution. The opacities
c     are computed for each frequency found in the file frequency.inp.
c     The dust species, size distributions and other relevant data
c     can be selected in the file dust.inp. The program computes the
c     opacity for each grain type separately and independently. 
c     Each grain type has various grain sizes. The opacity for each
c     of the discrete size-samples is computed PER unit TOTAL density
c     of the grain species (i.e. the density averaged over all grain
c     sizes!). So if you have -say- 7 g/cm^3 of olivine dust grains 
c     having a certain size distribution, then the total opacity 
c     (in [cm^-1]) would be: alpha = 7 * sum_j kappa_j, where kappa_j
c     are the values produced by this program, which are written
c     in files named "dustkappa_1.inp" etc. The "1" etc stand for the
c     dust species. Each file contains the opacities of all the sizes
c     of that species.
c
c     This program looks for the files:
c       - frequency.inp
c       - dust.inp
c       - the data files specified in dust.inp
c
c     The program produces the files:
c       - dustopac.inp
c       - dustkappa_1.inp .... dustkappa_*.inp
c     --------------------------------------------------------------
      program dustopacities
      implicit real*8 (a-h,o-z)
      parameter (naa=DUST_SIZE_MAX)
      parameter (ngg=DUST_SPECIES_MAX)
      parameter (nff=FRSIZE_FREQ)
      parameter (pi=PICONST)
      logical fexists
      character*99 ch,filename
c
      common/freq/waves(nff),freqs(nff)
      common/ifreq/nf
c
      common /opac/ opa_wgt(naa,ngg)       ,! weight opacity integration
C                                          ,!    kappa = cxxx(ia,ig,k)*opa_wgt(ia,ig)*hrod(ir)
     +              cext(naa,ngg,nff)      ,! opacity    C=Qext * surface
     +              cabs(naa,ngg,nff)      ,! opacity    C=Qabs * surface
     +              csca(naa,ngg,nff)      ,! opacity    C=Qsca * surface
     +              qalb(naa,ngg,nff)      ,! dust albedo  Qsca/(Qsca+abs)
     +              gsym(naa,ngg,nff)      ,! assymetry indicator
     +              dust_rho(ngg)          ,! density of (single) dust grain
     +              dust_alpha(ngg)        ,! powerlaw index, Q extrapolation at large lambda
     +              a_dust(naa,ngg)        ,! dust sizes
     +              a_min(ngg)             ,! minimum grain size
     +              a_max(ngg)             ,! maximum grain size
     +              a_delta(ngg)           ,! logarithmic stepsize grain size
     +              a_power(ngg)           ,! power index grainsize distribution function
     +              a_dist(naa,ngg)        ,! grain size distribution function (non-normalized) 
     +              a_norm(ngg)             ! normalization grainsize distribution / rho(ir)
      common/iopac/ na(ngg)                ,! number of grain sizes 
     +              ng                      ! number of grain types
      character*20  dust_name
      character*3   dust_shape
      common/sopac/ dust_name(ngg)         ,! dust name for output purposes
     +              dust_shape(ngg)         ! dust geometrical shape
c
c     Open file containing the frequency grid
c
      inquire(file='frequency.inp',EXIST=fexists)
      if(.not.fexists) then
          write(*,*) 'ERROR: Could not find the file: frequency.inp'
          stop 13
      endif      
c
c     Read the frequency in [Hz]
c
      open(unit=1,file='frequency.inp',status='old')
      read(1,*) nf
      if(nf.gt.FRSIZE_FREQ) then
          write(*,*) 'ERROR: Too many frequencies'
          stop 13
      endif
      do i=1,nf
          read(1,*) freqs(i)
      enddo
      close(1)
c
c     Now compute the wavelength in [microns]
c
      do i=1,nf
          waves(i) = 2.9979d14 / freqs(i)
      enddo
c
c     Read the dust.inp file, which contains the specification of
c     the dust grain species and size distributions
c
      call input_dust_properties()
c
c     Now compute the weights for each species/size, such that 
c     the opacities are normalized as the cross-section for 
c     1 g/cm^3 of that dust species
c
      call normalize_grainsize_distr()
c
c     Now we must write these kappas to files for the transfer
c     program. First the dustopac.inp file
c
      open(unit=1,file='dustopac.inp',status='unknown')
      write(1,705) 1
  705 format(I4,'             Format number of this file')
      write(1,701) ng
  701 format(I4,'             Nr of dust species')
      write(1,710) 
  710 format('===================================================')
      do ig = 1,ng
          write(1,704) 1
  704     format(I4,'             Input style for species ',
     %                '(1=dust kappa file)')
          write(1,702) ig
  702     format(I4,'             Extension of file dustkappa_*.inp')
          write(1,711) 
  711     format('---------------------------------------------------')
      enddo
      write(1,*) ' '
      write(1,*) ' '
      close(1)
c
c     Now the data files
c
      do ig = 1,ng
         if(na(ig).ne.1) stop 888
          if(ig.lt.10) then
              write(ch,'(I1)') ig
              filename = 'dustkappa_'//ch(1:1)//'.inp'
          elseif(ig.lt.100) then
              write(ch,'(I2)') ig
              filename = 'dustkappa_'//ch(1:2)//'.inp'
          elseif(ig.lt.1000) then
              write(ch,'(I3)') ig
              filename = 'dustkappa_'//ch(1:3)//'.inp'
          endif
          open(unit=1,file=filename,status='unknown')
          write(1,*) 2            ! Format number
          write(1,*) nf
          do k = 1,nf
              do ia = 1,na(ig)
                  if(cabs(ia,ig,k)*opa_wgt(ia,ig).le.1.d-20) then
                     cabs(ia,ig,k)=1d-20/opa_wgt(ia,ig)
                  endif
                  if(csca(ia,ig,k)*opa_wgt(ia,ig).le.1.d-20) then
                     csca(ia,ig,k)=1d-20/opa_wgt(ia,ig)
                  endif
                  write(1,110) waves(k),
     %                         cabs(ia,ig,k)*opa_wgt(ia,ig),
     %                         csca(ia,ig,k)*opa_wgt(ia,ig)
              enddo
          enddo
          close(1)
  110     format(E13.6,1x,E13.6,1x,E13.6)
      enddo
c
c     Done
c
      end




c     --------------------------------------------------------------
c       READ THE DUST OPACITY INFORMATION FROM THE DUST DATA FILES
c
c     This routine reads the information about the dust particles 
c     from the file called dust.inp, and computes the opacity data
c     for each frequency. The things that are computed are:
c
c       cext(ia,ig,k) = extinction cross section [cm^2]
c       cabs(ia,ig,k) = absorption cross section [cm^2]
c       csca(ia,ig,k) = scattering cross section [cm^2]
c       qalb(ia,ig,k) = albedo, qsca/(qsca+qabs) = qsca/qext
c       gsym(ia,ig,k) = assymetry  parameter
c        
c     For each dusttype ig is read:
c
c       dust_name()  - name of dust component
c       dust_shape() - geometrical shape of dust component
c       dust_abun()  - dust mass fraction
c       a_min()      - minimum grain size [log micron]
c       a_max()      - maximum grain size [log micron]
c       a_delta()    - log stepsize grain size distribution [log]
c       a_power()    - power law index of grain size distribution function
c       dust_rho()   - density of individual grain 
c       dust_alpha() - powerlaw index, Q extrapolation at large wavelengths
c
c     --------------------------------------------------------------
      subroutine input_dust_properties
      implicit real*8 (a-h,o-z)
c     
      parameter (nlmax = 3000) ! maximum size laboratory wavelength grid
      parameter (naa=DUST_SIZE_MAX)
      parameter (ngg=DUST_SPECIES_MAX)
      parameter (nff=FRSIZE_FREQ)
      parameter (pi=PICONST)
c
      common/freq/waves(nff),freqs(nff)
      common/ifreq/nf
c
      common /opac/ opa_wgt(naa,ngg)       ,! weight opacity integration
C                                          ,!    kappa = cxxx(ia,ig,k)*opa_wgt(ia,ig)*hrod(ir)
     +              cext(naa,ngg,nff)      ,! opacity    C=Qext * surface
     +              cabs(naa,ngg,nff)      ,! opacity    C=Qabs * surface
     +              csca(naa,ngg,nff)      ,! opacity    C=Qsca * surface
     +              qalb(naa,ngg,nff)      ,! dust albedo  Qsca/(Qsca+abs)
     +              gsym(naa,ngg,nff)      ,! assymetry indicator
     +              dust_rho(ngg)          ,! density of (single) dust grain
     +              dust_alpha(ngg)        ,! powerlaw index, Q extrapolation at large lambda
     +              a_dust(naa,ngg)        ,! dust sizes
     +              a_min(ngg)             ,! minimum grain size
     +              a_max(ngg)             ,! maximum grain size
     +              a_delta(ngg)           ,! logarithmic stepsize grain size
     +              a_power(ngg)           ,! power index grainsize distribution function
     +              a_dist(naa,ngg)        ,! grain size distribution function (non-normalized) 
     +              a_norm(ngg)             ! normalization grainsize distribution / rho(ir)
      common/iopac/ na(ngg)                ,! number of grain sizes 
     +              ng                      ! number of grain types
      character*20  dust_name
      character*3   dust_shape
      common/sopac/ dust_name(ngg)         ,! dust name for output purposes
     +              dust_shape(ngg)         ! dust geometrical shape
c
c
      dimension  xlab(nlmax),xn_opt(nlmax),xk_opt(nlmax)
      character*199 cname           ,! dust component name
     +              stripped_cname   
      character*199 dust_file
      character*3   dust_shape0
      character*199 path_dcl
      external len_trim
      double precision :: yinterp
c
c     Set the standard path
c
      path_dcl = 
     +DUST_DATABASE_PATH//
     +DUST_DATABASE_PATHE
c              = path to the dust component library, called in 
c                subroutine input_dust_properties
      n_dcl    = len_trim(path_dcl)
c              = nr of characters in path_dcl
c     
c     Open the main input file
c
      open (unit=30,file='dust.inp',status='unknown')
c
      ig = 0
c        = init number of dustgrain types
c
      write(6,*) 'Making the dust opacities'
c      write (6,*) "selected dust type(s): "
c
 1    read (30,'(a)',end=100) cname
      if (cname(1:1).eq.  '*') goto 1
      if (cname(1:3).eq.'EOF') goto 100
C
      read (30,*) dust_shape0
      read (30,*) status,dust_abun0,a_min0,a_max0,
     +            a_delta0,a_power0,dust_rho0,dust_alpha0
C
      lstr = len_trim(cname)
      if (cname(1:n_dcl).eq.path_dcl(1:n_dcl)) then
C        Standard dustcomponent library
         dust_file = cname(1:lstr)
         stripped_cname = cname(n_dcl+1:lstr)
      else
C        Local directory
c         write (6,*) "  using optical-constants file from local ",
c     %                 "directory"
         dust_file = cname(1:lstr)
         stripped_cname = cname(1:lstr)
      end if
C
      if (status.ne.0) then
C
C        This grain type will be selected
         ig = ig+1
         if(ig.gt.ngg) then
             write(*,*) 'Sorry: compiled for max nr of species = ',ngg
             stop 13
         endif
C
c         write (6,*) '    ',stripped_cname(1:lstr)
         dust_name(ig)    = stripped_cname
C                         = name of dust component
         dust_shape(ig)   = dust_shape0
C                         = geometrical shape of dust component
c         dust_abun(ig)    = dust_abun0
cC                         = dust mass fraction
         a_min(ig)        = a_min0
C                         = minimum grain size [log micron]
         a_max(ig)        = a_max0
C                         = maximum grain size [log micron]
         a_delta(ig)      = a_delta0
C                         = log stepsize grain size distribution [log]
         a_power(ig)      = a_power0
C                         = power law index of grain size distribution function
         dust_rho(ig)     = dust_rho0
C                         = density of individual grain 
         dust_alpha(ig)   = dust_alpha0
C                         = powerlaw index, Q extrapolation at large wavelengths
         if (a_max0.ne.a_min0) then
            if (a_delta(ig).lt.2.d0) then
C              we assume a_delta contains log stepsize grain size distr.
               na(ig)        = (a_max(ig)-a_min(ig))/a_delta(ig) + 1
C                            = number of grain sizes for dust type ig
            else
C              we assume a_delta contains nr of grain sizes for dust type ig
               na(ig)        = a_delta(ig)
               a_delta(ig)   = (a_max(ig)-a_min(ig))/(na(ig)-1)
            end if
         else 
            na(ig)        = 1
         end if
         if (na(ig).gt.naa) stop 'dust_properties: na(ig).gt.naa'
C
         do ia = 1,na(ig)
            a_dust(ia,ig) = 10**(a_min(ig)+(ia-1)*a_delta(ig))
C                         = grain sizes in micron
            a_dist(ia,ig) = (a_dust(ia,ig)/10**a_min0)**a_power(ig)
C                         = grainsize distribution function (non-normalized)
         end do              
C
C        -----------------------------------------------------------------
         if (status.eq.1) then ! q's are calculated from optical constants'
C        -----------------------------------------------------------------
C
C        Read optical constants for current grain component
C
C        ---------------- start reading optical data file
C
         open (unit=31,file=dust_file,status='old',iostat=ierr)
         if (ierr.ne.0) then
            write (6,*) ' Error opening file ',dust_file
            stop 'input_dust_properties: error opening file'
         end if
C
         klab = 1
 2       read (31,*,end=10) xlab(klab),xn_opt(klab),xk_opt(klab)
         klab = klab+1
         goto 2
 10      nlab = klab-1
         xmin = xlab(   1)
         xmax = xlab(nlab)
C
         close (unit=31)
C
C        ---------------- end   reading optical data file
C
C        Interpolate/extrapolate optical constants for all grain 
C        sizes from laboratory wavelength grid to continuum wavelength 
C        grid on which emergent spectrum will be calculated
C
C        To facilitate extrapolation:
         constant_n = (1.d0-xn_opt(1))/xmin/xmin
         constant_k =  xk_opt(1)/xmin/xmin/xmin
         if (waves(1).lt.xmin) write (6,200) xmin
 200     format (9x,'extrapolation optical constants <',
     +           0pf6.2,' micron')
C
         if(waves(nf).gt.waves(1)) then
             kkmin  = 1
             kkmax  = nf
             kkstep = 1
         else
             kkmin  = nf
             kkmax  = 1
             kkstep = -1
         endif 
         do k = kkmin,kkmax,kkstep
            waves0 = waves(k)       ! [micron]
C 
            if      (waves0.lt.xmax) then
               if (waves0.lt.xmin) then
C                 extrapolate k,n for short-wavelengths following
C                 Bohren & Huffman, p. 234 & p. 254
                  xns0 = 1.d0-constant_n*waves0*waves0
                  xks0 = constant_k*waves0*waves0*waves0
               else
C                 interpolate k,n from laboratory wavelength grid            
C                  call interp (xlab,xn_opt,nlab,3,waves0,xns0)
C                  call interp (xlab,xk_opt,nlab,3,waves0,xks0)
                  xns0 = yinterp(waves0,xlab,xn_opt,nlab,nlmax,0)
                  xks0 = yinterp(waves0,xlab,xk_opt,nlab,nlmax,0)
               end if
C
               do ia = 1,na(ig)
                  a0_cm         = 1.d-4*a_dust(ia,ig)
C                               = grain radius in cm
                  a_surface     = pi*a0_cm*a0_cm
C                               = grain disk surface in cm^2
C
C                 Select geometrical shape mode
                  if(dust_shape(ig).eq."MIE") then
C                       Mie theory
                        ttt = 0.d0
                        call q_mie(xns0,xks0,waves0,a_dust(ia,ig),
     %                             ttt,qe,qs,qa,ge)
                        ce = abs(qe) * a_surface
                        ca =     qa  * a_surface
                        cs =     qs  * a_surface
                    else if(dust_shape(ig).eq."CDE") then
C                       Continuous distribution of elipsoids at constant volume
                        call q_cde(xns0,xks0,waves0,a_dust(ia,ig),
     %                             ce,cs,ca,ge)
C                       convert cross section from micron^2 -> cm^2
                        ce = ce * 1.d-8
                        cs = cs * 1.d-8
                        ca = ca * 1.d-8
                    else
                        write (6,*) " non-valid geometrical mode"
                        stop
cccc                  end select
                    endif
C
                  cext(ia,ig,k) = ce
C                               = extinction cross section [cm^2]
                  cabs(ia,ig,k) = ca
C                               = absorption cross section [cm^2]
                  csca(ia,ig,k) = cs
C                               = scattering cross section [cm^2]
                  qalb(ia,ig,k) = cs/ce
C                               = albedo, qsca/(qsca+qabs) = qsca/qext
                  gsym(ia,ig,k) = abs(ge)
C                               = assymetry  parameter
               end do
            else
C              extrapolate k,n for long-wavelengths
               factor = (waves0/waves(k-kkstep))**dust_alpha(ig)
               do ia = 1,na(ig)
                  cext(ia,ig,k) = cext(ia,ig,k-kkstep)*factor
                  qalb(ia,ig,k) = qalb(ia,ig,k-kkstep)*factor
                  cabs(ia,ig,k) = (1.d0-qalb(ia,ig,k))*cext(ia,ig,k)
                  csca(ia,ig,k) = qalb(ia,ig,k) *cext(ia,ig,k)
                  gsym(ia,ig,k) = gsym(ia,ig,k-kkstep)*factor
               end do
            end if
         end do
C
C        -------------------------------------------
         else if (status.eq.-2) then ! benchmark q's'
C        -------------------------------------------
C
C        Benchmark absorption/scattering cross sections as given in 
C        Ivezic et al. 1997, MNRAS, 291, 121
C
         dust_name(ig) = 'benchmark-2'
C
         x1 = one
C           = normalizing wavelength [micron]
         do ia = 1,na(ig)
            a0_cm     = 1.d-4*a_dust(ia,ig)
C                     = grain radius in cm
            a_surface = pi*a0_cm*a0_cm
C                     = grain disk surface in cm^2
            x1 = 2.d0*pi*a_dust(ia,ig)
C              = limit of a >> lambda is chosen at           a = lambda
C                note: one may also define the edge at  2*pi*a = lambda
            do k = 1,nf
               waves0 = waves(k)       ! [micron]
               if (waves0.le.x1) then
                  cabs(ia,ig,k) = one * a_surface
                  csca(ia,ig,k) = one * a_surface
                  cext(ia,ig,k) = two * a_surface
               else
                  cabs(ia,ig,k) = x1/waves0 * a_surface
                  csca(ia,ig,k) = x1/waves0**4 * a_surface
                  cext(ia,ig,k) = (x1/waves0+x1/waves0**4) * 
     %                            a_surface 
               end if
            end do
         end do
C
C        ----------------------------------------------------
         end if ! end read options
C        ----------------------------------------------------
C
      end if
C
      goto 1
 100  ng = ig
C        = total number of grain components
C
      return
      end




c     --------------------------------------------------------------
c                 COMPUTE THE OPACITY INTEGRATION WEIGHTS
c     --------------------------------------------------------------
      subroutine normalize_grainsize_distr
      implicit real*8 (a-h,o-z)
c
      parameter (naa=DUST_SIZE_MAX)
      parameter (ngg=DUST_SPECIES_MAX)
      parameter (nff=FRSIZE_FREQ)
      parameter (pi=PICONST)
c
      common /opac/ opa_wgt(naa,ngg)       ,! weight opacity integration
C                                          ,!    kappa = cxxx(ia,ig,k)*opa_wgt(ia,ig)*hrod(ir)
     +              cext(naa,ngg,nff)      ,! opacity    C=Qext * surface
     +              cabs(naa,ngg,nff)      ,! opacity    C=Qabs * surface
     +              csca(naa,ngg,nff)      ,! opacity    C=Qsca * surface
     +              qalb(naa,ngg,nff)      ,! dust albedo  Qsca/(Qsca+abs)
     +              gsym(naa,ngg,nff)      ,! assymetry indicator
     +              dust_rho(ngg)          ,! density of (single) dust grain
     +              dust_alpha(ngg)        ,! powerlaw index, Q extrapolation at large lambda
     +              a_dust(naa,ngg)        ,! dust sizes
     +              a_min(ngg)             ,! minimum grain size
     +              a_max(ngg)             ,! maximum grain size
     +              a_delta(ngg)           ,! logarithmic stepsize grain size
     +              a_power(ngg)           ,! power index grainsize distribution function
     +              a_dist(naa,ngg)        ,! grain size distribution function (non-normalized) 
     +              a_norm(ngg)             ! normalization grainsize distribution / rho(ir)
      common/iopac/ na(ngg)                ,! number of grain sizes 
     +              ng                      ! number of grain types
      character*20  dust_name
      character*3   dust_shape
      common/sopac/ dust_name(ngg)         ,! dust name for output purposes
     +              dust_shape(ngg)         ! dust geometrical shape
c
      dimension b_dust(naa),b_wgt(naa),c_wgt(naa),d_wgt(naa)
C
C     Normalization grainsize distribution function
C
      do ig = 1,ng
         if (na(ig).eq.1) then
C           only one dust grainsize
            ia            = 1
            a0_cm         = 1.d-4*a_dust(ia,ig)
C                         = grainsize in cm
            arg           = 4.d0/3.d0*pi*dust_rho(ig)*a_dust(ia,ig)**3
            a_norm(ig)    = 1.d0/arg
C                         = normalization grainsize distribution [/rhod(ir)]
            arg0          = 4.d0/3.d0*pi*dust_rho(ig)*a0_cm*a0_cm*a0_cm
C                         = mass of single dust grain [gr]
            opa_wgt(ia,ig)= 1.d0/arg0
C                         = weight opacity integration
C                           kappa = c(ia,ig,k)*opa_wgt(ia,ig)*rhod(ir)
         else
C           multiple grain sizes
            do ia = 1,na(ig)
               b_dust(ia) = 1.d-4*a_dust(ia,ig)
C                         = store grainsizes in cm in single dimension array
            end do
C
            pas    = a_power(ig)
C                  = power-law index of grainsize distribution function
            pasneg = -pas
            pam    = pas + 3.d0
C                  = power-law index of grainmass distribution function
            pao    = pas + 2.d0
C                  = power-law index of opacity grainsize integration 
C
C        -  weights needed to determine the normalization of the grainsize
C           distribution (such that the total dust density is recovered)
C
            call momentweight (b_dust,b_wgt, pam,na(ig),sumb_wgt)
cC
cC        -  weights needed to determine the total number of grain particles
cC          
c            call momentweight (b_dust,c_wgt, pas,na(ig),sumc_wgt)
C
C        -  weights needed in opacity integration over all grainsizes
C
            call momentweight (b_dust,d_wgt, pao,na(ig),sumd_wgt)
C
            arg0 = 4.d0/3.d0*pi*dust_rho(ig)*b_dust(1)**pasneg*sumb_wgt
            a_norm(ig) = 1.d0/arg0
C                      = normalization grainsize distribution [/rhod(ir)]
C                        correct normalization is:  a_norm(ig)*rhod(ir)  
C
            do ia = 1,na(ig)
               a0_cm         = 1.d-4*a_dust(ia,ig)
C                            = grainsize in cm
               a_surface     = pi*a0_cm*a0_cm
C                            = grain disk surface in cm^2
               opa_wgt(ia,ig)= pi*a_norm(ig)*b_dust(1)**pasneg*
     +                         d_wgt(ia)/a_surface
C                            = weights needed in integrating opacity
C                              contributions of all grainsizes
C                              kappa = c(ia,ig,k)*opa_wgt(ia,ig)*rhod(ir)
            end do                       
         end if
      end do
C
      return
      end





      function len_trim (card)
      integer len_trim
C     =========================
C
C     Trims a string of trailing blanks
C        Function works only if card is character*80 or larger in 
C        calling routine!
C
      character*180 card
C      
      ncard = len(card)
      do icard = len(card),1,-1
         if (card(icard:icard).eq." " .or.
     +       card(icard:icard).eq.""       ) then
            ncard = ncard-1
         else
            len_trim = ncard
            return
         endif
      end do
C
      return
      end
C
C    
C     ******************************************************************
C
C
      subroutine momentweight (x,wgt,m,n,s)
      implicit real*8 (a-h,o-z)
C     =======================
C
C     Setup of weights of the m'th moment integral using trapezium
C     rule. Weights are in general defined as:
C
C        xmax                n
C        Int  f(x)x^m dx =  Sum  f(i)*wgt(i)   
C        xmin               i=1
C
C     where,  wgt(i) = d(x)              if m ==  0
C                    = ln(x)                m == -1
C                    = d(x^[m+1]/[m+1])     m <>  0,-1
C
C     The sum of the weights is returned in s
C
c      include './include/implic.dust'
      real*8 m,mp
      dimension x(n),wgt(n),y(32000)
C
C     Error checking
      if (n.gt.32000) stop 'momentweight: n.gt.32000'
C
      mp = m+1
      do i = 1,n
         if (x(i).eq.0.d0.and.m.lt.0) stop 'momentweight: div by 0'
         if      (m.eq. 0) then
            y(i) = x(i)
         else if (m.eq.-1) then
            y(i) = log(x(i))
         else
            y(i) = (x(i)**mp)/mp
         end if
      end do
C
      s = 0.d0
      if (n.eq.1) then
         wgt(1)    = 0.0d0
      else
         wgt(1)    = 0.5d0*dabs(y(  2)-y(  1))
         s         = s+wgt(1)
         do i = 2,n-1
            wgt(i) = 0.5d0*dabs(y(i+1)-y(i-1))
            s      = s+wgt(i)
         end do
         wgt(n)    = 0.5d0*dabs(y(  n)-y(n-1))
         s         = s+wgt(n)
      end if
C
      return
      end
C
C
C     ******************************************************************
C
C
      real*8 function yinterp (xint,x,y,n,ndim,mode)
      implicit real*8 (a-h,o-z)
C     =======================
C
C     Interpolation routine. Determines YINT = Y(XINT) from 
C     grid Y(X) with N points and dimension NDIM. 
C
C     mode = 0  -  interpolation: linear      in x, linear      in y
C            1  -                 linear      in x, logarithmic in y
C            2  -                 logarithmic in x, logarithmic in y
C
c      include './include/implic.dust'
c      include './include/consts.dust'
      dimension x(ndim),y(ndim)
C
C     bisection (see Numerical Recipes par 3.4 page 90)
      jl = 0 
      ju = n+1
10    if (ju-jl.gt.1) then
         jm = (ju+jl)/2
         if ((x(n).gt.x(1)).eqv.(xint.gt.x(jm))) then
            jl = jm
         else
            ju = jm
         end if
         go to 10
      end if
      j = min(jl,n-1)
      j = max(1,j)
C
C     interpolation
C
      if      (mode.eq.0) then
         slope   = (y(j+1)-y(j))/(x(j+1)-x(j))
         yinterp = slope*(xint-x(j))+y(j)
         return
      else if (mode.eq.1) then
         if (y(j+1).le.zero.or.y(j).le.zero) then
            yinterp = zero
            return
         else
            slope   = (log(y(j+1))-log(y(j)))/(x(j+1)-x(j))
            yinterp = slope*(xint-x(j))+log(y(j))
            if      (yinterp.gt. expmax) then
               yinterp = exp( expmax)
            else if (yinterp.lt.-expmax) then
               yinterp = exp(-expmax)
            else
               yinterp = exp(yinterp)
            end if
         end if
      else if (mode.eq.2) then
         if (y(j+1).le.zero.or.y(j).le.zero.or.
     +       x(j+1).le.zero.or.x(j).le.zero) then
            yinterp = zero
            return
         else
            slope   = (log(y(j+1))-log(y(j)))/(log(x(j+1))-log(x(j)))
            yinterp = slope*(log(xint)-log(x(j)))+log(y(j))
            if      (yinterp.gt. expmax) then
               yinterp = exp( expmax)
            else if (yinterp.lt.-expmax) then
               yinterp = exp(-expmax)
            else
               yinterp = exp(yinterp)
            end if
         end if
      else
         write (6,*) 'yinterp: unknown mode = ',mode
         stop
      end if
C
      return
      end
