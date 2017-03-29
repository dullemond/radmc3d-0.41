!======================================================================
!             PAH / VSG MODULE: QUANTUM-HEATED GRAINS
!
! This is a port of a code that was originally made for use with the
! original RADMC code. That's why there are comments from before the
! RADMC-3D code in there. 
!
! A main difference is that for the UV photons (the ones that excite
! the PAHs) we here use an independent frequency array, while the
! cooling is done using the standard frequency array freq_nu(). 
! This allows us to fine-tune the wavelengths of the PAH-exciting
! photons. This UV-photon array is quantum_frequencies(). 
!======================================================================
module quantum_module
use rtglobal_module
use dust_module
use constants_module

doubleprecision,allocatable :: emisquant(:,:)
doubleprecision,allocatable :: emisquant_loccum(:,:)
doubleprecision,allocatable :: emisquant_loctot(:)
doubleprecision,allocatable :: emisquant_cum(:)
doubleprecision :: emisquanttot
doubleprecision,allocatable :: quantum_frequencies(:)
doubleprecision,allocatable :: quantum_dnu(:)
integer :: quantum_nf
doubleprecision,allocatable :: quantum_temp_grid(:)
doubleprecision,allocatable :: quantum_temp_distr(:,:,:)
doubleprecision,allocatable :: quantum_mat(:,:)
doubleprecision,allocatable :: quantum_cooltime(:,:)
doubleprecision,allocatable :: quantum_pahdes_time(:)
doubleprecision,allocatable :: quantum_pahdes_pcrit(:)
integer,allocatable :: quantum_peak_itemp(:,:,:)
doubleprecision,allocatable :: quantum_peak_eps(:,:,:)
doubleprecision,allocatable :: quantum_peak_temp(:,:,:)
integer,allocatable :: quantum_ispec(:)
integer :: quantum_nrquantum=0
integer :: quantum_ntemp=100
doubleprecision :: quantum_temp0=1d0
doubleprecision :: quantum_temp1=1d4
doubleprecision,allocatable :: quantum_kappa(:,:)

contains

!---------------------------------------------------------------------
!                 Initialize the quantum mode
!---------------------------------------------------------------------
subroutine quantum_init(set_opacities,set_tempgrid,incl_emisquant)
implicit none
integer :: itemp,ispec,isq,inu
logical :: set_opacities,incl_emisquant,set_tempgrid
!
! Message and check
!
write(stdo,*) 'Initializing quantum module'
if(nrcellsmax.le.0) then
   write(stdo,*) 'INTERNAL ERROR IN QUANTUM_MODULE: MONTE CARLO NOT INITIALIZED'
   stop
endif
if(dust_nr_species.le.0) then
   write(stdo,*) 'Error in quantum-heating module: Dust opacities not yet read...'
   stop
endif
!
! Find which of the dust species are quantum
!
if(allocated(quantum_ispec)) deallocate(quantum_ispec)
quantum_nrquantum = 0
do ispec=1,dust_nr_species
   if(dust_quantum(ispec).ne.0) then
      quantum_nrquantum = quantum_nrquantum + 1
   endif
enddo
allocate(quantum_ispec(quantum_nrquantum))
isq = 0
do ispec=1,dust_nr_species
   if(dust_quantum(ispec).ne.0) then
      isq = isq + 1
      quantum_ispec(isq) = ispec
   endif
enddo
if(isq.ne.quantum_nrquantum) stop 123
!
! Get the opacities for the quantum wavelength grid
!
if(set_opacities) then
   !
   ! Allocate the opacity arrays
   !
   if(allocated(quantum_kappa)) deallocate(quantum_kappa)
   allocate(quantum_kappa(quantum_nf,quantum_nrquantum))
   !
   ! Fill these opacity arrays
   ! 
   do isq=1,quantum_nrquantum
      ispec = quantum_ispec(isq)
      do inu=1,quantum_nf
         quantum_kappa(inu,isq) = find_dust_kappa_interpol(  &
              quantum_frequencies(inu),ispec,100.d0,1,0,0)
      enddo
   enddo
endif
!
! Allocate the temperature grid. Note that the parameters can be
! fine-tuned in radmc3d.inp
!
if(set_tempgrid) then
   !
   ! Cleanup stuff before allocating
   !
   if(allocated(quantum_temp_grid)) deallocate(quantum_temp_grid)
   if(allocated(quantum_temp_distr)) deallocate(quantum_temp_distr)
   if(allocated(quantum_mat)) deallocate(quantum_mat)
   if(allocated(quantum_cooltime)) deallocate(quantum_cooltime)
   if(allocated(quantum_pahdes_time)) deallocate(quantum_pahdes_time)
   if(allocated(quantum_pahdes_pcrit)) deallocate(quantum_pahdes_pcrit)
   if(allocated(quantum_peak_itemp)) deallocate(quantum_peak_itemp)
   if(allocated(quantum_peak_eps)) deallocate(quantum_peak_eps)
   if(allocated(quantum_peak_temp)) deallocate(quantum_peak_temp)
   !
   ! The temperature grid itself
   !
   allocate(quantum_temp_grid(1:quantum_ntemp))
   do itemp=1,quantum_ntemp
      quantum_temp_grid(itemp) = quantum_temp0 *     &
           (quantum_temp1/quantum_temp0)**      &
           ((itemp-1.d0)/(quantum_ntemp-1.d0))
   enddo
   !
   ! Allocate the matrix required for the temperature
   ! distribution solution
   !
   allocate(quantum_mat(quantum_ntemp,quantum_ntemp))
   !
   ! Allocate the cooling time arrays
   !
   allocate(quantum_cooltime(quantum_ntemp,quantum_nrquantum))
   !
   ! Allocate the array for the temperature distribution
   ! function, but only for the quantum grains...
   !
   allocate(quantum_temp_distr(quantum_ntemp,quantum_nrquantum,nrcellsmax))
   !
   ! Allocate the PAH destruction arrays
   !
   allocate(quantum_pahdes_time(quantum_nrquantum))
   allocate(quantum_pahdes_pcrit(quantum_nrquantum))
   !
   ! Allocate the peak temperature arrays
   !
   allocate(quantum_peak_itemp(quantum_ntemp,quantum_nf,quantum_nrquantum))
   allocate(quantum_peak_eps(quantum_ntemp,quantum_nf,quantum_nrquantum))
   allocate(quantum_peak_temp(quantum_ntemp,quantum_nf,quantum_nrquantum))
endif
!
! Allocate the arrays necessary to include the quantum
! emission as a source in the thermal Monte Carlo.
! This is more self-consistent than simply computing
! the quantum emission and assuming the quantum emission
! does not affect the other (thermal) grains. But it costs
! a lot of memory, especially if a lot of frequencies are
! used in the global frequency array.
!
if(incl_emisquant) then
   allocate(emisquant(1:freq_nr,1:nrcellsmax))
   allocate(emisquant_loccum(1:freq_nr+1,1:nrcellsmax))
   allocate(emisquant_loctot(1:nrcellsmax))
   allocate(emisquant_cum(nrcellsmax+1))
endif
end subroutine quantum_init


!---------------------------------------------------------------------
!                   Cleanup the quantum mode
!---------------------------------------------------------------------
subroutine quantum_cleanup()
implicit none
if(allocated(emisquant)) deallocate(emisquant)
if(allocated(emisquant_loccum)) deallocate(emisquant_loccum)
if(allocated(emisquant_loctot)) deallocate(emisquant_loctot)
if(allocated(emisquant_cum)) deallocate(emisquant_cum)
if(allocated(quantum_frequencies)) deallocate(quantum_frequencies)
if(allocated(quantum_dnu)) deallocate(quantum_dnu)
if(allocated(quantum_temp_grid)) deallocate(quantum_temp_grid)
if(allocated(quantum_temp_distr)) deallocate(quantum_temp_distr)
if(allocated(quantum_mat)) deallocate(quantum_mat)
if(allocated(quantum_cooltime)) deallocate(quantum_cooltime)
if(allocated(quantum_ispec)) deallocate(quantum_ispec)
if(allocated(quantum_kappa)) deallocate(quantum_kappa)
if(allocated(quantum_pahdes_time)) deallocate(quantum_pahdes_time)
if(allocated(quantum_pahdes_pcrit)) deallocate(quantum_pahdes_pcrit)
if(allocated(quantum_peak_itemp)) deallocate(quantum_peak_itemp)
if(allocated(quantum_peak_eps)) deallocate(quantum_peak_eps)
if(allocated(quantum_peak_temp)) deallocate(quantum_peak_temp)
end subroutine quantum_cleanup



!---------------------------------------------------------------------
!        READ THE WAVELENGTH BINS FOR QUANTUM-CAPABLE PHOTONS
!
! Since quantum-heating of PAHs or other small carbonaceous grains
! is driven by UV photons, but not by infrared photons, it is 
! useful to have a separate wavelength grid for these UV photons.
!---------------------------------------------------------------------
subroutine quantum_read_wavelengths(action)
  implicit none
  integer :: action,inu
  logical :: fex
  doubleprecision, allocatable :: freq(:)
  !
  ! Action=0 means do nothing, action=1 means read if not yet read,
  ! action=2 means re-read.
  !
  if(action.eq.0) then
     return
  elseif(action.eq.1) then
     if(allocated(quantum_frequencies)) return
  elseif(action.eq.2) then
     if(allocated(quantum_frequencies)) deallocate(quantum_frequencies)
     if(allocated(quantum_dnu)) deallocate(quantum_dnu)
  endif
  !
  ! Check which file is present
  !
  inquire(file='quantum_wavelength_micron.inp',exist=fex)
  if(.not.fex) then
     write(stdo,*) 'ERROR: Could not find quantum_wavelength_micron.inp ', &
          'necessary for quantum heated grains.'
     stop
  endif
  !
  ! Message
  !
  write(stdo,*) 'Reading quantum frequencies/wavelengths...'
  call flush(stdo)
  !
  ! Read this file
  !
  open(unit=1,file='quantum_wavelength_micron.inp')
  read(1,*) quantum_nf
  if(quantum_nf.le.1) then
     write(stdo,*) 'ERROR in reading quantum_wavelength_micron.inp.'
     write(stdo,*) '   You must specify at least 2 wavelengths (in micron), so that '
     write(stdo,*) '   the size of the frequency bin can be determined. This is necessary '
     write(stdo,*) '   to be able to determine the total luminosity of PAH-exciting radiation.'
     stop
  endif
  !
  ! Allocate a temporary array
  !
  allocate(freq(quantum_nf))
  !
  ! Now read the wavelengths-in-microns
  !
  do inu=1,quantum_nf
     read(1,*) freq(inu)
     freq(inu) = 1d4*cc / freq(inu) 
  enddo
  close(1)
  !
  ! The actually used frequency array is 1 bin less
  !
  quantum_nf = quantum_nf - 1
  !
  ! Allocate the frequency arrays
  !
  allocate(quantum_frequencies(1:quantum_nf))
  allocate(quantum_dnu(1:quantum_nf))
  !
  ! Compute the average frequencies and the bin widths
  !
  do inu=1,quantum_nf
     quantum_frequencies(inu) = 0.5d0*(freq(inu+1)+freq(inu))
     quantum_dnu(inu)         = abs(freq(inu+1)-freq(inu))
  enddo
  !
  ! Deallocate the temporary array
  !
  deallocate(freq)
end subroutine quantum_read_wavelengths


!-----------------------------------------------------------------
!            SOLVE FOR THE TEMPERATURE DISTRIBUTION
!-----------------------------------------------------------------
subroutine quantum_solve_temp_distrib(isq,nf,ntemp,freq,meanint,tdist)
  implicit none
  integer :: nf,ntemp
  doubleprecision :: freq(nf),meanint(nf),dummy,tdist(ntemp)
  doubleprecision :: matbk(quantum_ntemp)
  doubleprecision :: matwork(quantum_ntemp,quantum_ntemp)
  doubleprecision :: tdtmp(quantum_ntemp),res(quantum_ntemp)
  integer :: isq,ispec,ievap,itemp,it
  integer :: indx(quantum_ntemp)
  !
  ! Check
  !
  if(ntemp.ne.quantum_ntemp) stop 2317
  !
  ! Get the dust ispec from isq
  !
  if(isq.gt.quantum_nrquantum) stop 8310
  ispec = quantum_ispec(isq)
  !
  ! Fill the matrix
  !
  call quantum_fill_matrix(isq,quantum_nf,quantum_ntemp,meanint)
  !
  ! This matrix obviously has one zero eigenvalue because all the
  ! temperature bins are converted into each other. This cannot be
  ! inverted.
  !
  ! So replace last line with normalization
  ! But remember the original line for later in this routine.
  !
  do itemp=1,quantum_ntemp
     matbk(itemp)                     = quantum_mat(quantum_ntemp,itemp)
     quantum_mat(quantum_ntemp,itemp) = 1.d0
  enddo
  !
  ! Put (0,0,...,0,1) as right-hand size
  !
  do itemp=1,quantum_ntemp-1
     tdist(itemp) = 0.d0
  enddo
  tdist(quantum_ntemp) = 1.d0
  !
  ! Precondition this matrix 
  !
  do itemp=1,quantum_ntemp-1
     dummy = 0.d0
     do it=1,quantum_ntemp
        if(abs(quantum_mat(itemp,it)).gt.dummy) dummy=abs(quantum_mat(itemp,it))
     enddo
     do it=1,quantum_ntemp
        quantum_mat(itemp,it) = quantum_mat(itemp,it) / dummy
     enddo
  enddo
  !
  ! Solve this matrix equation
  !
  dummy = 1.d0
  do itemp=1,quantum_ntemp
     do it=1,quantum_ntemp
        matwork(itemp,it) = quantum_mat(itemp,it)
     enddo
  enddo
  call ludcmp(matwork,quantum_ntemp,quantum_ntemp,indx,dummy)
  call lubksb(matwork,quantum_ntemp,quantum_ntemp,indx,tdist)
  !
  ! Solution has been obtained...
  !
  ! Check if it is indeed a solution
  !
  dummy = 0.d0
  do itemp=1,quantum_ntemp
     tdtmp(itemp)  = 0.d0
     do it=1,quantum_ntemp
        tdtmp(itemp) = tdtmp(itemp) + quantum_mat(itemp,it) * tdist(it)
     enddo
     if(itemp.lt.quantum_ntemp) then
        dummy = dummy + abs(tdtmp(itemp))
     else
        dummy = dummy + abs(tdtmp(itemp)-1.d0)
     endif
  enddo
  if(dummy.gt.1d-10) then
     write(stdo,*) 'ERROR in matix inversion quantum: Error = ',dummy
     stop
  endif
  !
  ! Check if the solution makes any sense...
  !
  dummy = 0.d0
  do itemp=1,quantum_ntemp
     dummy = dummy + tdist(itemp)
  enddo
  if(abs(dummy-1.d0).gt.1d-7) then 
     write(stdo,*) 'ERROR: Distribution does not add up to 1,'
     write(stdo,*) (tdist(itemp),itemp=1,quantum_ntemp)
     write(stdo,*) dummy
     stop 9564
  endif
  !
  ! Some values may be <0 due to numerical errors
  ! Put to 0 those values
  !
  do itemp=1,quantum_ntemp
     if(tdist(itemp).lt.0.d0) then
        tdist(itemp) = 0.d0
     endif
  enddo
  !
  ! Check again if the norm is not too much compromised
  !
  dummy = 0.d0
  do itemp=1,quantum_ntemp
     dummy = dummy + tdist(itemp)
  enddo
  if(abs(dummy-1.d0).gt.1d-2) then 
     write(stdo,*) 'ERROR: Distribution does not add up to 1,'
     write(stdo,*) (tdist(itemp),itemp=1,quantum_ntemp)
     write(stdo,*) dummy
     open(unit=2,file='matrix.dat')
     write(2,*) quantum_ntemp
     do itemp=1,quantum_ntemp
        do it=1,quantum_ntemp
           write(2,*) quantum_mat(itemp,it)
        enddo
     enddo
     close(2)
     stop 9565
  endif
  !
  ! Renormalize, to make sure that the norm is exactly 1
  !
  do itemp=1,quantum_ntemp
     tdist(itemp) = tdist(itemp) / dummy
  enddo
  !    
  ! The next part of this routine is only for computing
  ! the typical life time of the PAHs
  !
  ! Compute the typical life time of the grains
  !   
  if(dust_tmax(ispec).gt.0.d0) then
     if(quantum_pahdes_pcrit(isq).eq.0.d0) then
        !
        ! The original method of Kees
        !
        ! Find the index of the temperature grid 
        !
        call hunt(quantum_temp_grid,quantum_ntemp,dust_tmax(ispec),ievap)
        if(ievap.ge.quantum_ntemp) stop 90934
        if(ievap.lt.1) then
           write(stdo,*) dust_tmax(ispec)
           stop 90935
        endif
        !      
        ! Select all stuff below this temperature
        !     
        do itemp=1,ievap
           tdtmp(itemp) = tdist(itemp)
        enddo
        do itemp=ievap+1,quantum_ntemp
           tdtmp(itemp) = 0.d0
        enddo
        !
        ! Normalize
        !
        dummy = 0.d0
        do itemp=1,quantum_ntemp
           dummy = dummy + tdtmp(itemp)
        enddo
        do itemp=1,quantum_ntemp
           tdtmp(itemp) = tdtmp(itemp) / dummy
        enddo
        !
        ! Restore last line of the matrix
        !             
        do itemp=1,quantum_ntemp
           quantum_mat(quantum_ntemp,itemp) = matbk(itemp)
        enddo
        !             
        ! Do matrix multip
        !
        do itemp=1,quantum_ntemp
           res(itemp) = 0.d0
           do it=1,quantum_ntemp
              res(itemp) = res(itemp) + quantum_mat(itemp,it) * tdtmp(it)
           enddo
        enddo
        !
        ! Check out how much has exceeded the tevap
        ! BUGFIX-13-07-5: use res(itemp) instead of tdtmp(itemp)
        !
        dummy = 0.d0
        do itemp=ievap+1,quantum_ntemp
           dummy = dummy + res(itemp)
        enddo
        quantum_pahdes_time(isq) = 1.d0 / ( abs(dummy) + 1d-99 )
     else
        !
        ! The Habart method, slightly generalized
        ! [ADDED: 25.11.06]
        !
        ! Find the index of the temperature grid 
        !
        call hunt(quantum_temp_grid,quantum_ntemp,dust_tmax(ispec),ievap)
        if(ievap.ge.quantum_ntemp) stop 90934
        if(ievap.lt.1) then
           write(stdo,*) dust_tmax(ispec)
           stop 90935
        endif
        !
        ! Do the integral from T_crit to infty
        !
        dummy = 0.d0
        do itemp=ievap+1,quantum_ntemp
           dummy = dummy + tdist(itemp)
        enddo
        !
        ! For Habart, if this exceeds 1e-8, she destroys it
        ! Here I assume this 1e-8 to correspond to 1 Myr 
        ! (any correction to this assumtion is reflected in 
        ! another value of pahdes_pcrit read in from the dust
        ! opacity files).
        !
        quantum_pahdes_time(isq) = 3.1536d+13 * ( quantum_pahdes_pcrit(isq) / dummy )
        !
     endif
  else
     quantum_pahdes_time(isq) = 1d99
  endif
  !
end subroutine quantum_solve_temp_distrib


!-----------------------------------------------------------------
!          FILL THE MATRIX, FOR A GIVEN RADIATION FIELD
!-----------------------------------------------------------------
subroutine quantum_fill_matrix(isq,nf,ntemp,meanint)
  implicit none
  integer :: isq,nf,ntemp
  doubleprecision :: meanint(nf)
  doubleprecision :: epspeak,nphot,dum,dt
  integer :: inu,itemp,itpeak,ispec,it
  !
  ! Get the dust ispec from isq
  !
  if(isq.gt.quantum_nrquantum) stop 8310
  ispec = quantum_ispec(isq)
  !
  ! First clear the matrix
  !
  do itemp=1,quantum_ntemp
     do it=1,quantum_ntemp
        quantum_mat(itemp,it) = 0.d0
     enddo
  enddo
  !
  ! Fill matrix for the cooling
  !
  do itemp=quantum_ntemp,2,-1
     !
     ! Compute the time it takes to cool from itemp to itemp-1
     !
     ! We do this in an explicit way, i.e. the time scale for
     ! cooling from itemp to itemp-1 is computed using the 
     ! infrared-emission cooling rate at itemp. In this way
     ! we make sure that the shift of the material to lower 
     ! temperature over an explicit time step delta t corresponds
     ! to the cooling rate at the beginning of the time step
     ! multiplied by delta t.
     !
     dt = quantum_cooltime(itemp-1,isq) - quantum_cooltime(itemp,isq)
     !
     ! Compute the matrix contribution
     !
     dum = 1.d0 / dt
     !
     ! Add this to the matrix, such that mass is conserved
     !
     quantum_mat(itemp,itemp)   = quantum_mat(itemp,itemp)   - dum
     quantum_mat(itemp-1,itemp) = quantum_mat(itemp-1,itemp) + dum
  enddo
  !
  ! Fill the matrix for the heating
  !
  do itemp=1,quantum_ntemp-1
     !
     ! Excite from temperature level itemp to the temperature
     ! level given by the input photon energy. So make a loop
     ! over all incoming photons
     !
     do inu=1,quantum_nf
        !
        ! Find where the energy is dumped
        !
        itpeak  = quantum_peak_itemp(itemp,inu,isq)
        epspeak = quantum_peak_eps(itemp,inu,isq)
        !
        ! Find out how many photons are absorbed per second
        ! by 1 gram of PAHs 
        !
        nphot = 4 * pi * meanint(inu) * quantum_dnu(inu) *   &
                quantum_kappa(inu,isq) /                     &
                ( hh * quantum_frequencies(inu) )
        !
        ! Convert this into how many photons are absorbed per
        ! second by 1 PAH molecule
        !
        nphot = nphot * dust_mgrain(ispec)
        !
        ! This is the excitation rate...
        !
        ! Put this in the matrix
        !
        if(itpeak.lt.quantum_ntemp) then
           !
           ! Normal case
           !
           quantum_mat(itemp,itemp)    = quantum_mat(itemp,itemp) - nphot
           quantum_mat(itpeak,itemp)   = quantum_mat(itpeak,itemp) + (1.d0-epspeak)*nphot
           quantum_mat(itpeak+1,itemp) = quantum_mat(itpeak+1,itemp) + epspeak*nphot
        elseif(itpeak.eq.quantum_ntemp) then
           !
           ! Excited to upper bin
           !
           quantum_mat(itemp,itemp)    = quantum_mat(itemp,itemp) - nphot
           quantum_mat(itpeak,itemp)   = quantum_mat(itpeak,itemp) + nphot
        else
           !
           ! Internal error...
           !
           stop 93971
        endif
     enddo
  enddo
  !
end subroutine quantum_fill_matrix

!-----------------------------------------------------------------
!    COMPUTE PEAK TEMPERATURE FOR EXCITED PAH FROM SOME T-LEVEL 
!
! This subroutine computes the peak temperature obtained by a
! single photon with frequency nu, starting from a given temp level 
!-----------------------------------------------------------------
subroutine compute_peak_temperatures(ntemp,isq)
  implicit none
  integer :: ntemp
  doubleprecision :: cv_pah_gram,energy
  integer :: isq,itemp,ispec,inu,itstart
  doubleprecision :: heatcontent(ntemp)
  !
  ! Check
  !
  if(ntemp.ne.quantum_ntemp) stop 3676
  !
  ! Get the dust ispec from isq
  !
  if(isq.gt.quantum_nrquantum) stop 8310
  ispec = quantum_ispec(isq)
  !
  ! Compute a table for the heat content
  !
  ! Now compute the rest
  ! 
  do itemp=1,quantum_ntemp
     heatcontent(itemp) = dust_mgrain(ispec) * enthalpy_pah_gram(quantum_temp_grid(itemp))
  enddo
  !
  ! First a loop over starting temperatures
  ! 
  do itstart=1,quantum_ntemp-1 
     !
     ! Now do loop over frequency of the quantum-heating mean intensity
     ! array, and compute the peak temperature at every frequency (when a
     ! photon of that frequency hits a PAH).
     !     
     do inu=1,quantum_nf
        !     
        ! Compute the energy of a single photon
        !
        energy = hh * quantum_frequencies(inu)
        !
        ! Find the index of the heat content corresponding
        ! to that energy. This gives the peak temperature
        !
        if(energy.lt.heatcontent(quantum_ntemp-1)-heatcontent(itstart)) then
           !
           ! Find the index of the temperature
           !
           call hunt(heatcontent,quantum_ntemp,energy+heatcontent(itstart),itemp)
           !
           ! Find the actual temperature through linear interpolation
           !
           quantum_peak_itemp(itstart,inu,isq) = itemp
           quantum_peak_eps(itstart,inu,isq)   =                           &
                ( energy + heatcontent(itstart) - heatcontent(itemp) ) /   &
                ( heatcontent(itemp+1) - heatcontent(itemp) )
           if((quantum_peak_eps(itstart,inu,isq).lt.0.d0).or.              &
                (quantum_peak_eps(itstart,inu,isq).gt.1.d0)) then
              stop 82701
           endif
           quantum_peak_temp(itstart,inu,isq)  =                             &
                (1.d0-quantum_peak_eps(itstart,inu,isq)) * quantum_temp_grid(itemp) + &
                      quantum_peak_eps(itstart,inu,isq) * quantum_temp_grid(itemp+1)
           !
           ! If this temperature exceeds the evaporation temperature
           ! then put peak_itemp to quantum_ntemp as a signal
           !
           if((quantum_peak_temp(itstart,inu,isq).gt.dust_tmax(ispec)).and.  &
                (dust_tmax(ispec).ne.0.d0)) then 
              quantum_peak_itemp(itstart,inu,isq) = quantum_ntemp
           endif
        else
           quantum_peak_itemp(itstart,inu,isq) = quantum_ntemp
        endif
     enddo
     !
  end do
end subroutine compute_peak_temperatures


!-----------------------------------------------------------------
!               COMPUTE THE COOLING TIME LINE
!
! This routine returns the time the PAH reaches each of the 
! temperatures of the given temperature grid, starting from the
! highest temperature (i.e. the last temperature in the array),
! and cooling down to the lowest.
! For this, the routine uses dust opacity nr ispec.
!-----------------------------------------------------------------
subroutine quantum_cooling(isq)
  implicit none
  integer :: isq,ispec
  doubleprecision :: temp,dtime,cv_pah_gram
  doubleprecision :: dum,qcool
  doubleprecision :: kappa
  integer :: itemp,inu
  !
  ! Get the dust ispec from isq
  !
  if(isq.gt.quantum_nrquantum) stop 8310
  ispec = quantum_ispec(isq)
  !
  ! Test
  !
  if(quantum_temp_grid(quantum_ntemp).le.quantum_temp_grid(1)) then
     write(stdo,*) 'Quantum Tempeature grid must be in ascending order'
     stop
  endif
  !
  ! Time at max temp = 0
  ! 
  quantum_cooltime(quantum_ntemp,isq) = 0.d0
  !
  ! Now let's cool!
  !
  do itemp=quantum_ntemp-1,1,-1 
     !
     ! Use the temperature of itemp+1, which is important for
     ! the energy conservation of the multi-photon algorithm...
     !
     temp = quantum_temp_grid(itemp+1)
     !
     ! Compute cooling rate per gram of dust
     !
     dum = 0.d0
     do inu=1,freq_nr
        kappa = find_dust_kappa_interpol(                   &
                        freq_nu(inu),ispec,100.d0,1,0,0)
        dum   = dum + bplanck(temp,freq_nu(inu)) * kappa *  &
                        freq_dnu(inu)
     enddo
     qcool = 4 * pi * dum
     !
     ! Compute the time it takes to cool...
     !
     dtime = ( enthalpy_pah_gram(quantum_temp_grid(itemp+1)) -    &
               enthalpy_pah_gram(quantum_temp_grid(itemp)) ) /    &
             qcool
     !
     ! Compute the new time
     !
     quantum_cooltime(itemp,isq) = quantum_cooltime(itemp+1,isq) + dtime 
  enddo
  !
end subroutine quantum_cooling

!-----------------------------------------------------------------
!               ENTHALPY FOR PAH MOLECULES PER GRAM 
!
! According to Siebenmorgen & Kruegel (1992) the enthalpy of 
! PAH molecules in erg/atom is (formula by Chase et al 1985):
!-----------------------------------------------------------------
function enthalpy_pah_gram(temp)
  implicit none
  doubleprecision :: temp,enthalpy_pah_gram
  !
  ! The formula of Chase
  !
  enthalpy_pah_gram = 4.15D-22 * (temp**3.3) /            &
              ( 1.D0 + 6.51D-3*temp + 1.5D-6*temp**2 +    &
                      8.3D-7*(temp**2.3) )
  !
  ! Normalize to gram^(-1)
  !
  enthalpy_pah_gram = enthalpy_pah_gram / (12*mp)      
  !
  return
end function enthalpy_pah_gram

!-----------------------------------------------------------------
!                       COMPUTE EMISSION
!
! Note: This routine *adds* to the emissivity. So you have to
!       reset yourself.
!-----------------------------------------------------------------
subroutine quantum_add_emission(nf,ntemp,freq,kappa,tdist,emissivity)
  implicit none
  integer :: nf,ntemp
  doubleprecision :: freq(nf),kappa(nf),emissivity(nf),tdist(ntemp)
  doubleprecision :: absorb
  integer :: inu,itemp
  !
  ! Add emissivity of each of the temperature points,
  ! except for the itemp=1, because that is the zero-point
  ! 
  do itemp=2,quantum_ntemp
     do inu=1,nf
        emissivity(inu) = emissivity(inu) + tdist(itemp) *  &
             kappa(inu)*bplanck(quantum_temp_grid(itemp),freq(inu))
     enddo
  enddo
end subroutine quantum_add_emission

end module quantum_module
