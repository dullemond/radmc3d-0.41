module quantum_module
use rtglobal_module
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
integer,allocatable :: quantum_ispec(:)
integer :: quantum_nrquantum=0
integer :: quantum_ntemp=100
doubleprecision :: quantum_temp0=1d0
doubleprecision :: quantum_temp1=1d4

contains

!---------------------------------------------------------------------
!                 Initialize the quantum mode
!---------------------------------------------------------------------
subroutine quantum_init(set_tempgrid,incl_emisquant)
implicit none
integer :: ierr,itemp,ispec,isq
logical :: incl_emisquant,set_tempgrid
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
! Allocate the temperature grid. Note that the parameters can be
! fine-tuned in radmc3d.inp
!
if(set_tempgrid) then
   !
   ! The temperature grid itself
   !
   allocate(quantum_temp_grid(1:quantum_ntemp))
   do itemp=1,quantum_ntemp
      quantum_temp(itemp) = quantum_temp0 *     &
           (quantum_temp1/quantum_temp0)**      &
           ((itemp-1.d0)/(quantum_ntemp-1.d0))
   enddo
   !
   ! Find which of the dust species are quantum
   !
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
   !
   ! Allocate the array for the temperature distribution
   ! function, but only for the quantum grains...
   !
   allocate(quantum_temp_distr(quantum_ntemp,quantum_nrquantum,nrcellsmax))
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
   allocate(emisquant(1:freq_nr,1:nrcellsmax),STAT=ierr)
   if(ierr.ne.0) then
      write(stdo,*) 'ERROR: Could not allocate emisquant.'
      stop
   endif
   allocate(emisquant_loccum(1:freq_nr+1,1:nrcellsmax),STAT=ierr)
   if(ierr.ne.0) then
      write(stdo,*) 'ERROR: Could not allocate emisquant_loccum.'
      stop
   endif
   allocate(emisquant_loctot(1:nrcellsmax),STAT=ierr)
   if(ierr.ne.0) then
      write(stdo,*) 'ERROR: Could not allocate emisquant_loctot.'
      stop
   endif
   allocate(emisquant_cum(nrcellsmax+1),STAT=ierr)
   if(ierr.ne.0) then
      write(stdo,*) 'ERROR: Could not allocate emisquant_cum.'
      stop
   endif
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
if(allocated(quantum_ispec)) deallocate(quantum_ispec)
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
  allocate(quantum_frequencies(1:quantum_nf),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in Montecarlo Module: Could not allocate quantum_frequencies'
     stop 
  endif
  allocate(quantum_dnu(1:quantum_nf),STAT=ierr)
  if(ierr.ne.0) then
     write(stdo,*) 'ERROR in Montecarlo Module: Could not allocate quantum_dlgnu'
     stop 
  endif
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

end module quantum_module
