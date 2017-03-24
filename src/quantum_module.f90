module quantum_module
use rtglobal_module
use constants_module

doubleprecision,allocatable :: emisquant(:,:)
doubleprecision,allocatable :: emisquant_loccum(:,:)
doubleprecision,allocatable :: emisquant_loctot(:)
doubleprecision,allocatable :: emisquant_cum(:)
doubleprecision,allocatable :: miquant(:,:)
doubleprecision :: emisquanttot
doubleprecision,allocatable :: quantum_frequencies(:)
doubleprecision,allocatable :: quantum_dnu(:)
integer :: quantum_nf

contains

!---------------------------------------------------------------------
!                 Initialize the quantum mode
!---------------------------------------------------------------------
subroutine quantum_init()
implicit none
integer :: ierr
write(stdo,*) 'Initializing quantum module'
if(nrcellsmax.le.0) then
   write(stdo,*) 'INTERNAL ERROR IN QUANTUM_MODULE: MONTE CARLO NOT INITIALIZED'
   stop
endif
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
allocate(miquant(1:freq_nr,1:nrcellsmax),STAT=ierr)
if(ierr.ne.0) then
   write(stdo,*) 'ERROR: Could not allocate miquant.'
   stop
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
if(allocated(miquant)) deallocate(miquant)
if(allocated(quantum_frequencies)) deallocate(quantum_frequencies)
if(allocated(quantum_dnu)) deallocate(quantum_dnu)
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
