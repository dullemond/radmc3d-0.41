!======================================================================
!             PAH / VSG MODULE: QUANTUM-HEATED GRAINS
!
! This is a port of a code that was originally made for use with the
! original RADMC code. That's why there are comments from before the
! RADMC-3D code in there. 
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
integer :: ierr,itemp,ispec,isq
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
if(allocated(quantum_mat)) deallocate(quantum_mat)
if(allocated(quantum_cooltime)) deallocate(quantum_cooltime)
if(allocated(quantum_ispec)) deallocate(quantum_ispec)
if(allocated(quantum_kappa)) deallocate(quantum_kappa)
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
  call quantum_fill_matrix(ispec,meanint)
  !
  ! This matrix obviously has one zero eigenvalue because all the
  ! temperature bins are converted into each other. This cannot be
  ! inverted.
  !
  ! So replace last line with normalization
  ! But remember the original line for later in this routine.
  !
  do itemp=1,quantum_ntemp
     matbk(itemp)             = mat(quantum_ntemp,itemp)
     mat(quantum_ntemp,itemp) = 1.d0
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
        if(abs(mat(itemp,it)).gt.dummy) dummy=abs(mat(itemp,it))
     enddo
     do it=1,quantum_ntemp
        mat(itemp,it) = mat(itemp,it) / dummy
     enddo
  enddo
  !
  ! Solve this matrix equation
  !
  dummy = 1.d0
  do itemp=1,quantum_ntemp
     do it=1,quantum_ntemp
        matwork(itemp,it) = mat(itemp,it)
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
        tdtmp(itemp) = tdtmp(itemp) + mat(itemp,it) * tdist(it)
     enddo
     if(itemp.lt.quantum_ntemp) then
        dummy = dummy + abs(tdtmp(itemp))
     else
        dummy = dummy + abs(tdtmp(itemp)-1.d0)
     endif
  enddo
  if(dummy.gt.1d-10) then
     write(*,*) 'ERROR in matix inversion quantum: Error = ',dummy
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
     write(*,*) 'ERROR: Distribution does not add up to 1,'
     write(*,*) '   at ir,it = ',irr,itt
     write(*,*) (tdist(itemp),itemp=1,quantum_ntemp)
     write(*,*) dummy
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
     write(*,*) 'ERROR: Distribution does not add up to 1,'
     write(*,*) (tdist(itemp),itemp=1,quantum_ntemp)
     write(*,*) dummy
     open(unit=2,file='matrix.dat')
     write(2,*) quantum_ntemp
     do itemp=1,quantum_ntemp
        do it=1,quantum_ntemp
           write(2,*) mat(itemp,it)
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
  if(pahdes_temp(ispec).gt.0.d0) then
     if(pahdes_pcrit(ispec).eq.0.d0) then
        !
        ! The original method of Kees
        !
        ! Find the index of the temperature grid 
        !
        call hunt(quantum_temp,quantum_ntemp,pahdes_temp(ispec),ievap)
        if(ievap.ge.quantum_ntemp) stop 90934
        if(ievap.lt.1) then
           write(*,*) pahdes_temp(ispec)
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
           mat(quantum_ntemp,itemp) = matbk(itemp)
        enddo
        !             
        ! Do matrix multip
        !
        do itemp=1,quantum_ntemp
           res(itemp) = 0.d0
           do it=1,quantum_ntemp
              res(itemp) = res(itemp) + mat(itemp,it) * tdtmp(it)
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
        pahdes_time(ispec) = 1.d0 / ( abs(dummy) + 1d-99 )
     else
        !
        ! The Habart method, slightly generalized
        ! [ADDED: 25.11.06]
        !
        ! Find the index of the temperature grid 
        !
        call hunt(quantum_temp,quantum_ntemp,pahdes_temp(ispec),ievap)
        if(ievap.ge.quantum_ntemp) stop 90934
        if(ievap.lt.1) then
           write(*,*) pahdes_temp(ispec)
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
        pahdes_time(ispec) = 3.1536d+13 * ( pahdes_pcrit(ispec) / dummy )
        !
     endif
  else
     pahdes_time(ispec) = 1d99
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
        itpeak  = peak_itemp(itemp,inu)
        epspeak = peak_eps(itemp,inu)
        !
        ! Find out how many photons are absorbed per second
        ! by 1 gram of PAHs 
        !
        nphot = 4 * pi * meanint(inu) * dnu(inu) *   &
                quatum_kappa(inu) /                  &
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


end module quantum_module
