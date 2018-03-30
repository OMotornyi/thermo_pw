!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------
SUBROUTINE lanczos_write_restart(LR_iteration)
  !---------------------------------------------------------------------
  ! 
  ! This subroutine reads in and stores vectors necessary to
  ! restart the Lanczos recursion.
  !
  ! Modified by Osman Baris Malcioglu (2009)
  ! Modified by Xiaochuan Ge (May. 2013) to adapt pseudo-hermitian
  !
 USE kinds,                ONLY : DP
  USE io_files,             ONLY : tmp_dir, prefix, diropn
  
  USE lr_lanczos,   ONLY : evc1, evc1_new, evc1_old, sevc1, beta_store, &
                           gamma_store, zeta_store,iulanczos,iunrestart,nwordrestart
  USE lr_global,    ONLY :  rpert,size_evc1,pseudo_hermitian

  USE wvfct,                ONLY : nbnd, npwx
  USE fft_base,             ONLY : dfftp
  USE io_global,            ONLY : ionode, stdout
  USE klist,                ONLY : nks, nelec
  USE noncollin_module,     ONLY : nspin_mag, noncolin, npol
  use lsda_mod,             ONLY : nspin
  USE cell_base,            ONLY : alat, omega
  USE qpoint,               ONLY : nksq,xq
  !
  IMPLICIT NONE
  CHARACTER(len=6), EXTERNAL :: int_to_char
  !
  ! local variables
REAL(kind=dp) :: norm0(3)
  !
  INTEGER, EXTERNAL :: find_free_unit
  INTEGER, INTENT(IN) :: LR_iteration
  INTEGER :: i, j, k, pol_index,ibnd_occ,ibnd_virt
  CHARACTER (len=24) :: bgz_suffix
  CHARACTER(len=256) :: tempfile, filename
  LOGICAL :: exst
  real(kind=dp) :: degspin
  !
  !IF (lr_verbosity > 5) THEN
   WRITE(stdout,'("<lr_write_restart>")')
!  ENDIF
  !
  ! Note: ionode only operations are carried out in tmp_dir not wfc_dir
  !
  ! If there is only one polarization dir, storage is one rank less.
  !
  pol_index = 1
  ! 
! IF ( n_ipol /= 1 ) pol_index = LR_polarization
! !
! IF (eels) tmp_dir = tmp_dir_lr
! !
#if defined(__MPI)
  IF (ionode) THEN
#endif
  !
  !PRINT *, "lanczos_write_restart"
  ! Writing beta, gamma and zeta coefficients.
   bgz_suffix = TRIM ( ".beta_gamma_z." )
  !
     filename = trim(prefix) // trim(bgz_suffix) // trim("dat")
  tempfile = trim(tmp_dir) // trim(filename)
  !
  OPEN (158, file = tempfile, form = 'formatted', status = 'unknown')
  WRITE(158,*) LR_iteration
  !
  norm0(pol_index) = beta_store(1)
  WRITE(158,*) norm0(pol_index)
  !
  IF (nspin==2) THEN
        degspin = 1.0d0
  ELSE
        degspin = 2.0d0
  ENDIF
  IF (noncolin) degspin = 1.0d0
  !
  ! Write the degenaracy wrt spin
  !
  WRITE(158,*) degspin
  !
  ! ------ Needed for EELS ----------
  !
  ! Write the lattice parameter
  !
  WRITE(158,*) alat
  !
  ! Write the unit-cell volume
  !
  WRITE(158,*) omega
  !
  ! Write the number of valence (and semicore electrons) in the unit
  ! cell
  !
  WRITE(158,*) nelec
  !
  ! Write the components of the transferred momentum
  !
  WRITE(158,*) xq(1)
  WRITE(158,*) xq(2)
  WRITE(158,*) xq(3)
  !
  !-----------------------------------
  !
  DO i=1,LR_iteration-1
     !
     WRITE(158,*) beta_store(i+1)
     WRITE(158,*) gamma_store(i+1)
     !
     ! This is absolutely necessary for cross platform compatibility
     !
     DO j=1,rpert
      DO k=1,rpert
       WRITE(158,*) zeta_store (j,k,i)
      ENDDO
     ENDDO
     !
  ENDDO
  !
  ! X. Ge: Compatable with the old version. The beta & gamma will not be
  ! used in 
  ! the spectrum calculation.
  !
  WRITE(158,*) beta_store(LR_iteration)             
  WRITE(158,*) gamma_store(LR_iteration)             
  DO j=1,rpert                        
      DO k=1,rpert
       WRITE(158,*) zeta_store (j,k,i)
      ENDDO
  ENDDO
  !
  CLOSE(158)
! !
#if defined(__MPI)
  ENDIF
#endif
    !
    ! Parallel writing operations
    !
    ! Note: Restart files are writen in outdir.
    ! If you do not want them to be written,
    ! just disable restart saving completely.
    !
    ! Writing wavefuncion files for restart
    !
    nwordrestart = 2 * nbnd * npwx * npol * nksq
    nwordrestart = 2 * nbnd * npwx * npol * nksq*rpert
    IF (.NOT. pseudo_hermitian) nwordrestart=nwordrestart*2
!    nwordrestart = size_evc1
    !
              iunrestart = find_free_unit()
    CALL diropn ( iunrestart,'restart_lanczos.'//trim(int_to_char(1)),nwordrestart, exst)
    !
    !
!  IF (pseudo_hermitian) THEN
     CALL davcio(evc1(:,:,:,:),nwordrestart,iunrestart,1,1)
     CALL davcio(evc1_old(:,:,:,:),nwordrestart,iunrestart,2,1)
!  ELSE
!    CALL davcio(evc1(:,:,:,1),nwordrestart,iunrestart,1,1)
!    CALL davcio(evc1(:,:,:,2),nwordrestart,iunrestart,2,1)
!    CALL davcio(evc1_old(:,:,:,1),nwordrestart,iunrestart,3,1)
!    CALL davcio(evc1_old(:,:,:,2),nwordrestart,iunrestart,4,1)
!  ENDIF
    CLOSE( unit = iunrestart)
!   !
    RETURN
    !
!-----------------------------------------------------------------------

END SUBROUTINE lanczos_write_restart
