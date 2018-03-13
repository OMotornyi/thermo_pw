!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------
SUBROUTINE lr_restart(iter_restart,rflag)
  !---------------------------------------------------------------------
  !
  ! Restart the Lanczos recursion
  !
  ! Modified by Osman Baris Malcioglu (2009)
  ! Modified by Iurii Timrov (2013)
  !
  USE kinds,                ONLY : DP 
  USE io_global,            ONLY : stdout, ionode_id
  USE control_flags,        ONLY : gamma_only
  USE klist,                ONLY : nks, xk, ngk, igk_k
  USE io_files,             ONLY : tmp_dir, prefix, diropn, wfc_dir
  USE lr_lanczos,   ONLY : evc1, evc1_new, evc1_old, sevc1, beta_store, &
                             gamma_store, zeta_store,iulanczos,iunrestart,nwordrestart,lanczos_steps
  USE lr_global,    ONLY :  rpert, size_evc1
  USE wvfct,                ONLY : nbnd, npwx
  USE io_global,            ONLY : ionode
  USE mp,                   ONLY : mp_bcast
  USE mp_world,             ONLY : world_comm
  USE noncollin_module,     ONLY : nspin_mag, npol
  USE qpoint,               ONLY : nksq,xq

  IMPLICIT NONE
  !
  INTEGER, INTENT(OUT) :: iter_restart
  LOGICAL, INTENT(OUT) :: rflag
  !
  INTEGER, EXTERNAL :: find_free_unit
  CHARACTER(len=6), EXTERNAL :: int_to_char
  !
  CHARACTER (len=24) :: bgz_suffix
  ! local variables
  !
  INTEGER :: i,ibnd,ibnd_occ,ibnd_virt,temp
  INTEGER :: ik, ig, ip
  LOGICAL :: exst
  CHARACTER(len=256) :: tempfile, filename, tmp_dir_saved
  INTEGER :: pol_index
REAL(kind=dp) :: norm0(3)
  !
!  IF (lr_verbosity > 5) THEN
!    WRITE(stdout,'("<lr_restart>")')
!  ENDIF
  !
  pol_index = 1
  !
  !
 ! IF (.not.restart) RETURN
  !
  rflag = .false.
  !
  ! Optical case: recompute the kintic-energy g2kin and 
  ! beta functions vkb (needed only in the US case).
  ! Note, this is done only in the gamma_only case,
  ! because in the k-points version all is recomputed
  ! on-the-fly for every k point.
  !
  !
  ! Reading Lanczos coefficients
   bgz_suffix = TRIM ( ".beta_gamma_z." )
  !
!  IF (eels) THEN
    filename = trim(prefix) // trim(bgz_suffix) // trim("dat")
!  ELSE
!    filename = trim(prefix) // trim(bgz_suffix) // trim(int_to_char(LR_polarization))
!  ENDIF
  tempfile = trim(tmp_dir) // trim(filename)
  !
  INQUIRE (file = tempfile, exist = exst)
  !
  IF (.not.exst) THEN
     !
     WRITE( stdout,*) "WARNING: " // trim(filename) // " does not exist"
     rflag = .true.
     RETURN
     !
  ENDIF
  !
  ! Ionode only reads
  ! Note: ionode file I/O is done in tmp_dir
  !
#if defined(__MPI)
  IF (ionode) THEN
#endif
  !
  ! Read and broadcast beta, gamma, and zeta.
  !
  OPEN (158, file = tempfile, form = 'formatted', status = 'old')
  !
  READ(158,*,end=301,err=303) iter_restart
  !
  IF ( iter_restart >= lanczos_steps ) iter_restart = lanczos_steps
  !
  READ(158,*,end=301,err=303) norm0(pol_index)
  !
  ! X. Ge: New structure for the pseudo-Hermitian algorithm.
  !
  beta_store(1) = norm0(pol_index)
  !
  DO i=1,7
     READ(158,*,end=301,err=303)
  ENDDO
  ! 
  DO i=1,iter_restart-1
     READ(158,*,end=301,err=303) beta_store(i+1)
     READ(158,*,end=301,err=303) gamma_store(i+1)
     READ(158,*,end=301,err=303) zeta_store (:,:,i)
  ENDDO
  !
  READ(158,*,end=301,err=303) beta_store(iter_restart)
  READ(158,*,end=301,err=303) gamma_store(iter_restart)
  READ(158,*,end=301,err=303) zeta_store (:,:,iter_restart)
  !
  CLOSE(158)
  !
#if defined(__MPI)
  ENDIF
  CALL mp_bcast (iter_restart, ionode_id, world_comm)
  CALL mp_bcast (norm0(pol_index), ionode_id, world_comm)
  CALL mp_bcast (beta_store(:), ionode_id, world_comm)
  CALL mp_bcast (gamma_store(:), ionode_id, world_comm)
  CALL mp_bcast (zeta_store(:,:,:), ionode_id, world_comm)
#endif
  !
  ! Optical case: read projection
  !
  ! 
  iter_restart = iter_restart + 1 ! OM: NOT SURE IF IT IS NEEDED
  !
  ! Parallel reading
  ! Note: Restart files are always in outdir
  ! Reading Lanczos vectors
  !
  nwordrestart = 2 * nbnd * npwx * npol * nksq
!  nwordrestart = size_evc1
  !
  iunrestart = find_free_unit()
  CALL diropn ( iunrestart, 'restart_lanczos.'//trim(int_to_char(1)), nwordrestart, exst)
  !
 CALL davcio(evc1(:,:,:,1),nwordrestart,iunrestart,1,-1)
 CALL davcio(evc1(:,:,:,2),nwordrestart,iunrestart,2,-1)
 CALL davcio(evc1_old(:,:,:,1),nwordrestart,iunrestart,3,-1)
 CALL davcio(evc1_old(:,:,:,2),nwordrestart,iunrestart,4,-1)
  !
  CLOSE( unit = iunrestart)
  !
  ! Optical case: read the response charge density
  !
  !
  ! End of all file I/O for restart.
  PRINT *, "RESTART READ CORRECTLY?"
  !
  RETURN
  !
  301 CALL errore ('restart', 'A File is corrupted, file ended unexpectedly', 1 )
  303 CALL errore ('restart', 'A File is corrupted, error in reading data', 1)
  !
END SUBROUTINE lr_restart
