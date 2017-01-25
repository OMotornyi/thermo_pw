!
! Copyright (C) 2016 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE pheqscf()
  !-----------------------------------------------------------------------
  !
  ! ... This is the main driver for the calculation of the
  ! ... response to an electric field and related quantities.
  !
  USE kinds,           ONLY : DP
  USE io_global,       ONLY : stdout, meta_ionode
  USE paw_variables,   ONLY : okpaw
  USE klist,           ONLY : lgauss
  USE qpoint,          ONLY : xq
  USE uspp,            ONLY : okvan
  USE fft_base,        ONLY : dfftp
  USE wvfct,           ONLY : nbnd, npwx
  USE uspp_param,      ONLY : nhm
  USE ions_base,       ONLY : nat
  USE noncollin_module,ONLY : noncolin, nspin_mag, npol
  USE io_files,        ONLY : tmp_dir
  USE lsda_mod,        ONLY : nspin, lsda
  USE control_ph,      ONLY : convt, zeu, rec_code, rec_code_read, lnoloc, &
                              where_rec, done_epsil, done_zeu, epsil
  USE control_lr,      ONLY : lrpa
  USE control_flags,   ONLY : io_level
  USE output,          ONLY : fildrho
  USE ph_restart,      ONLY : ph_writefile
  USE lrus,            ONLY : int3, int3_nc, int3_paw
  USE freq_ph,         ONLY : fiu, fpol, nfs
  USE optical,         ONLY : current_w, fru, polarc, epsilonc, epsilonm1c,   &
                              intq, intq_nc, dmuxc_tran, chirr, chirz, chizz, &
                              chipm, chimp, chixx, chixy, chizr, lmagnon,     &
                              lall_tensor, lcharge, lchimag, epsm1,           &
                              lr1dwf, iu1dwf
  USE ramanm,          ONLY : ramtns, lraman, elop, done_lraman, done_elop
  USE buffers,         ONLY : close_buffer, open_buffer
  USE images_omega,    ONLY : comp_f
  !
  IMPLICIT NONE
  !
  INTEGER :: iu, ipol, jpol, iu_epsil, ierr
  COMPLEX(DP) :: epsi
  !
  LOGICAL :: exst_mem, exst
  !
  IF (okvan) THEN
     ALLOCATE (int3 ( nhm, nhm, nat, nspin_mag, 1))
     ALLOCATE (intq (nhm, nhm, nat) )
     IF (okpaw) ALLOCATE (int3_paw ( nhm, nhm, nat, nspin_mag, 1))
     IF (noncolin) THEN
        ALLOCATE(int3_nc( nhm, nhm, nat, nspin, 1))
        ALLOCATE(intq_nc( nhm, nhm, nat, nspin))
     ENDIF
     IF (fpol) CALL compute_intq()
  ENDIF

  IF (fpol.AND.lsda) THEN
     ALLOCATE (dmuxc_tran ( dfftp%nnr ) )
     CALL set_fxc_tran()
  END IF
  !
  IF (fpol) THEN    ! calculate freq. dependent polarizability
     !
     WRITE( stdout, '(/,5X,"Frequency Dependent TDDFT Dielectric &
                  &constants and susceptibilities",/)' )
     !
     ALLOCATE(chirr(nfs))
     ALLOCATE(chirz(nfs))
     ALLOCATE(chizr(nfs))
     ALLOCATE(chizz(nfs))
     ALLOCATE(chipm(nfs))
     ALLOCATE(chimp(nfs))
     ALLOCATE(chixx(nfs))
     ALLOCATE(chixy(nfs))
     ALLOCATE(epsm1(nfs))
     chirr=(0.0_DP,0.0_DP)
     chizr=(0.0_DP,0.0_DP)
     chirz=(0.0_DP,0.0_DP)
     chizz=(0.0_DP,0.0_DP)
     chipm=(0.0_DP,0.0_DP)
     chimp=(0.0_DP,0.0_DP)
     chixx=(0.0_DP,0.0_DP)
     chixy=(0.0_DP,0.0_DP)
     epsm1=(0.0_DP,0.0_DP)

     iu1dwf = 42
     lr1dwf = nbnd * npwx * npol
     CALL open_buffer (iu1dwf, 'mwf', lr1dwf, io_level, exst_mem, &
                                                           exst, tmp_dir)
     !
     DO iu = 1, nfs
        !
        IF (.NOT. comp_f(iu) ) CYCLE
        current_w=CMPLX(fru(iu), fiu(iu))
        WRITE( stdout, '(/,5x,70("-"))') 
        WRITE( stdout, '(10x,"q=(",3f15.5,")")') xq(1:3)
                              
        WRITE( stdout, '(/,10x,"Susceptibility at &
              &frequency",f9.4," +",f9.4," i Ry",i6," / ", i6)') current_w, &
                                                         iu, nfs
        !
        IF (.NOT.lsda) THEN
           CALL solve_eq( iu, 1 )
        ELSE
!
!   lchange computes the charge-charge and z-magnetization-charge response
!   the charge-charge response has the plasmon peaks
!
           IF (lcharge) THEN
              WRITE( stdout, '(/10x,"Applying V and computing n and m_z")')
              CALL solve_eq( iu, 1)
           END IF
!
!   lchimag computes the z-magnetization-z-magnetization response
!
           IF (lchimag) THEN
              WRITE( stdout, '(/,5x,70("_"))') 
              WRITE( stdout, '(/10x,"Applying B_z and computing n and m_z")')
              CALL solve_eq( iu, 2)
           END IF
!
!   lmagnon computes the \chi+- susceptibility whose peaks give the
!   magnons and the stoner spin-flip exitations (by default only this
!   is computed)
!
           IF (lmagnon) THEN
              WRITE( stdout, '(/,5x,70("_"))') 
              WRITE( stdout, '(/10x,"Applying B+ and computing m+")')
              CALL solve_eq_tran( iu, 1 )
           ENDIF
!
!   lall_tensor computes the \chi-+ susceptibility and gives all
!   the nonzero (within LSDA) cartesian components of the response.
!   Note that this is equivalent to compute \chi+- at -w* and to take
!   the complex conjugate.
!
           IF (lall_tensor) THEN
              WRITE( stdout, '(/,5x,70("_"))') 
              WRITE( stdout, '(/,10x,"Applying B- and computing m-")')
              CALL solve_eq_tran( iu, 2 )
           ENDIF
        ENDIF
        CALL write_chi_on_disk(iu)
        !
     END DO 

     DEALLOCATE(chirr)
     DEALLOCATE(chirz)
     DEALLOCATE(chizr)
     DEALLOCATE(chizz)
     DEALLOCATE(chipm)
     DEALLOCATE(chimp)
     DEALLOCATE(chixx)
     DEALLOCATE(chixy)
     DEALLOCATE(epsm1)
     !
     WRITE( stdout, '(/,6X,"End of Frequency Dependent Susceptibility Calculation")' )
     !
  ENDIF
  !
  CALL close_buffer(iu1dwf, 'delete')
  IF (okvan) THEN
     DEALLOCATE (intq)
     DEALLOCATE (int3)
     IF (okpaw) DEALLOCATE (int3_paw)
     IF (noncolin) THEN
        DEALLOCATE(int3_nc)
        DEALLOCATE(intq_nc)
     END IF
  ENDIF
  IF (fpol.AND.lsda) DEALLOCATE (dmuxc_tran)
  !
  RETURN
  !
END SUBROUTINE pheqscf

SUBROUTINE write_chi_on_disk(iu)
USE kinds, ONLY : DP
USE optical,          ONLY : current_w, fru, lcharge, &
                             intq, intq_nc, dmuxc_tran, chirr, chirz, chizz, &
                             chipm, chimp, chixx, chixy, chizr, epsm1
USE lsda_mod,         ONLY : nspin, lsda
USE freq_ph,          ONLY : fiu
USE mp_images,        ONLY : my_image_id
USE io_global,        ONLY : stdout, ionode
USE noncollin_module, ONLY : noncolin, nspin_mag

IMPLICIT NONE
INTEGER, INTENT(IN) :: iu
INTEGER :: iu_epsil
COMPLEX(DP) :: epsi
CHARACTER(LEN=6) :: int_to_char
CHARACTER(LEN=256) :: filename

LOGICAL :: exst

IF (ionode) THEN
   IF (nspin_mag==1) THEN
!
!   nonmagnetic case or noncollinear with time reversal
!
       iu_epsil=2
       filename='dynamical_matrices/chirr'
       IF (my_image_id>0) filename=TRIM(filename)//'_'//int_to_char(my_image_id)
       INQUIRE(FILE=TRIM(filename), exist=exst)
       IF (exst) THEN
          OPEN (UNIT=iu_epsil, FILE=TRIM(filename), STATUS='old', &
                               POSITION='append', FORM='formatted')
       ELSE
          OPEN (UNIT=iu_epsil, FILE=TRIM(filename), &
                                 STATUS='unknown', FORM='formatted')
          WRITE(iu_epsil,'("#    Re(w)     Im(w)     Re(chirr)       Im(chirr)&
                             &      Re (1/chirr)      Im(1/chirr) ")')
       END IF

       epsi = (0.0_DP,0.0_DP)
       IF (ABS(chirr(iu))>1.D-10) epsi = CMPLX(1.0_DP,0.0_DP) / chirr(iu)
       WRITE(iu_epsil,'(2f10.5, 4e15.7)') fru(iu), fiu(iu),   &
               DREAL(chirr(iu)), DIMAG(chirr(iu)), &
               DREAL(epsi), DIMAG(epsi)
       CLOSE(iu_epsil)

   ELSEIF (nspin_mag==2) THEN
!
!  lsda case
!
      iu_epsil=2
      filename='dynamical_matrices/chimag_re'
      IF (my_image_id>0) filename=TRIM(filename)//'_'//int_to_char(my_image_id)
      INQUIRE(FILE=TRIM(filename), exist=exst)
      IF (exst) THEN
         OPEN (UNIT=iu_epsil, FILE=TRIM(filename), &
                          STATUS='old', POSITION='append', FORM='formatted')
      ELSE
         OPEN (UNIT=iu_epsil, FILE=TRIM(filename), STATUS='unknown', &
                                             FORM='formatted')
         WRITE(iu_epsil,'("#  Re(w)     Im(w)     rr       rz        zz&
                           &         +-        -+       xx        xy")')
      END IF

      WRITE(iu_epsil,'(2f10.5,7e15.7)') fru(iu), fiu(iu),   &
                  DREAL(chirr(iu)), DREAL(chirz(iu)), &
                  DREAL(chizz(iu)), DREAL(chipm(iu)), &
                  DREAL(chimp(iu)), DREAL(chixx(iu)), &
                  DREAL(chixy(iu))

      CLOSE(iu_epsil)

      iu_epsil=2
      filename='dynamical_matrices/chimag_im'
      IF (my_image_id>0) filename=TRIM(filename)//'_'//int_to_char(my_image_id)
      INQUIRE(FILE=TRIM(filename), exist=exst)
      IF (exst) THEN
         OPEN (UNIT=iu_epsil, FILE=TRIM(filename), &
                        STATUS='old', POSITION='append', FORM='formatted')
      ELSE
         OPEN (UNIT=iu_epsil, FILE=TRIM(filename), &
                            STATUS='unknown', FORM='formatted')
         WRITE(iu_epsil,'("#  Re(w)     Im(w)     rr       rz        zz&
                           &         +-        -+       xx        xy")')
      END IF

      WRITE(iu_epsil,'(2f10.5,7e15.7)') fru(iu), fiu(iu),   &
                    DIMAG(chirr(iu)), DIMAG(chirz(iu)), &
                    DIMAG(chizz(iu)), DIMAG(chipm(iu)), &
                    DIMAG(chimp(iu)), DIMAG(chixx(iu)), &
                    DIMAG(chixy(iu))

      CLOSE(iu_epsil)
   ELSE
!
!  Noncollinear case: not yet implemented
!
   ENDIF

   IF (.NOT.lsda.OR.lcharge) THEN
      iu_epsil=2
      filename='dynamical_matrices/epsilon'
      IF (my_image_id>0) filename=TRIM(filename)//'_'//int_to_char(my_image_id)
      INQUIRE(FILE=TRIM(filename), exist=exst)
      IF (exst) THEN
         OPEN (UNIT=iu_epsil, FILE=TRIM(filename), STATUS='old', &
                               POSITION='append', FORM='formatted')
      ELSE
         OPEN (UNIT=iu_epsil, FILE=TRIM(filename), STATUS='unknown', &
                                                            FORM='formatted')
         WRITE(iu_epsil,'("#    Re(w)     Im(w)     Re(1/eps)       Im(1/eps)&
                        &      Re (eps)      Im(eps) ")')
      END IF

      epsi = (0.0_DP,0.0_DP)
      IF (ABS(epsm1(iu))>1.D-10) epsi = CMPLX(1.0_DP,0.0_DP) / epsm1(iu)
         WRITE(iu_epsil,'(2f10.5, 4e15.7)') fru(iu), fiu(iu),   &
                  DREAL(epsm1(iu)), DIMAG(epsm1(iu)), &
                  DREAL(epsi), DIMAG(epsi)
     
      CLOSE(iu_epsil)
   END IF
END IF

RETURN
END SUBROUTINE write_chi_on_disk

SUBROUTINE collect_all_chi()

USE kinds, ONLY : DP
USE noncollin_module,  ONLY : nspin_mag
USE lsda_mod,  ONLY : lsda
USE freq_ph,   ONLY : fiu, nfs
USE optical,   ONLY : lcharge, fru, chirr, epsm1
USE optical,   ONLY : fru, lcharge, chirr, chirz, chizz, &
                      chipm, chimp, chixx, chixy, chizr, epsm1
USE io_global, ONLY : meta_ionode, ionode, ionode_id
USE mp_images, ONLY : my_image_id, intra_image_comm, inter_image_comm
USE mp,        ONLY : mp_bcast, mp_sum


IMPLICIT NONE
INTEGER :: iu_epsil, iuf, iu
LOGICAL :: exst
REAL(DP) :: chirrr, chirri, epsir, epsii, epsm1r, epsm1i, ro, ri
REAL(DP) :: chirzr, chizzr, chipmr, chimpr, chixxr, chixyr, chirzi, chizzi, &
            chipmi, chimpi, chixxi, chixyi
REAL(DP), ALLOCATABLE :: chirrb(:), chirzb(:), chizzb(:), chipmb(:), &
                         chimpb(:), chixxb(:), chixyb(:), chirrc(:), &
                         chirzc(:), chizzc(:), chipmc(:), chimpc(:), &
                         chixxc(:), chixyc(:)
CHARACTER(LEN=256) :: filename
CHARACTER(LEN=6)   :: int_to_char


ALLOCATE(chirr(nfs))
ALLOCATE(chirz(nfs))
ALLOCATE(chizz(nfs))
ALLOCATE(chipm(nfs))
ALLOCATE(chimp(nfs))
ALLOCATE(chixx(nfs))
ALLOCATE(chixy(nfs))
ALLOCATE(epsm1(nfs))

ALLOCATE(chirrb(nfs))
ALLOCATE(chirzb(nfs))
ALLOCATE(chizzb(nfs))
ALLOCATE(chipmb(nfs))
ALLOCATE(chimpb(nfs))
ALLOCATE(chixxb(nfs))
ALLOCATE(chixyb(nfs))

ALLOCATE(chirrc(nfs))
ALLOCATE(chirzc(nfs))
ALLOCATE(chizzc(nfs))
ALLOCATE(chipmc(nfs))
ALLOCATE(chimpc(nfs))
ALLOCATE(chixxc(nfs))
ALLOCATE(chixyc(nfs))

chirr=0.0_DP
epsm1=0.0_DP

IF (ionode) THEN
   iu_epsil=2
   IF (nspin_mag==1) THEN
!
!  nonmagnetic case of noncollinear with time reversal
!
      filename='dynamical_matrices/chirr'
      IF (my_image_id>0) filename=TRIM(filename)//'_'//int_to_char(my_image_id)
      INQUIRE(FILE=TRIM(filename), exist=exst)
      IF (exst) THEN
         OPEN (UNIT=iu_epsil, FILE=TRIM(filename), STATUS='old', &
                                                              FORM='formatted')
         READ(iu_epsil, *)
      END IF

      DO iuf=1,nfs
         READ(iu_epsil,'(2f10.5, 4e15.7)',END=100) ro, ri,   &
                                                chirrr, chirri, epsir, epsii
         DO iu=1,nfs
            IF ((ABS(ro-fru(iu))+ABS(ri-fiu(iu)))<1.D-4) THEN
               chirr(iu)=CMPLX(chirrr,chirri)
               EXIT
            END IF
         END DO
      END DO
100   CONTINUE
      CLOSE(iu_epsil)

   ELSEIF (nspin_mag==2) THEN

      filename='dynamical_matrices/chimag_re'
      IF (my_image_id>0) filename=TRIM(filename)//'_'//int_to_char(my_image_id)
      INQUIRE(FILE=TRIM(filename), exist=exst)
      IF (exst) THEN
         OPEN (UNIT=iu_epsil, FILE=TRIM(filename), STATUS='old', &
                                                              FORM='formatted')
         READ(iu_epsil, *)
      END IF

      DO iuf=1,nfs
         READ(iu_epsil,'(2f10.5,7e15.7)',END=200) ro, ri, chirrr, chirzr, & 
                                  chizzr, chipmr, chimpr, chixxr, chixyr
         DO iu=1,nfs
            IF ((ABS(ro-fru(iu))+ABS(ri-fiu(iu)))<1.D-4) THEN
               chirrb(iu)=chirrr
               chirzb(iu)=chirzr
               chizz(iu)=chizzr
               chipmb(iu)=chipmr
               chimpb(iu)=chimpr
               chixxb(iu)=chixxr
               chixyb(iu)=chixyr
               EXIT
            ENDIF
         ENDDO
      ENDDO
200   CONTINUE
      CLOSE(iu_epsil)

      filename='dynamical_matrices/chimag_im'
      IF (my_image_id>0) filename=TRIM(filename)//'_'//int_to_char(my_image_id)
      INQUIRE(FILE=TRIM(filename), exist=exst)
      IF (exst) THEN
         OPEN (UNIT=iu_epsil, FILE=TRIM(filename), STATUS='old', &
                                                              FORM='formatted')
         READ(iu_epsil, *)
      END IF

      DO iuf=1,nfs
         READ(iu_epsil,'(2f10.5,7e15.7)',END=300) ro, ri, chirri, chirzi, & 
                                  chizzi, chipmi, chimpi, chixxi, chixyi
         DO iu=1,nfs
            IF ((ABS(ro-fru(iu))+ABS(ri-fiu(iu)))<1.D-4) THEN
               chirrc(iu)=chirri
               chirzc(iu)=chirzi
               chizzc(iu)=chizzi
               chipmc(iu)=chipmi
               chimpc(iu)=chimpi
               chixxc(iu)=chixxi
               chixyc(iu)=chixyi
               EXIT
            ENDIF
         ENDDO
      ENDDO
300   CONTINUE
      CLOSE(iu_epsil)
!
!  Finally collect all the results
!
      chirr(:)=CMPLX(chirrb(:),chirrc(:))
      chirz(:)=CMPLX(chirzb(:),chirzc(:))
      chizz(:)=CMPLX(chizzb(:),chizzc(:))
      chipm(:)=CMPLX(chipmb(:),chipmc(:))
      chimp(:)=CMPLX(chimpb(:),chimpc(:))
      chixx(:)=CMPLX(chixxb(:),chixxc(:))
      chixy(:)=CMPLX(chixyb(:),chixyc(:))

   ELSE
!
!   Noncollinear magnetic case, not yet implemented
!
   ENDIF

   IF (.NOT.lsda.OR.lcharge) THEN
      iu_epsil=2
      filename='dynamical_matrices/epsilon'
      IF (my_image_id>0) filename=TRIM(filename)//'_'//int_to_char(my_image_id)
      INQUIRE(FILE=TRIM(filename), exist=exst)
      IF (exst) THEN
         OPEN (UNIT=iu_epsil, FILE=TRIM(filename), STATUS='old', &
                                                            FORM='formatted')
         READ(iu_epsil, *)
      END IF

      DO iuf=1,nfs
         READ(iu_epsil,'(2f10.5, 4e15.7)',END=400) ro, ri, epsm1r, epsm1i, &
                                                                   epsir, epsii
         DO iu=1,nfs
            IF ((ABS(ro-fru(iu))+ABS(ri-fiu(iu)))<1.D-4) THEN
               epsm1(iu)=CMPLX(epsm1r,epsm1i)
               EXIT
            END IF
         END DO
      END DO
400   CONTINUE
      CLOSE(iu_epsil)
   ENDIF
ENDIF
!
! Now transfer the susceptibility and epsm1 to all nodes. First collect
! all the frequencies and then send to all the nodes of each image
!
IF (nspin_mag==1) THEN
   CALL mp_sum(chirr, inter_image_comm)
   CALL mp_bcast(chirr, ionode_id, intra_image_comm)
ELSEIF(nspin_mag==2) THEN
   CALL mp_sum(chirr, inter_image_comm)
   CALL mp_bcast(chirr, ionode_id, intra_image_comm)
   CALL mp_sum(chirz, inter_image_comm)
   CALL mp_bcast(chirz, ionode_id, intra_image_comm)
   CALL mp_sum(chizz, inter_image_comm)
   CALL mp_bcast(chizz, ionode_id, intra_image_comm)
   CALL mp_sum(chipm, inter_image_comm)
   CALL mp_bcast(chipm, ionode_id, intra_image_comm)
   CALL mp_sum(chimp, inter_image_comm)
   CALL mp_bcast(chimp, ionode_id, intra_image_comm)
   CALL mp_sum(chixx, inter_image_comm)
   CALL mp_bcast(chixx, ionode_id, intra_image_comm)
   CALL mp_sum(chixy, inter_image_comm)
   CALL mp_bcast(chixy, ionode_id, intra_image_comm)
ELSE

ENDIF

IF (.NOT.lsda.OR.lcharge) THEN
   CALL mp_sum(epsm1, inter_image_comm)
   CALL mp_bcast(epsm1, ionode_id, intra_image_comm)
ENDIF
!
!  Now remove all the files created by the images
!
IF (ionode) THEN
   IF (nspin_mag==1) THEN
      filename='dynamical_matrices/chirr'
      IF (my_image_id>0) filename=TRIM(filename)//'_'//int_to_char(my_image_id)
      OPEN (UNIT=iu_epsil, FILE=TRIM(filename), STATUS='old', FORM='formatted')
      CLOSE (UNIT=iu_epsil, STATUS='DELETE')
   ELSEIF (nspin_mag==2) THEN
      filename='dynamical_matrices/chimag_re'
      IF (my_image_id>0) filename=TRIM(filename)//'_'//int_to_char(my_image_id)
      OPEN (UNIT=iu_epsil, FILE=TRIM(filename), STATUS='old', FORM='formatted')
      CLOSE (UNIT=iu_epsil, STATUS='DELETE')
      filename='dynamical_matrices/chimag_im'
      IF (my_image_id>0) filename=TRIM(filename)//'_'//int_to_char(my_image_id)
      OPEN (UNIT=iu_epsil, FILE=TRIM(filename), STATUS='old', FORM='formatted')
      CLOSE (UNIT=iu_epsil, STATUS='DELETE')
   ELSE
   ENDIF

   IF (.NOT.lsda.OR.lcharge) THEN
      filename='dynamical_matrices/epsilon'
      IF (my_image_id>0) filename=TRIM(filename)//'_'//int_to_char(my_image_id)
      OPEN (UNIT=iu_epsil, FILE=TRIM(filename), STATUS='old', FORM='formatted')
      CLOSE (UNIT=iu_epsil, STATUS='DELETE')
   ENDIF
ENDIF
!
!  now rewrite on a single file all the collected chi.
!  Only one processor writes on disk
!
IF (meta_ionode) THEN
   DO iu=1,nfs
      CALL write_chi_on_disk(iu)
   END DO
END IF

DEALLOCATE(chirr)
DEALLOCATE(chirz)
DEALLOCATE(chizz)
DEALLOCATE(chipm)
DEALLOCATE(chimp)
DEALLOCATE(chixx)
DEALLOCATE(chixy)
DEALLOCATE(epsm1)

DEALLOCATE(chirrb)
DEALLOCATE(chirzb)
DEALLOCATE(chizzb)
DEALLOCATE(chipmb)
DEALLOCATE(chimpb)
DEALLOCATE(chixxb)
DEALLOCATE(chixyb)

DEALLOCATE(chirrc)
DEALLOCATE(chirzc)
DEALLOCATE(chizzc)
DEALLOCATE(chipmc)
DEALLOCATE(chimpc)
DEALLOCATE(chixxc)
DEALLOCATE(chixyc)

RETURN
END SUBROUTINE collect_all_chi
