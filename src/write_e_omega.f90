!
! Copyright (C) 2014 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE write_e_omega()
!
! This routine receives as input the coefficients of the quadratic
! or quartic polynomial that fit the energy as a function of the volume
! and gives as output the energy as a function of the volume and
! the pressure as a function of the volume. In the anisotropic case
! it gives also the crystal parameters as a function of the pressure.
!
! in output omega in (a.u.)**3, p in kbar, e in Ry
! 
!
USE kinds,            ONLY : DP
USE constants,        ONLY : ry_kbar
USE data_files,       ONLY : flevdat
USE cell_base,        ONLY : ibrav
USE thermo_mod,       ONLY : omega_geo, celldm_geo, energy_geo
USE control_pressure, ONLY : pressure_kb
USE control_mur,      ONLY : vmin_input, vmax_input
USE control_quartic_energy, ONLY : lquartic, lsolve
USE quadratic_surfaces, ONLY : fit_multi_quadratic, find_fit_extremum
USE quartic_surfaces, ONLY : fit_multi_quartic, compute_quartic_var, &
                             find_quartic_extremum
USE mp_images,        ONLY : root_image, my_image_id
USE io_global,        ONLY : ionode

IMPLICIT NONE

CHARACTER(LEN=256) :: filename, filename1
CHARACTER(LEN=8)   :: float_to_char
INTEGER  :: i, iu_mur, npress, ipress, idata, degree, nvar, nvar4, ndata
INTEGER  :: find_free_unit, compute_nwork
REAL(DP) :: press_min, press_max, deltap, ymin, ymin4
REAL(DP) :: compute_omega_geo
REAL(DP), ALLOCATABLE :: coeff(:), f(:), x(:,:), x_pos_min(:), x_min_4(:), &
                         coeff4(:), p(:), e(:), omega(:), celldmp(:,:)

IF (my_image_id /= root_image) RETURN

filename="energy_files/"//TRIM(flevdat)//'_mur'
filename1="energy_files/"//TRIM(flevdat)//'_mur_celldm'
IF (pressure_kb /= 0.0_DP) &
   filename=TRIM(filename)//'.'//TRIM(float_to_char(pressure_kb,1))
IF (pressure_kb /= 0.0_DP) &
   filename1=TRIM(filename1)//'.'//TRIM(float_to_char(pressure_kb,1))

npress=51
ndata=compute_nwork()
press_min=-50.0_DP / ry_kbar
press_max=80.0_DP / ry_kbar
deltap= ( press_max - press_min ) / npress

CALL compute_degree(ibrav,degree,nvar)

ALLOCATE(x(degree,ndata))
ALLOCATE(x_pos_min(degree))
ALLOCATE(f(ndata))
ALLOCATE(coeff(nvar))
ALLOCATE(p(npress))
ALLOCATE(e(npress))
ALLOCATE(omega(npress))
ALLOCATE(celldmp(6,npress))

IF (lquartic) THEN
   nvar4=compute_quartic_var(degree)
   ALLOCATE(x_min_4(degree))
   ALLOCATE(coeff4(nvar4))
ENDIF

DO ipress=1, npress
   p(ipress) = press_min + ( ipress - 1 ) * deltap
   DO idata=1, ndata
     f(idata)=energy_geo(idata) + p(ipress) * omega_geo(idata)
   END DO
   CALL set_x_from_celldm(ibrav, degree, ndata, x, celldm_geo)
   CALL fit_multi_quadratic(ndata,degree,nvar,x,f,coeff)
   CALL find_fit_extremum(degree,nvar,x_pos_min,ymin,coeff)
   IF (lquartic) THEN
      CALL fit_multi_quartic(ndata,degree,nvar4,lsolve,x,f,coeff4)
      x_min_4=x_pos_min
      CALL find_quartic_extremum(degree,nvar4,x_min_4,ymin4,coeff4)
      CALL set_celldm_from_xmin(ibrav, degree, x_min_4, celldmp(1,ipress))
      omega(ipress)=compute_omega_geo(ibrav,celldmp(1,ipress))
      e(ipress)=ymin4 - p(ipress) * omega(ipress)
   ELSE
      CALL set_celldm_from_xmin(ibrav, degree, x_pos_min, celldmp(1,ipress))
      omega(ipress)=compute_omega_geo(ibrav,celldmp(1,ipress))
      e(ipress) = ymin - p(ipress) * omega(ipress)
   ENDIF
ENDDO

IF (vmin_input == 0.0_DP) vmin_input=omega(npress) * 0.98_DP
IF (vmax_input == 0.0_DP) vmax_input=omega(1) * 1.02_DP

IF (ionode) THEN
!
!  Print the volume and the energy as a function of volume
!
   iu_mur=find_free_unit()
   OPEN(UNIT=iu_mur, FILE=TRIM(filename), STATUS='UNKNOWN', FORM='FORMATTED')
   IF (pressure_kb /= 0.0_DP) THEN
      WRITE(iu_mur,'( "# omega (a.u.)**3      enthalpy (Ry)   pressure (kbar)" )')
   ELSE
      WRITE(iu_mur,'( "# omega (a.u.)**3       energy (Ry)      pressure (kbar)" )')
   END IF
   DO ipress=1,npress
      WRITE(iu_mur,'(3f20.10)') omega(ipress), e(ipress), p(ipress) * ry_kbar
   ENDDO 
   CLOSE(UNIT=iu_mur, STATUS='KEEP')
!
!  print the crystal parameter as a function of pressure
!
   OPEN(UNIT=iu_mur, FILE=TRIM(filename1), STATUS='UNKNOWN', FORM='FORMATTED')
   WRITE(iu_mur,'( "# pressure (kbar)   celldm(1)   celldm(2)  celldm(3) &
                      celldm(4)   celldm(5)    celldm(6) " )')
   DO ipress=1,npress
      WRITE(iu_mur,'(7f15.8)') p(ipress)*ry_kbar, celldmp(1, ipress), &
               celldmp(2, ipress), celldmp(3, ipress), celldmp(4, ipress), &
               celldmp(5, ipress), celldmp(6, ipress)
   ENDDO 
   CLOSE(UNIT=iu_mur, STATUS='KEEP')
END IF

DEALLOCATE(x)
DEALLOCATE(x_pos_min)
DEALLOCATE(f)
DEALLOCATE(coeff)
DEALLOCATE(p)
DEALLOCATE(e)
DEALLOCATE(omega)
DEALLOCATE(celldmp)

IF (lquartic) THEN
   DEALLOCATE(x_min_4)
   DEALLOCATE(coeff4)
ENDIF

RETURN
END SUBROUTINE write_e_omega
