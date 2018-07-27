subroutine intersect_elsewhere( &
     poly, &
     xyz, &
     separable, &
     randomize )
  use mod_math
  use mod_polynomial
  use mod_geometry
  use mod_separation
  implicit none
  LOGICAL, PARAMETER :: DEBUG = ( GLOBALDEBUG .AND. .false. )
  type(ptr_polynomial), intent(in)  :: poly(2)
  real(kind=fp),        intent(in)  :: xyz(3)
  logical,              intent(out) :: separable
  logical, optional,    intent(in)  :: randomize
  type(type_matrix)                 :: sep(2)
  integer                           :: nbcp, nsep(2)
  real(kind=fp)                     :: rot(3,3), vec(3)
  integer                           :: i

  ! rearrange Bezier control points
  do i = 1,2 ! <-----------------------------------------------------------+
     !                                                                     !                                                         
     if ( poly(i)%ptr%dim < 3 ) STOP 'intersect_elsewhere : poly%dim < 3'  !   
     !                                                                     !                                                         
     select case (poly(i)%ptr%nvar) ! <--------------------------------+   !
        !                                                              !   !
     case (1) ! -------------------------------------------------------+   !
        nbcp = poly(i)%ptr%degr(1) + 1                                 !   !
        allocate( sep(i)%mat(nbcp,3) )                                 !   !
        call rearrange_for_separability_test( &                        !   !
             poly(i)%ptr%coef(1:nbcp,1:3,1), &                         !   !
             nbcp, &                                                   !   !
             xyz, &                                                    !   !
             sep(i)%mat, &                                             !   !
             nsep(i), &                                                !   !
             .false. )                                                 !   !
        !                                                              !   !
     case (2) ! -------------------------------------------------------+   !
        nbcp = (poly(i)%ptr%degr(1) + 1) * (poly(i)%ptr%degr(2) + 1)   !   !
        allocate( sep(i)%mat(nbcp,3) )                                 !   !
        call rearrange_for_separability_test( &                        !   !
             reshape( poly(i)%ptr%coef(&                               !   !
             1:poly(i)%ptr%degr(1)+1,&                                 !   !
             1:poly(i)%ptr%degr(2)+1,&                                 !   !
             1:3), [nbcp,3] ), &                                       !   !
             nbcp, &                                                   !   !
             xyz, &                                                    !   !
             sep(i)%mat, &                                             !   !
             nsep(i), &                                                !   !
             .false. )                                                 !   !
        !                                                              !   !
     case default ! ---------------------------------------------------+   !
        STOP 'intersect_elsewhere : poly%nvar /= 1,2'                  !   !
     end select ! <----------------------------------------------------+   !
     !                                                                     !
     if ( nsep(i) >= nbcp ) then
        PRINT *,'******* intersect_elsewhere : NSEP > NBCP *******'
     end if
  end do ! <---------------------------------------------------------------+


  ! perform spatial separability test
  if ( any(nsep < 1) ) then ! <--------------------------------------------+
     !                                                                     !
     separable = .true.                                                    !
     !                                                                     !
  else ! ------------------------------------------------------------------+
     if ( present(randomize) ) then ! <-----------------------------+      !
        if ( randomize ) then ! <-------------------------------+   !      !
           call random_rotation_matrix3d( rot )                 !   !      !
           do i = 1,2 ! <-----------------------------------+   !   !      !
              sep(i)%mat(1:nsep(i),1:3) = &                 !   !   !      !
                   matmul( sep(i)%mat(1:nsep(i),1:3), rot ) !   !   !      !
           end do ! <---------------------------------------+   !   !      !
        end if ! <----------------------------------------------+   !      !
     end if ! <-----------------------------------------------------+      !
     !                                                                     !
     call separating_plane( &                                              !
          sep(1)%mat(1:nsep(1),1:3), &                                     !
          sep(2)%mat(1:nsep(2),1:3), &                                     !
          nsep(1), &                                                       !
          nsep(2), &                                                       !
          vec, &                                                           !
          separable )                                                      !
     !                                                                     !
  end if ! <---------------------------------------------------------------+
  
  IF ( DEBUG ) THEN
     !IF (.NOT.SEPARABLE) THEN
     PRINT *,'XYZ =',xyz
     IF (SEPARABLE) THEN
        PRINT *,'VEC = ',matmul(rot, vec)
     END IF
     CALL WRITE_POLYNOMIAL(POLY(1)%PTR, 'dev_intersection/sep1_poly.bern')
     CALL WRITE_POLYNOMIAL(POLY(2)%PTR, 'dev_intersection/sep2_poly.bern')
     CALL WRITE_MATRIX(SEP(1)%MAT(1:NSEP(1),1:3), NSEP(1), 3, 'dev_intersection/sep1.dat')
     CALL WRITE_MATRIX(SEP(2)%MAT(1:NSEP(2),1:3), NSEP(2), 3, 'dev_intersection/sep2.dat')
     do i = 1,2
        IF ( NSEP(I) >= SIZE(SEP(I)%MAT,1) ) THEN
           PRINT *,'I =',I
           PRINT *,'XYZ =',xyz
           STOP 'DEBUG SEPARATION'
        END IF
     END do
     !END IF
  END IF

  deallocate( sep(1)%mat, sep(2)%mat )

end subroutine intersect_elsewhere
