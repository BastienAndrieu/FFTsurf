subroutine loop_detection_criterion( &
     bpn, &
     stat, &
     param_vector, &
     known_to_intersect, &
     n_collineal, &
     randomize )
  USE MOD_UTIL
  use mod_math
  use mod_geometry
  use mod_polynomial
  use mod_separation
  ! Returns:  stat = 0 if Hohmeyer's criterion is satisfied
  !                = 1 if Gauss map #1 must be subdivided
  !                = 2 if Gauss map #2 must be subdivided
  !                = 3 if both Gauss maps must be subdivided
  implicit none
  LOGICAL, PARAMETER :: DEBUG = ( GLOBALDEBUG .AND. .false. )
  CHARACTER :: STRNUM
  type(ptr_polynomial), intent(in)  :: bpn(2)
  integer,              intent(out) :: stat
  real(kind=fp),        intent(out) :: param_vector(3)
  logical,              intent(in)  :: known_to_intersect
  real(kind=fp),        intent(in)  :: n_collineal(3)
  logical, optional,    intent(in)  :: randomize
  type(type_matrix)                 :: sep(2)
  integer                           :: nbcp, nsep(2)
  real(kind=fp)                     :: rot(3,3)
  logical                           :: separable
  real(kind=fp)                     :: vec1(3), vec2(3), vecsqr
  real(kind=fp)                     :: gaussmapsize(2), coneaxe(3)

 IF ( DEBUG ) PRINT *,'KNOWN TO INTERSECT?',known_to_intersect

 if ( present(randomize) ) then ! <-----------------------------+
    if ( randomize ) then ! <-------------------------------+   !
       call random_rotation_matrix3d( rot )                 !   !
    end if ! <----------------------------------------------+   !
 end if ! <-----------------------------------------------------+

  do isurf = 1,2
     nbcp = (bpn(isurf)%ptr%degr(1) + 1) * (bpn(isurf)%ptr%degr(2) + 1)
     allocate( sep(isurf)%mat(nbcp,3) )
     call rearrange_for_separability_test( &
          reshape( bpn(isurf)%ptr%coef(&
          1:bpn(isurf)%ptr%degr(1)+1,&
          1:bpn(isurf)%ptr%degr(2)+1,&
          1:3), [nbcp,3] ), &
          nbcp, &
          n_collineal, &
          sep(isurf)%mat, &
          nsep(isurf), &
          .true. )
     if ( known_to_intersect ) then
        nsep(isurf) = nsep(isurf) + 1
        sep(isurf)%mat(nsep(isurf),1:3) = n_collineal
     end if
     if ( present(randomize) ) then ! <------------------------------+
        if ( randomize ) then ! <--------------------------------+   !
           sep(isurf)%mat(1:nsep(isurf),1:3) = &                 !   !
                matmul( sep(isurf)%mat(1:nsep(isurf),1:3), rot ) !   !
        end if ! <-----------------------------------------------+   !
     end if ! <------------------------------------------------------+

     IF ( DEBUG ) THEN
        WRITE (STRNUM,'(I1)') ISURF
        CALL WRITE_MATRIX( SEP(ISURF)%MAT(1:NSEP(ISURF),1:3), NSEP(ISURF), 3, &
             'dev_intersection/sep' // strnum // '.dat' )
     END IF
  end do

  ! first, we search for a vector vec1 such that vec1.n1 < 0 and vec1.n2 > 0
  ! for all ni in the Gauss map of surface #i
  if ( known_to_intersect ) then ! <-------------------------------------------+
     ! if the two gauss maps are known to intersect at a point                 !
     ! (with common normal direction n_collineal), then vec1                   !
     ! only has 1 degree of freedom                                            !
     call intersect_gaussmaps_elsewhere( &                                     !
          sep, &                                                               !
          nsep-1, &                                                              !
          n_collineal, &                                                       !
          separable, &                                                         !
          stat, &                                                              !
          vec1 )                                                               !
  else ! ----------------------------------------------------------------------+
     call separating_plane( &                                                  !
          sep(1)%mat(1:nsep(1),1:3), &                                         !
          sep(2)%mat(1:nsep(2),1:3), &                                         !
          nsep(1), &                                                           !
          nsep(2), &                                                           !
          vec1, &                                                              !
          separable )                                                          !
  end if ! <-------------------------------------------------------------------+
  
  IF ( DEBUG ) THEN
     PRINT *,'P1?',SEPARABLE
     IF ( SEPARABLE ) PRINT *,vec1
  END IF

  if ( separable ) then ! <----------------------------------------------------+
     ! then we search for a vector vec2 such that vec2.n1 > 0 and vec2.n2 > 0  !
     call separating_plane( &                                                  !
          sep(1)%mat(1:nsep(1),1:3), &                                         !
          -sep(2)%mat(1:nsep(2),1:3), &                                        !
          nsep(1), &                                                           !
          nsep(2), &                                                           !
          vec2, &                                                              !
          separable )                                                          !
     IF ( DEBUG ) THEN
        PRINT *,'P2?',SEPARABLE
        IF ( SEPARABLE ) PRINT *,vec2
     END IF
  end if ! <-------------------------------------------------------------------+
     

  if ( separable ) then ! <----------------------------------------------------+
     param_vector = cross( vec1, vec2 )                                        !
     vecsqr = sum( param_vector**2 )                                           !
     if ( vecsqr > EPSfpsqr ) then ! <----------------+                        !
        stat = 0                                      !                        !
        param_vector = param_vector / sqrt( vecsqr )  !                        !
        if ( present(randomize) ) then ! <--------+
           if ( randomize ) then ! <-----------+  !
              param_vector = &                 !  !
                   matmul( rot, param_vector ) !  !
           end if ! <--------------------------+  !
        end if ! <--------------------------------+
        IF ( DEBUG ) THEN
           PRINT *,'P = ',param_vector
           DO ISURF = 1,2
              WRITE (STRNUM,'(I1)') ISURF
              CALL WRITE_MATRIX( SEP(ISURF)%MAT(1:NSEP(ISURF),1:3), NSEP(ISURF), 3, &
                   'dev_intersection/sep' // strnum // '.dat' )
              CALL WRITE_POLYNOMIAL( BPN(ISURF)%PTR, 'dev_intersection/debugld_reg' // strnum // '.bern' )
           END DO
           STOP
        END IF
     else ! ------------------------------------------+                        !
        ! an error has occured                        !                        !
        stat = -1                                     !                        !
        separable = .false.                           !                        !
     end if ! <---------------------------------------+                        !
  end if ! <-------------------------------------------------------------------+


  if ( .not.separable ) then ! <-----------------------------------------------+
     ! which surface has the "largest" Gauss map?                              !
     ! we use an approximate solid angle as a surrogate for Gauss map size     !
     do isurf = 1,2 ! <--------------------------------------------------+     !
        call bounding_cone( &                                            !     !
          sep(isurf)%mat(1:nsep(isurf),1:3), &                           !     !
          coneaxe, &                                                     !     !
          gaussmapsize(isurf) )                                          !     !
     end do ! <----------------------------------------------------------+     !
     !                                                                         !
     if ( all(gaussmapsize > 0.4_fp * CSTpi) ) then ! <------------------+     !
        stat = 3                                                         !     !
     else ! -------------------------------------------------------------+     !
        stat = maxloc(gaussmapsize, dim=1)                               !     !
     end if ! <----------------------------------------------------------+     !
  end if ! <-------------------------------------------------------------------+

  deallocate( sep(1)%mat, sep(2)%mat )

end subroutine loop_detection_criterion
