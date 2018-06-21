subroutine intersect_gaussmaps_elsewhere( &
     bpn1, &
     bpn2, &
     xyzinter, &
     vec, &
     separable )
  use mod_math
  use mod_polynomial
  use mod_separation
  implicit none
  LOGICAL, PARAMETER :: DEBUG = .false.
  type(type_polynomial), intent(in)  :: bpn1, bpn2
  real(kind=fp),         intent(in)  :: xyzinter(3)
  logical,               intent(out) :: separable
  real(kind=fp),         intent(out) :: vec(3)
  real(kind=fp)                      :: sep1((bpn1%degr(1)+1)*(bpn1%degr(2)+1),3)
  real(kind=fp)                      :: sep2((bpn2%degr(1)+1)*(bpn2%degr(2)+1),3)
  integer                            :: n1, n2
  integer                            :: nbcp
  real(kind=fp)                      :: rot(3,3), wedge(2)

  nbcp = ( bpn1%degr(1) + 1 ) * ( bpn1%degr(2) + 1 )
  call rearrange_for_separability_test( &
       reshape( bpn1%coef(1:bpn1%degr(1)+1,1:bpn1%degr(2)+1,1:bpn1%dim), [nbcp,3] ), &
       nbcp, &
       xyzinter, &              ! pivot
       [0._fp, 0._fp, 0._fp], & ! origin
       sep1, &
       n1 )

  nbcp = ( bpn2%degr(1) + 1 ) * ( bpn2%degr(2) + 1 )
  call rearrange_for_separability_test( &
       reshape( bpn2%coef(1:bpn2%degr(1)+1,1:bpn2%degr(2)+1,1:bpn2%dim), [nbcp,3] ), &
       nbcp, &
       xyzinter, &              ! pivot
       [0._fp, 0._fp, 0._fp], & ! origin
       sep2, &
       n2 )

  rot(:,1) = xyzinter! / norm2( xyzinter )
  rot(:,2) = 0._fp
  rot(minloc(abs(rot(:,1))),2) = 1._fp
  rot(:,2) = rot(:,2) - dot_product( rot(:,2), rot(:,1) ) * rot(:,1)
  rot(:,3) = cross( rot(:,1), rot(:,2) )

  sep1(1:n1,1:3) = matmul( sep1(1:n1,1:3), rot )
  sep2(1:n2,1:3) = matmul( sep2(1:n2,1:3), rot )


  if ( n1 < 1 .and. n2 < 1 ) then ! <---------------------------------------+
     !                                                                      !
     separable = .false.                                                    !
     !                                                                      !
  elseif ( n1 < 1 ) then ! -------------------------------------------------+
     !                                                                      !
     call minimal_bounding_sector( &                                        !
          sep2(1:n2,2:3), &                                                 !
          n2, &                                                             !
          wedge, &                                                          !
          separable )                                                       !
     if ( separable ) vec = [ 0._fp, cos(wedge(1)), sin(wedge(1)) ]         !
     !                                                                      !
  elseif ( n2 < 1 ) then ! -------------------------------------------------+ 
     !                                                                      !
     call minimal_bounding_sector( &                                        !
          sep1(1:n1,2:3), &                                                 !
          n1, &                                                             !
          wedge, &                                                          !
          separable )                                                       !
     if ( separable ) vec = [ 0._fp, cos(wedge(1)), sin(wedge(1)) ]         !
     !                                                                      !
  else ! -------------------------------------------------------------------+
     !                                                                      !
     call separate_spherical_bounding_boxes( &                              !
          sep1(1:n1,1:3), &                                                 !
          sep2(1:n2,1:3), &                                                 !
          n1, &                                                             !
          n2, &                                                             !
          vec, &                                                            !
          separable, &                                                      !
          mask_axes=[ .true., .false., .false. ] )                          !
     !                                                                      !
  end if ! <----------------------------------------------------------------+

  if ( separable ) vec = matmul( rot, vec ) 

  IF ( DEBUG ) THEN
     PRINT *,'MAY INTERSECT ELSEWHERE?',.NOT.SEPARABLE
     IF ( SEPARABLE ) PRINT *,' VEC =',VEC
  END IF

  !IF ( DEBUG .and. separable ) THEN
  IF ( .FALSE. ) THEN
     CALL WRITE_POLYNOMIAL( BPN1, 'dev_intersection_surface_surface/bpn1.bern' )
     CALL WRITE_POLYNOMIAL( BPN2, 'dev_intersection_surface_surface/bpn2.bern' )
     CALL WRITE_MATRIX( SEP1(1:N1,1:3), N1, 3, 'dev_intersection_surface_surface/sep1.dat' )
     CALL WRITE_MATRIX( SEP2(1:N2,1:3), N2, 3, 'dev_intersection_surface_surface/sep2.dat' )
     CALL WRITE_MATRIX( RESHAPE(VEC,[1,3]), 1, 3, 'dev_intersection_surface_surface/vec.dat' )
     CALL WRITE_MATRIX( RESHAPE(XYZINTER,[1,3]), 1, 3, 'dev_intersection_surface_surface/normal_collineal.dat' )
     PRINT *,'SEPARABLE =',SEPARABLE
     STOP '---> DEBUG_SEPARATION'
  END IF

end subroutine intersect_gaussmaps_elsewhere
