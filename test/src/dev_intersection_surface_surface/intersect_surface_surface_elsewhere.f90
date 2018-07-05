subroutine intersect_surface_surface_elsewhere( &
     b1, &
     b2, &
     xyzinter, &
     separable, &
     randomize )
  use mod_math
  use mod_polynomial
  use mod_geometry
  use mod_separation
  implicit none
  type(type_polynomial), intent(in)  :: b1
  type(type_polynomial), intent(in)  :: b2
  real(kind=fp),         intent(in)  :: xyzinter(3)
  logical,               intent(out) :: separable
  logical, optional,     intent(in)  :: randomize
  real(kind=fp)                      :: sep1((b1%degr(1)+1)*(b1%degr(2)+1),3)
  real(kind=fp)                      :: sep2((b2%degr(1)+1)*(b2%degr(2)+1),3)
  real(kind=fp)                      :: vec(3), rot(3,3)
  integer                            :: nbcps, n1, n2

  nbcps = ( b1%degr(1) + 1 ) * ( b1%degr(2) + 1 )
  call rearrange_for_separability_test( &
       reshape( b1%coef(1:b1%degr(1)+1,1:b1%degr(2)+1,1:3), [nbcps,3] ), &
       nbcps, &
       xyzinter, & ! pivot
       xyzinter, & ! origin
       sep1, &
       n1 )

  nbcps = ( b2%degr(1) + 1 ) * ( b2%degr(2) + 1 )
  call rearrange_for_separability_test( &
       reshape( b2%coef(1:b2%degr(1)+1,1:b2%degr(2)+1,1:3), [nbcps,3] ), &
       nbcps, &
       xyzinter, & ! pivot
       xyzinter, & ! origin
       sep2, &
       n2 )

  if ( n1 < 1 .or. n2 < 1 ) then
     separable = .true.
  else
     if ( present(randomize) ) then
        if ( randomize ) then
           call random_rotation_matrix3d( rot )
           IF ( .false. ) THEN
              PRINT *,'';PRINT *,'';PRINT *,''
              PRINT *,'RANDOM ROTATION MATRIX ='
              CALL PRINT_MAT( ROT )
           END IF
           sep1(1:n1,:) = matmul( sep1(1:n1,:), rot )
           sep2(1:n2,:) = matmul( sep2(1:n2,:), rot )
        end if
     end if

     call separating_plane( &
          sep1(1:n1,1:3), &
          sep2(1:n2,1:3), &
          n1, &
          n2, &
          vec, &
          separable )

     IF ( .NOT.SEPARABLE ) THEN
        !PRINT *,'';PRINT *,'';PRINT *,''
        !PRINT *,'-----------------------------'
        !PRINT *,'XYZ_SEP (1) ='
        !CALL PRINT_MAT( SEP1(1:N1,:) )
        !PRINT *,'XYZ_SEP (2) ='
        !CALL PRINT_MAT( SEP2(1:N2,:) )     
        !PRINT *,'-----------------------------'
        CALL WRITE_POLYNOMIAL( B1, 'dev_intersection_surface_surface/bpn1.bern' )
        CALL WRITE_POLYNOMIAL( B2, 'dev_intersection_surface_surface/bpn2.bern' )
        CALL WRITE_MATRIX( SEP1(1:N1,1:3), N1, 3, 'dev_intersection_surface_surface/sep1.dat' )
        CALL WRITE_MATRIX( SEP2(1:N2,1:3), N2, 3, 'dev_intersection_surface_surface/sep2.dat' )
        CALL WRITE_MATRIX( ROT, 3, 3, 'dev_intersection_simple_surface/rot.dat' )
     END IF
  end if

end subroutine intersect_surface_surface_elsewhere
