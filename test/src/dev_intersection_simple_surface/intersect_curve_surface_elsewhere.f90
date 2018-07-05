subroutine intersect_curve_surface_elsewhere( &
     b_c, &
     b_s, &
     xyzinter, &
     separable, &
     randomize )
  use mod_math
  use mod_polynomial
  use mod_geometry
  use mod_separation
  implicit none
  type(type_polynomial), intent(in)  :: b_c
  type(type_polynomial), intent(in)  :: b_s
  real(kind=fp),         intent(in)  :: xyzinter(3)
  logical,               intent(out) :: separable
  logical, optional,     intent(in)  :: randomize
  real(kind=fp)                      :: sep_c(b_c%degr(1)+1,3)
  real(kind=fp)                      :: sep_s((b_s%degr(1)+1)*(b_s%degr(2)+1),3)
  real(kind=fp)                      :: vec(3), rot(3,3)
  integer                            :: nbcps, nc, ns

  call rearrange_for_separability_test( &
       b_c%coef(1:b_c%degr(1)+1,1:3,1), &
       b_c%degr(1)+1, &
       xyzinter, & ! pivot
       xyzinter, & ! origin
       sep_c, &
       nc )

  nbcps = ( b_s%degr(1) + 1 ) * ( b_s%degr(2) + 1 )
  call rearrange_for_separability_test( &
       reshape( b_s%coef(1:b_s%degr(1)+1,1:b_s%degr(2)+1,1:3), [nbcps,3] ), &
       nbcps, &
       xyzinter, & ! pivot
       xyzinter, & ! origin
       sep_s, &
       ns )

  if ( nc < 1 .or. ns < 1 ) then
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
           sep_c(1:nc,:) = matmul( sep_c(1:nc,:), rot )
           sep_s(1:ns,:) = matmul( sep_s(1:ns,:), rot )
        end if
     end if

     call separating_plane( &
          sep_c(1:nc,1:3), &
          sep_s(1:ns,1:3), &
          nc, &
          ns, &
          vec, &
          separable )

     IF ( .NOT.SEPARABLE ) THEN
        !PRINT *,'';PRINT *,'';PRINT *,''
        !PRINT *,'-----------------------------'
        !PRINT *,'XYZ_SEP (C) ='
        !CALL PRINT_MAT( SEP_C(1:NC,:) )
        !PRINT *,'XYZ_SEP (S) ='
        !CALL PRINT_MAT( SEP_S(1:NS,:) )     
        !PRINT *,'-----------------------------'
        CALL WRITE_MATRIX( SEP_C(1:NC,1:3), NC, 3, 'dev_intersection_simple_surface/sepc.dat' )
        CALL WRITE_MATRIX( SEP_S(1:NS,1:3), NS, 3, 'dev_intersection_simple_surface/seps.dat' )
        CALL WRITE_MATRIX( ROT, 3, 3, 'dev_intersection_simple_surface/rot.dat' )
     END IF
  end if

end subroutine intersect_curve_surface_elsewhere
