subroutine hohmeyer_loop_detection_bezier( &
    bpn1, &
    n1, &
    bpn2, &
    n2, &
    vec, &
    stat ) 
  use mod_obb
  use mod_separation
  ! Returns :
  !  'stat' = 0 if Hohmeyer's criterion is satisfied
  !         = 1 if Gauss map #1 must be subdivided
  !         = 2 if Gauss map #2 must be subdivided
  !         = 3 if both Gauss maps must be subdivided
  implicit none
  logical, parameter             :: VERBOSE = .false.
  integer                        :: n1, n2
  real(kind=MATHpr), intent(in)  :: bpn1(n1,3), bpn2(n2,3)
  real(kind=MATHpr), intent(out) :: vec(3)
  integer,           intent(out) :: stat
  logical                        :: separable
  real(kind=MATHpr)              :: vec1(3), vec2(3), normvec
  real(kind=MATHpr)              :: gaussmapsize(2), coneaxe(3)

  !PRINT *,'BPN1='
  !CALL PRINT_MAT(REAL(BPN1))
  !PRINT *,'BPN2='
  !CALL PRINT_MAT(REAL(BPN2))
  
  ! search for a vector vec1 such that vec1.n1 < 0 and vec1.n2 > 0
  ! for all ni in pnbox(i)
  call separating_plane( &
       bpn1, &
       bpn2, &
       vec1, &
       separable )
  IF (VERBOSE) PRINT *,separable

  if ( separable ) then
     IF (VERBOSE) PRINT *,'VEC1 =',REAL(VEC1)
     ! search for a vector vec2 such that vec2.n1 > 0 and vec2.n2 > 0 
     call separating_plane( &
          bpn1, &
          -bpn2, &
          vec2, &
          separable )
     PRINT *,separable

  end if

  if ( separable ) then
     IF (VERBOSE) PRINT *,'VEC2 =',REAL(VEC2)
     vec = cross( vec1, vec2 )
     normvec = norm2( vec )
     if ( normvec > MATHeps ) then
        vec = vec / normvec
        return
     end if
  end if

  ! which surface has the "largest" Gauss map?
  ! we use an approximate solid angle as a surrogate for Gauss map size
  call bounding_cone( &
       bpn1, &
       coneaxe, &
       gaussmapsize(1) )
  call bounding_cone( &
       bpn2, &
       coneaxe, &
       gaussmapsize(2) )
  IF (VERBOSE) PRINT *,'   GAUSS MAP SIZE / PI =',REAL( GAUSSMAPSIZE / MATHPI )

  if ( all( gaussmapsize > 0.4_MATHpr * MATHpi ) ) then
     stat = 3
  else
     stat = maxloc( gaussmapsize, 1 )
  end if

end subroutine hohmeyer_loop_detection_bezier

