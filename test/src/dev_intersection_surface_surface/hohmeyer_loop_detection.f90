subroutine hohmeyer_loop_detection( &
     bpn1, &
     bpn2, &
     vec, &
     stat, &
     known_to_intersect, &
     normal_collineal ) 
  use mod_obb
  use mod_separation
  ! Returns :
  !  'stat' = 0 if Hohmeyer's criterion is satisfied
  !         = 1 if Gauss map #1 must be subdivided
  !         = 2 if Gauss map #2 must be subdivided
  !         = 3 if both Gauss maps must be subdivided
  implicit none
  logical, parameter         :: DEBUG = .true.
  logical, parameter         :: VERBOSE = .false.
  type(type_polynomial), intent(in)  :: bpn1, bpn2
  real(kind=fp),         intent(out) :: vec(3)
  integer,               intent(out) :: stat
  logical,               intent(in)  :: known_to_intersect
  real(kind=fp),         intent(in)  :: normal_collineal(3)
  real(kind=fp),         allocatable :: sep1(:,:), sep2(:,:)
  integer                            :: n1, n2
  logical                            :: separable
  real(kind=fp)                      :: vec1(3), vec2(3), normvec
  real(kind=fp)                      :: gaussmapsize(2), coneaxe(3)

  IF ( DEBUG ) PRINT *,'KNOWN TO INTERSECT?',known_to_intersect

  n1 = product( bpn1%degr(1:2) + 1 )
  n2 = product( bpn2%degr(1:2) + 1 )
  allocate( sep1(n1,3), sep2(n2,3) )
  sep1 = reshape( bpn1%coef(1:bpn1%degr(1)+1,1:bpn1%degr(2)+1,1:3), [n1,3] )
  sep2 = reshape( bpn2%coef(1:bpn2%degr(1)+1,1:bpn2%degr(2)+1,1:3), [n2,3] )
  
  sep1 = sep1 / spread( sqrt( sum( sep1**2, 2 ) ), dim=2, ncopies=3 )
  sep2 = sep2 / spread( sqrt( sum( sep2**2, 2 ) ), dim=2, ncopies=3 )

  ! search for a vector vec1 such that vec1.n1 < 0 and vec1.n2 > 0
  ! for all ni in pnbox(i)
  if ( known_to_intersect ) then
     call intersect_gaussmaps_elsewhere( &
          bpn1, &
          bpn2, &
          normal_collineal, &
          vec1, &
          separable )
     !PRINT *,'HKTI SEPARABLE?',separable
     !IF ( .NOT.SEPARABLE ) THEN
        !CALL WRITE_POLYNOMIAL( BPN1, 'dev_intersection_surface_surface/bpn1.bern' )
        !CALL WRITE_POLYNOMIAL( BPN2, 'dev_intersection_surface_surface/bpn2.bern' )
        !CALL WRITE_MATRIX( normal_collineal, 1, 3, 'dev_intersection_surface_surface/norcol.dat' )
     !   STOP '>>> DEBUG'
     !END IF
  else
     call separating_plane( &
          sep1(1:n1,1:3), &
          sep2(1:n2,1:3), &
          n1, &
          n2, &
          vec1, &
          separable )

     IF ( .false. ) THEN
        PRINT *,'SEPARABLE?',separable
        CALL WRITE_MATRIX( SEP1(1:N1,1:3), N1, 3, 'trace_intersection_polyline/sep1.dat' )
        CALL WRITE_MATRIX( SEP2(1:N2,1:3), N2, 3, 'trace_intersection_polyline/sep2.dat' )
        PRINT *,'***'
        STOP
     END IF

  end if
  !IF (VERBOSE) PRINT *,separable

  IF ( DEBUG ) THEN
     PRINT *,'P1?',SEPARABLE
     IF ( SEPARABLE ) PRINT *,vec1
  END IF

  if ( separable ) then
     IF (VERBOSE) PRINT *,'VEC1 =',REAL(VEC1)
     ! search for a vector vec2 such that vec2.n1 > 0 and vec2.n2 > 0 
     call separating_plane( &
          sep1(1:n1,1:3), &
          -sep2(1:n2,1:3), &
          n1, &
          n2, &
          vec2, &
          separable )
     IF ( DEBUG ) THEN
        PRINT *,'P2?',SEPARABLE
        IF ( SEPARABLE ) PRINT *,vec2
     END IF
  end if

  if ( separable ) then
     IF (VERBOSE) PRINT *,'VEC2 =',REAL(VEC2)

     vec = cross( vec1, vec2 )
     normvec = norm2( vec )
     if ( normvec < EPSmath ) then
        PRINT *,'VEC1 =', VEC1
        PRINT *,'VEC2 =', VEC2
        stat = -1
        return
        !STOP 'hohmeyer_loop_detection : normvec << 1'
     else
        vec = vec / normvec
        stat = 0
        IF ( DEBUG ) PRINT *,'P = ',vec
        IF ( VERBOSE ) THEN
           PRINT *,'KNOWN TO INTERSECT?',known_to_intersect
           PRINT *,'VEC1 =', VEC1
           PRINT *,'VEC2 =', VEC2
           PRINT *,'   P =', VEC
        END IF
     end if
  else
     ! which surface has the "largest" Gauss map?
     ! we use an approximate solid angle as a surrogate for Gauss map size
     call bounding_cone( &
          sep1(1:n1,1:3), &
          coneaxe, &
          gaussmapsize(1) )
     call bounding_cone( &
          sep2(1:n2,1:3), &
          coneaxe, &
          gaussmapsize(2) )
     IF (VERBOSE) PRINT *,'   GAUSS MAP SIZE / PI =',REAL( GAUSSMAPSIZE / CSTPI )

     if ( all( gaussmapsize > 0.4_fp * MATHpi ) ) then
        stat = 3
     else
        stat = maxloc( gaussmapsize, 1 )
     end if
  end if

  deallocate( sep1, sep2 )

end subroutine hohmeyer_loop_detection
