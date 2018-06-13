subroutine hohmeyer_loop_detection( &
    pnbox, &
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
  type(type_obb),    intent(in)  :: pnbox(2)
  real(kind=MATHpr), intent(out) :: vec(3)
  integer,           intent(out) :: stat
  real(kind=MATHpr)              :: xyzsep(16,3)
  logical                        :: separable
  real(kind=MATHpr)              :: vec1(3), vec2(3), normvec
  real(kind=MATHpr)              :: gaussmapsize(2), coneaxe(3)
  integer                        :: isurf, i

  stat = 0

  ! check if Gauss map fits in one single hemisphere
  ! (we check if the bounding box contains the origin)
  do isurf = 1,2
     if ( is_inside_OBB( &
          real( [0,0,0], kind=MATHpr ), &
          pnbox(isurf) ) ) then
        IF (VERBOSE) PRINT *,'   GAUSS MAP #',ISURF,' CONTAINS ORIGIN'
        stat = stat + isurf
     end if
  end do

  !if ( any(subdivsurf) ) return
  if ( stat > 0 ) return

  ! get xyz coordinates of the vertices of both bounding boxes
  do isurf = 1,2
     call OBB_vertices( &
          pnbox(isurf), &
          xyzsep(8*(isurf-1)+1:8*isurf,:) )
  end do

  ! search for a vector vec1 such that vec1.n1 < 0 and vec1.n2 > 0
  ! for all ni in pnbox(i)
  call separating_plane( &
       xyzsep(1:8,:), &
       xyzsep(9:16,:), &
       8, &
       8, &
       vec1, &
       separable )
  IF (VERBOSE) PRINT *,separable

  if ( separable ) then
     IF (VERBOSE) PRINT *,'VEC1 =',REAL(VEC1)
     ! search for a vector vec2 such that vec2.n1 > 0 and vec2.n2 > 0 
     call separating_plane( &
          xyzsep(1:8,:), &
          -xyzsep(9:16,:), &
          8, &
          8, &
          vec2, &
          separable )
     !PRINT *,separable

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
  xyzsep = xyzsep / spread( sqrt(sum(xyzsep**2,2)), dim=2, ncopies=3 )
  gaussmapsize(:) = 0._MATHpr
  do isurf = 1,2
     coneaxe = pnbox(isurf)%ctr
     coneaxe = coneaxe / norm2( coneaxe )
     do i = (isurf-1)*8+1,isurf*8
        gaussmapsize(isurf) = max( gaussmapsize(isurf), &
             acos( &
             min( 1._MATHpr, max( -1._MATHpr, dot_product( coneaxe, xyzsep(i,:) ) ) ) &
             ) )
        
     end do
  end do
  IF (VERBOSE) PRINT *,'   GAUSS MAP SIZE / PI =',REAL( GAUSSMAPSIZE / MATHPI )

  if ( all( gaussmapsize > 0.4_MATHpr * MATHpi ) ) then
     stat = 3
  else
     stat = maxloc( gaussmapsize, 1 )
  end if

end subroutine hohmeyer_loop_detection
