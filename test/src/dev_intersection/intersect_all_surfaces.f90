subroutine intersect_all_surfaces( &
     surf, &
     nsurf, &
     interdata_global, &
     mask )
  USE MOD_UTIL
  use mod_math
  use mod_polynomial
  use mod_diffgeom
  use mod_regiontree
  use mod_types_intersection
  implicit none
  LOGICAL, PARAMETER :: DEBUG = ( GLOBALDEBUG .AND. .TRUE. )
  CHARACTER :: STRNUM
  integer,                      intent(in)    :: nsurf            ! number of surfaces
  type(type_surface), target,   intent(in)    :: surf(nsurf)      ! surfaces
  type(type_intersection_data), intent(inout) :: interdata_global ! global intersection data collection
  logical, optional,            intent(in)    :: mask(nsurf)      ! used to skip some surfaces when computing intersections
  type(ptr_surface)                           :: surfpair(2)   
  type(type_region), target                   :: root(nsurf)
  type(ptr_region)                            :: region(2)
  type(type_intersection_data)                :: interdata_local
  real(kind=fp), allocatable                  :: uvxyz(:,:)
  integer                                     :: nuvxyz
  integer                                     :: stat_degeneracy
  integer                                     :: isurf, jsurf, ic, i

  ! initialize all region trees
  do isurf = 1,nsurf
     call init_region( &
       root(isurf), &
       2, &
       [( [-1._fp, 1._fp], i=1,2 )] ) 

     ! compute Bezier control points for position and pseudonormal vectors
     allocate(root(isurf)%poly(2))
     allocate(root(isurf)%poly(1)%ptr, root(isurf)%poly(2)%ptr)
     call cheb2bern( &
          surf(isurf)%x, &
          root(isurf)%poly(1)%ptr )
     call cheb2bern( &
          surf(isurf)%pn, &
          root(isurf)%poly(2)%ptr )

     !IF ( DEBUG ) THEN
     !   WRITE (STRNUM,'(I1)') ISURF
     !   CALL WRITE_POLYNOMIAL( root(isurf)%poly(1)%ptr, 'dev_intersection/surfroot' // strnum // '_x.bern'  )
     !   CALL WRITE_POLYNOMIAL( root(isurf)%poly(2)%ptr, 'dev_intersection/surfroot' // strnum // '_pn.bern' )
     !END IF
  end do


  nuvxyz = 0
  allocate(uvxyz(7,10))

  ! loop over all pairs of DISTINCT surfaces and compute their intersection
  outer : do isurf = 1,nsurf-1
     if ( present(mask) ) then
        if ( .not.mask(isurf) ) cycle outer
     end if

     inner : do jsurf = isurf+1,nsurf
        if ( present(mask) ) then
           if ( .not.mask(jsurf) ) cycle inner
        end if

        ! check if there is a previously discovered tangential intersection curve between 
        ! the two intersected surfaces
        do ic = 1,interdata_global%nc
           if ( .not.interdata_global%curves(ic)%smooth ) cycle
           if ( ( &
                associated( interdata_global%curves(ic)%surf(1)%ptr, surf(isurf) ) .and. &
                associated( interdata_global%curves(ic)%surf(2)%ptr, surf(jsurf) ) ) .or. &
                ( &
                associated( interdata_global%curves(ic)%surf(1)%ptr, surf(jsurf) ) .and. &
                associated( interdata_global%curves(ic)%surf(2)%ptr, surf(isurf) ) ) &
                ) cycle inner ! we skip this pair of surfaces
        end do

        IF ( DEBUG ) THEN
           PRINT *,''; PRINT *,''; PRINT *,''
           PRINT *,'PAIR :',ISURF,JSURF
        END IF

        ! initialize pointers to surfaces and region trees
        region(1)%ptr   => root(isurf)
        region(2)%ptr   => root(jsurf)
        surfpair(1)%ptr => surf(isurf)
        surfpair(2)%ptr => surf(jsurf)
        stat_degeneracy =  0

        ! if some isolated or branch singular intersection points between 
        ! the two intersected surfaces have already been discovered, copy them into
        ! the local intersection data collection
        ! (...)        

        allocate(interdata_local%points(10), interdata_local%curves(10))
        interdata_local%np = 0
        interdata_local%nc = 0
        
        ! compute topology of the intersection between the current pair of surfaces
        nuvxyz = 0
        call intersect_surface_pair( &
             surfpair, &
             region, &
             interdata_local, &
             uvxyz, &
             nuvxyz, &
             stat_degeneracy ) 
        
        IF ( DEBUG ) THEN
           PRINT *,'STAT_DEGENERACY =',stat_degeneracy
           IF ( NUVXYZ > 0 ) THEN
              CALL WRITE_MATRIX( TRANSPOSE(UVXYZ(1:7,1:NUVXYZ)), NUVXYZ, 7, &
                   'dev_intersection/uvxyz.dat' )
           END IF
           CALL WRITE_INTERSECTION_DATA( INTERDATA_LOCAL, &
                'dev_intersection/interdata_points.dat', &
                'dev_intersection/interdata_curves.dat' )

           PRINT *,NUVXYZ,            ' INTERSECTION POINT(S)'
           PRINT *,INTERDATA_LOCAL%NC,' INTERSECTION CURVE(S)'

           !DO i = 1,2
           !   WRITE (STRNUM,'(I1)') i
           !   CALL EXPORT_REGION_TREE( REGION(i)%PTR, 'dev_intersection/treessi_' // strnum // '.dat' )
           !END DO
        END IF

        ! trace intersection curves and append them to the global intersection data collection
        call merge_intersection_data( &
             surfpair, &
             uvxyz, &
             nuvxyz, &
             interdata_local, &
             interdata_global )

        ! reset local intersection data collection
        call free_intersection_data(interdata_local)
        call free_ipts(region(1)%ptr)
        call free_ipts(region(2)%ptr)

        ! if a degeneracy has been encountered, report it
        if ( stat_degeneracy > 0 ) exit outer

     end do inner
  end do outer

  if ( allocated(uvxyz) ) deallocate(uvxyz)


  ! free all region trees
  do isurf = 1,nsurf
     IF ( DEBUG ) THEN
        WRITE (STRNUM,'(I1)') isurf
        CALL EXPORT_REGION_TREE( root(isurf), 'dev_intersection/treessi_' // strnum // '.dat' )
     END IF

     call free_polynomial(root(isurf)%poly(1)%ptr)
     call free_polynomial(root(isurf)%poly(2)%ptr)
     deallocate(root(isurf)%poly(1)%ptr, root(isurf)%poly(2)%ptr)
     deallocate(root(isurf)%poly)
     call free_region_tree(root(isurf)) 
  end do

  nullify(region(1)%ptr, region(2)%ptr)

end subroutine intersect_all_surfaces
