subroutine intersect_border_surface( &
     root_s, &
     root_c, &
     region, &
     icurv, &
     ivar, &
     ival, &
     uvxyz, &
     nuvxyz, &
     iptsbs, &
     nptsbs, &
     stat_degeneracy )
  use mod_math
  use mod_polynomial
  use mod_diffgeom2
  use mod_regiontree
  use mod_types_intersection
  use mod_tolerances
  implicit none
  LOGICAL, PARAMETER :: DEBUG = ( GLOBALDEBUG .AND. .false. )
  type(ptr_surface),          intent(in)    :: root_s(2)
  type(type_curve),           intent(in)    :: root_c
  type(ptr_region),           intent(inout) :: region(2)
  integer,                    intent(in)    :: icurv, ivar, ival
  real(kind=fp), allocatable, intent(inout) :: uvxyz(:,:)
  integer,                    intent(inout) :: nuvxyz
  integer, allocatable,       intent(inout) :: iptsbs(:)
  integer,                    intent(inout) :: nptsbs
  integer,                    intent(inout) :: stat_degeneracy
  integer                                   :: isurf, jvar
  type(type_region)                         :: region_c
  type(type_region)                         :: region_s
  real(kind=fp), allocatable                :: tuvxyz(:,:)
  integer                                   :: ntuvxyz, ntuvxyz_tmp
  integer, allocatable                      :: sharedpts(:)
  integer                                   :: n_sharedpts
  real(kind=fp)                             :: tmp(7), uv(2,2)
  integer                                   :: ipt, jpt

  IF ( DEBUG ) THEN
     PRINT *,''; PRINT *,'';
     PRINT *,'INTERSECT_BORDER_SURFACE'
     PRINT *,'ICURV, IVAR, IVAL =',ICURV,IVAR,IVAL
  END IF

  isurf = 1 + mod(icurv,2)
  jvar  = 1 + mod(ivar ,2)

  ! initialize a temporary curve region tree
  call init_region( &
       region_c, &
       1, &
       [ -1._fp, 1._fp ] ) 

  ! initialize a temporary surface region tree
  call init_region( &
       region_s, &
       2, &
       region(isurf)%ptr%uvbox ) 

  ! compute the curve's Bezier control points
  allocate(region_c%poly(1))
  allocate(region_c%poly(1)%ptr)
  call cheb2bern( &
       root_c%x, &
       region_c%poly(1)%ptr )

  ! copy the surface's Bezier control points
  allocate(region_s%poly(1))
  region_s%poly(1)%ptr => region(isurf)%ptr%poly(1)%ptr
  
  ntuvxyz = 0
  allocate(tuvxyz(6,10))
  ! get the list of already discovered intersection points contained in both surface regions *
  n_sharedpts = 0                                                                          ! *
  if ( region(1)%ptr%npts > 0 .and. region(2)%ptr%npts > 0 ) then ! <-----------------------+ PASSER
     call intersection_arrays( &                                                            ! DANS 
          region(1)%ptr%ipts(1:region(1)%ptr%npts), &                                       ! INTERSECT
          region(2)%ptr%ipts(1:region(2)%ptr%npts), &                                       ! SURFACE-
          sharedpts )                                                                       ! SURFACE
     if ( allocated(sharedpts) ) n_sharedpts = size(sharedpts)                              !*
  end if ! <--------------------------------------------------------------------------------+*
  IF ( DEBUG ) THEN
     PRINT *,'N_SHAREDPTS =',n_sharedpts
     IF ( N_SHAREDPTS > 0 ) THEN
        PRINT *,sharedpts(1:n_sharedpts)
        CALL PRINT_MAT( TRANSPOSE(UVXYZ(1:4,sharedpts(1:n_sharedpts))) )
     END IF
  END IF

  ! check if some of these points are located on the current border
  do jpt = 1,n_sharedpts ! <----------------------------------------------------------------+
     ipt = sharedpts(jpt)                                                                   !
     tmp(1) = uvxyz(2*(icurv-1)+ivar,ipt)                                                   !
     !                                                                                      !
     if ( abs( tmp(1) - region(icurv)%ptr%uvbox(2*(ivar-1)+ival) ) < EPSuv ) then ! <---+   !
        tmp(1) = ab2n1p1( &                                                             !   !
             uvxyz(2*(icurv-1)+jvar,ipt), &                                             !   !
             region(icurv)%ptr%uvbox(2*jvar-1), &                                       !   !
             region(icurv)%ptr%uvbox(2*jvar  ) )                                        !   !
        tmp(2:3) = uvxyz(2*isurf-1:2*isurf,ipt)                                         !   !
        tmp(4:6) = uvxyz(5:7,ipt)                                                       !   !
        IF ( DEBUG ) THEN
           PRINT *,'+ 1PT ON BORDER :'
           PRINT *,'    UV =',UVXYZ(1:4,IPT)
           PRINT *,'   XYZ =',UVXYZ(5:7,IPT)
           PRINT *,'TUVXYZ =',TMP(1:6)
        END IF
        !                                                                               !   !
        call append_vector( &                                                           !   !
             tmp(1:6), &                                                                !   !
             6, &                                                                       !   !
             tuvxyz, &                                                                  !   !
             ntuvxyz )                                                                  !   !
        !                                                                               !   !
        call append_n( &                                                                !   !
             iptsbs, &                                                                  !   !
             nptsbs, &                                                                  !   !
             [ipt], &                                                                   !   !
             1, &                                                                       !   !
             unique=.true. )                                                            !   !
     end if ! <-------------------------------------------------------------------------+   !
  end do ! <--------------------------------------------------------------------------------+
  ntuvxyz_tmp = ntuvxyz

  IF ( DEBUG ) THEN
     PRINT *,'NTUVXYZ_TMP =',NTUVXYZ!_TMP
     PRINT *,'TUVXYZ ='
     CALL PRINT_MAT( TRANSPOSE(TUVXYZ(1:6,1:NTUVXYZ)) )
  END IF


  ! compute the curve-surface intersection
  stat_degeneracy = 0
  !ntuvxyz = 0
  call intersect_curve_surface( &
       root_c, &
       root_s(isurf)%ptr, &
       region_c, &
       region_s, &
       tuvxyz, &
       ntuvxyz, &
       stat_degeneracy )

  IF ( DEBUG ) THEN
     PRINT *,'NTUVXYZ =',NTUVXYZ
     IF ( NTUVXYZ > 0 ) CALL PRINT_MAT( TRANSPOSE(TUVXYZ(:,1:NTUVXYZ)) )

     !CALL EXPORT_REGION_TREE( REGION_C, 'dev_intersection/treebsi_c.dat' )
     !CALL EXPORT_REGION_TREE( REGION_S, 'dev_intersection/treebsi_s.dat' )

  END IF

  ! free the curve region tree
  call free_polynomial(region_c%poly(1)%ptr)
  deallocate(region_c%poly(1)%ptr)
  deallocate(region_c%poly)
  call free_region_tree(region_c)

  ! free the (temporary) surface region tree
  nullify(region_s%poly(1)%ptr)
  call free_region_tree(region_s)


  ! manage the newly discovered intersection points
  do ipt = ntuvxyz_tmp+1,ntuvxyz ! <--------------------------------+
     ! convert from (t,u,v) to (u1,v1,u2,v2)                        !
     uv(ivar,icurv) = region(icurv)%ptr%uvbox(2*(ivar-1)+ival)      !
     uv(jvar,icurv) = n1p12ab( &                                    !
          tuvxyz(1,ipt), &                                          !
          region(icurv)%ptr%uvbox(2*jvar-1), &                      !
          region(icurv)%ptr%uvbox(2*jvar  ) )                       !
     uv(:,isurf) = tuvxyz(2:3,ipt)                                  !
     !                                                              !
     tmp(1:2) = uv(:,1)                                             !
     tmp(3:4) = uv(:,2)                                             !
     tmp(5:7) = tuvxyz(4:6,ipt)                                     !
     !                                                              !
     ! check unicity                                                !
     if ( nuvxyz > 0 ) then ! <-------+                             !
        call check_unicity( &         !                             !
             tmp(5:7), &              !                             !
             3, &                     !                             !
             uvxyz(5:7,1:nuvxyz), &   !                             !
             nuvxyz, &                !                             !
             EPSxyz, &                !                             !
             jpt )                    !                             !
     else ! --------------------------+                             !
        jpt = 1                       !                             !
     end if ! <-----------------------+                             !
     !                                                              !
     if ( jpt > nuvxyz ) then ! <-----+                             !
        call append_vector( &         !                             !
             tmp(1:7), &              !                             !
             7, &                     !                             !
             uvxyz, &                 !                             !
             nuvxyz )                 !                             !
        jpt = nuvxyz                  !                             !
     end if ! <-----------------------+                             !
     !                                                              !
     call append_n( &                                               !
          iptsbs, &                                                 !
          nptsbs, &                                                 !
          [jpt], &                                                  !
          1, &                                                      !
          unique=.true. )                                           !
     !                                                              !
     IF ( DEBUG ) THEN
        PRINT *,'+1 UVXYZ :',TMP
        PRINT *,'IPTSBS <---',NUVXYZ
     END IF
  end do ! <--------------------------------------------------------+
  
  if ( allocated(sharedpts) ) deallocate(sharedpts)
  if ( allocated(tuvxyz)    ) deallocate(tuvxyz   )


end subroutine intersect_border_surface
