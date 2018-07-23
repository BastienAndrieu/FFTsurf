subroutine intersect_border_surface( &
     root_s, &
     root_c, &
     region, &
     icurv, &
     ivar, &
     ival, &
     ipts_ss, &
     npts_ss, &
     uvxyz, &
     nuvxyz, &
     ipts_bs, &
     npts_bs, &
     stat_degeneracy )
  use mod_util
  use mod_math
  use mod_polynomial
  use mod_diffgeom
  use mod_regiontree
  use mod_types_intersection
  use mod_tolerances
  implicit none
  LOGICAL, PARAMETER :: DEBUG = ( GLOBALDEBUG .AND. .true. )
  type(ptr_surface),          intent(in)    :: root_s(2)
  type(type_curve),           intent(in)    :: root_c
  type(ptr_region),           intent(inout) :: region(2)
  integer,                    intent(in)    :: icurv, ivar, ival
  integer, allocatable,       intent(in)    :: ipts_ss(:)
  integer,                    intent(in)    :: npts_ss
  real(kind=fp), allocatable, intent(inout) :: uvxyz(:,:)
  integer,                    intent(inout) :: nuvxyz
  integer, allocatable,       intent(inout) :: ipts_bs(:)
  integer,                    intent(inout) :: npts_bs
  integer,                    intent(inout) :: stat_degeneracy
  integer                                   :: isurf, jvar
  type(type_region)                         :: region_c
  type(type_region)                         :: region_s
  real(kind=fp), allocatable                :: tuvxyz(:,:)
  integer                                   :: ntuvxyz, ntuvxyz_tmp
  real(kind=fp)                             :: tmp(7), uv(2,2)
  integer                                   :: ipt, jpt

  IF ( DEBUG ) THEN
     PRINT *,''; PRINT *,'';
     PRINT *,'INTERSECT_BORDER_SURFACE'
     PRINT *,'ICURV, IVAR, IVAL =',ICURV,IVAR,IVAL
     PRINT *,'UVBOXES ='
     DO ISURF = 1,2 ! <-----------------+
        PRINT *,REGION(ISURF)%PTR%UVBOX !
     END DO ! <-------------------------+
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
  IF ( DEBUG ) THEN
     PRINT *,'NPTS_SS =',npts_ss
     IF ( NPTS_SS > 0 ) THEN
        PRINT *,ipts_ss(1:npts_ss)
        CALL PRINT_MAT( TRANSPOSE(UVXYZ(1:4,ipts_ss(1:npts_ss))) )
     END IF
  END IF

  ! check if some of these points are located on the current border
  do jpt = 1,npts_ss ! <--------------------------------------------------------------------+
     ipt = ipts_ss(jpt)                                                                     !
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
             ipts_bs, &                                                                 !   !
             npts_bs, &                                                                 !   !
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
  if ( stat_degeneracy == 50 ) stat_degeneracy = 0
  IF ( DEBUG ) THEN
     PRINT *,'NTUVXYZ =',NTUVXYZ
     IF ( NTUVXYZ > 0 ) CALL PRINT_MAT( TRANSPOSE(TUVXYZ(:,1:NTUVXYZ)) )
     !CALL EXPORT_REGION_TREE( REGION_C, 'dev_intersection/treebsi_c.dat' )
     !CALL EXPORT_REGION_TREE( REGION_S, 'dev_intersection/treebsi_s.dat' )
  END IF

  IF ( stat_degeneracy > 0 ) THEN
     CALL WRITE_POLYNOMIAL( ROOT_C%X, 'dev_intersection/debugbsi_c.cheb' )
     CALL WRITE_POLYNOMIAL( ROOT_S(isurf)%ptr%X, 'dev_intersection/debugbsi_s.cheb' )
     CALL WRITE_POLYNOMIAL( REGION_C%POLY(1)%PTR, 'dev_intersection/debugbsi_c.bern' )
     CALL WRITE_POLYNOMIAL( REGION_S%POLY(1)%PTR, 'dev_intersection/debugbsi_s.bern' )
     CALL EXPORT_REGION_TREE( REGION_C, 'dev_intersection/treebsi_c.dat' )
     CALL EXPORT_REGION_TREE( REGION_S, 'dev_intersection/treebsi_s.dat' )
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
          ipts_bs, &                                                !
          npts_bs, &                                                !
          [jpt], &                                                  !
          1, &                                                      !
          unique=.true. )                                           !
     !                                                              !
     IF ( DEBUG ) THEN
        PRINT *,'+1 UVXYZ :',TMP
        PRINT *,'IPTS_BS <---',NUVXYZ
     END IF
  end do ! <--------------------------------------------------------+
  
  if ( allocated(tuvxyz)    ) deallocate(tuvxyz   )


end subroutine intersect_border_surface
