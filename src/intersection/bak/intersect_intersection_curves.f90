subroutine intersect_intersection_curves( &
     interdata, &
     curvpair, &
     stat )
  USE MOD_UTIL
  use mod_math
  use mod_diffgeom
  use mod_tolerances
  use mod_types_intersection
  implicit none
  LOGICAL, PARAMETER :: DEBUG = .true.!( GLOBALDEBUG .AND. .false. )
  type(type_intersection_data), target, intent(inout) :: interdata
  integer,                              intent(in)    :: curvpair(2)
  type(ptr_intersection_curve)                        :: curv(2)
  logical                                             :: may_intersect
  real(kind=fp)                                       :: uvbox(2,2)
  integer                                             :: numsurf(2)
  type(type_matrix)                                   :: polylineuv(2)
  integer                                             :: np(2)
  integer, allocatable                                :: ipls(:,:)
  real(kind=fp), allocatable                          :: lambda(:,:)
  integer                                             :: npts
  type(ptr_surface)                                   :: surf(3)
  real(kind=fp)                                       :: uv3(2,3), xyz(3), uv2(2,2)
  real(kind=fp)                                       :: w, wi, wip1
  real(kind=fp), dimension(6)                         :: lowerb, upperb
  integer                                             :: stat
  integer                                             :: idpt
  integer                                             :: isurf, jsurf, ivar, icurv, ipt, jpt
  INTEGER :: FID = 13, I
  CHARACTER :: STRNUM

  do icurv = 1,2
     curv(icurv)%ptr => interdata%curves(curvpair(icurv))
  end do

  ! first, check if the two curves can intersect, i.e. they share a common incident surface,
  ! and their uv-boxes on this surface overlap
  may_intersect = .false.
  jloop : do jsurf = 1,2 !  <------------------------------------------------------------------------+
     iloop : do isurf = 1,2 ! <-------------------------------------------------------------------+  !
        if ( associated(curv(1)%ptr%surf(isurf)%ptr, curv(2)%ptr%surf(jsurf)%ptr) ) then ! <---+  !  !
           ! the two intersection curve share a common incident surface                        !  !  !
           ! compute the intersection of their bounding boxes in the uv-space of that surface  !  !  !
           do ivar = 1,2 ! <----------------------------------------------------------------+  !  !  !
              uvbox(:,ivar) = [ &                                                           !  !  !  !
                   max(curv(1)%ptr%uvbox(1,ivar,isurf), curv(2)%ptr%uvbox(1,ivar,jsurf)), & !  !  !  !
                   min(curv(1)%ptr%uvbox(2,ivar,isurf), curv(2)%ptr%uvbox(2,ivar,jsurf)) ]  !  !  !  !
              if ( uvbox(2,ivar) - uvbox(1,ivar) < EPSuv ) then ! <------------+            !  !  !  !
                 ! the uv-boxes do not overlap                                 !            !  !  !  !
                 cycle iloop                                                   !            !  !  !  !
              end if ! <-------------------------------------------------------+            !  !  !  !
           end do ! <-----------------------------------------------------------------------+  !  !  !
           ! if we got this far, the two curves may intersect                                  !  !  !
           numsurf = [isurf, jsurf]                                                            !  !  !
           may_intersect = .true.                                                              !  !  !
           exit jloop                                                                          !  !  !
        end if ! <-----------------------------------------------------------------------------+  !  !
     end do iloop ! <-----------------------------------------------------------------------------+  !
  end do jloop ! <-----------------------------------------------------------------------------------+

  if (.not.may_intersect) return ! the curves cannot intersect

  ! find intersection points between the polylines in the uv-space of the common incident surface
  do icurv = 1,2
     np(icurv) = curv(icurv)%ptr%polyline%np
     allocate(polylineuv(icurv)%mat(2,np(icurv)))
     polylineuv(icurv)%mat(1:2,1:np(icurv)) = curv(icurv)%ptr%polyline%uv(1:2,numsurf(icurv),1:np(icurv))
  end do
  npts = 0
  call intersect_2Dpolylines( &
       polylineuv, &
       [1,1], &
       np, &
       npts, &
       ipls, &
       lambda )
  IF ( DEBUG ) THEN
     DO ICURV = 1,2
        WRITE (STRNUM,'(I1)') ICURV
        OPEN(UNIT=FID, FILE='polylines/xy'//strnum//'.dat', ACTION='WRITE')
        WRITE (FID,*) NP(ICURV)
        DO IPT = 1,NP(ICURV)
           WRITE (FID,*) POLYLINEUV(ICURV)%MAT(1:2,IPT)
        END DO
        CLOSE(FID)
     END DO
     OPEN(UNIT=FID, FILE='polylines/result.dat', ACTION='WRITE')
     WRITE (FID,*) NPTS
     DO IPT = 1,NPTS
        WRITE (FID,*) IPLS(:,IPT)
     END DO
     DO IPT = 1,NPTS
        WRITE (FID,*) LAMBDA(:,IPT)
     END DO
     CLOSE(FID)
  END IF
  deallocate(polylineuv(1)%mat, polylineuv(2)%mat)

  if ( npts < 1 ) return ! the polylines do not intersect, we assume the curves do not either

  !! refine all intersection points using Newton-Raphson algorithm and add them to the collections  
  ! form the triplet of surfaces (one common incident surface, plus one other for each curve)
  surf(1)%ptr => curv(1)%ptr%surf(numsurf(1))%ptr
  do icurv = 1,2
     surf(1+icurv)%ptr => curv(icurv)%ptr%surf(1+mod(numsurf(icurv),2))%ptr
  end do
  ! lower and upper bounds of feasible domain
  lowerb(1:2) = uvbox(1,1:2)
  upperb(1:2) = uvbox(2,1:2)
  do icurv = 1,2
     isurf = 1 + mod(numsurf(icurv),2)
     lowerb(2*(icurv+1)-1:2*(icurv+1)) = curv(icurv)%ptr%uvbox(1,1:2,isurf)
     upperb(2*(icurv+1)-1:2*(icurv+1)) = curv(icurv)%ptr%uvbox(2,1:2,isurf)
  end do
  lowerb = lowerb - EPSuv
  upperb = upperb + EPSuv

  do ipt = 1,npts ! <----------------------------------------------------------------------------------+
     PRINT *,'IPLS   =',IPLS(:,IPT)
     PRINT *,'LAMBDA =',LAMBDA(:,IPT)
     PRINT *,'SEGMENTS :'
     DO ICURV = 1,2
        PRINT *,CURV(ICURV)%PTR%POLYLINE%UV(:,numsurf(ICURV),ipls(ICURV,ipt))
        PRINT *,CURV(ICURV)%PTR%POLYLINE%UV(:,numsurf(ICURV),ipls(ICURV,ipt)+1)
     END DO
     PRINT *,'UV ='
     DO ICURV = 1,2
        PRINT *,(1._FP - LAMBDA(ICURV,IPT))*CURV(ICURV)%PTR%POLYLINE%UV(:,numsurf(ICURV),ipls(ICURV,ipt)) + &
             LAMBDA(ICURV,IPT)*CURV(ICURV)%PTR%POLYLINE%UV(:,numsurf(ICURV),ipls(ICURV,ipt)+1)
     END DO
     ! set initial iterate (result of polyline intersection)                                           !
     uv3(:,1) = &                                                                                      !
          (1._fp - lambda(1,ipt)) * curv(1)%ptr%polyline%uv(:,numsurf(1),ipls(1,ipt)) + &              !
          lambda(1,ipt) * curv(1)%ptr%polyline%uv(:,numsurf(1),ipls(1,ipt)+1)                          !
     do icurv = 1,2 ! <-----------------------------------------------------------------------------+  !
        isurf = 1 + mod(numsurf(icurv),2)                                                           !  !
        uv3(:,1+icurv) = &                                                                          !  !
             (1._fp - lambda(icurv,ipt)) * curv(icurv)%ptr%polyline%uv(:,isurf,ipls(icurv,ipt)) + & !  !
             lambda(icurv,ipt) * curv(icurv)%ptr%polyline%uv(:,isurf,ipls(icurv,ipt)+1)             !  !
     end do ! <-------------------------------------------------------------------------------------+  !
     !                                                                                                 !
     ! run Newton-Raphson algorithm                                                                    !
     PRINT *,'SMOOTH? ', curv(1)%ptr%smooth, curv(2)%ptr%smooth
     PRINT *,'UV0 = ',UV3
     if ( curv(2)%ptr%smooth ) then
        call newton_three_surfaces_1tangential( &
             surf, &
             lowerb, &
             upperb, &
             stat, &
             uv3, &
             xyz )
     else
        call newton_three_surfaces( &                                                                     !
             surf, &                                                                                      !
             lowerb, &                                                                                    !
             upperb, &                                                                                    !
             stat, &                                                                                      !
             uv3, &                                                                                       !
             xyz )                                                                                        !
     end if
     !                                                                                                 !
     IF ( DEBUG ) THEN
        IF ( STAT == 0 ) THEN
           PRINT *,'NEWTON 3 SURFACES CONVERGED'
           PRINT *,' UV =',UV3
           PRINT *,'XYZ =',XYZ
        ELSE
           CALL WRITE_POLYNOMIAL( SURF(1)%PTR%X, 'dev_intersection/debug3si_surf1.cheb' )
           CALL WRITE_POLYNOMIAL( SURF(2)%PTR%X, 'dev_intersection/debug3si_surf2.cheb' )
           CALL WRITE_POLYNOMIAL( SURF(3)%PTR%X, 'dev_intersection/debug3si_surf3.cheb' )
           CALL GET_FREE_UNIT(FID)
           OPEN(UNIT=FID, FILE='dev_intersection/debug3si.dat', ACTION='WRITE')
           WRITE(FID,*) LOWERB
           WRITE(FID,*) UPPERB
           DO ICURV = 1,2
              WRITE(FID,*) NUMSURF(ICURV)
              WRITE(FID,*) CURV(ICURV)%PTR%POLYLINE%NP
              DO I = 1,CURV(ICURV)%PTR%POLYLINE%NP
                 WRITE(FID,*) CURV(ICURV)%PTR%POLYLINE%UV(:,:,I)
              END DO
           END DO
           CLOSE(FID)
           STOP '----> DEBUG 3SI'
        END IF
     END IF
     if ( stat > 0 ) then ! <---------------+                                                          !
        ! Newton failed to converge         !                                                          !
        STOP 'NEWTON_3_SURFACES FAILED TO CONVERGE'
        return                              !                                                          !
     end if ! <-----------------------------+                                                          !
     !                                                                                                 !
     ! add the split point the intersection points collection                                          !
     call add_intersection_point( &                                                                    !
          uv3, &                                                                                       !
          xyz, &                                                                                       !
          surf, &                                                                                      !
          3, &                                                                                         !
          interdata, &                                                                                 !
          idpt )                                                                                       !
     !                                                                                                 !
     ! add the split point to each intersected curve's polyline                                        !
     curv_loop : do icurv = 1,2 ! <---------------------------------------------------+                !
        ! compare the new split point to already discovered split points              !                !
        do jpt = 1,curv(icurv)%ptr%nsplit ! <---------------------------------+       !                !
           if ( idpt == curv(icurv)%ptr%isplit(1,jpt) ) cycle curv_loop       !       !                !
        end do ! <------------------------------------------------------------+       !                !
        !                                                                             !                !
        ! locate the segement in which the new split point must be inserted           !                !
        if ( curv(icurv)%ptr%smooth ) then
           ! ...
        else
           w    = dot_product(curv(icurv)%ptr%param_vector, xyz)                         !                !
           wi   = dot_product(curv(icurv)%ptr%param_vector, &                            !                !
                curv(icurv)%ptr%polyline%xyz(:,ipls(icurv,ipt)))                   !                !
           wip1 = dot_product(curv(icurv)%ptr%param_vector, &                            !                !
                curv(icurv)%ptr%polyline%xyz(:,ipls(icurv,ipt)+1))                 !                !
           while_loop : do ! <---------------------------------------------------+       !                !
              if ( ipls(icurv,ipt) < 1 .or. &                                    !       !                !
                   ipls(icurv,ipt) >= curv(icurv)%ptr%polyline%np ) then ! <--+  !       !                !
                 PRINT *,'W =',W
                 PRINT *,&
                      dot_product(curv(icurv)%ptr%param_vector, curv(icurv)%ptr%polyline%xyz(:,1)), &
                      dot_product(curv(icurv)%ptr%param_vector, curv(icurv)%ptr%polyline%xyz(:,curv(icurv)%ptr%polyline%np))
                 STOP 'COULD NOT FIND CORRECT SEGMENT'                        !  !       !                !
              end if ! <------------------------------------------------------+  !       !                !
              if ( w > wip1 ) then ! <-----------------------------------+       !       !                !
                 ! the split point is located downstream                 !       !       !                !
                 ipls(icurv,ipt) = ipls(icurv,ipt) + 1                   !       !       !                !
                 wi = wip1                                               !       !       !                !
                 wip1 = dot_product(curv(icurv)%ptr%param_vector, &      !       !       !                !
                      curv(icurv)%ptr%polyline%xyz(:,ipls(icurv,ipt)+1)) !       !       !                !
              elseif ( w < wi ) then ! ----------------------------------+       !       !                !
                 ! the split point is located upstream                   !       !       !                !
                 ipls(icurv,ipt) = ipls(icurv,ipt) - 1                   !       !       !                !
                 IF ( ipls(icurv,ipt) < 1 ) THEN
                    PRINT *,idpt,curv(icurv)%ptr%isplit(1,1:curv(icurv)%ptr%nsplit)
                    PRINT *,XYZ
                    DO JPT = 1,curv(icurv)%ptr%nsplit
                       PRINT *,INTERDATA%POINTS(curv(icurv)%ptr%isplit(1,JPT))%XYZ
                    END DO
                 END IF
                 wi = dot_product(curv(icurv)%ptr%param_vector, &        !       !       !                !
                      curv(icurv)%ptr%polyline%xyz(:,ipls(icurv,ipt)))   !       !       !                !
                 wip1 = wi                                               !       !       !                !
              else ! ----------------------------------------------------+       !       !                !
                 ! we found the correct segment                          !       !       !                !
                 exit while_loop                                         !       !       !                !
              end if ! <-------------------------------------------------+       !       !                !
           end do while_loop ! <-------------------------------------------------+       !                !
        end if
        !                                                                             !                !
        ! insert a new polyline point                                                 !                !
        uv2(:,numsurf(icurv)) = uv3(:,1)                                              !                !
        uv2(:,1+mod(numsurf(icurv),2)) = uv3(:,1+icurv)                               !                !
        call insert_polyline_point( &                                                 !                !
             uv2, &                                                                   !                !
             xyz, &                                                                   !                !
             stat, &                                                                  !                !
             curv(icurv)%ptr%polyline, &                                              !                !
             i=ipls(icurv,ipt) )                                                      !                !
        !                                                                             !                !
        if ( stat > 0 ) then ! <----------------------+                               !                !
           ! the polyine point could no be inserted   !                               !                !
           stat = stat + 2                            !                               !                !
           return                                     !                               !                !
        end if ! <------------------------------------+                               !                !
        !                                                                             !                !
        ! locate in the curve's split points collection where the new point           !                !
        ! must be inserted (split points are sorted in ascending order of w)          !                !
        ! at the same time, update indices of split points located downstream         !                !
        do jpt = curv(icurv)%ptr%nsplit,1,-1 ! <-----------------------------------+  !                !
           if ( curv(icurv)%ptr%isplit(2,jpt) > ipls(icurv,ipt) ) then ! <------+  !  !                !
              curv(icurv)%ptr%isplit(2,jpt) = curv(icurv)%ptr%isplit(2,jpt) + 1 !  !  !                !
           else ! --------------------------------------------------------------+  !  !                !
              exit                                                              !  !  !                !
           end if ! <-----------------------------------------------------------+  !  !                !
        end do ! <-----------------------------------------------------------------+  !                !
        !                                                                             !                !
        ! finally, add the split point the curve's split points list                  !                !
        call insert_column_after( &                                                   !                !
             curv(icurv)%ptr%isplit, &                                                !                !
             2, &                                                                     !                !
             curv(icurv)%ptr%nsplit, &                                                !                !
             [idpt, ipls(icurv,ipt)+1], &                                             !                !
             jpt )                                                                    !                !
     end do curv_loop ! <-------------------------------------------------------------+                !
     !                                                                                                 !
  end do ! <-------------------------------------------------------------------------------------------+

end subroutine intersect_intersection_curves
