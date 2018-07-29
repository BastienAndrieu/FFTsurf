recursive subroutine intersect_simple_surfaces( &
     surfroot, &
     region, &
     param_vector, &
     interdata, &
     uvxyz, &
     nuvxyz, &
     stat_degeneracy )
  use mod_util
  use mod_math
  use mod_polynomial
  use mod_bernstein
  use mod_diffgeom
  use mod_regiontree
  use mod_types_intersection
  implicit none
  LOGICAL, PARAMETER :: DEBUG = ( GLOBALDEBUG .AND. .true. )
  type(ptr_surface),            intent(in)    :: surfroot(2)
  type(ptr_region),             intent(inout) :: region(2)
  real(kind=fp),                intent(in)    :: param_vector(3)
  integer,                      intent(inout) :: stat_degeneracy
  type(type_intersection_data), intent(inout) :: interdata
  real(kind=fp), allocatable,   intent(inout) :: uvxyz(:,:)
  integer,                      intent(inout) :: nuvxyz
  integer, allocatable                        :: ipts_ss(:)
  integer                                     :: npts_ss
  type(type_polynomial)                       :: regc
  type(type_curve)                            :: root_c
  integer, allocatable                        :: ipts_bs(:)
  integer                                     :: npts_bs
  real(kind=fp)                               :: uv_subdiv(2,2)
  integer                                     :: stat_subdiv(2), nchild(2)
  type(ptr_region)                            :: newregion(2)
  integer                                     :: stat_point(2)
  integer                                     :: order(2)
  integer                                     :: icurv, ivar, ival, ipt, ichild, jchild
  INTEGER :: ISURF

  if ( stat_degeneracy /= 0 ) return ! a degeneracy has been encountered

  IF ( DEBUG ) THEN
     !PAUSE
     PRINT *,''; PRINT *,'';
     PRINT *,'INTERSECT_SIMPLE_SURFACES'
     PRINT *,'UVBOXES ='
     DO ISURF = 1,2 ! <-----------------+
        PRINT *,REGION(ISURF)%PTR%UVBOX !
     END DO ! <-------------------------+
  END IF

  ! inherit from parent regions all intersection points contained in current regions
  if ( nuvxyz > 0 ) then ! <----------------------------+
     do isurf = 1,2 ! <---------------------------+     !
        call inherit_points( &                    !     !
             region(isurf)%ptr, &                 !     !
             uvxyz(2*isurf-1:2*isurf,1:nuvxyz), & !     !
             uvxyz(8,1:nuvxyz), &                 !     !
             nuvxyz )                             !     !
     end do ! <-----------------------------------+     !
  end if ! <--------------------------------------------+

  IF ( DEBUG ) THEN
     PRINT *,'IPTS ='
     DO ISURF = 1,2
        IF ( REGION(ISURF)%PTR%NPTS < 1 ) THEN
           PRINT *,'N/A'
        ELSE
           PRINT *,REGION(ISURF)%PTR%IPTS(1:REGION(ISURF)%PTR%NPTS)
        END IF
     END DO
  END IF

  ! get the list of already discovered intersection points contained in both surface regions
  npts_ss = 0
  if ( region(1)%ptr%npts > 0 .and. region(2)%ptr%npts > 0 ) then ! <-----------------------+
     call intersection_arrays( &                                                            !
          region(1)%ptr%ipts(1:region(1)%ptr%npts), &                                       !
          region(2)%ptr%ipts(1:region(2)%ptr%npts), &                                       !
          ipts_ss )                                                                         !
     if ( allocated(ipts_ss) ) npts_ss = size(ipts_ss)                                      !
  end if ! <--------------------------------------------------------------------------------+

  ! intersect the 4*2 pairs of border-surface
  npts_bs = 0
  outer : do icurv = 1,2 ! <------------------------------------------------+
     ! compute the parameterization of the surface region                   !
     ! of which we consider the border (Chebyshev polynomial basis)         !
     call chgvar2( &                                                        !
          surfroot(icurv)%ptr%x, &                                          !
          regc, &                                                           !
          region(icurv)%ptr%uvbox([1,3]), &                                 !
          region(icurv)%ptr%uvbox([2,4]) )                                  !
     !                                                                      !
     do ivar = 1,2 ! <---------------------------------------------------+  !
        !                                                                !  !
        do ival = 1,2 ! <---------------------------------------------+  !  !
           ! compute the parameterization of the border curve         !  !  !
           ! (still expressed in Chebyshev polynomial basis)          !  !  !
           call bivar2univar( &                                       !  !  !
                regc, &                                               !  !  !
                root_c%x, &                                           !  !  !
                ivar, &                                               !  !  !
                ival )                                                !  !  !
           ! perform degree reduction (without loss of accuracy)      !  !  !
           call economize1(root_c%x, EPSmath)                         !  !  !
           ! compute first and second derivatives of the curve        !  !  !
           call compute_deriv1(root_c)                                !  !  !
           call compute_deriv2(root_c)                                !  !  !
           !                                                          !  !  !
           ! compute border-surface intersection                      !  !  !
           call intersect_border_surface( &                           !  !  !
                surfroot, &                                           !  !  !
                root_c, &                                             !  !  !
                region, &                                             !  !  !
                icurv, &                                              !  !  !
                ivar, &                                               !  !  !
                ival, &                                               !  !  !
                ipts_ss, &                                            !  !  !
                npts_ss, &                                            !  !  !
                uvxyz, &                                              !  !  !
                nuvxyz, &                                             !  !  !
                ipts_bs, &                                            !  !  !
                npts_bs, &                                            !  !  !
                stat_degeneracy )                                     !  !  !
           !                                                          !  !  !
           if ( stat_degeneracy > 1 ) then ! <-------------------+    !  !  !
              ! a degeneracy has been encountered                !    !  !  !
              exit outer                                         !    !  !  !
           end if ! <--------------------------------------------+    !  !  !
           !                                                          !  !  !
        end do ! <----------------------------------------------------+  !  !
     end do ! <----------------------------------------------------------+  !
  end do outer ! <----------------------------------------------------------+
  if ( allocated(ipts_ss) ) deallocate(ipts_ss)

  ! free allocated polynomials
  call free_polynomial(regc      )
  call free_polynomial(root_c%x  )
  call free_polynomial(root_c%xt )
  call free_polynomial(root_c%xtt)

  ! append the new border-surface intersection points to each region's collection
  if ( npts_bs > 0 ) then ! <-----------------------------+
     do isurf = 1,2 ! <------------------------+          !
        call add_points_bottom_up( &           !          !
             region(isurf)%ptr, &              !          !
             ipts_bs(1:npts_bs), &             !          !
             npts_bs )                         !          !
     end do ! <--------------------------------+          !
  end if ! <----------------------------------------------+

  if ( stat_degeneracy > 1 ) then ! <---------------------+
     ! a degeneracy has been encountered                  !
     if ( allocated(ipts_bs) ) deallocate(ipts_bs)        !
     return                                               !
  end if ! <----------------------------------------------+

  IF ( DEBUG ) THEN
     PRINT *,'NPTS_BS =',NPTS_BS
     IF ( NPTS_BS > 0 ) THEN
        PRINT *,'IPTS_BS =',IPTS_BS
        CALL PRINT_MAT( TRANSPOSE(UVXYZ(1:7,IPTS_BS(1:NPTS_BS))) )
     END IF
  END IF

  if ( npts_bs == 0 ) then
     ! no border-surface intersection points
     if ( allocated(ipts_bs) ) deallocate(ipts_bs)
     return
  end if

  if ( npts_bs <= 2 ) then ! <------------------------------------------------+
     ! classify the border-surface intersection points (entering, exiting,    !
     ! isolated)                                                              !
     IF ( DEBUG ) THEN ! <---------------------------------------+            !
        PRINT *,'BSI POINTS: UV ='                               !            !
        CALL PRINT_MAT(TRANSPOSE(UVXYZ(1:4,IPTS_BS(1:NPTS_BS)))) !            !
     END IF ! <--------------------------------------------------+            !
     call classify_border_surface_intersection_point( &                       !
          surfroot, &                                                         !
          region, &                                                           !
          reshape( uvxyz(1:4,ipts_bs(1:npts_bs)), [2,2,npts_bs] ), &          !
          npts_bs, &                                                          !
          stat_point(1:npts_bs) )                                             !
     !                                                                        !
     !                                                                        !
     if ( npts_bs == 2 ) then ! <-------------------------------------+       !
        if ( stat_point(1)*stat_point(2) < 0 ) then ! <---------+     !       !
           ! 1 entering point and 1 exiting point.              !     !       !
           ! Re-order the points from entering to exiting       !     !       !
           if ( stat_point(1) < 0 ) then ! <---------------+    !     !       !
              order = [2,1]                                !    !     !       !
           else ! -----------------------------------------+    !     !       !
              order = [1,2]                                !    !     !       !
           end if ! <--------------------------------------+    !     !       !
           !                                                    !     !       !
           ! trace curve without ambiguity                      !     !       !
           call add_intersection_curve( &                       !     !       !
                interdata, &                                    !     !       !
                param_vector, &                                 !     !       !
                ipts_bs(order), &                               !     !       !
                reshape([region(1)%ptr%uvbox,&                  !     !       !
                region(2)%ptr%uvbox], [2,2,2]) )                !     !       !
           IF ( DEBUG ) PRINT *,'+1 INTERSECTION CURVE'
           !                                                    !     !       !
        elseif ( all(stat_point == 0) ) then ! -----------------+     !       !
           ! 2 isolated points                                  !     !       !
        else ! -------------------------------------------------+     !       !
           ! incorrect configuration                            !     !       !
           PRINT *,'------------------------------------------'
           PRINT *,'2 BSI POINTS - INCORRECT CONFIGURATION'
           PRINT *,'UVBOXES ='
           DO ISURF = 1,2 ! <-----------------+
              PRINT *,REGION(ISURF)%PTR%UVBOX !
           END DO ! <-------------------------+
           PRINT *,'STAT_POINT =',STAT_POINT
           PRINT *,'UVXYZ ='
           CALL PRINT_MAT( TRANSPOSE(UVXYZ(1:7,IPTS_BS(1:2))) )
           PRINT *,'TOLUV =',UVXYZ(8,IPTS_BS(1:2))
           CALL WRITE_POLYNOMIAL( SURFROOT(1)%PTR%X, 'dev_intersection/debugssi_surf1.cheb' )
           CALL WRITE_POLYNOMIAL( SURFROOT(2)%PTR%X, 'dev_intersection/debugssi_surf2.cheb' )
           CALL WRITE_POLYNOMIAL( REGION(1)%PTR%POLY(1)%PTR, 'dev_intersection/debugbsi_reg1.bern' )
           CALL WRITE_POLYNOMIAL( REGION(2)%PTR%POLY(1)%PTR, 'dev_intersection/debugbsi_reg2.bern' )
           PRINT *,'------------------------------------------'
           stat_degeneracy = 35                                 !     !       !
        end if ! <----------------------------------------------+     !       !
     end if ! <-------------------------------------------------------+       !
     !                                                                        !
     !                                                                        !
  elseif ( npts_bs > 2 ) then ! ----------------------------------------------+
     ! More than 2 border-surface intersection points have been found,        !
     ! the situation is ambiguous so we need to carry on the recursion.       !
     ! first, try to subdivide at an bsi point interior to at least one of    !
     ! the intersected surface regions.                                       !
     do ipt = 1,npts_bs ! <----------------------------------------------+    !
        do isurf = 1,2 ! <--------------------------------------------+  !    !
           uv_subdiv(:,isurf) = uvxyz(2*isurf-1:2*isurf,ipts_bs(ipt)) !  !    !
           call subdiv_region( &                                      !  !    !
                region(isurf)%ptr, &                                  !  !    !
                uv_subdiv(:,isurf), &                                 !  !    !
                stat_subdiv(isurf) )                                  !  !    !
        end do ! <----------------------------------------------------+  !    !
        if ( any(stat_subdiv /= 2) ) exit                                !    !
     end do ! <----------------------------------------------------------+    !
     !                                                                        !
     if ( all(stat_subdiv == 2) ) then ! <-------------------------------+    !
        ! all bsi points are at corners in both regions, subdivide at    !    !
        ! parametric midpoint                                            !    !
        do isurf = 1,2 ! <--------------------------------------------+  !    !
           uv_subdiv(:,isurf) = 0.5_fp * ( &                          !  !    !
                region(isurf)%ptr%uvbox([1,3]) + &                    !  !    !
                region(isurf)%ptr%uvbox([2,4]) )                      !  !    !
           call subdiv_region( &                                      !  !    !
                region(isurf)%ptr, &                                  !  !    !
                uv_subdiv(:,isurf), &                                 !  !    !
                stat_subdiv(isurf) )                                  !  !    !
        end do ! <----------------------------------------------------+  !    !
     end if ! <----------------------------------------------------------+    !
     !
     if ( all(stat_subdiv == 2) ) then ! <-------------------------------+    !
        ! error, at least one region should be subdivided                !    !
        stat_degeneracy = 37                                             !    !
        return                                                           !    !
     end if ! <----------------------------------------------------------+    !
     !                                                                        !
     ! compute Bezier control points for children regions                     !
     do isurf = 1,2 ! <--------------------------------------------------+    !
        if ( stat_subdiv(isurf) == 2 ) then ! <-----------------------+  !    !
           ! the subdivision point is at one of the surface's corners !  !    !
           nchild(isurf) = 1                                          !  !    !
           cycle                                                      !  !    !
        end if ! <----------------------------------------------------+  !    !
        ! the subdivision point is interior to the surface               !    !
        nchild(isurf) = size(region(isurf)%ptr%child)                    !    !
        if ( stat_subdiv(isurf) < 0 ) then ! <-----------+               !    !
           ! the region's children are already allocated !               !    !
           cycle                                         !               !    !
        end if ! <---------------------------------------+               !    !
        !                                                                !    !
        ! map uv_subdiv to the region's local frame ([0,1]^2)            !    !
        do ivar = 1,2 ! <--------------------------------+               !    !
           uv_subdiv(ivar,isurf) = 0.5_fp * ( 1._fp + &  !               !    !
                ab2n1p1( &                               !               !    !
                uv_subdiv(ivar,isurf), &                 !               !    !
                region(isurf)%ptr%uvbox(2*ivar-1), &     !               !    !
                region(isurf)%ptr%uvbox(2*ivar) ) )      !               !    !
        end do ! <---------------------------------------+               !    !
        !                                                                !    !
        ! allocate polynomials                                           !    !
        do ichild = 1,nchild(isurf) ! <--------------------------+       !    !
           allocate(region(isurf)%ptr%child(ichild)%poly(1)    ) !       !    !
           allocate(region(isurf)%ptr%child(ichild)%poly(1)%ptr) !       !    !
        end do ! <-----------------------------------------------+       !    !
        !                                                                !    !
        if ( stat_subdiv(isurf) == 0 ) then ! <----------------------+   !    !
           ! 4 children                                              !   !    !
           call subdiv_bezier2( &                                    !   !    !
                region(isurf)%ptr%poly(1)%ptr, &                     !   !    !
                uv_subdiv(:,isurf), &                                !   !    !
                bsw=region(isurf)%ptr%child(1)%poly(1)%ptr, &        !   !    !
                bse=region(isurf)%ptr%child(2)%poly(1)%ptr, &        !   !    !
                bnw=region(isurf)%ptr%child(3)%poly(1)%ptr, &        !   !    !
                bne=region(isurf)%ptr%child(4)%poly(1)%ptr )         !   !    !
        elseif ( stat_subdiv(isurf) == 1 ) then ! -------------------+   !    !
           ! 2 children                                              !   !    !
           if ( region(isurf)%ptr%child(2)%uvbox(1) <= &             !   !    !
                region(isurf)%ptr%uvbox(1) + EPSregion ) then ! <--+ !   !    !
              ! the subdivision point is on a iso-u boundary       ! !   !    !
              call subdiv_bezier2_only_v( &                        ! !   !    !
                   region(isurf)%ptr%poly(1)%ptr, &                ! !   !    !
                   v=uv_subdiv(2,isurf), &                         ! !   !    !
                   bs=region(isurf)%ptr%child(1)%poly(1)%ptr, &    ! !   !    !
                   bn=region(isurf)%ptr%child(2)%poly(1)%ptr )     ! !   !    !
           else ! -------------------------------------------------+ !   !    !
              ! the subdivision point is on a iso-v boundary       ! !   !    !
              call subdiv_bezier2_only_u( &                        ! !   !    !
                   region(isurf)%ptr%poly(1)%ptr, &                ! !   !    !
                   u=uv_subdiv(1,isurf), &                         ! !   !    !
                   bw=region(isurf)%ptr%child(1)%poly(1)%ptr, &    ! !   !    !
                   be=region(isurf)%ptr%child(2)%poly(1)%ptr )     ! !   !    !
           end if ! <----------------------------------------------+ !   !    !
        end if ! <---------------------------------------------------+   !    !
     end do ! -----------------------------------------------------------+    !
     !                                                                        !
     ! carry on the recursion with pairs of children regions                  !
     do jchild = 1,nchild(2) ! <--------------------------------------+       !
        if ( nchild(2) == 1 ) then ! <-------------------------+      !       !
           newregion(2)%ptr => region(2)%ptr                   !      !       !
        else ! ------------------------------------------------+      !       !
           newregion(2)%ptr => region(2)%ptr%child(jchild)     !      !       !
        end if ! <---------------------------------------------+      !       !
        !                                                             !       !
        do ichild = 1,nchild(1) ! <-------------------------------+   !       !
           if ( nchild(1) == 1 ) then ! <----------------------+  !   !       !
              newregion(1)%ptr => region(1)%ptr                !  !   !       !
           else ! ---------------------------------------------+  !   !       !
              newregion(1)%ptr => region(1)%ptr%child(ichild)  !  !   !       !
           end if ! <------------------------------------------+  !   !       !
           !                                                      !   !       !
           call intersect_simple_surfaces( &                      !   !       !
                surfroot, &                                       !   !       !
                newregion, &                                      !   !       !
                param_vector, &                                   !   !       !
                interdata, &                                      !   !       !
                uvxyz, &                                          !   !       !
                nuvxyz, &                                         !   !       !
                stat_degeneracy )                                 !   !       !
           !                                                      !   !       !
        end do ! <------------------------------------------------+   !       !
     end do ! <-------------------------------------------------------+       !
     !                                                                        !
  end if ! <------------------------------------------------------------------+

  if ( allocated(ipts_bs) ) deallocate(ipts_bs)
  
end subroutine intersect_simple_surfaces
