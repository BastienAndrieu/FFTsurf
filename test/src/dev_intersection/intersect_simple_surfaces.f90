recursive subroutine intersect_simple_surfaces( &
     surfroot, &
     region, &
     param_vector, &
     interdata, &
     uvxyz, &
     nuvxyz, &
     stat_degeneracy )
  use mod_math
  use mod_polynomial
  use mod_bernstein2
  use mod_diffgeom2
  use mod_regiontree
  use mod_types_intersection
  implicit none
  LOGICAL, PARAMETER :: DEBUG = ( GLOBALDEBUG .AND. .TRUE. )
  type(ptr_surface),            intent(in)    :: surfroot(2)
  type(ptr_region),             intent(inout) :: region(2)
  real(kind=fp),                intent(in)    :: param_vector(3)
  integer,                      intent(inout) :: stat_degeneracy
  type(type_intersection_data), intent(inout) :: interdata
  real(kind=fp), allocatable,   intent(inout) :: uvxyz(:,:)
  integer,                      intent(inout) :: nuvxyz
  type(type_polynomial)                       :: regc
  type(type_curve)                            :: root_c
  integer, allocatable                        :: iptsbs(:)
  integer                                     :: nptsbs
  real(kind=fp)                               :: uv_subdiv(2)
  integer                                     :: stat_subdiv
  type(ptr_region)                            :: newregion(2)
  integer                                     :: stat_point(2)
  integer                                     :: order(2)
  integer                                     :: icurv, ivar, ival, ichild, jchild

  if ( stat_degeneracy /= 0 ) return ! a degeneracy has been encountered

  IF ( DEBUG ) THEN
     PRINT *,''; PRINT *,'';
     PRINT *,'INTERSECT_SIMPLE_SURFACES'
     PRINT *,'UVBOXES ='
     DO ISURF = 1,2 ! <-----------------+
        PRINT *,REGION(ISURF)%PTR%UVBOX !
     END DO ! <-------------------------+
     PRINT *,'IPTS ='
     DO ISURF = 1,2
        IF ( REGION(ISURF)%PTR%NPTS < 1 ) THEN
           PRINT *,'N/A'
        ELSE
           PRINT *,REGION(ISURF)%PTR%IPTS(1:REGION(ISURF)%PTR%NPTS)
        END IF
     END DO
  END IF

  ! inherit from parent regions all intersection points contained in current regions
  if ( nuvxyz > 0 ) then ! <-----------------+
     do isurf = 1,2 ! <----------------+     !
        call inherit_points( &         !     !
             region(isurf)%ptr, &      !     !
             uvxyz, &                  !     !
             nuvxyz )                  !     !
     end do ! <------------------------+     !
  end if ! <---------------------------------+


  ! intersect the 4*2 pairs of border-surface
  nptsbs = 0
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
                uvxyz, &                                              !  !  !
                nuvxyz, &                                             !  !  !
                iptsbs, &                                             !  !  !
                nptsbs, &                                             !  !  !
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

  ! free allocated polynomials
  call free_polynomial(regc      )
  call free_polynomial(root_c%x  )
  call free_polynomial(root_c%xt )
  call free_polynomial(root_c%xtt)
  
  if ( nptsbs > 0 ) then ! <------------------------------+
     do isurf = 1,2 ! <------------------------+          !
        call add_points_bottom_up( &           !          !
             region(isurf)%ptr, &              !          !
             iptsbs(1:nptsbs), &               !          !
             nptsbs )                          !          !
     end do ! <--------------------------------+          !
  end if ! <----------------------------------------------+


  if ( stat_degeneracy > 1 ) then ! <---------------------+
     if ( allocated(iptsbs) ) deallocate(iptsbs)          !
     return                                               !
  end if ! <----------------------------------------------+

  IF ( DEBUG ) THEN
     PRINT *,'NPTSBS =',NPTSBS
     IF ( NPTSBS > 0 ) THEN
        PRINT *,'IPTSBS =',IPTSBS
        CALL PRINT_MAT( TRANSPOSE(UVXYZ(:,IPTSBS(1:NPTSBS))) )
     END IF

     PRINT *,'REGION%IPTS ='
     DO ISURF = 1,2
        IF ( REGION(ISURF)%PTR%NPTS < 1 ) THEN
           PRINT *,'N/A'
        ELSE
           PRINT *,REGION(ISURF)%PTR%IPTS(1:REGION(ISURF)%PTR%NPTS)
        END IF
     END DO

  END IF

  if ( nptsbs > 2 ) then ! <--------------------------------------------------+
     ! More than 2 border-surface intersection points have been found,        !
     ! the situation is ambiguous so we need to carry on the recursion.       !
     ! Both regions are subdivided at their parametric center point           !
     do isurf = 1,2 ! <----------------------------------------------------+  !
        uv_subdiv = 0.5_fp * ( &                                           !  !
             region(isurf)%ptr%uvbox([1,3]) + &                            !  !
             region(isurf)%ptr%uvbox([2,4]) )                              !  !
        !                                                                  !  !
        call subdiv_region( &                                              !  !
             region(isurf)%ptr, &                                          !  !
             uv_subdiv, &                                                  !  !
             stat_subdiv )                                                 !  !
        !                                                                  !  !
        if ( stat_subdiv > 0 ) then ! <---------------------------------+  !  !
           ! the region is subdivided into less than 4 children,        !  !  !
           ! this should not happen                                     !  !  !
           stat_degeneracy = 33                                         !  !  !
           return                                                       !  !  !
        elseif ( stat_subdiv < 0 ) then ! ------------------------------+  !  !
           ! the region already has children                            !  !  !
           if ( size(region(isurf)%ptr%child) < 4 ) then ! <---+        !  !  !
              ! the region has less than 4 children,           !        !  !  !
              ! this should not happen                         !        !  !  !
              stat_degeneracy = 34                             !        !  !  !
           end if ! <------------------------------------------+        !  !  !
        elseif ( stat_subdiv == 0 ) then ! -----------------------------+  !  !
           ! the region does not have children yet                      !  !  !
           do ichild = 1,4 ! <----------------------------------------+ !  !  !
              allocate(region(isurf)%ptr%child(ichild)%poly(1)    )   ! !  !  !
              allocate(region(isurf)%ptr%child(ichild)%poly(1)%ptr)   ! !  !  !
           end do ! <-------------------------------------------------+ !  !  !
           call subdiv_bezier2( &                                       !  !  !
                region(isurf)%ptr%poly(1)%ptr, &                        !  !  !
                [0.5_fp, 0.5_fp], &                                     !  !  !
                bsw=region(isurf)%ptr%child(1)%poly(1)%ptr, &           !  !  !
                bse=region(isurf)%ptr%child(2)%poly(1)%ptr, &           !  !  !
                bnw=region(isurf)%ptr%child(3)%poly(1)%ptr, &           !  !  !
                bne=region(isurf)%ptr%child(4)%poly(1)%ptr )            !  !  !
        end if ! <------------------------------------------------------+  !  !
     end do ! <------------------------------------------------------------+  !
     !                                                                        !
     ! carry on the recursion with the 4*4 new pairs of regions               !
     do jchild = 1,4 ! <---------------------------------------------------+  !
        newregion(2)%ptr => region(2)%ptr%child(jchild)                    !  !
        do ichild = 1,4 ! <---------------------------------------------+  !  !
           newregion(1)%ptr => region(1)%ptr%child(ichild)              !  !  !
           call intersect_simple_surfaces( &                            !  !  !
                surfroot, &                                             !  !  !
                newregion, &                                            !  !  !
                param_vector, &                                         !  !  !
                interdata, &                                            !  !  !
                uvxyz, &                                                !  !  !
                nuvxyz, &                                               !  !  !
                stat_degeneracy )                                       !  !  !
        end do ! <------------------------------------------------------+  !  !
     end do ! <------------------------------------------------------------+  !
     !                                                                        !
  elseif ( nptsbs > 0 ) then ! -----------------------------------------------+
     ! 0 < npts <= 2                                                          !
     ! classify the border-surface intersection points (entering, exiting,    !
     ! isolated)                                                              !
     call classify_border_surface_intersection_point( &                       !
          surfroot, &                                                         !
          region, &                                                           !
          reshape( uvxyz(1:4,iptsbs(1:nptsbs)), [2,2,nptsbs] ), &             !
          nuvxyz, &                                                           !
          stat_point(1:nptsbs) )                                              !
     !                                                                        !
     if ( nptsbs == 2 ) then ! <--------------------------------------+       !
        IF ( DEBUG ) PRINT *,'HERE'
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
                iptsbs(order), &                                !     !       !
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
           PRINT *,'STAT_POINT =',STAT_POINT
           PRINT *,'UVXYZ ='
           CALL PRINT_MAT( TRANSPOSE(UVXYZ(:,IPTSBS(1:2))) )
           PRINT *,'------------------------------------------'
           stat_degeneracy = 35                                 !     !       !
        end if ! <----------------------------------------------+     !       !
        !                                                             !       !
     elseif ( nptsbs == 1 ) then ! -----------------------------------+       !
        if ( stat_point(1) == 0 ) then ! <----------------------+     !       !
           ! 1 isolated point                                   !     !       !
        else ! -------------------------------------------------+     !       !
           ! incorrect configuration                            !     !       !
           PRINT *,'------------------------------------------'
           PRINT *,'1 BSI POINT - INCORRECT CONFIGURATION'
           PRINT *,'STAT_POINT =',STAT_POINT
           PRINT *,'UVXYZ ='
           PRINT *,UVXYZ(:,IPTSBS(1))
           PRINT *,'------------------------------------------'
           stat_degeneracy = 36                                 !     !       !
        end if ! <----------------------------------------------+     !       !
        !                                                             !       !
     end if ! <-------------------------------------------------------+       !
     !                                                                        !
  end if ! <------------------------------------------------------------------+

  IF ( DEBUG ) PRINT *,'FREE IPTSBS...'
  if ( allocated(iptsbs) ) deallocate(iptsbs)
  IF ( DEBUG ) PRINT *,'           ...OK'

end subroutine intersect_simple_surfaces
