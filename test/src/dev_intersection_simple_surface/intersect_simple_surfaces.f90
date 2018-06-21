recursive subroutine intersect_simple_surfaces( &
     surfroot, &
     region, &
     param_vector, &
     interdat, &
     !listcurv, &
     uvxyz, &
     nuvxyz, &
     stat_degeneracy )
  use mod_math
  use mod_diffgeom2
  use mod_regiontree
  implicit none
  integer, parameter                          :: nuvxyz_init = 10
  type(ptr_surface),            intent(in)    :: surfroot(2)
  type(ptr_region),             intent(inout) :: region(2)
  real(kind=fp),                intent(in)    :: param_vector(3)
  integer,                      intent(inout) :: stat_degeneracy
  type(type_intersection_data), intent(inout) :: interdat
  !type(type_listcurves),        intent(inout) :: listcurv
  real(kind=fp), allocatable,   intent(inout) :: uvxyz(:,:)
  integer,                      intent(inout) :: nuvxyz
  type(type_polynomial)                       :: regc
  type(type_curve)                            :: root_c
  integer                                     :: stat_point(2)
  real(kind=fp)                               :: uv_subdiv(2)
  integer                                     :: stat_subdiv
  type(ptr_region)                            :: newregion(2)
  logical, allocatable                        :: mask(:,:)
  real(kind=fp), allocatable                  :: newuvxyz(:,:)
  integer                                     :: nnewuvxyz
  integer                                     :: idnewpoint(2), order(2)
  real(kind=fp)                               :: uv_endpoints(2,2,2), xyz_endpoints(3,2)
  integer                                     :: icurv, ivar, ival, ipt
  integer                                     :: ichild, jchild
  CHARACTER                                   :: STR1

  if ( stat_degeneracy > 1 ) return

  !PRINT *,'';PRINT *,'';PRINT *,'';
  !PRINT *,'UVBOXES ='
  !DO ISURF = 1,2
  !   PRINT *,REGION(ISURF)%PTR%UVBOX
  !END DO

  ! intersect the 4 borders of each surface with the other surface
  !nuvxyz = 0
  !allocate( uvxyz(7,nuvxyz_init) )
  if ( .not.allocated(uvxyz) ) allocate( uvxyz(7,nuvxyz_init) )
  !PRINT *,'BEFORE, ',NUVXYZ,' POINTS'
  do icurv = 1,2 ! <--------------------------------------------------------+
     call chgvar2( &                                                        !
       surfroot(icurv)%ptr%x, &                                             !
       regc, &                                                              !
       region(icurv)%ptr%uvbox([1,3]), &                                    !
       region(icurv)%ptr%uvbox([2,4]) )                                     !
     !                                                                      !
     do ivar = 1,2 ! <---------------------------------------------------+  !
        do ival = 1,2 ! <---------------------------------------------+  !  !
           !                                                          !  !  !
           call bivar2univar( &                                       !  !  !
                regc, &                                               !  !  !
                root_c%x, &                                           !  !  !
                ivar, &                                               !  !  !
                ival )                                                !  !  !
           call economize1( root_c%x, EPSmath )                       !  !  !
           call compute_deriv1( root_c )                              !  !  !
           call compute_deriv2( root_c )                              !  !  !
           !                                                          !  !  !
           call intersect_border_surface( &                           !  !  !
                surfroot, &                                           !  !  !
                root_c, &                                             !  !  !  
                region, &                                             !  !  !
                icurv, &                                              !  !  !
                ivar, &                                               !  !  !
                ival, &                                               !  !  !
                uvxyz, &                                              !  !  !
                nuvxyz, &                                             !  !  !
                stat_degeneracy )                                     !  !  !
           !                                                          !  !  !
           if ( stat_degeneracy > 1 ) then ! <-------------------+    !  !  !
              ! propagation erreur vers routine appelante ...    !    !  !  !
              PRINT *,'STAT_DEGENERACY =',STAT_DEGENERACY
              return                                             !    !  !  !
           end if ! <--------------------------------------------+    !  !  !
           !                                                          !  !  !
        end do ! <----------------------------------------------------+  !  !
     end do ! <----------------------------------------------------------+  !
  end do ! <----------------------------------------------------------------+
  call free_polynomial( regc )
  call free_polynomial( root_c%x )
  call free_polynomial( root_c%xt )
  call free_polynomial( root_c%xtt )
  !PRINT *,'AFTER, ',NUVXYZ,' POINTS'


  if ( nuvxyz > 0 ) then ! <---------------------------------------------------------------+
     !PRINT *,NUVXYZ,' POINTS:'
     !CALL PRINT_MAT( TRANSPOSE( UVXYZ(:,1:NUVXYZ) ) )
     if ( nuvxyz > 2 ) then ! <------------------------------------------------------+     !
        ! subdivide both regions at their center point                               !     !
        do isurf = 1,2 ! <-------------------------------------+                     !     !
           !                                                   !                     !     !
           uv_subdiv = 0.5_fp * ( &                            !                     !     !
                region(isurf)%ptr%uvbox([1,3]) + &             !                     !     !
                region(isurf)%ptr%uvbox([2,4]) )               !                     !     !
           !                                                   !                     !     !
           call subdiv_region( &                               !                     !     !
                region(isurf)%ptr, &                           !                     !     !
                uv_subdiv, &                                   !                     !     !
                stat_subdiv )                                  !                     !     !
           !                                                   !                     !     !
           if ( stat_subdiv > 0 ) then ! <-----+               !                     !     !
              stat_degeneracy = 33             !               !                     !     !
              PRINT *,'STAT_DEGENERACY =',STAT_DEGENERACY
              return                           !               !                     !     !
           end if ! <--------------------------+               !                     !     !
           !                                                   !                     !     !
           if ( stat_subdiv == 0 ) then ! ----------------+    !                     !     !
              allocate( &                                 !    !                     !     !
                   region(isurf)%ptr%child(1)%poly(1), &     !    !                     !     !
                   region(isurf)%ptr%child(2)%poly(1), &     !    !                     !     !
                   region(isurf)%ptr%child(3)%poly(1), &     !    !                     !     !
                   region(isurf)%ptr%child(4)%poly(1) )      !    !                     !     !
              allocate( &                                 !    !                     !     !
                   region(isurf)%ptr%child(1)%poly(1)%ptr, &     !    !                     !     !
                   region(isurf)%ptr%child(2)%poly(1)%ptr, &     !    !                     !     !
                   region(isurf)%ptr%child(3)%poly(1)%ptr, &     !    !                     !     !
                   region(isurf)%ptr%child(4)%poly(1)%ptr )      !    !                     !     !
              call subdiv_bezier2( &                      !    !                     !     ! 
                   region(isurf)%ptr%poly(1)%ptr, &              !    !                     !     !
                   [0.5_fp, 0.5_fp], &                    !    !                     !     !
                   bsw=region(isurf)%ptr%child(1)%poly(1)%ptr, & !    !                     !     !
                   bse=region(isurf)%ptr%child(2)%poly(1)%ptr, & !    !                     !     !
                   bnw=region(isurf)%ptr%child(3)%poly(1)%ptr, & !    !                     !     !
                   bne=region(isurf)%ptr%child(4)%poly(1)%ptr )  !    !                     !     !
           end if ! <-------------------------------------+    !                     !     !
           !                                                   !                     !     !
        end do ! <---------------------------------------------+                     !     !
        !                                                                            !     !
        ! carry on the recursion with the 4*4 new pairs of regions                   !     !
        allocate( mask(nuvxyz,2) )                                                   !     !
        mask(:,:) = .true.                                                           !     !
        do jchild = 1,size(region(2)%ptr%child) ! <-----------------------------+    !     !
           newregion(2)%ptr => region(2)%ptr%child(jchild)                      !    !     ! 
           !                                                                    !    !     !     
           do ipt = 1,nuvxyz ! <---------------------------------------------+  !    !     !
              if ( .not.is_in_interval( &                                    !  !    !     !
                   uvxyz(3,ipt), &                                           !  !    !     !
                   newregion(2)%ptr%uvbox(1), &                              !  !    !     !
                   newregion(2)%ptr%uvbox(2) ) .or. &                        !  !    !     !
                   .not.is_in_interval( &                                    !  !    !     !
                   uvxyz(4,ipt), &                                           !  !    !     !
                   newregion(2)%ptr%uvbox(3), &                              !  !    !     !
                   newregion(2)%ptr%uvbox(4) ) ) mask(ipt,2) = .false.       !  !    !     !
           end do ! <--------------------------------------------------------+  !    !     !
           !PRINT *,'JCHILD =',JCHILD,', UVBOX=',NEWREGION(2)%PTR%UVBOX
           !PRINT *,'MASK =',MASK(:,2)
           !                                                                    !    !     !     
           do ichild = 1,size(region(1)%ptr%child) ! <-----------------------+  !    !     !    
              newregion(1)%ptr => region(1)%ptr%child(ichild)                !  !    !     !
              !                                                              !  !    !     !
              do ipt = 1,nuvxyz ! <---------------------------------------+  !  !    !     !
                 if ( .not.is_in_interval( &                              !  !  !    !     !
                      uvxyz(1,ipt), &                                     !  !  !    !     !
                      newregion(1)%ptr%uvbox(1), &                        !  !  !    !     !
                      newregion(1)%ptr%uvbox(2) ) .or. &                  !  !  !    !     !
                      .not.is_in_interval( &                              !  !  !    !     !
                      uvxyz(2,ipt), &                                     !  !  !    !     !
                      newregion(1)%ptr%uvbox(3), &                        !  !  !    !     !
                      newregion(1)%ptr%uvbox(4) ) ) mask(ipt,1) = .false. !  !  !    !     !
              end do ! <--------------------------------------------------+  !  !    !     !
              !PRINT *,'ICHILD =',ICHILD,', UVBOX=',NEWREGION(1)%PTR%UVBOX
              !PRINT *,'MASK =',MASK(:,1)
              !                                                              !  !    !     !
              nnewuvxyz = count( mask(:,1) .and. mask(:,2) )                 !  !    !     !
              !PRINT *,'NNEWUVXYZ =',NNEWUVXYZ
              if ( nnewuvxyz > 0 ) then ! <-------------------+              !  !    !     !
                 allocate( newuvxyz(7,nnewuvxyz) )            !              !  !    !     !
                 newuvxyz(1:7,1:nnewuvxyz) = reshape( &       !              !  !    !     !
                      pack( &                                 !              !  !    !     !
                      uvxyz(1:7,1:nuvxyz), &                  !              !  !    !     !
                      spread( mask(:,1).and.mask(:,2), &      !              !  !    !     !
                      dim=1, ncopies=7 ) ), [7,nnewuvxyz] )   !              !  !    !     !
                 !PRINT *,'NEWUVXYZ ='
                 !CALL PRINT_MAT( TRANSPOSE( NEWUVXYZ(:,1:NNEWUVXYZ) ) )
              end if ! <--------------------------------------+              !  !    !     !
              !                                                              !  !    !     !
              call intersect_simple_surfaces( &                              !  !    !     !
                   surfroot, &                                               !  !    !     !
                   newregion, &                                              !  !    !     !
                   param_vector, &                                           !  !    !     !
                   interdat, &
                   !listcurv, &                                               !  !    !     !
                   newuvxyz, &                                               !  !    !     !
                   nnewuvxyz, &                                              !  !    !     !
                   stat_degeneracy )                                         !  !    !     !
              if ( allocated(newuvxyz) ) deallocate( newuvxyz )              !  !    !     !
              !                                                              !  !    !     !  
           end do ! <--------------------------------------------------------+  !    !     !
        end do ! <--------------------------------------------------------------+    !     !
        deallocate( mask )                                                           !     ! 
        !                                                                            !     !
     else ! -------------------------------------------------------------------------+     !
        ! (nuvxyz <= 2)                                                              !     !
        ! classifier points : entering, exiting, isolated                            !     !
        call classify_border_surface_intersection_point( &                           !     !
             surfroot, &                                                             !     !
             region, &                                                               !     !
             reshape( uvxyz(1:4,1:nuvxyz), [2,2,nuvxyz] ), &                         !     !
             nuvxyz, &                                                               !     !
             stat_point(1:nuvxyz) )                                                  !     !
        !                                                                            !     !
        if ( nuvxyz == 2 ) then ! <----------------------------+                     !     !
           if ( product(stat_point) < 0 ) then ! <----------+  !                     !     !
              ! trace courbe...                             !  !                     !     !
              PRINT *,'1 CURVE :)' 
              CALL PRINT_MAT( TRANSPOSE(UVXYZ(:,1:2)) )

              !listcurv%nc = listcurv%nc + 1                    !** à supprimer
              !listcurv%curve(listcurv%nc)%uvxyz = uvxyz(:,1:2) !** à supprimer

              interdat%nc = interdat%nc + 1              
              do isurf = 1,2
                 interdat%curves(interdat%nc)%surf(isurf)%ptr => surfroot(isurf)%ptr
                 !interdat%curves(interdat%nc)%region(isurf)%ptr => region(isurf)%ptr
                 interdat%curves(interdat%nc)%uvbox(:,:,isurf) = reshape( region(isurf)%ptr%uvbox, [2,2] )
              end do

              ! add the two new intersection points to the intersection data structure
              do ipt = 1,2
                 call add_intersection_point( &
                      reshape( uvxyz(1:4,ipt), [2,2] ), &
                      uvxyz(5:7,ipt), &
                      surfroot, &
                      interdat, &
                      idnewpoint(ipt) )
              end do
              do isurf = 1,2
                 call append_n( &
                      region(isurf)%ptr%ipts, &
                      region(isurf)%ptr%npts, &
                      idnewpoint(1:2), &
                      2, &
                      unique=.true. )
              end do

              ! re-order the endpoints from entering to exiting
              if ( stat_point(1) < 0 ) then
                 order = [2,1]
              else
                 order = [1,2]
              end if
              interdat%curves(interdat%nc)%root%endpoints = idnewpoint(order)
              uv_endpoints = reshape( uvxyz(1:4,order), [2,2,2] )
              xyz_endpoints = uvxyz(5:7,order)

              interdat%curves(interdat%nc)%param_vector = param_vector
              !INTERDAT%CURVES(INTERDAT%NC)%PARAM_VECTOR = UVXYZ(5:7,2) - UVXYZ(5:7,1)

              ! ****************************************
              IF ( .FALSE. ) THEN
                 DO ISURF = 1,2
                    WRITE (STR1,'(I1)') ISURF
                    CALL WRITE_POLYNOMIAL( SURFROOT(ISURF)%PTR%X, 'trace_intersection_polyline/c' // STR1 // '.cheb' )
                 END DO
                 OPEN(UNIT=13, FILE='trace_intersection_polyline/data.dat', ACTION='WRITE')
                 WRITE (13,*) 'UV_ENDPOINTS'
                 DO IPT = 1,2
                    DO ISURF = 1,2
                       WRITE (13,*) uv_endpoints(:,isurf,ipt)
                    END DO
                 END DO
                 WRITE (13,*) 'XYZ_ENDPOINTS'
                 DO IPT = 1,2
                    WRITE (13,*) xyz_endpoints(:,ipt)
                 END DO
                 WRITE (13,*) 'PARAM_VECTOR'
                 WRITE (13,*) interdat%curves(interdat%nc)%param_vector
                 CLOSE(13)
              END IF
              ! ****************************************

              ALLOCATE( &
                   INTERDAT%CURVES(INTERDAT%NC)%POLYLINE%UV(2,2,100), &
                   INTERDAT%CURVES(INTERDAT%NC)%POLYLINE%XYZ(3,100) )
              call trace_intersection_polyline( &
                   surfroot, &
                   uv_endpoints, &
                   xyz_endpoints, &
                   interdat%curves(interdat%nc)%param_vector, &
                   interdat%curves(interdat%nc)%polyline, &
                   HMIN=REAL(1.E-3,KIND=FP), &
                   HMAX=REAL(2.E-1,KIND=FP) )  

           elseif ( all(stat_point == 0) ) then ! ----------+  !                     !     !
              ! 2 points isolés...                          !  !                     !     !
              PRINT *,'2 ISOLATED POINTS ...'
              CALL PRINT_MAT( TRANSPOSE(UVXYZ(:,1:2)) )
              STOP
           else ! ------------------------------------------+  !                     !     !
              ! erreur                                      !  !                     !     !
              stat_degeneracy = 555                         !  !                     !     !
              PRINT *,'***********'
              PRINT *,'2 POINTS, INCORRECT CONFIGURATION ='
              CALL PRINT_MAT( TRANSPOSE(UVXYZ(:,1:2)) )
              PRINT *,'STAT_POINT =',STAT_POINT(1:2)
              DO ISURF = 1,2
                 WRITE (STR1,'(I1)') ISURF
                 CALL WRITE_POLYNOMIAL( &
                      REGION(ISURF)%PTR%POLY(1)%ptr, &
                      'dev_intersection_simple_surface/region_' // STR1 // '.bern' )
              END DO
              PRINT *,'***********'
              return                                        !  !                     !     !
           end if ! <---------------------------------------+  !                     !     !
        else ! ------------------------------------------------+                     !     !
           ! nuvxyz = 1                                        !                     !     !
           if ( stat_point(1) == 0 ) then ! <---------------+  !                     !     !
              ! 1 point isolé...                            !  !                     !     !
              PRINT *,'1 ISOLATED POINT ...'
              PRINT *, UVXYZ(:,1)
              STOP
           else ! ------------------------------------------+  !                     !     !
              ! erreur                                      !  !                     !     !
              stat_degeneracy = 666                         !  !                     !     !
              PRINT *,'***********'
              PRINT *,'1 POINT, INCORRECT CONFIGURATION ='
              PRINT *, UVXYZ(:,1)
              PRINT *,'STAT_POINT =',STAT_POINT(1)
              PRINT *,'***********'
              return                                        !  !                     !     !
           end if ! <---------------------------------------+  !                     !     !
        end if ! <---------------------------------------------+                     !     ! 
        !                                                                            !     !
     end if ! <----------------------------------------------------------------------+     !
  end if ! <-------------------------------------------------------------------------------+

  if ( allocated(uvxyz) ) deallocate( uvxyz )

end subroutine intersect_simple_surfaces
