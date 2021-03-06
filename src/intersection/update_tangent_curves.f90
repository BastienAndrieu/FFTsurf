subroutine update_intersection_curves( &
     interdata )
  use mod_types_intersection
  use mod_diffgeom
  use mod_tolerances
  implicit none
  type(type_intersection_data), intent(inout), target :: interdata
  type(type_point_on_surface), pointer                :: pos => null()
  type(type_intersection_curve), pointer              :: curve => null()
  real(kind=fp)                                       :: uv(2,2), xyz(3,2), err
  real(kind=fp), dimension(4)                         :: lowerb, upperb
  integer                                             :: stat
  integer                                             :: ipoint, ipos, icurv, i, j
 
  upperb(:) = 1._fp + EPSuv
  lowerb = -upperb

  do ipoint = 1,interdata%np
     interdata%points(ipoint)%xyz(1:3) = 0._fp
     pos => interdata%points(ipoint)%pos
     do ipos = 1,interdata%points(ipoint)%npos
        call eval( &
             xyz(1:3,1), &
             pos%surf, &
             pos%uv )
        interdata%points(ipoint)%xyz = interdata%points(ipoint)%xyz + xyz(1:3,1)
        pos => pos%next
     end do
     interdata%points(ipoint)%xyz = interdata%points(ipoint)%xyz / real(interdata%points(ipoint)%npos, kind=fp)
  end do
  
  do icurv = 1,interdata%nc ! <---------------------------------+
     !PRINT *,'ICURV =',ICURV
     curve => interdata%curves(icurv)                           !
     if ( .not.curve%smooth ) cycle                             !
     !                                                          !
     do i = 1,curve%polyline%np ! <--------------------------+  !
        uv = curve%polyline%uv(1:2,1:2,i)                    !  !
        do j = 1,2 ! <----------------+                      !  !
           call eval( &               !                      !  !
                xyz(:,j), &           !                      !  !
                curve%surf(j)%ptr, &  !                      !  !
                uv(:,j) )             !                      !  !
        end do ! <--------------------+                      !  !
        curve%polyline%xyz(1:3,i) = 0.5_fp * sum(xyz, 2)     !  !
        !                                                    !  !
        err = sum((xyz(:,1) - xyz(:,2))**2)                  !  !
        if ( err > EPSxyzsqr ) then ! <-----------------+    !  !
           PRINT *,'CURVE #',ICURV
           PRINT *,I,sqrt(err)
           call simultaneous_point_inversions( &        !    !  !
                curve%surf, &                           !    !  !
                lowerb, &                               !    !  !
                upperb, &                               !    !  !
                stat, &                                 !    !  !
                uv, &                                   !    !  !
                xyz(:,1) )                              !    !  !
           if ( stat > 0 ) then ! <------------------+  !    !  !
              ! keep current coordinates             !  !    !  !
              PRINT *,'update_intersection_curves: *** /!\ FAILED TO RECOVER POINT ON TANGENT INTERSECTION, ERR=',SQRT(ERR)
           else ! -----------------------------------+  !    !  !
              curve%polyline%uv(1:2,1:2,i) = uv      !  !    !  !
              PRINT *,I,NORM2(XYZ(:,1) - CURVE%POLYLINE%XYZ(:,I))
              PRINT *,I,uv
              curve%polyline%xyz(1:3,i) = xyz(:,1)   !  !    !  !
           end if ! <--------------------------------+  !    !  !
        end if ! <--------------------------------------+    !  !
        !                                                    !  !
     end do ! <----------------------------------------------+  !
  end do ! <----------------------------------------------------+
  
end subroutine update_intersection_curves
