! ********* INCOMPLET *********
subroutine classify_border_surface_intersection_point( &
     surf, &
     region, &
     uv, &
     npts, &
     stat_point ) 
  use mod_math
  use mod_diffgeom2
  use mod_regiontree
  !stat_point = -1 : exiting
  !              0 : isolated
  !              1 : entering
  implicit none
  type(ptr_surface), intent(in)  :: surf(2)
  type(ptr_region),  intent(in)  :: region(2)
  integer,           intent(in)  :: npts
  real(kind=fp),     intent(in)  :: uv(2,2,npts)
  integer,           intent(out) :: stat_point(npts)
  real(kind=fp)                  :: duv_ds(2,2,2), dxyz_ds(3,2)
  integer                        :: stat_contactpoint
  integer                        :: stat(2), sumstat
  integer                        :: isurf, ipt

  do ipt = 1,npts ! <---------------------------------------------+
     call diffgeom_intersection_curve( &                          !
          surf, &                                                 !
          uv(:,:,ipt), &                                          !
          duv_ds, &                                               !
          dxyz_ds, &                                              !
          stat_contactpoint )                                     !
     !                                                            !
     if ( stat_contactpoint > 1 ) then ! <---------------------+  !
        PRINT *,'classify_border_surface_intersection_point : &
             &stat_contactpoint = ',stat_contactpoint          !  !
        STOP                                                   !  !
        if ( stat_contactpoint == 2 ) then ! <-------------+   !  !
           ! two distinct tangential directions            !   !  !
           ! (...)                                         !   !  !
        else ! --------------------------------------------+   !  !                                               
           ! isolated point                                !   !  !
           ! (...)                                         !   !  !
        end if ! <-----------------------------------------+   !  !
     end if ! <------------------------------------------------+  !
     !                                                            !
     do isurf = 1,2                                               !
        call classify_endpoint( &                                 !
             uv(:,isurf,ipt), &                                   !
             duv_ds(:,1,isurf), &                                 !
             region(isurf)%ptr%uvbox, &                           !
             stat(isurf) )                                        !
     end do                                                       !
     !                                                            !
     sumstat = sum(stat)                                          !
     if ( sumstat == 0 ) then ! <---------+                       !
        stat_point(ipt) = 0               !                       !
     elseif ( sumstat < 0 ) then ! -------+                       !
        stat_point(ipt) = -1              !                       !
     else ! ------------------------------+                       !
        stat_point(ipt) = 1               !                       !
     end if ! <---------------------------+                       !
  end do ! <------------------------------------------------------+

end subroutine classify_border_surface_intersection_point
