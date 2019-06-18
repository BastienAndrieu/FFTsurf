subroutine remesh_intersection_polyline( &
     surf, &
     lowerb, &
     upperb, &
     polyline, &
     head, &
     tail, &
     tolchord, &
     hmin, &
     hmax, &
     stat )
  use mod_math
  use mod_diffgeom
  use mod_types_intersection
  implicit none
  type(ptr_surface),                intent(in)    :: surf(2)
  real(kind=fp),                    intent(in)    :: lowerb(4)
  real(kind=fp),                    intent(in)    :: upperb(4)
  type(type_intersection_polyline), intent(inout) :: polyline
  integer,                          intent(in)    :: head
  integer,                          intent(inout) :: tail
  real(kind=fp),                    intent(in)    :: tolchord
  real(kind=fp),                    intent(in)    :: hmin
  real(kind=fp),                    intent(in)    :: hmax
  integer,                          intent(out)   :: stat
  real(kind=fp)                                   :: FRACcurvature_radius
  integer                                         :: stat_tangent
  real(kind=fp)                                   :: duv_ds(2,2,2), dxyz_ds(3,2)
  real(kind=fp)                                   :: curvature(2)
  real(kind=fp)                                   :: h(2)
  real(kind=fp)                                   :: hmintol, hmaxtol
  real(kind=fp)                                   :: delta_s, lambda
  real(kind=fp)                                   :: uv(2,2), xyz(3)
  integer                                         :: stat_relax, stat_insert, stat_remove
  integer                                         :: ipt

  FRACcurvature_radius = 2._fp*sqrt(tolchord*(2._fp - tolchord))

  stat = 0
  ipt = head

  call diffgeom_intersection( &
       surf, &
       polyline%uv(1:2,1:2,ipt), &
       duv_ds, &
       dxyz_ds, &
       stat_tangent, &
       curvature )
  h(1) = FRACcurvature_radius/max(EPSfp, curvature(1)) ! /!\ if smooth curve
  h(1) = min(hmax, max(hmin, h(1)))


  outer : do 
     PRINT *,'   IPT =',IPT,', HEAD =', HEAD,', TAIL =',TAIL
     if ( ipt == tail ) exit outer
     
     call diffgeom_intersection( &
          surf, &
          polyline%uv(1:2,1:2,ipt+1), &
          duv_ds, &
          dxyz_ds, &
          stat_tangent, &
          curvature )
     h(2) = FRACcurvature_radius/max(EPSfp, curvature(1)) ! /!\ if smooth curve
     h(2) = min(hmax, max(hmin, h(2)))

     hmintol = (0.25_fp*sum(h))**2
     hmaxtol = sum(h)**2
     PRINT *,'   HMINTOL, HMAXTOL =', SQRT(HMINTOL), SQRT(HMAXTOL)

     delta_s = sum((polyline%xyz(1:3,ipt+1) - polyline%xyz(1:3,ipt))**2)
     PRINT *,'   DELTA_S =', SQRT(DELTA_S)
     if ( delta_s > hmaxtol ) then ! <-------------------------------------------------------+
        PRINT *,'--> TOO LONG'
        ! segment too long --> split it                                                      !
        ! get weighted midpoint                                                              !
        lambda = 0.5_fp!h(1)/sum(h)                                                                 !
        ! initial guess                                                                      !
        uv = (1._fp - lambda)*polyline%uv(1:2,1:2,ipt) + lambda*polyline%uv(1:2,1:2,ipt+1)   !
        ! relax to exact intersection                                                        !
        call simultaneous_point_inversions( &                                                !
             surf, &                                                                         !
             lowerb, &                                                                       !
             upperb, &                                                                       !
             stat_relax, &                                                                   !
             uv, &                                                                           !
             xyz )                                                                           !
        if ( stat_relax > 0 ) then ! <---------------------------------------------------+   !
           ! error                                                                       !   !
           PRINT *,'!!! remesh_intersection_polyline: FAILED TO RELAX POINT'             !   !
        else ! --------------------------------------------------------------------------+   !
           ! ok :)                                                                       !   !
           call insert_polyline_point( &                                                 !   !
                uv, &                                                                    !   !
                xyz, &                                                                   !   !
                stat_insert, &                                                           !   !
                polyline, &                                                              !   !
                ipt )                                                                    !   !
           if ( stat_insert == 0 ) then ! <------------------+                           !   !
              tail = tail + 1                                !                           !   !
              cycle outer                                    !                           !   !
           else ! -------------------------------------------+                           !   !
              ! error                                        !                           !   !
              PRINT *,'!!! remesh_intersection_polyline: FAILED TO INSERT POLYLINE POINT'!   !
           end if ! <----------------------------------------+                           !   !
        end if ! <-----------------------------------------------------------------------+   !
        !                                                                                    !
     elseif ( delta_s < hmintol ) then ! <---------------------------------------------------+
        PRINT *,'--> TOO SHORT'
        ! too short --> collapse it                                                          !
        if ( ipt+1 == tail) then ! <-------------------------+                               !
           ! remove point #ipt                               !                               !
           call remove_polyline_point( &                     !                               !
                polyline, &                                  !                               !
                ipt, &                                       !                               !
                stat_remove )                                !                               !
        else ! ----------------------------------------------+                               !
           ! remove point #(ipt+1)                           !                               !
           call remove_polyline_point( &                     !                               !
                polyline, &                                  !                               !
                ipt+1, &                                     !                               !
                stat_remove )                                !                               !
        end if ! <-------------------------------------------+                               !
        !                                                                                    !
        if ( stat_remove == 0 ) then ! <---------------------+                               !
           tail = tail - 1                                   !                               !
           cycle outer                                       !                               !
        else ! ----------------------------------------------+                               !
           ! error                                           !                               !
           PRINT *,'!!! remesh_intersection_polyline: FAILED TO REMOVE POLYLINE POINT'       !
        end if ! <-------------------------------------------+                               !
        !                                                                                    !
     end if ! <------------------------------------------------------------------------------+

     ipt = ipt + 1

  end do outer

end subroutine remesh_intersection_polyline
