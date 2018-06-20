! *** INCOMPLET ***
subroutine tangent_to_intersection_curve( &
     surf, &
     uv, &
     uv_s, &
     xyz_s, &
     stat )
  use mod_math
  use mod_diffgeom2
  ! ( see "Shape interrogation for computer aided design and manufacturing", &
  ! Patrikalakis et al. (2009), pp.166-175)
  implicit none
  type(ptr_surface), intent(in)  :: surf(2)
  real(kind=fp),     intent(in)  :: uv(2,2)
  real(kind=fp),     intent(out) :: uv_s(2,2,2)
  real(kind=fp),     intent(out) :: xyz_s(3,2)
  integer,           intent(out) :: stat
  real(kind=fp)                  :: xyz_uv(3,2,2)
  real(kind=fp)                  :: n(3,2), nsqr(2), aux(3), normxyz_s
  integer                        :: isurf, ivar, i

  ! compute normal vectors to both surfaces
  do isurf = 1,2 ! <--------------------------------------------------+
     do ivar = 1,2 ! <-------------------------+                      !
        call evald1( &                         !                      !
             xyz_uv(:,ivar,isurf), &           !                      !
             surf(isurf)%ptr, &                !                      !
             uv(:,isurf), &                    !                      !
             ivar )                            !                      !
     end do ! <--------------------------------+                      !
     n(:,isurf) = cross( xyz_uv(:,1,isurf), xyz_uv(:,2,isurf) )       !
  end do ! <----------------------------------------------------------+

  xyz_s(:,1) = cross( n(:,1), n(:,2) )
  normxyz_s = norm2( xyz_s(:,1) )
  if ( normxyz_s > EPSmath ) then ! <------------------------------------------------+
     ! the normals are not collinear, this is a transversal intersection point       !
     stat = 0                                                                        !
     xyz_s(:,1) = xyz_s(:,1) / norm2( xyz_s(:,1) )                                   !
  else ! ----------------------------------------------------------------------------+
     ! this is a tangential intersection point, the tangent direction(s) to          !
     ! the intersection banch(es) is (are) computed using 2nd derivatives            !
     call characterize_tangential_intersection_point( &                              !
          surf, &                                                                    !
          uv, &                                                                      !
          xyz_uv, &                                                                  !
          n, &                                                                       !
          xyz_s, &                                                                   !
          stat )                                                                     !
  end if ! <-------------------------------------------------------------------------+

  if ( stat < 2 ) then ! <---------------------------------------------------------------------+
     ! compute tangential direction(s) in the uv-space of each surface                         !
     nsqr = sum( n**2, 1 )                                                                     !
     do isurf = 1,2 ! loop over surfaces ! <---------------------------------------------+     !
        do i = 1,stat+1 ! loop over tangential directions ! <-------------------------+  !     ! 
           aux = cross( xyz_s(:,i), n(:,isurf) )                                      !  !     !
           do ivar = 1,2 ! loop over parameters u,v ! <----------------------------+  !  !     !
              uv_s(ivar,i,isurf) = real( (-1)**ivar, kind=fp) * &                  !  !  !     !
                   dot_product( xyz_uv(:,1+mod(ivar,2),isurf), aux ) / nsqr(isurf) !  !  !     !
           end do ! <--------------------------------------------------------------+  !  !     !
        end do ! <--------------------------------------------------------------------+  !     !
     end do ! <--------------------------------------------------------------------------+     !
  else ! --------------------------------------------------------------------------------------+
     ! the tangential direction is undefined                                                   !
     ! ...                                                                                     !
  end if! <------------------------------------------------------------------------------------+

end subroutine tangent_to_intersection_curve
