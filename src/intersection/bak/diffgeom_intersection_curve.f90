subroutine diffgeom_intersection_curve( &
    surf, &
    uv, &
    duv_ds, &
    dxyz_ds, &
    stat, &
    curvature )
  ! Computes tangent direction(s) and curvature of the intersection curve between
  ! two surfaces at a given intersection point.
  ! Returns: stat =-1 if one surface is singular at that point (normal = zero vector)
  !               = 0 if the point is on a single transveral intersection curve
  !               = 1 if the point is on a single tangential intersection curve
  !               = 2 if the point is a branch point (junction of two intersection curves)
  !               = 3 if the point is an isolated tangential contact point (degenerate curve)
  !               = 4 if the point is an high-order tangential contact point (possible degeneracy)
  use mod_math
  use mod_diffgeom
  ! ( see "Shape interrogation for computer aided design and manufacturing", &
  ! Patrikalakis et al. (2009), pp.166-175, 
  ! and "Tracing surface intersections with validated ODE system solver", Mukundan et al (2004)
  implicit none
  type(ptr_surface),       intent(in)  :: surf(2)
  real(kind=fp),           intent(in)  :: uv(2,2)
  real(kind=fp),           intent(out) :: duv_ds(2,2,2) ! u/v, #branch, #surf
  real(kind=fp),           intent(out) :: dxyz_ds(3,2)
  integer,                 intent(out) :: stat
  real(kind=fp), optional, intent(out) :: curvature(2)
  real(kind=fp)                        :: dxyz_duv(3,2,2), n(3,2), d2xyz_duv2(3)
  real(kind=fp)                        :: cosn, nsqr(2), kn(2), denom
  real(kind=fp), dimension(3,2)        :: EFG, LMN
  integer                              :: ndir
  logical                              :: singular
  integer                              :: isurf, ivar, jvar, i

  ! compute tangent and normal vectors to both surfaces at the intersection point
  do isurf = 1,2 ! <--------------------------------------------------+
     ! tangent vectors                                                !
     do ivar = 1,2 ! <-------------------------+                      !
        call evald1( &                         !                      !
             dxyz_duv(:,ivar,isurf), &         !                      !
             surf(isurf)%ptr, &                !                      !
             uv(:,isurf), &                    !                      !
             ivar )                            !                      !
     end do ! <--------------------------------+                      !
     ! (pseudo)normal vector                                          !
     n(:,isurf) = cross( dxyz_duv(:,1,isurf), dxyz_duv(:,2,isurf) )   !
     ! unit normal vector                                             !
     nsqr(isurf) = sum(n(:,isurf)**2)                                 !
     if ( nsqr(isurf) < 3._fp*EPSfpsqr ) then ! <----+                !
        ! singular surface                           !                !
        stat = -1                                    !                !
        return                                       !                !
     else ! -----------------------------------------+                !
        n(:,isurf) = n(:,isurf) / sqrt(nsqr(isurf))  !                !
     end if ! <--------------------------------------+                !
  end do ! <----------------------------------------------------------+

  ! cosine of the angle between the two normal vectors
  cosn = dot_product( n(:,1), n(:,2) )

  if ( abs(cosn) < 1._fp - EPSmath ) then ! <--------------------------------------------+
     ! the normals are not parallel, this is a transversal intersection point            !
     stat = 0                                                                            !
     ndir = 1 ! one single tangential direction                                          !
     dxyz_ds(:,1) = cross( n(:,1), n(:,2) )                                              !
     dxyz_ds(:,1) = dxyz_ds(:,1) / norm2( dxyz_ds(:,1) )                                 !
  else ! <-------------------------------------------------------------------------------+
     ! this is a tangential intersection point, the tangent direction(s) to              !
     ! the intersection banch(es) is (are) computed using 2nd derivatives                !
     call characterize_tangential_intersection_point( &                                  !
          surf, &                                                                        !
          uv, &                                                                          !
          dxyz_duv, &                                                                    !
          n, &                                                                           !
          stat, &                                                                        !
          dxyz_ds )                                                                      !
     ndir = stat ! 1 or 2 tangential direction(s)                                        !
  end if ! <-----------------------------------------------------------------------------+

  if ( stat > 2 ) then ! <-----------------------+
     ! isolated or high-order contact point      !
     ! => the tangential direction is undefined  !
     return                                      !
  end if ! <-------------------------------------+
  
  ! (at this stage, the number of intersection branches is equal to ndir)

  ! compute tangential direction(s) in the uv-space of each surface
  do isurf = 1,2 ! loop over surfaces ! <------------------------------------------------+
     ! 1st fundamental form coefficients                                                 !
     do ivar = 1,2 ! <-----------------------------------------------+                   !
        do jvar = ivar,2 ! <--------------------------------------+  !                   !
           EFG(ivar + jvar - 1,isurf) = dot_product( &            !  !                   !
                dxyz_duv(:,ivar,isurf), dxyz_duv(:,jvar,isurf) )  !  !                   !
        end do ! <------------------------------------------------+  !                   !
     end do ! <------------------------------------------------------+                   !
     !                                                                                   !
     do i = 1,ndir ! loop over tangential directions ! <-------------+                   !
        call solve_2x2( &                                            !                   !
             duv_ds(:,i,isurf), &                                    !                   !
             EFG(1,isurf), &                                         !                   !
             EFG(2,isurf), &                                         !                   !
             EFG(2,isurf), &                                         !                   !
             EFG(3,isurf), &                                         !                   !
             [dot_product(dxyz_ds(:,i), dxyz_duv(:,1,isurf)), &      !                   !
             dot_product(dxyz_ds(:,i), dxyz_duv(:,2,isurf))], &      !                   !
             singular )                                              !                   !
     end do ! <------------------------------------------------------+                   !
  end do ! <-----------------------------------------------------------------------------+

  if ( present(curvature) ) then ! <-----------------------------------------------------+
     !! compute curvature of the intersection branch(es) at the current point            !
     do isurf = 1,2 ! <--------------------------------------------------------------+   !
        ! 2nd fundamental form coefficients                                          !   !
        do ivar = 1,3 ! <-----------------------------------------------+            !   !
           call evald2( &                                               !            !   !
                d2xyz_duv2, &                                           !            !   !
                surf(isurf)%ptr, &                                      !            !   !
                uv(:,isurf), &                                          !            !   !
                ivar )                                                  !            !   !
           LMN(ivar,isurf) = dot_product(n(:,isurf), d2xyz_duv2)        !            !   !
        end do ! <------------------------------------------------------+            !   !
     end do ! <----------------------------------------------------------------------+   !
     !                                                                                   !
     do i = 1,ndir ! <---------------------------------------------------------------+   !
        do isurf = 1,2 ! <--------------------------------------------------------+  !   !
           kn(isurf) = ( &                                                        !  !   !
                LMN(1,isurf) * duv_ds(1,i,isurf)**2 + &                           !  !   !
                2._fp * LMN(2,isurf) * duv_ds(1,i,isurf) * duv_ds(2,i,isurf) + &  !  !   !
                LMN(3,isurf) * duv_ds(2,i,isurf)**2 ) / ( &                       !  !   !
                EFG(1,isurf) * duv_ds(1,i,isurf)**2 + &                           !  !   !
                2._fp * EFG(2,isurf) * duv_ds(1,i,isurf) * duv_ds(2,i,isurf) + &  !  !   !
                EFG(3,isurf) * duv_ds(2,i,isurf)**2 )                             !  !   !
        end do ! <----------------------------------------------------------------+  !   !
        denom = max(1._fp - cosn**2, EPSfp)                                          !   !
        curvature(i) = (kn(1)**2 - 2._fp*kn(1)*kn(2)*cosn + kn(2)**2) / denom        !   !
        curvature(i) = sqrt(max(0._fp, curvature(i)))                                !   !
     end do ! <----------------------------------------------------------------------+   !
  end if  ! <----------------------------------------------------------------------------+

end subroutine diffgeom_intersection_curve
