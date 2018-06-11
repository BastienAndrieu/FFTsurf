subroutine compute_tangent_to_intersection_curve( &
     surf, &
     uv, &
     uv_s, &
     xyz_s, &
     stat )
  use mod_math
  use mod_diffgeom
  use mod_errors_intersection
  ! ( see "Shape interrogation for computer aided design and manufacturing", &
  ! Patrikalakis et al. (2009), pp.166-175)
  implicit none
  type(type_parametric_surface), intent(in)  :: surf(2)
  real(kind=MATHpr),             intent(in)  :: uv(2,2)
  real(kind=MATHpr),             intent(out) :: uv_s(2,2,2)
  real(kind=MATHpr),             intent(out) :: xyz_s(3,2)
  integer,                       intent(out) :: stat
  real(kind=MATHpr)                          :: xyz_uv(3,2,2)
  real(kind=MATHpr)                          :: n(3,2), nsqr(2), aux(3)
  real(kind=MATHpr)                          :: errorxyz(3)
  character(200)                             :: errormsg, errorlogfile
  integer                                    :: isurf, ivar, i

  ! compute normal vectors to both surfaces
  do isurf = 1,2
     do ivar = 1,2
        call evald( xyz_uv(:,ivar,isurf), surf(isurf), uv(:,isurf), ivar )
     end do
     n(:,isurf) = cross( xyz_uv(:,1,isurf), xyz_uv(:,2,isurf) )
  end do
  nsqr = sum( n**2, 1 )

  ! if the normals are not collinear, this is a transversal intersection point ...
  stat = 0
  xyz_s(:,1) = cross( n(:,1), n(:,2) )
  xyz_s(:,1) = xyz_s(:,1) / norm2( xyz_s(:,1) )
  
  ! ...else this is a tangential intersection point, the tangent direction(s) to 
  ! the intersection banch(es) is computed from second derivatives of the surfaces
  call characterize_tangential_intersection_point( &
       surf, &
       uv, &
       xyz_s, &
       stat, &
       xyz_uv, &
       n )

  if ( stat < 2 ) then
     ! compute tangential direction(s) in the uv-space of each surface
     do isurf = 1,2
        do i = 1,stat+1
           aux = cross( xyz_s(:,i), n(:,isurf) )
           do ivar = 1,2
              uv_s(ivar,i,isurf) = real( (-1)**ivar, kind=MATHpr) * &
                   dot_product( xyz_uv(:,1+mod(ivar,2),isurf), aux ) / nsqr(isurf)
           end do
        end do
     end do

  else
     ! the tangential direction is undefined, throws an error
     call eval( errorxyz, surf(1), uv(:,1) )

     select case (stat)
     case (2) ! isolated tangential contact point
        errormsg = 'compute_tangent_to_intersection_curve : &
             & isolated tangential contact point --> potential non-manifoldness'
     case (3) ! higher-order contact point
        errormsg = 'compute_tangent_to_intersection_curve : &
             & high-order contact point --> cannot determine nature of intersection'
     end select

     call error_intersection_point( &
          errormsg, &
          surf, &
          uv, &
          errorxyz, &
          errorlogfile )
     print *, 'compute_tangent_to_intersection_curve : see &
          &', trim(errorlogfile)
     STOP

  end if
  


end subroutine compute_tangent_to_intersection_curve
