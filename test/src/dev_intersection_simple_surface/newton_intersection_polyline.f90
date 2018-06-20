subroutine newton_intersection_polyline( &
     surf, &
     !uvbox, &
     xyz_prev, &
     htargetsqr, &
     tolhsqr, &
     uv, &
     xyz, &
     stat )     
  use mod_math
  use mod_diffgeom2
  use mod_tolerances
  implicit none
  integer, parameter               :: nitmax = 10
  type(ptr_surface), intent(in)    :: surf(2)
  !real(kind=fp),     intent(in)    :: uvbox(4,2)
  real(kind=fp),     intent(in)    :: xyz_prev(3)
  real(kind=fp),     intent(in)    :: htargetsqr
  real(kind=fp),     intent(in)    :: tolhsqr
  real(kind=fp),     intent(inout) :: uv(2,2)
  real(kind=fp),     intent(out)   :: xyz(3)
  integer,           intent(out)   :: stat
  !real(kind=fp), dimension(4)      :: lowerb, upperb, rng
  real(kind=fp)                    :: s(3,2), s1(3,2,2)
  real(kind=fp), dimension(3)      :: r1, r2
  real(kind=fp)                    :: resx, resh
  real(kind=fp)                    :: jac(4,4), duv(4)!, lambda
  logical                          :: singular
  integer                          :: it, isurf, ivar

  stat = 1
  !lowerb = reshape( uvbox(1:2,1:2), [4] )
  !upperb = reshape( uvbox(3:4,1:2), [4] )
  !rng = upperb - lowerb
  !lowerb = lowerb - EPsuv*rng
  !upperb = upperb + EPsuv*rng
  jac(:,:) = 0._fp

  do it = 1,nitmax

     do isurf = 1,2
        call eval( &
             s(:,isurf), &
             surf(isurf)%ptr, &
             uv(:,isurf) )
     end do

     r1 = s(:,1) - s(:,2)
     r2 = s(:,1) - xyz_prev

     resx = sum( r1**2 )
     resh = sum( r2**2 ) - htargetsqr

     ! convergence criterion
     if ( resx < EPSxyzsqr .and. abs(resh) < tolhsqr ) then
        stat = 0
        xyz = 0.5_fp * sum( s, 2 )
        return     
     end if

     ! Jacobian matrix
     do isurf = 1,2
        do ivar = 1,2
           call evald1( &
                s1(:,ivar,isurf), &
                surf(isurf)%ptr, &
                uv(:,isurf), &
                ivar )
        end do
        jac(1:3,2*isurf-[1,0]) = real( (-1)**(isurf+1), kind=fp ) * s1(1:3,1:2,isurf)
        if ( isurf == 1 ) then
           do ivar = 1,2
              jac(4,ivar) = 2._fp * dot_product( s1(1:3,ivar,1), r2 )
           end do
        end if
     end do
     
     ! solve for Newton step
     call linsolve_QR( &
          duv, &
          jac, &
          -[r1,resh], &
          4, &
          4, &
          singular )
     
     if ( singular ) then
        ! singular Jacobian matrix (should not happen)
        stat = 2
        return
     end if

     !! scale down Newton step to keep the solution inside feasible region
     !call nd_box_constraint( &
     !     reshape( uv, [4] ), &
     !     lowerb, &
     !     upperb, &
     !     duv, &
     !     lambda )
     !if ( lambda < -EPSmath ) return ! negative damped factor
     !duv = lambda * duv
     if ( sum(duv(1:2)**2) < EPSuvsqr .or. sum(duv(3:4)**2) < EPSuvsqr ) then
        ! damped Newton step is too small
        return
     end if

     ! update solution
     uv(:,1) = uv(:,1) + duv(1:2)
     uv(:,2) = uv(:,2) + duv(3:4)
     
  end do


end subroutine newton_intersection_polyline
