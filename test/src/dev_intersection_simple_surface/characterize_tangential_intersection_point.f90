subroutine characterize_tangential_intersection_point( &
     surf, &
     uv, &
     s1, &
     n, &
     xyz_s, &
     stat )
  use mod_math
  use mod_diffgeom2
  ! ( see "Shape interrogation for computer aided design and manufacturing", &
  ! Patrikalakis et al. (2009), pp.172-175 )
  ! stat = 0 : single tangential intersection curve
  ! stat = 1 : two tangential intersection branches
  ! stat = 2 : isolated tangential contact point
  ! stat = 3 : higher-order contact point
  implicit none
  real(kind=fp), parameter       :: EPS = EPSmath
  type(ptr_surface), intent(in)  :: surf(2)
  real(kind=fp),     intent(in)  :: uv(2,2)
  real(kind=fp),     intent(in)  :: s1(3,2,2)
  real(kind=fp),     intent(in)  :: n(3,2)
  real(kind=fp),     intent(out) :: xyz_s(3,2)
  integer,           intent(out) :: stat
  real(kind=fp)                  :: s2(3)
  real(kind=fp)                  :: aux(3,2), A(2,2), LMN(3,2), lambda, B(3)
  real(kind=fp)                  :: discr, w(2)
  integer                        :: isurf, ivar

  ! compute A
  do ivar = 1,2
     aux(:,ivar) = cross( n(:,1), s1(:,1+mod(ivar,2),2) )
  end do
  aux(:,1) = -aux(:,1)
  A = matmul( transpose(aux), s1(:,:,1) )

  ! compute second fundamental forms
  do isurf = 1,2
     do ivar = 1,3
        call evald2( &
             s2, &
             surf(isurf)%ptr, &
             uv(:,isurf), &
             ivar )
        LMN(ivar,isurf) = dot_product( n(:,isurf), s2 )
     end do
  end do

  ! compute B = [ b11 ; b12 ; b22 ]
  B(1) = LMN(1,2)*A(1,1)**2 + &
       LMN(2,2)*2*A(1,1)*A(2,1) + &
       LMN(3,2)*A(2,1)**2

  B(2) = LMN(1,2)*A(1,1)*A(1,2) + &
       LMN(2,2)*(A(1,1)*A(2,2) + A(1,2)*A(2,1)) + &
       LMN(3,2)*A(2,1)*A(2,2)

  B(3) = LMN(1,2)*A(1,2)**2 + &
       LMN(2,2)*2*A(1,2)*A(2,2) + &
       LMN(3,2)*A(2,2)**2

  lambda = 1._fp / &
       ( sqrt( product(sum(n**2, 1)) ) * dot_product( n(:,1), n(:,2) ) )

  B = LMN(:,1) - lambda * B

  ! analyze B
  discr = B(2)**2 - B(1)*B(3)

  if ( all( abs(B) < EPS ) ) then
     ! higher order contact point (xyz_s cannot be determined)
     stat = 3
     return

  elseif ( abs( discr ) < EPS ) then
     ! tangential intersection curve (one single xyz_s)
     stat = 0
     if ( abs(B(1)) > abs(B(3)) ) then
        w(1) = -B(2) / B(1)
        xyz_s(:,1) = w(1)*s1(:,1,1) + s1(:,2,1)
     else
        w(1) = -B(2) / B(3)
        xyz_s(:,1) = s1(:,1,1) + w(1)*s1(:,2,1)
     end if

  elseif ( discr < 0._fp ) then
     ! isolated tangential contact point (xyz_s undefined)
     stat = 2
     return

  else
     ! branch point (two disctinct xyz_s)
     stat = 1
     w = real( [-1,1], kind=fp ) * sqrt( discr )
     if ( abs(B(1)) > abs(B(3)) ) then
        w = (w - B(2)) / B(1)
        xyz_s = outer_product( s1(:,1,1), w ) + &
             spread( s1(:,2,1), dim=2, ncopies=2 )
     else
        w = (w - B(2)) / B(3)
        xyz_s = spread( s1(:,1,1), dim=2, ncopies=2 ) + &
             outer_product( s1(:,2,1), w )
     end if

  end if

  ! unitize xyz_s
  do ivar = 1,stat+1
     xyz_s(:,ivar) = xyz_s(:,ivar) / norm2( xyz_s(:,ivar) )
  end do

end subroutine characterize_tangential_intersection_point
