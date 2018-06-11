subroutine find_collineal_corners( &
     surf, &
     uv, &
     stat )
  use mod_math
  use mod_chebyshev
  use mod_diffgeom
  use mod_tolerances
  ! Searches for a pair of collineal points among the 4x4 pairs of corners 
  ! of two rectangular parametric surfaces.
  ! Returns 'stat' > 0 if no such pair has been found;
  !         'stat' =-1 if a tangential contact point has been found;
  !         'stat' = 0 else (pair of non-coincident collineal points);
  ! as well as the uv-coordniates of those potential points.
  implicit none
  type(ptr_parametric_surface), intent(in)  :: surf(2)
  real(kind=MATHpr),            intent(out) :: uv(2,2)
  integer,                      intent(out) :: stat
  real(kind=MATHpr)                         :: s(3,4,2), s1(3,4,2,2), f(4)
  real(kind=MATHpr), dimension(3)           :: n, r
  integer                                   :: i, j

  stat = 1

  do i = 1,2
     call chebval2_corners( surf(i)%ptr%s,  s(:,:,i)    )
     call chebval2_corners( surf(i)%ptr%su, s1(:,:,1,i) )
     call chebval2_corners( surf(i)%ptr%sv, s1(:,:,2,i) )
  end do

  do i = 1,4
     n = cross( s1(:,i,1,1), s1(:,i,2,1) )
     do j = 1,4
        r = s(:,i,1) - s(:,j,2)
        f(1) = dot_product( n, s1(:,j,1,2) )
        f(2) = dot_product( n, s1(:,j,2,2) )
        f(3) = dot_product( r, s1(:,i,1,1) )
        f(4) = dot_product( r, s1(:,i,2,1) )
        if ( sum(f**2) < EPScollinealsqr ) then
           uv(:,1) = real( (-1) ** [ i, 1 + i/3 ], kind=MATHpr )
           uv(:,2) = real( (-1) ** [ j, 1 + j/3 ], kind=MATHpr )
           if ( sum(r**2) < EPSxyzsqr ) then
              stat = -1
           else
              stat = 0
           end if
           return
        end if
     end do
  end do

end subroutine find_collineal_corners
