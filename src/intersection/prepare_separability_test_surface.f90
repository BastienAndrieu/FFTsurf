subroutine prepare_separability_test_surface( &
    surfc, &
    xyzinter, &
    transfmat, &
    xyzsep, &
    n )
  use mod_math
  use mod_chebyshev
  use mod_types_intersection
  use mod_tolerances
  ! Computes data necessary for performing Hohmeyer's separability test
  implicit none
  type(type_chebyshev_series2),  intent(in)    :: surfc
  real(kind=MATHpr),             intent(in)    :: xyzinter(3)
  type(type_matrix),             intent(inout) :: transfmat(:)
  real(kind=MATHpr),             intent(out)   :: xyzsep(:,:)
  integer,                       intent(out)   :: n
  real(kind=MATHpr)                            :: c2bu(surfc%degr(1)+1,surfc%degr(1)+1)
  real(kind=MATHpr)                            :: c2bv(surfc%degr(2)+1,surfc%degr(2)+1)
  real(kind=MATHpr)                            :: bcp(surfc%degr(1)+1,surfc%degr(2)+1,3)
  real(kind=MATHpr)                            :: xyzi(3)
  integer                                      :: i, j

  ! get adequate Chebyshev to Bernstein transformation matrices
  call get_ch2be_from_collection( &
       transfmat, &
       surfc%degr(1)+1, &
       c2bu )
  call get_ch2be_from_collection( &
       transfmat, &
       surfc%degr(2)+1, &
       c2bv )
  c2bv = transpose( c2bv )

  ! compute Bezier control points
  do i = 1,3
     bcp(:,:,i) = matmul( matmul( c2bu, surfc%coef(:,:,i) ), c2bv )
  end do

  ! rearrange and translate control points for separability test
  n = 0
  do j = 1,surfc%degr(2)+1
     do i = 1,surfc%degr(1)+1
        xyzi = bcp(i,j,:) - xyzinter
        if ( sum(xyzi**2) > EPSxyzsqr ) then
           n = n + 1
           if ( n > size(xyzsep,1) ) STOP 'prepare_separability_test_surface :&
                & n > size(xyzsep,1)'
           xyzsep(n,:) = xyzi
        end if
     end do
  end do

end subroutine prepare_separability_test_surface
