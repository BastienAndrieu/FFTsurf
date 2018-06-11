subroutine prepare_separability_test_curve( &
    curvc, &
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
  type(type_chebyshev_series1),  intent(in)    :: curvc
  real(kind=MATHpr),             intent(in)    :: xyzinter(3)
  type(type_matrix),             intent(inout) :: transfmat(:)
  real(kind=MATHpr),             intent(out)   :: xyzsep(:,:)
  integer,                       intent(out)   :: n
  real(kind=MATHpr)                            :: c2b(curvc%degr+1,curvc%degr+1)
  real(kind=MATHpr)                            :: bcp(curvc%degr+1,3)
  real(kind=MATHpr)                            :: xyzi(3)
  integer                                      :: i

  ! get adequate Chebyshev to Bernstein transformation matrix
  call get_ch2be_from_collection( &
       transfmat, &
       curvc%degr + 1, &
       c2b )

  ! compute Bezier control points
  bcp = matmul( c2b, curvc%coef )

  ! rearrange and translate control points for separability test
  n = 0
  do i = 1,curvc%degr+1
     xyzi = bcp(i,:) - xyzinter
     if ( sum(xyzi**2) > EPSxyzsqr ) then
        n = n + 1
        if ( n > size(xyzsep,1) ) STOP 'prepare_separability_test_curve :&
             & n > size(xyzsep,1)'
        xyzsep(n,:) = xyzi
     end if
  end do

end subroutine prepare_separability_test_curve
