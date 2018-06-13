! *** A PRIORI NOT USED ***
subroutine intersect_elsewhere_surface_surface( &
     surfc, &
     xyzinter, &
     transfmat, &
     separable )
  use mod_math
  use mod_chebyshev
  use mod_separation
  use mod_types_intersection
  ! Hohmeyer's spherical separability test between two surfaces that 
  ! are know to intersect at corners.
  ! ( see "Robust and Efﬁcient Surface Intersection for Solid Modeling", 
  ! M.E. Hohmeyer (1986), pp.61-64. ).
  ! Returns 'separable'=.true. if the surfaces are guaranteed not to intersect 
  ! at any other point.
  implicit none
  type(type_chebyshev_series2),  intent(in)    :: surfc(2)
  real(kind=MATHpr),             intent(in)    :: xyzinter(3)
  type(type_matrix),             intent(inout) :: transfmat(:)
  logical,                       intent(out)   :: separable
  real(kind=MATHpr)                            :: xyzsep1((surfc(1)%degr(1)+1)*(surfc(1)%degr(2)+1)-1,3)
  real(kind=MATHpr)                            :: xyzsep2((surfc(2)%degr(1)+1)*(surfc(2)%degr(2)+1)-1,3)
  integer                                      :: n1, n2
  real(kind=MATHpr)                            :: separation_vector(3)
  
  call prepare_separability_test_surface( &
       surfc(1), &
       xyzinter, &
       transfmat, &
       xyzsep1, &
       n1 )

  call prepare_separability_test_surface( &
       surfc(2), &
       xyzinter, &
       transfmat, &
       xyzsep2, &
       n2 )

  ! run separability test
  call separating_plane( &
       xyzsep1(1:n1,:), &
       xyzsep2(1:n2,:), &
       n1, &
       n2, &
       separation_vector, &
       separable )

end subroutine intersect_elsewhere_surface_surface
