subroutine intersect_gaussmaps_elsewhere( &
     surroot, &
     surfpn1, &
     surfpn2, &
     xyzinter, &
     region1, &
     region2, &
     interdat )
  use mod_math
  use mod_chebyshev
  use mod_separation
  use mod_types_intersection
  implicit none
  type(ptr_parametric_surface), intent(in)    :: surfroot(2)
  type(type_chebyshev_series2), intent(in)    :: surfpn1(:), surfpn2(:)
  real(kind=MATHpr),            intent(in)    :: xyzinter
  type(ptr_surface_region),     intent(inout) :: region1(:), region2(:)
  type(type_intersection_data), intent(inout) :: interdat
  real(kind=MATHpr)                           :: rotation_matrix(3,3)
  real(kind=MATHpr)                           :: xyzsep1((surfpn1(1)%degr(1)+1)*(surfpn1(1)%degr(2)+1),3,size(surfpn1))
  real(kind=MATHpr)                           :: xyzsep2((surfpn2(1)%degr(1)+1)*(surfpn2(1)%degr(2)+1),3,size(surfpn2))
  integer                                     :: n1(size(surfpn1)), n2(size(surfpn2))
  real(kind=MATHpr)                           :: param_vector(3)
  logical                                     :: separable
  integer                                     :: ichild, jchild

  
  rotation_matrix(:,3) = xyzinter / norm2( xyzinter )
  rotation_matrix(:,1) = 0._MATHpr
  rotation_matrix( minloc( abs(xyzinter), 1 ),1) = 1._MATHpr
  rotation_matrix(:,1) = rotation_matrix(:,1) - &
       dot_product( rotation_matrix(:,1), rotation_matrix(:,3) ) * rotation_matrix(:,3)
  rotation_matrix(:,2) = cross( rotation_matrix(:,3), rotation_matrix(:,1) )

  
  do ichild = 1,size(surfpn1)
     call prepare_separability_test_surface( &
          surfpn1(ichild), &
          xyzinter, &
          ch2be_matrices, &
          xyzsep1(:,:,1), &
          n1(ichild) )
     xyzsep1(1:n1(ichild),:,ichild) = matmul( xyzsep1(1:n1(ichild),:,ichild), rotation_matrix )
  end do

  do ichild = 1,size(surfpn2)
     call prepare_separability_test_surface( &
          surfpn2(ichild), &
          xyzinter, &
          ch2be_matrices, &
          xyzsep2(:,:,2), &
          n2(ichild) )
     xyzsep2(1:n2(ichild),:,ichild) = matmul( xyzsep2(1:n2(ichild),:,ichild), rotation_matrix )
  end do

  do jchild = 1,size(surfpn2)
     do ichild = 1,size(surfpn1)
        call separate_spherical_bounding_boxes( &
             xyzsep1(1:n1(ichild),:,ichild), &
             xyzsep2(1:n2(ichild),:,jchild), &
             param_vector, &
             separable, &
             [.false., .false., .true.] )

        if ( separable ) then
           ! call intersect_simple_surfaces ...
           PRINT *,'PARAM_VECTOR =',REAL(PARAM_VECTOR)
           call intersect_simple_surfaces( &
                surfroot, &
                surfc, &
                region, &
                param_vector, &
                interdat )

        else
           ! subdivide and call intersect_surface_surface ...
        end if

     end do
  end do


end subroutine intersect_gaussmaps_elsewhere
