module mod_diffgeomnew
  
  use mod_polynomial

  type type_curve
     type(type_polynomial) :: x, xt, xtt
  end type type_curve

  type type_surface
     type(type_polynomial) :: x, xu, xv, pn, xuu, xuv, xvv
  end type type_surface

  type ptr_surface
     type(type_surface), pointer :: ptr => null()
  end type ptr_surface
  


  

  interface compute_deriv1
     module procedure compute_deriv1_curve, compute_deriv1_surface
  end interface compute_deriv1
     
  interface compute_deriv2
     module procedure compute_deriv2_curve, compute_deriv2_surface
  end interface compute_deriv2




contains


  subroutine compute_deriv1_curve( curv )
    implicit none
    type(type_curve), intent(inout) :: curv

    call diff1( curv%x, curv%xt )

  end subroutine compute_deriv1_curve



  
  subroutine compute_deriv1_surface( surf )
    implicit none
    type(type_surface), intent(inout) :: surf

    call diff2( surf%x, du=surf%xu, dv=surf%xv )

  end subroutine compute_deriv1_surface


  
  
  subroutine compute_deriv2_curve( curv )
    implicit none
    type(type_curve), intent(inout) :: curv

    call diff1( curv%xt, curv%xtt )

  end subroutine compute_deriv2_curve




  subroutine compute_deriv2_surface( surf )
    implicit none
    type(type_surface), intent(inout) :: surf

    call diff2( surf%xu, du=surf%xuu, dv=surf%xuv )
    call diff2( surf%xv, dv=surf%xvv )

  end subroutine compute_deriv2_surface

  

  



end module mod_diffgeomnew
