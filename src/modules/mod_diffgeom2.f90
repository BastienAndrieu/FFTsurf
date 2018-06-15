module mod_diffgeom2
  
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

  interface eval
     module procedure eval_curve, eval_surface
  end interface eval
  
  interface evald1
     module procedure evald1_curve, evald1_surface
  end interface evald1

  interface evald2
     module procedure evald2_curve, evald2_surface
  end interface evald2


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

  

  subroutine eval_curve( &
       xyz, &
       curv, &
       t )
    implicit none
    type(type_curve), intent(in)  :: curv
    real(kind=fp),    intent(in)  :: t
    real(kind=fp),    intent(out) :: xyz(3)
    
    call polyval1( xyz, curv%x, [t], 1 )

  end subroutine eval_curve




  subroutine evald1_curve( &
       xyz, &
       curv, &
       t )
    implicit none
    type(type_curve), intent(in)  :: curv
    real(kind=fp),    intent(in)  :: t
    real(kind=fp),    intent(out) :: xyz(3)

    call polyval1( xyz, curv%xt, [t], 1 )

  end subroutine evald1_curve



  
  subroutine evald2_curve( &
       xyz, &
       curv, &
       t )
    implicit none
    type(type_curve), intent(in)  :: curv
    real(kind=fp),    intent(in)  :: t
    real(kind=fp),    intent(out) :: xyz(3)

    call polyval1( xyz, curv%xtt, [t], 1 )

  end subroutine evald2_curve





  subroutine eval_surface( &
       xyz, &
       surf, &
       uv )
    implicit none
    type(type_surface), intent(in)  :: surf
    real(kind=fp),      intent(in)  :: uv(2)
    real(kind=fp),      intent(out) :: xyz(3)
    
    call polyval2( xyz, surf%x, uv, 1 )

  end subroutine eval_surface


  
  subroutine evald1_surface( &
       xyz, &
       surf, &
       uv, &
       ivar )
    implicit none
    type(type_surface), intent(in)  :: surf
    real(kind=fp),      intent(in)  :: uv(2)
    real(kind=fp),      intent(out) :: xyz(3)
    integer,            intent(in)  :: ivar

    select case (ivar)
       case (1)
       call polyval2( xyz, surf%xu, uv, 1 )
    case (2)
       call polyval2( xyz, surf%xv, uv, 1 )
    case default
       STOP 'evald1_surface : ivar =/ 1,2'
    end select
    
  end subroutine evald1_surface



  
  subroutine evald2_surface( &
       xyz, &
       surf, &
       uv, &
       ivar )
    implicit none
    type(type_surface), intent(in)  :: surf
    real(kind=fp),      intent(in)  :: uv(2)
    real(kind=fp),      intent(out) :: xyz(3)
    integer,            intent(in)  :: ivar

    select case (ivar)
    case (1)
       call polyval2( xyz, surf%xuu, uv, 1 )
    case (2)
       call polyval2( xyz, surf%xuv, uv, 1 )
    case (3)
       call polyval2( xyz, surf%xvv, uv, 1 )
    case default
       STOP 'evald1_surface : ivar =/ 1,2,3'
    end select

  end subroutine evald2_surface


end module mod_diffgeom2
