module mod_diffgeom
  
  use mod_math
  use mod_chebyshev

  implicit none
  
  type type_parametric_surface
     type(type_chebyshev_series2) :: s, su, sv, suu, suv, svv, pn
  end type type_parametric_surface

  type ptr_parametric_surface
     type(type_parametric_surface), pointer :: ptr => null()
  end type ptr_parametric_surface

  type type_parametric_curve
     type(type_chebyshev_series1) :: c, ct, ctt
  end type type_parametric_curve

  interface eval
     module procedure ceval_point, seval_point
  end interface eval

  interface evald
     module procedure cevald_point, sevald_point
  end interface evald

  interface evald2
     module procedure evald2_point
  end interface evald2

  interface compute_first_derivatives
     module procedure curv_compute_first_derivative, surf_compute_first_derivatives
  end interface compute_first_derivatives

  interface compute_second_derivatives
     module procedure surf_compute_second_derivatives
  end interface compute_second_derivatives

contains

  subroutine curv_compute_first_derivative( curv )
    implicit none
    type(type_parametric_curve), intent(inout) :: curv

    call chebdiff1( curv%c, curv%ct )
    
  end subroutine curv_compute_first_derivative


  subroutine surf_compute_first_derivatives( surf )
    implicit none
    type(type_parametric_surface), intent(inout) :: surf
    
    call chebdiff2( surf%s, surf%su, surf%sv )

  end subroutine surf_compute_first_derivatives



  subroutine surf_compute_second_derivatives( surf )
    implicit none
    type(type_parametric_surface), intent(inout) :: surf
    
    call chebdiff2( surf%su, surf%suu )
    call chebdiff2( surf%sv, surf%suv, surf%svv )

  end subroutine surf_compute_second_derivatives


  
  subroutine compute_pseudonormal( surf )
    implicit none
    type(type_parametric_surface), intent(inout)                      :: surf
    real(kind=MATHpr), dimension(2*surf%s%degr(1),2*surf%s%degr(2),3) :: xu, xv, xn
    integer                                                           :: k

    if ( surf%s%dim < 3 ) STOP 'compute_pseudonormal : dim < 3'
    
    call ifcht2_padded( &
         surf%su, &
         xu, &
         2*surf%s%degr(1), &
         2*surf%s%degr(2) )
    call ifcht2_padded( surf%sv, &
         xv, &
         2*surf%s%degr(1), &
         2*surf%s%degr(2) )
    
    do k = 1,3
       xn(:,:,k) = xu(:,:,1+mod(k,3)) * xv(:,:,1+mod(k+1,3)) - &
            xu(:,:,1+mod(k+1,3)) * xv(:,:,1+mod(k,3))
    end do
    
    call fcht2( &
         xn, &
         surf%pn, &
         2*surf%s%degr(1), &
         2*surf%s%degr(2), &
         3 )
    
  end subroutine compute_pseudonormal







  subroutine ceval_point( &
       xyz, &
       curv, &
       t )
    implicit none
    type(type_parametric_curve), intent(in)  :: curv
    real(kind=MATHpr),           intent(in)  :: t
    real(kind=MATHpr),           intent(out) :: xyz(3)

    call chebval( xyz, curv%c, t )

  end subroutine ceval_point

  
  subroutine seval_point( &
       xyz, &
       surf, &
       uv )
    implicit none
    type(type_parametric_surface), intent(in)  :: surf
    real(kind=MATHpr),             intent(in)  :: uv(2)
    real(kind=MATHpr),             intent(out) :: xyz(3)

    call chebval( xyz, surf%s, uv )

  end subroutine seval_point



  subroutine cevald_point( &
       xyz, &
       curv, &
       t )
    implicit none
    type(type_parametric_curve), intent(in)  :: curv
    real(kind=MATHpr),           intent(in)  :: t
    real(kind=MATHpr),           intent(out) :: xyz(3)

    call chebval( xyz, curv%ct, t )

  end subroutine cevald_point



  subroutine sevald_point( &
       xyz, &
       surf, &
       uv, &
       ivar )
    implicit none
    type(type_parametric_surface), intent(in)  :: surf
    real(kind=MATHpr),             intent(in)  :: uv(2)
    integer,                       intent(in)  :: ivar
    real(kind=MATHpr),             intent(out) :: xyz(3)

    if ( ivar == 2 ) then
       call chebval( xyz, surf%sv, uv )
    else
       call chebval( xyz, surf%su, uv )
    end if

  end subroutine sevald_point


  subroutine evald2_point( &
       xyz, &
       surf, &
       uv, &
       ivar )
    implicit none
    type(type_parametric_surface), intent(in)  :: surf
    real(kind=MATHpr),             intent(in)  :: uv(2)
    integer,                       intent(in)  :: ivar
    real(kind=MATHpr),             intent(out) :: xyz(3)

    if ( ivar == 3 ) then
       call chebval( xyz, surf%svv, uv )
    elseif ( ivar == 2 ) then
       call chebval( xyz, surf%suv, uv )
    else
       call chebval( xyz, surf%suu, uv )
    end if

  end subroutine evald2_point

end module mod_diffgeom
