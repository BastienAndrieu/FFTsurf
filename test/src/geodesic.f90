program geodesic

  use mod_util
  use mod_math
  use mod_polynomial
  use mod_diffgeom

  implicit none
  
  integer, parameter :: nmax = 10000
  real(kind=fp), parameter :: tolchord = 5.d-4
  real(kind=fp), parameter :: frac_chord = 2._fp*sqrt(tolchord*(2._fp - tolchord))
  real(kind=fp), parameter :: hmin = 1.d-4
  real(kind=fp), parameter :: hmax = 1.d-1
  
  
  
  integer            :: narg
  character(100)     :: arg
  
  type(type_surface) :: surf
  real(kind=fp)      :: uv(2,nmax)
  real(kind=fp)      :: angle, duv_ds(2), d2uv_ds2(2)
  real(kind=fp)      :: dxyz_duv(3,2), dxyz_ds(3), scale
  real(kind=fp)      :: d2xyz_duv2(3,3), EFG(3), LMN(3), kn, h
  logical            :: singular
  real(kind=fp)      :: step = 1.d-3
  integer            :: ipoint, n
  integer            :: fileunit

  call read_polynomial( &
       surf%x, &
       'geodesic/c.cheb', &
       nvar=2, &
       base=1 )
  call economize2(surf%x, EPSmath)
  call compute_deriv1(surf)
  call compute_deriv2(surf)

  ! intial point and direction
   narg = command_argument_count()
  if ( narg < 2 ) then
     uv(1:2,1) = 0._fp
  else
     call get_command_argument(1, arg)
     read (arg,*) uv(1,1)
     call get_command_argument(2, arg)
     read (arg,*) uv(2,1)
  end if
  
  if ( narg < 3 ) then
     call random_number(duv_ds)
     duv_ds = 2._fp*duv_ds - 1._fp
  else
     call get_command_argument(3, arg)
     read (arg,*) angle
     angle = angle*CSTpi/180._fp
     duv_ds(1) = cos(angle)
     duv_ds(2) = sin(angle)
  end if
  call normalize_duv_ds(surf, uv(1:2,1), duv_ds)

  ! integrate geodesic
  do ipoint = 1,nmax-1     
     if ( .false. ) then
        call eval_d2uv_ds2(surf, uv(1:2,ipoint), duv_ds, singular, d2uv_ds2)
        if ( singular ) exit

        duv_ds = duv_ds + step*d2uv_ds2
        call normalize_duv_ds(surf, uv(1:2,ipoint), duv_ds)
        uv(1:2,ipoint+1) = uv(1:2,ipoint) + step*duv_ds
     else
        
        call evald1_all(surf, uv(1:2,ipoint), dxyz_duv)
        call first_fundamental_form_coeffs(dxyz_duv, EFG)
        dxyz_ds = matmul(dxyz_duv, duv_ds)
        scale = 1._fp/norm2(dxyz_ds)
        dxyz_ds = dxyz_ds*scale
        duv_ds = duv_ds*scale

        call evald2_all(surf, uv(1:2,ipoint), d2xyz_duv2)
        call second_fundamental_form_coeffs(dxyz_duv, d2xyz_duv2, LMN)
        call normal_curvature(duv_ds, EFG, LMN, kn)

        h = frac_chord/abs(kn)
        h = min(hmax, max(hmin, h))
        
        uv(1:2,ipoint+1) = uv(1:2,ipoint) + h*duv_ds

        call evald1_all(surf, uv(1:2,ipoint+1), dxyz_duv)
        call project_onto_uv(dxyz_duv, dxyz_ds, singular, duv_ds)
        
        if ( singular ) exit
        
     end if

     if ( any(abs(uv) > 1._fp) ) exit
  end do
  n = ipoint

  call get_free_unit(fileunit)
  open(unit=fileunit, file='geodesic/geodesic_uv.dat', action='write')
  do ipoint = 1,n
     write (fileunit, *) uv(1:2,ipoint)
  end do
  close(fileunit)
  
  
contains

  subroutine normal_curvature(duv, EFG, LMN, kn)
    implicit none
    real(kind=fp), intent(in)  :: duv(2)
    real(kind=fp), intent(in)  :: EFG(3)
    real(kind=fp), intent(in)  :: LMN(3)
    real(kind=fp), intent(out) :: kn

    kn = (LMN(1)*duv(1)**2 + 2._fp*LMN(2)*duv(1)*duv(2) + LMN(3)*duv(2)**2) / &
         (EFG(1)*duv(1)**2 + 2._fp*EFG(2)*duv(1)*duv(2) + EFG(3)*duv(2)**2)
    
  end subroutine normal_curvature
  

  subroutine project_onto_uv(dxyz_duv, vec, singular, duv)
    implicit none
    real(kind=fp), intent(in)  :: dxyz_duv(3,2)
    real(kind=fp), intent(in)  :: vec(3)
    logical,       intent(out) :: singular
    real(kind=fp), intent(out) :: duv(2)
    real(kind=fp)              :: EFG(3)
    
    call first_fundamental_form_coeffs(dxyz_duv, EFG)
    call solve_2x2( &
         duv, &
         EFG(1), &
         EFG(2), &
         EFG(2), &
         EFG(3), &
         matmul(transpose(dxyz_duv), vec), &
         singular )
        
  end subroutine project_onto_uv
  

  subroutine normalize_duv_ds(surf, uv, duv_ds)
    implicit none
    type(type_surface), intent(in)  :: surf
    real(kind=fp),      intent(in)  :: uv(2)
    real(kind=fp),      intent(out) :: duv_ds(2)
    real(kind=fp)                   :: dxyz_duv(3,2), dxyz_ds(3)
    
    call evald1_all(surf, uv, dxyz_duv)
    dxyz_ds = matmul(dxyz_duv, duv_ds)
    duv_ds = duv_ds/norm2(dxyz_ds)
    
  end subroutine normalize_duv_ds
  
  

  subroutine evald1_all(surf, uv, dxyz_duv)
    implicit none
    type(type_surface), intent(in)  :: surf
    real(kind=fp),      intent(in)  :: uv(2)
    real(kind=fp),      intent(out) :: dxyz_duv(3,2)
    integer                         :: ivar

    do ivar = 1,2
       call evald1( &
            dxyz_duv(1:3,ivar), &
            surf, &
            uv, &
            ivar )
    end do
  end subroutine evald1_all



  subroutine evald2_all(surf, uv, d2xyz_duv2)
    implicit none
    type(type_surface), intent(in)  :: surf
    real(kind=fp),      intent(in)  :: uv(2)
    real(kind=fp),      intent(out) :: d2xyz_duv2(3,3)
    integer                         :: ivar

    do ivar = 1,3
       call evald2( &
            d2xyz_duv2(1:3,ivar), &
            surf, &
            uv, &
            ivar )
    end do
  end subroutine evald2_all



  subroutine first_fundamental_form_coeffs(dxyz_duv, EFG)
    implicit none
    real(kind=fp), intent(in)  :: dxyz_duv(3,2)
    real(kind=fp), intent(out) :: EFG(3)

    EFG(1) = dot_product(dxyz_duv(1:3,1), dxyz_duv(1:3,1))
    EFG(2) = dot_product(dxyz_duv(1:3,1), dxyz_duv(1:3,2))
    EFG(3) = dot_product(dxyz_duv(1:3,2), dxyz_duv(1:3,2))
  end subroutine first_fundamental_form_coeffs
  

  
  subroutine second_fundamental_form_coeffs(dxyz_duv, d2xyz_duv2, LMN)
    implicit none
    real(kind=fp), intent(in)  :: dxyz_duv(3,2)
    real(kind=fp), intent(in)  :: d2xyz_duv2(3,3)
    real(kind=fp), intent(out) :: LMN(3)
    real(kind=fp)              :: n(3)
    integer                    :: i

    n = cross(dxyz_duv(1:3,1), dxyz_duv(1:3,2))
    n = n/norm2(n)

    do i = 1,3
       LMN(i) = dot_product(n, d2xyz_duv2(1:3,i))
    end do
  end subroutine second_fundamental_form_coeffs

  

  subroutine christoffel_symbols(surf, uv, singular, csymb)
    implicit none
    type(type_surface), intent(in)  :: surf
    real(kind=fp),      intent(in)  :: uv(2)
    logical,            intent(out) :: singular
    real(kind=fp),      intent(out) :: csymb(3,2)
    real(kind=fp)                   :: dxyz_duv(3,2), d2xyz_duv2(3,3)
    real(kind=fp), dimension(3)     :: EFG, dEFG_du, dEFG_dv
    real(kind=fp)                   :: detEFG
    
    ! first derivatives
    call evald1_all(surf, uv, dxyz_duv)

    ! second derivatives
    call evald2_all(surf, uv, d2xyz_duv2)

    ! first fundamental form coefficients
    call first_fundamental_form_coeffs(dxyz_duv, EFG)

    ! first fundamental form determinant
    detEFG = EFG(1)*EFG(3) - EFG(2)**2
    singular = detEFG < EPSmath
    if ( singular ) return

    ! u-derivatives of first fundamental form coefficients
    dEFG_du(1) = 2._fp*dot_product(dxyz_duv(1:3,1), d2xyz_duv2(1:3,1))
    dEFG_du(2) = dot_product(d2xyz_duv2(1:3,1), dxyz_duv(1:3,2)) + &
         dot_product(dxyz_duv(1:3,1), d2xyz_duv2(1:3,2))
    dEFG_du(3) = 2._fp*dot_product(dxyz_duv(1:3,2), d2xyz_duv2(1:3,2))

    ! v-derivatives of first fundamental form coefficients
    dEFG_dv(1) = 2._fp*dot_product(dxyz_duv(1:3,1), d2xyz_duv2(1:3,2))
    dEFG_dv(2) = dot_product(d2xyz_duv2(1:3,2), dxyz_duv(1:3,2)) + &
         dot_product(dxyz_duv(1:3,1), d2xyz_duv2(1:3,3))
    dEFG_dv(3) = 2._fp*dot_product(dxyz_duv(1:3,2), d2xyz_duv2(1:3,3))
   
    ! Christoffel symbols
    csymb(1,1) = EFG(3)*dEFG_du(1) - 2._fp*EFG(2)*dEFG_du(2) + EFG(2)*dEFG_dv(1)
    csymb(2,1) = EFG(3)*dEFG_dv(1) - EFG(2)*dEFG_du(3)
    csymb(3,1) = 2._fp*EFG(3)*dEFG_dv(2) - EFG(3)*dEFG_du(3) + EFG(2)*dEFG_dv(3)
    csymb(1,2) = 2._fp*EFG(1)*dEFG_du(2) - EFG(1)*dEFG_dv(1) + EFG(2)*dEFG_du(1)
    csymb(2,2) = EFG(1)*dEFG_du(3) - EFG(2)*dEFG_dv(1)
    csymb(3,2) = EFG(1)*dEFG_dv(3) - 2._fp*EFG(2)*dEFG_dv(2) + EFG(2)*dEFG_du(3)

    csymb = 0.5_fp*csymb/detEFG
    
  end subroutine christoffel_symbols




  subroutine eval_d2uv_ds2(surf, uv, duv_ds, singular, d2uv_ds2)
    implicit none
    type(type_surface), intent(in)  :: surf
    real(kind=fp),      intent(in)  :: uv(2)
    real(kind=fp),      intent(in)  :: duv_ds(2)
    logical,            intent(out) :: singular
    real(kind=fp),      intent(out) :: d2uv_ds2(2)
    real(kind=fp)                   :: csymb(3,2)
    integer                         :: ivar

    call christoffel_symbols(surf, uv, singular, csymb)

    do ivar = 1,2
       d2uv_ds2(ivar) = -( &
            csymb(1,ivar)*duv_ds(1)**2 + &
            2._fp*csymb(2,ivar)*duv_ds(1)*duv_ds(2) + &
            csymb(2,ivar)*duv_ds(2)**2 )
    end do
    
  end subroutine eval_d2uv_ds2
  
end program geodesic
