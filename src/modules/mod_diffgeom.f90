module mod_diffgeom
  
  use mod_math
  use mod_polynomial

  implicit none

  type type_curve
     type(type_polynomial) :: x, xt, xtt
  end type type_curve

  type type_surface
     type(type_polynomial) :: x, xu, xv, pn, xuu, xuv, xvv
     integer               :: tag = 1
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

  


  subroutine compute_pseudonormal( surf )
    use mod_chebyshev
    implicit none
    type(type_surface), intent(inout) :: surf
    real(kind=fp), allocatable        :: xu(:,:,:), xv(:,:,:), pn(:,:,:)
    integer                           :: mu, nu, mv, nv, p, q, k

    mu = surf%xu%degr(1)+1
    nu = surf%xu%degr(2)+1
    mv = surf%xv%degr(1)+1
    nv = surf%xv%degr(2)+1

    p = mu + mv
    q = nu + nv

    allocate( xu(p,q,3), xv(p,q,3), pn(p,q,3) )
    
    call ifcht2( &
         surf%xu%coef(1:surf%xu%degr(1)+1,1:surf%xu%degr(2)+1,1:3), &
         xu, &
         p, &
         q, &
         3 )
    call ifcht2( &
         surf%xv%coef(1:surf%xv%degr(1)+1,1:surf%xv%degr(2)+1,1:3), &
         xv, &
         p, &
         q, &
         3 )
    
    do k = 1,3
       pn(:,:,k) = &
            xu(:,:,1+mod(k,3))   * xv(:,:,1+mod(k+1,3)) - &
            xu(:,:,1+mod(k+1,3)) * xv(:,:,1+mod(k,3))
    end do

    call reset_polynomial( &
         poly=surf%pn, &
         nvar=2, &
         base=1, &
         degr=[p-1,q-1], &
         dim=3 )

    call fcht2( &
         pn, &
         surf%pn%coef(1:p,1:q,1:3), &
         p, &
         q, &
         3 )

    deallocate( xu, xv, pn )

  end subroutine compute_pseudonormal








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


  subroutine eval_minimum_curvature_radius( &
       surf, &
       uv, &
       r )
    implicit none
    type(type_surface), intent(in)  :: surf
    real(kind=fp),      intent(in)  :: uv(2)
    real(kind=fp),      intent(out) :: r
    real(kind=fp)                   :: kmean, kgauss, k1, k2

    call eval_curvature_surface( &
         surf, &
         uv, &
         kmean, &
         kgauss, &
         k1, &
         k2 )

    r = min(abs(1._fp/k1), abs(1._fp/k2))
    
  end subroutine eval_minimum_curvature_radius

  

  subroutine eval_curvature_surface( &
       surf, &
       uv, &
       kmean, &
       kgauss, &
       k1, &
       k2 )
    implicit none
    type(type_surface), intent(in)  :: surf
    real(kind=fp),      intent(in)  :: uv(2)
    real(kind=fp),      intent(out) :: kmean
    real(kind=fp),      intent(out) :: kgauss
    real(kind=fp),      intent(out) :: k1
    real(kind=fp),      intent(out) :: k2
    real(kind=fp)                   :: dxyz_duv(3,2), EFG(3), invdetI, n(3)
    real(kind=fp)                   :: d2xyz_duv2(3), LMN(3), s
    integer                         :: ivar, jvar

    ! First Fundamental Form (metric tensor)
    do ivar = 1,2
       call evald1( &
            dxyz_duv(:,ivar), &
            surf, &
            uv, &
            ivar )
    end do
    do ivar = 1,2
       do jvar = ivar,2 
          EFG(ivar + jvar - 1) = dot_product(dxyz_duv(:,ivar), dxyz_duv(:,jvar))
       end do
    end do
    invdetI = EFG(1)*EFG(3) - EFG(2)**2
    if ( invdetI < EPSmath ) then
       ! the surface is singular at that point
       kmean = huge(1._fp)
       kgauss = kmean
       k1 = kmean
       k2 = kmean
       return
    else
       ! regular surface
       invdetI = 1._fp / invdetI
    end if

    ! unit normal
    n = cross(dxyz_duv(:,1), dxyz_duv(:,2)) * sqrt(invdetI)
    
    ! Second Fundamental Form (shape operator)
    do ivar = 1,3
       call evald2( &
            d2xyz_duv2, &
            surf, &
            uv, &
            ivar )
       LMN(ivar) = dot_product(n, d2xyz_duv2)
    end do

    ! mean curvature
    kmean = 0.5_fp * invdetI * (EFG(1)*LMN(3) + EFG(3)*LMN(1) - 2._fp*EFG(2)*LMN(2))

    ! Gaussian curvature
    kgauss = (LMN(1)*LMN(3) - LMN(2)**2) * invdetI

    ! principal curvatures
    s = kmean**2 - kgauss
    s = sqrt(max(0._fp, s))

    k1 = kmean + s
    k2 = kmean - s
    
  end subroutine eval_curvature_surface



  subroutine eval_curvature_curve( &
       curv, &
       w, &
       k, &
       d1opt, &
       d2opt )
    implicit none
    type(type_curve),        intent(in)  :: curv
    real(kind=fp),           intent(in)  :: w
    real(kind=fp),           intent(out) :: k
    real(kind=fp), optional, intent(in)  :: d1opt(curv%x%dim)
    real(kind=fp), optional, intent(in)  :: d2opt(curv%x%dim)
    real(kind=fp), dimension(curv%x%dim) :: d1, d2
    real(kind=fp)                        :: d1sqr

    if ( present(d1opt) ) then
       d1 = d1opt
    else
       call evald1( &
            d1, &
            curv, &
            w )
    end if

    if ( present(d2opt) ) then
       d2 = d2opt
    else
       call evald2( &
            d2, &
            curv, &
            w )
    end if

    d1sqr = sum(d1**2)

    select case (size(d1))
    case (2)
       k = (d1(1)*d2(2) - d1(2)*d2(1))/(d1sqr*sqrt(d1sqr))
    case (3)
       k = norm2(cross(d1, d2))/(d1sqr*sqrt(d1sqr))
    case default
       STOP 'eval_curvature_curve: curv%x%dim must equal 2 or 3'
    end select

  end subroutine eval_curvature_curve

end module mod_diffgeom
