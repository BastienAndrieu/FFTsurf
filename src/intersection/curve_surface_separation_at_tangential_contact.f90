subroutine curve_surface_separation_at_tangential_contact( &
     surf, &
     uv, &
     xyz, &
     bc, &
     bs, &
     separable )
  use mod_math
  use mod_polynomial
  use mod_diffgeom
  use mod_tolerances
  implicit none
  type(type_surface),    intent(in)  :: surf
  real(kind=fp),         intent(in)  :: uv(2)
  real(kind=fp),         intent(in)  :: xyz(3)
  type(type_polynomial), intent(in)  :: bc, bs
  logical,               intent(out) :: separable
  real(kind=fp), dimension(3)        :: su, sv, n
  real(kind=fp)                      :: dotn, mn(2), mx(2)
  integer                            :: i, j

  !CALL WRITE_POLYNOMIAL(BC, 'dev_intersection/debugcsstc_bc.bern')
  !CALL WRITE_POLYNOMIAL(BS, 'dev_intersection/debugcsstc_bs.bern')

  call evald1(su, surf, uv, 1)
  call evald1(sv, surf, uv, 2)
  n = cross(su, sv)
  !PRINT *,' UV =',UV
  !PRINT *,'XYZ =',XYZ
  !PRINT *,'  N =',N

  mn(:) =  huge(1._fp)
  mx(:) = -huge(1._fp)
  
  do i = 1,bc%degr(1)+1
     dotn = dot_product(bc%coef(i,1:3,1) - xyz, n)
     !PRINT *,I,DOTN
     mn(1) = min(mn(1), dotn)
     mx(1) = max(mx(1), dotn)
  end do

  do j = 1,bs%degr(2)+1
     do i = 1,bs%degr(1)+1
        !PRINT *,I,J,DOTN
        dotn = dot_product(bs%coef(i,j,1:3) - xyz, n)
        mn(2) = min(mn(2), dotn)
        mx(2) = max(mx(2), dotn)
     end do
  end do

  !PRINT *,'MAX(MIN) - MIN(MAX) =', maxval(mn) - minval(mx)
  
  separable = ( maxval(mn) < minval(mx) + EPSxyz )
  
end subroutine curve_surface_separation_at_tangential_contact
