subroutine check_curve_surface_intersection_point( &
   curv, &
   surf, &
   t, &
   uv, &
   stat )
  use mod_math
  use mod_diffgeom
  use mod_tolerances
  implicit none
  LOGICAL, PARAMETER :: DEBUG = ( GLOBALDEBUG .AND. .false. )
  type(type_curve),   intent(in)  :: curv
  type(type_surface), intent(in)  :: surf
  real(kind=fp),      intent(in)  :: t
  real(kind=fp),      intent(in)  :: uv(2)
  integer,            intent(out) :: stat
  real(kind=fp)                   :: dxyz_dt(3), dxyz_duv(3,2), n(3)
  real(kind=fp)                   :: EFG(3), duv_dt(2)
  logical                         :: singular
  real(kind=fp)                   :: d2xyz_dt2(3), d2xyz_duv2(3,3), y(3)
  integer                         :: ivar, jvar

  ! curve tangent
  call evald1( &
       dxyz_dt, &
       curv, &
       t )
  
  ! surface tangents
  do ivar = 1,2
     call evald1( &
          dxyz_duv(:,ivar), &
          surf, &
          uv, &
          ivar)
  end do

  ! surface normal
  n = cross(dxyz_duv(:,1), dxyz_duv(:,2))
  n = n / norm2(n)
  
  IF ( DEBUG ) PRINT *,'TUV =',T,UV
  IF ( DEBUG ) PRINT *,'CT.N =',DOT_PRODUCT(DXYZ_DT, N)

  if ( abs(dot_product(dxyz_dt, n)) > EPSxyz ) then
     stat = 0
  else
     do ivar = 1,2
        do jvar = ivar,2
           EFG(ivar+jvar-1) = dot_product(dxyz_duv(:,ivar), dxyz_duv(:,jvar))
        end do
     end do
     call solve_2x2( &
          duv_dt, &
          EFG(1), &
          EFG(2), &
          EFG(2), &
          EFG(3), &
          [ dot_product(dxyz_dt, dxyz_duv(:,1)), &
          dot_product(dxyz_dt, dxyz_duv(:,2)) ], &
          singular )
     
     ! curve 2nd derivative
     call evald2( &
          d2xyz_dt2, &
          curv, &
          t )

     ! surface 2nd derivatives
     do ivar = 1,3
        call evald2( &
             d2xyz_duv2(:,ivar), &
             surf, &
             uv, &
             ivar)
     end do

     y = d2xyz_dt2 - ( &
          d2xyz_duv2(:,1)*duv_dt(1)**2 + &
          d2xyz_duv2(:,2)*duv_dt(1)*duv_dt(2) + &
          d2xyz_duv2(:,3)*duv_dt(2)**2 )
     IF ( DEBUG ) PRINT *,'Y.N =',DOT_PRODUCT(Y, N)
     if ( dot_product(y, n) > EPSxyz ) then
        stat = 1
     else
        stat = 2
     end if
  end if

end subroutine check_curve_surface_intersection_point
