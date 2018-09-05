program diffgeom
  use mod_math
  use mod_polynomial
  use mod_diffgeom
  use mod_intersection

  implicit none

  logical, parameter :: ECONOMIZE = .false.
  type(ptr_surface)  :: surf(2)
  real(kind=fp)      :: uv(2,2)
  real(kind=fp)      :: duv_ds(2,2,2) ! u/v, #branch, #surf
  real(kind=fp)      :: dxyz_ds(3,2)
  real(kind=fp)      :: curvature(2)
  integer            :: stat, i
  
  allocate(surf(1)%ptr, surf(2)%ptr)
  
  call read_polynomial( &
       surf(1)%ptr%x, &
       'diffgeom/C1.cheb', &
       !'Jouke/propergol/C_005.cheb', &
       nvar=2, &
       base=1 )
  call read_polynomial( &
       surf(2)%ptr%x, &
       'diffgeom/C2.cheb', &
       !'Jouke/propergol/C_006.cheb', &
       nvar=2, &
       base=1 )

  do i = 1,2
     if ( ECONOMIZE ) call economize2(surf(i)%ptr%x, EPSmath)
     call compute_deriv1(surf(i)%ptr)
     call compute_deriv2(surf(i)%ptr)
  end do

  !call random_number(r)
  !r = 2._fp*r - 1._fp
  !call myrand(r)
  !print *,r
  !uv(:,1) = [r, 1._fp]
  !uv(:,2) = [-r, 1._fp]
  uv(:,:) = 0._fp

  call diffgeom_intersection( &
       surf, &
       uv, &
       duv_ds, &
       dxyz_ds, &
       stat, &
       curvature )
  PRINT *,'STAT =',STAT
  PRINT *,'DUV_DS ='
  DO I = 1,max(1,STAT)
     PRINT *,DUV_DS(:,I,:)
  END DO
  PRINT *,'DXYZ_DS ='
  DO I = 1,max(1,STAT)
     PRINT *,DXYZ_DS(:,I)
  END DO
  if ( stat == 0 ) print *,'curvature =',curvature(1)
  !PRINT *,''; PRINT *,''; PRINT *,''

contains

   subroutine myrand( u )
    implicit none
    real(kind=fp), intent(out) :: u
    integer                    :: time(8)
    real(kind=fp)              :: x
    call date_and_time(values=time)
    x = dble( time(4) * (360000*time(5) + 6000*time(6) + 100*time(7) + time(8)) )
    u = sin( exp(1._fp)*x )
  end subroutine myrand
  
end program diffgeom
