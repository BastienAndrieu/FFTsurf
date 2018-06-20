subroutine intersect_border_surface( &
     root_s, &
     root_c, &
     region, &
     icurv, &
     ivar, &
     ival, &
     uvxyz, &
     nuvxyz, &
     stat_degeneracy )
  use mod_util
  use mod_math
  use mod_chebyshev2
  use mod_polynomial
  use mod_diffgeom2
  use mod_regiontree
  use mod_tolerances
  implicit none
  LOGICAL, PARAMETER :: DEBUG = .false.
  integer, parameter                        :: npts_init = 10
  type(ptr_surface),          intent(in)    :: root_s(2)
  type(ptr_region),           intent(inout) :: region(2)
  integer,                    intent(in)    :: icurv, ivar, ival
  real(kind=fp), allocatable, intent(inout) :: uvxyz(:,:)
  integer,                    intent(inout) :: nuvxyz
  integer,                    intent(inout) :: stat_degeneracy
  !type(type_polynomial)                     :: regc
  type(type_curve),           intent(in)    :: root_c
  type(type_region)                         :: region_c
  type(type_region)                         :: region_s
  real(kind=fp), allocatable                :: tuvxyz(:,:)
  integer                                   :: ntuvxyz
  real(kind=fp)                             :: uv(2,2)
  integer                                   :: isurf, ipt, jpt

  !IF ( DEBUG ) THEN
  !PRINT *,''; PRINT *,''; PRINT *,''; PRINT *,''
  !PRINT *,'   ICURV, IVAR, IVAL =', ICURV, IVAR, IVAL
  !END IF

  isurf = 1 + mod(icurv,2)

  !! convert surface border to parametric curve
  !! <-------------------------------------------+
  !call chgvar2( &                               !
  !     root_s(icurv)%ptr%x, &                    !
  !     regc, &                                   !-- METTRE DANS INTERSECT_SIMPLE_SURFACES
  !     region(icurv)%ptr%uvbox([1,3]), &         !
  !     region(icurv)%ptr%uvbox([2,4]) )          !
  !! <-------------------------------------------+
  !call bivar2univar( &
  !     regc, &!root_s(icurv)%ptr%x, &
  !     root_c%x, &
  !     ivar, &
  !     ival )       
  !all economize1( root_c%x, EPSmath )
  !all compute_deriv1( root_c )
  !all compute_deriv2( root_c )

  IF ( DEBUG ) THEN
     call write_polynomial( root_c%x, 'dev_intersection_simple_surface/root_c_x.cheb' )
     call write_polynomial( root_c%xt, 'dev_intersection_simple_surface/root_c_xt.cheb' )
     call write_polynomial( root_s(isurf)%ptr%x, 'dev_intersection_simple_surface/root_s_x.cheb' )
     call write_polynomial( root_s(isurf)%ptr%xu, 'dev_intersection_simple_surface/root_s_xu.cheb' )
     call write_polynomial( root_s(isurf)%ptr%xv, 'dev_intersection_simple_surface/root_s_xv.cheb' )
  END IF


  ! initialize curve region tree
  call init_region( &
       region_c, &
       1, &
       [ -1._fp, 1._fp ] ) 

  ! initialize surface region tree
  call init_region( &
       region_s, &
       2, &
       region(isurf)%ptr%uvbox ) 

  ! compute curve Bezier control points
  allocate( region_c%poly )
  call cheb2bern_poly( &
       root_c%x, &
       region_c%poly )

  ! copy surface Bezier control points
  region_s%poly => region(isurf)%ptr%poly

  IF ( DEBUG ) THEN
     CALL WRITE_POLYNOMIAL( &
          REGION_C%POLY, &
          'dev_intersection_simple_surface/root_c_bezier.bern' )
     CALL WRITE_POLYNOMIAL( &
          REGION_S%POLY, &
          'dev_intersection_simple_surface/root_s_bezier.bern' )
  END IF






  allocate( tuvxyz(6,npts_init) )
  ntuvxyz = 0
  stat_degeneracy = 0
  call intersect_curve_surface( &
       root_c, &
       root_s(isurf)%ptr, &
       region_c, &
       region_s, &
       tuvxyz, &
       ntuvxyz, &
       stat_degeneracy )

  IF ( DEBUG ) THEN
     open( &
          unit=13, &
          file='dev_intersection_simple_surface/tuv_xyz.dat', &
          action='write' )
     if ( ntuvxyz < 1 ) then
        write (13,*) ''
     else
        do ipt = 1,ntuvxyz
           write (13,*) tuvxyz(:,ipt)
        end do
     end if
     close(13)

     call export_region_tree( &
          region_s, &
          'dev_intersection_simple_surface/tree_s.dat' )
     call export_region_tree( &
          region_c, &
          'dev_intersection_simple_surface/tree_c.dat' )
  END IF

  IF ( STAT_DEGENERACY > 10 ) THEN
     PRINT *,'STAT_DEGENERACY =',STAT_DEGENERACY
     RETURN
     !STOP '********************'
  END IF

  !PRINT  *,'';PRINT  *,'';PRINT  *,'';PRINT  *,'';PRINT  *,'';PRINT  *,'';
  !PRINT *,'     BEFORE :',NUVXYZ,' POINTS'
  !IF ( NUVXYZ > 0 ) CALL PRINT_MAT( TRANSPOSE(UVXYZ(:,1:NUVXYZ)) )
  outer : do ipt = 1,ntuvxyz
     do jpt = 1,nuvxyz
        if ( sum( ( tuvxyz(4:6,ipt) - uvxyz(5:7,jpt) )**2 ) < EPSxyzsqr ) cycle outer
     end do

     uv(ivar,icurv) = region(icurv)%ptr%uvbox(2*(ivar-1)+ival)
     uv(1+mod(ivar,2),icurv) = n1p12ab( &
          tuvxyz(1,ipt), &
          region(icurv)%ptr%uvbox(2*(1+mod(ivar,2))-1), &
          region(icurv)%ptr%uvbox(2*(1+mod(ivar,2))) )
     uv(:,isurf) = tuvxyz(2:3,ipt)
     !PRINT *,uv,tuvxyz(4:6,ipt)
     !PRINT *,'       XYZ  =',TUVXYZ(4:6,IPT)
     !PRINT *,'       T    =',TUVXYZ(1,IPT)
     !PRINT *,'       UVBOX=',REGION(ICURV)%PTR%UVBOX(2*IVAR+[-1,0])
     !PRINT *,'       UV   =',UV(:,ICURV)

     !CALL WRITE_POLYNOMIAL( &
     !     REGION_C%POLY, '&
     !     dev_intersection_simple_surface/root_c_bezier.bern' )

     call append_vector( &
          [ uv(:,1), uv(:,2), tuvxyz(4:6,ipt) ], &
          7, &
          uvxyz, &
          nuvxyz )

  end do outer
  !PRINT *,'     AFTER  :',NUVXYZ,' POINTS'
  !IF ( NUVXYZ > 0 ) CALL PRINT_MAT( TRANSPOSE(UVXYZ(:,1:NUVXYZ)) )

  call free_polynomial( region_c%poly )
  deallocate( region_c%poly )

  nullify( region_s%poly )

  call free_region_tree( region_s )
  call free_region_tree( region_c )

  !call free_polynomial( root_c%x )
  !call free_polynomial( root_c%xt )
  !call free_polynomial( root_c%xtt )

end subroutine intersect_border_surface
