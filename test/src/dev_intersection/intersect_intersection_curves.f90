subroutine intersect_intersection_curves( &
     curv )
  USE MOD_UTIL
  use mod_math
  use mod_types_intersection
  implicit none
  type(ptr_intersection_curve), intent(inout) :: curv(2)
  logical                                     :: may_intersect
  real(kind=fp)                               :: uvbox(2,2)
  integer                                     :: numsurf(2)
  type(type_matrix)                           :: polylineuv(2)
  integer                                     :: np(2)
  integer, allocatable                        :: ipls(:,:)
  real(kind=fp), allocatable                  :: lambda(:,:)
  integer                                     :: npts
  type(ptr_surface)                           :: surf(3)
  real(kind=fp)                               :: uv(2,2), xyz(3), w
  real(kind=fp), dimension(6)                 :: lowerb, upperb
  integer                                     :: stat
  integer                                     :: isurf, jsurf, ivar, icurv, ipt
  INTEGER :: FID, I

  may_intersect = .false.
  jloop : do jsurf = 1,2
     iloop : do isurf = 1,2
        if ( associated(curv(1)%ptr%surf(isurf)%ptr, curv(2)%ptr%surf(jsurf)%ptr) ) then
           ! the two intersection curve share one incident surface
           ! compute the intersection of their bounding boxes in the uv-space of that surface
           do ivar = 1,2
              uvbox(:,ivar) = [ &
                   max(curv(1)%ptr%uvbox(1,ivar,isurf), curv(2)%ptr%uvbox(1,ivar,jsurf)), &
                   min(curv(1)%ptr%uvbox(2,ivar,isurf), curv(2)%ptr%uvbox(2,ivar,jsurf)) ]
              if ( uvbox(2,ivar) - uvbox(1,ivar) < epsilon(1._fp) ) cycle iloop ! empty intersection
           end do
           ! if we got this far, the two curves may intersect
           numsurf = [isurf, jsurf]
           may_intersect = .true.
           exit jloop
        end if
     end do iloop
  end do jloop

  if (.not.may_intersect) return ! the curves cannot intersect

  do icurv = 1,2
     np(icurv) = curv(icurv)%ptr%polyline%np
     allocate(polylineuv(icurv)%mat(2,np(icurv)))
     polylineuv(icurv)%mat(1:2,1:np(icurv)) = curv(icurv)%ptr%polyline%uv(1:2,numsurf(icurv),1:np(icurv))
  end do

  ! intersect the polylines in the uv-space of the shared incident surface
  npts = 0
  call intersect_2Dpolylines( &
       polylineuv, &
       [1,1], &
       np, &
       npts, &
       ipls, &
       lambda )
  deallocate(polylineuv(1)%mat, polylineuv(2)%mat)

  if ( npts < 1 ) return ! the polylines do not intersect, we assume the curves do not either

  do ipt = 1,npts
     !! refine all intersection points using Newton algorithm
     surf(1)%ptr => curv(1)%ptr%surf(numsurf(1))%ptr
     do icurv = 1,2
        surf(1+icurv)%ptr => curv(icurv)%ptr%surf(1+mod(numsurf(icurv),2))%ptr
     end do

     ! set initial iterate
     uv(:,1) = &
          (1._fp - lambda(1,ipt)) * curv(1)%ptr%polyline%uv(:,numsurf(1),ipls(1,ipt)) + &
          lambda(1,ipt) * curv(1)%ptr%polyline%uv(:,numsurf(1),ipls(1,ipt)+1)
     do icurv = 1,2
        isurf = 1 + mod(numsurf(icurv),2)
        uv(:,1+icurv) = &
             (1._fp - lambda(icurv,ipt)) * curv(icurv)%ptr%polyline%uv(:,isurf,ipls(icurv,ipt)) + &
             lambda(icurv,ipt) * curv(icurv)%ptr%polyline%uv(:,isurf,ipls(icurv,ipt)+1)
     end do

     ! lower and upper bounds of feasible domain
     lowerb(1:2) = uvbox(1,1:2)
     upperb(1:2) = uvbox(2,1:2)
     do icurv = 1,2
        isurf = 1 + mod(numsurf(icurv),2)
        lowerb(2*(icurv+1)-1:2*(icurv+1)) = curv(icurv)%ptr%uvbox(1,1:2,isurf)
        upperb(2*(icurv+1)-1:2*(icurv+1)) = curv(icurv)%ptr%uvbox(2,1:2,isurf)
     end do

     call newton_three_surfaces( &
          surf, &
          lowerb, &
          upperb, &
          stat, &
          uv, &
          xyz )

     IF ( STAT == 0 ) THEN
        PRINT *,'NEWTON 3 SURFACES CONVERGED'
        PRINT *,' UV =',UV
        PRINT *,'XYZ =',XYZ
     ELSE
        CALL WRITE_POLYNOMIAL( SURF(1)%PTR%X, 'dev_intersection/debug3si_surf1.cheb' )
        CALL WRITE_POLYNOMIAL( SURF(2)%PTR%X, 'dev_intersection/debug3si_surf2.cheb' )
        CALL WRITE_POLYNOMIAL( SURF(3)%PTR%X, 'dev_intersection/debug3si_surf3.cheb' )
        CALL GET_FREE_UNIT(FID)
        OPEN(UNIT=FID, FILE='dev_intersection/debug3si.dat', ACTION='WRITE')
        WRITE(FID,*) LOWERB
        WRITE(FID,*) UPPERB
        DO ICURV = 1,2
           WRITE(FID,*) NUMSURF(ICURV)
           WRITE(FID,*) CURV(ICURV)%PTR%POLYLINE%NP
           DO I = 1,CURV(ICURV)%PTR%POLYLINE%NP
              WRITE(FID,*) CURV(ICURV)%PTR%POLYLINE%UV(:,:,I)
           END DO
        END DO
        CLOSE(FID)
        STOP '----> DEBUG 3SI'
     END IF
     
     
     do icurv = 1,2
        w = dot_product(xyz, curv(icurv)%ptr%param_vector)
        
     end do

  end do

end subroutine intersect_intersection_curves
