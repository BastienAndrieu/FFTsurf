module mod_intersection
implicit none
LOGICAL, PARAMETER :: GLOBALDEBUG = .false.
contains



subroutine inherit_points( &
     region, &
     coords, &
     tol, &
     npts )
  use mod_util
  use mod_math
  use mod_regiontree
  use mod_tolerances
  implicit none
  LOGICAL, PARAMETER :: DEBUG = ( GLOBALDEBUG .AND. .false. )
  type(type_region), intent(inout) :: region
  integer,           intent(in)    :: npts
  real(kind=fp),     intent(in)    :: coords(region%dim,npts)
  real(kind=fp),     intent(in)    :: tol(npts)
  logical, allocatable             :: mask(:)
  integer                          :: idim, ipt, jpt

  if ( .not.associated(region%parent) ) return
  if ( region%parent%npts < 1 ) return
  
  allocate( mask(region%parent%npts) )
  mask(:) = .true.
  outer : do jpt = 1,region%parent%npts ! <-------------+
     ipt = region%parent%ipts(jpt)                      !
     IF ( IPT < 0 .OR. IPT > NPTS ) THEN
        PRINT *,'*** REGION%PARENT%IPTS =',REGION%PARENT%IPTS(1:REGION%PARENT%NPTS)
     END IF
     do idim = 1,region%dim ! <--------------------+    !
        if ( .not.is_in_closed_interval( &         !    !
             coords(idim,ipt), &                   !    !
             region%uvbox(2*idim-1), &             !    !
             region%uvbox(2*idim), &               !    !
             !tolerance=EPSregion) ) then ! <--+    !    !
             tolerance=max(EPSuv, tol(ipt))) ) then
           mask(jpt) = .false.                !    !    !
           cycle outer                        !    !    !
        end if ! <----------------------------+    !    !
     end do ! <------------------------------------+    !
  end do outer ! <--------------------------------------+

  IF ( DEBUG ) THEN
     PRINT *,'   REGION%UVBOX =',REGION%UVBOX
     DO IPT = 1,region%parent%npts
        PRINT *,REGION%PARENT%IPTS(IPT), COORDS(:,REGION%PARENT%IPTS(IPT)), TOL(REGION%PARENT%IPTS(IPT)), MASK(IPT)
     END DO
     PRINT *,'   INHERITS',pack(region%parent%ipts(1:region%parent%npts), mask)
     !DO IPT = 1,region%parent%npts
     !   IF ( MASK(IPT) ) PRINT *,COORDS(:,REGION%PARENT%IPTS(IPT))
     !END DO
  END IF
  call append_n( &
       region%ipts, &
       region%npts, &
       pack(region%parent%ipts(1:region%parent%npts), mask), &
       count(mask), &
       unique=.true. )

  deallocate( mask )

end subroutine inherit_points



subroutine add_intersection_point( &
    uv, &
    xyz, &
    surf, &
    nsurf, &
    interdata, &
    id ) 
  use mod_math
  use mod_diffgeom
  use mod_types_intersection
  use mod_tolerances
  implicit none
  integer, parameter                          :: PARAM_xtra_np = 10
  integer,                      intent(in)    :: nsurf
  real(kind=fp),                intent(in)    :: uv(2,nsurf)
  real(kind=fp),                intent(in)    :: xyz(3)
  type(ptr_surface),            intent(in)    :: surf(nsurf)
  type(type_intersection_data), intent(inout) :: interdata
  integer,                      intent(out)   :: id
  type(type_point_on_surface), pointer        :: pos
  logical                                     :: newsurf(nsurf)
  type(type_intersection_data)                :: tmp
  integer                                     :: isurf, ipos

  ! compare the point to be inserted with already collected points
  do id = 1,interdata%np ! <-------------------------------------------------------+
     if ( sum((xyz - interdata%points(id)%xyz)**2) < EPSxyzsqr ) then ! <-------+  !
        pos => interdata%points(id)%pos                                         !  !
        newsurf(1:nsurf) = .true.                                               !  !
        do ipos = 1,interdata%points(id)%npos ! <----------------------------+  !  !
           if ( all(.not.newsurf) ) exit                                     !  !  !
           do isurf = 1,nsurf ! <-----------------------------------------+  !  !  !
              if ( .not.newsurf(isurf) ) cycle                            !  !  !  !
              if ( associated(pos%surf, surf(isurf)%ptr) ) then ! <----+  !  !  !  !
                 newsurf(isurf)= .false.                               !  !  !  !  !
                 exit                                                  !  !  !  !  !
              end if ! <-----------------------------------------------+  !  !  !  !
           end do ! <-----------------------------------------------------+  !  !  !
           if ( ipos < interdata%points(id)%npos ) pos => pos%next           !  !  !
        end do ! <-----------------------------------------------------------+  !  !
        !                                                                       !  !
        do isurf = 1,nsurf ! <-----------------------------------------------+  !  !
           if ( newsurf(isurf) ) then ! <---------------------------------+  !  !  !
              interdata%points(id)%npos = interdata%points(id)%npos + 1   !  !  !  !
              allocate(pos%next)                                          !  !  !  !
              pos => pos%next                                             !  !  !  !
              pos%uv = uv(:,isurf)                                        !  !  !  !
              pos%surf => surf(isurf)%ptr                                 !  !  !  !
           end if ! <-----------------------------------------------------+  !  !  !
        end do ! <-----------------------------------------------------------+  !  !
        nullify(pos)                                                            !  !
        return                                                                  !  !
        !                                                                       !  !
     end if ! <-----------------------------------------------------------------+  !
  end do ! <-----------------------------------------------------------------------+

  
  if ( .not.allocated(interdata%points) ) allocate(interdata%points(PARAM_xtra_np) )

  ! the point to be inserted is not a duplicate, we insert it

  ! reallocate if necessary
  if ( interdata%np + 1 > size(interdata%points) ) then
     allocate(tmp%points(interdata%np))
     call transfer_intersection_points( &
          from=interdata, &
          to=tmp )

     allocate(interdata%points(tmp%np + PARAM_xtra_np))
     call transfer_intersection_points( &
          from=tmp, &
          to=interdata )
  end if
  
  ! insert new point
  interdata%np = interdata%np + 1
  id = interdata%np
  interdata%points(id)%xyz = xyz
  interdata%points(id)%npos = nsurf
  allocate(interdata%points(id)%pos)
  pos => interdata%points(id)%pos
  do isurf = 1,nsurf ! <-----------------------------------------------+
     pos%uv = uv(:,isurf)                                              !
     pos%surf => surf(isurf)%ptr                                       !
     if ( isurf < nsurf ) allocate(pos%next)                           !
     pos => pos%next                                                   !
  end do ! <-----------------------------------------------------------+

  nullify(pos)

end subroutine add_intersection_point



! ********* INCOMPLET *********
subroutine classify_border_surface_intersection_point( &
     surf, &
     region, &
     uv, &
     npts, &
     stat_point ) 
  use mod_math
  use mod_diffgeom
  use mod_regiontree
  !stat_point = -1 : exiting
  !              0 : isolated
  !              1 : entering
  implicit none
  LOGICAL, PARAMETER :: DEBUG = ( GLOBALDEBUG .AND. .false. )
  type(ptr_surface), intent(in)  :: surf(2)
  type(ptr_region),  intent(in)  :: region(2)
  integer,           intent(in)  :: npts
  real(kind=fp),     intent(in)  :: uv(2,2,npts)
  integer,           intent(out) :: stat_point(npts)
  real(kind=fp)                  :: duv_ds(2,2,2), dxyz_ds(3,2)
  integer                        :: stat_contactpoint
  integer                        :: stat(2), sumstat
  integer                        :: isurf, ipt

  do ipt = 1,npts ! <---------------------------------------------+
     call diffgeom_intersection( &                                !
          surf, &                                                 !
          uv(:,:,ipt), &                                          !
          duv_ds, &                                               !
          dxyz_ds, &                                              !
          stat_contactpoint )                                     !
     !call diffgeom_intersection_curve( &                          !
     !     surf, &                                                 !
     !     uv(:,:,ipt), &                                          !
     !     duv_ds, &                                               !
     !     dxyz_ds, &                                              !
     !     stat_contactpoint )                                     !
     IF (DEBUG ) PRINT *,'    UV =',uv(:,:,ipt)
     IF (DEBUG ) PRINT *,'DUV_DS =',duv_ds(:,1,:)
     !                                                            !
     if ( stat_contactpoint > 1 ) then ! <---------------------+  !
        PRINT *,'classify_border_surface_intersection_point : &
             &stat_contactpoint = ',stat_contactpoint          !  !
        STOP                                                   !  !
        if ( stat_contactpoint == 2 ) then ! <-------------+   !  !
           ! two distinct tangential directions            !   !  !
           ! (...)                                         !   !  !
        else ! --------------------------------------------+   !  !                                               
           ! isolated point                                !   !  !
           ! (...)                                         !   !  !
        end if ! <-----------------------------------------+   !  !
     end if ! <------------------------------------------------+  !
     !                                                            !
     do isurf = 1,2                                               !
        call classify_endpoint( &                                 !
             uv(:,isurf,ipt), &                                   !
             duv_ds(:,1,isurf), &                                 !
             region(isurf)%ptr%uvbox, &                           !
             stat(isurf) )                                        !
     end do                                                       !
     IF (DEBUG ) PRINT *,'STATSURF =',STAT
     !                                                            !
     sumstat = sum(stat)                                          !
     if ( any(stat == 2) ) then ! <-------+                       !
        stat_point(ipt) = 2
        !STOP 'AMBIGUOUS POINT'
     elseif ( sumstat == 0 ) then ! <-----+                       !
        stat_point(ipt) = 0               !                       !
     elseif ( sumstat < 0 ) then ! -------+                       !
        stat_point(ipt) = -1              !                       !
     else ! ------------------------------+                       !
        stat_point(ipt) = 1               !                       !
     end if ! <---------------------------+                       !
  end do ! <------------------------------------------------------+

  IF ( ANY(stat_point == 2) ) STOP 'AMBIGUOUS POINT'

end subroutine classify_border_surface_intersection_point



subroutine characterize_tangential_intersection_point( &
     surf, &
     uv, &
     dxyz_duv, &
     n, &
     stat, &
     dxyz_ds )
  use mod_math
  use mod_diffgeom
  ! NOT USED ANYMORE
  ! ( see "Shape interrogation for computer aided design and manufacturing", &
  ! Patrikalakis et al. (2009), pp.172-175 )
  ! Returns stat = 1 : single tangential intersection curve
  !              = 2 : two tangential intersection branches
  !              = 3 : isolated tangential contact point
  !              = 4 : higher-order contact point
  implicit none
  real(kind=fp), parameter       :: EPS = EPSmath
  type(ptr_surface), intent(in)  :: surf(2)
  real(kind=fp),     intent(in)  :: uv(2,2)
  real(kind=fp),     intent(in)  :: dxyz_duv(3,2,2)
  real(kind=fp),     intent(in)  :: n(3,2)
  integer,           intent(out) :: stat
  real(kind=fp),     intent(out) :: dxyz_ds(3,2)
  real(kind=fp)                  :: s2(3)
  real(kind=fp)                  :: aux(3,2), A(2,2), LMN(3,2), lambda, B(3)
  real(kind=fp)                  :: discr, w(2)
  integer                        :: isurf, ivar

  ! compute A
  do ivar = 1,2
     aux(:,ivar) = cross( n(:,1), dxyz_duv(:,1+mod(ivar,2),2) )
  end do
  aux(:,1) = -aux(:,1)
  A = matmul( transpose(aux), dxyz_duv(:,:,1) )

  ! compute second fundamental forms
  do isurf = 1,2 ! <-----------------------------------------+
     do ivar = 1,3 ! <------------------------------------+  !
        call evald2( &                                    !  !
             s2, &                                        !  !
             surf(isurf)%ptr, &                           !  !
             uv(:,isurf), &                               !  !
             ivar )                                       !  !
        LMN(ivar,isurf) = dot_product( n(:,isurf), s2 )   !  !
     end do ! <-------------------------------------------+  !
  end do ! <-------------------------------------------------+

  ! compute B = [ b11 ; b12 ; b22 ]
  B(1) = LMN(1,2)*A(1,1)**2 + &
       LMN(2,2)*2*A(1,1)*A(2,1) + &
       LMN(3,2)*A(2,1)**2

  B(2) = LMN(1,2)*A(1,1)*A(1,2) + &
       LMN(2,2)*(A(1,1)*A(2,2) + A(1,2)*A(2,1)) + &
       LMN(3,2)*A(2,1)*A(2,2)

  B(3) = LMN(1,2)*A(1,2)**2 + &
       LMN(2,2)*2*A(1,2)*A(2,2) + &
       LMN(3,2)*A(2,2)**2

  lambda = 1._fp / &
       ( sqrt( product(sum(n**2, 1)) ) * dot_product( n(:,1), n(:,2) ) )

  B = LMN(:,1) - lambda * B

  ! analyze B
  discr = B(2)**2 - B(1)*B(3)

  if ( all( abs(B) < EPS ) ) then ! <--------------------------------+
     ! higher order contact point (dxyz_ds cannot be determined)     !
     stat = 4                                                        !
     return                                                          !
     !                                                               !
  elseif ( abs( discr ) < EPS ) then ! ------------------------------+
     ! tangential intersection curve (one single dxyz_ds)            !
     stat = 1                                                        !
     if ( abs(B(1)) > abs(B(3)) ) then ! <------------------------+  !
        w(1) = -B(2) / B(1)                                       !  !
        dxyz_ds(:,1) = w(1)*dxyz_duv(:,1,1) + dxyz_duv(:,2,1)     !  !
     else ! ------------------------------------------------------+  !
        w(1) = -B(2) / B(3)                                       !  !
        dxyz_ds(:,1) = dxyz_duv(:,1,1) + w(1)*dxyz_duv(:,2,1)     !  !
     end if ! <---------------------------------------------------+  !
     !                                                               !
  elseif ( discr < 0._fp ) then ! -----------------------------------+
     ! isolated tangential contact point (dxyz_ds undefined)         !
     stat = 3                                                        !
     return                                                          !
     !                                                               !
  else ! ------------------------------------------------------------+
     ! branch point (two disctinct dxyz_ds)                          !
     stat = 2                                                        !
     w = real( [-1,1], kind=fp ) * sqrt( discr )                     !
     if ( abs(B(1)) > abs(B(3)) ) then ! <------------------------+  !
        w = (w - B(2)) / B(1)                                     !  !
        dxyz_ds = outer_product( dxyz_duv(:,1,1), w ) + &         !  !
             spread( dxyz_duv(:,2,1), dim=2, ncopies=2 )          !  !
     else ! ------------------------------------------------------+  !
        w = (w - B(2)) / B(3)                                     !  !
        dxyz_ds = spread( dxyz_duv(:,1,1), dim=2, ncopies=2 ) + & !  !
             outer_product( dxyz_duv(:,2,1), w )                  !  !
     end if ! <---------------------------------------------------+  !
     !                                                               !
  end if ! <---------------------------------------------------------+

  ! unitize dxyz_ds (the number of distinct tangential directions
  ! is equal to stat)
  do ivar = 1,stat
     dxyz_ds(:,ivar) = dxyz_ds(:,ivar) / norm2( dxyz_ds(:,ivar) )
  end do

end subroutine characterize_tangential_intersection_point



subroutine newton_curve_surface_dynshift( &
     curv, &
     surf, &
     lowerb, &
     upperb, &
     stat, &
     tuv, &
     toltuv, &
     xyz )
  use mod_math
  use mod_linalg
  use mod_diffgeom
  use mod_tolerances
  implicit none
  LOGICAL, PARAMETER :: DEBUG = ( GLOBALDEBUG .AND. .false. )
  integer,       parameter          :: itmax = 20!2 + ceiling(-log10(EPSuv))
  real(kind=fp), parameter          :: shiftd = real(1.d-1, kind=fp)
  type(type_curve),   intent(in)    :: curv
  type(type_surface), intent(in)    :: surf
  real(kind=fp),      intent(in)    :: lowerb(3)
  real(kind=fp),      intent(in)    :: upperb(3)
  integer,            intent(out)   :: stat
  real(kind=fp),      intent(inout) :: tuv(3)
  real(kind=fp),      intent(out)   :: toltuv
  real(kind=fp),      intent(out)   :: xyz(3)
  real(kind=fp)                     :: shiftv(3), sgn
  real(kind=fp), dimension(3)       :: xyz_c, xyz_s, r
  real(kind=fp), dimension(3)       :: ct, ctt, su, sv, suu, suv, svv, n
  real(kind=fp)                     :: resxyz
  real(kind=fp)                     :: grad(3), hess(3,3), dtuv(3), cond, errtuv
  integer                           :: it
  integer                           :: stat_refl
  
  call evald2(ctt, curv, tuv(1))

  stat = 1
  errtuv = 0._fp
  cond = 1._fp
  do it = 1,itmax
     call evald1(ct, curv, tuv(1))
     call evald1(su, surf, tuv(2:3), 1)
     call evald1(sv, surf, tuv(2:3), 2)
     n = cross(su, sv)
     
     sgn = sign(1._fp, dot_product(n, ctt))
     shiftv = sgn * shiftd * n
     
     !! compute residual
     call eval(xyz_c, curv, tuv(1)  ) ! curve's position vector
     call eval(xyz_s, surf, tuv(2:3)) ! surface's position vector
     r = xyz_s - xyz_c
     resxyz = sum(r**2)
     r = r - shiftv
     IF ( DEBUG ) PRINT *,SQRT(RESXYZ), SQRT(ERRTUV), EPSFP*COND
     
     !! compute gradient of distance
     grad(1) = -dot_product(ct, r)
     grad(2) =  dot_product(su, r)
     grad(3) =  dot_product(sv, r)

     !! compute Hessian matrix
     if ( it > 1 ) call evald2(ctt, curv, tuv(1))
     call evald2(suu, surf, tuv(2:3), 1)
     call evald2(suv, surf, tuv(2:3), 2)
     call evald2(svv, surf, tuv(2:3), 3)

     hess(1,1) = -dot_product(ctt,r) + dot_product(ct,ct)
     hess(1,2) = -dot_product(ct,su)
     hess(1,3) = -dot_product(ct,sv)
     hess(2,2) =  dot_product(suu,r) + dot_product(su,su)
     hess(2,3) =  dot_product(suv,r) + dot_product(su,sv)
     hess(3,3) =  dot_product(svv,r) + dot_product(sv,sv)
     hess(2,1) = hess(1,2) ! ...
     hess(3,1) = hess(1,3) ! symmetry
     hess(3,2) = hess(2,3) ! ...

     !! solve for Newton step
     call linsolve_svd( &
          dtuv, &
          hess, &
          -grad, &
          3, &
          3, &
          1, &
          cond )

     errtuv = sum(dtuv**2)

     ! correct Newton step to keep the iterate inside feasible region
     call nd_box_reflexions( &
          tuv, &
          lowerb, &
          upperb, &
          dtuv, &
          3, &
          stat_refl )
     if ( stat_refl > 0 ) return
     
     ! update solution
     tuv = tuv + dtuv

     !! termination criterion
     if ( errtuv < max(EPSuvsqr, EPSfpsqr*cond**2) ) then
        if ( resxyz < EPSxyzsqr ) then
           if ( errtuv > EPSuvsqr ) then
              PRINT *,'newton_curve_surface_dynshift : /!\ toltuv > EPSuv'
           end if
           ! converged to a curve-surface intersection point
           stat = 0
           toltuv = max(EPSmath, sqrt(errtuv))
           call eval(xyz_c, curv, tuv(1)  )
           call eval(xyz_s, surf, tuv(2:3))
           xyz = 0.5_fp * (xyz_c + xyz_s)
           IF ( DEBUG ) PRINT *,SQRT(RESXYZ), SQRT(ERRTUV), EPSFP*COND
           IF ( DEBUG ) PRINT *,'CONVERGED, TUV =',TUV,', XYZ =',XYZ,', TOL =',TOLTUV
           !PAUSE
        else
           IF ( DEBUG ) PRINT *,'STAGNATION'
        end if
        return
     end if
  end do
end subroutine newton_curve_surface_dynshift



subroutine diffgeom_intersection_curve( &
    surf, &
    uv, &
    duv_ds, &
    dxyz_ds, &
    stat, &
    curvature )
  ! Computes tangent direction(s) and curvature of the intersection curve between
  ! two surfaces at a given intersection point.
  ! Returns: stat =-1 if one surface is singular at that point (normal = zero vector)
  !               = 0 if the point is on a single transveral intersection curve
  !               = 1 if the point is on a single tangential intersection curve
  !               = 2 if the point is a branch point (junction of two intersection curves)
  !               = 3 if the point is an isolated tangential contact point (degenerate curve)
  !               = 4 if the point is an high-order tangential contact point (possible degeneracy)
  use mod_math
  use mod_diffgeom
  ! ( see "Shape interrogation for computer aided design and manufacturing", &
  ! Patrikalakis et al. (2009), pp.166-175, 
  ! and "Tracing surface intersections with validated ODE system solver", Mukundan et al (2004)
  implicit none
  type(ptr_surface),       intent(in)  :: surf(2)
  real(kind=fp),           intent(in)  :: uv(2,2)
  real(kind=fp),           intent(out) :: duv_ds(2,2,2) ! u/v, #branch, #surf
  real(kind=fp),           intent(out) :: dxyz_ds(3,2)
  integer,                 intent(out) :: stat
  real(kind=fp), optional, intent(out) :: curvature(2)
  real(kind=fp)                        :: dxyz_duv(3,2,2), n(3,2), d2xyz_duv2(3)
  real(kind=fp)                        :: cosn, nsqr(2), kn(2), denom
  real(kind=fp), dimension(3,2)        :: EFG, LMN
  integer                              :: ndir
  logical                              :: singular
  integer                              :: isurf, ivar, jvar, i

  ! compute tangent and normal vectors to both surfaces at the intersection point
  do isurf = 1,2 ! <--------------------------------------------------+
     ! tangent vectors                                                !
     do ivar = 1,2 ! <-------------------------+                      !
        call evald1( &                         !                      !
             dxyz_duv(:,ivar,isurf), &         !                      !
             surf(isurf)%ptr, &                !                      !
             uv(:,isurf), &                    !                      !
             ivar )                            !                      !
     end do ! <--------------------------------+                      !
     ! (pseudo)normal vector                                          !
     n(:,isurf) = cross( dxyz_duv(:,1,isurf), dxyz_duv(:,2,isurf) )   !
     ! unit normal vector                                             !
     nsqr(isurf) = sum(n(:,isurf)**2)                                 !
     if ( nsqr(isurf) < 3._fp*EPSfpsqr ) then ! <----+                !
        ! singular surface                           !                !
        stat = -1                                    !                !
        return                                       !                !
     else ! -----------------------------------------+                !
        n(:,isurf) = n(:,isurf) / sqrt(nsqr(isurf))  !                !
     end if ! <--------------------------------------+                !
  end do ! <----------------------------------------------------------+

  ! cosine of the angle between the two normal vectors
  cosn = dot_product( n(:,1), n(:,2) )

  if ( abs(cosn) < 1._fp - EPSmath ) then ! <--------------------------------------------+
     ! the normals are not parallel, this is a transversal intersection point            !
     stat = 0                                                                            !
     ndir = 1 ! one single tangential direction                                          !
     dxyz_ds(:,1) = cross( n(:,1), n(:,2) )                                              !
     dxyz_ds(:,1) = dxyz_ds(:,1) / norm2( dxyz_ds(:,1) )                                 !
  else ! <-------------------------------------------------------------------------------+
     ! this is a tangential intersection point, the tangent direction(s) to              !
     ! the intersection banch(es) is (are) computed using 2nd derivatives                !
     call characterize_tangential_intersection_point( &                                  !
          surf, &                                                                        !
          uv, &                                                                          !
          dxyz_duv, &                                                                    !
          n, &                                                                           !
          stat, &                                                                        !
          dxyz_ds )                                                                      !
     ndir = stat ! 1 or 2 tangential direction(s)                                        !
  end if ! <-----------------------------------------------------------------------------+

  if ( stat > 2 ) then ! <-----------------------+
     ! isolated or high-order contact point      !
     ! => the tangential direction is undefined  !
     return                                      !
  end if ! <-------------------------------------+
  
  ! (at this stage, the number of intersection branches is equal to ndir)

  ! compute tangential direction(s) in the uv-space of each surface
  do isurf = 1,2 ! loop over surfaces ! <------------------------------------------------+
     ! 1st fundamental form coefficients                                                 !
     do ivar = 1,2 ! <-----------------------------------------------+                   !
        do jvar = ivar,2 ! <--------------------------------------+  !                   !
           EFG(ivar + jvar - 1,isurf) = dot_product( &            !  !                   !
                dxyz_duv(:,ivar,isurf), dxyz_duv(:,jvar,isurf) )  !  !                   !
        end do ! <------------------------------------------------+  !                   !
     end do ! <------------------------------------------------------+                   !
     !                                                                                   !
     do i = 1,ndir ! loop over tangential directions ! <-------------+                   !
        call solve_2x2( &                                            !                   !
             duv_ds(:,i,isurf), &                                    !                   !
             EFG(1,isurf), &                                         !                   !
             EFG(2,isurf), &                                         !                   !
             EFG(2,isurf), &                                         !                   !
             EFG(3,isurf), &                                         !                   !
             [dot_product(dxyz_ds(:,i), dxyz_duv(:,1,isurf)), &      !                   !
             dot_product(dxyz_ds(:,i), dxyz_duv(:,2,isurf))], &      !                   !
             singular )                                              !                   !
     end do ! <------------------------------------------------------+                   !
  end do ! <-----------------------------------------------------------------------------+

  if ( present(curvature) ) then ! <-----------------------------------------------------+
     !! compute curvature of the intersection branch(es) at the current point            !
     do isurf = 1,2 ! <--------------------------------------------------------------+   !
        ! 2nd fundamental form coefficients                                          !   !
        do ivar = 1,3 ! <-----------------------------------------------+            !   !
           call evald2( &                                               !            !   !
                d2xyz_duv2, &                                           !            !   !
                surf(isurf)%ptr, &                                      !            !   !
                uv(:,isurf), &                                          !            !   !
                ivar )                                                  !            !   !
           LMN(ivar,isurf) = dot_product(n(:,isurf), d2xyz_duv2)        !            !   !
        end do ! <------------------------------------------------------+            !   !
     end do ! <----------------------------------------------------------------------+   !
     !                                                                                   !
     do i = 1,ndir ! <---------------------------------------------------------------+   !
        do isurf = 1,2 ! <--------------------------------------------------------+  !   !
           kn(isurf) = ( &                                                        !  !   !
                LMN(1,isurf) * duv_ds(1,i,isurf)**2 + &                           !  !   !
                2._fp * LMN(2,isurf) * duv_ds(1,i,isurf) * duv_ds(2,i,isurf) + &  !  !   !
                LMN(3,isurf) * duv_ds(2,i,isurf)**2 ) / ( &                       !  !   !
                EFG(1,isurf) * duv_ds(1,i,isurf)**2 + &                           !  !   !
                2._fp * EFG(2,isurf) * duv_ds(1,i,isurf) * duv_ds(2,i,isurf) + &  !  !   !
                EFG(3,isurf) * duv_ds(2,i,isurf)**2 )                             !  !   !
        end do ! <----------------------------------------------------------------+  !   !
        denom = max(1._fp - cosn**2, EPSfp)                                          !   !
        curvature(i) = (kn(1)**2 - 2._fp*kn(1)*kn(2)*cosn + kn(2)**2) / denom        !   !
        curvature(i) = sqrt(max(0._fp, curvature(i)))                                !   !
     end do ! <----------------------------------------------------------------------+   !
  end if  ! <----------------------------------------------------------------------------+

end subroutine diffgeom_intersection_curve



subroutine simultaneous_point_inversions( &
     surf, &
     lowerb, &
     upperb, &
     stat, &
     uv, &
     xyz )
  use mod_math
  use mod_linalg
  use mod_diffgeom
  use mod_tolerances
  implicit none
  LOGICAL, PARAMETER :: DEBUG = .false. !( GLOBALDEBUG .AND. .true. )
  integer,           parameter     :: itmax = 2 + ceiling(-log10(EPSuv))
  type(ptr_surface), intent(in)    :: surf(2)
  real(kind=fp),     intent(in)    :: lowerb(4)
  real(kind=fp),     intent(in)    :: upperb(4)
  integer,           intent(out)   :: stat
  real(kind=fp),     intent(inout) :: uv(2,2)
  real(kind=fp),     intent(out)   :: xyz(3)
  real(kind=fp)                    :: xyzs(3,2), r(3)
  real(kind=fp)                    :: dxyz_duv(3,2), mat(2,2), rhs(2), duv(2)
  real(kind=fp), dimension(2)      :: cond, erruv
  integer                          :: it, isurf, ivar
  integer                          :: stat_refl
  
  !do isurf = 1,2
  !   uv(1:2,isurf) = max(lowerb(2*isurf-1:2*isurf), min(upperb(2*isurf-1:2*isurf), uv(1:2,isurf)))
  !end do
  
  IF ( DEBUG ) THEN
     PRINT *,''
     PRINT *,'SIMULTANEOUS_POINT_INVERSION'
     PRINT *,'UV =',uv
  END IF
  
  stat = 1
  erruv(1:2) = 0._fp
  cond(1:2) = 1._fp

  IF ( DEBUG ) PRINT *,'|X1 - X2| , |DUV1| , |DUV2| , eps*cond1 , eps*cond2'
  do it = 1,itmax ! <-------------------------------------------------------+
     !IF ( DEBUG ) PRINT *,'IT. #',IT
     !IF ( DEBUG ) THEN
     !   PRINT *,'UV='
     !   CALL PRINT_MAT(transpose(UV))
     !END IF
     do isurf = 1,2 ! <--------------+                                      !
        call eval( &                 !                                      !
             xyzs(:,isurf), &        !                                      !
             surf(isurf)%ptr, &      !                                      !
             uv(:,isurf) )           !                                      !
        !PRINT *,'xyz =',xyzs(:,isurf)
     end do ! <----------------------+                                      !
     r = 0.5_fp * (xyzs(:,1) - xyzs(:,2))                                   !
     !                                                                      !
     do isurf = 1,2 ! <--------------------------------------------------+  !
        do ivar = 1,2 ! <---------------+                                !  !
           call evald1( &               !                                !  !
                dxyz_duv(:,ivar), &     !                                !  !
                surf(isurf)%ptr, &      !                                !  !
                uv(:,isurf), &          !                                !  !
                ivar )                  !                                !  !
        end do ! <----------------------+                                !  !
        mat = matmul(transpose(dxyz_duv), dxyz_duv)                      !  !
        rhs = matmul(transpose(dxyz_duv), real((-1)**isurf, kind=fp)*r)  !  !
        call linsolve_svd( &                                             !  !
             duv, &                                                      !  !
             mat, &                                                      !  !
             rhs, &                                                      !  !
             2, &                                                        !  !
             2, &                                                        !  !
             1, &                                                        !  !
             cond(isurf), &                                              !  !
             tol=EPSmath )                                               !  !
        !                                                                !  !
        erruv(isurf) = sum(duv**2)                                       !  !
        !                                                                !  !
        ! modify Newton step => keep iterate inside feasible region      !  !
        call nd_box_reflexions( &                                        !  !
             uv(:,isurf), &                                              !  !
             lowerb(2*isurf-1:2*isurf), &                                !  !
             upperb(2*isurf-1:2*isurf), &                                !  !
             duv, &                                                      !  !
             2, &                                                        !  !
             stat_refl )                                                 !  !
        if ( stat_refl > 0 ) return                                      !  !
        !                                                                !  !
        ! update solution                                                !  !
        uv(:,isurf) = uv(:,isurf) + duv                                  !  !
     end do ! <----------------------------------------------------------+  !
     !                                                                      !
     IF ( DEBUG ) PRINT *,norm2(r), sqrt(erruv), EPSfp*cond
     ! termination criteria                                                 !
     if ( all(erruv < max(EPSuvsqr, EPSfpsqr*cond**2)) ) then ! <--------+  !
        if ( sum(r**2) < 0.25_fp*EPSxyzsqr ) then ! <----+               !  !
           stat = 0                                      !               !  !
           do isurf = 1,2 ! <--------------+             !               !  !
              call eval( &                 !             !               !  !
                   xyzs(:,isurf), &        !             !               !  !
                   surf(isurf)%ptr, &      !             !               !  !
                   uv(:,isurf) )           !             !               !  !
           end do ! <----------------------+             !               !  !
           xyz = 0.5_fp * sum(xyzs, 2)                   !               !  !
           IF ( DEBUG ) PRINT *,'CONVERGED, UV =',UV,', XYZ =',XYZ
        else ! ------------------------------------------+               !  !
           IF ( DEBUG ) PRINT *,'STAGNATION'
        end if ! <---------------------------------------+               !  !
        return                                                           !  !
     end if ! <----------------------------------------------------------+  !
  end do ! <----------------------------------------------------------------+

  if ( stat > 0 ) then
     if ( all(erruv < max(100._fp * EPSuvsqr, EPSfpsqr*cond**2)) ) then
        if ( sum(r**2) < 25_fp*EPSxyzsqr ) then
           stat = 0
           do isurf = 1,2 ! <--------------+
              call eval( &                 !
                   xyzs(:,isurf), &        !
                   surf(isurf)%ptr, &      !
                   uv(:,isurf) )           !
           end do ! <----------------------+
           xyz = 0.5_fp * sum(xyzs, 2)
           IF ( DEBUG ) PRINT *,'CONVERGED*, UV =',UV,', XYZ =',XYZ
        end if
     end if
  end if

end subroutine simultaneous_point_inversions



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



subroutine update_intersection_curves( &
     interdata )
  use mod_types_intersection
  use mod_diffgeom
  use mod_tolerances
  implicit none
  type(type_intersection_data), intent(inout), target :: interdata
  type(type_point_on_surface), pointer                :: pos => null()
  type(type_intersection_curve), pointer              :: curve => null()
  real(kind=fp)                                       :: uv(2,2), xyz(3,2), err
  real(kind=fp), dimension(4)                         :: lowerb, upperb
  integer                                             :: stat
  integer                                             :: ipoint, ipos, icurv, i, j
 
  upperb(:) = 1._fp + EPSuv
  lowerb = -upperb

  do ipoint = 1,interdata%np
     interdata%points(ipoint)%xyz(1:3) = 0._fp
     pos => interdata%points(ipoint)%pos
     do ipos = 1,interdata%points(ipoint)%npos
        call eval( &
             xyz(1:3,1), &
             pos%surf, &
             pos%uv )
        interdata%points(ipoint)%xyz = interdata%points(ipoint)%xyz + xyz(1:3,1)
        pos => pos%next
     end do
     interdata%points(ipoint)%xyz = interdata%points(ipoint)%xyz / real(interdata%points(ipoint)%npos, kind=fp)
  end do
  
  do icurv = 1,interdata%nc ! <---------------------------------+
     !PRINT *,'ICURV =',ICURV
     curve => interdata%curves(icurv)                           !
     if ( .not.curve%smooth ) cycle                             !
     !                                                          !
     do i = 1,curve%polyline%np ! <--------------------------+  !
        uv = curve%polyline%uv(1:2,1:2,i)                    !  !
        do j = 1,2 ! <----------------+                      !  !
           call eval( &               !                      !  !
                xyz(:,j), &           !                      !  !
                curve%surf(j)%ptr, &  !                      !  !
                uv(:,j) )             !                      !  !
        end do ! <--------------------+                      !  !
        curve%polyline%xyz(1:3,i) = 0.5_fp * sum(xyz, 2)     !  !
        !                                                    !  !
        err = sum((xyz(:,1) - xyz(:,2))**2)                  !  !
        if ( err > EPSxyzsqr ) then ! <-----------------+    !  !
           PRINT *,'CURVE #',ICURV
           PRINT *,I,sqrt(err)
           call simultaneous_point_inversions( &        !    !  !
                curve%surf, &                           !    !  !
                lowerb, &                               !    !  !
                upperb, &                               !    !  !
                stat, &                                 !    !  !
                uv, &                                   !    !  !
                xyz(:,1) )                              !    !  !
           if ( stat > 0 ) then ! <------------------+  !    !  !
              ! keep current coordinates             !  !    !  !
              PRINT *,'update_intersection_curves: *** /!\ FAILED TO RECOVER POINT ON TANGENT INTERSECTION, ERR=',SQRT(ERR)
           else ! -----------------------------------+  !    !  !
              curve%polyline%uv(1:2,1:2,i) = uv      !  !    !  !
              PRINT *,I,NORM2(XYZ(:,1) - CURVE%POLYLINE%XYZ(:,I))
              PRINT *,I,uv
              curve%polyline%xyz(1:3,i) = xyz(:,1)   !  !    !  !
           end if ! <--------------------------------+  !    !  !
        end if ! <--------------------------------------+    !  !
        !                                                    !  !
     end do ! <----------------------------------------------+  !
  end do ! <----------------------------------------------------+
  
end subroutine update_intersection_curves



subroutine diffgeom_intersection( &
     surf, &
     uv, &
     duv_ds, &
     dxyz_ds, &
     stat, &
     curvature )
  ! Computes tangent direction(s) and curvature of the intersection curve between
  ! two surfaces at a given intersection point.
  ! Returns: stat =-1 if one surface is singular at that point (normal = zero vector)
  !               = 0 if the point is on a single transveral intersection curve
  !               = 1 if the point is on a single tangential intersection curve
  !               = 2 if the point is a branch point (junction of two intersection curves)
  !               = 3 if the point is an isolated tangential contact point (degenerate curve)
  !               = 4 if the point is an high-order tangential contact point (possible degeneracy)
  use mod_math
  use mod_diffgeom
  ! ( see "Shape interrogation for computer aided design and manufacturing", &
  ! Patrikalakis et al. (2009), pp.166-175, 
  ! and "Tracing surface intersections with validated ODE system solver", Mukundan et al (2004)
  implicit none
  real(kind=fp), parameter             :: EPS = 1.d-9
  type(ptr_surface),       intent(in)  :: surf(2)
  real(kind=fp),           intent(in)  :: uv(2,2)
  real(kind=fp),           intent(out) :: duv_ds(2,2,2) ! u/v, #branch, #surf
  real(kind=fp),           intent(out) :: dxyz_ds(3,2)
  integer,                 intent(out) :: stat
  real(kind=fp), optional, intent(out) :: curvature(2)
  real(kind=fp)                        :: dxyz_duv(3,2,2), n(3,2), d2xyz_duv2(3)
  real(kind=fp)                        :: dotn, sgndotn, aux(3), w
  real(kind=fp), dimension(3,2)        :: EFG, LMN
  real(kind=fp), dimension(2)          :: detEFG, invsqrtdetEFG
  real(kind=fp)                        :: A(2,2), B(3), discr, invdenom, kn(2)
  integer                              :: isurf, jsurf, ivar, i, j

  do isurf = 1,2 ! <------------------------------------------------------------------+
     ! tangent vectors                                                                !
     do ivar = 1,2 ! <-----------------------------+                                  !
        call evald1( dxyz_duv(:,ivar,isurf), &     !                                  !
             surf(isurf)%ptr, &                    !                                  !
             uv(:,isurf), &                        !                                  !
             ivar )                                !                                  !
     end do ! <------------------------------------+                                  !
     ! First Fundamental Form coefficients                                            !
     EFG(1,isurf) = dot_product(dxyz_duv(:,1,isurf), dxyz_duv(:,1,isurf))             !
     EFG(2,isurf) = dot_product(dxyz_duv(:,1,isurf), dxyz_duv(:,2,isurf))             !
     EFG(3,isurf) = dot_product(dxyz_duv(:,2,isurf), dxyz_duv(:,2,isurf))             !
     detEFG(isurf) = EFG(1,isurf)*EFG(3,isurf) - EFG(2,isurf)**2                      !
     invsqrtdetEFG(isurf) = sqrt(1._fp / detEFG(isurf))                               !
     !                                                                                !
     ! unit normal vector                                                             !
     n(:,isurf) = cross(dxyz_duv(:,1,isurf), dxyz_duv(:,2,isurf))                     !
     !PRINT *,'N =',N(:,ISURF)
     !n(:,isurf) = n(:,isurf) * invsqrtdetEFG(isurf)
     !n(:,isurf) = n(:,isurf) / norm2(n(:,isurf))
  end do ! <--------------------------------------------------------------------------+
  !PRINT *,'DXYZ_DUV(1) ='
  !CALL PRINT_MAT(DXYZ_DUV(:,:,1))
  !PRINT *,'EFG ='
  !CALL PRINT_MAT(EFG)
  !PRINT *,'DETEFG =',DETEFG
  !PRINT *,'INVSQRTDETEFG =',INVSQRTDETEFG
  !PRINT *,'N ='
  !CALL PRINT_MAT(N)
  dotn = dot_product(n(:,1), n(:,2))
  !PRINT *,'DOTN =',DOTN

  if ( abs(dotn) > (1._fp - EPS)*product(sqrt(detEFG)) ) then ! <---------------------+
     !! tangential intersection                                                       !
     sgndotn = sign(1._fp, dotn)                                                      !
     ! Second Fundamental Form coefficients                                           !
     do isurf = 1,2 ! <-----------------------------------------------+               !
        do ivar = 1,3 ! <-----------------------------------------+   !               !
           call evald2( d2xyz_duv2, &                             !   !               !
                surf(isurf)%ptr, &                                !   !               !
                uv(:,isurf), &                                    !   !               !
                ivar )                                            !   !               !
           LMN(ivar,isurf) = dot_product(n(:,isurf), d2xyz_duv2)  !   !               !
        end do ! <------------------------------------------------+   !               !
     end do ! <-------------------------------------------------------+               !
     where( abs(LMN) < EPSmath ) LMN = 0._fp
     do isurf = 1,2
        LMN(:,isurf) = LMN(:,isurf) * invsqrtdetEFG(isurf)
     end do
     !PRINT *,'DETLMN =',LMN(1,:)*LMN(3,:) - LMN(2,:)**2
     !PRINT *,'LMN ='
     !CALL PRINT_MAT(LMN)
     jsurf = minloc(detEFG, 1)
     isurf = 1+mod(jsurf,2)
     ! linear realtion duvj_ds = A * duvi_ds                                          !
     do i = 1,2 ! <---------------------------------------------------+               !
        aux = real((-1)**i, kind=fp) * &                              !               !
             cross(n(:,jsurf), dxyz_duv(:,1+mod(i,2),isurf))          !               !
        !PRINT *,'AUX =',AUX
        do j = 1,2 ! <-------------------------------------+          !               !
           A(i,j) = dot_product(aux, dxyz_duv(:,j,jsurf))  !          !               !
        end do ! <-----------------------------------------+          !               !
     end do ! <-------------------------------------------------------+               !
     !A = sign(invsqrtdetEFG(jsurf), dotn)*A
     where( abs(A) < EPSmath ) A = 0._fp
     A = A / dotn
     !PRINT *,'A ='
     !CALL PRINT_MAT(A)
     ! quadratic equation for duv1_ds                                                 !
     B(1) = LMN(1,isurf)*A(1,1)**2 + 2._fp*LMN(2,isurf)*A(1,1)*A(2,1) + LMN(3,isurf)*A(2,1)**2
     B(2) = LMN(1,isurf)*A(1,1)*A(1,2) + LMN(2,isurf)*(A(1,2)*A(2,1) + A(1,1)*A(2,2)) + LMN(3,isurf)*A(2,1)*A(2,2)
     B(3) = LMN(1,isurf)*A(1,2)**2 + 2._fp*LMN(2,isurf)*A(1,2)*A(2,2) + LMN(3,isurf)*A(2,2)**2
     !B = B / detEFG(jsurf)
     B = LMN(:,jsurf) - sign(1._fp, dotn)*B
     where( abs(B) < EPSmath ) B = 0._fp
     !PRINT *,'B =',B
     discr = B(2)**2 - B(1)*B(3)                                                      !
     !A = sign(invsqrtdetEFG(jsurf), dotn)*A
     !PRINT *,'DISCR =',DISCR
     !                                                                                !
     if ( all(abs(B) < EPS) ) then ! <-------------------------------------+          !
        ! high-order contact point (dxyz_ds cannot be determined)          !          !
        stat = 4                                                           !          !
        return                                                             !          !
        !                                                                  !          !
     elseif ( discr < -EPS ) then ! ---------------------------------------+          !
        ! isolated tangential contact point (dxyz_ds undefined)            !          !
        stat = 3                                                           !          !
        return                                                             !          !
        !                                                                  !          !
     else ! ---------------------------------------------------------------+          !
        if ( discr < EPS ) then ! <------------------------------+         !          !
           ! tangential intersection curve (one single dxyz_ds)  !         !          !
           stat = 1                                              !         !          !
           discr = 0._fp                                         !         !          !
        else ! --------------------------------------------------+         !          !
           ! branch point (two disctinct dxyz_ds)                !         !          !
           stat = 2                                              !         !          !
           discr = sqrt(discr)                                   !         !          !
        end if ! <-----------------------------------------------+         !          !
        !                                                                  !          !
        ! tangent direction(s) to the intersection branch(es)              !          !
        ! (the number of branches is equal to stat)                        !          !
        if ( abs(B(1)) > abs(B(3)) ) then ! <---------------------------+  !          !
           do i = 1,stat ! <-----------------------------------------+  !  !          !
              w = (discr*real((-1)**i, kind=fp) - B(2)) / B(1)       !  !  !          !
              duv_ds(:,i,jsurf) = [w, 1._fp] / &                     !  !  !          !
                   sqrt(EFG(1,jsurf)*w**2 + &                        !  !  !          !
                   2._fp*EFG(2,jsurf)*w + &                          !  !  !          !
                   EFG(3,jsurf))                                     !  !  !          !
           end do ! <------------------------------------------------+  !  !          !
        else ! ---------------------------------------------------------+  !          !
           do i = 1,stat ! <-----------------------------------------+  !  !          !
              w = (discr*real((-1)**i, kind=fp) - B(2)) / B(3)       !  !  !          !
              duv_ds(:,i,jsurf) = [1._fp, w] / &                     !  !  !          !
                   sqrt(EFG(1,jsurf) + &                             !  !  !          !
                   2._fp*EFG(2,jsurf)*w + &                          !  !  !          !
                   EFG(3,jsurf)*w**2)                                !  !  !          !
           end do ! <------------------------------------------------+  !  !          !
        end if ! <------------------------------------------------------+  !          !
        do i = 1,stat ! <--------------------------------------------+     !          !
           duv_ds(:,i,isurf) = matmul(A, duv_ds(:,i,jsurf))          !     !          !
           dxyz_ds(:,i) = matmul(dxyz_duv(:,:,jsurf), &              !     !          !
                duv_ds(:,i,jsurf))                                   !     !          !
        end do ! <---------------------------------------------------+     !          !
        !                                                                  !          !
     end if ! <------------------------------------------------------------+          !
     !                                                                                !
     !                                                                                !
  else ! -----------------------------------------------------------------------------+
     !! transversal intersection                                                      !
     stat = 0                                                                         !
     invdenom = 1._fp * sqrt(1._fp - dotn**2)                                         !
     dxyz_ds(:,1) = cross(n(:,1), n(:,2))                                             !
     dxyz_ds(:,1) = dxyz_ds(:,1) / norm2(dxyz_ds(:,1))                                !
     !                                                                                !
     do isurf = 1,2 ! <-----------------------------------------------+               !
        aux(1) = dot_product(dxyz_ds(:,1), dxyz_duv(:,1,isurf))       !               !
        aux(2) = dot_product(dxyz_ds(:,1), dxyz_duv(:,2,isurf))       !               !
        duv_ds(1,1,isurf) = aux(1)*EFG(3,isurf) - aux(2)*EFG(2,isurf) !               !
        duv_ds(2,1,isurf) = aux(2)*EFG(1,isurf) - aux(1)*EFG(2,isurf) !               !
        duv_ds(:,1,isurf) = duv_ds(:,1,isurf) / detEFG(isurf)         !               !
     end do ! <-------------------------------------------------------+               !
     !                                                                                !
     if ( present(curvature) ) then ! <--------------------------------------------+  !
        do isurf = 1,2 ! <--------------------------------------------------+      !  !
           ! Second Fundamental Form coefficients                           !      !  !
           do ivar = 1,3 ! <-----------------------------------------+      !      !  !
              call evald2( d2xyz_duv2, &                             !      !      !  !
                   surf(isurf)%ptr, &                                !      !      !  !
                   uv(:,isurf), &                                    !      !      !  !
                   ivar )                                            !      !      !  !
              LMN(ivar,isurf) = dot_product(n(:,isurf), d2xyz_duv2)  !      !      !  !
           end do ! <------------------------------------------------+      !      !  !
           LMN(:,isurf) = LMN(:,isurf) * invsqrtdetEFG(isurf)               !      !  !
           !                                                                !      !  !
           ! normal curvature                                               !      !  !
           kn(isurf) =  &                                                   !      !  !
                LMN(1,isurf)*duv_ds(1,1,isurf)**2 + &                       !      !  !
                2._fp*LMN(2,isurf)*duv_ds(1,1,isurf)*duv_ds(2,1,isurf) + &  !      !  !
                LMN(3,isurf)*duv_ds(2,1,isurf)**2                           !      !  !
        end do ! <----------------------------------------------------------+      !  !
        !                                                                          !  !
        curvature(1) = sqrt(kn(1)**2 + kn(2)**2 - 2._fp*kn(1)*kn(2)*dotn)*invdenom !  !
     end if ! <--------------------------------------------------------------------+  !
  end if ! <--------------------------------------------------------------------------+
  
end subroutine diffgeom_intersection



subroutine check_unicity( &
     vec, &
     dim, &
     array, &
     n, &
     tol, &
     id )
  use mod_math
  implicit none
  LOGICAL, PARAMETER :: DEBUG = ( GLOBALDEBUG .AND. .false. )
  integer,       intent(in)  :: dim
  real(kind=fp), intent(in)  :: vec(dim)
  integer,       intent(in)  :: n
  real(kind=fp), intent(in)  :: array(dim,n)
  real(kind=fp), intent(in)  :: tol
  integer,       intent(out) :: id
  real(kind=fp)              :: tolsqr

  IF ( DEBUG ) PRINT *,'CHECK UNICITY OF',VEC(1:DIM)
  
  if ( n < 1 ) then
     id = 1
  else
     tolsqr = tol**2
     do id = 1,n
        IF ( DEBUG ) PRINT *,array(1:dim,id), NORM2(vec(1:dim) - array(1:dim,id))
        if ( sum((vec(1:dim) - array(1:dim,id))**2) < tolsqr ) then
           IF ( DEBUG ) PRINT '(A17,1x,I0,1x,A1)','DUPLICATE POINT (',ID,')'
           exit
        end if
     end do
  end if
  IF ( DEBUG ) PRINT '(A4,1X,I0)','ID =',ID
  
end subroutine check_unicity



subroutine intersect_border_surface( &
     root_s, &
     root_c, &
     region, &
     icurv, &
     ivar, &
     ival, &
     ipts_ss, &
     npts_ss, &
     uvxyz, &
     nuvxyz, &
     ipts_bs, &
     npts_bs, &
     stat_degeneracy )
  use mod_util
  use mod_math
  use mod_polynomial
  use mod_diffgeom
  use mod_regiontree
  use mod_types_intersection
  use mod_tolerances
  implicit none
  LOGICAL, PARAMETER :: DEBUG = ( GLOBALDEBUG .AND. .true. )
  type(ptr_surface),          intent(in)    :: root_s(2)
  type(type_curve),           intent(in)    :: root_c
  type(ptr_region),           intent(inout) :: region(2)
  integer,                    intent(in)    :: icurv, ivar, ival
  integer, allocatable,       intent(in)    :: ipts_ss(:)
  integer,                    intent(in)    :: npts_ss
  real(kind=fp), allocatable, intent(inout) :: uvxyz(:,:)
  integer,                    intent(inout) :: nuvxyz
  integer, allocatable,       intent(inout) :: ipts_bs(:)
  integer,                    intent(inout) :: npts_bs
  integer,                    intent(inout) :: stat_degeneracy
  integer                                   :: isurf, jvar
  type(type_region)                         :: region_c
  type(type_region)                         :: region_s
  real(kind=fp), allocatable                :: tuvxyz(:,:)
  integer                                   :: ntuvxyz, ntuvxyz_tmp
  real(kind=fp)                             :: tmp(8), uv(2,2)
  integer                                   :: ipt, jpt

  IF ( DEBUG ) THEN
     PRINT *,''; PRINT *,'';
     PRINT *,'INTERSECT_BORDER_SURFACE'
     PRINT *,'ICURV, IVAR, IVAL =',ICURV,IVAR,IVAL
     PRINT *,'UVBOXES ='
     DO ISURF = 1,2 ! <-----------------+
        PRINT *,REGION(ISURF)%PTR%UVBOX !
     END DO ! <-------------------------+
  END IF

  isurf = 1 + mod(icurv,2)
  jvar  = 1 + mod(ivar ,2)

  ! initialize a temporary curve region tree
  call init_region( &
       region_c, &
       1, &
       [ -1._fp, 1._fp ] ) 

  ! initialize a temporary surface region tree
  call init_region( &
       region_s, &
       2, &
       region(isurf)%ptr%uvbox ) 

  ! compute the curve's Bezier control points
  allocate(region_c%poly(1))
  allocate(region_c%poly(1)%ptr)
  call cheb2bern( &
       root_c%x, &
       region_c%poly(1)%ptr )

  ! copy the surface's Bezier control points
  allocate(region_s%poly(1))
  region_s%poly(1)%ptr => region(isurf)%ptr%poly(1)%ptr
  
  ntuvxyz = 0
  allocate(tuvxyz(7,10))
  IF ( DEBUG ) THEN
     PRINT *,'NPTS_SS =',npts_ss
     IF ( NPTS_SS > 0 ) THEN
        PRINT *,ipts_ss(1:npts_ss)
        CALL PRINT_MAT( TRANSPOSE(UVXYZ(1:4,ipts_ss(1:npts_ss))) )
     END IF
  END IF

  ! check if some of these points are located on the current border
  do jpt = 1,npts_ss ! <--------------------------------------------------------------------+
     ipt = ipts_ss(jpt)                                                                     !
     tmp(1) = uvxyz(2*(icurv-1)+ivar,ipt)                                                   !
     !                                                                                      !
     if ( abs( tmp(1) - region(icurv)%ptr%uvbox(2*(ivar-1)+ival) ) < EPSuv ) then ! <---+   !
        tmp(1) = ab2n1p1( &                                                             !   !
             uvxyz(2*(icurv-1)+jvar,ipt), &                                             !   !
             region(icurv)%ptr%uvbox(2*jvar-1), &                                       !   !
             region(icurv)%ptr%uvbox(2*jvar  ) )                                        !   !
        tmp(2:3) = uvxyz(2*isurf-1:2*isurf,ipt)                                         !   !
        tmp(4:6) = uvxyz(5:7,ipt)                                                       !   !
        tmp(7)   = uvxyz(8,ipt)                                                         !   !
        IF ( DEBUG ) THEN
           PRINT *,'+ 1PT ON BORDER :'
           PRINT *,'    UV =',UVXYZ(1:4,IPT)
           PRINT *,'   XYZ =',UVXYZ(5:7,IPT)
           PRINT *,' TOLUV =',UVXYZ(8,IPT)
           PRINT *,'TUVXYZ =',TMP(1:6)
        END IF
        !                                                                               !   !
        call append_vector( &                                                           !   !
             tmp(1:7), &                                                                !   !
             7, &                                                                       !   !
             tuvxyz, &                                                                  !   !
             ntuvxyz )                                                                  !   !
        !                                                                               !   !
        call append_n( &                                                                !   !
             ipts_bs, &                                                                 !   !
             npts_bs, &                                                                 !   !
             [ipt], &                                                                   !   !
             1, &                                                                       !   !
             unique=.true. )                                                            !   !
        !                                                                               !   !
        call add_point_bottom_up( region_c, ntuvxyz )
        call add_point_bottom_up( region_s, ntuvxyz )
     end if ! <-------------------------------------------------------------------------+   !
  end do ! <--------------------------------------------------------------------------------+
  ntuvxyz_tmp = ntuvxyz

  IF ( DEBUG ) THEN
     PRINT *,'NTUVXYZ_TMP =',NTUVXYZ!_TMP
     PRINT *,'TUVXYZ ='
     CALL PRINT_MAT( TRANSPOSE(TUVXYZ(1:6,1:NTUVXYZ)) )
     PRINT *,'TOLTUV =',TUVXYZ(7,1:NTUVXYZ)
  END IF


  ! compute the curve-surface intersection
  stat_degeneracy = 0
  !ntuvxyz = 0
  call intersect_curve_surface( &
       root_c, &
       root_s(isurf)%ptr, &
       region_c, &
       region_s, &
       tuvxyz, &
       ntuvxyz, &
       stat_degeneracy )
  if ( stat_degeneracy == 50 ) stat_degeneracy = 0
  IF ( DEBUG ) THEN
     PRINT *,'NTUVXYZ =',NTUVXYZ
     IF ( NTUVXYZ > 0 ) THEN
        CALL PRINT_MAT( TRANSPOSE(TUVXYZ(:,1:NTUVXYZ)) )
        PRINT *,'TOLTUV =',TUVXYZ(7,1:NTUVXYZ)
     END IF
     !CALL EXPORT_REGION_TREE( REGION_C, 'dev_intersection/treebsi_c.dat' )
     !CALL EXPORT_REGION_TREE( REGION_S, 'dev_intersection/treebsi_s.dat' )
  END IF

  IF ( stat_degeneracy > 0 ) THEN
     CALL WRITE_POLYNOMIAL( ROOT_C%X, 'dev_intersection/debugbsi_c.cheb' )
     CALL WRITE_POLYNOMIAL( ROOT_S(icurv)%ptr%X, 'dev_intersection/debugbsi_scurv.cheb' )
     CALL WRITE_POLYNOMIAL( ROOT_S(isurf)%ptr%X, 'dev_intersection/debugbsi_s.cheb' )
     CALL WRITE_POLYNOMIAL( REGION_C%POLY(1)%PTR, 'dev_intersection/debugbsi_c.bern' )
     CALL WRITE_POLYNOMIAL( REGION_S%POLY(1)%PTR, 'dev_intersection/debugbsi_s.bern' )
     CALL EXPORT_REGION_TREE( REGION_C, 'dev_intersection/treebsi_c.dat' )
     CALL EXPORT_REGION_TREE( REGION_S, 'dev_intersection/treebsi_s.dat' )
  END IF

  ! free the curve region tree
  call free_polynomial(region_c%poly(1)%ptr)
  deallocate(region_c%poly(1)%ptr)
  deallocate(region_c%poly)
  call free_region_tree(region_c)

  ! free the (temporary) surface region tree
  nullify(region_s%poly(1)%ptr)
  call free_region_tree(region_s)


  ! manage the newly discovered intersection points
  do ipt = ntuvxyz_tmp+1,ntuvxyz ! <--------------------------------+
     ! convert from (t,u,v) to (u1,v1,u2,v2)                        !
     uv(ivar,icurv) = region(icurv)%ptr%uvbox(2*(ivar-1)+ival)      !
     uv(jvar,icurv) = n1p12ab( &                                    !
          tuvxyz(1,ipt), &                                          !
          region(icurv)%ptr%uvbox(2*jvar-1), &                      !
          region(icurv)%ptr%uvbox(2*jvar  ) )                       !
     uv(:,isurf) = tuvxyz(2:3,ipt)                                  !
     !                                                              !
     tmp(1:2) = uv(:,1)                                             !
     tmp(3:4) = uv(:,2)                                             !
     tmp(5:7) = tuvxyz(4:6,ipt)                                     !
     tmp(8)   = tuvxyz(7,ipt)                                       !
     !                                                              !
     ! check unicity                                                !
     if ( nuvxyz > 0 ) then ! <---------------+                     !
        call check_unicity_with_tolerance( &  !                     !
             tmp(1:4), &                      !                     !
             tmp(5:7), &                      !                     !
             tmp(8), &                        !                     !
             4, &                             !                     !
             uvxyz(1:4,1:nuvxyz), &           !                     !
             uvxyz(5:7,1:nuvxyz), &           !                     !
             uvxyz(8,1:nuvxyz), &             !                     !
             nuvxyz, &                        !                     !
             jpt )                            !                     !
     else ! ----------------------------------+                     !
        jpt = 1                               !                     !
     end if ! <-------------------------------+                     !
     !                                                              !
     if ( jpt > nuvxyz ) then ! <-----+                             !
        call append_vector( &         !                             !
             tmp(1:8), &              !                             !
             8, &                     !                             !
             uvxyz, &                 !                             !
             nuvxyz )                 !                             !
        jpt = nuvxyz                  !                             !
     end if ! <-----------------------+                             !
     !                                                              !
     call append_n( &                                               !
          ipts_bs, &                                                !
          npts_bs, &                                                !
          [jpt], &                                                  !
          1, &                                                      !
          unique=.true. )                                           !
     !                                                              !
     IF ( DEBUG ) THEN
        PRINT *,'+1 UVXYZ :',TMP
        PRINT *,'IPTS_BS <---',NUVXYZ
     END IF
  end do ! <--------------------------------------------------------+
  
  if ( allocated(tuvxyz) ) deallocate(tuvxyz)


end subroutine intersect_border_surface



subroutine add_intersection_curve( &
     interdata, &
     param_vector, &
     iendpoints, &
     uvbox )
  use mod_math
  use mod_types_intersection
  implicit none
  LOGICAL, PARAMETER :: DEBUG = ( GLOBALDEBUG .AND. .false. )
  integer, parameter                          :: PARAM_xtra_nc = 10
  type(type_intersection_data), intent(inout) :: interdata
  real(kind=fp),                intent(in)    :: param_vector(3)
  integer,                      intent(in)    :: iendpoints(2)
  real(kind=fp),                intent(in)    :: uvbox(2,2,2)
  type(type_intersection_data)                :: tmp
  integer                                     :: ic

  IF ( DEBUG ) THEN
     PRINT *,''; PRINT *,'';
     PRINT *,'ADD_INTERSECTION_CURVE'
     PRINT *,'INTERDATA%NC =',interdata%nc
     PRINT *,'ALLOCATED(INTERDATA%CURVES)?',allocated(interdata%curves)
     IF ( allocated(interdata%curves) ) PRINT *,'SIZE =',size(interdata%curves)
  END IF

  if ( .not.allocated(interdata%curves) ) allocate(interdata%curves(PARAM_xtra_nc) )

  ! reallocate if necessary
  if ( interdata%nc + 1 > size(interdata%curves) ) then
     allocate(tmp%curves(interdata%nc))
     IF ( DEBUG ) PRINT *,'INTERDATA%CURVES --> TMP'
     call transfer_intersection_curves( &
          from=interdata, &
          to=tmp )

     IF ( DEBUG ) PRINT *,'INTERDATA%CURVES <-- TMP'
     allocate(interdata%curves(tmp%nc + PARAM_xtra_nc))
     call transfer_intersection_curves( &
          from=tmp, &
          to=interdata )
  end if

  !! insert new curve
  interdata%nc = interdata%nc + 1
  ic = interdata%nc
  interdata%curves(ic)%uvbox          = uvbox
  interdata%curves(ic)%param_vector   = param_vector
  ! add endpoints as split points
  interdata%curves(ic)%nsplit = 2
  allocate(interdata%curves(ic)%isplit(2,2))
  interdata%curves(ic)%isplit(1,1:2) = iendpoints
  interdata%curves(ic)%isplit(2,1:2) = -1 ! do not forget to fill in when tracing the polyline!
  IF ( DEBUG ) PRINT *,'ALL IS OK, NC =',interdata%nc

end subroutine add_intersection_curve



subroutine append_vector( &
     vec, &
     dim, &
     array, &
     n )
  use mod_math
  implicit none
  integer,                    intent(in)    :: dim
  real(kind=fp),              intent(in)    :: vec(dim)
  real(kind=fp), allocatable, intent(inout) :: array(:,:)
  integer,                    intent(inout) :: n
  real(kind=fp), allocatable                :: tmp(:,:)

  if ( .not.allocated(array) ) allocate( array(dim,1) )
  if ( dim > size(array,1) ) STOP 'append_vector : dim > size(array,1)'

  if ( n + 1 > size(array,2) ) then
     call move_alloc( from=array, to=tmp )
     allocate( array(dim,n+1) )
     array(:,1:n) = tmp(:,1:n)
     deallocate(tmp)
  end if
  n = n + 1
  array(1:dim,n) = vec(1:dim)

end subroutine append_vector



subroutine intersect_all_surfaces( &
     surf, &
     nsurf, &
     interdata_global, &
     mask, &
     tolchord, &
     hmin, &
     hmax )
  USE MOD_UTIL
  use mod_math
  use mod_polynomial
  use mod_diffgeom
  use mod_regiontree
  use mod_types_intersection
  implicit none
  LOGICAL, PARAMETER :: DEBUG = ( GLOBALDEBUG .AND. .TRUE. )
  CHARACTER(2) :: STRNUM
  integer,                      intent(in)    :: nsurf            ! number of surfaces
  type(type_surface), target,   intent(in)    :: surf(nsurf)      ! surfaces
  type(type_intersection_data), intent(inout) :: interdata_global ! global intersection data collection
  logical, optional,            intent(in)    :: mask(nsurf)      ! used to skip some surfaces when computing intersections
  real(kind=fp),                intent(in)    :: tolchord, hmin, hmax
  type(ptr_surface)                           :: surfpair(2)   
  type(type_region), target                   :: root(nsurf)
  type(ptr_region)                            :: region(2)
  type(type_intersection_data)                :: interdata_local
  real(kind=fp), allocatable                  :: uvxyz(:,:)
  integer                                     :: nuvxyz
  integer                                     :: stat_degeneracy
  integer                                     :: isurf, jsurf, ic, i

  ! initialize all region trees
  do isurf = 1,nsurf
     call init_region( &
       root(isurf), &
       2, &
       [( [-1._fp, 1._fp], i=1,2 )] ) 

     ! compute Bezier control points for position and pseudonormal vectors
     allocate(root(isurf)%poly(2))
     allocate(root(isurf)%poly(1)%ptr, root(isurf)%poly(2)%ptr)
     call cheb2bern( &
          surf(isurf)%x, &
          root(isurf)%poly(1)%ptr )
     call cheb2bern( &
          surf(isurf)%pn, &
          root(isurf)%poly(2)%ptr )

     !IF ( DEBUG ) THEN
     !   WRITE (STRNUM,'(I1)') ISURF
     !   CALL WRITE_POLYNOMIAL( root(isurf)%poly(1)%ptr, 'dev_intersection/surfroot' // strnum // '_x.bern'  )
     !   CALL WRITE_POLYNOMIAL( root(isurf)%poly(2)%ptr, 'dev_intersection/surfroot' // strnum // '_pn.bern' )
     !END IF
  end do


  nuvxyz = 0
  allocate(uvxyz(8,10))

  ! loop over all pairs of DISTINCT surfaces and compute their intersection
  outer : do isurf = 1,nsurf-1
     if ( present(mask) ) then
        if ( .not.mask(isurf) ) cycle outer
     end if

     inner : do jsurf = isurf+1,nsurf
        if ( present(mask) ) then
           if ( .not.mask(jsurf) ) cycle inner
        end if

        ! check if there is a previously discovered tangential intersection curve between 
        ! the two intersected surfaces
        do ic = 1,interdata_global%nc
           if ( .not.interdata_global%curves(ic)%smooth ) cycle
           if ( ( &
                associated( interdata_global%curves(ic)%surf(1)%ptr, surf(isurf) ) .and. &
                associated( interdata_global%curves(ic)%surf(2)%ptr, surf(jsurf) ) ) .or. &
                ( &
                associated( interdata_global%curves(ic)%surf(1)%ptr, surf(jsurf) ) .and. &
                associated( interdata_global%curves(ic)%surf(2)%ptr, surf(isurf) ) ) &
                ) cycle inner ! we skip this pair of surfaces
        end do

        IF ( isurf == 3 .or. jsurf == 3 ) then
           PRINT *,''; PRINT *,'PAIR :',ISURF,JSURF
        end IF
        ! initialize pointers to surfaces and region trees
        region(1)%ptr   => root(isurf)
        region(2)%ptr   => root(jsurf)
        surfpair(1)%ptr => surf(isurf)
        surfpair(2)%ptr => surf(jsurf)
        stat_degeneracy =  0

        ! if some isolated or branch singular intersection points between 
        ! the two intersected surfaces have already been discovered, copy them into
        ! the local intersection data collection
        ! (...)        

        allocate(interdata_local%points(10), interdata_local%curves(10))
        interdata_local%np = 0
        interdata_local%nc = 0
        
        ! compute topology of the intersection between the current pair of surfaces
        nuvxyz = 0
        call intersect_surface_pair( &
             surfpair, &
             region, &
             interdata_local, &
             uvxyz, &
             nuvxyz, &
             stat_degeneracy ) 
        IF ( isurf == 3 .or. jsurf == 3 ) then!.true. ) THEN
           PRINT *,'STAT_DEGENERACY =',stat_degeneracy
           !IF ( NUVXYZ > 0 ) THEN
           !   CALL WRITE_MATRIX( TRANSPOSE(UVXYZ(1:7,1:NUVXYZ)), NUVXYZ, 7, &
           !        'dev_intersection/uvxyz.dat' )
           !END IF
           !CALL WRITE_INTERSECTION_DATA( INTERDATA_LOCAL, &
           !     'dev_intersection/interdata_points.dat', &
           !     'dev_intersection/interdata_curves.dat' )

           PRINT *,NUVXYZ,            ' INTERSECTION POINT(S)'
           PRINT *,INTERDATA_LOCAL%NC,' INTERSECTION CURVE(S)'
           IF ( stat_degeneracy /= 0 ) PAUSE

           !DO i = 1,2
           !   WRITE (STRNUM,'(I1)') i
           !   CALL EXPORT_REGION_TREE( REGION(i)%PTR, 'dev_intersection/treessi_' // strnum // '.dat' )
           !END DO
        END IF

        ! trace intersection curves and append them to the global intersection data collection
        call merge_intersection_data( &
             surfpair, &
             uvxyz, &
             nuvxyz, &
             interdata_local, &
             interdata_global, &
             tolchord, &
             hmin, &
             hmax )

        ! reset local intersection data collection
        call free_intersection_data(interdata_local)
        call free_ipts(region(1)%ptr)
        call free_ipts(region(2)%ptr)

        ! if a degeneracy has been encountered, report it
        if ( stat_degeneracy > 0 ) exit outer

     end do inner
  end do outer

  if ( allocated(uvxyz) ) deallocate(uvxyz)


  ! free all region trees
  do isurf = 1,nsurf
     IF ( DEBUG ) THEN
        WRITE (STRNUM,'(I2.2)') isurf
        CALL EXPORT_REGION_TREE( root(isurf), 'dev_intersection/treessi_' // strnum // '.dat' )
     END IF

     call free_polynomial(root(isurf)%poly(1)%ptr)
     call free_polynomial(root(isurf)%poly(2)%ptr)
     deallocate(root(isurf)%poly(1)%ptr, root(isurf)%poly(2)%ptr)
     deallocate(root(isurf)%poly)
     call free_region_tree(root(isurf)) 
  end do

  nullify(region(1)%ptr, region(2)%ptr)

  do ic = 1,interdata_global%nc
     if ( .not.allocated(interdata_global%curves(ic)%iedge) ) then
        allocate(interdata_global%curves(ic)%iedge(interdata_global%curves(ic)%nsplit-1))
        interdata_global%curves(ic)%iedge(:) = 0
     end if

     ! compute polyline arclength
     if ( allocated(interdata_global%curves(ic)%polyline%s) ) deallocate(interdata_global%curves(ic)%polyline%s)
     allocate(interdata_global%curves(ic)%polyline%s(interdata_global%curves(ic)%polyline%np))
     interdata_global%curves(ic)%polyline%s(1) = 0._fp
     do i = 2,interdata_global%curves(ic)%polyline%np
        interdata_global%curves(ic)%polyline%s(i) = interdata_global%curves(ic)%polyline%s(i-1) + &
             norm2(interdata_global%curves(ic)%polyline%xyz(:,i) - interdata_global%curves(ic)%polyline%xyz(:,i-1))
     end do
  end do

end subroutine intersect_all_surfaces



subroutine write_intersection_data( &
     interdata, &
     filepoints, &
     filecurves )
  use mod_util
  use mod_types_intersection
  implicit none
  type(type_intersection_data), intent(in) :: interdata
  character(*),                 intent(in) :: filepoints, filecurves
  integer                                  :: fileunit
  integer                                  :: ip, ic, is

  call get_free_unit( fileunit )

  ! intersection points
  open( &
       unit = fileunit, &
       file = filepoints, &
       action = 'write' )
  do ip = 1,interdata%np
     write ( fileunit, * ) interdata%points(ip)%xyz
  end do
  close( fileunit )

  ! intersection curves
  open( &
       unit = fileunit, &
       file = filecurves, &
       action = 'write' )
  write ( fileunit, * ) interdata%nc
  do ic = 1,interdata%nc
     write ( fileunit, * ) logic2int(interdata%curves(ic)%dummy)
     write ( fileunit, * ) logic2int(interdata%curves(ic)%smooth)
     write ( fileunit, * ) interdata%curves(ic)%uvbox(:,:,1)
     write ( fileunit, * ) interdata%curves(ic)%uvbox(:,:,2)
     write ( fileunit, * ) interdata%curves(ic)%nsplit
     do ip = 1,interdata%curves(ic)%nsplit
        write ( fileunit, * ) interdata%curves(ic)%isplit(:,ip)
     end do
     do is = 1,interdata%curves(ic)%nsplit-1
        write ( fileunit, * ) 1!class(is)
     end do
     if ( associated(interdata%curves(ic)%polyline) ) then
        write ( fileunit, * ) interdata%curves(ic)%polyline%np
        do ip = 1,interdata%curves(ic)%polyline%np
           write ( fileunit, * ) interdata%curves(ic)%polyline%uv(:,:,ip), interdata%curves(ic)%polyline%xyz(:,ip)
        end do
     else
        write ( fileunit, * ) 0
     end if
  end do
  close( fileunit )    

end subroutine write_intersection_data



subroutine insert_polyline_point( &
    uv, &
    xyz, &
    stat, &
    polyline, &
    i )
  use mod_math
  use mod_types_intersection
  ! Inserts a uv-xyz point in an intersection_polyline after the i-th point.
  ! If 'i' is not provided, the point is inserted after the current last point.
  implicit none
  integer, parameter                              :: PARAM_xtra_np = 10
  real(kind=fp),                    intent(in)    :: uv(2,2)
  real(kind=fp),                    intent(in)    :: xyz(3)
  integer,                          intent(out)   :: stat
  type(type_intersection_polyline), intent(inout) :: polyline
  integer, optional,                intent(in)    :: i
  integer                                         :: iprev

  if ( present(i) ) then
     iprev = i
     iprev = min(iprev,polyline%np)
     iprev = max(iprev,0)
  else
     iprev = polyline%np
  end if

  if ( .not.allocated(polyline%uv)  ) allocate(polyline%uv(2,2,PARAM_xtra_np))
  if ( .not.allocated(polyline%xyz) ) allocate(polyline%xyz(3,PARAM_xtra_np) )

  if ( polyline%np + 1 > size(polyline%xyz,2) ) then
     call reallocate_polyline( &
          polyline, &
          polyline%np + PARAM_xtra_np, &
          stat )
     if ( stat > 0 ) return
  end if

  stat = 0
  if ( iprev < polyline%np ) then
     polyline%uv(:,:,iprev+2:polyline%np+1) = polyline%uv(:,:,iprev+1:polyline%np)
     polyline%xyz(:,iprev+2:polyline%np+1)  = polyline%xyz(:,iprev+1:polyline%np)
  end if
  polyline%uv(:,:,iprev+1) = uv
  polyline%xyz(:,iprev+1)  = xyz
  polyline%np = polyline%np + 1

end subroutine insert_polyline_point



subroutine transfer_intersection_curves( &
     from, &
     to )
  use mod_types_intersection
  implicit none
  type(type_intersection_data), intent(inout) :: from
  type(type_intersection_data), intent(inout) :: to
  integer                                     :: ic

  to%nc = from%nc

  do ic = 1,from%nc
     to%curves(ic)%dummy        =  from%curves(ic)%dummy
     to%curves(ic)%smooth       =  from%curves(ic)%smooth
     to%curves(ic)%surf(1)%ptr  => from%curves(ic)%surf(1)%ptr
     to%curves(ic)%surf(2)%ptr  => from%curves(ic)%surf(2)%ptr
     to%curves(ic)%uvbox        =  from%curves(ic)%uvbox
     to%curves(ic)%param_vector =  from%curves(ic)%param_vector
     to%curves(ic)%w0           =  from%curves(ic)%w0
     if ( allocated(from%curves(ic)%isplit) ) then
        call move_alloc(from=from%curves(ic)%isplit, to=to%curves(ic)%isplit)
     end if
     to%curves(ic)%nsplit       =  from%curves(ic)%nsplit
     if ( allocated(from%curves(ic)%iedge) ) then
        call move_alloc(from=from%curves(ic)%iedge, to=to%curves(ic)%iedge)
     end if
     to%curves(ic)%polyline     => from%curves(ic)%polyline
     nullify( &
          from%curves(ic)%surf(1)%ptr, &
          from%curves(ic)%surf(2)%ptr, &
          from%curves(ic)%polyline)
  end do

  deallocate(from%curves)
  from%nc = 0

end subroutine transfer_intersection_curves



subroutine merge_intersection_data( &
     surf, &
     uvxyz, &
     nuvxyz, &
     interdata_local, &
     interdata_global, &
     tolchord, &
     hmin, &
     hmax )
  USE MOD_UTIL
  use mod_math
  use mod_diffgeom
  use mod_types_intersection
  use mod_tolerances
  ! Trace all intersection curves, intersect them and subidivide them accordingly and 
  implicit none
  LOGICAL, PARAMETER :: DEBUG = ( GLOBALDEBUG .AND. .false. )
  type(ptr_surface),            intent(in)    :: surf(2)
  integer,                      intent(in)    :: nuvxyz
  real(kind=fp),                intent(in)    :: uvxyz(8,nuvxyz)
  type(type_intersection_data), intent(in)    :: interdata_local
  type(type_intersection_data), intent(inout) :: interdata_global
  real(kind=fp),                intent(in)    :: tolchord, hmin, hmax
  integer                                     :: id_global(nuvxyz)
  integer                                     :: stat
  integer                                     :: nc
  integer                                     :: ip, ic, jc, isurf

  IF ( DEBUG ) THEN
     PRINT *,'UVXYZ ='
     CALL PRINT_MAT(TRANSPOSE(UVXYZ(1:7,1:NUVXYZ)))
     PRINT *,'TOLUV =',UVXYZ(8,1:NUVXYZ)
  END IF
  ! add new intersection points
  do ip = 1,nuvxyz
     call add_intersection_point( &
          reshape(uvxyz(1:4,ip), [2,2]), &
          uvxyz(5:7,ip), &
          surf, &
          2, &
          interdata_global, &
          id_global(ip) ) 
  end do

  nc = interdata_global%nc ! number of curves before adding new ones
  do ic = 1,interdata_local%nc
     ! add curve
     call add_intersection_curve( &
          interdata_global, &
          interdata_local%curves(ic)%param_vector, &
          id_global(interdata_local%curves(ic)%isplit(1,1:2)), &
          interdata_local%curves(ic)%uvbox )
     do isurf = 1,2
        interdata_global%curves(interdata_global%nc)%surf(isurf)%ptr => surf(isurf)%ptr
     end do

     ! trace polyline
     IF ( DEBUG ) PRINT *,'TRACE POLYLINE...'
     allocate(interdata_global%curves(interdata_global%nc)%polyline)
     call trace_intersection_polyline( &
          surf, &
          interdata_global%curves(nc+ic)%uvbox, &
          interdata_global%curves(nc+ic)%param_vector, &
          reshape(uvxyz(1:4,interdata_local%curves(ic)%isplit(1,1:2)),[2,2,2]), &
          uvxyz(5:7,interdata_local%curves(ic)%isplit(1,1:2)), &
          stat, &
          interdata_global%curves(nc+ic)%polyline, &
          interdata_global%curves(nc+ic)%w0, &
          tolchord, &
          hmin, &
          hmax )
     IF ( DEBUG ) PRINT *,'...OK'
     IF ( interdata_global%curves(nc+ic)%polyline%np < 2 ) THEN
        PRINT *,'CURVE #', nc+ic,', POLYLINE%NP =', interdata_global%curves(nc+ic)%polyline%np
        PAUSE
     END IF
     if ( stat > 0 ) then
        PRINT *,'STAT_TRACE_INTERSECTION_POLYLINE = ',STAT
        PRINT *,'PARAM_VECTOR =',interdata_global%curves(nc+ic)%param_vector
        PAUSE
        return
     end if

     ! add the polyline's endpoints in isplit
     interdata_global%curves(nc+ic)%isplit(2,1) = 1
     interdata_global%curves(nc+ic)%isplit(2,2) = interdata_global%curves(nc+ic)%polyline%np

     ! intersection with other curves
     IF ( DEBUG ) PRINT *,'INTERSECT WITH OTHER CURVES...'
     do jc = 1,nc
        IF ( DEBUG ) PRINT *,'CURVE #',JC
        if ( interdata_global%curves(jc)%dummy ) cycle
        call intersect_intersection_curves( &
             interdata_global, &
             [nc+ic,jc], &
             stat )
        if ( stat > 0 ) then
           PRINT *,'STAT_INTERSECT_INTERSECTION_CURVES = ',STAT
           PAUSE
           return
        end if
     end do
     IF ( DEBUG ) PRINT *,'...OK'
  end do

end subroutine merge_intersection_data



subroutine newton_three_surfaces( &
     surf, &
     lowerb, &
     upperb, &
     stat, &
     uv, &
     xyz)
  use mod_math
  use mod_linalg
  use mod_diffgeom
  use mod_tolerances  
  implicit none
  LOGICAL, PARAMETER :: DEBUG = .false.!( GLOBALDEBUG .AND. .true. )
  integer,           parameter     :: itmax = 2 + ceiling(-log10(EPSuv))
  type(ptr_surface), intent(in)    :: surf(3)
  real(kind=fp),     intent(in)    :: lowerb(6)
  real(kind=fp),     intent(in)    :: upperb(6)
  integer,           intent(out)   :: stat
  real(kind=fp),     intent(inout) :: uv(2,3)
  real(kind=fp),     intent(out)   :: xyz(3)
  real(kind=fp)                    :: xyzs(3,3)
  real(kind=fp)                    :: r(6), jac(6,6), duv(6)
  real(kind=fp)                    :: cond, erruv
  integer                          :: it, isurf
  integer                          :: stat_refl
  
  IF ( DEBUG ) THEN
     PRINT *,''; PRINT *,'';
     PRINT *,'NEWTON_THREE_SURFACES'
     PRINT *,'LOWERB =',LOWERB
     PRINT *,'UPPERB =',UPPERB
  END IF
  stat = 1
  erruv = 0._fp
  cond = 1._fp

  do it = 1,itmax
     !IF ( DEBUG ) PRINT *,'UV =',UV
     !! compute residual
     do isurf = 1,3
        call eval( &
             xyzs(:,isurf), &
             surf(isurf)%ptr, &
             uv(:,isurf) )
     end do
     r(1:3) = xyzs(:,1) - xyzs(:,2)
     r(4:6) = xyzs(:,1) - xyzs(:,3)

     !! compute Jacobian matrix
     call evald1(jac(1:3,1), surf(1)%ptr, uv(:,1), 1)
     call evald1(jac(1:3,2), surf(1)%ptr, uv(:,1), 2)
     call evald1(jac(1:3,3), surf(2)%ptr, uv(:,2), 1)
     call evald1(jac(1:3,4), surf(2)%ptr, uv(:,2), 2)
     call evald1(jac(4:6,5), surf(3)%ptr, uv(:,3), 1)
     call evald1(jac(4:6,6), surf(3)%ptr, uv(:,3), 2)
     jac(4:6,1:2) = jac(1:3,1:2)
     jac(1:6,3:6) = -jac(1:6,3:6)
     jac(1:3,5:6) = 0._fp
     jac(4:6,3:4) = 0._fp

     IF ( .false. ) THEN
        PRINT *,'JAC ='
        CALL PRINT_MAT(JAC)

        PRINT *,'RHS =',-r
     END IF

     !! solve for Newton step
     call linsolve_svd( &
          duv, &
          jac, &
          -r, &
          6, &
          6, &
          1, &
          cond, &
          tol=0.1_fp*EPSuv )

     IF ( COND > EPSUV/EPSFP ) THEN
        if ( max(sum(r(1:3)**2), sum(r(4:6)**2)) < EPSxyzsqr ) then
           stat = 0
           erruv = 10._fp*EPSuv !*******
           IF ( DEBUG ) PRINT *,'newton_three_surfaces : /!\ singular Jacobian at solution'
           xyz = sum(xyzs, 2)/3._fp
           return
        end if
     END IF

     erruv = maxval([sum(duv(1:2)**2), sum(duv(3:4)**2), sum(duv(5:6))**2])

     !! correct Newton step to keep the iterate inside feasible region
     call nd_box_reflexions( &
          reshape(uv, [6]), &
          lowerb, &
          upperb, &
          duv, &
          6, &
          stat_refl )
     if ( stat_refl > 0 ) return

     !! update solution
     uv(:,1) = uv(:,1) + duv(1:2)
     uv(:,2) = uv(:,2) + duv(3:4)
     uv(:,3) = uv(:,3) + duv(5:6)

     IF ( DEBUG ) PRINT *,norm2(r(1:3)), norm2(r(4:6)), sqrt(erruv), EPSfp*cond
     !! termination criteria
     if ( erruv < max(EPSuvsqr, EPSfpsqr*cond**2) ) then
        if ( max(sum(r(1:3)**2), sum(r(4:6)**2)) < EPSxyzsqr ) then
           if ( erruv > EPSuvsqr ) then
              IF ( DEBUG ) PRINT *,'newton_three_surfaces : /!\ toluv > EPSuv'
              pause
           end if
           stat = 0
           do isurf = 1,3
              call eval( &
                   xyzs(:,isurf), &
                   surf(isurf)%ptr, &
                   uv(:,isurf) )
           end do
           xyz = sum(xyzs, 2)/3._fp
           IF ( DEBUG ) PRINT *,'CONVERGED, UV =',UV,', XYZ =',XYZ
        else
           IF ( DEBUG ) PRINT *,'STAGNATION'
        end if
        return
     end if

  end do

end subroutine newton_three_surfaces



subroutine copy_tangent_intersections( &
     from, &
     to )
  use mod_types_intersection
  implicit none
  type(type_intersection_data), intent(in), target :: from
  type(type_intersection_data), intent(inout)      :: to
  type(type_point_on_surface), pointer             :: pos => null()
  real(kind=fp), allocatable                       :: uv(:,:)
  type(ptr_surface), allocatable                   :: surf(:)
  integer                                          :: icurv, j, ipoint, iendpoints(2), isurf, np

  do icurv = 1,from%nc
     if ( .not.from%curves(icurv)%smooth ) cycle
     
     do j = 1,2
        !if ( j == 1 ) then
        !   ipoint = from%curves(icurv)%isplit(1,1)
        !else
        !   ipoint = from%curves(icurv)%isplit(1,from%curves(icurv)%nsplit)
        !end if
        ipoint = from%curves(icurv)%isplit(1,1 + (j-1)*(from%curves(icurv)%nsplit-1))

        allocate(uv(2,from%points(ipoint)%npos), surf(from%points(ipoint)%npos))
        pos => from%points(ipoint)%pos
        do isurf = 1,from%points(ipoint)%npos
           uv(:,isurf) = pos%uv 
           surf(isurf)%ptr => pos%surf
           pos => pos%next
        end do        
        
        call add_intersection_point( &
             uv(1:2,1:from%points(ipoint)%npos), &
             from%points(ipoint)%xyz, &
             surf(1:from%points(ipoint)%npos), &
             from%points(ipoint)%npos, &
             to, &
             iendpoints(j) )
        deallocate(uv, surf)
        
     end do

     call add_intersection_curve( &
          to, &
          from%curves(icurv)%param_vector, &
          iendpoints, &
          from%curves(icurv)%uvbox )

     to%curves(to%nc)%dummy = from%curves(icurv)%dummy
     to%curves(to%nc)%smooth = from%curves(icurv)%smooth
     do isurf = 1,2
        to%curves(to%nc)%surf(isurf)%ptr => from%curves(icurv)%surf(isurf)%ptr
     end do
     !to%curves(to%nc)%polyline => from%curves(icurv)%polyline
     allocate(to%curves(to%nc)%polyline)
     np = from%curves(icurv)%polyline%np
     to%curves(to%nc)%polyline%np = np
     allocate(to%curves(to%nc)%polyline%uv(2,2,np), to%curves(to%nc)%polyline%xyz(3,np))
     to%curves(to%nc)%polyline%uv(1:2,1:2,1:np) = from%curves(icurv)%polyline%uv(1:2,1:2,1:np)
     to%curves(to%nc)%polyline%xyz(1:3,1:np) = from%curves(icurv)%polyline%xyz(1:3,1:np)
     
     if ( allocated(to%curves(to%nc)%iedge) ) deallocate(to%curves(to%nc)%iedge)
     to%curves(to%nc)%isplit(2,1:2) = [1,np]!from%curves(icurv)%isplit(2,[1,from%curves(icurv)%nsplit])
     
  end do
  
end subroutine copy_tangent_intersections



subroutine check_unicity_with_tolerance( &
     unew, &
     xnew, &
     tolunew, &
     dim, &
     u, &
     x, &
     tolu, &
     n, &
     id )
  use mod_math
  use mod_tolerances
  implicit none
  integer,       intent(in)  :: dim
  real(kind=fp), intent(in)  :: unew(dim)
  real(kind=fp), intent(in)  :: xnew(3)
  real(kind=fp), intent(in)  :: tolunew
  integer,       intent(in)  :: n
  real(kind=fp), intent(in)  :: u(dim,n)
  real(kind=fp), intent(in)  :: x(3,n)
  real(kind=fp), intent(in)  :: tolu(n)
  integer,       intent(out) :: id

  if ( n < 1 ) then
     id = 1
  else
     do id = 1,n
        if ( maxval(abs(unew - u(:,id))) < tolunew + tolu(id) .or. & 
           !sum((unew - u(:,id))**2) < (tolunew + tolu(id))**2 .or. &
             sum((xnew - x(:,id))**2) < EPSxyzsqr ) then
           exit
        end if
     end do
  end if
 
end subroutine check_unicity_with_tolerance



subroutine rearrange_for_separability_test( &
     bcp, &
     nbcp, &
     vec, &
     sep, &
     nsep, &
     gaussmap )
  use mod_math
  use mod_tolerances
  implicit none
  integer,       intent(in)  :: nbcp
  real(kind=fp), intent(in)  :: bcp(nbcp,3)
  real(kind=fp), intent(in)  :: vec(3)
  logical,       intent(in)  :: gaussmap
  real(kind=fp), intent(out) :: sep(nbcp,3)
  integer,       intent(out) :: nsep
  real(kind=fp)              :: tmp(3), tmpsqr
  integer                    :: i

  if ( gaussmap ) tmpsqr = sum(vec**2)
  nsep = 0
  do i = 1,nbcp

     if ( gaussmap ) then
        tmp = bcp(i,:) / norm2(bcp(i,:))
        if ( tmpsqr > EPSxyzsqr .and. abs(dot_product(tmp, vec)) >= 1._fp - EPSxyz ) cycle
     else
        tmp = bcp(i,:) - vec
        tmpsqr = sum( tmp**2 )
        if ( tmpsqr <= 400.D0 * EPSxyzsqr ) cycle
        tmp = tmp / sqrt(tmpsqr)
     end if

     nsep = nsep + 1
     sep(nsep,:) = tmp

  end do

end subroutine rearrange_for_separability_test



subroutine newton_curve_surface( &
     curv, &
     surf, &
     lowerb, &
     upperb, &
     stat, &
     tuv, &
     toltuv, &
     xyz )
  use mod_math
  use mod_linalg
  use mod_diffgeom
  use mod_tolerances    
  ! stat = 0 : converged
  !        1 : not converged
  !        2 : degeneracy
  implicit none
  LOGICAL, PARAMETER :: DEBUG = ( GLOBALDEBUG .AND. .TRUE. )
  logical,       parameter          :: acceleration = .false.
  real(kind=fp), parameter          :: THRESHOLD = real(1.d-2, kind=fp)**2
  integer,       parameter          :: itmax = 30!2 + ceiling(-log10(EPSuv))
  integer,       parameter          :: itconv = 5
  type(type_curve),   intent(in)    :: curv
  type(type_surface), intent(in)    :: surf
  real(kind=fp),      intent(in)    :: lowerb(3)
  real(kind=fp),      intent(in)    :: upperb(3)
  integer,            intent(out)   :: stat
  real(kind=fp),      intent(inout) :: tuv(3)
  real(kind=fp),      intent(out)   :: toltuv
  real(kind=fp),      intent(out)   :: xyz(3)
  real(kind=fp), dimension(3)       :: xyz_c, xyz_s, r
  real(kind=fp)                     :: resxyz, resconv
  real(kind=fp)                     :: jac(3,3), dtuv(3), cond, errtuv, errtuvprev
  integer                           :: it, rank
  logical                           :: linear_conv
  real(kind=fp)                     :: fac
  integer                           :: stat_refl

  IF ( DEBUG ) THEN
     PRINT *,''; PRINT *,'';
     PRINT *,'NEWTON_CURVE_SURFACE'
     PRINT *,'LOWERB =',LOWERB
     PRINT *,'UPPERB =',UPPERB
  END IF

  stat = 1
  errtuv = 0._fp
  errtuvprev = 0._fp
  cond = 1._fp
  fac = 1._fp
  linear_conv = .false.
  !IF ( DEBUG ) PRINT *,'|XS - XC|, |DTUV|, EPS*COND(J)'
  do it = 1,itmax
     !! compute residual
     call eval(xyz_c, curv, tuv(1)  ) ! curve's position vector
     call eval(xyz_s, surf, tuv(2:3)) ! surface's position vector

     r = xyz_s - xyz_c

     ! check signs of convergence
     resxyz = sum(r**2)
     IF ( DEBUG ) PRINT *,SQRT(RESXYZ), SQRT(ERRTUV), EPSFP*COND
     if ( it == 1 ) resconv = THRESHOLD * resxyz
     if ( it > itconv ) then
        if ( resxyz > resconv ) then
           ! Newton sequence not likely to converge, presumably no solution
           IF ( DEBUG ) PRINT *,'NO SIGN OF CONVERGENCE'
           return
        elseif ( errtuv > THRESHOLD * errtuvprev ) then
           ! linear convergence
           IF ( DEBUG ) PRINT *,'LINEAR CONVERGENCE, MU =',sqrt(errtuv/errtuvprev)
           linear_conv = .true.
           fac = 2._fp
        else
           ! superlinear convergence
           fac = 1._fp
           linear_conv = .false.
        end if
     end if

     !! compute Jacobian matrix
     call evald1(jac(:,1), curv, tuv(1))
     jac(:,1) = -jac(:,1)
     call evald1(jac(:,2), surf, tuv(2:3), 1)
     call evald1(jac(:,3), surf, tuv(2:3), 2)
     IF ( DEBUG ) THEN
        PRINT *,'JAC ='
        CALL PRINT_MAT(JAC)
        PRINT *,'RHS =',-R
     END IF

     !! solve for Newton step
     call linsolve_svd( &
          dtuv, &
          jac, &
          -r, &
          3, &
          3, &
          1, &
          cond, &
          rank, &
          tol=EPSmath )
     !PRINT *,'RANK =',RANK
     if ( acceleration ) then
        !PRINT *,'FAC =',FAC
        if ( mod(it,3) == 0 ) dtuv = fac * dtuv
     end if
     errtuvprev = errtuv
     errtuv = sum(dtuv**2)

     ! correct Newton step to keep the iterate inside feasible region
     IF ( DEBUG ) THEN
        PRINT *,' TUV =',TUV
        PRINT *,'DTUV =',DTUV
     END IF
     if ( .true. ) then
        if ( errtuv > sum((upperb - lowerb)**2) ) then
           dtuv = sqrt(sum((upperb - lowerb)**2) / errtuv) * dtuv
           IF ( DEBUG ) PRINT *,'DTUV*=',DTUV
        end if
     end if
     call nd_box_reflexions( &
          tuv, &
          lowerb, &
          upperb, &
          dtuv, &
          3, &
          stat_refl )
     if ( stat_refl > 0 ) return
     IF ( .false. ) PRINT *,'...OK'

     ! update solution
     tuv = tuv + dtuv

     !! termination criterion
     if ( errtuv < max(EPSuvsqr, EPSfpsqr*cond**2) ) then
        if ( resxyz < EPSxyzsqr ) then
           if ( errtuv > EPSuvsqr ) then
              IF ( DEBUG ) PRINT *,'newton_curve_surface : /!\ toltuv > EPSuv'
              if ( linear_conv ) then
                 IF ( DEBUG ) THEN
                    PRINT *,'--> NEWTON_CS_DYNSHIFT'
                    PRINT *,'TUV =',TUV
                 END IF
                 call newton_curve_surface_dynshift( &
                      curv, &
                      surf, &
                      lowerb, &
                      upperb, &
                      stat, &
                      tuv, &
                      toltuv, &
                      xyz )
                 return
              end if
           end if
           ! converged to a curve-surface intersection point
           stat = 0
           toltuv = max(EPSmath, sqrt(errtuv))
           call eval(xyz_c, curv, tuv(1)  )
           call eval(xyz_s, surf, tuv(2:3))
           xyz = 0.5_fp * (xyz_c + xyz_s)
           IF ( DEBUG ) PRINT *,SQRT(RESXYZ), SQRT(ERRTUV), EPSFP*COND
           IF ( DEBUG ) PRINT *,'CONVERGED, TUV =',TUV,', XYZ =',XYZ
           !PAUSE
        else
           IF ( DEBUG ) PRINT *,'STAGNATION'
        end if
        return
     end if

  end do

end subroutine newton_curve_surface



subroutine find_collineal_points( &
     surf, &
     lowerb, &
     upperb, &
     stat, &
     uv_collineal, &
     n_collineal, &
     toluv, &
     xyz_collineal )
  use mod_math
  use mod_linalg
  use mod_diffgeom
  use mod_tolerances
  ! Searches for a pair of collineal points on two rectangular parametric surfaces
  ! using a box-constrained Newton-Raphson algorithm. 
  ! The lower (resp. upper) bounds of the 4-dimensional feasible domain are stored 
  ! in 'lowerb' (resp. 'upperb').
  ! Returns: stat > 0 if no such pair has been found;
  !               =-1 if a tangential contact point has been found;
  !               = 0 else (pair of non-coincident collineal points).
  !          uv_collineal  : uv-coordinates of the collineal points
  !          xyz_collineal : xyz-coordinates (relevent only if stat =-1)
  !          n_collineal   : the common (unit) normal direction at the collineal points.
  implicit none
  LOGICAL, PARAMETER :: DEBUG = ( GLOBALDEBUG .AND. .false. )
  integer,           parameter     :: itmax = 2 + ceiling(-log10(EPScollineal))
  type(ptr_surface), intent(in)    :: surf(2)
  integer,           intent(out)   :: stat
  real(kind=fp),     intent(inout) :: uv_collineal(2,2)
  real(kind=fp),     intent(out)   :: n_collineal(3)
  real(kind=fp),     intent(out)   :: toluv
  real(kind=fp),     intent(out)   :: xyz_collineal(3)
  real(kind=fp),     intent(in)    :: lowerb(4)
  real(kind=fp),     intent(in)    :: upperb(4)
  real(kind=fp)                    :: xyz(3,2), dxyz_duv(3,2,2), d2xyz_duv2(3,3,2), n(3), r(3)
  real(kind=fp)                    :: f(4), jac(4,4), duv(4)
  real(kind=fp)                    :: cond, erruv, lambda
  integer                          :: it, isurf, ivar
  integer                           :: stat_refl

  IF ( DEBUG ) THEN
     PRINT *,''; PRINT *,'';
     PRINT *,'FIND_COLLINEAL_POINTS'
     PRINT *,'LOWERB =',LOWERB
     PRINT *,'UPPERB =',UPPERB
  END IF

  stat = 1
  erruv = 0._fp
  cond = 1._fp 
  
  do it = 1,itmax 
     !IF ( DEBUG ) PRINT *,'IT #',IT
     !! compute residual
     do isurf = 1,2 ! <--------------------------+
        ! position vector                        !
        call eval( &                             !
             xyz(:,isurf), &                     !
             surf(isurf)%ptr, &                  !
             uv_collineal(:,isurf) )             !
        !                                        !
        ! tangent vectors                        !
        do ivar = 1,2 ! <--------------------+   !
           call evald1( &                    !   !
                dxyz_duv(:,ivar,isurf), &    !   !
                surf(isurf)%ptr, &           !   !
                uv_collineal(:,isurf), &     !   !
                ivar )                       !   !
        end do ! <---------------------------+   !
     end do ! <----------------------------------+
     
     n = cross( dxyz_duv(:,1,1), dxyz_duv(:,2,1) ) ! (pseudo-)normal to surface 1
     r = xyz(:,1) - xyz(:,2)

     do ivar = 1,2 ! <-------------------------------------+
        f(ivar)   = dot_product( dxyz_duv(:,ivar,2), n )   !
        f(2+ivar) = dot_product( dxyz_duv(:,ivar,1), r )   !
     end do ! <--------------------------------------------+

     !IF ( DEBUG ) PRINT *,'F =',F
     IF ( DEBUG ) PRINT *, NORM2(R), NORM2(F), SQRT(ERRUV), EPSFP*COND

     !! compute Jacobian matrix
     do isurf = 1,2 ! <-----------------------------+
        do ivar = 1,3 ! <-----------------------+   !
           call evald2( &                       !   !
                d2xyz_duv2(:,ivar,isurf), &     !   !
                surf(isurf)%ptr, &              !   !
                uv_collineal(:,isurf), &        !   !
                ivar )                          !   !
        end do ! <------------------------------+   !
     end do ! <-------------------------------------+

     do ivar = 1,2
        jac(1:2,ivar) = matmul( &
             transpose(dxyz_duv(:,:,2)), &
             cross( d2xyz_duv2(:,ivar,1), dxyz_duv(:,2,1) ) + &
             cross( dxyz_duv(:,1,1), d2xyz_duv2(:,ivar+1,1) ) &
             )
        jac(3:4,ivar) = matmul( &
             transpose(d2xyz_duv2(:,ivar:ivar+1,1)), r ) + &
             matmul( transpose(dxyz_duv(:,:,1)), dxyz_duv(:,ivar,1) &
             )
        jac(1:2,2+ivar) = matmul( transpose(d2xyz_duv2(:,ivar:ivar+1,2)), n )
        jac(3:4,2+ivar) = -matmul( transpose(dxyz_duv(:,:,1)), dxyz_duv(:,ivar,2) )
     end do
     

     !! solve for Newton step
     !PRINT *,'FIND_COLLINEAL_POINTS, IT#',IT
     !CALL PRINT_MAT(JAC)
     call linsolve_svd( &
          duv, &
          jac, &
          -f, &
          4, &
          4, &
          1, &
          cond )
     erruv = max(sum(duv(1:2)**2), sum(duv(3:4)**2))
     IF (.true.) THEN
        ! (seems to be faster)
        ! scale down Newton step to keep the solution inside feasible region
        call nd_box_constraint( &
             reshape( uv_collineal, [4] ), &
             lowerb, &
             upperb, &
             duv, &
             lambda )
        if ( lambda < EPSfp ) then ! <---+
           ! non-positive scaling factor !
           return                        !
        end if ! <-----------------------+
        duv = lambda * duv
     ELSE
        ! correct Newton step to keep the iterate inside feasible region
        call nd_box_reflexions( &
             reshape( uv_collineal, [4] ), &
             lowerb, &
             upperb, &
             duv, &
             4, &
             stat_refl )
        if ( stat_refl > 0 ) return
     END IF
     ! update solution
     uv_collineal(:,1) = uv_collineal(:,1) + duv(1:2)
     uv_collineal(:,2) = uv_collineal(:,2) + duv(3:4)
     
     !! termination criteria
     if ( erruv < max(EPSuvsqr, EPSfpsqr*cond**2) ) then ! <------------+
        if ( sum(f**2) < EPScollinealsqr ) then ! <----------------+    !
           IF ( DEBUG ) PRINT *, NORM2(R), NORM2(F), SQRT(ERRUV), EPSFP*COND
           if ( erruv > EPSuvsqr ) then
              PRINT *,'find_collineal_points : /!\ toluv > EPSuv'
           end if
           toluv = sqrt(erruv)                                     !    !
           if ( sum(r**2) < EPSxyzsqr ) then ! <--------------+    !    !
              ! converged to a tangential intersection point  !    !    !
              stat = -1                                       !    !    !
              do isurf = 1,2 ! <----------------+             !    !    !
                 call eval( &                   !             !    !    !
                      xyz(:,isurf), &           !             !    !    !
                      surf(isurf)%ptr, &        !             !    !    !
                      uv_collineal(:,isurf) )   !             !    !    !
              end do ! <------------------------+             !    !    !
              xyz_collineal = 0.5_fp * sum(xyz, 2)            !    !    !
           else ! --------------------------------------------+    !    !
              ! converged to a pair of collineal points       !    !    !
              stat = 0                                        !    !    !
           end if ! <-----------------------------------------+    !    !
           n_collineal = n / norm2( n )                            !    !
        else ! ----------------------------------------------------+    !
           IF ( DEBUG ) PRINT *,'STAGNATION'                       !    !
        end if ! <-------------------------------------------------+    !
        return                                                          !
     end if ! <---------------------------------------------------------+

  end do

end subroutine find_collineal_points



subroutine newton_three_surfaces_1tangential( &
     surf, &
     lowerb, &
     upperb, &
     stat, &
     uv, &
     xyz )
  use mod_math
  use mod_linalg
  use mod_diffgeom
  use mod_tolerances
  USE MOD_UTIL
  USE MOD_POLYNOMIAL
  implicit none
  LOGICAL, PARAMETER :: DEBUG = .false.!( GLOBALDEBUG .AND. .true. )
  integer,           parameter     :: itmax = 2 + ceiling(-log10(EPSuv))
  type(ptr_surface), intent(in)    :: surf(3)
  real(kind=fp),     intent(in)    :: lowerb(6)
  real(kind=fp),     intent(in)    :: upperb(6)
  integer,           intent(out)   :: stat
  real(kind=fp),     intent(inout) :: uv(2,3)
  real(kind=fp),     intent(out)   :: xyz(3)
  real(kind=fp)                    :: xyzs(3,3), r(6)
  type(ptr_surface)                :: surftng(2)
  real(kind=fp)                    :: uvtng(2,2)
  real(kind=fp)                    :: duv_ds(2,2,2), dxyz_ds(3,2), curvature(2)
  real(kind=fp)                    :: mat(3,3), dtuv(3), duv(6)
  integer                          :: stat_contactpoint
  real(kind=fp)                    :: cond, erruv
  integer                          :: it, isurf, ivar
  integer                          :: stat_refl
  
  IF ( DEBUG ) THEN
     PRINT *,''
     PRINT *,'NEWTON_THREE_SURFACES_1TANGENTIAL'
     PRINT *,'UV0 ='
     CALL PRINT_MAT(transpose(UV))
  END IF
  
  stat = 1
  erruv = 0._fp
  cond = 1._fp

  surftng(1)%ptr => surf(1)%ptr
  surftng(2)%ptr => surf(3)%ptr

  IF ( DEBUG ) PRINT *,'|X1 - X2| , |X1 - X3| , |DUV| , eps*cond'
  do it = 1,itmax
     IF ( DEBUG ) PRINT *,'IT. #',IT
     uvtng = uv(:,[1,3])
     !! relax tangential intersection
     call simultaneous_point_inversions( &
          surftng, &
          lowerb([1,2,5,6]) - 100._fp*EPSuv, &
          upperb([1,2,5,6]) + 100._fp*EPSuv, &
          stat, &
          uvtng, &
          xyz )
     if ( stat > 0 ) then
        CALL WRITE_POLYNOMIAL(surf(1)%ptr%x, '../debug/debug_n3s1t_surf1.cheb')
        CALL WRITE_POLYNOMIAL(surf(2)%ptr%x, '../debug/debug_n3s1t_surf2.cheb')
        CALL WRITE_POLYNOMIAL(surf(3)%ptr%x, '../debug/debug_n3s1t_surf3.cheb')
        CALL write_matrix(uv, 2, 3, '../debug/debug_n3s1t_uv.dat')
        STOP 'newton_three_surfaces_1tangential : failed to relax onto tangential intersection :('
     end if
     uv(:,[1,3]) = uvtng
     IF ( DEBUG ) THEN
        PRINT *,'UV='
        CALL PRINT_MAT(transpose(UV))
     END IF
     
     !! compute residual
     do isurf = 1,3
        call eval( &
             xyzs(:,isurf), &
             surf(isurf)%ptr, &
             uv(:,isurf) )
     end do
     IF ( DEBUG ) THEN
        PRINT *,'XYZ='
        CALL PRINT_MAT(transpose(XYZS))
     END IF
     
     r(1:3) = xyzs(:,1) - xyzs(:,2)
     r(4:6) = xyzs(:,1) - xyzs(:,3)
     IF ( .true. ) THEN
        if ( max(sum(r(1:3)**2), sum(r(4:6)**2)) < EPSxyzsqr ) then
           stat = 0
           xyz = sum(xyzs, 2)/3._fp
           IF ( DEBUG ) PRINT *,'CONVERGED, UV =',UV,', XYZ =',XYZ
           return
        end if
     END IF

     !! compute derivatives
     if ( sum(r(4:6)**2) > EPSxyzsqr ) then
        PRINT *,'|X1 - X3| = ',norm2(r(4:6))
        STOP 'newton_three_surfaces_1tangential : not an intersection point :('
     else
        ! tangent to tangential intersection curve
        call diffgeom_intersection( &
             surftng, &
             uv(:,[1,3]), &
             duv_ds, &
             dxyz_ds, &
             stat_contactpoint, &
             curvature )
        !IF ( DEBUG ) PRINT *,'STAT_CONTACTPOINT =',stat_contactpoint
        IF ( DEBUG ) PRINT *,'DUV_DS  =',duv_ds(:,1,:)
        !IF ( DEBUG ) PRINT *,'DXYZ_DS =',dxyz_ds(:,1)
        if ( stat_contactpoint /= 1 ) then
           PRINT *,'STAT_CONTACTPOINT =',STAT_CONTACTPOINT
           STOP 'newton_three_surfaces_1tangential : not a tangential intersection curve'
        end if
     end if
     mat(:,1) = -dxyz_ds(:,1)

     ! tangents to surface
     do ivar = 1,2
        call evald1( &
             mat(:,ivar+1), &
             surf(2)%ptr, &
             uv(:,2), &
             ivar )
     end do

     !! solve for Newton step
     IF ( DEBUG ) PRINT *,'MAT ='
     IF ( DEBUG ) CALL PRINT_MAT(MAT)
     IF ( DEBUG ) PRINT *,'RHS =',0.5_fp*(xyzs(:,1) + xyzs(:,3)) - xyzs(:,2)
     call linsolve_svd( &
          dtuv, &
          mat, &
          0.5_fp*(xyzs(:,1) + xyzs(:,3)) - xyzs(:,2), &
          3, &
          3, &
          1, &
          cond, &
          tol=EPSmath )
     PRINT *,'DTUV =',DTUV

     duv(1:2) = dtuv(1)*duv_ds(:,1,1)
     duv(3:4) = dtuv(2:3)
     duv(5:6) = dtuv(1)*duv_ds(:,1,2)
     
     erruv = sum(duv**2) / 3._fp

     !! correct Newton step to keep the iterate inside feasible region
     call nd_box_reflexions( &
          reshape(uv, [6]), &
          lowerb, &
          upperb, &
          duv, &
          6, &
          stat_refl )
     if ( stat_refl > 0 ) return
     !IF ( DEBUG ) PRINT *,'DUV =',DUV
     IF ( DEBUG ) THEN
        PRINT *,'DUV ='
        DO IVAR = 1,SIZE(DUV)
           PRINT *,DUV(IVAR)
        END DO
     END IF

     !! update solution
     uv(:,1) = uv(:,1) + duv(1:2)
     uv(:,2) = uv(:,2) + duv(3:4)
     uv(:,3) = uv(:,3) + duv(5:6)

     IF ( DEBUG ) PRINT *,norm2(r(1:3)), norm2(r(4:6)), sqrt(erruv), EPSfp*cond

     !! termination criteria
     if ( erruv < max(EPSuvsqr, EPSfpsqr*cond**2) ) then
        if ( max(sum(r(1:3)**2), sum(r(4:6)**2)) < EPSxyzsqr ) then
           stat = 0
           do isurf = 1,3
              call eval( &
                   xyzs(:,isurf), &
                   surf(isurf)%ptr, &
                   uv(:,isurf) )
           end do
           xyz = sum(xyzs, 2)/3._fp
           IF ( DEBUG ) PRINT *,'CONVERGED, UV =',UV,', XYZ =',XYZ
        else
           IF ( DEBUG ) PRINT *,'STAGNATION'
        end if
        return
     end if
  end do


end subroutine newton_three_surfaces_1tangential



subroutine trace_intersection_polyline( &
     surf, &
     uvbox, &
     param_vector, &
     uv_endpoints, &
     xyz_endpoints, &
     stat, &
     polyline, &
     w0, &
     tolchord, &
     hmin, &
     hmax )
  use mod_math
  use mod_diffgeom
  use mod_types_intersection
  use mod_tolerances
  ! Returns: stat =-1 if one surface is singular at a polyline point (normal = zero vector)
  !               = 0 if everything went OK :)
  !               = 1 if an error ocurred when reallocating polyline%uv
  !               = 2 if an error ocurred when reallocating polyline%xyz
  !               = 3 if an isolated tangential contact point was encountered
  !               = 4 if an high-order tangential contact point was encountered
  !               = 5 if backtracking did not terminate successfully
  implicit none
  LOGICAL, PARAMETER :: DEBUG = ( GLOBALDEBUG .AND. .true. )
  real(kind=fp), parameter                        :: FRACbacktrack = 0.5_fp
  real(kind=fp), parameter                        :: EPSbacktrack = real(1d-2, kind=fp)
  type(ptr_surface),                intent(in)    :: surf(2)
  real(kind=fp),                    intent(in)    :: uvbox(2,2,2)
  real(kind=fp),                    intent(inout) :: param_vector(3)
  real(kind=fp),                    intent(in)    :: uv_endpoints(2,2,2) ! u/v, #surf, #point
  real(kind=fp),                    intent(in)    :: xyz_endpoints(3,2)
  integer,                          intent(out)   :: stat
  type(type_intersection_polyline), intent(inout) :: polyline
  real(kind=fp),                    intent(out)   :: w0
  real(kind=fp),                    intent(in)    :: tolchord
  real(kind=fp),                    intent(in)    :: hmin, hmax
  real(kind=fp)                                   :: FRACcurvature_radius 
  integer                                         :: stat_tangent, stat_insertion, stat_newton
  real(kind=fp), dimension(4)                     :: lowerb, upperb
  real(kind=fp)                                   :: Dw, w, wprev, dist_from_end
  real(kind=fp)                                   :: duv_ds(2,2,2), dxyz_ds(3,2), curvature(2)
  real(kind=fp)                                   :: h_endpoints(2), h, EPSh, hnext
  real(kind=fp)                                   :: uv(2,2), xyz(3), lambda
  integer                                         :: ipt

  FRACcurvature_radius = 2._fp*sqrt(tolchord*(2._fp - tolchord))
  
  stat = 0
  lowerb = reshape(uvbox(1,1:2,1:2), [4]) - EPSuv
  upperb = reshape(uvbox(2,1:2,1:2), [4]) + EPSuv
  IF ( DEBUG ) THEN
     PRINT *,''; PRINT *,'';
     PRINT *,'TRACE_INTERSECTION_POLYLINE'
     PRINT *,' PARAM_VECTOR =', PARAM_VECTOR
     PRINT *,' UV_ENDPOINTS ='; PRINT *,uv_endpoints(:,:,1); PRINT *,uv_endpoints(:,:,2)
     PRINT *,'XYZ_ENDPOINTS ='; PRINT *,xyz_endpoints(:,1);  PRINT *,xyz_endpoints(:,2)
     PRINT *,'LOWERB =',LOWERB; PRINT *,'UPPERB =',UPPERB
  END IF
  w0 = dot_product(param_vector, xyz_endpoints(:,1))
  Dw = dot_product(param_vector, xyz_endpoints(:,2)) - w0
  param_vector = param_vector/Dw
  w0 = w0/Dw

  !PRINT *,'CURVATURE RADIUS AT ENDPOINTS ='
  do ipt = 2,1,-1
     ! compute tangent direction and curvature of the intersection curve at endpoints
     ! (in reverse order, so data at first endpoint is kept in memory)
     call diffgeom_intersection( &
          surf, &
          uv_endpoints(:,:,ipt), &
          duv_ds, &
          dxyz_ds, &
          stat_tangent, &
          curvature )
     !call diffgeom_intersection_curve( &
     !     surf, &
     !     uv_endpoints(:,:,ipt), &
     !     duv_ds, &
     !     dxyz_ds, &
     !     stat_tangent, &
     !     curvature )
     IF ( DEBUG ) PRINT *,'CURVATURE RADIUS = ',1._FP/CURVATURE(1)
     if ( stat_tangent < 0 .or. stat_tangent > 2 ) then ! <----+
        IF ( DEBUG ) PRINT *,'stat_tangent = ',stat_tangent    !
        stat = stat_tangent                                    !
        return                                                 !
     end if ! <------------------------------------------------+
     curvature(1) = max(EPSfp, curvature(1))
     h_endpoints(ipt) = FRACcurvature_radius / curvature(1)
  end do
  h_endpoints = max(h_endpoints, hmin) 
  h_endpoints = min(h_endpoints, hmax)

  
  ! first point
  call insert_polyline_point( &
       uv_endpoints(1:2,1:2,1), &
       xyz_endpoints(1:3,1), &
       stat_insertion, &
       polyline, &
       i=0 )
  if ( stat_insertion > 0 ) then ! <--------+
     stat = stat_insertion                  !
     return                                 !
  end if ! <--------------------------------+

  wprev = 0._fp
  h = h_endpoints(1)
  outer : do ! <-----------------------------------------------------------------+
     IF ( DEBUG ) PRINT *,'HTARGET =',H
     ! check if the current is close enough to the end of the polyline           !
     dist_from_end = sum((polyline%xyz(:,polyline%np) - xyz_endpoints(:,2))**2)  !
     if ( dist_from_end < h**2 ) then ! <-----------------------+                !
        if ( dist_from_end < h_endpoints(2)**2 ) then ! <----+  !                !
           if ( dist_from_end < (TOLh * min(h, h_endpoints(2)))**2 ) polyline%np = polyline%np - 1
           exit outer                                        !  !                !
        else ! ----------------------------------------------+  !                !
           IF ( DEBUG ) PRINT *,'SHORTEN H : H =',h,', DIST =',sqrt(dist_from_end),', HEND =',h_endpoints(2)
           h = sqrt(dist_from_end)                           !  !                !
        end if ! <-------------------------------------------+  !                !
     end if ! <-------------------------------------------------+                !
     !                                                                           !
     ! compute next point                                                        !
     EPSh = EPSbacktrack * h                                                     !
     inner : do ! <-----------------------------------------------------------+  !
        ! set initial iterate                                                 !  !
        call nd_box_constraint( &                                             !  !
             reshape(polyline%uv(:,:,polyline%np), [4]), &                    !  !
             lowerb, &                                                        !  !
             upperb, &                                                        !  !
             h*reshape(duv_ds(:,1,:), [4]), &                                 !  !
             lambda )                                                         !  !
        !                                                                     !  !
        IF (DEBUG .AND. LAMBDA < 1._FP) THEN
           PRINT *,' UV =',polyline%uv(:,:,polyline%np)
           PRINT *,'DUV =',h*duv_ds(:,1,:)
           PRINT *,'LAMBDA =',LAMBDA
        END IF
        uv = polyline%uv(:,:,polyline%np) + lambda*h*duv_ds(:,1,:)            !  !
        !                                                                     !  !
        ! refine using Newton-Raphson algorithm                               !  !
        IF ( DEBUG ) PRINT *,'UV0 =',UV
        call newton_intersection_polyline( &                                  !  !
             surf, &                                                          !  !
             lowerb, &                                                        !  !
             upperb, &                                                        !  !
             polyline%xyz(1:3,polyline%np), &                                 !  !
             (lambda*h)**2, &                                                 !  !
             stat_newton, &                                                   !  !
             uv, &                                                            !  !
             xyz )                                                            !  !
        !                                                                     !  !
        IF ( DEBUG ) PRINT *,'STAT_NEWTON =',STAT_NEWTON
        if ( stat_newton == 0 ) then ! <----------------------------------+   !  !
           ! Newton has converged, check if w is monotonic                !   !  !
           w = dot_product(param_vector, xyz) - w0                        !   !  !
           if ( is_in_open_interval(w, wprev, 1._fp) ) then ! <--------+  !   !  !
              ! compute tangent direction and curvature                !  !   !  !
              ! of the intersection curve at the current point         !  !   !  !
              call diffgeom_intersection( &                            !  !   !  !
                   surf, &                                             !  !   !  !
                   uv, &                                               !  !   !  !
                   duv_ds, &                                           !  !   !  !
                   dxyz_ds, &                                          !  !   !  !
                   stat_tangent, &                                     !  !   !  !
                   curvature )                                         !  !   !  !
              IF ( DEBUG ) PRINT *,'CURVATURE RADIUS = ',1._FP/CURVATURE(1)
              if ( stat_tangent < 0 ) then ! <--------------------+    !  !   !  !
                 stat = stat_tangent                              !    !  !   !  !
                 return ! singular surface (undefined normal)     !    !  !   !  !
              elseif ( stat > 2 ) then ! -------------------------+    !  !   !  !
                 stat = stat_tangent                              !    !  !   !  !
                 return ! undefined tangent to intersection curve !    !  !   !  !
              end if ! <------------------------------------------+    !  !   !  !
              !                                                        !  !   !  !
              curvature(1) = max(EPSfp, curvature(1))                  !  !   !  !
              hnext = FRACcurvature_radius / curvature(1)              !  !   !  !
              hnext = max(hnext, hmin)                                 !  !   !  !
              hnext = min(hnext, hmax)                                 !  !   !  !
              if ( lambda*h <= hnext ) then ! <---------+              !  !   !  !
                 h = hnext                              !              !  !   !  !
                 exit inner                             !              !  !   !  !
              end if ! <--------------------------------+              !  !   !  !
           end if ! <--------------------------------------------------+  !   !  !
        end if ! <--------------------------------------------------------+   !  !
        !                                                                     !  !
        ! Newton either failed to converge or produced                        !  !
        ! an unsatisfactory solution => backtrack                             !  !
        h = FRACbacktrack * h                                                 !  !
        if ( h < EPSh ) then ! <---------------+                              !  !
           ! h << h0 (indefinite backtracking) !                              !  !
           stat = 5                            !                              !  !
           return                              !                              !  !
        end if ! <-----------------------------+                              !  !
        !                                                                     !  !
     end do inner ! <---------------------------------------------------------+  !
     !                                                                           !
     ! insert the new point                                                      !
     call insert_polyline_point( &                                               !
          uv, &                                                                  !
          xyz, &                                                                 !
          stat_insertion, &                                                      !
          polyline, &                                                            !
          i=polyline%np )                                                        !
     if ( stat_insertion > 0 ) then ! <--------+                                 !
        stat = stat_insertion                  !                                 !
        return                                 !                                 !
     end if ! <--------------------------------+                                 !
     !                                                                           !
     wprev = w                                                                   !
     !                                                                           !
  end do outer ! <---------------------------------------------------------------+

  ! last point
  polyline%np = max(1, polyline%np)
  call insert_polyline_point( &
       uv_endpoints(1:2,1:2,2), &
       xyz_endpoints(1:3,2), &
       stat_insertion, &
       polyline, &
       i=polyline%np ) 
  if ( stat_insertion > 0 ) then ! <--------+
     stat = stat_insertion                  !
     return                                 !
  end if ! <--------------------------------+

  IF ( DEBUG ) PRINT *,'TRACING OK, NP =',polyline%np

end subroutine trace_intersection_polyline



subroutine loop_detection_criterion( &
     bpn, &
     stat, &
     param_vector, &
     known_to_intersect, &
     n_collineal, &
     randomize )
  USE MOD_UTIL
  use mod_math
  use mod_geometry
  use mod_polynomial
  use mod_separation
  ! Returns:  stat = 0 if Hohmeyer's criterion is satisfied
  !                = 1 if Gauss map #1 must be subdivided
  !                = 2 if Gauss map #2 must be subdivided
  !                = 3 if both Gauss maps must be subdivided
  implicit none
  LOGICAL, PARAMETER :: DEBUG = ( GLOBALDEBUG .AND. .false. )
  CHARACTER :: STRNUM
  type(ptr_polynomial), intent(in)  :: bpn(2)
  integer,              intent(out) :: stat
  real(kind=fp),        intent(out) :: param_vector(3)
  logical,              intent(in)  :: known_to_intersect
  real(kind=fp),        intent(in)  :: n_collineal(3)
  logical, optional,    intent(in)  :: randomize
  type(type_matrix)                 :: sep(2)
  integer                           :: nbcp, nsep(2)
  real(kind=fp)                     :: rot(3,3)
  logical                           :: separable
  real(kind=fp)                     :: vec1(3), vec2(3), vecsqr
  real(kind=fp)                     :: gaussmapsize(2), coneaxe(3)
  integer                           :: isurf

 IF ( DEBUG ) PRINT *,'KNOWN TO INTERSECT?',known_to_intersect

 if ( present(randomize) ) then ! <-----------------------------+
    if ( randomize ) then ! <-------------------------------+   !
       call random_rotation_matrix3d( rot )                 !   !
    end if ! <----------------------------------------------+   !
 end if ! <-----------------------------------------------------+

  do isurf = 1,2
     nbcp = (bpn(isurf)%ptr%degr(1) + 1) * (bpn(isurf)%ptr%degr(2) + 1)
     allocate( sep(isurf)%mat(nbcp,3) )
     call rearrange_for_separability_test( &
          reshape( bpn(isurf)%ptr%coef(&
          1:bpn(isurf)%ptr%degr(1)+1,&
          1:bpn(isurf)%ptr%degr(2)+1,&
          1:3), [nbcp,3] ), &
          nbcp, &
          n_collineal, &
          sep(isurf)%mat, &
          nsep(isurf), &
          .true. )
     if ( known_to_intersect ) then
        nsep(isurf) = nsep(isurf) + 1
        sep(isurf)%mat(nsep(isurf),1:3) = n_collineal
     end if
     if ( present(randomize) ) then ! <------------------------------+
        if ( randomize ) then ! <--------------------------------+   !
           sep(isurf)%mat(1:nsep(isurf),1:3) = &                 !   !
                matmul( sep(isurf)%mat(1:nsep(isurf),1:3), rot ) !   !
        end if ! <-----------------------------------------------+   !
     end if ! <------------------------------------------------------+

     IF ( DEBUG ) THEN
        WRITE (STRNUM,'(I1)') ISURF
        CALL WRITE_MATRIX( SEP(ISURF)%MAT(1:NSEP(ISURF),1:3), NSEP(ISURF), 3, &
             'dev_intersection/sep' // strnum // '.dat' )
     END IF
  end do

  ! first, we search for a vector vec1 such that vec1.n1 < 0 and vec1.n2 > 0
  ! for all ni in the Gauss map of surface #i
  if ( known_to_intersect ) then ! <-------------------------------------------+
     ! if the two gauss maps are known to intersect at a point                 !
     ! (with common normal direction n_collineal), then vec1                   !
     ! only has 1 degree of freedom                                            !
     call intersect_gaussmaps_elsewhere( &                                     !
          sep, &                                                               !
          nsep-1, &                                                              !
          n_collineal, &                                                       !
          separable, &                                                         !
          stat, &                                                              !
          vec1 )                                                               !
  else ! ----------------------------------------------------------------------+
     call separating_plane( &                                                  !
          sep(1)%mat(1:nsep(1),1:3), &                                         !
          sep(2)%mat(1:nsep(2),1:3), &                                         !
          nsep(1), &                                                           !
          nsep(2), &                                                           !
          vec1, &                                                              !
          separable )                                                          !
  end if ! <-------------------------------------------------------------------+
  
  IF ( DEBUG ) THEN
     PRINT *,'P1?',SEPARABLE
     IF ( SEPARABLE ) PRINT *,vec1
  END IF

  if ( separable ) then ! <----------------------------------------------------+
     ! then we search for a vector vec2 such that vec2.n1 > 0 and vec2.n2 > 0  !
     call separating_plane( &                                                  !
          sep(1)%mat(1:nsep(1),1:3), &                                         !
          -sep(2)%mat(1:nsep(2),1:3), &                                        !
          nsep(1), &                                                           !
          nsep(2), &                                                           !
          vec2, &                                                              !
          separable )                                                          !
     IF ( DEBUG ) THEN
        PRINT *,'P2?',SEPARABLE
        IF ( SEPARABLE ) PRINT *,vec2
     END IF
  end if ! <-------------------------------------------------------------------+
     

  if ( separable ) then ! <----------------------------------------------------+
     param_vector = cross( vec1, vec2 )                                        !
     vecsqr = sum( param_vector**2 )                                           !
     if ( vecsqr > EPSfpsqr ) then ! <----------------+                        !
        stat = 0                                      !                        !
        param_vector = param_vector / sqrt( vecsqr )  !                        !
        if ( present(randomize) ) then ! <--------+
           if ( randomize ) then ! <-----------+  !
              param_vector = &                 !  !
                   matmul( rot, param_vector ) !  !
           end if ! <--------------------------+  !
        end if ! <--------------------------------+
        IF ( DEBUG ) THEN
           PRINT *,'P = ',param_vector
           DO ISURF = 1,2
              WRITE (STRNUM,'(I1)') ISURF
              CALL WRITE_MATRIX( SEP(ISURF)%MAT(1:NSEP(ISURF),1:3), NSEP(ISURF), 3, &
                   'dev_intersection/sep' // strnum // '.dat' )
              CALL WRITE_POLYNOMIAL( BPN(ISURF)%PTR, 'dev_intersection/debugld_reg' // strnum // '.bern' )
           END DO
           STOP
        END IF
     else ! ------------------------------------------+                        !
        ! an error has occured                        !                        !
        stat = -1                                     !                        !
        separable = .false.                           !                        !
     end if ! <---------------------------------------+                        !
  end if ! <-------------------------------------------------------------------+


  if ( .not.separable ) then ! <-----------------------------------------------+
     ! which surface has the "largest" Gauss map?                              !
     ! we use an approximate solid angle as a surrogate for Gauss map size     !
     do isurf = 1,2 ! <--------------------------------------------------+     !
        call bounding_cone( &                                            !     !
          sep(isurf)%mat(1:nsep(isurf),1:3), &                           !     !
          coneaxe, &                                                     !     !
          gaussmapsize(isurf) )                                          !     !
     end do ! <----------------------------------------------------------+     !
     !                                                                         !
     if ( all(gaussmapsize > 0.4_fp * CSTpi) ) then ! <------------------+     !
        stat = 3                                                         !     !
     else ! -------------------------------------------------------------+     !
        stat = maxloc(gaussmapsize, dim=1)                               !     !
     end if ! <----------------------------------------------------------+     !
  end if ! <-------------------------------------------------------------------+

  deallocate( sep(1)%mat, sep(2)%mat )

end subroutine loop_detection_criterion



subroutine classify_endpoint( &
     uv, &
     tng, &
     uvbox, &
     stat )
  use mod_math
  use mod_tolerances
  !stat = -1 : exiting
  !        0 : interior
  !        1 : entering
  !        2 : ambiguous
  implicit none
  real(kind=fp), intent(in)  :: uv(2)
  real(kind=fp), intent(in)  :: tng(2) ! tangential direction
  real(kind=fp), intent(in)  :: uvbox(4)
  integer,       intent(out) :: stat
  real(kind=fp)              :: uvloc, tngmax
  integer                    :: statvar(2)
  integer                    :: ivar

  tngmax = max(maxval(abs(tng)), 1._fp)
  do ivar = 1,2
     uvloc = ab2n1p1(uv(ivar), uvbox(2*ivar-1), uvbox(2*ivar))
     if ( is_in_open_interval(uvloc, -1._fp, 1._fp, tolerance=EPSuv) ) then
        statvar(ivar) = 0 ! interior
     else
        if ( abs(tng(ivar)) < 100._fp*EPSuv*tngmax ) then
           statvar(ivar) = 2 ! ambiguous
        elseif ( tng(ivar)*uvloc > 0._fp ) then
           statvar(ivar) = -1 ! exiting
        else
           statvar(ivar) =  1 ! entering
        end if
     end if
  end do
  
  if ( any(statvar == -1) ) then
     stat = -1 ! exiting
  else
     if ( any(statvar == 1) ) then
        stat = 1 ! entering
     else
        !stat = sum(statvar) ! interior or ambiguous
        stat = 0 ! interior
     end if
  end if

end subroutine classify_endpoint



subroutine intersect_intersection_curves( &
     interdata, &
     curvpair, &
     stat )
  USE MOD_UTIL
  use mod_math
  use mod_diffgeom
  use mod_tolerances
  use mod_types_intersection
  implicit none
  LOGICAL, PARAMETER :: DEBUG = ( GLOBALDEBUG .AND. .false. )
  type(type_intersection_data), target, intent(inout) :: interdata
  integer,                              intent(in)    :: curvpair(2)
  type(ptr_intersection_curve)                        :: curv(2)
  logical                                             :: may_intersect
  real(kind=fp)                                       :: uvbox(2,2)
  integer                                             :: numsurf(2)
  type(type_matrix)                                   :: polylineuv(2)
  integer                                             :: np(2)
  integer, allocatable                                :: ipls(:,:)
  real(kind=fp), allocatable                          :: lambda(:,:)
  integer                                             :: npts
  type(ptr_surface)                                   :: surf(3)
  real(kind=fp)                                       :: uv3(2,3), xyz(3), uv2(2,2)
  real(kind=fp)                                       :: w, wi, wip1
  real(kind=fp), dimension(6)                         :: lowerb, upperb
  real(kind=fp), allocatable                          :: uv3tmp(:,:,:)
  integer                                             :: stat
  integer                                             :: idpt
  integer                                             :: isurf, jsurf, ivar, icurv, ipt, jpt
  INTEGER :: FID = 13, I
  CHARACTER :: STRNUM

  do icurv = 1,2
     curv(icurv)%ptr => interdata%curves(curvpair(icurv))
  end do

  ! first, check if the two curves can intersect, i.e. they share a common incident surface,
  ! and their uv-boxes on this surface overlap
  may_intersect = .false.
  jloop : do jsurf = 1,2 !  <------------------------------------------------------------------------+
     iloop : do isurf = 1,2 ! <-------------------------------------------------------------------+  !
        if ( associated(curv(1)%ptr%surf(isurf)%ptr, curv(2)%ptr%surf(jsurf)%ptr) ) then ! <---+  !  !
           ! the two intersection curve share a common incident surface                        !  !  !
           ! compute the intersection of their bounding boxes in the uv-space of that surface  !  !  !
           do ivar = 1,2 ! <----------------------------------------------------------------+  !  !  !
              uvbox(:,ivar) = [ &                                                           !  !  !  !
                   max(curv(1)%ptr%uvbox(1,ivar,isurf), curv(2)%ptr%uvbox(1,ivar,jsurf)), & !  !  !  !
                   min(curv(1)%ptr%uvbox(2,ivar,isurf), curv(2)%ptr%uvbox(2,ivar,jsurf)) ]  !  !  !  !
              if ( uvbox(2,ivar) - uvbox(1,ivar) < EPSuv ) then ! <------------+            !  !  !  !
                 ! the uv-boxes do not overlap                                 !            !  !  !  !
                 cycle iloop                                                   !            !  !  !  !
              end if ! <-------------------------------------------------------+            !  !  !  !
           end do ! <-----------------------------------------------------------------------+  !  !  !
           ! if we got this far, the two curves may intersect                                  !  !  !
           numsurf = [isurf, jsurf]                                                            !  !  !
           may_intersect = .true.                                                              !  !  !
           exit jloop                                                                          !  !  !
        end if ! <-----------------------------------------------------------------------------+  !  !
     end do iloop ! <-----------------------------------------------------------------------------+  !
  end do jloop ! <-----------------------------------------------------------------------------------+

  if (.not.may_intersect) return ! the curves cannot intersect

  ! find intersection points between the polylines in the uv-space of the common incident surface
  do icurv = 1,2
     np(icurv) = curv(icurv)%ptr%polyline%np
     allocate(polylineuv(icurv)%mat(2,np(icurv)))
     polylineuv(icurv)%mat(1:2,1:np(icurv)) = curv(icurv)%ptr%polyline%uv(1:2,numsurf(icurv),1:np(icurv))
  end do
  npts = 0
  call intersect_2Dpolylines( &
       polylineuv, &
       [1,1], &
       np, &
       npts, &
       ipls, &
       lambda )
  IF ( .false. ) THEN
     DO ICURV = 1,2
        WRITE (STRNUM,'(I1)') ICURV
        OPEN(UNIT=FID, FILE='polylines/xy'//strnum//'.dat', ACTION='WRITE')
        WRITE (FID,*) NP(ICURV)
        DO IPT = 1,NP(ICURV)
           WRITE (FID,*) POLYLINEUV(ICURV)%MAT(1:2,IPT)
        END DO
        CLOSE(FID)
     END DO
     OPEN(UNIT=FID, FILE='polylines/result.dat', ACTION='WRITE')
     WRITE (FID,*) NPTS
     DO IPT = 1,NPTS
        WRITE (FID,*) IPLS(:,IPT)
     END DO
     DO IPT = 1,NPTS
        WRITE (FID,*) LAMBDA(:,IPT)
     END DO
     CLOSE(FID)
  END IF

  if ( npts < 1 ) then
     ! the polylines do not intersect, we assume the curves do not either
     deallocate(polylineuv(1)%mat, polylineuv(2)%mat)
     return
  end if

  !! refine all intersection points using Newton-Raphson algorithm and add them to the collections  
  ! form the triplet of surfaces (one common incident surface, plus one other for each curve)
  surf(1)%ptr => curv(1)%ptr%surf(numsurf(1))%ptr
  do icurv = 1,2
     surf(1+icurv)%ptr => curv(icurv)%ptr%surf(1+mod(numsurf(icurv),2))%ptr
  end do
  ! lower and upper bounds of feasible domain
  lowerb(1:2) = uvbox(1,1:2)
  upperb(1:2) = uvbox(2,1:2)
  do icurv = 1,2
     isurf = 1 + mod(numsurf(icurv),2)
     lowerb(2*(icurv+1)-1:2*(icurv+1)) = curv(icurv)%ptr%uvbox(1,1:2,isurf)
     upperb(2*(icurv+1)-1:2*(icurv+1)) = curv(icurv)%ptr%uvbox(2,1:2,isurf)
  end do
  lowerb = lowerb - EPSuv
  upperb = upperb + EPSuv

  allocate(uv3tmp(2,3,npts))
  do ipt = 1,npts
     uv3tmp(:,1,ipt) = &
          (1._fp - lambda(1,ipt)) * curv(1)%ptr%polyline%uv(:,numsurf(1),ipls(1,ipt)) + &
          lambda(1,ipt) * curv(1)%ptr%polyline%uv(:,numsurf(1),ipls(1,ipt)+1)
     do icurv = 1,2 ! <-----------------------------------------------------------------------------+
        isurf = 1 + mod(numsurf(icurv),2)                                                           !
        uv3tmp(:,1+icurv,ipt) = &                                                                   !
             (1._fp - lambda(icurv,ipt)) * curv(icurv)%ptr%polyline%uv(:,isurf,ipls(icurv,ipt)) + & !
             lambda(icurv,ipt) * curv(icurv)%ptr%polyline%uv(:,isurf,ipls(icurv,ipt)+1)             !
     end do ! <-------------------------------------------------------------------------------------+
  end do

  do ipt = 1,npts ! <----------------------------------------------------------------------------------+
     ! set initial iterate (result of polyline intersection)                                           !
     uv3 = uv3tmp(:,:,ipt)                                                                             !
     !uv3(:,1) = &                                                                                      !
     !     (1._fp - lambda(1,ipt)) * curv(1)%ptr%polyline%uv(:,numsurf(1),ipls(1,ipt)) + &              !
     !     lambda(1,ipt) * curv(1)%ptr%polyline%uv(:,numsurf(1),ipls(1,ipt)+1)                          !
     !do icurv = 1,2 ! <-----------------------------------------------------------------------------+  !
     !   isurf = 1 + mod(numsurf(icurv),2)                                                           !  !
     !   uv3(:,1+icurv) = &                                                                          !  !
     !        (1._fp - lambda(icurv,ipt)) * curv(icurv)%ptr%polyline%uv(:,isurf,ipls(icurv,ipt)) + & !  !
     !        lambda(icurv,ipt) * curv(icurv)%ptr%polyline%uv(:,isurf,ipls(icurv,ipt)+1)             !  !
     !end do ! <-------------------------------------------------------------------------------------+  !
     !                                                                                                 !
     ! run Newton-Raphson algorithm                                                                    !
     IF ( DEBUG ) THEN
        PRINT *,'SMOOTH? ', curv(1)%ptr%smooth, curv(2)%ptr%smooth
        PRINT *,'UV0 = ',UV3
     END IF
     if ( curv(2)%ptr%smooth ) then
        call newton_three_surfaces_1tangential( &
             surf, &
             lowerb, &
             upperb, &
             stat, &
             uv3, &
             xyz )
     else
        call newton_three_surfaces( &                                                                     !
             surf, &                                                                                      !
             lowerb, &                                                                                    !
             upperb, &                                                                                    !
             stat, &                                                                                      !
             uv3, &                                                                                       !
             xyz )                                                                                        !
     end if
     !                                                                                                 !
     IF ( STAT /= 0 ) THEN
        PRINT *,'******************************************'
        PRINT *,' NEWTON 3 SURFACES DID NOT CONVERGE'
        PRINT *,'******************************************'
        PAUSE
        CYCLE
     END IF
     IF ( DEBUG ) THEN
        IF ( STAT == 0 ) THEN
           PRINT *,'NEWTON 3 SURFACES CONVERGED'
           PRINT *,' UV =',UV3
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
     END IF
     if ( stat > 0 ) then ! <---------------+                                                          !
        ! Newton failed to converge         !                                                          !
        STOP 'NEWTON_3_SURFACES FAILED TO CONVERGE'
        return                              !                                                          !
     end if ! <-----------------------------+                                                          !
     !                                                                                                 !
     ! add the split point the intersection points collection                                          !
     call add_intersection_point( &                                                                    !
          uv3, &                                                                                       !
          xyz, &                                                                                       !
          surf, &                                                                                      !
          3, &                                                                                         !
          interdata, &                                                                                 !
          idpt )                                                                                       !
     !                                                                                                 !
     ! add the split point to each intersected curve's polyline                                        !
     curv_loop : do icurv = 1,2 ! <---------------------------------------------------+                !
        ! compare the new split point to already discovered split points              !                !
        do jpt = 1,curv(icurv)%ptr%nsplit ! <---------------------------------+       !                !
           if ( idpt == curv(icurv)%ptr%isplit(1,jpt) ) cycle curv_loop       !       !                !
        end do ! <------------------------------------------------------------+       !                !
        !                                                                             !                !
        ! locate the segement in which the new split point must be inserted           !                !
        if ( curv(icurv)%ptr%smooth ) then
           ! ...
        else
           w    = dot_product(curv(icurv)%ptr%param_vector, xyz)                         !                !
           wi   = dot_product(curv(icurv)%ptr%param_vector, &                            !                !
                curv(icurv)%ptr%polyline%xyz(:,ipls(icurv,ipt)))                   !                !
           wip1 = dot_product(curv(icurv)%ptr%param_vector, &                            !                !
                curv(icurv)%ptr%polyline%xyz(:,ipls(icurv,ipt)+1))                 !                !
           while_loop : do ! <---------------------------------------------------+       !                !
              if ( ipls(icurv,ipt) < 1 .or. &                                    !       !                !
                   ipls(icurv,ipt) >= curv(icurv)%ptr%polyline%np ) then ! <--+  !       !                !
                 PRINT *,'W =',W
                 PRINT *,&
                      dot_product(curv(icurv)%ptr%param_vector, curv(icurv)%ptr%polyline%xyz(:,1)), &
                      dot_product(curv(icurv)%ptr%param_vector, curv(icurv)%ptr%polyline%xyz(:,curv(icurv)%ptr%polyline%np))
                 STOP 'COULD NOT FIND CORRECT SEGMENT'                        !  !       !                !
              end if ! <------------------------------------------------------+  !       !                !
              if ( w > wip1 ) then ! <-----------------------------------+       !       !                !
                 ! the split point is located downstream                 !       !       !                !
                 ipls(icurv,ipt) = ipls(icurv,ipt) + 1                   !       !       !                !
                 wi = wip1                                               !       !       !                !
                 wip1 = dot_product(curv(icurv)%ptr%param_vector, &      !       !       !                !
                      curv(icurv)%ptr%polyline%xyz(:,ipls(icurv,ipt)+1)) !       !       !                !
              elseif ( w < wi ) then ! ----------------------------------+       !       !                !
                 ! the split point is located upstream                   !       !       !                !
                 ipls(icurv,ipt) = ipls(icurv,ipt) - 1                   !       !       !                !
                 IF ( ipls(icurv,ipt) < 1 ) THEN
                    PRINT *,idpt,curv(icurv)%ptr%isplit(1,1:curv(icurv)%ptr%nsplit)
                    PRINT *,XYZ
                    DO JPT = 1,curv(icurv)%ptr%nsplit
                       PRINT *,INTERDATA%POINTS(curv(icurv)%ptr%isplit(1,JPT))%XYZ
                    END DO
                 END IF
                 wi = dot_product(curv(icurv)%ptr%param_vector, &        !       !       !                !
                      curv(icurv)%ptr%polyline%xyz(:,ipls(icurv,ipt)))   !       !       !                !
                 wip1 = wi                                               !       !       !                !
              else ! ----------------------------------------------------+       !       !                !
                 ! we found the correct segment                          !       !       !                !
                 exit while_loop                                         !       !       !                !
              end if ! <-------------------------------------------------+       !       !                !
           end do while_loop ! <-------------------------------------------------+       !                !
        end if
        !                                                                             !                !
        ! insert a new polyline point                                                 !                !
        uv2(:,numsurf(icurv)) = uv3(:,1)                                              !                !
        uv2(:,1+mod(numsurf(icurv),2)) = uv3(:,1+icurv)                               !                !
        call insert_polyline_point( &                                                 !                !
             uv2, &                                                                   !                !
             xyz, &                                                                   !                !
             stat, &                                                                  !                !
             curv(icurv)%ptr%polyline, &                                              !                !
             i=ipls(icurv,ipt) )                                                      !                !
        !                                                                             !                !
        if ( stat > 0 ) then ! <----------------------+                               !                !
           ! the polyine point could no be inserted   !                               !                !
           stat = stat + 2                            !                               !                !
           return                                     !                               !                !
        end if ! <------------------------------------+                               !                !
        !                                                                             !                !
        ! locate in the curve's split points collection where the new point           !                !
        ! must be inserted (split points are sorted in ascending order of w)          !                !
        ! at the same time, update indices of split points located downstream         !                !
        do jpt = curv(icurv)%ptr%nsplit,1,-1 ! <-----------------------------------+  !                !
           if ( curv(icurv)%ptr%isplit(2,jpt) > ipls(icurv,ipt) ) then ! <------+  !  !                !
              curv(icurv)%ptr%isplit(2,jpt) = curv(icurv)%ptr%isplit(2,jpt) + 1 !  !  !                !
           else ! --------------------------------------------------------------+  !  !                !
              exit                                                              !  !  !                !
           end if ! <-----------------------------------------------------------+  !  !                !
        end do ! <-----------------------------------------------------------------+  !                !
        !                                                                             !                !
        ! finally, add the split point the curve's split points list                  !                !
        call insert_column_after( &                                                   !                !
             curv(icurv)%ptr%isplit, &                                                !                !
             2, &                                                                     !                !
             curv(icurv)%ptr%nsplit, &                                                !                !
             [idpt, ipls(icurv,ipt)+1], &                                             !                !
             jpt )                                                                    !                !
     end do curv_loop ! <-------------------------------------------------------------+                !
     !                                                                                                 !
  end do ! <-------------------------------------------------------------------------------------------+
  
  deallocate(polylineuv(1)%mat, polylineuv(2)%mat)
  
end subroutine intersect_intersection_curves



recursive subroutine intersect_curve_surface( &
     root_c, &
     root_s, &
     region_c, &
     region_s, &
     tuvxyz, &
     ntuvxyz, &
     stat_degeneracy )
  use mod_util
  use mod_math
  use mod_bernstein
  use mod_polynomial
  use mod_diffgeom
  use mod_obb
  use mod_regiontree
  use mod_tolerances
  implicit none
  LOGICAL, PARAMETER :: DEBUG = ( GLOBALDEBUG .AND. .TRUE. )
  type(type_curve),           intent(in)    :: root_c
  type(type_surface),         intent(in)    :: root_s
  type(type_region), target,  intent(inout) :: region_c
  type(type_region), target,  intent(inout) :: region_s
  real(kind=fp), allocatable, intent(inout) :: tuvxyz(:,:)
  integer,                    intent(inout) :: ntuvxyz
  integer,                    intent(inout) :: stat_degeneracy
  integer, allocatable                      :: ipts(:)
  integer                                   :: npts
  logical                                   :: overlap
  logical                                   :: interior(2)
  type(ptr_polynomial)                      :: poly(2)
  logical                                   :: separable
  real(kind=fp), dimension(3)               :: xyz_c, xyz_s
  real(kind=fp), dimension(3)               :: tuv, xyz, tuv_subdiv
  real(kind=fp)                             :: toltuv
  integer                                   :: stat_newpoint, stat_contact
  integer                                   :: stat_subdiv
  integer                                   :: nchild(2)
  type(type_region), pointer                :: newregion_c
  type(type_region), pointer                :: newregion_s
  integer                                   :: i, j, k, ipt, jpt, ichild, jchild

  ! return if a degeneracy has been encountered previously
  if ( stat_degeneracy /= 0 ) return

  IF ( DEBUG ) THEN
     !PAUSE
     PRINT *,''; PRINT *,'';
     PRINT *,'INTERSECT_CURVE_SURFACE'
     PRINT *,' TBOX =',REGION_C%UVBOX
     PRINT *,'UVBOX =',REGION_S%UVBOX
  END IF

  ! inherit from parents all the  points contained in the current regions
  if ( ntuvxyz > 0 ) then ! <-------------------------------------------------+
     call inherit_points( &                                                   !
          region_c, &                                                         !
          tuvxyz(1,1:ntuvxyz), &                                              !
          tuvxyz(7,1:ntuvxyz), &                                              !
          ntuvxyz )                                                           !
     call inherit_points( &                                                   !
          region_s, &                                                         !
          tuvxyz(2:3,1:ntuvxyz), &                                            !
          tuvxyz(7,1:ntuvxyz), &                                              !
          ntuvxyz )                                                           !
  end if ! <------------------------------------------------------------------+

  IF ( DEBUG ) THEN
     PRINT *,'T, TOL ='
     IF ( REGION_C%NPTS > 0) THEN
        CALL PRINT_MAT(TRANSPOSE( &
             TUVXYZ([1,7],REGION_C%IPTS(1:REGION_C%NPTS)) ))
     END IF
     PRINT *,'UV, TOL ='
     IF ( REGION_S%NPTS > 0) THEN
        CALL PRINT_MAT(TRANSPOSE( &
             TUVXYZ([2,3,7],REGION_S%IPTS(1:REGION_S%NPTS)) ))
     END IF
  END IF
 

  ! get list of already discovered points contained in both the current regions  
  npts = 0
  if ( region_c%npts > 0 .and. region_s%npts > 0 ) then ! <-------------------+
     call intersection_arrays( &                                              !
          region_c%ipts(1:region_c%npts), &                                   !
          region_s%ipts(1:region_s%npts), &                                   !
          ipts )                                                              !
     if ( allocated(ipts) ) npts = size(ipts)                                 !
  end if ! <------------------------------------------------------------------+


  ! if there are no already discovered points contained in both the current regions, check if 
  ! they intersect at endpoints/corners (2*4 tests)
  if ( npts == 0 ) then ! <-------------------------------------------------------------------+
     do k = 1,2 ! <----------------------------------------------------------------+          !
        do j = 1,2 ! <----------------------------------------------------------+  !          !
           xyz_s = region_s%poly(1)%ptr%coef( &                                 !  !          !
                1 + (j-1)*region_s%poly(1)%ptr%degr(1), &                       !  !          !
                1 + (k-1)*region_s%poly(1)%ptr%degr(2), &                       !  !          !
                1:3 )                                                           !  !          !
           do i = 1,2 ! <----------------------------------------------------+  !  !          !
              xyz_c = region_c%poly(1)%ptr%coef( &                           !  !  !          !
                   1 + (i-1)*region_c%poly(1)%ptr%degr(1), &                 !  !  !          !
                   1:3, &                                                    !  !  !          !
                   1 )                                                       !  !  !          !
              if ( sum( (xyz_s - xyz_c)**2 ) < EPSxyzsqr ) then ! <----+     !  !  !          !
                 xyz = 0.5_fp * ( xyz_c + xyz_s )                      !     !  !  !          !
                 tuv = real([i,j,k]-1, kind=fp)                        !     !  !  !          !
                 tuv(1) = (1._fp - tuv(1))*region_c%uvbox(1) + &       !     !  !  !          !
                      tuv(1)*region_c%uvbox(2)                         !     !  !  !          !
                 tuv(2) = (1._fp - tuv(2))*region_s%uvbox(1) + &       !     !  !  !          !
                      tuv(2)*region_s%uvbox(2)                         !     !  !  !          !
                 tuv(3) = (1._fp - tuv(3))*region_s%uvbox(3) + &       !     !  !  !          !
                      tuv(3)*region_s%uvbox(4)                         !     !  !  !          !
                 IF ( DEBUG ) PRINT *,'TUV =',TUV
                 IF ( DEBUG ) PRINT *,'XYZ =',XYZ
                 !                                                     !     !  !  !          !
                 call check_curve_surface_intersection_point( &        !     !  !  !          !
                      root_c, &                                        !     !  !  !          !
                      root_s, &                                        !     !  !  !          !
                      tuv(1), &                                        !     !  !  !          !
                      tuv(2:3), &                                      !     !  !  !          !
                      stat_contact )                                   !     !  !  !          !
                 IF ( DEBUG ) PRINT *,'CURVE-SURFACE STATPOINT =',stat_contact
                 if ( stat_contact == 2 ) then ! <----------+          !     !  !  !          !
                    ! high-order tangential contact point   !          !     !  !  !          !
                    ! => the curve is presumably a subset   !          !     !  !  !          !
                    ! of the surface                        !          !     !  !  !          !
                    stat_degeneracy = 50                    !          !     !  !  !          !
                    return                                  !          !     !  !  !          !
                 end if ! <---------------------------------+          !     !  !  !          !
                 !                                                     !     !  !  !          !
                 call append_vector( &                                 !     !  !  !          !
                      [tuv,xyz,EPSuv], &                               !     !  !  !          !
                      7, &                                             !     !  !  !          !
                      tuvxyz, &                                        !     !  !  !          !
                      ntuvxyz )                                        !     !  !  !          !
                 !                                                     !     !  !  !          !
                 call append_n( &                                      !     !  !  !          !
                      ipts, &                                          !     !  !  !          !
                      npts, &                                          !     !  !  !          !
                      [ntuvxyz], &                                     !     !  !  !          !
                      1, &                                             !     !  !  !          !
                      unique=.true. )                                  !     !  !  !          !
                 !                                                     !     !  !  !          !
                 call add_point_bottom_up( region_c, ntuvxyz )         !     !  !  !          !
                 call add_point_bottom_up( region_s, ntuvxyz )         !     !  !  !          !
                 !                                                     !     !  !  !          !
              end if ! <-----------------------------------------------+     !  !  !          !
           end do ! <--------------------------------------------------------+  !  !          !
        end do ! <--------------------------------------------------------------+  !          !
     end do ! <--------------------------------------------------------------------+          !
  end if ! <----------------------------------------------------------------------------------+

  if ( npts == 0 ) then ! <-------------------------------------------------------------------+
     ! if there are still no intersection points, check if the bounding boxes overlap         !
     if ( .not.associated(region_c%xyzbox) ) then ! <--------+                                !
        allocate( region_c%xyzbox )                          !                                !
        call bernOBB1( &                                     !                                !
             region_c%poly(1)%ptr%coef( &                    !                                !
             1:region_c%poly(1)%ptr%degr(1)+1, &             !                                !
             1:3, &                                          !                                !
             1), &                                           !                                !
             region_c%poly(1)%ptr%degr(1), &                 !                                !
             region_c%xyzbox )                               !                                !
     end if ! <----------------------------------------------+                                !
     !                                                                                        !
     if ( .not.associated(region_s%xyzbox) ) then ! <--------+                                !
        allocate( region_s%xyzbox )                          !                                !
        call bernOBB2( &                                     !                                !
             region_s%poly(1)%ptr%coef( &                    !                                !
             1:region_s%poly(1)%ptr%degr(1)+1, &             !                                !
             1:region_s%poly(1)%ptr%degr(2)+1, &             !                                !
             1:3), &                                         !                                !
             region_s%poly(1)%ptr%degr, &                    !                                !
             region_s%xyzbox )                               !                                !
     end if ! <----------------------------------------------+                                !
     !                                                                                        !
     call overlap_OBBs( &                                                                     !
          region_c%xyzbox, &                                                                  !
          region_s%xyzbox, &                                                                  !
          overlap )                                                                           !
     !                                                                                        !
     if ( .not.overlap ) return ! disjoint bounding boxes => empty intersection               !
     !                                                                                        !
  elseif ( npts > 0 ) then ! <----------------------------------------------------------------+
     ! are some of the already discovered points interior to the curve or surface region?...  !
     do jpt = 1,npts  ! <------------------------------------------+                          !
        ipt = ipts(jpt)                                            !                          !
        interior(1) = is_in_open_interval( &                       !                          !
             tuvxyz(1,ipt), &                                      !                          !
             region_c%uvbox(1), region_c%uvbox(2), &               !                          !
             tolerance=EPSregion )                                 !                          !
        interior(2) = ( &                                          !                          !
             is_in_open_interval( &                                !                          !
             tuvxyz(2,ipt), &                                      !                          !
             region_s%uvbox(1), region_s%uvbox(2), &               !                          !
             tolerance=EPSregion ) .or. &                          !                          !
             is_in_open_interval( &                                !                          !
             tuvxyz(3,ipt), &                                      !                          !
             region_s%uvbox(3), region_s%uvbox(4), &               !                          !
             tolerance=EPSregion ) )                               !                          !
        !                                                          !                          !
        ! ... if so, subdivide both regions at that point          !                          !
        if ( any(interior) ) then ! <-------------------+          !                          !
           tuv_subdiv = tuvxyz(1:3,ipt)                 !          !                          !
           exit                                         !          !                          !
        end if ! <--------------------------------------+          !                          !
     end do ! <----------------------------------------------------+                          !
     !                                                                                        !
     if ( all(.not.interior) ) then  ! <--------------------------------------------------+   !
        ! if only one intersection point, check if the pair can intersect at other points !   !
        if ( npts == 1 ) then ! <-----------------------------------------------+         !   !
           IF ( DEBUG ) THEN
              PRINT *,'TUV =',TUVXYZ(1:3,IPT)
              PRINT *,'ALL NOT INTERIOR'
           END IF
           poly(1)%ptr => region_c%poly(1)%ptr                                  !         !   !
           poly(2)%ptr => region_s%poly(1)%ptr                                  !         !   !
           call intersect_elsewhere( &                                          !         !   !
                poly, &                                                         !         !   !
                tuvxyz(4:6,ipts(1)), &                                          !         !   !
                separable, &                                                    !         !   !
                randomize=.true. )                                              !         !   !
           nullify(poly(1)%ptr, poly(2)%ptr)                                    !         !   !
           IF ( DEBUG ) PRINT *,'CURVE-SURFACE SEPARABLE?',SEPARABLE
           !                                                                    !         !   !
           if ( separable ) return ! the pair cannot intersect at other points  !         !   !
           !                                                                    !         !   !
           ! special case: tangential contact                                   !         !   !
           call check_curve_surface_intersection_point( &                       !         !   !
                root_c, &                                                       !         !   !
                root_s, &                                                       !         !   !
                tuvxyz(1,ipts(1)), &                                            !         !   !
                tuvxyz(2:3,ipts(1)), &                                          !         !   !
                stat_contact )                                                  !         !   !
           if ( stat_contact == 1 ) then ! <----------------------------+       !         !   !
              call curve_surface_separation_at_tangential_contact( &    !       !         !   !
                   root_s, &                                            !       !         !   !
                   tuvxyz(2:3,ipts(1)), &                               !       !         !   !
                   tuvxyz(4:6,ipts(1)), &                               !       !         !   !
                   region_c%poly(1)%ptr, &                              !       !         !   !
                   region_s%poly(1)%ptr, &                              !       !         !   !
                   separable )                                          !       !         !   !
              IF ( DEBUG ) PRINT *,'CS SEPARABLE AT TANGENTIAL CONTACT POINT?',separable
              if ( separable ) return                                   !       !         !   !
           end if ! <---------------------------------------------------+       !         !   !
        end if ! <--------------------------------------------------------------+         !   !
     end if ! <---------------------------------------------------------------------------+   !
     !                                                                                        !
  end if ! <----------------------------------------------------------------------------------+


  ! search for a new intersection point
  if ( npts < 1 .or. all(.not.interior) ) then ! <--------------------------------------------+
     ! use box-constrained Newton-Raphson algorithm                                           !
     ! set initial iterate to the regions' parametric center point                            !
     tuv(1)   = 0.5_fp * (region_c%uvbox(1)     + region_c%uvbox(2)    )                      !
     tuv(2:3) = 0.5_fp * (region_s%uvbox([1,3]) + region_s%uvbox([2,4]))                      !
     call newton_curve_surface( &                                                             !
          root_c, &                                                                           !
          root_s, &                                                                           !
          [region_c%uvbox(1), region_s%uvbox([1,3])] - EPSuv, &                               !
          [region_c%uvbox(2), region_s%uvbox([2,4])] + EPSuv, &                               !
          ![region_c%uvbox(1), region_s%uvbox([1,3])], &                                       !
          ![region_c%uvbox(2), region_s%uvbox([2,4])], &                                       !
          stat_newpoint, &                                                                    !
          tuv, &                                                                              !
          toltuv, &                                                                           !
          xyz )                                                                               !
     !                                                                                        !
     ! if a degeneracy has been encountered, report it                                        !
     if ( stat_newpoint > 1 ) then ! <---------------+                                        !
        stat_degeneracy = stat_newpoint              !                                        !
        CALL WRITE_POLYNOMIAL( REGION_C%POLY(1)%PTR, 'dev_intersection/debugncs_regc.bern' )
        CALL WRITE_POLYNOMIAL( REGION_S%POLY(1)%PTR, 'dev_intersection/debugncs_regs.bern' )
        return                                       !                                        !
     end if ! <--------------------------------------+                                        !
     !                                                                                        !
     if ( stat_newpoint <= 0 ) then ! <---------------------------------------------------+   !
        IF ( DEBUG ) PRINT *,'NEWTON -> TUVXYZ=',TUV,XYZ
        call check_curve_surface_intersection_point( &                                    !   !
             root_c, &                                                                    !   !
             root_s, &                                                                    !   !
             tuv(1), &                                                                    !   !
             tuv(2:3), &                                                                  !   !
             stat_contact )                                                               !   !
        IF ( DEBUG ) PRINT *,'NEWTON, CURVE-SURFACE STATPOINT =',stat_contact
        if ( stat_contact == 2 ) then ! <-----------------------+                         !   !
           ! high-order tangential contact point                !                         !   !
           ! => the curve is presumably a subset of the surface !                         !   !
           stat_degeneracy = 50                                 !                         !   !
           return                                               !                         !   !
        end if ! <----------------------------------------------+                         !   !
        !                                                                                 !   !
        ! if a point has been found, check whether it is a duplicate                      !   !
        if ( ntuvxyz > 0 ) then ! <-----------------------+                               !   !
           call check_unicity_with_tolerance( &           !                               !   !
                tuv, &                                    !                               !   !
                xyz, &                                    !                               !   !
                toltuv, &                                 !                               !   !
                3, &                                      !                               !   !
                tuvxyz(1:3,1:ntuvxyz), &                  !                               !   !
                tuvxyz(4:6,1:ntuvxyz), &                  !                               !   !
                tuvxyz(7,1:ntuvxyz), &                    !                               !   !
                ntuvxyz, &                                !                               !   !
                ipt )                                     !                               !   !
        else ! -------------------------------------------+                               !   !
           ipt = 1                                        !                               !   !
        end if ! <----------------------------------------+                               !   !
        !                                                                                 !   !
        if ( ipt > ntuvxyz ) then ! <---------------------+                               !   !
           ! if this is actually a new intersection point !                               !   !
           call append_vector( &                          !                               !   !
                [tuv,xyz,toltuv], &                       !                               !   !
                7, &                                      !                               !   !
                tuvxyz, &                                 !                               !   !
                ntuvxyz )                                 !                               !   !
           IF ( DEBUG ) THEN
              PRINT *,'TUVXYZ ='
              CALL PRINT_MAT( TRANSPOSE(TUVXYZ(:,1:NTUVXYZ)) )
              PRINT *,'TOLTUV =',TUVXYZ(7,1:NTUVXYZ)
           END IF
        else ! -------------------------------------------+                               !   !
           stat_newpoint = 1                              !                               !   !
        end if ! <----------------------------------------+                               !   !
        ! append that point to the current regions' lists and to that of their ascendants !   !
        if ( is_in_closed_interval( &
             tuv(1), &
             region_c%uvbox(1), &
             region_c%uvbox(2), &
             tolerance=EPSuv) ) then
           !PRINT *,'DIFF T:',tuv(1)-region_c%uvbox(1), region_c%uvbox(2)-tuv(1)
           call add_point_bottom_up( region_c, ipt )                                      !   !
        else
           PRINT *,'*** WARNING : POINT NOT INSIDE TBOX'
        end if
        if ( is_in_closed_interval( &
             tuv(2), &
             region_s%uvbox(1), &
             region_s%uvbox(2), &
             tolerance=EPSuv) .or. &
             is_in_closed_interval( &
             tuv(3), &
             region_s%uvbox(3), &
             region_s%uvbox(4), &
             tolerance=EPSuv) ) then
           !PRINT *,'DIFF UV:',tuv(2)-region_s%uvbox(1), region_s%uvbox(2)-tuv(2), tuv(3)-region_s%uvbox(3), region_s%uvbox(4)-tuv(3)
           call add_point_bottom_up( region_s, ipt )                                      !   !
        else
           PRINT *,'*** WARNING : POINT NOT INSIDE UVBOX'
        end if
        !                                                                                 !   !
     end if  ! <--------------------------------------------------------------------------+   !
     !                                                                                        !
     if ( stat_newpoint <= 0 ) then ! <---------------------------------------------------+   !
        ! subdivide both regions at that point                                            !   !
        tuv_subdiv = tuv                                                                  !   !
     else ! ------------------------------------------------------------------------------+   !
        ! if no new point has been discovered, subdivide both regions at the centerpoint  !   !
        tuv_subdiv(1)   = 0.5_fp * (region_c%uvbox(1)     + region_c%uvbox(2)    )        !   !
        tuv_subdiv(2:3) = 0.5_fp * (region_s%uvbox([1,3]) + region_s%uvbox([2,4]))        !   !
     end if  ! <--------------------------------------------------------------------------+   !
     !                                                                                        !
  end if ! <----------------------------------------------------------------------------------+

  IF ( DEBUG ) PRINT *,'TUV_SUBDIV =',tuv_subdiv

  !! subdivide the curve region
  call subdiv_region( &
       region_c, &
       tuv_subdiv(1), &
       stat_subdiv )
  ! linear change of variable --> local frame of Bézier subcurve
  tuv_subdiv(1) = 0.5_fp * ( ab2n1p1(tuv_subdiv(1), region_c%uvbox(1), region_c%uvbox(2)) + 1._fp )

  if ( stat_subdiv == 1 ) then ! <--------------------------------------------------+
     ! the subdivision point is at one of the curve's endpoints                     !
     nchild(1) = 1                                                                  !
  else ! ---------------------------------------------------------------------------+
     ! the subdivision point is interior to the curve                               !
     nchild(1) = size(region_c%child)                                               !
     if ( stat_subdiv == 0 ) then ! <-------------------------+                     !
        ! the curve region has no children yet                !                     !
        do ichild = 1,size(region_c%child) ! <-------------+  !                     !
           allocate(region_c%child(ichild)%poly(1))        !  !                     !
           allocate(region_c%child(ichild)%poly(1)%ptr)    !  !                     !
        end do ! <-----------------------------------------+  !                     !
        call subdiv_bezier1( &                                !                     !
             region_c%poly(1)%ptr, &                          !                     !
             tuv_subdiv(1), &                                 !                     !
             bl=region_c%child(1)%poly(1)%ptr, &              !                     !
             br=region_c%child(2)%poly(1)%ptr )               !                     !
     end if ! <-----------------------------------------------+                     !
  end if ! <------------------------------------------------------------------------+


  !! subdivide the surface region
  call subdiv_region( &
       region_s, &
       tuv_subdiv(2:3), &
       stat_subdiv )
  ! linear change of variables --> local frame of Bézier subsurface
  tuv_subdiv(2) = 0.5_fp * ( ab2n1p1(tuv_subdiv(2), region_s%uvbox(1), region_s%uvbox(2)) + 1._fp )
  tuv_subdiv(3) = 0.5_fp * ( ab2n1p1(tuv_subdiv(3), region_s%uvbox(3), region_s%uvbox(4)) + 1._fp )

  if ( stat_subdiv == 2 ) then ! <--------------------------------------------------+
     ! the subdivision point is at one of the surface's corners                     !
     nchild(2) = 1                                                                  !
  else ! ---------------------------------------------------------------------------+
     ! the subdivision point is interior to the surface                             !
     nchild(2) = size(region_s%child)                                               !
     if ( stat_subdiv >= 0 ) then ! <-------------------------+                     !
        ! the surface region has no children yet              !                     !
        do ichild = 1,size(region_s%child) ! <-------------+  !                     !
           allocate(region_s%child(ichild)%poly(1))        !  !                     !
           allocate(region_s%child(ichild)%poly(1)%ptr)    !  !                     !
        end do ! <-----------------------------------------+  !                     !
     end if ! <-----------------------------------------------+                     !
  end if ! <------------------------------------------------------------------------+

  if ( stat_subdiv == 0 ) then ! <--------------------------------------------------+
     ! subdivide into 4 children                                                    !
     call subdiv_bezier2( &                                                         !
          region_s%poly(1)%ptr, &                                                   !
          tuv_subdiv(2:3), &                                                        !
          bsw=region_s%child(1)%poly(1)%ptr, &                                      !
          bse=region_s%child(2)%poly(1)%ptr, &                                      !
          bnw=region_s%child(3)%poly(1)%ptr, &                                      !
          bne=region_s%child(4)%poly(1)%ptr )                                       !
  elseif ( stat_subdiv == 1 ) then ! -----------------------------------------------+
     ! subdivide into 2 children                                                    !
     if ( region_s%child(2)%uvbox(1) <= region_s%uvbox(1) + EPSregion ) then ! <--+ !
        call subdiv_bezier2_only_v( &                                             ! !
             region_s%poly(1)%ptr, &                                              ! !
             v=tuv_subdiv(3), &                                                   ! !
             bs=region_s%child(1)%poly(1)%ptr, &                                  ! !
             bn=region_s%child(2)%poly(1)%ptr )                                   ! !
     else ! ----------------------------------------------------------------------+ !
        call subdiv_bezier2_only_u( &                                             ! !
             region_s%poly(1)%ptr, &                                              ! !
             u=tuv_subdiv(2), &                                                   ! !
             bw=region_s%child(1)%poly(1)%ptr, &                                  ! !
             be=region_s%child(2)%poly(1)%ptr )                                   ! !
     end if ! <-------------------------------------------------------------------+ !
  end if ! <------------------------------------------------------------------------+

  ! carry on the recursion with new pairs of regions
  if ( all(nchild < 2) ) then
     ! erreur
     PRINT *,'************************************************'
     PRINT *,' TBOX =',REGION_C%UVBOX
     PRINT *,'UVBOX =',REGION_S%UVBOX
     PRINT *,'intersect_curve_surface : NO MORE SUBDIVISIONS !'
     PRINT *,'************************************************'
     STAT_DEGENERACY = 51
     RETURN
     !STOP
  end if

  do jchild = 1,nchild(2) ! <--------------------------------+
     if ( nchild(2) == 1 ) then ! <-----------------+        !
        newregion_s => region_s                     !        !
     else ! ----------------------------------------+        !
        newregion_s => region_s%child(jchild)       !        !
     end if ! <-------------------------------------+        !
     !                                                       !
     do ichild = 1,nchild(1) ! <-------------------------+   !
        if ( nchild(1) == 1 ) then ! <--------------+    !   !
           newregion_c => region_c                  !    !   !
        else ! -------------------------------------+    !   ! 
           newregion_c => region_c%child(ichild)    !    !   !
        end if ! <----------------------------------+    !   !
        !                                                !   !
        call intersect_curve_surface( &                  !   !
             root_c, &                                   !   !
             root_s, &                                   !   !
             newregion_c, &                              !   !
             newregion_s, &                              !   !
             tuvxyz, &                                   !   !
             ntuvxyz, &                                  !   !
             stat_degeneracy )                           !   !
        !                                                !   !
     end do ! <------------------------------------------+   !
     !                                                       !
  end do ! <-------------------------------------------------+

  nullify(newregion_s, newregion_c)

end subroutine intersect_curve_surface



subroutine newton_intersection_polyline( &
     surf, &
     lowerb, &
     upperb, &
     xyzp_prev, &
     htargetsqr, &
     stat, &
     uv, &
     xyzp )
  use mod_math
  use mod_linalg
  use mod_diffgeom
  use mod_tolerances  
  implicit none
  LOGICAL, PARAMETER :: DEBUG = ( GLOBALDEBUG .AND. .false. )
  integer,           parameter     :: itmax = 2 + ceiling(-log10(EPSuv))
  type(ptr_surface), intent(in)    :: surf(2)
  real(kind=fp),     intent(in)    :: lowerb(4)
  real(kind=fp),     intent(in)    :: upperb(4)
  real(kind=fp),     intent(in)    :: xyzp_prev(3)
  real(kind=fp),     intent(in)    :: htargetsqr
  integer,           intent(out)   :: stat
  real(kind=fp),     intent(inout) :: uv(2,2)
  real(kind=fp),     intent(out)   :: xyzp(3)
  real(kind=fp)                    :: xyz(3,2), resh
  real(kind=fp)                    :: r1(3), r2(3), jac(4,4), duv(4)
  real(kind=fp)                    :: cond, erruv
  integer                          :: it, isurf, ivar
  integer                           :: stat_refl
  
  IF ( DEBUG ) THEN
     PRINT *,''; PRINT *,'';
     PRINT *,'NEWTON_INTERSECTION_POLYLINE'
     PRINT *,'HTARGETSQR =',htargetsqr
     !PRINT *,'LOWERB =',LOWERB
     !PRINT *,'UPPERB =',UPPERB
  END IF

  stat = 1
  erruv = 0._fp
  cond = 1._fp

  do it = 1,itmax
     !IF ( DEBUG ) PRINT *,'NEWTON_INTERSECTION_POLYLINE, IT#',IT
     !IF ( DEBUG ) PRINT *,'UV =',UV
     !! compute residual
     do isurf = 1,2
        call eval( &
             xyz(:,isurf), &
             surf(isurf)%ptr, &
             uv(:,isurf) )
     end do
     IF ( IT == 1 ) THEN
        !PRINT *,'XYZ1 =',xyz(:,1)
        !PRINT *,'XYZ2 =',xyz(:,2)
     END IF
     r1 = xyz(:,1) - xyz(:,2)
     r2 = xyz(:,1) - xyzp_prev
     resh = sum(r2**2) - htargetsqr

     !! compute Jacobian matrix
     do isurf = 1,2
        do ivar = 1,2
           call evald1( &
                jac(1:3,2*(isurf-1)+ivar), &
                surf(isurf)%ptr, &
                uv(:,isurf), &
                ivar)
        end do
     end do
     jac(1:3,3:4) = -jac(1:3,3:4)
     do ivar = 1,2
        jac(4,ivar) = 2._fp*dot_product(jac(1:3,ivar), r2)
     end do
     jac(4,3:4) = 0._fp

     !! solve for Newton step
     !PRINT *,'RHS=',-R1,-RESH
     !PRINT *,'JAC='
     !CALL PRINT_MAT(JAC)
     call linsolve_svd( &
          duv, &
          jac, &
          -[r1,resh], &
          4, &
          4, &
          1, &
          cond )
     erruv = max(sum(duv(1:2)**2), sum(duv(3:4)**2))

     !PRINT *,'DUV =',DUV

     !! correct Newton step to keep the iterate inside feasible region
     call nd_box_reflexions( &
          reshape(uv, [4]), &
          lowerb, &
          upperb, &
          duv, &
          4, &
          stat_refl )
     if ( stat_refl > 0 ) return
          
     !PRINT *,'DUV* =',DUV
     
     !! update solution
     uv(:,1) = uv(:,1) + duv(1:2)
     uv(:,2) = uv(:,2) + duv(3:4)
     
     !! termination criteria
     IF ( DEBUG ) PRINT *,norm2(r1), sqrt(abs(resh)), sqrt(erruv), EPSfp*cond
     if ( erruv < max(EPSuvsqr, EPSfpsqr*cond**2) ) then ! <--------+
        if ( sum(r1**2) < EPSxyzsqr .and. &                         !
             abs(resh)  < tolhsqr ) then ! <---------------------+  !
           stat = 0                                              !  !
           do isurf = 1,2 ! <------------+                       !  !
              call eval( &               !                       !  !
                   xyz(:,isurf), &       !                       !  !
                   surf(isurf)%ptr, &    !                       !  !
                   uv(:,isurf) )         !                       !  !
           end do ! <--------------------+                       !  !
           xyzp = 0.5_fp * sum(xyz, 2)                           !  !
           IF ( DEBUG ) PRINT *,'CONVERGED, UV =',UV,', XYZ =',XYZP
        else ! --------------------------------------------------+  !
           IF ( DEBUG ) PRINT *,'STAGNATION'                     !  !
        end if ! ------------------------------------------------+  !
        return                                                      !
     end if ! <-----------------------------------------------------+

  end do

end subroutine newton_intersection_polyline



recursive subroutine intersect_simple_surfaces( &
     surfroot, &
     region, &
     param_vector, &
     interdata, &
     uvxyz, &
     nuvxyz, &
     stat_degeneracy )
  use mod_util
  use mod_math
  use mod_polynomial
  use mod_bernstein
  use mod_diffgeom
  use mod_regiontree
  use mod_types_intersection
  implicit none
  LOGICAL, PARAMETER :: DEBUG = ( GLOBALDEBUG .AND. .false. )
  type(ptr_surface),            intent(in)    :: surfroot(2)
  type(ptr_region),             intent(inout) :: region(2)
  real(kind=fp),                intent(in)    :: param_vector(3)
  integer,                      intent(inout) :: stat_degeneracy
  type(type_intersection_data), intent(inout) :: interdata
  real(kind=fp), allocatable,   intent(inout) :: uvxyz(:,:)
  integer,                      intent(inout) :: nuvxyz
  integer, allocatable                        :: ipts_ss(:)
  integer                                     :: npts_ss
  type(type_polynomial)                       :: regc
  type(type_curve)                            :: root_c
  integer, allocatable                        :: ipts_bs(:)
  integer                                     :: npts_bs
  real(kind=fp)                               :: uv_subdiv(2,2)
  integer                                     :: stat_subdiv(2), nchild(2)
  type(ptr_region)                            :: newregion(2)
  integer                                     :: stat_point(2)
  integer                                     :: order(2)
  integer                                     :: icurv, ivar, ival, ipt, ichild, jchild
  INTEGER :: ISURF

  if ( stat_degeneracy /= 0 ) return ! a degeneracy has been encountered

  IF ( DEBUG ) THEN
     !PAUSE
     PRINT *,''; PRINT *,'';
     PRINT *,'INTERSECT_SIMPLE_SURFACES'
     PRINT *,'UVBOXES ='
     DO ISURF = 1,2 ! <-----------------+
        PRINT *,REGION(ISURF)%PTR%UVBOX !
     END DO ! <-------------------------+
  END IF

  ! inherit from parent regions all intersection points contained in current regions
  if ( nuvxyz > 0 ) then ! <----------------------------+
     do isurf = 1,2 ! <---------------------------+     !
        call inherit_points( &                    !     !
             region(isurf)%ptr, &                 !     !
             uvxyz(2*isurf-1:2*isurf,1:nuvxyz), & !     !
             uvxyz(8,1:nuvxyz), &                 !     !
             nuvxyz )                             !     !
     end do ! <-----------------------------------+     !
  end if ! <--------------------------------------------+

  IF ( DEBUG ) THEN
     PRINT *,'IPTS ='
     DO ISURF = 1,2
        IF ( REGION(ISURF)%PTR%NPTS < 1 ) THEN
           PRINT *,'N/A'
        ELSE
           PRINT *,REGION(ISURF)%PTR%IPTS(1:REGION(ISURF)%PTR%NPTS)
        END IF
     END DO
  END IF

  ! get the list of already discovered intersection points contained in both surface regions
  npts_ss = 0
  if ( region(1)%ptr%npts > 0 .and. region(2)%ptr%npts > 0 ) then ! <-----------------------+
     call intersection_arrays( &                                                            !
          region(1)%ptr%ipts(1:region(1)%ptr%npts), &                                       !
          region(2)%ptr%ipts(1:region(2)%ptr%npts), &                                       !
          ipts_ss )                                                                         !
     if ( allocated(ipts_ss) ) npts_ss = size(ipts_ss)                                      !
  end if ! <--------------------------------------------------------------------------------+

  ! intersect the 4*2 pairs of border-surface
  npts_bs = 0
  outer : do icurv = 1,2 ! <------------------------------------------------+
     ! compute the parameterization of the surface region                   !
     ! of which we consider the border (Chebyshev polynomial basis)         !
     call chgvar2( &                                                        !
          surfroot(icurv)%ptr%x, &                                          !
          regc, &                                                           !
          region(icurv)%ptr%uvbox([1,3]), &                                 !
          region(icurv)%ptr%uvbox([2,4]) )                                  !
     !                                                                      !
     do ivar = 1,2 ! <---------------------------------------------------+  !
        !                                                                !  !
        do ival = 1,2 ! <---------------------------------------------+  !  !
           ! compute the parameterization of the border curve         !  !  !
           ! (still expressed in Chebyshev polynomial basis)          !  !  !
           call bivar2univar( &                                       !  !  !
                regc, &                                               !  !  !
                root_c%x, &                                           !  !  !
                ivar, &                                               !  !  !
                ival )                                                !  !  !
           ! perform degree reduction (without loss of accuracy)      !  !  !
           call economize1(root_c%x, EPSmath)                         !  !  !
           ! compute first and second derivatives of the curve        !  !  !
           call compute_deriv1(root_c)                                !  !  !
           call compute_deriv2(root_c)                                !  !  !
           !                                                          !  !  !
           ! compute border-surface intersection                      !  !  !
           call intersect_border_surface( &                           !  !  !
                surfroot, &                                           !  !  !
                root_c, &                                             !  !  !
                region, &                                             !  !  !
                icurv, &                                              !  !  !
                ivar, &                                               !  !  !
                ival, &                                               !  !  !
                ipts_ss, &                                            !  !  !
                npts_ss, &                                            !  !  !
                uvxyz, &                                              !  !  !
                nuvxyz, &                                             !  !  !
                ipts_bs, &                                            !  !  !
                npts_bs, &                                            !  !  !
                stat_degeneracy )                                     !  !  !
           !                                                          !  !  !
           if ( stat_degeneracy > 1 ) then ! <-------------------+    !  !  !
              ! a degeneracy has been encountered                !    !  !  !
              exit outer                                         !    !  !  !
           end if ! <--------------------------------------------+    !  !  !
           !                                                          !  !  !
        end do ! <----------------------------------------------------+  !  !
     end do ! <----------------------------------------------------------+  !
  end do outer ! <----------------------------------------------------------+
  if ( allocated(ipts_ss) ) deallocate(ipts_ss)

  ! free allocated polynomials
  call free_polynomial(regc      )
  call free_polynomial(root_c%x  )
  call free_polynomial(root_c%xt )
  call free_polynomial(root_c%xtt)

  ! append the new border-surface intersection points to each region's collection
  if ( npts_bs > 0 ) then ! <-----------------------------+
     do isurf = 1,2 ! <------------------------+          !
        call add_points_bottom_up( &           !          !
             region(isurf)%ptr, &              !          !
             ipts_bs(1:npts_bs), &             !          !
             npts_bs )                         !          !
     end do ! <--------------------------------+          !
  end if ! <----------------------------------------------+

  if ( stat_degeneracy > 1 ) then ! <---------------------+
     ! a degeneracy has been encountered                  !
     if ( allocated(ipts_bs) ) deallocate(ipts_bs)        !
     return                                               !
  end if ! <----------------------------------------------+

  IF ( DEBUG ) THEN
     PRINT *,'NPTS_BS =',NPTS_BS
     IF ( NPTS_BS > 0 ) THEN
        PRINT *,'IPTS_BS =',IPTS_BS
        CALL PRINT_MAT( TRANSPOSE(UVXYZ(1:7,IPTS_BS(1:NPTS_BS))) )
     END IF
  END IF

  if ( npts_bs == 0 ) then
     ! no border-surface intersection points
     if ( allocated(ipts_bs) ) deallocate(ipts_bs)
     return
  end if

  if ( npts_bs <= 2 ) then ! <------------------------------------------------+
     ! classify the border-surface intersection points (entering, exiting,    !
     ! isolated)                                                              !
     IF ( DEBUG ) THEN ! <---------------------------------------+            !
        PRINT *,'BSI POINTS: UV ='                               !            !
        CALL PRINT_MAT(TRANSPOSE(UVXYZ(1:4,IPTS_BS(1:NPTS_BS)))) !            !
     END IF ! <--------------------------------------------------+            !
     call classify_border_surface_intersection_point( &                       !
          surfroot, &                                                         !
          region, &                                                           !
          reshape( uvxyz(1:4,ipts_bs(1:npts_bs)), [2,2,npts_bs] ), &          !
          npts_bs, &                                                          !
          stat_point(1:npts_bs) )                                             !
     !                                                                        !
     !                                                                        !
     if ( npts_bs == 2 ) then ! <-------------------------------------+       !
        if ( stat_point(1)*stat_point(2) < 0 ) then ! <---------+     !       !
           ! 1 entering point and 1 exiting point.              !     !       !
           ! Re-order the points from entering to exiting       !     !       !
           if ( stat_point(1) < 0 ) then ! <---------------+    !     !       !
              order = [2,1]                                !    !     !       !
           else ! -----------------------------------------+    !     !       !
              order = [1,2]                                !    !     !       !
           end if ! <--------------------------------------+    !     !       !
           !                                                    !     !       !
           ! trace curve without ambiguity                      !     !       !
           call add_intersection_curve( &                       !     !       !
                interdata, &                                    !     !       !
                param_vector, &                                 !     !       !
                ipts_bs(order), &                               !     !       !
                reshape([region(1)%ptr%uvbox,&                  !     !       !
                region(2)%ptr%uvbox], [2,2,2]) )                !     !       !
           IF ( DEBUG ) PRINT *,'+1 INTERSECTION CURVE'
           !                                                    !     !       !
        elseif ( all(stat_point == 0) ) then ! -----------------+     !       !
           ! 2 isolated points                                  !     !       !
        else ! -------------------------------------------------+     !       !
           ! incorrect configuration                            !     !       !
           PRINT *,'------------------------------------------'
           PRINT *,'2 BSI POINTS - INCORRECT CONFIGURATION'
           PRINT *,'UVBOXES ='
           DO ISURF = 1,2 ! <-----------------+
              PRINT *,REGION(ISURF)%PTR%UVBOX !
           END DO ! <-------------------------+
           PRINT *,'STAT_POINT =',STAT_POINT
           PRINT *,'UVXYZ ='
           CALL PRINT_MAT( TRANSPOSE(UVXYZ(1:7,IPTS_BS(1:2))) )
           PRINT *,'TOLUV =',UVXYZ(8,IPTS_BS(1:2))
           CALL WRITE_POLYNOMIAL( SURFROOT(1)%PTR%X, 'dev_intersection/debugssi_surf1.cheb' )
           CALL WRITE_POLYNOMIAL( SURFROOT(2)%PTR%X, 'dev_intersection/debugssi_surf2.cheb' )
           CALL WRITE_POLYNOMIAL( REGION(1)%PTR%POLY(1)%PTR, 'dev_intersection/debugbsi_reg1.bern' )
           CALL WRITE_POLYNOMIAL( REGION(2)%PTR%POLY(1)%PTR, 'dev_intersection/debugbsi_reg2.bern' )
           PRINT *,'------------------------------------------'
           stat_degeneracy = 35                                 !     !       !
        end if ! <----------------------------------------------+     !       !
     end if ! <-------------------------------------------------------+       !
     !                                                                        !
     !                                                                        !
  elseif ( npts_bs > 2 ) then ! ----------------------------------------------+
     ! More than 2 border-surface intersection points have been found,        !
     ! the situation is ambiguous so we need to carry on the recursion.       !
     ! first, try to subdivide at an bsi point interior to at least one of    !
     ! the intersected surface regions.                                       !
     do ipt = 1,npts_bs ! <----------------------------------------------+    !
        do isurf = 1,2 ! <--------------------------------------------+  !    !
           uv_subdiv(:,isurf) = uvxyz(2*isurf-1:2*isurf,ipts_bs(ipt)) !  !    !
           call subdiv_region( &                                      !  !    !
                region(isurf)%ptr, &                                  !  !    !
                uv_subdiv(:,isurf), &                                 !  !    !
                stat_subdiv(isurf) )                                  !  !    !
        end do ! <----------------------------------------------------+  !    !
        if ( any(stat_subdiv /= 2) ) exit                                !    !
     end do ! <----------------------------------------------------------+    !
     !                                                                        !
     if ( all(stat_subdiv == 2) ) then ! <-------------------------------+    !
        ! all bsi points are at corners in both regions, subdivide at    !    !
        ! parametric midpoint                                            !    !
        do isurf = 1,2 ! <--------------------------------------------+  !    !
           uv_subdiv(:,isurf) = 0.5_fp * ( &                          !  !    !
                region(isurf)%ptr%uvbox([1,3]) + &                    !  !    !
                region(isurf)%ptr%uvbox([2,4]) )                      !  !    !
           call subdiv_region( &                                      !  !    !
                region(isurf)%ptr, &                                  !  !    !
                uv_subdiv(:,isurf), &                                 !  !    !
                stat_subdiv(isurf) )                                  !  !    !
        end do ! <----------------------------------------------------+  !    !
     end if ! <----------------------------------------------------------+    !
     !
     if ( all(stat_subdiv == 2) ) then ! <-------------------------------+    !
        ! error, at least one region should be subdivided                !    !
        stat_degeneracy = 37                                             !    !
        return                                                           !    !
     end if ! <----------------------------------------------------------+    !
     !                                                                        !
     ! compute Bezier control points for children regions                     !
     do isurf = 1,2 ! <--------------------------------------------------+    !
        if ( stat_subdiv(isurf) == 2 ) then ! <-----------------------+  !    !
           ! the subdivision point is at one of the surface's corners !  !    !
           nchild(isurf) = 1                                          !  !    !
           cycle                                                      !  !    !
        end if ! <----------------------------------------------------+  !    !
        ! the subdivision point is interior to the surface               !    !
        nchild(isurf) = size(region(isurf)%ptr%child)                    !    !
        if ( stat_subdiv(isurf) < 0 ) then ! <-----------+               !    !
           ! the region's children are already allocated !               !    !
           cycle                                         !               !    !
        end if ! <---------------------------------------+               !    !
        !                                                                !    !
        ! map uv_subdiv to the region's local frame ([0,1]^2)            !    !
        do ivar = 1,2 ! <--------------------------------+               !    !
           uv_subdiv(ivar,isurf) = 0.5_fp * ( 1._fp + &  !               !    !
                ab2n1p1( &                               !               !    !
                uv_subdiv(ivar,isurf), &                 !               !    !
                region(isurf)%ptr%uvbox(2*ivar-1), &     !               !    !
                region(isurf)%ptr%uvbox(2*ivar) ) )      !               !    !
        end do ! <---------------------------------------+               !    !
        !                                                                !    !
        ! allocate polynomials                                           !    !
        do ichild = 1,nchild(isurf) ! <--------------------------+       !    !
           allocate(region(isurf)%ptr%child(ichild)%poly(1)    ) !       !    !
           allocate(region(isurf)%ptr%child(ichild)%poly(1)%ptr) !       !    !
        end do ! <-----------------------------------------------+       !    !
        !                                                                !    !
        if ( stat_subdiv(isurf) == 0 ) then ! <----------------------+   !    !
           ! 4 children                                              !   !    !
           call subdiv_bezier2( &                                    !   !    !
                region(isurf)%ptr%poly(1)%ptr, &                     !   !    !
                uv_subdiv(:,isurf), &                                !   !    !
                bsw=region(isurf)%ptr%child(1)%poly(1)%ptr, &        !   !    !
                bse=region(isurf)%ptr%child(2)%poly(1)%ptr, &        !   !    !
                bnw=region(isurf)%ptr%child(3)%poly(1)%ptr, &        !   !    !
                bne=region(isurf)%ptr%child(4)%poly(1)%ptr )         !   !    !
        elseif ( stat_subdiv(isurf) == 1 ) then ! -------------------+   !    !
           ! 2 children                                              !   !    !
           if ( region(isurf)%ptr%child(2)%uvbox(1) <= &             !   !    !
                region(isurf)%ptr%uvbox(1) + EPSregion ) then ! <--+ !   !    !
              ! the subdivision point is on a iso-u boundary       ! !   !    !
              call subdiv_bezier2_only_v( &                        ! !   !    !
                   region(isurf)%ptr%poly(1)%ptr, &                ! !   !    !
                   v=uv_subdiv(2,isurf), &                         ! !   !    !
                   bs=region(isurf)%ptr%child(1)%poly(1)%ptr, &    ! !   !    !
                   bn=region(isurf)%ptr%child(2)%poly(1)%ptr )     ! !   !    !
           else ! -------------------------------------------------+ !   !    !
              ! the subdivision point is on a iso-v boundary       ! !   !    !
              call subdiv_bezier2_only_u( &                        ! !   !    !
                   region(isurf)%ptr%poly(1)%ptr, &                ! !   !    !
                   u=uv_subdiv(1,isurf), &                         ! !   !    !
                   bw=region(isurf)%ptr%child(1)%poly(1)%ptr, &    ! !   !    !
                   be=region(isurf)%ptr%child(2)%poly(1)%ptr )     ! !   !    !
           end if ! <----------------------------------------------+ !   !    !
        end if ! <---------------------------------------------------+   !    !
     end do ! -----------------------------------------------------------+    !
     !                                                                        !
     ! carry on the recursion with pairs of children regions                  !
     do jchild = 1,nchild(2) ! <--------------------------------------+       !
        if ( nchild(2) == 1 ) then ! <-------------------------+      !       !
           newregion(2)%ptr => region(2)%ptr                   !      !       !
        else ! ------------------------------------------------+      !       !
           newregion(2)%ptr => region(2)%ptr%child(jchild)     !      !       !
        end if ! <---------------------------------------------+      !       !
        !                                                             !       !
        do ichild = 1,nchild(1) ! <-------------------------------+   !       !
           if ( nchild(1) == 1 ) then ! <----------------------+  !   !       !
              newregion(1)%ptr => region(1)%ptr                !  !   !       !
           else ! ---------------------------------------------+  !   !       !
              newregion(1)%ptr => region(1)%ptr%child(ichild)  !  !   !       !
           end if ! <------------------------------------------+  !   !       !
           !                                                      !   !       !
           call intersect_simple_surfaces( &                      !   !       !
                surfroot, &                                       !   !       !
                newregion, &                                      !   !       !
                param_vector, &                                   !   !       !
                interdata, &                                      !   !       !
                uvxyz, &                                          !   !       !
                nuvxyz, &                                         !   !       !
                stat_degeneracy )                                 !   !       !
           !                                                      !   !       !
        end do ! <------------------------------------------------+   !       !
     end do ! <-------------------------------------------------------+       !
     !                                                                        !
  end if ! <------------------------------------------------------------------+

  if ( allocated(ipts_bs) ) deallocate(ipts_bs)
  
end subroutine intersect_simple_surfaces



subroutine intersect_elsewhere( &
     poly, &
     xyz, &
     separable, &
     randomize )
  use mod_math
  use mod_polynomial
  use mod_geometry
  use mod_separation
  implicit none
  LOGICAL, PARAMETER :: DEBUG = ( GLOBALDEBUG .AND. .false. )
  type(ptr_polynomial), intent(in)  :: poly(2)
  real(kind=fp),        intent(in)  :: xyz(3)
  logical,              intent(out) :: separable
  logical, optional,    intent(in)  :: randomize
  type(type_matrix)                 :: sep(2)
  integer                           :: nbcp, nsep(2)
  real(kind=fp)                     :: rot(3,3), vec(3)
  integer                           :: i

  ! rearrange Bezier control points
  do i = 1,2 ! <-----------------------------------------------------------+
     !                                                                     !                                                         
     if ( poly(i)%ptr%dim < 3 ) STOP 'intersect_elsewhere : poly%dim < 3'  !   
     !                                                                     !                                                         
     select case (poly(i)%ptr%nvar) ! <--------------------------------+   !
        !                                                              !   !
     case (1) ! -------------------------------------------------------+   !
        nbcp = poly(i)%ptr%degr(1) + 1                                 !   !
        allocate( sep(i)%mat(nbcp,3) )                                 !   !
        call rearrange_for_separability_test( &                        !   !
             poly(i)%ptr%coef(1:nbcp,1:3,1), &                         !   !
             nbcp, &                                                   !   !
             xyz, &                                                    !   !
             sep(i)%mat, &                                             !   !
             nsep(i), &                                                !   !
             .false. )                                                 !   !
        !                                                              !   !
     case (2) ! -------------------------------------------------------+   !
        nbcp = (poly(i)%ptr%degr(1) + 1) * (poly(i)%ptr%degr(2) + 1)   !   !
        allocate( sep(i)%mat(nbcp,3) )                                 !   !
        call rearrange_for_separability_test( &                        !   !
             reshape( poly(i)%ptr%coef(&                               !   !
             1:poly(i)%ptr%degr(1)+1,&                                 !   !
             1:poly(i)%ptr%degr(2)+1,&                                 !   !
             1:3), [nbcp,3] ), &                                       !   !
             nbcp, &                                                   !   !
             xyz, &                                                    !   !
             sep(i)%mat, &                                             !   !
             nsep(i), &                                                !   !
             .false. )                                                 !   !
        !                                                              !   !
     case default ! ---------------------------------------------------+   !
        STOP 'intersect_elsewhere : poly%nvar /= 1,2'                  !   !
     end select ! <----------------------------------------------------+   !
     !                                                                     !
     if ( nsep(i) >= nbcp ) then
        PRINT *,'******* intersect_elsewhere : NSEP > NBCP *******'
     end if
  end do ! <---------------------------------------------------------------+


  ! perform spatial separability test
  if ( any(nsep < 1) ) then ! <--------------------------------------------+
     !                                                                     !
     separable = .true.                                                    !
     !                                                                     !
  else ! ------------------------------------------------------------------+
     if ( present(randomize) ) then ! <-----------------------------+      !
        if ( randomize ) then ! <-------------------------------+   !      !
           call random_rotation_matrix3d( rot )                 !   !      !
           do i = 1,2 ! <-----------------------------------+   !   !      !
              sep(i)%mat(1:nsep(i),1:3) = &                 !   !   !      !
                   matmul( sep(i)%mat(1:nsep(i),1:3), rot ) !   !   !      !
           end do ! <---------------------------------------+   !   !      !
        end if ! <----------------------------------------------+   !      !
     end if ! <-----------------------------------------------------+      !
     !                                                                     !
     call separating_plane( &                                              !
          sep(1)%mat(1:nsep(1),1:3), &                                     !
          sep(2)%mat(1:nsep(2),1:3), &                                     !
          nsep(1), &                                                       !
          nsep(2), &                                                       !
          vec, &                                                           !
          separable )                                                      !
     !                                                                     !
  end if ! <---------------------------------------------------------------+
  
  IF ( DEBUG ) THEN
     !IF (.NOT.SEPARABLE) THEN
     PRINT *,'XYZ =',xyz
     IF (SEPARABLE) THEN
        PRINT *,'VEC = ',matmul(rot, vec)
     END IF
     CALL WRITE_POLYNOMIAL(POLY(1)%PTR, 'dev_intersection/sep1_poly.bern')
     CALL WRITE_POLYNOMIAL(POLY(2)%PTR, 'dev_intersection/sep2_poly.bern')
     CALL WRITE_MATRIX(SEP(1)%MAT(1:NSEP(1),1:3), NSEP(1), 3, 'dev_intersection/sep1.dat')
     CALL WRITE_MATRIX(SEP(2)%MAT(1:NSEP(2),1:3), NSEP(2), 3, 'dev_intersection/sep2.dat')
     do i = 1,2
        IF ( NSEP(I) >= SIZE(SEP(I)%MAT,1) ) THEN
           PRINT *,'I =',I
           PRINT *,'XYZ =',xyz
           STOP 'DEBUG SEPARATION'
        END IF
     END do
     !END IF
  END IF

  deallocate( sep(1)%mat, sep(2)%mat )

end subroutine intersect_elsewhere



subroutine read_intersection_curves( &
     filename, &
     surf, &
     interdata )
  use mod_util
  use mod_diffgeom
  use mod_types_intersection
  use mod_tolerances
  implicit none
  character(*),                 intent(in)    :: filename
  type(type_surface), target,   intent(in)    :: surf(:)
  type(type_intersection_data), intent(inout) :: interdata
  integer                                     :: fid, ncurves, np
  real(kind=fp), allocatable                  :: xyz(:,:), uv(:,:,:)
  type(ptr_surface)                           :: surfpair(2)
  integer                                     :: ic, ip, pair(2), endpt(2)

  call get_free_unit(fid)
  open( &
       unit=fid, &
       file=filename, &
       action='read' )
  read (fid,*) ncurves
  do ic = 1,ncurves
     read (fid,*) pair
     surfpair(1)%ptr => surf(pair(1))
     surfpair(2)%ptr => surf(pair(2))
     read (fid,*) np

     allocate(xyz(3,np), uv(2,2,np))
     do ip = 1,np
        read (fid,*) xyz(:,ip)
     end do
     do ip = 1,np
        read (fid,*) uv(:,1,ip), uv(:,2,ip)
     end do

     call add_intersection_point( &
          uv(:,:,1), &
          xyz(:,1), &
          surfpair, &
          2, &
          interdata, &
          endpt(1) ) 
     call add_intersection_point( &
          uv(:,:,np), &
          xyz(:,np), &
          surfpair, &
          2, &
          interdata, &
          endpt(2) )

     call add_intersection_curve( &
          interdata, &
          [0._fp, 0._fp, 0._fp], &
          endpt, &
          spread(spread([-1._fp-EPSuv, 1._fp+EPSuv], 2, 2), 3, 2) )
     interdata%curves(interdata%nc)%surf(1)%ptr => surf(pair(1))
     interdata%curves(interdata%nc)%surf(2)%ptr => surf(pair(2))
     interdata%curves(interdata%nc)%isplit(2,1:2) = [1,np]

     allocate(interdata%curves(interdata%nc)%polyline)
     interdata%curves(interdata%nc)%polyline%np = np
     call move_alloc(from=xyz, to=interdata%curves(interdata%nc)%polyline%xyz)
     call move_alloc(from=uv , to=interdata%curves(interdata%nc)%polyline%uv )

  end do
  close(fid)

end subroutine read_intersection_curves



subroutine intersect_gaussmaps_elsewhere( &
    sep, &
    nsep, &
    n_collineal, &
    separable, &
    stat, &
    vec )
  use mod_math
  use mod_geometry
  use mod_separation
  implicit none
  type(type_matrix), intent(inout) :: sep(2)
  integer,           intent(in)    :: nsep(2)
  real(kind=fp),     intent(in)    :: n_collineal(3)
  logical,           intent(out)   :: separable
  integer,           intent(out)   :: stat
  real(kind=fp),     intent(out)   :: vec(3)
  real(kind=fp)                    :: rot(3,3), wedge(2)
  integer                          :: isurf
  
  ! first rotate the frame to align the collineal normal with the x-axis
  rot(:,1) = n_collineal
  call complete_orthonormal_basis( rot(:,1), rot(:,2), rot(:,3) )
  do isurf = 1,2
     !PRINT *,'BEFORE, #',ISURF
     !CALL PRINT_MAT( sep(isurf)%mat(1:nsep(isurf),1:3) )
     sep(isurf)%mat(1:nsep(isurf),1:3) = &
          matmul( sep(isurf)%mat(1:nsep(isurf),1:3), rot )
  end do

  stat = 0
  if ( all(nsep < 1) ) then ! <----------------------------------------+
     stat = -2 ! degeneracy                                            !
     separable = .false.                                               !
     !                                                                 !
  elseif ( any(nsep < 1) ) then ! -------------------------------------+
     !                                                                 !
     isurf = maxloc( nsep, dim=1 )                                     !
     call minimal_bounding_sector( &                                   !
          sep(isurf)%mat(1:nsep(isurf),2:3), &                         !
          nsep(isurf), &                                               !
          wedge, &                                                     !
          separable )                                                  !
     if ( separable ) vec = [ 0._fp, cos(wedge(1)), sin(wedge(1)) ]    !
     !                                                                 !
  else ! --------------------------------------------------------------+
     !                                                                 !
     call separate_spherical_bounding_boxes( &                         !
          sep(1)%mat(1:nsep(1),1:3), &                                 !
          sep(2)%mat(1:nsep(2),1:3), &                                 !
          nsep(1), &                                                   !
          nsep(2), &                                                   !
          vec, &                                                       !
          separable, &                                                 !
          mask_axes=[.true., .false., .false.] )                       !
     !                                                                 !
  end if ! <-----------------------------------------------------------+

  ! rotate back to original frame
  do isurf = 1,2
     sep(isurf)%mat(1:nsep(isurf),1:3) = &
          matmul( sep(isurf)%mat(1:nsep(isurf),1:3), transpose(rot) )
     !PRINT *,'AFTER, #',ISURF
     !CALL PRINT_MAT( sep(isurf)%mat(1:nsep(isurf),1:3) )
  end do
  if ( separable ) vec = matmul( rot, vec )

end subroutine intersect_gaussmaps_elsewhere



recursive subroutine intersect_2Dpolylines( &
     xy, &
     head, &
     tail, &
     npts, &
     isegments, &
     lambda )
  use mod_util
  use mod_math
  implicit none
  real(kind=fp), parameter                  :: tol = 1.d-7
  real(kind=fp), parameter                  :: tolsqr = tol**2
  type(type_matrix),          intent(in)    :: xy(2)
  integer,                    intent(in)    :: head(2)
  integer,                    intent(in)    :: tail(2)
  integer,                    intent(inout) :: npts
  integer, allocatable,       intent(inout) :: isegments(:,:)
  real(kind=fp), allocatable, intent(inout) :: lambda(:,:)
  integer                                   :: nchild(2), ind(3,2)
  real(kind=fp)                             :: mat(2,2), rhs(2), t(2)
  logical                                   :: singular
  real(kind=fp)                             :: xybox(2,2,2,2) ! x/y, min/max, #child, #polyline
  real(kind=fp), dimension(2)               :: p, q, vec
  real(kind=fp)                             :: invvecsqr, dist
  integer                                   :: ntmp
  integer                                   :: i, j, ichild, jchild

  !if ( npts > 0 ) return

  do i = 1,2
     if ( tail(i) > head(i) + 1 ) then
        nchild(i) = 2
        ind(1:3,i) = [head(i), (head(i) + tail(i))/2, tail(i)]
     else
        nchild(i) = 1
        ind(1:2,i) = [head(i), tail(i)]
     end if
  end do

  if ( all(nchild == 1) ) then
     ! solve intersection of two segments
     mat(:,1) = xy(1)%mat(:,tail(1)) - xy(1)%mat(:,head(1))
     mat(:,2) = xy(2)%mat(:,head(2)) - xy(2)%mat(:,tail(2))
     rhs      = xy(2)%mat(:,head(2)) - xy(1)%mat(:,head(1))

     call solve_2x2( &
          t, &
          mat(1,1), &
          mat(1,2), &
          mat(2,1), &
          mat(2,2), &
          rhs, &
          singular )

     if ( .not.singular ) then
        !PRINT *,'intersect_2Dpolylines: SEGMENTS =',HEAD,', T =',T
        if ( minval(t) > -tol .and. maxval(t) - 1._fp < tol ) then
           t = min(1._fp, max(0._fp, t))
           ntmp = npts
           call append_vec( &
                head, &
                2, &
                isegments, &
                ntmp )
           call append_vec( &
                t, &
                2, &
                lambda, &
                npts )
           return
        end if
     end if

     IF ( .FALSE. ) THEN
        do i = 1,2 ! <------------------------------------+
           if ( i == 1 ) then ! <-------+                 !
              p = xy(i)%mat(:,head(i))  !                 !
           else ! ----------------------+                 !
              p = xy(i)%mat(:,tail(i))  !                 !
           end if ! <-------------------+                 !
           j = 1 + mod(i,2)                               !
           q = xy(j)%mat(:,head(j))                       !
           vec = mat(:,j)                                 !
           invvecsqr = 1._fp / sum(vec**2)                !
           do j = 1,2 ! <-----------------------------+   !
              t = dot_product(vec,p - q) * invvecsqr  !   !
              dist = sum((q + t*vec - p)**2)          !   !
              if ( dist < tolsqr ) then ! <---+       !   !
                 ntmp = npts                  !       !   !
                 call append_vec( &           !       !   !
                      head, &                 !       !   !
                      2, &                    !       !   !
                      isegments, &            !       !   !
                      ntmp )                  !       !   !
                 call append_vec( &           !       !   !
                      t, &                    !       !   !
                      2, &                    !       !   !
                      lambda, &               !       !   !
                      npts )                  !       !   !
                 return                       !       !   !
              end if ! <----------------------+       !   !
           end do ! <---------------------------------+   !
        end do ! <----------------------------------------+
     END IF

     return
  end if

  ! compute axis-aligned bounding boxes for children
  do i = 1,2
     do ichild = 1,nchild(i)
        xybox(:,1,ichild,i) = minval(xy(i)%mat(:,ind(ichild,i):ind(ichild+1,i)), dim=2)
        xybox(:,2,ichild,i) = maxval(xy(i)%mat(:,ind(ichild,i):ind(ichild+1,i)), dim=2)
     end do
  end do

  ! recurse with pairs of children
  do jchild = 1,nchild(2)
     do ichild = 1,nchild(1)
        if ( overlap_intervals(xybox(1,:,ichild,1), xybox(1,:,jchild,2), tol) .and. &
             overlap_intervals(xybox(2,:,ichild,1), xybox(2,:,jchild,2), tol) ) then
           call intersect_2Dpolylines( &
                xy, &
                [ind(ichild,1),   ind(jchild,2)  ], &
                [ind(ichild+1,1), ind(jchild+1,2)], &
                npts, &
                isegments, &
                lambda )
        end if
     end do
  end do

end subroutine intersect_2Dpolylines



recursive subroutine intersect_surface_pair( &
     surfroot, &
     region, &
     interdata, &
     uvxyz, &
     nuvxyz, &
     stat_degeneracy )
  use mod_util
  use mod_math
  use mod_bernstein
  use mod_diffgeom
  use mod_regiontree
  use mod_tolerances
  use mod_types_intersection
  implicit none
  LOGICAL, PARAMETER :: DEBUG = ( GLOBALDEBUG .AND. .true. )
  type(ptr_surface),            intent(in)    :: surfroot(2)
  type(ptr_region),             intent(inout) :: region(2)
  type(type_intersection_data), intent(inout) :: interdata
  real(kind=fp), allocatable,   intent(inout) :: uvxyz(:,:)
  integer,                      intent(inout) :: nuvxyz
  integer,                      intent(inout) :: stat_degeneracy
  logical                                     :: overlap, separable
  integer                                     :: stat_collineal
  real(kind=fp)                               :: uv_collineal(2,2), toluv
  real(kind=fp), dimension(3)                 :: xyz_collineal, n_collineal
  integer                                     :: stat_loopdetection
  real(kind=fp)                               :: tmp(8)
  real(kind=fp)                               :: uv_subdiv(2,2)
  type(ptr_polynomial)                        :: poly(2)
  real(kind=fp)                               :: param_vector(3)
  integer                                     :: stat_singularpoint
  real(kind=fp)                               :: dxyz_duv(3,2,2), n(3,2), dxyz_ds(3,2), DUV_DS(2,2,2), CURVATURE(2)
  type(ptr_region)                            :: newregion(2)
  integer, dimension(2)                       :: stat_subdiv, nchild
  integer                                     :: isurf, ipt, ivar, ichild, jchild, ipoly

  if ( stat_degeneracy /= 0 ) return ! a degeneracy has been encountered

  IF ( DEBUG ) THEN
     !PAUSE
     PRINT *,''; PRINT *,'';
     PRINT *,'INTERSECT_SURFACE_PAIR'
     PRINT *,'UVBOXES ='
     DO ISURF = 1,2 ! <-----------------+
        PRINT *,REGION(ISURF)%PTR%UVBOX !
     END DO ! <-------------------------+
  END IF

  ! inherit from parent regions all intersection points contained in current regions
  if ( nuvxyz > 0 ) then ! <----------------------------+
     do isurf = 1,2 ! <---------------------------+     !
        call inherit_points( &                    !     !
             region(isurf)%ptr, &                 !     !
             uvxyz(2*isurf-1:2*isurf,1:nuvxyz), & !     !
             uvxyz(8,1:nuvxyz), &                 !     !
             nuvxyz )                             !     !
     end do ! <-----------------------------------+     !
  end if ! <--------------------------------------------+

  IF ( DEBUG ) THEN
     PRINT *,'UV, TOL ='
     DO ISURF = 1,2
        PRINT *,'REGION #',ISURF
        IF ( REGION(ISURF)%PTR%NPTS > 0) THEN
           CALL PRINT_MAT(TRANSPOSE( &
                UVXYZ([2*ISURF-1,2*ISURF,8],REGION(ISURF)%PTR%IPTS(1:REGION(ISURF)%PTR%NPTS)) ))
        END IF
     END DO
  END IF


  ! compute bounding boxes for each region...
  do isurf = 1,2 ! <------------------------------------------------------+
     if ( .not.associated(region(isurf)%ptr%xyzbox) ) then ! <----+       !
        allocate( region(isurf)%ptr%xyzbox )                      !       !
        call bernOBB2( &                                          !       !
             region(isurf)%ptr%poly(1)%ptr%coef, &                !       !
             region(isurf)%ptr%poly(1)%ptr%degr, &                !       !
             region(isurf)%ptr%xyzbox )                           !       !
     end if ! <---------------------------------------------------+       !
  end do ! <--------------------------------------------------------------+

  ! ... and check if these bounding boxes overlap
  call overlap_OBBs( &
       region(1)%ptr%xyzbox, &
       region(2)%ptr%xyzbox, &
       overlap )

  if ( .not.overlap ) return ! disjoint bounding boxes => empty intersection

  ! check if there is a pair of collineal corner points
  call find_collineal_corners( &
       region, &
       stat_collineal, &
       uv_collineal, &
       n_collineal, &
       xyz_collineal )
  IF ( DEBUG ) PRINT *,'STAT_COLLINEAL(CORNERS) =',stat_collineal

  if ( stat_collineal <= 0 ) then ! <----------------------------+
     if ( stat_collineal < 0 ) then !<-----------------+         !
        if ( nuvxyz > 0) then ! <----------------+     !         !
           call check_unicity_with_tolerance( &  !     !         !
                reshape(uv_collineal, [4]), &    !     !         !
                xyz_collineal, &                 !     !         !
                toluv, &                         !     !         !
                4, &                             !     !         !
                uvxyz(1:4,1:nuvxyz), &           !     !         !
                uvxyz(5:7,1:nuvxyz), &           !     !         !
                uvxyz(8,1:nuvxyz), &             !     !         !
                nuvxyz, &                        !     !         !
                ipt )                            !     !         !
        else ! ----------------------------------+     !         !
           ipt = 1                               !     !         !
        end if ! <-------------------------------+     !         !
        !                                              !         !
        if ( ipt > nuvxyz ) then ! <-------------+     !         !
           tmp(1:2) = uv_collineal(:,1)          !     !         !
           tmp(3:4) = uv_collineal(:,2)          !     !         !
           tmp(5:7) = xyz_collineal              !     !         !
           tmp(8) = toluv                        !     !         !
           call append_vector( &                 !     !         !
                tmp(1:8), &                      !     !         !
                8, &                             !     !         !
                uvxyz, &                         !     !         !
                nuvxyz )                         !     !         !
        end if ! <-------------------------------+     !         !
        !                                              !         !
        do isurf = 1,2 ! <-----------------------+     !         !
           call add_points_bottom_up( &          !     !         !
                region(isurf)%ptr, &             !     !         !
                [nuvxyz], &                      !     !         !
                1 )                              !     !         !
        end do ! <-------------------------------+     !         !
        !                                              !         !
        do isurf = 1,2 ! <-----------------------+     !         !
           poly(isurf)%ptr => &                  !     !         !
                region(isurf)%ptr%poly(1)%ptr    !     !         !
        end do ! <-------------------------------+     !         !
        call intersect_elsewhere( &                    !         !
             poly, &                                   !         !
             xyz_collineal, &                          !         !
             separable )                               !         !
        IF ( DEBUG ) PRINT *,'SEPARABLE?',separable    !         !
        if ( separable ) return                        !         !
        !                                              !         !
     end if ! <----------------------------------------+         !
     toluv = EPSuv                                               !
     ! there is a pair of collineal corner points, we will       !
     ! subdivide both region at their parametric center point    !
     do isurf = 1,2 ! <----------------------------+             !
        uv_subdiv(:,isurf) = 0.5_fp * ( &          !             !
             region(isurf)%ptr%uvbox([1,3]) + &    !             !
             region(isurf)%ptr%uvbox([2,4]) )      !             !
     end do ! <------------------------------------+             !
  elseif ( stat_collineal < 2 ) then ! --------------------------+
     n_collineal(:) = 0._fp                                      !
  else ! --------------------------------------------------------+
     ! one surface has a singular corner                         !
     stat_degeneracy = 20                                        !
     return                                                      !
  end if ! <-----------------------------------------------------+


  ! check Hohmeyer's loop detection criterion
  do isurf = 1,2 ! <----------------------------------+
     poly(isurf)%ptr => region(isurf)%ptr%poly(2)%ptr !
  end do ! <------------------------------------------+
  call loop_detection_criterion( &
       poly, &
       stat_loopdetection, &
       param_vector, &
       ( stat_collineal <= 0 ), &
       n_collineal, &
       randomize=.false. )
  nullify(poly(1)%ptr, poly(2)%ptr)
  IF ( DEBUG ) PRINT *,'STAT_LOOPDETECTION =',stat_loopdetection

  if ( stat_loopdetection < 0 ) then ! <--------------+
     ! a degeneracy has been encountered              !
     stat_degeneracy = stat_loopdetection             !
     return                                           !
  end if ! <------------------------------------------+

  ! If the loop detection criterion IS satisfied, then the intersection of the 
  ! current surface regions is a (possibly empty) set of isolated points and/or 
  ! open curve branches. All the isolated points as well as the curves' endpoints
  ! are located on the boundary of a surface region.
  if ( stat_loopdetection == 0 ) then ! <---------------------------------------+
     ! first we start new temporary region trees, rooted at copies of current   !
     ! regions (all the data contained in the regions is copied, except the     !
     ! pointers to children)                                                    !
     do isurf = 1,2 ! <-------------------+                                     !
        allocate( newregion(isurf)%ptr )  !                                     !
        call copy_region( &               !                                     !
             region(isurf)%ptr, &         !                                     !
             newregion(isurf)%ptr )       !                                     !
     end do ! <---------------------------+                                     !
     !                                                                          !
     call intersect_simple_surfaces( &                                          !
          surfroot, &                                                           !
          newregion, &                                                          !
          param_vector, &                                                       !
          interdata, &                                                          !
          uvxyz, &                                                              !
          nuvxyz, &                                                             !
          stat_degeneracy )                                                     !
     !                                                                          !
     do isurf = 1,2 ! <------------------------------------------------------+  !
        ! copy the data in temporay regions back to the current regions      !  !
        ! (essentially indices of newly discovered intersection points)...   !  !
        if ( newregion(isurf)%ptr%npts > 0 ) then ! <---------------------+  !  !
           call append_n( &                                               !  !  !
                region(isurf)%ptr%ipts, &                                 !  !  !
                region(isurf)%ptr%npts, &                                 !  !  !
                newregion(isurf)%ptr%ipts(1:newregion(isurf)%ptr%npts), & !  !  !
                newregion(isurf)%ptr%npts, &                              !  !  !
                unique=.true. )                                           !  !  !
        end if ! <--------------------------------------------------------+  !  !
        ! ...then free the temporay region trees...                          !  !
        nullify( &                                                           !  !
             newregion(isurf)%ptr%xyzbox,      &                             !  !
             newregion(isurf)%ptr%poly(1)%ptr, &                             !  !
             newregion(isurf)%ptr%poly(2)%ptr )                              !  !
        deallocate( newregion(isurf)%ptr%poly )                              !  !
        call free_region_tree( newregion(isurf)%ptr )                        !  !
        deallocate( newregion(isurf)%ptr )                                   !  !
     end do ! <--------------------------------------------------------------+  !
     !                                                                          !
     ! ... and finally return (there cannot be other intersection points/curves ! 
     ! between the current surface regions)                                     !
     return                                                                     !
     !                                                                          !
  end if ! <--------------------------------------------------------------------+


  if ( stat_collineal > 0 ) then ! <--------------------------------------------+
     ! There are no pair of collineal corners have, so we now                   !
     ! search for a pair of interior collineal points.                          !
     ! First, set initial iterate to the parametric center of both regions ...  !
     do isurf = 1,2 ! <------------------------------------------+              !
        uv_collineal(:,isurf) = 0.5_fp * ( &                     !              !
             region(isurf)%ptr%uvbox([1,3]) + &                  !              !
             region(isurf)%ptr%uvbox([2,4]) )                    !              !
     end do ! <--------------------------------------------------+              !
     ! ... then run box-constrained Newton-Raphson algorithm                    !
     call find_collineal_points( &                                              !
          surfroot, &                                                           !
          [region(1)%ptr%uvbox([1,3]), region(2)%ptr%uvbox([1,3])] - EPSuv, &   !
          [region(1)%ptr%uvbox([2,4]), region(2)%ptr%uvbox([2,4])] + EPSuv, &   !
          stat_collineal, &                                                     !
          uv_collineal, &                                                       !
          n_collineal, &                                                        !
          toluv, &
          xyz_collineal )                                                       !
     IF ( DEBUG ) THEN
        PRINT *,'STAT_COLLINEAL =',stat_collineal
        IF ( stat_collineal <= 0 ) PRINT *,'UV_COLLINEAL =',UV_COLLINEAL
     END IF
     !                                                                          !
     ! if a pair of collineal points have been discovered, we will subdivide    !
     ! the surface regions at that point                                        !
     if ( stat_collineal <= 0 ) uv_subdiv = uv_collineal                        !
  end if ! <--------------------------------------------------------------------+


  if ( stat_collineal < 0 ) then  ! <-------------------------------------------+
     ! The collineal points are coincident: this is a tangential contact point. !
     ! first, add the point to the collection...                                !
     if ( nuvxyz > 0) then ! <-------------------------+                        !
        call check_unicity_with_tolerance( &           !                        !
             reshape(uv_collineal, [4]), &             !                        !
             xyz_collineal, &                          !                        !
             toluv, &                                  !                        !
             4, &                                      !                        !
             uvxyz(1:4,1:nuvxyz), &                    !                        !
             uvxyz(5:7,1:nuvxyz), &                    !                        !
             uvxyz(8,1:nuvxyz), &                      !                        !
             nuvxyz, &                                 !                        !
             ipt )                                     !                        !
     else ! -------------------------------------------+                        !
        ipt = 1                                        !                        !
     end if ! <----------------------------------------+                        !
     !                                                                          !
     if ( ipt > nuvxyz ) then ! <----------------------+                        !
        tmp(1:2) = uv_collineal(:,1)                   !                        !
        tmp(3:4) = uv_collineal(:,2)                   !                        !
        tmp(5:7) = xyz_collineal                       !                        !
        tmp(8) = toluv                                 !                        !
        call append_vector( &                          !                        !
             tmp(1:8), &                               !                        !
             8, &                                      !                        !
             uvxyz, &                                  !                        !
             nuvxyz )                                  !                        !
     end if ! <----------------------------------------+                        !
     !                                                                          !
     do isurf = 1,2 ! <--------------------------------+                        !
        call add_points_bottom_up( &                   !                        !
             region(isurf)%ptr, &                      !                        !
             [nuvxyz], &                               !                        !
             1 )                                       !                        !
     end do ! <----------------------------------------+                        !
     !                                                                          !
     ! ... then determine the nature of that singular point                     !
     ! compute tangent and normal vectors                                       !
     do isurf = 1,2 ! <--------------------------------+                        !
        do ivar = 1,2 ! <------------------+           !                        !
           call evald1( &                  !           !                        !
                dxyz_duv(:,ivar,isurf), &  !           !                        !
                surfroot(isurf)%ptr, &     !           !                        !
                uv_collineal(:,isurf), &   !           !                        !
                ivar )                     !           !                        !
        end do ! <-------------------------+           !                        !
        n(:,isurf) = cross(&                           !                        !
             dxyz_duv(:,1,isurf), &                    !                        !
             dxyz_duv(:,2,isurf))                      !                        !
     end do ! <----------------------------------------+                        !
     call characterize_tangential_intersection_point( &                         !
          surfroot, &                                                           !
          uv_collineal, &                                                       !
          dxyz_duv, &                                                           !
          n, &                                                                  !
          stat_singularpoint, &                                                 !
          dxyz_ds )                                                             !
     !                                                                          !
     select case (stat_singularpoint) ! <-------------------------------+       !
        case (1) ! -----------------------------------------------------+       !
           ! point on tangential intersection curve                     !       !
           PRINT *,'>>> POINT ON TAGENTIAL INTSERSECTION CURVE'         !       !
           call diffgeom_intersection( &
                surfroot, &
                uv_collineal, &
                duv_ds, &
                dxyz_ds, &
                stat_singularpoint, &
                curvature )
           !call diffgeom_intersection_curve( &
           !     surfroot, &
           !     uv_collineal, &
           !     duv_ds, &
           !     dxyz_ds, &
           !     stat_singularpoint, &
           !     curvature )
           PRINT *,'CURVATURE =',CURVATURE(1)
           PRINT *,'DXYZ_DS =',dxyz_ds(:,1)
           PRINT *,'DUV_DS='
           CALL PRINT_MAT(DUV_DS(:,1,:))
        case (2) ! -----------------------------------------------------+       !
           ! branch point (two curves meet at that point)               !       !
           stat_loopdetection = 3 ! force subdivision of both surfaces  !       !
           PRINT *,'>>> BRANCH POINT'                                   !       !
           call diffgeom_intersection( &
                surfroot, &
                uv_collineal, &
                duv_ds, &
                dxyz_ds, &
                stat_singularpoint, &
                curvature )
           !call diffgeom_intersection_curve( &
           !     surfroot, &
           !     uv_collineal, &
           !     duv_ds, &
           !     dxyz_ds, &
           !     stat_singularpoint, &
           !     curvature )
           PRINT *,'CURVATURE =',CURVATURE
           PRINT *,'DXYZ_DS ='
           CALL PRINT_MAT(dxyz_ds)
           PRINT *,'DUV_DS,1='
           CALL PRINT_MAT(DUV_DS(:,1,:))
           PRINT *,'DUV_DS,2='
           CALL PRINT_MAT(DUV_DS(:,2,:))
        case (3) ! -----------------------------------------------------+       !
           ! isolated tangential contact point                          !       !
           stat_loopdetection = 3 ! force subdivision of both surfaces  !       !
           PRINT *,'>>> ISOLATED TANGENTIAL CONTACT POINT'              !       !
        case (4) ! -----------------------------------------------------+       !
           ! high-order contact point                                   !       !
           PRINT *,'>>> HIGH-ORDER CONTACT POINT'                       !       !
     end select ! <-----------------------------------------------------+       !
     !                                                                          !
  elseif ( stat_collineal > 0 ) then ! -----------------------------------------+
     ! no pair of collineal points has been found, we subdivide the surface     !
     ! with largest Gauss map at its center point                               !
     !                                                                          !
     do isurf = 1,2 ! <------------------------------------------+              !
        uv_subdiv(:,isurf) = 0.5_fp * ( &                        !              !
             region(isurf)%ptr%uvbox([1,3]) + &                  !              !
             region(isurf)%ptr%uvbox([2,4]) )                    !              !
     end do ! <--------------------------------------------------+              !
     !                                                                          !
  end if ! <--------------------------------------------------------------------+

  ! subdivide the surface regions
  IF ( DEBUG ) PRINT *,'UV_SUBDIV =',UV_SUBDIV
  nchild(:) = 1
  stat_subdiv(:) = 2
  do isurf = 1,2 ! <------------------------------------------------------------+
     if ( stat_loopdetection == isurf .or. stat_loopdetection == 3 ) then ! <-+ !
        call subdiv_region( &                                                 ! !
             region(isurf)%ptr, &                                             ! !
             uv_subdiv(:,isurf), &                                            ! !
             stat_subdiv(isurf) )                                             ! !
        !                                                                     ! !
        if ( stat_subdiv(isurf) == 2 ) then ! <------------------------+      ! !
           ! the subdivision point is at one of the surface's corners  !      ! !
           nchild(isurf) = 1                                           !      ! !
        else ! --------------------------------------------------------+      ! !
           ! the subdivision point is interior to the surface          !      ! !
           nchild(isurf) = size(region(isurf)%ptr%child)               !      ! !
        end if ! <-----------------------------------------------------+      ! !
     end if ! <---------------------------------------------------------------+ !
  end do ! <--------------------------------------------------------------------+

  if ( all(nchild < 2) ) then  ! <---------------------------------+
     IF ( .FALSE. ) THEN
        PRINT *,'STAT_COLLINEAL =',STAT_COLLINEAL
        PRINT *,'UVBOXES ='
        DO ISURF = 1,2 ! <-----------------+
           PRINT *,REGION(ISURF)%PTR%UVBOX !
        END DO ! <-------------------------+
        PRINT *,'UV_COLLINEAL =',UV_COLLINEAL
        STOP '%%%%%%%%%%%%%%%%%%%'
     END IF
     if ( stat_loopdetection == 3 ) then ! <------------------+    !
        ! recursion terminates prematurely because both       !    !
        ! regions have ranges smaller than 2*EPSregion        !    !
        stat_degeneracy = 21                                  !    !
        return                                                !    !
     else ! --------------------------------------------------+    !
        isurf = 3 - stat_loopdetection                        !    !
        call subdiv_region( &                                 !    !
             region(isurf)%ptr, &                             !    !
             uv_subdiv(:,isurf), &                            !    !
             stat_subdiv(isurf) )                             !    !
        if ( stat_subdiv(isurf) == 2 ) then ! <-----------+   !    !
           ! recursion terminates prematurely because     !   !    !
           ! both regions are too small                   !   !    !
           stat_degeneracy = 22                           !   !    !
           return                                         !   !    !
        else ! -------------------------------------------+   !    !
           nchild(isurf) = size(region(isurf)%ptr%child)  !   !    !
        end if ! <----------------------------------------+   !    !
     end if ! <-----------------------------------------------+    !
  end if ! <-------------------------------------------------------+


  ! compute Bezier control points for children regions
  do isurf = 1,2 ! <------------------------------------------------------------+
     if ( nchild(isurf) < 2 ) cycle ! no children                               !
     if ( stat_subdiv(isurf) < 0 ) cycle ! the children are already allocated   !
     !                                                                          !
     ! map uv_subdiv to the local frame ([0,1]^2) of the region                 !
     do ivar = 1,2 ! <--------------------------------+                         !
        uv_subdiv(ivar,isurf) = 0.5_fp * ( 1._fp + &  !                         !
             ab2n1p1( &                               !                         !
             uv_subdiv(ivar,isurf), &                 !                         !
             region(isurf)%ptr%uvbox(2*ivar-1), &     !                         !
             region(isurf)%ptr%uvbox(2*ivar) ) )      !                         !
     end do ! <---------------------------------------+                         !
     !                                                                          !
     ! allocate polynomials                                                     !
     do ichild = 1,nchild(isurf) ! <---------------------------+                !
        allocate(region(isurf)%ptr%child(ichild)%poly(2)    )  !                !
        allocate(region(isurf)%ptr%child(ichild)%poly(1)%ptr)  !                !
        allocate(region(isurf)%ptr%child(ichild)%poly(2)%ptr)  !                !
     end do ! <------------------------------------------------+                !
     !                                                                          !
     if ( stat_subdiv(isurf) == 0 ) then ! <------------------------------+     !
        ! 4 children                                                      !     !
        do ipoly = 1,2 ! <------------------------------------------+     !     !
           call subdiv_bezier2( &                                   !     !     !
                region(isurf)%ptr%poly(ipoly)%ptr, &                !     !     !
                uv_subdiv(:,isurf), &                               !     !     !
                bsw=region(isurf)%ptr%child(1)%poly(ipoly)%ptr, &   !     !     !
                bse=region(isurf)%ptr%child(2)%poly(ipoly)%ptr, &   !     !     !
                bnw=region(isurf)%ptr%child(3)%poly(ipoly)%ptr, &   !     !     !
                bne=region(isurf)%ptr%child(4)%poly(ipoly)%ptr )    !     !     !
        end do ! <--------------------------------------------------+     !     !
        !                                                                 !     !
     elseif ( stat_subdiv(isurf) == 1 ) then ! ---------------------------+     !
        ! 2 children                                                      !     !
        if ( region(isurf)%ptr%child(2)%uvbox(1) <= &                     !     !
             region(isurf)%ptr%uvbox(1) + EPSregion ) then ! <---------+  !     !
           ! the subdivision point is on a iso-u boundary              !  !     !
           do ipoly = 1,2 ! <---------------------------------------+  !  !     !
              call subdiv_bezier2_only_v( &                         !  !  !     !
                   region(isurf)%ptr%poly(ipoly)%ptr, &             !  !  !     !
                   v=uv_subdiv(2,isurf), &                          !  !  !     !
                   bs=region(isurf)%ptr%child(1)%poly(ipoly)%ptr, & !  !  !     !
                   bn=region(isurf)%ptr%child(2)%poly(ipoly)%ptr )  !  !  !     !
           end do ! <-----------------------------------------------+  !  !     !
           !                                                           !  !     !
        else ! --------------------------------------------------------+  !     !
           ! the subdivision point is on a iso-v boundary              !  !     !
           do ipoly = 1,2 ! <---------------------------------------+  !  !     !
              call subdiv_bezier2_only_u( &                         !  !  !     !
                   region(isurf)%ptr%poly(ipoly)%ptr, &             !  !  !     !
                   u=uv_subdiv(1,isurf), &                          !  !  !     !
                   bw=region(isurf)%ptr%child(1)%poly(ipoly)%ptr, & !  !  !     !
                   be=region(isurf)%ptr%child(2)%poly(ipoly)%ptr )  !  !  !     !
           end do ! <-----------------------------------------------+  !  !     !
           !                                                           !  !     !
        end if ! <-----------------------------------------------------+  !     !
        !                                                                 !     !
     end if ! <-----------------------------------------------------------+     !
     !                                                                          !
  end do ! <--------------------------------------------------------------------+


  ! carry on the recursion with pairs of children regions
  do jchild = 1,nchild(2) ! <--------------------------------------+
     if ( nchild(2) == 1 ) then ! <-------------------------+      !
        newregion(2)%ptr => region(2)%ptr                   !      !
     else ! ------------------------------------------------+      !
        newregion(2)%ptr => region(2)%ptr%child(jchild)     !      !
     end if ! <---------------------------------------------+      !
     !                                                             !
     do ichild = 1,nchild(1) ! <-------------------------------+   !
        if ( nchild(1) == 1 ) then ! <----------------------+  !   !
           newregion(1)%ptr => region(1)%ptr                !  !   !
        else ! ---------------------------------------------+  !   !
           newregion(1)%ptr => region(1)%ptr%child(ichild)  !  !   !
        end if ! <------------------------------------------+  !   !
        !                                                      !   !
        call intersect_surface_pair( &                         !   !
             surfroot, &                                       !   !
             newregion, &                                      !   !
             interdata, &                                      !   !
             uvxyz, &                                          !   !
             nuvxyz, &                                         !   !
             stat_degeneracy )                                 !   !
        !                                                      !   !
     end do ! <------------------------------------------------+   !
  end do ! <-------------------------------------------------------+

end subroutine intersect_surface_pair



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



subroutine transfer_intersection_points( &
     from, &
     to )
  use mod_types_intersection
  implicit none
  type(type_intersection_data), intent(inout) :: from
  type(type_intersection_data), intent(inout) :: to
  integer                                     :: ip
  
  to%np = from%np

  do ip = 1,from%np
     to%points(ip)%xyz   =  from%points(ip)%xyz
     to%points(ip)%pos   => from%points(ip)%pos
     to%points(ip)%npos  =  from%points(ip)%npos
     to%points(ip)%ivert =  from%points(ip)%ivert
     nullify(from%points(ip)%pos)
  end do
  deallocate(from%points)
  from%np = 0

end subroutine transfer_intersection_points



subroutine refine_intersection_polyline( &
     curve, &
     tol2, &
     tol3 )
  use mod_math
  use mod_types_intersection
  use mod_diffgeom
  use mod_geometry
  USE MOD_UTIL
  implicit none
  LOGICAL, PARAMETER :: DEBUG = .false.
  type(type_intersection_curve), intent(inout) :: curve
  real(kind=fp),                 intent(in)    :: tol2
  real(kind=fp),                 intent(in)    :: tol3
  real(kind=fp)                                :: fractol2sqr, fractol3sqr
  type(type_intersection_polyline), pointer    :: polyline => null()
  real(kind=fp), dimension(4)                  :: lowerb, upperb
  real(kind=fp)                                :: uv(2,2), xyz(3), ctr(3), radsqr
  integer                                      :: stat
  logical                                      :: crit2, crit3
  integer                                      :: ivert, isurf
  INTEGER :: FID

  polyline => curve%polyline
  lowerb = reshape(curve%uvbox(1,1:2,1:2), [4])
  upperb = reshape(curve%uvbox(2,1:2,1:2), [4])

  fractol3sqr = 4._fp*(2._fp - tol3)*tol3**2
  fractol2sqr = 4._fp*(2._fp - tol2)*tol2**2
  
  IF ( DEBUG ) THEN
     CALL GET_FREE_UNIT(FID)
     OPEN(UNIT=FID, FILE='../debug/debug_rip_xyz0.dat', ACTION='WRITE')
     DO IVERT = 1,POLYLINE%NP
        WRITE (FID,*) POLYLINE%XYZ(1:3,IVERT)
     END DO
     CLOSE(FID)
     OPEN(UNIT=FID, FILE='../debug/debug_rip_uv0.dat', ACTION='WRITE')
     DO IVERT = 1,POLYLINE%NP
        WRITE (FID,*) POLYLINE%UV(1:2,1:2,IVERT)
     END DO
     CLOSE(FID)
  END IF
  
  
  ivert = 1
  outer_loop : do while ( ivert < polyline%np )
     IF ( DEBUG ) PRINT *,'IVERT =', IVERT, ', NP =', POLYLINE%NP
     !! get exact midpoint
     ! initial guess
     uv = 0.5_fp*(polyline%uv(1:2,1:2,ivert) + polyline%uv(1:2,1:2,ivert+1))
     ! relax to exact intersection
     call simultaneous_point_inversions( &
          curve%surf, &
          lowerb, &
          upperb, &
          stat, &
          uv, &
          xyz )
     if ( stat > 0 ) then
        PRINT *,'refine_intersection_polyline: failed to relax midpoint on segment #', ivert
        ivert = ivert + 1
        cycle outer_loop
     else
        IF ( DEBUG ) THEN
           PRINT *,'XYZ VERTICES ='
           CALL PRINT_MAT(TRANSPOSE(POLYLINE%XYZ(1:3,[IVERT,IVERT+1])))
           PRINT *,'XYZ MIDPOINT =',XYZ
           PRINT *,'UV VERTICES ='
           PRINT *,POLYLINE%UV(1:2,1:2,IVERT)
           PRINT *,POLYLINE%UV(1:2,1:2,IVERT+1)
           PRINT *,'UV MIDPOINT =',UV
        END IF
     end if
     !
     !! check chordal error
     ! in xyz space
     call circumcircle( &
       polyline%xyz(1:3,ivert), &
       xyz, &
       polyline%xyz(1:3,ivert+1), &
       ctr, &
       radsqr )
     !
     if ( radsqr < EPSfp ) radsqr = huge(1._fp)
     IF ( DEBUG ) PRINT *,'RAD3 =',SQRT(RADSQR)
     IF ( DEBUG ) PRINT *,SQRT(sum((polyline%xyz(1:3,ivert+1) - polyline%xyz(1:3,ivert))**2) / radsqr)
     crit3 = ( sum((polyline%xyz(1:3,ivert+1) - polyline%xyz(1:3,ivert))**2) < radsqr*fractol3sqr ) 
     !
     if ( crit3 ) then
        ! in uv space
        do isurf = 1,2
           call circumcircle( &
                [polyline%uv(1:2,isurf,ivert), 0._fp], &
                [uv(1:2,isurf), 0._fp], &
                [polyline%uv(1:2,isurf,ivert+1), 0._fp], &
                ctr, &
                radsqr )
           if ( radsqr < EPSfp ) radsqr = huge(1._fp)
           IF ( DEBUG ) PRINT *,'RAD2 =',SQRT(RADSQR)
           IF ( DEBUG ) PRINT *,SQRT(sum((polyline%uv(1:2,isurf,ivert+1) - polyline%uv(1:2,isurf,ivert))**2) / radsqr)
           crit2 = ( sum((polyline%uv(1:2,isurf,ivert+1) - polyline%uv(1:2,isurf,ivert))**2) < radsqr*fractol2sqr )
           if ( .not.crit2 ) exit
        end do
     else
        crit2 = .false.
     end if
     !
     IF ( DEBUG ) THEN
        PRINT *,'CRIT2, CRIT3 =', CRIT2, CRIT3
        PRINT *,''
     END IF
     if ( crit3 .and. crit2 ) then
        ivert = ivert + 1
        cycle outer_loop
     end if
     !
     !! insert midpoint
     call insert_polyline_point( &
          uv, &
          xyz, &
          stat, &
          polyline, &
          ivert )
     if ( stat == 0 ) then
        ! shift downstream polyline split points
        where ( curve%isplit(2,:) > ivert ) curve%isplit(2,:) = curve%isplit(2,:) + 1
     else
        PRINT *,'refine_intersection_polyline: failed to insert midpoint on segment #', ivert
        ivert = ivert + 1
        cycle outer_loop
     end if
     !
  end do outer_loop


  IF ( DEBUG ) THEN
     CALL GET_FREE_UNIT(FID)
     OPEN(UNIT=FID, FILE='../debug/debug_rip_xyz1.dat', ACTION='WRITE')
     DO IVERT = 1,POLYLINE%NP
        WRITE (FID,*) POLYLINE%XYZ(1:3,IVERT)
     END DO
     CLOSE(FID)
     OPEN(UNIT=FID, FILE='../debug/debug_rip_uv1.dat', ACTION='WRITE')
     DO IVERT = 1,POLYLINE%NP
        WRITE (FID,*) POLYLINE%UV(1:2,1:2,IVERT)
     END DO
     CLOSE(FID)

     PAUSE
  END IF
  
end subroutine refine_intersection_polyline



subroutine reallocate_polyline( &
     polyline, &
     np, &
     stat_alloc )
  use mod_math
  use mod_types_intersection
  ! Reallocates the uv and xyz arrays of an intersection_polyline
  ! to size (2,2,np) and (3,np), respectively.
  ! The actual length of the polyline (pline%np) is kept unchanged.
  implicit none
  type(type_intersection_polyline), intent(inout) :: polyline
  integer,                          intent(in)    :: np
  integer,                          intent(out)   :: stat_alloc
  real(kind=fp), allocatable                      :: uvtmp(:,:,:), xyztmp(:,:)

  if ( allocated(polyline%uv) ) then ! <--------------------------------------------+
     if ( size(polyline%uv,3) < np ) then ! <-----------------------------------+   !
        call move_alloc( polyline%uv, uvtmp )                                   !   !
        allocate( polyline%uv(2,2,np), stat=stat_alloc )                        !   !
        if ( stat_alloc == 0 ) polyline%uv(1:2,1:2,1:size(uvtmp,3)) = uvtmp     !   !
     else ! --------------------------------------------------------------------+   !
        stat_alloc = 0                                                          !   !
     end if ! <-----------------------------------------------------------------+   !
  else ! ---------------------------------------------------------------------------+
     allocate( polyline%uv(2,2,np), stat=stat_alloc )                               !
  end if ! <------------------------------------------------------------------------+

  if ( stat_alloc /= 0 ) then ! <-----------+
     stat_alloc = 1                         !
     return                                 !
  end if ! <--------------------------------+

  
  if ( allocated(polyline%xyz) ) then ! <-------------------------------------------+
     if ( size(polyline%xyz,2) < np ) then ! <----------------------------------+   !
        call move_alloc( polyline%xyz, xyztmp )                                 !   !
        allocate( polyline%xyz(3,np), stat=stat_alloc )                         !   !
        if ( stat_alloc == 0 ) polyline%xyz(1:3,1:size(xyztmp,2)) = xyztmp      !   !
     else ! --------------------------------------------------------------------+   !
        stat_alloc = 0                                                          !   !
     end if ! <-----------------------------------------------------------------+   !
  else ! ---------------------------------------------------------------------------+
     allocate( polyline%xyz(3,np), stat=stat_alloc )                                !
  end if ! <------------------------------------------------------------------------+

  if ( stat_alloc /= 0 ) then ! <-----------+
     stat_alloc = 2                         !
     return                                 !
  end if ! <--------------------------------+

end subroutine reallocate_polyline



subroutine find_collineal_corners( &
     region, &
     stat, &
     uv_collineal, &
     n_collineal, &
     xyz_collineal )
  use mod_math
  use mod_regiontree
  use mod_tolerances
  ! Searches for a pair of collineal points among the 4x4 pairs of corners 
  ! of two rectangular surface region, each one being described by a tensor-
  ! product grid of Bezier control points.
  ! Returns: stat > 0 if no such pair has been found;
  !               =-1 if a tangential contact point has been found;
  !               = 0 else (pair of non-coincident collineal points).
  !          uv_collineal  : uv-coordinates of the collineal points
  !          xyz_collineal : xyz-coordinates (relevent only if stat =-1)
  !          n_collineal   : the common (unit) normal direction at the collineal points.
  implicit none
  LOGICAL, PARAMETER :: DEBUG = ( GLOBALDEBUG .AND. .false. )
  type(ptr_region), intent(in)  :: region(2)
  integer,          intent(out) :: stat
  real(kind=fp),    intent(out) :: uv_collineal(2,2)
  real(kind=fp),    intent(out) :: n_collineal(3)
  real(kind=fp),    intent(out) :: xyz_collineal(3)
  real(kind=fp)                 :: s(3,2), n(3,2), r(3), nsqr
  integer                       :: i, j, k, l

  stat = 1

  do l = 1,2 ! <----------------------------------------------------------------------------------+
     do k = 1,2 ! <----------------------------------------------------------------------------+  !
        !                                                                                      !  !
        s(:,2) = region(2)%ptr%poly(1)%ptr%coef( &                                             !  !
             1 + (k-1)*region(2)%ptr%poly(1)%ptr%degr(1), &                                    !  !
             1 + (l-1)*region(2)%ptr%poly(1)%ptr%degr(2), &                                    !  !
             1:3)                                                                              !  !
        n(:,2) = region(2)%ptr%poly(2)%ptr%coef( &                                             !  !
             1 + (k-1)*region(2)%ptr%poly(2)%ptr%degr(1), &                                    !  !
             1 + (l-1)*region(2)%ptr%poly(2)%ptr%degr(2), &                                    !  !
             1:3)                                                                              !  !
        nsqr = sum(n(:,2)**2)                                                                  !  !
        if ( nsqr > EPSxyzsqr ) then ! <-----------------+                                     !  !
           n(:,2) = n(:,2) / sqrt(nsqr)                  !                                     !  !
        else ! ------------------------------------------+                                     !  !
           stat = 2                                      !                                     !  !
           return                                        !                                     !  !
        end if ! <---------------------------------------+                                     !  !
        !                                                                                      !  !
        do j = 1,2 ! <----------------------------------------------------------------------+  !  !
           do i = 1,2 ! <----------------------------------------------------------------+  !  !  !
              !                                                                          !  !  !  !
              s(:,1) = region(1)%ptr%poly(1)%ptr%coef( &                                 !  !  !  !
                   1 + (i-1)*region(1)%ptr%poly(1)%ptr%degr(1), &                        !  !  !  !
                   1 + (j-1)*region(1)%ptr%poly(1)%ptr%degr(2), &                        !  !  !  !
                   1:3)                                                                  !  !  !  !
              n(:,1) = region(1)%ptr%poly(2)%ptr%coef( &                                 !  !  !  !
                   1 + (i-1)*region(1)%ptr%poly(2)%ptr%degr(1), &                        !  !  !  !
                   1 + (j-1)*region(1)%ptr%poly(2)%ptr%degr(2), &                        !  !  !  !
                   1:3)                                                                  !  !  !  !
              nsqr = sum(n(:,1)**2)                                                      !  !  !  !
              if ( nsqr > EPSxyzsqr ) then ! <-----------------+                         !  !  !  !
                 n(:,1) = n(:,1) / sqrt(nsqr)                  !                         !  !  !  !
              else ! ------------------------------------------+                         !  !  !  !
                 stat = 2                                      !                         !  !  !  !
                 return                                        !                         !  !  !  !
              end if ! <---------------------------------------+                         !  !  !  !
              !                                                                          !  !  !  !
              r = s(:,1) - s(:,2)                                                        !  !  !  !
              !                                                                          !  !  !  !
              IF ( DEBUG ) THEN
                 PRINT '(I1,1X,I1,1X,I1,1X,I1,1X,E22.15,1X,E22.15)',I,J,K,L,&
                      NORM2(cross( n(:,1), r )**2), NORM2(cross( n(:,1), n(:,2) )**2)
              END IF
              if ( sum( cross( n(:,1), r )**2 ) + &                                      !  !  !  !
                   sum( cross( n(:,1), n(:,2) )**2 ) < EPScollinealsqr ) then ! <---+    !  !  !  !
                 if ( sum(r**2) < EPSxyzsqr ) then ! <-------------------------+    !    !  !  !  !
                    stat = -1                                                  !    !    !  !  !  !
                 else ! -------------------------------------------------------+    !    !  !  !  !
                    stat = 0                                                   !    !    !  !  !  !
                 end if ! <----------------------------------------------------+    !    !  !  !  !
                 uv_collineal(:,1) = region(1)%ptr%uvbox([i,2+j])                   !    !  !  !  !
                 uv_collineal(:,2) = region(2)%ptr%uvbox([k,2+l])                   !    !  !  !  !
                 xyz_collineal = 0.5_fp * sum( s, 2 )                               !    !  !  !  !
                 n_collineal = n(:,1) / norm2( n(:,1) )                             !    !  !  !  !
                 return                                                             !    !  !  !  !
              end if ! <------------------------------------------------------------+    !  !  !  !
              !                                                                          !  !  !  !
           end do ! <--------------------------------------------------------------------+  !  !  !
        end do ! <--------------------------------------------------------------------------+  !  !
     end do ! <--------------------------------------------------------------------------------+  !
  end do ! <--------------------------------------------------------------------------------------+

end subroutine find_collineal_corners
end module mod_intersection
