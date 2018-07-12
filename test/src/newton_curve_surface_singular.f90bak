program newton_curve_surface_singular
  
  use mod_util
  use mod_math
  use mod_polynomial
  use mod_diffgeom2

  implicit none
  
  type(type_curve)            :: curv
  type(type_surface)          :: surf
  real(kind=fp), dimension(3) :: tuv0, tuv, xyz, dxyz_dt, n, xyz_c, xyz_s
  real(kind=fp)               :: dxyz_duv(3,2), mat(3,3), cond, wt(3)
  real(kind=fp)               :: d2xyz_duv2(3,3), d2xyz_dt2(3), y(3)
  integer                     :: stat, ivar, fileunit, rank

  call get_free_unit( fileunit )
  open( unit=fileunit, file='newton_curve_surface_singular/tuv.dat', action='read' )
  read ( fileunit, * ) tuv0(1)
  read ( fileunit, * ) tuv0(2:3)
  close( fileunit )
  
  call read_polynomial( &
       curv%x, &
       'newton_curve_surface_singular/curv.cheb', &
       1, &
       1 )

  !CALL PRINT_MAT( CURV%X%COEF(:,:,1) )
  !STOP




  call read_polynomial( &
       surf%x, &
       'newton_curve_surface_singular/surf.cheb', &
       2, &
       1 )


  call write_polynomial( &
       curv%x, &
       'newton_curve_surface_singular/curv_verif.cheb' )
  call write_polynomial( &
       surf%x, &
       'newton_curve_surface_singular/surf_verif.cheb' )

  call compute_deriv1( curv )
  call compute_deriv2( curv )

  call compute_deriv1( surf )
  call compute_deriv2( surf )

  !tuv(:) = 0._fp
  call random_number( xyz )
  tuv = max( -1._fp, min( 1._fp, tuv0 + 0.8_fp * (2._fp * xyz - 1._fp ) ) )
  call newton_curve_surface( &
       curv, &
       surf, &
       [-1._fp, 1._fp], &
       [-1._fp, 1._fp, -1._fp, 1._fp], &
       tuv, &
       stat, &
       xyz )

  if ( stat == 0 ) then
     PRINT *,''
     PRINT *,'   TUV =',tuv
     PRINT *,'ERRTUV =',NORM2( tuv - tuv0 )
     call get_free_unit( fileunit )
     open( unit=fileunit, file='newton_curve_surface_singular/solution.dat', action='write' )
     write ( fileunit, * ) tuv
     close( fileunit )

     call eval( xyz_c, curv, tuv(1) )   ! curve
     call eval( xyz_s, surf, tuv(2:3) ) ! surface

     PRINT *,'   |R| =', NORM2( XYZ_C - XYZ_S )

     call evald1( &
          dxyz_dt, &
          curv, &
          tuv(1) )

     do ivar = 1,2
        call evald1( &
             dxyz_duv(:,ivar), &
             surf, &
             tuv(2:3), &
             ivar )
     end do
     n = cross( dxyz_duv(:,1), dxyz_duv(:,2) )
     n = n / norm2( n )

     PRINT *,'   XYZ =',XYZ
     PRINT *,'|GT.N| =',abs( dot_product( dxyz_dt, n ) ) / norm2( dxyz_dt )

     
     !print *,'J='
     !call print_mat( dxyz_duv )
     !print *,'n='
     !print *,n
     !print *,'[J,n]='
     !call print_mat( [dxyz_duv, n] )
     !mat = reshape( [dxyz_duv(:,1), dxyz_duv(:,2), n ], [3,3] )
     mat(:,1:2) = dxyz_duv
     mat(:,3) = n
     !call print_mat( mat )

     call linsolve( wt, mat, dxyz_dt, 3, 3, 1, rank, cond )
     print *,cond
     print *,wt


     call evald2( &
          d2xyz_dt2, &
          curv, &
          tuv(1) )

     do ivar = 1,3
        call evald2( &
             d2xyz_duv2(:,ivar), &
             surf, &
             tuv(2:3), &
             ivar )
     end do

     !y = d2xyz_dt2 - matmul( wt(1)*d2xyz_duv2(:,1:2) + wt(2)*d2xyz_duv2(:,2:3), wt(1:2) )
     y = d2xyz_dt2 - ( &
          d2xyz_duv2(:,1)*wt(1)**2 + &
          2._fp * d2xyz_duv2(:,2)*wt(1)*wt(2) + &
          d2xyz_duv2(:,3)*wt(2)**2 )
     PRINT *,' |Y.N| =',abs( dot_product( y, n ) ) / norm2( y )
     

  end if
  

contains

  subroutine newton_curve_surface( &
       curv, &
       surf, &
       tbox, &
       uvbox, &
       tuv, &
       stat, &
       xyz )
    use mod_util
    use mod_math
    use mod_diffgeom2
    use mod_tolerances    
    ! stat = 0 : converged
    !        1 : not converged
    !        2 : degeneracy
    implicit none
    real(kind=fp), parameter          :: THRESHOLD = real(1.e-2, kind=fp)
    real(kind=fp), parameter          :: EPS = EPSmath !real(1.e-12, kind=fp)
    real(kind=fp), parameter          :: TOL = 1.d-13!EPSxyz
    real(kind=fp), parameter          :: EPSsqr = EPS**2
    real(kind=fp), parameter          :: TOLsqr = TOL**2
    integer, parameter                :: nitmax = 50!ceiling(-log10(EPS))
    integer, parameter                :: nitcheck = 5

    type(type_curve),   intent(in)    :: curv
    type(type_surface), intent(in)    :: surf
    real(kind=fp),      intent(in)    :: tbox(2)
    real(kind=fp),      intent(in)    :: uvbox(4)
    real(kind=fp),      intent(inout) :: tuv(3)
    integer,            intent(out)   :: stat
    real(kind=fp),      intent(out)   :: xyz(3)
    real(kind=fp), dimension(3)       :: lowerb, upperb, rng
    real(kind=fp), dimension(3)       :: xyz_c, xyz_s, r
    real(kind=fp), dimension(3)       :: tuvtmp, dtuv
    real(kind=fp)                     :: rescheck, res, restmp
    real(kind=fp)                     :: jac(3,3), lambda
    !logical                           :: singular
    integer                           :: rank
    integer                           :: it

    integer                           :: fileunit
    real(kind=fp)                     :: hess(3,3), grad(3)
    real(kind=fp)                     :: cond, condH

    ! Feasible t,u,v-domain
    lowerb = [ tbox(1), uvbox([1,3]) ]
    upperb = [ tbox(2), uvbox([2,4]) ]
    rng = upperb - lowerb
    lowerb = lowerb - EPsuv*rng
    upperb = upperb + EPsuv*rng

    stat = 1
    restmp = huge(1._fp)
    
    call get_free_unit( fileunit )
    open( unit=fileunit, file='newton_curve_surface_singular/residu.dat', action='write' )

    PRINT *,'      |R|       |    ~cond(J)     |    ~cond(H)    '
    PRINT *,'----------------+-----------------+----------------'
    newton_iteration : do it = 1,nitmax

       ! position vector
       call eval( xyz_c, curv, tuv(1) )   ! curve
       call eval( xyz_s, surf, tuv(2:3) ) ! surface

       r = xyz_s - xyz_c ! residual vector
       res = sum( r**2 ) ! squared norm of residual vector
       !PRINT *,'NEWTON, IT#',IT,', RES =',REAL(NORM2(R))
       
       ! check signs of convergence
       if ( it == 1 ) rescheck = THRESHOLD * res
       if ( it > nitcheck .and. res > rescheck ) then
          PRINT *,'NO SIGN OF CONVERGENCE, STOP NEWTON ITERATION'
          close(fileunit)
          return
       end if

       
       ! Jacobian matrix
       call evald1( jac(:,1), curv, tuv(1) )
       jac(:,1) = -jac(:,1)
       call evald1( jac(:,2), surf, tuv(2:3), 1 )
       call evald1( jac(:,3), surf, tuv(2:3), 2 )

       
       ! *********
       call grad_hessian_curve_surface( &
            curv, &
            surf, &
            tuv, &
            hess, &
            grad )
       call linsolve( &
            dtuv, &
            hess, &
            -grad, &
            3, &
            3, &
            1, &
            rank, &
            condH )
       ! *********


       call linsolve( &
            dtuv, &
            jac, &
            -r, &
            3, &
            3, &
            1, &
            rank, &
            cond )

       !PRINT *,'DET(J) =',triple_product( jac(:,1),jac(:,2),jac(:,3) )
       !call QR_colpiv(jac, Qj, Rj, Pj, rank)
       !PRINT *,NORM2(R),triple_product( jac(:,1),jac(:,2),jac(:,3) ), abs(Rj(3,3))/abs(Rj(1,1)), rank
       !cond = abs( Rj(1,1) / Rj(3,3) )
       PRINT '(ES16.9,1x,A1,ES16.9,1x,A1,ES16.9)', NORM2(R), '|', cond, '|', condH
       !call print_mat( jac )
       !call print_mat( Rj )
       !print *,'cond =',abs( Rj(1,1) / Rj(3,3) )
       !print *,''
       if ( it == 1 ) then
          !PRINT *,NORM2(R), 1.0, 1.0
          !PRINT *,NORM2(R), 1.0
          write (fileunit,*) norm2(r), 1.d0, cond, tuv
       else
          !PRINT *, NORM2(R), abs(dtuv(1)), norm2(dtuv(2:3))
          !PRINT *, NORM2(R), norm2(dtuv)
          write (fileunit,*) norm2(r), norm2(dtuv), cond, tuv
       end if

       ! convergence criterion
       if ( res < TOLsqr ) then
          stat = 0
          if ( res < restmp ) then
             restmp = res
             tuvtmp = tuv
             xyz = 0.5_fp * ( xyz_s + xyz_c )
          end if
          if ( restmp < EPSsqr .and. sum(dtuv**2) < (2._fp*cond*epsilon(dtuv))**2 ) then
             tuv = tuvtmp
             !PRINT *,'CONVERGED'
             !PRINT *,restmp,' < ',EPSsqr, '  &  ', sum(dtuv**2),' < ',(condprev*epsilon(dtuv))**2
             close(fileunit)
             return
          end if
       end if

       ! solve for Newton step
       !call linsolve_QR( &
       !     dtuv, &
       !     jac, &
       !     -r, &
       !     3, &
       !     3, &
       !     singular, &
       !     rank )
       !condprev = cond

       if ( rank < 3 ) then ! degeneracy
          PRINT *,''
          PRINT *,'IT #',IT
          PRINT *,'TUV =',TUV
          PRINT *,'RES =',SQRT(RES)
          PRINT *,'JACOBIAN ='
          CALL PRINT_MAT( JAC )
          PRINT *,'RANK =',RANK
          PRINT *,'DTUV=',DTUV
          if ( stat == 0 ) then
             stat = 2
             return
          end if
       end if

       ! scale down Newton step to keep the solution inside feasible region
       call nd_box_constraint( &
            tuv, &
            lowerb, &
            upperb, &
            dtuv, &
            lambda )

       if ( lambda < -EPSmath ) return ! negative damped factor

       dtuv = lambda * dtuv
       if ( abs(dtuv(1)) < EPSuv .or. sum(dtuv(2:3)**2) < EPSuvsqr ) then
          ! damped Newton step is too small
          !PRINT *,'|DT| =',ABS(DTUV(1)),' |DUV| =',NORM2( DTUV(2:3) )
          !PRINT *,'NEWTON STEP TOO SMALL'
          !close(fileunit)
          !return
       end if

       ! update solution
       !if ( mod(it,3) == 0 ) lambda = 2._fp * lambda
       tuv = tuv + lambda * dtuv

    end do newton_iteration

    if ( stat == 0 ) tuv = tuvtmp

    close(fileunit)
    
  end subroutine newton_curve_surface







  subroutine grad_hessian_curve_surface( &
       curv, &
       surf, &
       tuv, &
       hess, &
       grad )
    type(type_curve),   intent(in)  :: curv
    type(type_surface), intent(in)  :: surf
    real(kind=fp),      intent(in)  :: tuv(3)
    real(kind=fp),      intent(out) :: hess(3,3)
    real(kind=fp),      intent(out) :: grad(3)
    real(kind=fp), dimension(3)     :: xyz_c, xyz_s, r
    real(kind=fp)                   :: dxyz_dt(3), dxyz_duv(3,2)
    real(kind=fp)                   :: d2xyz_dt2(3), d2xyz_duv2(3,3)
    integer                         :: ivar, jvar, i, j

    ! position vector
    call eval( xyz_c, curv, tuv(1) )   ! curve
    call eval( xyz_s, surf, tuv(2:3) ) ! surface

    r = xyz_s - xyz_c ! residual vector

    ! first derivatives
    call evald1( dxyz_dt, curv, tuv(1) )
    do ivar = 1,2
       call evald1( dxyz_duv(:,ivar), surf, tuv(2:3), ivar )
    end do

    ! gradient
    grad(1) = - dot_product( dxyz_dt, r )
    do ivar = 1,2
       grad(1+ivar) = dot_product( dxyz_duv(:,ivar), r )
    end do

    ! second derivatives
    call evald2( d2xyz_dt2, curv, tuv(1) )
    do ivar = 1,3
       call evald2( d2xyz_duv2(:,ivar), surf, tuv(2:3), ivar )
    end do

    ! Hessian matrix
    hess(1,1) = dot_product( dxyz_dt, dxyz_dt ) - dot_product( d2xyz_dt2, r )
    do ivar = 1,2
       hess(1+ivar,1) = - dot_product( dxyz_dt, dxyz_duv(:,ivar) )
       do jvar = ivar,2
          hess(1+jvar,ivar) = dot_product( dxyz_duv(:,ivar), dxyz_duv(:,jvar) ) + &
               dot_product( d2xyz_duv2(:,ivar+jvar-1), r )
       end do
    end do
    do j = 2,3
       do i = 1,j-1
          hess(i,j) = hess(j,i)
       end do
    end do

  end subroutine grad_hessian_curve_surface

end program newton_curve_surface_singular
