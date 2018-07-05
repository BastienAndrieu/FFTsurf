program newton_raphson



contains


  subroutine newton_curve_surface( &
       curv, &
       surf, &
       x0, &
       lowerb, &
       upperb, &
       EPSf, &
       EPSx, &
       stat, &
       xsol, &
       fsol )
    implicit none
    
    integer, parameter :: itcheckconv = 3

    type(type_curve),   intent(in)  :: curv
    type(type_surface), intent(in)  :: surf
    real(kind=fp),      intent(in)  :: uv0(3)
    real(kind=fp),      intent(in)  :: lowerb(3)
    real(kind=fp),      intent(in)  :: upperb(dim)
    real(kind=fp),      intent(in)  :: EPSxyz
    real(kind=fp),      intent(in)  :: EPSuv
    integer,            intent(out) :: stat
    real(kind=fp),      intent(out) :: uvsol(3)
    real(kind=fp),      intent(out) :: xyzsol(3)
    integer                         :: nitmax
    
    stat = 1
    convergence_order = 2
    it = 0
    nitmax = huge(1)
    do while( it < nitmax )
       !! compute "f"
       ! position vector
       call eval( xyz_c, curv, x(1) )   ! curve
       call eval( xyz_s, surf, x(2:3) ) ! surface
       r = xyz_s - xyz_c ! residual vector
       res = norm2( r )

       !! check convergence criteria
       if ( res < EPSxyz ) stat = 0
       if ( it > 1 ) then
          if ( erruv > EPSuv * condJ ) stat = 1
       end if
       if ( stat == 0 ) then
          xyzsol = 0.5_fp * ( xyz_s + xyz_c )
          return
       end if

       !! compute Jacobian matrix
       call evald1( jac(:,1), curv, tuv(1) )
       jac(:,1) = -jac(:,1)
       call evald1( jac(:,2), surf, tuv(2:3), 1 )
       call evald1( jac(:,3), surf, tuv(2:3), 2 )

       !! solve for Newton step using QR-factorization
       call linsolve( &
            dtuv, &
            jac, &
            -r, &
            3, &
            3, &
            1, &
            rankJ, &
            condJ )

       if ( it == itcheckconv ) then
          ! assess the type of convergence
          rate = normdtuv / normdtuvprev
          !if ( rate > 0.33_fp ) convergence_order = 1
       end if

       if ( convergence_order == 2 ) then
          erruv = normdtuv
       else
          nullspacedim = nint( mu / ( 1._fp - mu ) )
          erruv = ( mu normdtuv
       end if

    end do

    


  end subroutine newton_curve_surface


end program newton_raphson
