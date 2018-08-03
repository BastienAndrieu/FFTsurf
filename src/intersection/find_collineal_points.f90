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
             4 )
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
