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
  real(kind=fp)                 :: s(3,2), n(3,2), r(3)
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
