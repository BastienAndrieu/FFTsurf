subroutine find_collineal_corners( &
     region, &
     uv, &
     ninter, &
     stat )
  use mod_math
  use mod_regiontree
  use mod_tolerances
  ! Searches for a pair of collineal points among the 4x4 pairs of corners 
  ! of two rectangular surface region, each one being described by a Bezier patch.
  ! Returns 'stat' > 0 if no such pair has been found;
  !         'stat' =-1 if a tangential contact point has been found;
  !         'stat' = 0 else (pair of non-coincident collineal points);
  ! as well as the uv-coordinates of those potential points.
  implicit none
  type(ptr_region), intent(in)  :: region(2)
  real(kind=fp),    intent(out) :: uv(2,2)
  real(kind=fp),    intent(out) :: ninter(3)
  integer,          intent(out) :: stat
  real(kind=fp)                 :: s(3,2), n(3,2), r(3)
  integer                       :: i, j, k, l

  stat = 1
  
  do l = 1,2 ! <----------------------------------------------------------------------------------+
     do k = 1,2 ! <----------------------------------------------------------------------------+  !
        !                                                                                      !  !
        !PRINT *,1 + (k-1)*region(2)%ptr%poly(1)%ptr%degr(1), 1 + (l-1)*region(2)%ptr%poly(1)%ptr%degr(2), &
        !     SIZE(region(2)%ptr%poly(1)%ptr%COEF,1), SIZE(region(2)%ptr%poly(1)%ptr%COEF,2)
        !PRINT *,region(2)%ptr%poly(1)%ptr%coef(1,1,1:3) 

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
              !PRINT *,'--'
              !PRINT *,'I,J,K,L ',I,J,K,L
              !PRINT *,'UV =',region(1)%ptr%uvbox([i,2+j]),region(2)%ptr%uvbox([k,2+l])
              !PRINT *,NORM2( CROSS(N(:,1),R) ), NORM2( CROSS(N(:,1),N(:,2)) )
              !PRINT *,I,J,K,L,NORM2( CROSS(N(:,1),R) ), NORM2( CROSS(N(:,1),N(:,2)) )
              if ( sum( cross( n(:,1), r )**2 ) + &                                      !  !  !  !
                   sum( cross( n(:,1), n(:,2) )**2 ) < EPScollinealsqr ) then ! <---+    !  !  !  !
                 if ( sum(r**2) < EPSxyzsqr ) then ! <-------------------------+    !    !  !  !  !
                    stat = -1                                                  !    !    !  !  !  !
                 else ! -------------------------------------------------------+    !    !  !  !  !
                    stat = 0                                                   !    !    !  !  !  !
                 end if ! <----------------------------------------------------+    !    !  !  !  !
                 uv(:,1) = region(1)%ptr%uvbox([i,2+j])                             !    !  !  !  !
                 uv(:,2) = region(2)%ptr%uvbox([k,2+l])                             !    !  !  !  !
                 ninter = n(:,1) / norm2( n(:,1) )                                  !    !  !  !  !
                 return                                                             !    !  !  !  !
              end if ! <------------------------------------------------------------+    !  !  !  !
              !                                                                          !  !  !  !
           end do ! <--------------------------------------------------------------------+  !  !  !
        end do ! <--------------------------------------------------------------------------+  !  !
     end do ! <--------------------------------------------------------------------------------+  !
  end do ! <--------------------------------------------------------------------------------------+

end subroutine find_collineal_corners
