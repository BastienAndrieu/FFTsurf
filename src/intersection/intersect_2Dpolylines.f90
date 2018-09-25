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
