module mod_graph

  implicit none

contains

  subroutine disable_dangling_branches( &
       arc2nod, &
       nod2arc_in, &
       nin, &
       nod2arc_out, &
       nout, &
       active_arc, &
       active_nod, &
       narc, &
       nnod )
    use mod_util
    implicit none
    integer, intent(in)    :: narc, nnod
    integer, intent(in)    :: arc2nod(2,narc)
    integer, intent(inout) :: nin(nnod)
    integer, intent(inout) :: nout(nnod)
    integer, intent(inout) :: nod2arc_in(narc,nnod)
    integer, intent(inout) :: nod2arc_out(narc,nnod)
    logical, intent(inout) :: active_arc(narc)
    logical, intent(inout) :: active_nod(nnod)
    logical                :: change
    integer                :: inod, jnod, jarc, iarc
    
    do ! <------------------------------------------------------------+
       change = .false.                                               !
       do inod = 1,nnod ! <---------------------------------------+   !
          if ( .not.active_nod(inod) ) cycle                      !   !
          if ( nin(inod) < 1 .or. nout(inod) < 1 ) then ! <---+   !   !
             change = .true.                                  !   !   !
             !                                                !   !   !
             do jarc = 1,nin(inod) ! <----------+             !   !   !
                iarc = nod2arc_in(jarc,inod)    !             !   !   !
                active_arc(iarc) = .false.      !             !   !   !
                jnod = arc2nod(1,iarc)          !             !   !   !
                call remove_from_list( &        !             !   !   !
                     iarc, &                    !             !   !   !
                     nod2arc_out(:,jnod), &     !             !   !   !
                     nout(jnod) )               !             !   !   !
             end do ! <-------------------------+             !   !   !
             !                                                !   !   !
             do jarc = 1,nout(inod) ! <---------+             !   !   !
                iarc = nod2arc_out(jarc,inod)   !             !   !   !
                active_arc(iarc) = .false.      !             !   !   !
                jnod = arc2nod(2,iarc)          !             !   !   !
                call remove_from_list( &        !             !   !   !
                     iarc, &                    !             !   !   !
                     nod2arc_in(:,jnod), &      !             !   !   !
                     nin(jnod) )                !             !   !   !
             end do ! <-------------------------+             !   !   !
             !                                                !   !   !
             active_nod(inod) = .false.                       !   !   !
             !                                                !   !   !
          end if ! <------------------------------------------+   !   !
       end do ! <-------------------------------------------------+   !
       !                                                              !
       if ( .not.change ) return                                      !
    end do ! <--------------------------------------------------------+
    
  end subroutine disable_dangling_branches




  
  subroutine make_loops( &
       arc2nod, &
       arc_angles, &
       narc, &
       nnod, &
       nloops, &
       looparc, &
       !loopnod, &
       lenloop )
    use mod_util
    use mod_math
    implicit none
    integer,       intent(in)     :: narc, nnod
    integer,       intent(in)     :: arc2nod(2,narc)
    real(kind=fp), intent(in)     :: arc_angles(2,narc)
    integer,       intent(out)    :: nloops
    integer,       intent(out)    :: looparc(min(narc,nnod),min(narc,nnod))!, loopnod(min(narc,nnod),min(narc,nnod))
    integer,       intent(out)    :: lenloop(min(narc,nnod))
    integer, dimension(nnod)      :: nin, nout
    integer, dimension(narc,nnod) :: nod2arc_in, nod2arc_out
    logical                       :: active_arc(narc), active_nod(nnod)
    real(kind=fp)                 :: angle, dangle, maxangle_in, maxangle_out
    logical                       :: validloop
    integer                       :: inod, iarc, jarc, karc

    ! compute initial node -> arc incidence
    nin(:)  = 0
    nout(:) = 0
    do iarc = 1,narc ! <---------------------+
       inod = arc2nod(2,iarc)                !
       nin(inod) = nin(inod) + 1             !
       nod2arc_in(nin(inod), inod) = iarc    !
       !                                     !
       inod = arc2nod(1,iarc)                !
       nout(inod) = nout(inod) + 1           !
       nod2arc_out(nout(inod), inod) = iarc  !
    end do ! <-------------------------------+

    active_arc(:) = .true.
    active_nod(:) = .true.

    nloops = 0
    do karc = 1,narc ! <-------------------------------------------------------------------------+
       ! remove dangling branches                                                                !
       call disable_dangling_branches( &                                                         !
            arc2nod, &                                                                           !
            nod2arc_in, &                                                                        !
            nin, &                                                                               !
            nod2arc_out, &                                                                       !
            nout, &                                                                              !
            active_arc, &                                                                        !
            active_nod, &                                                                        !
            narc, &                                                                              !
            nnod )                                                                               !
       !                                                                                         !
       if ( .not.active_arc(karc) ) cycle                                                        !
       !                                                                                         !
       lenloop(nloops+1) = 1                                                                     !
       looparc(1,nloops+1) = karc                                                                !
       !loopnod(1,nloops+1) = arc2nod(1,karc)                                                     !
       inod = arc2nod(2,karc)                                                                    !
       validloop = .true.                                                                        !
       while_loop : do ! <------------------------------------------------------------------+    !
          angle = arc_angles(2,looparc(lenloop(nloops+1),nloops+1))                         !    !
          ! find leftmost ingoing arc at current node                                       !    !
          maxangle_in = -CSTpi                                                              !    !
          do jarc = 1,nin(inod) ! <----------------------------------------------------+    !    !
             if ( nod2arc_in(jarc,inod) == looparc(lenloop(nloops+1),nloops+1) ) cycle !    !    !
             dangle = diff_angle(arc_angles(2,nod2arc_in(jarc,inod)) + CSTpi, angle)   !    !    !
             maxangle_in = max(maxangle_in, dangle)                                    !    !    !
          end do ! <-------------------------------------------------------------------+    !    !
          !                                                                                 !    !
          ! find leftmost outgoing arc at current node                                      !    !
          maxangle_out = -CSTpi                                                             !    !
          iarc = 0                                                                          !    !
          do jarc = 1,nout(inod) ! <------------------------------------------+             !    !
             dangle = diff_angle(arc_angles(1,nod2arc_out(jarc,inod)), angle) !             !    !
             if ( dangle > maxangle_out ) then ! <-----+                      !             !    !
                iarc = nod2arc_out(jarc,inod)          !                      !             !    !
                maxangle_out = dangle                  !                      !             !    !
             end if ! <--------------------------------+                      !             !    !
          end do ! <----------------------------------------------------------+             !    !
          !                                                                                 !    !
          if ( maxangle_out < maxangle_in ) then ! <------------+                           !    !
             validloop = .false.                                !                           !    !
             exit while_loop                                    !                           !    !
          else ! -----------------------------------------------+                           !    !
             if ( iarc == looparc(1,nloops+1) ) exit while_loop !                           !    !
             lenloop(nloops+1) = lenloop(nloops+1) + 1          !                           !    !
             looparc(lenloop(nloops+1),nloops+1) = iarc         !                           !    !
             !loopnod(lenloop(nloops+1),nloops+1) = inod         !                           !    !
             inod = arc2nod(2,iarc)                             !                           !    !
          end if ! <--------------------------------------------+                           !    !
       end do while_loop ! <----------------------------------------------------------------+    !
       !                                                                                         !
       if ( validloop ) then ! <----------------------------------+                              !
          do jarc = 1,lenloop(nloops+1) ! <-------------------+   !                              !
             iarc = looparc(jarc,nloops+1)                    !   !                              !
             active_arc(iarc) = .false.                       !   !                              !
             call remove_from_list( &                         !   !                              !
                  iarc, &                                     !   !                              !
                  nod2arc_in(:,arc2nod(2,iarc)), &            !   !                              !
                  nin(arc2nod(2,iarc)) )                      !   !                              !
             call remove_from_list( &                         !   !                              !
                  iarc, &                                     !   !                              !
                  nod2arc_out(:,arc2nod(1,iarc)), &           !   !                              !
                  nout(arc2nod(1,iarc)) )                     !   !                              !
          end do ! <------------------------------------------+   !                              !
          nloops = nloops + 1                                     !                              !
       end if ! <-------------------------------------------------+                              !
       !                                                                                         !
    end do ! <-----------------------------------------------------------------------------------+
    
  end subroutine make_loops

  



  subroutine make_faces( &
       nascendants, &
       ascendants, &
       nloop, &
       outer, &
       inner, &
       ninner, &
       nfaces )
    implicit none
    integer, intent(in)  :: nloop
    integer, intent(in)  :: nascendants(nloop)
    integer, intent(in)  :: ascendants(nloop,nloop)
    integer, intent(out) :: outer(nloop)
    integer, intent(out) :: inner(nloop,nloop)
    integer, intent(out) :: ninner(nloop)
    integer, intent(out) :: nfaces
    integer              :: nchildren(nloop), children(nloop,nloop)
    integer              :: maxlevel, level
    integer              :: iloop, jloop, kloop
    
    maxlevel = 0
    do iloop = 1,nloop
       maxlevel = max(maxlevel, nascendants(iloop))
    end do

    nchildren(:) = 0
    do level = 0,maxlevel ! <------------------------------------------------+
       do iloop = 1,nloop ! <---------------------------------------------+  !
          if ( nascendants(iloop) == level ) then ! <------------------+  !  !
             do kloop = 1,level ! <---------------------------------+  !  !  !
                jloop = ascendants(kloop,iloop)                     !  !  !  !
                if ( nascendants(jloop) == level - 1 ) then ! <--+  !  !  !  !
                   nchildren(jloop) = nchildren(jloop) + 1       !  !  !  !  !
                   children(nchildren(jloop),jloop) = iloop      !  !  !  !  !
                   exit                                          !  !  !  !  !
                end if ! <---------------------------------------+  !  !  !  !
             end do ! <---------------------------------------------+  !  !  !
          end if ! <---------------------------------------------------+  !  !
       end do ! <---------------------------------------------------------+  !
    end do ! <---------------------------------------------------------------+

    nfaces = 0
    do iloop = 1,nloop ! <----------------------------------------------------------+
       if ( mod(nascendants(iloop),2) == 0 ) then ! <----------------------------+  !
          nfaces = nfaces + 1                                                    !  !
          outer(nfaces) = iloop                                                  !  !
          ninner(nfaces) = nchildren(iloop)                                      !  !
          inner(1:nchildren(iloop),nfaces) = children(1:nchildren(iloop),iloop)  !  !
       end if ! <----------------------------------------------------------------+  !
    end do ! <----------------------------------------------------------------------+

  end subroutine make_faces

  
end module mod_graph
