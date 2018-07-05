subroutine add_intersection_point( &
    uv, &
    xyz, &
    surf, &
    interdat, &
    id ) 
  use mod_math
  use mod_diffgeom2
  use mod_tolerances
  implicit none
  integer, parameter                          :: PARAM_xtra_np = 10
  real(kind=fp),                intent(in)    :: uv(2,2)
  real(kind=fp),                intent(in)    :: xyz(3)
  type(ptr_surface),            intent(in)    :: surf(2)
  type(type_intersection_data), intent(inout) :: interdat
  integer,                      intent(out)   :: id
  type(type_point_on_surface), pointer        :: pos
  logical                                     :: newsurf(2)
  type(type_intersection_point), allocatable  :: tmp(:)
  integer                                     :: ipt, isurf, ipos

  do ipt = 1,interdat%np ! <-------------------------------------------------------------------------+
     if ( interdat%points(ipt)%npos > 0 ) then ! <-----------------------------------------------+   !
        if ( sum( (xyz - interdat%points(ipt)%xyz)**2 ) < EPSxyzsqr ) then ! <---------------+   !   !
           ! this is a duplicate point                                                       !   !   !
           pos => interdat%points(ipt)%pos                                                   !   !   !
           newsurf(:) = .true.                                                               !   !   !
           do ipos = 1,interdat%points(ipt)%npos ! <-------------------------------------+   !   !   !
              do isurf = 1,2 ! <------------------------------------------------------+  !   !   !   !
                 if ( associated(pos%surf,surf(isurf)%ptr) ) newsurf(isurf) = .false. !  !   !   !   !
              end do ! <--------------------------------------------------------------+  !   !   !   !
              if ( ipos < interdat%points(ipt)%npos ) pos => pos%next                    !   !   !   !
           end do ! <--------------------------------------------------------------------+   !   !   !
           !                                                                                 !   !   !
           do isurf = 1,2 ! <-------------------------------------------------+              !   !   !
              if ( newsurf(isurf) ) then ! <-------------------------------+  !              !   !   !
                 interdat%points(ipt)%npos = interdat%points(ipt)%npos + 1 !  !              !   !   !
                 allocate( pos%next )                                      !  !              !   !   !
                 pos => pos%next                                           !  !              !   !   !
                 pos%uv = uv(:,isurf)                                      !  !              !   !   !
                 pos%surf => surf(isurf)%ptr                               !  !              !   !   !
              end if ! <---------------------------------------------------+  !              !   !   !
           end do ! <---------------------------------------------------------+              !   !   !
           !                                                                                 !   !   !
           id = ipt                                                                          !   !   !
           nullify(pos)                                                                      !   !   !
           return                                                                            !   !   !
        end if ! <---------------------------------------------------------------------------+   !   !
     else ! -------------------------------------------------------------------------------------+   !
           STOP 'add_intersection_point : npos <= 0'                                             !   !
     end if ! <----------------------------------------------------------------------------------+   !
  end do ! <-----------------------------------------------------------------------------------------+

  if ( .not.allocated(interdat%points) ) allocate(interdat%points(PARAM_xtra_np) )

  ! this is a new point
  PRINT *,'ADDING INTERSECTION POINT :'
  PRINT *,'XYZ =',XYZ
  PRINT *,' UV =',UV
  interdat%np = interdat%np + 1
  id = interdat%np
  if ( id > size(interdat%points) ) then ! <---------------------------+
     ! reallocate interdat%points                                      !
     allocate( tmp(size(interdat%points)) )                            !
     do ipt = 1,size(interdat%points) ! <-----------+                  !
        tmp(ipt)%xyz = interdat%points(ipt)%xyz     !                  !
        tmp(ipt)%npos = interdat%points(ipt)%npos   !                  !
        tmp(ipt)%pos => interdat%points(ipt)%pos    !                  !
        nullify( interdat%points(ipt)%pos )         !                  !
     end do ! <-------------------------------------+                  !
     deallocate( interdat%points )                                     !
     !                                                                 !
     allocate( interdat%points(interdat%np + PARAM_xtra_np) )          !
     do ipt = 1,size(tmp) ! <-----------------------+                  !
        interdat%points(ipt)%xyz = tmp(ipt)%xyz     !                  !
        interdat%points(ipt)%npos = tmp(ipt)%npos   !                  !
        interdat%points(ipt)%pos => tmp(ipt)%pos    !                  !
        nullify( tmp(ipt)%pos )                     !                  !
     end do ! <-------------------------------------+                  !
     deallocate( tmp )                                                 !
  end if ! <-----------------------------------------------------------+

  interdat%points(id)%xyz = xyz
  allocate( interdat%points(id)%pos )
  pos => interdat%points(id)%pos
  do isurf = 1,2 ! <---------------------------------------------------+
     interdat%points(id)%npos = interdat%points(id)%npos + 1           !
     pos%uv = uv(:,isurf)                                              !
     pos%surf => surf(isurf)%ptr                                       !
     if ( isurf < 2 ) allocate( pos%next )                             !
     pos => pos%next                                                   !
  end do ! <-----------------------------------------------------------+

  nullify(pos)

end subroutine add_intersection_point
