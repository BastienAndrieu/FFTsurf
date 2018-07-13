subroutine add_intersection_point( &
    uv, &
    xyz, &
    surf, &
    interdata, &
    id ) 
  use mod_math
  use mod_diffgeom2
  use mod_types_intersection
  use mod_tolerances
  implicit none
  LOGICAL, PARAMETER :: DEBUG = ( GLOBALDEBUG .AND. .TRUE. )
  integer, parameter                          :: PARAM_xtra_np = 10
  real(kind=fp),                intent(in)    :: uv(2,2)
  real(kind=fp),                intent(in)    :: xyz(3)
  type(ptr_surface),            intent(in)    :: surf(2)
  type(type_intersection_data), intent(inout) :: interdata
  integer,                      intent(out)   :: id
  type(type_point_on_surface), pointer        :: pos
  logical                                     :: newsurf(2)
  type(type_intersection_point), allocatable  :: tmp(:)
  integer                                     :: ipt, isurf, ipos

  do ipt = 1,interdata%np ! <------------------------------------------------------------------------+
     if ( interdata%points(ipt)%npos > 0 ) then ! <----------------------------------------------+   !
        if ( sum( (xyz - interdata%points(ipt)%xyz)**2 ) < EPSxyzsqr ) then ! <--------------+   !   !
           ! this is a duplicate point                                                       !   !   !
           pos => interdata%points(ipt)%pos                                                  !   !   !
           newsurf(:) = .true.                                                               !   !   !
           do ipos = 1,interdata%points(ipt)%npos ! <------------------------------------+   !   !   !
              do isurf = 1,2 ! <------------------------------------------------------+  !   !   !   !
                 if ( associated(pos%surf,surf(isurf)%ptr) ) newsurf(isurf) = .false. !  !   !   !   !
              end do ! <--------------------------------------------------------------+  !   !   !   !
              if ( ipos < interdata%points(ipt)%npos ) pos => pos%next                   !   !   !   !
           end do ! <--------------------------------------------------------------------+   !   !   !
           !                                                                                 !   !   !
           do isurf = 1,2 ! <---------------------------------------------------+            !   !   !
              if ( newsurf(isurf) ) then ! <---------------------------------+  !            !   !   !
                 interdata%points(ipt)%npos = interdata%points(ipt)%npos + 1 !  !            !   !   !
                 allocate( pos%next )                                        !  !            !   !   !
                 pos => pos%next                                             !  !            !   !   !
                 pos%uv = uv(:,isurf)                                        !  !            !   !   !
                 pos%surf => surf(isurf)%ptr                                 !  !            !   !   !
              end if ! <-----------------------------------------------------+  !            !   !   !
           end do ! <-----------------------------------------------------------+            !   !   !
           !                                                                                 !   !   !
           id = ipt                                                                          !   !   !
           nullify(pos)                                                                      !   !   !
           return                                                                            !   !   !
        end if ! <---------------------------------------------------------------------------+   !   !
     else ! -------------------------------------------------------------------------------------+   !
           STOP 'add_intersection_point : npos <= 0'                                             !   !
     end if ! <----------------------------------------------------------------------------------+   !
  end do ! <-----------------------------------------------------------------------------------------+

  if ( .not.allocated(interdata%points) ) allocate(interdata%points(PARAM_xtra_np) )

  ! this is a new point
  IF ( DEBUG ) THEN
     PRINT *,'ADDING INTERSECTION POINT :'
     PRINT *,'XYZ =',XYZ
     PRINT *,' UV =',UV
  END IF

  interdata%np = interdata%np + 1
  id = interdata%np
  if ( id > size(interdata%points) ) then ! <--------------------------+
     ! reallocate interdata%points                                     !
     allocate( tmp(size(interdata%points)) )                           !
     do ipt = 1,size(interdata%points) ! <-----------+                 !
        tmp(ipt)%xyz = interdata%points(ipt)%xyz     !                 !
        tmp(ipt)%npos = interdata%points(ipt)%npos   !                 !
        tmp(ipt)%pos => interdata%points(ipt)%pos    !                 !
        nullify( interdata%points(ipt)%pos )         !                 !
     end do ! <-------------------------------------+                  !
     deallocate( interdata%points )                                    !
     !                                                                 !
     allocate( interdata%points(interdata%np + PARAM_xtra_np) )        !
     do ipt = 1,size(tmp) ! <-----------------------+                  !
        interdata%points(ipt)%xyz = tmp(ipt)%xyz     !                 !
        interdata%points(ipt)%npos = tmp(ipt)%npos   !                 !
        interdata%points(ipt)%pos => tmp(ipt)%pos    !                 !
        nullify( tmp(ipt)%pos )                     !                  !
     end do ! <-------------------------------------+                  !
     deallocate( tmp )                                                 !
  end if ! <-----------------------------------------------------------+

  interdata%points(id)%xyz = xyz
  allocate( interdata%points(id)%pos )
  pos => interdata%points(id)%pos
  do isurf = 1,2 ! <---------------------------------------------------+
     interdata%points(id)%npos = interdata%points(id)%npos + 1         !
     pos%uv = uv(:,isurf)                                              !
     pos%surf => surf(isurf)%ptr                                       !
     if ( isurf < 2 ) allocate( pos%next )                             !
     pos => pos%next                                                   !
  end do ! <-----------------------------------------------------------+

  nullify(pos)

end subroutine add_intersection_point
