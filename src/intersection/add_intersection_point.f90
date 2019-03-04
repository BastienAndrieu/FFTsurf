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
        !PRINT *,'POINT #',ID,' -> +1 INTERSECTION POS, NPOS =',interdata%points(id)%npos
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

  !PRINT *,'+1 INTERSECTION POINT, NP =',INTERDATA%NP
  
  nullify(pos)

end subroutine add_intersection_point
