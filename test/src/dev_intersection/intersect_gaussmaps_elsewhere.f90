subroutine intersect_gaussmaps_elsewhere( &
    sep, &
    nsep, &
    n_collineal, &
    separable, &
    stat, &
    vec )
  use mod_math
  use mod_geometry
  use mod_separation
  implicit none
  type(type_matrix), intent(inout) :: sep(2)
  integer,           intent(in)    :: nsep(2)
  real(kind=fp),     intent(in)    :: n_collineal(3)
  logical,           intent(out)   :: separable
  integer,           intent(out)   :: stat
  real(kind=fp),     intent(out)   :: vec(3)
  real(kind=fp)                    :: rot(3,3), wedge(2)
  integer                          :: isurf
  
  ! first rotate the frame to align the collineal normal with the x-axis
  rot(:,1) = n_collineal
  call complete_orthonormal_basis( rot(:,1), rot(:,2), rot(:,3) )
  do isurf = 1,2
     !PRINT *,'BEFORE, #',ISURF
     !CALL PRINT_MAT( sep(isurf)%mat(1:nsep(isurf),1:3) )
     sep(isurf)%mat(1:nsep(isurf),1:3) = &
          matmul( sep(isurf)%mat(1:nsep(isurf),1:3), rot )
  end do

  stat = 0
  if ( all(nsep < 1) ) then ! <----------------------------------------+
     stat = -2 ! degeneracy                                            !
     separable = .false.                                               !
     !                                                                 !
  elseif ( any(nsep < 1) ) then ! -------------------------------------+
     !                                                                 !
     isurf = maxloc( nsep, dim=1 )                                     !
     call minimal_bounding_sector( &                                   !
          sep(isurf)%mat(1:nsep(isurf),2:3), &                         !
          nsep(isurf), &                                               !
          wedge, &                                                     !
          separable )                                                  !
     if ( separable ) vec = [ 0._fp, cos(wedge(1)), sin(wedge(1)) ]    !
     !                                                                 !
  else ! --------------------------------------------------------------+
     !                                                                 !
     call separate_spherical_bounding_boxes( &                         !
          sep(1)%mat(1:nsep(1),1:3), &                                 !
          sep(2)%mat(1:nsep(2),1:3), &                                 !
          nsep(1), &                                                   !
          nsep(2), &                                                   !
          vec, &                                                       !
          separable, &                                                 !
          mask_axes=[.true., .false., .false.] )                       !
     !                                                                 !
  end if ! <-----------------------------------------------------------+

  ! rotate back to original frame
  do isurf = 1,2
     sep(isurf)%mat(1:nsep(isurf),1:3) = &
          matmul( sep(isurf)%mat(1:nsep(isurf),1:3), transpose(rot) )
     !PRINT *,'AFTER, #',ISURF
     !CALL PRINT_MAT( sep(isurf)%mat(1:nsep(isurf),1:3) )
  end do
  if ( separable ) vec = matmul( rot, vec )

end subroutine intersect_gaussmaps_elsewhere
