subroutine add_intersection_curve( &
     interdata, &
     param_vector, &
     iendpoints, &
     uvbox )
  use mod_math
  use mod_types_intersection
  implicit none
  LOGICAL, PARAMETER :: DEBUG = ( GLOBALDEBUG .AND. .false. )
  integer, parameter                          :: PARAM_xtra_nc = 10
  type(type_intersection_data), intent(inout) :: interdata
  real(kind=fp),                intent(in)    :: param_vector(3)
  integer,                      intent(in)    :: iendpoints(2)
  real(kind=fp),                intent(in)    :: uvbox(2,2,2)
  type(type_intersection_data)                :: tmp
  integer                                     :: ic

  IF ( DEBUG ) THEN
     PRINT *,''; PRINT *,'';
     PRINT *,'ADD_INTERSECTION_CURVE'
     PRINT *,'INTERDATA%NC =',interdata%nc
     PRINT *,'ALLOCATED(INTERDATA%CURVES)?',allocated(interdata%curves)
     IF ( allocated(interdata%curves) ) PRINT *,'SIZE =',size(interdata%curves)
  END IF

  if ( .not.allocated(interdata%curves) ) allocate(interdata%curves(PARAM_xtra_nc) )

  ! reallocate if necessary
  if ( interdata%nc + 1 > size(interdata%curves) ) then
     allocate(tmp%curves(interdata%nc))
     IF ( DEBUG ) PRINT *,'INTERDATA%CURVES --> TMP'
     call transfer_intersection_curves( &
          from=interdata, &
          to=tmp )

     IF ( DEBUG ) PRINT *,'INTERDATA%CURVES <-- TMP'
     allocate(interdata%curves(tmp%nc + PARAM_xtra_nc))
     call transfer_intersection_curves( &
          from=tmp, &
          to=interdata )
  end if

  !! insert new curve
  interdata%nc = interdata%nc + 1
  ic = interdata%nc
  interdata%curves(ic)%uvbox          = uvbox
  interdata%curves(ic)%param_vector   = param_vector
  ! add endpoints as split points
  interdata%curves(ic)%nsplit = 2
  allocate(interdata%curves(ic)%isplit(2,2))
  interdata%curves(ic)%isplit(1,1:2) = iendpoints
  interdata%curves(ic)%isplit(2,1:2) = -1 ! do not forget to fill in when tracing the polyline!
  IF ( DEBUG ) PRINT *,'ALL IS OK, NC =',interdata%nc
  
end subroutine add_intersection_curve
