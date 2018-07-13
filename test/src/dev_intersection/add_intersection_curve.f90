subroutine add_intersection_curve( &
     interdata, &
     param_vector, &
     iendpoints, &
     uvbox )
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

   tmp%nc = 0
   if ( allocated(interdata%curves) ) then
      if ( interdata%nc + 1 > size(interdata%curves) ) then
      IF ( DEBUG ) PRINT *,'INTERDATA%CURVES --> TMP'
         allocate( tmp%curves(interdata%nc) )
         call transfer_intersection_curves( &
              from=interdata, &
              to=tmp )
      end if
   end if

   if ( .not.allocated(interdata%curves) ) then
      allocate( interdata%curves(tmp%nc + PARAM_xtra_nc) )
      if ( tmp%nc > 0 ) then
         IF ( DEBUG ) PRINT *,'INTERDATA%CURVES <-- TMP'
         call transfer_intersection_curves( &
              from=tmp, &
              to=interdata )
      end if
   end if

   interdata%nc = interdata%nc + 1
   ic = interdata%nc
   IF ( DEBUG ) PRINT *,'BEFORE, ASSOCIATED(ROOT)?',ASSOCIATED(interdata%curves(ic)%root)
   allocate(interdata%curves(ic)%root)
   IF ( DEBUG ) PRINT *,'AFTER, ASSOCIATED(ROOT)?',ASSOCIATED(interdata%curves(ic)%root)
   interdata%curves(ic)%uvbox          = uvbox
   interdata%curves(ic)%param_vector   = param_vector
   interdata%curves(ic)%root%endpoints = iendpoints
   IF ( DEBUG ) PRINT *,'ALL IS OK, NC =',interdata%nc

end subroutine add_intersection_curve
