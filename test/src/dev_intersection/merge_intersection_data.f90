subroutine merge_intersection_data( &
     surf, &
     uvxyz, &
     nuvxyz, &
     interdata_local, &
     interdata_global )
  use mod_math
  use mod_types_intersection
  ! Trace all intersection curves, intersect them and subidivide them accordingly and 
  implicit none
  type(ptr_surface),            intent(in)    :: surf(2)
  integer,                      intent(in)    :: nuvxyz
  real(kind=fp),                intent(in)    :: uvxyz(7,nuvxyz)
  type(type_intersection_data), intent(in)    :: interdata_local
  type(type_intersection_data), intent(inout) :: interdata_global
  integer                                     :: id_global(nuvxyz)
  integer                                     :: stat
  integer                                     :: nc
  integer                                     :: ip, ic, jc, isurf, jsurf

  ! add new intersection points
  do ip = 1,nuvxyz
     call add_intersection_point( &
          reshape(uvxyz(1:4,ip), [2,2]), &
          uvxyz(5:7,ip), &
          surf, &
          interdata_global, &
          id_global(ip) ) 
  end do

  nc = interdata_global%nc ! number of curves before adding new ones
  do ic = 1,interdata_local%nc
     ! add curve
     call add_intersection_curve( &
          interdata_global, &
          interdata_local%curves(ic)%param_vector, &
          id_global(interdata_local%curves(ic)%root%endpoints), &
          interdata_local%curves(ic)%uvbox )
     do isurf = 1,2
        interdata_global%curves(interdata_global%nc)%surf(isurf)%ptr => surf(isurf)%ptr
     end do

     ! trace polyline
     allocate(interdata_global%curves(interdata_global%nc)%polyline)
     call trace_intersection_polyline( &
          surf, &
          interdata_local%curves(ic)%uvbox, &
          interdata_local%curves(ic)%param_vector, &
          reshape(uvxyz(1:4,interdata_local%curves(ic)%root%endpoints),[2,2,2]), &
          uvxyz(5:7,interdata_local%curves(ic)%root%endpoints), &
          stat, &
          interdata_global%curves(interdata_global%nc)%polyline, &
          HMIN=REAL(1.E-3,KIND=FP), &
          HMAX=REAL(1.E-1,KIND=FP) )

     if ( stat > 0 ) then
        PRINT *,'STAT = ',STAT
        return!STOP
     end if

     ! check intersection with other curves
     do jc = 1,nc
        jloop : do jsurf = 1,2
           do isurf = 1,2
              if ( associated(&
                   interdata_global%curves(jc)%surf(jsurf)%ptr, &
                   surf(isurf)%ptr) ) then
                 if ( overlap_intervals(&
                      interdata_global%curves(jc)%uvbox(:,1,jsurf), &
                      interdata_local%curves(ic)%uvbox(:,1,isurf)) .and. &
                      overlap_intervals(&
                      interdata_global%curves(jc)%uvbox(:,2,jsurf), &
                      interdata_local%curves(ic)%uvbox(:,2,isurf))) then!exit jloop
                     PRINT *,'BLA'
                 end if
              end if
           end do
        end do jloop
        if ( jsurf > 2 ) cycle
     end do

  end do





end subroutine merge_intersection_data
