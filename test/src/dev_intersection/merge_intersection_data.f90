subroutine merge_intersection_data( &
     surf, &
     uvxyz, &
     nuvxyz, &
     interdata_local, &
     interdata_global )
  USE MOD_UTIL
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
  real(kind=fp)                               :: uvbox(2,2)
  integer                                     :: np(2)
  type(type_matrix)                           :: polylineuv(2)
  integer, allocatable                        :: isegm(:,:)
  real(kind=fp), allocatable                  :: lambda(:,:)
  integer                                     :: npts
  integer                                     :: ip, ic, jc, isurf, jsurf, ivar
  INTEGER :: FID, I

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
           iloop : do isurf = 1,2
              if ( associated(&
                   interdata_global%curves(jc)%surf(jsurf)%ptr, &
                   surf(isurf)%ptr) ) then
                 do ivar = 1,2
                    uvbox(:,ivar) = [ &
                         max( &
                         interdata_global%curves(jc)%uvbox(1,ivar,jsurf), &
                         interdata_local%curves(ic)%uvbox(1,ivar,isurf) ) , &
                         min( &
                         interdata_global%curves(jc)%uvbox(2,ivar,jsurf), &
                         interdata_local%curves(ic)%uvbox(2,ivar,isurf) ) ]
                    if ( uvbox(2,ivar) - uvbox(1,ivar) < epsilon(1._fp) ) cycle iloop
                 end do

                 np(1) = interdata_global%curves(nc+ic)%polyline%np
                 np(2) = interdata_global%curves(jc   )%polyline%np
                 allocate(polylineuv(1)%mat(2,np(1)), polylineuv(2)%mat(2,np(2)))
                 polylineuv(1)%mat(1:2,1:np(1)) = interdata_global%curves(nc+ic)%polyline%uv(1:2,isurf,1:np(1))
                 polylineuv(2)%mat(1:2,1:np(2)) = interdata_global%curves(jc   )%polyline%uv(1:2,jsurf,1:np(2))
                 npts = 0
                 call intersect_2Dpolylines( &
                      polylineuv, &
                      [1,1], &
                      np, &
                      npts, &
                      isegm, &
                      lambda )
                 deallocate(polylineuv(1)%mat, polylineuv(2)%mat)
 
                 IF ( NPTS > 0 ) THEN
                    PRINT *,'NPTS =',NPTS
                    PRINT *,'ISEGM ='
                    CALL PRINT_MAT( TRANSPOSE(ISEGM(:,1:NPTS)) )
                    PRINT *,'LAMBDA ='
                    CALL PRINT_MAT( TRANSPOSE(LAMBDA(:,1:NPTS)) )
                    PRINT *,''

                    CALL WRITE_POLYNOMIAL(surf(isurf)%ptr%x, 'dev_intersection/debugmrg_surf1.cheb')
                    CALL WRITE_POLYNOMIAL(surf(1+mod(isurf,2))%ptr%x, 'dev_intersection/debugmrg_surf2.cheb')
                    CALL WRITE_POLYNOMIAL(interdata_global%curves(jc)%surf(1+mod(jsurf,2))%ptr%x, &
                         'dev_intersection/debugmrg_surf3.cheb')

                    CALL GET_FREE_UNIT(FID)
                    OPEN(UNIT=FID, FILE='dev_intersection/merge_interdata.dat', ACTION='WRITE')
                    WRITE(FID,*) interdata_local%curves(ic)%uvbox(:,:,isurf)
                    WRITE(FID,*) interdata_global%curves(jc)%uvbox(:,:,jsurf)
                    WRITE(FID,*) uvbox
                    WRITE(FID,*) interdata_local%curves(ic)%uvbox(:,:,1+mod(isurf,2))
                    WRITE(FID,*) interdata_global%curves(jc)%uvbox(:,:,1+mod(jsurf,2))
                    WRITE(FID,*) interdata_global%curves(nc+ic)%polyline%np
                    DO I = 1,interdata_global%curves(nc+ic)%polyline%np
                       WRITE(FID,*) interdata_global%curves(nc+ic)%polyline%uv(:,isurf,i), &
                            interdata_global%curves(nc+ic)%polyline%uv(:,1+mod(isurf,2),i), &
                            interdata_global%curves(nc+ic)%polyline%xyz(:,i)
                    END DO
                    WRITE(FID,*) interdata_global%curves(jc)%polyline%np
                    DO I = 1,interdata_global%curves(jc)%polyline%np
                       WRITE(FID,*) interdata_global%curves(jc)%polyline%uv(:,jsurf,i), &
                            interdata_global%curves(jc)%polyline%uv(:,1+mod(jsurf,2),i), &
                            interdata_global%curves(jc)%polyline%xyz(:,i)
                    END DO
                    CLOSE(FID)
                    IF ( NPTS > 1 ) STOP
                 END IF

                 if ( allocated(isegm ) ) deallocate(isegm )
                 if ( allocated(lambda) ) deallocate(lambda)
                 
              end if
           end do iloop
        end do jloop
        if ( jsurf > 2 ) cycle
     end do

  end do





end subroutine merge_intersection_data
