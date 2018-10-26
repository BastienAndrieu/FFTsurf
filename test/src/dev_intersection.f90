program dev_intersection

  use mod_util
  use mod_math
  use mod_polynomial
  use mod_diffgeom
  use mod_types_intersection
  use mod_intersection

  implicit none
  
  integer                    :: narg
  character(100)             :: arg
  integer, allocatable       :: numtest(:), numsurf(:)
  character                  :: strnum
  character(2)               :: strnum2
  integer*8                  :: tic, toc, count_rate

  logical, parameter                      :: ECONOMIZE = .true.
  integer                                 :: nsurf
  type(type_surface), target, allocatable :: surf(:)
  logical, allocatable                    :: mask(:)
  type(type_intersection_data)            :: interdata
  integer, allocatable                    :: surf2curv(:,:)
  integer                                 :: nsurf2curv
  integer                                 :: isurf, i, j, fid, ic, jsurf, ksurf

  ! =================================================================================
  ! Read command line arguments
  narg = command_argument_count()
  if ( narg < 1 ) then
     nsurf = 2
     allocate( numtest(nsurf), numsurf(nsurf) )
     numtest(:) = 1
     numsurf = [1,2]
  elseif ( mod(narg,2) /= 0 ) then
     nsurf = 2
     allocate( numtest(nsurf), numsurf(nsurf) )
     call get_command_argument(1, arg)
     read (arg,*) numtest(1)
     numtest(2) = numtest(1)
     numsurf = [1,2]
  else
     nsurf = narg / 2
     allocate( numtest(nsurf), numsurf(nsurf) )
     do isurf = 1,nsurf
        call get_command_argument(2*isurf-1, arg)
        read (arg,*) numtest(isurf)
        call get_command_argument(2*isurf, arg)
        read (arg,*) numsurf(isurf)
     end do
  end if


  PRINT *,'NUMTEST =',NUMTEST
  PRINT *,'NUMSURF =',NUMSURF

  
  ! =================================================================================
  ! Import surfaces
  allocate( surf(nsurf) )
  do isurf = 1,nsurf
     write (strnum2,'(I2.2)') numtest(isurf)
     write (strnum,'(I1)') numsurf(isurf)
     
     call read_polynomial( &
          surf(isurf)%x, &
          'coeffstest/C' // strnum // '_test' // strnum2 // '.txt', &
          nvar=2, &
          base=1 )

     if ( ECONOMIZE ) call economize2( surf(isurf)%x, EPSmath )

     call compute_deriv1( surf(isurf) )
     call compute_deriv2( surf(isurf) )
     call compute_pseudonormal( surf(isurf) )
     if ( ECONOMIZE ) call economize2( surf(isurf)%pn, EPSmath )

     write (strnum2,'(I2.2)') isurf
     CALL WRITE_POLYNOMIAL( surf(isurf)%x,  'dev_intersection/surfroot' // strnum2 // '_x.cheb'  )
     CALL WRITE_POLYNOMIAL( surf(isurf)%pn, 'dev_intersection/surfroot' // strnum2 // '_pn.cheb' )
  end do
  ! =================================================================================



  
  ! =================================================================================
  ! Compute all intersection
  allocate(mask(nsurf))
  mask(:) = .true.
  call system_clock( tic, count_rate )
  call intersect_all_surfaces( &
       surf, &
       nsurf, &
       interdata, &
       mask, &
       5.D-3, &
       1.D-3, &
       1.D-2 )
  call system_clock( toc )
  PRINT *,''; PRINT *,''; PRINT *,''
  PRINT *,'ELAPSED =',REAL( TOC - TIC ) / REAL( COUNT_RATE )

  call write_intersection_data( &
       interdata, &
       'dev_intersection/interdataglobal_points.dat', &
       'dev_intersection/interdataglobal_curves.dat' )

  PRINT *,'TOTAL :'
  PRINT *,INTERDATA%NP,' INTERSECTION POINT(S)'
  PRINT *,INTERDATA%NC,' INTERSECTION CURVE(S)'

  IF ( .false. ) THEN
     DO I = 1,INTERDATA%NC
        PRINT *,'- - - - - -'
        PRINT *,'CURV #',I
        PRINT *,'NP =',INTERDATA%CURVES(I)%POLYLINE%NP
        PRINT *,'NSPLIT =',INTERDATA%CURVES(I)%NSPLIT, ALLOCATED(INTERDATA%CURVES(I)%ISPLIT)
        DO J = 1,INTERDATA%CURVES(I)%NSPLIT
           PRINT *,INTERDATA%CURVES(I)%ISPLIT(:,J)
        END DO
     END DO
  END IF



  ! surface -> incident intersection curves (brute force)
  call get_free_unit( fid )
  open(unit=fid, file='dev_intersection/interdataglobal_surf2curv.dat', action='write')
  do isurf = 1,nsurf
     nsurf2curv = 0
     do ic = 1,interdata%nc
        do jsurf = 1,2
           if ( associated(interdata%curves(ic)%surf(jsurf)%ptr, surf(isurf)) ) then
              do ksurf = 1,nsurf
                 if ( associated(interdata%curves(ic)%surf(1+mod(jsurf,2))%ptr, surf(ksurf)) ) exit
              end do
              call insert_column_after( &
                   surf2curv, &
                   3, &
                   nsurf2curv, &
                   [ic, jsurf, ksurf], &
                   nsurf2curv )
           end if
        end do
     end do
     write (fid,*) nsurf2curv
     do ic = 1,nsurf2curv
        write (fid,*) surf2curv(:,ic)
     end do
     if ( allocated(surf2curv) ) deallocate(surf2curv)
  end do
  close(fid)


  

  call free_intersection_data(interdata)
  ! =================================================================================

  ! =================================================================================
  ! Free memory
  deallocate(numtest, numsurf)
  deallocate(surf, mask)
  call free_transformation_matrices()
  ! =================================================================================
  



contains
!end program dev_intersection


  








  subroutine classify_in_on_out( &
       interdata, &
       icurv, &
       class )
    ! class =  1 in
    !          0 on
    !         -1 out
    use mod_math
    use mod_diffgeom
    use mod_types_intersection
    use mod_intersection
    implicit none
    type(type_intersection_data), target, intent(in)    :: interdata
    integer,                              intent(in)    :: icurv
    integer, allocatable,                 intent(inout) :: class(:)
    type(type_point_on_surface), pointer                :: pos => null()
    real(kind=fp)                                       :: duv_ds(2,2,2), dxyz_ds(3,2)
    real(kind=fp)                                       :: dxyz_duv(3,2), n(3)
    integer                                             :: stat_point
    integer                                             :: is, ivar

    if ( allocated(class) ) deallocate(class)

    if ( interdata%curves(icurv)%nsplit < 3 ) return
    allocate(class(interdata%curves(icurv)%nsplit-1))
    class(:) = 1
    RETURN
    split_loop : do is = 2,interdata%curves(icurv)%nsplit-1
       ! get tangent to intersection curve at current split point
       call diffgeom_intersection_curve( &
            interdata%curves(icurv)%surf, &
            interdata%curves(icurv)%polyline%uv(:,:,interdata%curves(icurv)%isplit(2,is)), &
            duv_ds, &
            dxyz_ds, &
            stat_point )
       pos => interdata%points(interdata%curves(icurv)%isplit(1,is))%pos
       pos_loop : do while ( associated(pos) )
          if ( associated(pos%surf, interdata%curves(icurv)%surf(1)%ptr) .or. &
               associated(pos%surf, interdata%curves(icurv)%surf(2)%ptr) ) then
             pos => pos%next
             cycle pos_loop
          end if

          ! compute normal of current surface at curve-surface intersection point
          do ivar = 1,2
             call evald1(dxyz_duv(:,ivar), pos%surf, pos%uv, ivar)
          end do
          n = cross(dxyz_duv(:,1), dxyz_duv(:,2))

          ! check sign of dot product between the curve tangent and the surface normal
          if ( dot_product(dxyz_ds(:,1), n) > 0._fp ) then ! check small values...
             class(is) = -1
          else
             class(is-1) = -1
          end if

          pos => pos%next
       end do pos_loop

    end do split_loop
    
  end subroutine classify_in_on_out


  
  
end program dev_intersection
