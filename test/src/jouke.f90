program jouke

  use mod_util
  use mod_math
  use mod_polynomial
  use mod_diffgeom
  use mod_types_intersection
  use mod_intersection

  implicit none
  
  character(3)               :: strnum3
  integer                    :: fid
  integer*8                  :: tic, toc, count_rate

  logical, parameter                      :: ECONOMIZE = .true.
  integer                                 :: nsurf
  type(type_surface), target, allocatable :: surf(:)
  logical, allocatable                    :: mask(:)
  type(type_intersection_data)            :: interdata
  integer, allocatable                    :: surf2curv(:,:)
  integer                                 :: nsurf2curv
  integer                                 :: isurf, ic, jsurf, ksurf

  ! =================================================================================
  ! Read propergol brep
  call get_free_unit(fid)
  open(unit=fid, file='Jouke/propergol/brep.dat', action='read')
  read(fid,*) !'n_Pa'
  read(fid,*) nsurf
  !read(fid,*) !'G1 continuity'
  !read(fid,*) nsmooth
  !do ic = 1,nsmooth
  !   read(fid,*) isurf, jsurf
  !   call add_intersection_curve( &
  !        interdata, &
  !        [0._fp, 0._fp, 0._fp], &
  !        [0, 0], &
  !        spread(spread([-1._fp, 1._fp], 2, 2), 3, 2) )
  !   interdata%curves(ic)%dummy = .true.
  !   interdata%curves(ic)%smooth = .true.
  !end do
  close(fid)
  !PRINT *,'';PRINT *,'';PRINT *,'';
  
  ! =================================================================================
  ! Import surfaces
  allocate( surf(nsurf) )
  do isurf = 1,nsurf
     write (strnum3,'(I3.3)') isurf
     call read_polynomial( &
          surf(isurf)%x, &
          'Jouke/propergol/C_' // strnum3 // '.cheb', &
          nvar=2, &
          base=1 )

     if ( ECONOMIZE ) call economize2( surf(isurf)%x, EPSmath )

     call compute_deriv1( surf(isurf) )
     call compute_deriv2( surf(isurf) )
     call compute_pseudonormal( surf(isurf) )
     if ( ECONOMIZE ) call economize2( surf(isurf)%pn, EPSmath )
  end do
  ! =================================================================================

  call read_curves('Jouke/propergol/tangent_curves.dat', surf, interdata)
  do ic = 1,interdata%nc
     interdata%curves(ic)%smooth = .true.
     !PRINT *,'NP = ',INTERDATA%CURVES(IC)%POLYLINE%NP, ALLOCATED(INTERDATA%CURVES(IC)%POLYLINE%UV)
     !interdata%curves(ic)%dummy = .true.
  end do
  !PRINT *,INTERDATA%NC,' DUMMY CURVES'
  !STOP
  ! =================================================================================
  ! Compute all intersection curves
  allocate(mask(nsurf))
  if (.false.) then
     mask(:) = .false.
     mask(1:7) = .true.
     mask(28:44) = .true.
  else
     mask(:) = .true.
  end if
  call system_clock( tic, count_rate )
  call intersect_all_surfaces( &
       surf, &
       nsurf, &
       interdata, &
       mask )
  call system_clock( toc )
  PRINT *,''; PRINT *,''; PRINT *,''
  PRINT *,'ELAPSED =',REAL( TOC - TIC ) / REAL( COUNT_RATE )

  call write_intersection_data( &
       interdata, &
       'Jouke/result/interdata_points.dat', &
       'Jouke/result/interdata_curves.dat' )

  PRINT *,'TOTAL :'
  PRINT *,INTERDATA%NP,' INTERSECTION POINT(S)'
  PRINT *,INTERDATA%NC,' INTERSECTION CURVE(S)'


  ! surface -> incident intersection curves (brute force)
  call get_free_unit( fid )
  open(unit=fid, file='Jouke/result/interdata_surf2curv.dat', action='write')
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
  deallocate(surf, mask)
  call free_transformation_matrices()
  ! =================================================================================
  



contains
!end program dev_intersection


  subroutine write_intersection_data( &
       interdata, &
       filepoints, &
       filecurves )
    use mod_util
    use mod_types_intersection
    implicit none
    type(type_intersection_data), intent(in) :: interdata
    character(*),                 intent(in) :: filepoints, filecurves
    integer                                  :: fileunit
    !integer, allocatable                     :: class(:)
    integer                                  :: ip, ic, is

    call get_free_unit( fileunit )
    
    ! intersection points
    open( &
         unit = fileunit, &
         file = filepoints, &
         action = 'write' )
    do ip = 1,interdata%np
       write ( fileunit, * ) interdata%points(ip)%xyz
    end do
    close( fileunit )

    ! intersection curves
    open( &
         unit = fileunit, &
         file = filecurves, &
         action = 'write' )
    write ( fileunit, * ) interdata%nc
    do ic = 1,interdata%nc
       !call classify_in_on_out( &
       !     interdata, &
       !     ic, &
       !     class )
       write ( fileunit, * ) logic2int(interdata%curves(ic)%dummy)
       write ( fileunit, * ) interdata%curves(ic)%uvbox(:,:,1)
       write ( fileunit, * ) interdata%curves(ic)%uvbox(:,:,2)
       write ( fileunit, * ) interdata%curves(ic)%nsplit
       do ip = 1,interdata%curves(ic)%nsplit
          write ( fileunit, * ) interdata%curves(ic)%isplit(:,ip)
       end do
       do is = 1,interdata%curves(ic)%nsplit-1
          write ( fileunit, * ) 1!class(is)
       end do
       if ( associated(interdata%curves(ic)%polyline) ) then
          write ( fileunit, * ) interdata%curves(ic)%polyline%np
          do ip = 1,interdata%curves(ic)%polyline%np
             write ( fileunit, * ) interdata%curves(ic)%polyline%uv(:,:,ip), interdata%curves(ic)%polyline%xyz(:,ip)
          end do
       else
          write ( fileunit, * ) 0
       end if
    end do
    close( fileunit )    

  end subroutine write_intersection_data









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






  subroutine read_curves(filename, surf, interdata)
    character(*),                 intent(in)    :: filename
    type(type_surface), target,   intent(in)    :: surf(:)
    type(type_intersection_data), intent(inout) :: interdata
    integer                                     :: fid, ncurves, np
    real(kind=fp), allocatable                  :: xyz(:,:), uv(:,:,:)
    type(ptr_surface)                           :: surfpair(2)
    integer                                     :: ic, ip, pair(2), endpt(2)

    call get_free_unit(fid)
    open(unit=fid, file=filename, action='read')
    read (fid,*) ncurves
    do ic = 1,ncurves
       read (fid,*) pair
       surfpair(1)%ptr => surf(pair(1))
       surfpair(2)%ptr => surf(pair(2))
       read (fid,*) np
       
       allocate(xyz(3,np), uv(2,2,np))
       do ip = 1,np
          read (fid,*) xyz(:,ip)
       end do
       do ip = 1,np
          read (fid,*) uv(:,1,ip), uv(:,2,ip)
       end do

       call add_intersection_point( &
            uv(:,:,1), &
            xyz(:,1), &
            surfpair, &
            2, &
            interdata, &
            endpt(1) ) 
       call add_intersection_point( &
            uv(:,:,np), &
            xyz(:,np), &
            surfpair, &
            2, &
            interdata, &
            endpt(2) )

       call add_intersection_curve( &
            interdata, &
            [0._fp, 0._fp, 0._fp], &
            endpt, &
            spread(spread([-1._fp, 1._fp], 2, 2), 3, 2) )
       interdata%curves(interdata%nc)%surf(1)%ptr => surf(pair(1))
       interdata%curves(interdata%nc)%surf(2)%ptr => surf(pair(2))
       interdata%curves(interdata%nc)%isplit(2,:) = [1,np]

       allocate(interdata%curves(interdata%nc)%polyline)
       interdata%curves(interdata%nc)%polyline%np = np
       call move_alloc(from=xyz, to=interdata%curves(interdata%nc)%polyline%xyz)
       call move_alloc(from=uv , to=interdata%curves(interdata%nc)%polyline%uv )
       
    end do
    close(fid)
    
  end subroutine read_curves

  
  
end program jouke
