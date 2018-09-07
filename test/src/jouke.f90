program jouke

  use mod_util
  use mod_math
  use mod_polynomial
  use mod_diffgeom
  use mod_types_intersection
  use mod_intersection
  use mod_brep2
  use mod_hypergraph
  
  implicit none
  
  character(3)               :: strnum3
  integer                    :: fid
  integer*8                  :: tic, toc, count_rate

  logical, parameter                      :: ECONOMIZE = .true.
  integer                                 :: nsurf
  type(type_surface), target, allocatable :: surf(:)
  logical, allocatable                    :: mask(:)
  type(type_intersection_data)            :: interdata
  type(type_BREP)                         :: brep
  integer, allocatable                    :: surf2curv(:,:), faces(:), valence(:)
  logical, allocatable                    :: feat_edge(:), feat_vert(:)
  type(type_hyperface), allocatable       :: hyperfaces(:)
  type(type_hyperedge), allocatable       :: hyperedges(:)
  integer                                 :: nsurf2curv, nfaces, nhf, nhe
  integer                                 :: isurf, ic, jsurf, ksurf, iface, I, J, k

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
     !MASK(1:36) = .FALSE.
     !MASK(38:100) = .FALSE.
     !MASK(112:) = .FALSE.
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



  PRINT *,'TOTAL :'
  PRINT *,INTERDATA%NP,' INTERSECTION POINT(S)'
  PRINT *,INTERDATA%NC,' INTERSECTION CURVE(S)'


  !allocate(brep%dcel%faces(1:200), brep%dcel%vertices(1:200), brep%dcel%halfedges(1:700), brep%faces(1:200), brep%edges(1:350))
  !call brep_from_intersection_data( &
  !     surf, &
  !     nsurf, &
  !     interdata, &
  !     brep )
  call make_brep_from_intersection_data( &
       surf, &
       nsurf, &
       interdata, &
       brep )

  PRINT *,'BREP : NV =',BREP%NV,', NE =',BREP%NE,', NF =',BREP%NF
  PRINT *,'EULER-POINCARE =',BREP%NV - BREP%NE + BREP%NF

  call system_clock( toc )
  PRINT *,''; PRINT *,''; PRINT *,''
  PRINT *,'ELAPSED =',REAL( TOC - TIC ) / REAL( COUNT_RATE )


  call write_brep_files( &
       brep, &
       'Jouke/brep/verts.dat', &
       'Jouke/brep/edges.dat', &
       'Jouke/brep/faces.dat' )

  IF ( .FALSE. ) THEN ! HYPERGRAPHE
     ! get feature edges
     allocate(feat_edge(brep%ne))
     call get_feature_edges(brep, feat_edge)
     call get_free_unit( fid )
     open(unit=fid, file='Jouke/brep/feat_edge.dat', action='write')
     do i = 1,brep%ne
        if ( feat_edge(i) ) write (fid,*) i
     end do
     close(fid)

     ! get hyperfaces
     call get_hyperfaces(brep, feat_edge, hyperfaces, nhf)
     call get_free_unit( fid )
     open(unit=fid, file='Jouke/brep/hyperfaces.dat', action='write')
     write (fid,*) nhf
     do i = 1,nhf
        write (fid,*) hyperfaces(i)%nf
        write (fid,*) hyperfaces(i)%faces(1:hyperfaces(i)%nf)
     end do
     close(fid)

     ! get feature vertices
     allocate(feat_vert(brep%nv), valence(brep%nv))
     call get_feature_vertices( &
          brep, &
          feat_edge, &
          feat_vert, &
          valence )
     call get_free_unit( fid )
     open(unit=fid, file='Jouke/brep/feat_vert.dat', action='write')
     do i = 1,brep%nv
        if ( feat_vert(i) ) then
           write (fid,'(I1,1x)',advance='no') 1
        else
           write (fid,'(I1,1x)',advance='no') 0
        end if
        write (fid,*) valence(i)
     end do
     close(fid)


     ! get hyperedges
     call get_hyperedges( &
          brep, &
          feat_edge, &
          feat_vert, &
          valence, &
          hyperedges, &
          nhe )
     call get_free_unit( fid )
     open(unit=fid, file='Jouke/brep/hyperedges.dat', action='write')
     write (fid,*) nhe
     do i = 1,nhe
        write (fid,*) hyperedges(i)%ne
        write (fid,*) hyperedges(i)%verts
        do j = 1,hyperedges(i)%ne
           write (fid,*) hyperedges(i)%edges(1:2,j)
        end do
     end do
     close(fid)
  END IF

  

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

  call write_intersection_data( &
       interdata, &
       'Jouke/result/interdata_points.dat', &
       'Jouke/result/interdata_curves.dat' )



  

  open(unit=13, file='Jouke/brep/v2f.dat', action='write')
  do i = 1,brep%nv
     call get_v2f( &
          brep, &
          i, &
          faces, &
          nfaces )
     write (13,*) nfaces
     write (13,*) faces(1:nfaces)
  end do
  close(13)

  
  IF ( .true. ) THEN
     do iface = 1,brep%nf
        PRINT *,'meshgen face #',iface
        write (strnum3,'(I3.3)') iface
        call generate_face_mesh( &       
             brep, &
             iface, &
             'Jouke/meshgen/coeffs/c_'//strnum3//'.cheb', &
             'Jouke/meshgen/contours/uv_'//strnum3//'.dat', &
             'Jouke/meshgen/contours/edges_'//strnum3//'.dat', &
             'Jouke/meshgen/info.dat', &
             'Jouke/meshgen/tri_'//strnum3//'.dat', &
             'Jouke/meshgen/uv_'//strnum3//'.dat' )
     end do
  END IF

  
  IF ( .false.)  THEN
     DO I = 1,BREP%NE
        IF ( IS_BOUNDARY_EDGE(BREP, I) ) PRINT *,'EDGE #',I,' ON BOUNDARY'
     END DO
     
     DO I = 1,BREP%NE-1
        DO J = I+1,BREP%NE
           IF ( ASSOCIATED(BREP%EDGES(I)%CURVE, BREP%EDGES(J)%CURVE) .AND. &
                BREP%EDGES(I)%ISPLIT == BREP%EDGES(J)%ISPLIT ) THEN
              PRINT *,'EDGES',I,J,' -> SAME CURVE SEGMENT'
           END IF
        END DO
     END DO

     k = 0
     open(unit=13, file='Jouke/debug_curve_segments.dat', action='write')
     do ic = 1,interdata%nc
        do j = 1,interdata%curves(ic)%nsplit-1
           if ( interdata%curves(ic)%iedge(j) > 0 ) then
              k = k + 1
              write (13,*) interdata%curves(ic)%isplit(2,j+1) - interdata%curves(ic)%isplit(2,j) + 1
              do i = interdata%curves(ic)%isplit(2,j),interdata%curves(ic)%isplit(2,j+1)
                 write (13,*) interdata%curves(ic)%polyline%xyz(:,i)
              end do
           end if
        end do
     end do
     close(13)
     PRINT *,'K =',K
  END IF

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
       write ( fileunit, * ) logic2int(interdata%curves(ic)%smooth)
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
