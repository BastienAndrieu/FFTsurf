program jouke

  use mod_util
  use mod_math
  use mod_polynomial
  use mod_diffgeom
  use mod_types_intersection
  use mod_intersection
  use mod_types_brep
  use mod_brep
  use mod_hypergraph
  use mod_mesh
  use mod_tolerances
  use mod_optimmesh

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
  type(type_surface_mesh)                 :: mesh
  real(kind=fp)                           :: xyzverif(3,2)
  real(kind=fp), allocatable              :: wei(:), hve(:), ener(:), grad(:,:), hess(:,:,:)
  real(kind=fp)                           :: hminsqr
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

  !call read_curves('Jouke/propergol/tangent_curves.dat', surf, interdata)
  call read_intersection_curves( &
     'Jouke/propergol/tangent_curves.dat', &
     surf, &
     interdata )
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
       mask, &
       5.d-3, &
       1.d-3, &
       1.d-2 )
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

  IF ( .true. ) THEN ! HYPERGRAPHE
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
        call reverse_hyperedge(hyperedges(i))
        write (fid,*) hyperedges(i)%ne
        write (fid,*) hyperedges(i)%verts
        write (fid,*) hyperedges(i)%hyperfaces
        do j = 1,hyperedges(i)%ne
           !PRINT *,get_polyline_length(brep, hyperedges(i)%halfedges(1:2,j))
           write (fid,*) hyperedges(i)%halfedges(1:2,j)
        end do
     end do
     close(fid)
  END IF



  
  !PRINT *,'EDGE -> HYPEREDGE'
  !DO I = 1,BREP%NE
  !   IF ( FEAT_EDGE(I) ) PRINT *, I, BREP%EDGES(I)%HYPEREDGE
  !END DO
  !PAUSE
  

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



  PRINT *,'GENERATE BREP MESH...'
  call generate_brep_mesh( &
       brep, &
       1.d-3, &!PARAM_hmin, &
       1.d-2, &!PARAM_hmax, &
       5.d-3, &!TOLchord, &
       feat_edge, &
       feat_vert, &
       hyperedges(1:nhe), &
       nhe, &
       mesh )
  PRINT *,'...OK'

  PRINT *,'MAKE HALFEDGES...'
  call make_halfedges( &
       mesh )
  PRINT *,'...OK'
  
  open(unit=13, file='Jouke/meshgen/brepmesh/mv2h.dat', action='write')
  do i = 1,mesh%nv
     write (13,*) mesh%v2h(:,i)
  end do
  close(13)
  open(unit=13, file='Jouke/meshgen/brepmesh/mtwin.dat', action='write')
  do i = 1,mesh%nt
     write (13,*) mesh%twin(:,:,i)
  end do
  close(13)

  call write_inria_mesh( &
       mesh, &
       'Jouke/meshgen/brepmesh/brepmesh.mesh' )
  
  call write_mesh_files( &
       mesh, &
       'Jouke/meshgen/brepmesh/tri.dat', &
       'Jouke/meshgen/brepmesh/xyz.dat', &
       'Jouke/meshgen/brepmesh/uv.dat', &
       'Jouke/meshgen/brepmesh/idstyp.dat', &
       'Jouke/meshgen/brepmesh/paths.dat' )

  
  ! CHECK UVs
  IF ( .true. ) THEN
     PRINT *,'CHECK UVs'
     do i = 1,mesh%nv
        select case(mesh%typ(i))
        case (0)
        case (1)
           do j = 1,2
              iface = brep%edges(mesh%ids(i))%halfedges(1+mod(j,2))%face
              call eval( &
                   xyzverif(:,j), &
                   brep%faces(iface)%surface, &
                   mesh%uv(:,j,i) )
           end do
           !PRINT *,'I =', I,', ERR =', norm2(mesh%xyz(:,i) - xyzverif(:,1)), norm2(mesh%xyz(:,i) - xyzverif(:,2))
        case (2)
           call eval( &
                xyzverif(:,1), &
                brep%faces(mesh%ids(i))%surface, &
                mesh%uv(:,1,i) )
           IF ( norm2(mesh%xyz(:,i) - xyzverif(:,1)) > 1.D-13 ) THEN
              PRINT *,'I =', I,', ERR =', norm2(mesh%xyz(:,i) - xyzverif(:,1))
           END IF
        end select
     end do
  END IF


  ! CONTRACT SMALL EDGES
  IF ( .false. ) THEN
     STOP
  ELSE
     call contract_small_edges( &
          brep, &
          hyperedges(1:nhe), &
          nhe, &
          mesh, &
          1.d-3 )!PARAM_hmin )
     PRINT *,'(RE)MAKE HALFEDGES...'
     deallocate(mesh%v2h, mesh%twin)
     call make_halfedges( &
          mesh )
     PRINT *,'...OK'
  END IF

  call write_inria_mesh( &
       mesh, &
       'Jouke/meshgen/brepmesh/brepmesh_ec.mesh' )

  call write_mesh_files( &
       mesh, &
       'Jouke/meshgen/brepmesh/tri_ec.dat', &
       'Jouke/meshgen/brepmesh/xyz_ec.dat', &
       'Jouke/meshgen/brepmesh/uv_ec.dat', &
       'Jouke/meshgen/brepmesh/idstyp_ec.dat', &
       'Jouke/meshgen/brepmesh/paths_ec.dat' )
  
  open(unit=13, file='Jouke/meshgen/brepmesh/mv2h_ec.dat', action='write')
  do i = 1,mesh%nv
     write (13,*) mesh%v2h(:,i)
  end do
  close(13)
  open(unit=13, file='Jouke/meshgen/brepmesh/mtwin_ec.dat', action='write')
  do i = 1,mesh%nt
     write (13,*) mesh%twin(:,:,i)
  end do
  close(13)

  !STOP
  
  
  IF ( .TRUE. ) THEN
     call write_tecplot_mesh( &
          mesh, &
          'Jouke/meshgen/brepmesh/brepmesh_optim_00.dat', &
          'pass_00' )
     !call write_obj_mesh( &
     !     mesh, &
     !     'Jouke/meshgen/brepmesh/brepmesh_optim_00.obj' )
     call optim_jiao( &
          brep, &
          hyperedges(1:nhe), &
          nhe, &
          mesh, &
          1.0_fp, &
          0.75_fp, &
          5, &
          15, &
          100, &
          1.d-3, &!PARAM_hmin, &
          1.d-2 )!PARAM_hmax )
     call write_inria_mesh( &
          mesh, &
          'Jouke/meshgen/brepmesh/brepmesh_optim.mesh' )
     call write_mesh_files( &
       mesh, &
       'Jouke/meshgen/brepmesh/tri.dat', &
       'Jouke/meshgen/brepmesh/xyz_smooth.dat', &
       'Jouke/meshgen/brepmesh/uv_smooth.dat', &
       'Jouke/meshgen/brepmesh/idstyp_smooth.dat', &
       'Jouke/meshgen/brepmesh/paths.dat' )
  ELSE
     allocate(wei(mesh%nt), hve(mesh%nv), ener(mesh%nt), grad(3,mesh%nv), hess(3,3,mesh%nv))
     wei(:) = 1._fp
     hve(:) = 1._fp
     call energy_jiao( &
          mesh%nt, &
          mesh%tri(1:3,1:mesh%nt), &
          wei, &
          mesh%nv, &
          mesh%xyz(1:3,1:mesh%nv), &
          hve, &
          0.7_fp, &
          hminsqr, &
          ener, &
          grad, &
          hess )

     open(unit=13, file='Jouke/meshgen/brepmesh/ener.dat', action='write')
     do i = 1,mesh%nt
        write (13,*) ener(i)
     end do
     close(13)

     open(unit=13, file='Jouke/meshgen/brepmesh/grad.dat', action='write')
     do i = 1,mesh%nv
        write (13,*) grad(:,i)
     end do
     close(13)

     open(unit=13, file='Jouke/meshgen/brepmesh/hess.dat', action='write')
     do i = 1,mesh%nv
        write (13,*) hess(:,:,i)
     end do
     close(13)
  END IF

! CHECK UVs
  IF ( .true. ) THEN
     PRINT *,'CHECK UVs'
     do i = 1,mesh%nv
        select case(mesh%typ(i))
        case (0)
        case (1)
           do j = 1,2
              iface = brep%edges(mesh%ids(i))%halfedges(1+mod(j,2))%face
              call eval( &
                   xyzverif(:,j), &
                   brep%faces(iface)%surface, &
                   mesh%uv(:,j,i) )
           end do
           IF ( max(norm2(mesh%xyz(:,i) - xyzverif(:,1)), norm2(mesh%xyz(:,i) - xyzverif(:,2))) > 1.D-13 ) THEN
              PRINT *,'I =', I,', ERR =', norm2(mesh%xyz(:,i) - xyzverif(:,1)), norm2(mesh%xyz(:,i) - xyzverif(:,2))
           END IF
        case (2)
           call eval( &
                xyzverif(:,1), &
                brep%faces(mesh%ids(i))%surface, &
                mesh%uv(:,1,i) )
           IF ( norm2(mesh%xyz(:,i) - xyzverif(:,1)) > 1.D-13 ) THEN
              PRINT *,'I =', I,', ERR =', norm2(mesh%xyz(:,i) - xyzverif(:,1))
           END IF
        end select
     end do
  END IF
  
  
  STOP

  IF ( .true. ) THEN
     do iface = 1,brep%nf
        PRINT '(A14,I0,A1,I0)','meshgen face #',iface,'/',brep%nf
        write (strnum3,'(I3.3)') iface
        call generate_face_mesh( &       
             brep, &
             iface, &
             'Jouke/meshgen/coeffs/c_'//strnum3//'.cheb', &
             'Jouke/meshgen/contours/uv_'//strnum3//'.dat', &
             'Jouke/meshgen/contours/edges_'//strnum3//'.dat', &
             'Jouke/meshgen/info.dat', &
             'Jouke/meshgen/tri_'//strnum3//'.dat', &
             'Jouke/meshgen/uv_'//strnum3//'.dat', &
             'Jouke/meshgen/xyz_'//strnum3//'.dat' )
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










































  subroutine write_tecplot_mesh_iface( &
       mesh, &
       filename, &
       zonename, &
       iface )
    use mod_util
    implicit none
    type(type_surface_mesh), intent(in) :: mesh
    character(*),            intent(in) :: filename
    character(*),            intent(in) :: zonename
    integer,                 intent(in) :: iface(mesh%nt)
    integer                             :: fid, i, j

    call get_free_unit(fid)
    open(unit=fid, file=filename, action='write')

    write (fid,*) 'VARIABLES = "X" "Y" "Z" "IFACE"'

    write (fid,*) 'ZONE T="' // zonename // '"'
    write (fid,'(A2,I0,A3,I0)') 'N=',mesh%nv,' E=',mesh%nt
    write (fid,*) 'ZONETYPE=FETriangle'
    write (fid,*) 'DATAPACKING=BLOCK'

    write (fid,*) 'VARLOCATION = ([1-3]=NODAL, [4]=CELLCENTERED)'

    do i = 1,3
       write (fid,*) ''
       do j = 1,mesh%nv
          write (fid,'(ES22.15)') mesh%xyz(i,j)
       end do
    end do

    write (fid,*) ''

    do j = 1,mesh%nt
       write (fid,'(I0)') iface(j)
    end do

    write (fid,*) ''

    do j = 1,mesh%nt
       write (fid,'(I0,1x,I0,1x,I0)') mesh%tri(:,j)
    end do

    close(fid)
    
  end subroutine write_tecplot_mesh_iface


end program jouke
