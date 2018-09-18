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
       PARAM_hmin, &
       PARAM_hmax, &
       TOLchord, &
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
  
  
  IF ( .TRUE. ) THEN
     call write_tecplot_mesh( &
          mesh, &
          'Jouke/meshgen/brepmesh/brepmesh_optim_00.dat', &
          'pass_00' )
     call write_obj_mesh( &
          mesh, &
          'Jouke/meshgen/brepmesh/brepmesh_optim_00.obj' )
     call optim_jiao( &
          brep, &
          hyperedges(1:nhe), &
          nhe, &
          mesh, &
          20, &
          0.5_fp, &
          PARAM_hmin, &
          PARAM_hmax )
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

































  subroutine generate_brep_mesh( &
       brep, &
       hmin, &
       hmax, &
       tol, &
       feat_edge, &
       feat_vert, &
       hyperedges, &
       nhe, &
       mesh )
    use mod_brep
    use mod_hypergraph
    use mod_mesh
    use mod_util
    implicit none
    type(type_brep),         intent(in)    :: brep
    real(kind=fp),           intent(in)    :: hmin
    real(kind=fp),           intent(in)    :: hmax
    real(kind=fp),           intent(in)    :: tol
    logical,                 intent(in)    :: feat_edge(brep%ne)
    logical,                 intent(in)    :: feat_vert(brep%nv)
    integer,                 intent(in)    :: nhe
    type(type_hyperedge),    intent(inout) :: hyperedges(nhe)
    type(type_surface_mesh), intent(inout) :: mesh
    integer                                :: finf, fpts, fedg
    integer                                :: bv2mv(brep%nv), be2mv(2,brep%ne)
    integer, allocatable                   :: loc2glob(:), idsf(:), typf(:)
    real(kind=fp), allocatable             :: uvf(:,:,:), xyzf(:,:), uvi(:,:)
    integer, allocatable                   :: trif(:,:)
    integer                                :: npf, np, ntrif
    integer                                :: ifirstedge, ifirstpoint, sens, ifirst, ilast
    type(type_path)                        :: pathtmp
    integer                                :: iface, iwire, ihedg(2), ivert, ihype, iedge, i
    character(3)                           :: strnum
    type(type_surface_mesh)                :: meshvisu
    integer, allocatable                   :: t2bf(:)
    integer                                :: nt
    
    ! write 'info' file (common to all faces)
    call get_free_unit(finf)
    open(unit=finf, file='tmp/info.dat', action='write')
    write (finf,*) hmin
    write (finf,*) hmax
    write (finf,*) tol
    close(finf)

    allocate(loc2glob(100), idsf(100), typf(100))
    mesh%nv = 0
    mesh%nt = 0
    bv2mv(:) = 0
    be2mv(:,:) = 0
    NT = 0
    faces : do iface = 1,brep%nf
       ! write polynomial parametric surface 
       call write_polynomial(brep%faces(iface)%surface%x, 'tmp/c.cheb')
       !
       ! open boundary points & edges files
       call get_free_unit(fpts)
       open(unit=fpts, file='tmp/bpts.dat', action='write')
       call get_free_unit(fedg)
       open(unit=fedg, file='tmp/bedg.dat', action='write')
       !
       npf = 0
       wires : do iwire = 0,brep%faces(iface)%ninner ! <-----------------------------------+
          ! first halfedge of the wire                                                     !
          if ( iwire == 0 ) then ! <-------------------+                                   !
             ihedg = brep%faces(iface)%outer           !                                   !
          else ! --------------------------------------+                                   !
             ihedg = brep%faces(iface)%inner(:,iwire)  !                                   !
          end if ! <-----------------------------------+                                   !
          ifirstedge = ihedg(1)                                                            !
          !                                                                                !
          ifirstpoint = npf + 1                                                            !
          ! traverse the wire                                                              !
          halfedges : do ! <------------------------------------------------------------+  !
             ! get polyline points                                                      !  !
             call get_polyline_endpoints( &
                  brep, &
                  ihedg, &
                  ifirst, &
                  ilast, &
                  sens, &
                  np )
             !
             if ( .not.allocated(loc2glob) .or. &
                  size(loc2glob) < npf + np - 1 ) then ! <------+
                call reallocate_list(loc2glob, npf + np + 100)  !
             end if ! <-----------------------------------------+
             !
             if ( .not.allocated(idsf) .or. &
                  size(idsf) < npf + np - 1 ) then ! <----------+
                call reallocate_list(idsf, npf + np + 100)      !
             end if ! <-----------------------------------------+
             !
             if ( .not.allocated(typf) .or. &
                  size(typf) < npf + np - 1 ) then ! <----------+
                call reallocate_list(typf, npf + np + 100)      !
             end if ! <-----------------------------------------+
             !
             if ( allocated(uvf) .and. size(uvf,3) < np ) deallocate(uvf)
             if ( .not.allocated(uvf) ) allocate(uvf(2,2,np))
             if ( sens == 2 ) then
                uvf(1:2,1:2,1:np) = brep%edges(ihedg(1))%curve%polyline%uv&
                     (1:2,[2,1],ifirst:ilast)
             else
                uvf(1:2,1:2,1:np) = brep%edges(ihedg(1))%curve%polyline%uv&
                     (1:2,1:2,ifirst:ilast:-1)
             end if
             !
             if ( allocated(xyzf) .and. size(xyzf,2) < np ) deallocate(xyzf)
             if ( .not.allocated(xyzf) ) allocate(xyzf(3,np))
             xyzf(1:3,1:np) = brep%edges(ihedg(1))%curve%polyline%xyz&
                  (1:3,ifirst:ilast:(-1)**sens)
             !                                                                          !  !
             ! get brep vertex index of first point
             ivert = get_orig(brep, ihedg)
             if ( bv2mv(ivert) == 0 ) then ! <-----------------------------------+
                ! the brep vertex is not yet present in the mesh                 !
                bv2mv(ivert) = mesh%nv+1                                         !
                if ( feat_vert(ivert) ) then ! <------------------------------+  !
                   ! the mesh vertex is locked on a brep vertex               !  !
                   idsf(npf+1) = ivert                                        !  !
                   typf(npf+1) = 0                                            !  !
                else ! -------------------------------------------------------+  !
                   if ( feat_edge(ihedg(1)) ) then ! <---------------------+  !  !
                      ! the mesh vertex is locked on a hyperedge           !  !  !
                      idsf(npf+1) = ihedg(1)                               !  !  !
                      typf(npf+1) = 1                                      !  !  !
                   else ! -------------------------------------------------+  !  !
                      ! the mesh vertex is free to move along a hyperface  !  !  !
                      idsf(npf+1) = iface                                  !  !  !
                      typf(npf+1) = 2                                      !  !  !
                   end if ! <----------------------------------------------+  !  !
                end if ! <----------------------------------------------------+  !
                !                                                                !
                ! append new vertex to the mesh                                  !
                call append_vertices( &                                          !
                     mesh, &                                                     !
                     xyzf(1:3,1), &                                              !
                     uvf(1:2,1:2,1), &                                           !
                     idsf(npf+1), &                                              !
                     typf(npf+1), &                                              !
                     1 )                                                         !
             else
                if ( feat_edge(ihedg(1)) ) then
                   if ( mesh%typ(bv2mv(ivert)) == 2 ) then
                      mesh%ids(bv2mv(ivert)) = ihedg(1)
                      mesh%typ(bv2mv(ivert)) = 1
                      mesh%uv(1:2,1:2,bv2mv(ivert)) = uvf(1:2,1:2,1)
                   end if
                end if
             end if ! <----------------------------------------------------------+
             loc2glob(npf+1) = bv2mv(ivert)
             !
             if ( be2mv(1,ihedg(1)) == 0 ) then ! <-----------------------------------------+
                ! the brep edge is not yet present in the mesh                              !
                be2mv(1,ihedg(1)) = mesh%nv + 1                                             !
                be2mv(2,ihedg(1)) = mesh%nv + np - 2                                        !
                !                                                                           !
                if ( feat_edge(ihedg(1)) ) then ! <-----------------------+                 !
                   ! the mesh vertices are locked on a hyperedge          !                 !
                   idsf(npf+2:npf+np-1) = ihedg(1)                        !                 !
                   typf(npf+2:npf+np-1) = 1                               !                 !
                else ! ---------------------------------------------------+                 !
                   ! the mesh vertices are free to move along a hyperface !                 !
                   idsf(npf+2:npf+np-1) = iface                           !                 !
                   typf(npf+2:npf+np-1) = 2                               !                 !
                end if ! <------------------------------------------------+                 !
                !                                                                           !
                if ( sens == 1 ) then ! <-----------------------------------+               !
                   loc2glob(npf+2:npf+np-1) = mesh%nv + [(i, i=np-2,1,-1)]  !               !
                   call append_vertices( &                                  !               !
                     mesh, &                                                !               !
                     xyzf(1:3,np-1:2:-1), &                                 !               !
                     uvf(1:2,1:2,np-1:2:-1), &                              !               !
                     idsf(npf+2:npf+np-1), &                                !               !
                     typf(npf+2:npf+np-1), &                                !               !
                     np-2 )                                                 !               !
                else ! -----------------------------------------------------+               !
                   loc2glob(npf+2:npf+np-1) = mesh%nv + [(i, i=1,np-2)]     !               !
                   call append_vertices( &                                  !               !
                     mesh, &                                                !               !
                     xyzf(1:3,2:np-1), &                                    !               !
                     uvf(1:2,1:2,2:np-1), &                                 !               !
                     idsf(npf+2:npf+np-1), &                                !               !
                     typf(npf+2:npf+np-1), &                                !               !
                     np-2 )                                                 !               !
                end if ! <--------------------------------------------------+               !
             else ! ------------------------------------------------------------------------+
                if ( sens == 1 ) then ! <-----------------------------------+               !
                   loc2glob(npf+2:npf+np-1) = &                             !               !
                        [(i, i=be2mv(2,ihedg(1)),be2mv(1,ihedg(1)),-1)]     !               !
                else ! -----------------------------------------------------+               !
                   loc2glob(npf+2:npf+np-1) = &                             !               !
                        [(i, i=be2mv(1,ihedg(1)),be2mv(2,ihedg(1)))]        !               !
                end if ! <--------------------------------------------------+               !
             end if ! <---------------------------------------------------------------------+
             !
             ! write points...                            
             do i = 1,np-1 ! <---------------------------+
                write (fpts,*) uvf(:,1,i)             !
             end do ! <----------------------------------+
             ! ...and edges
             do i = 1,np-2 ! <---------------------------+
                write (fedg,*) npf + i, npf + i + 1      !
             end do ! <----------------------------------+
             !
             npf = npf + np - 1
             !
             ! move on to the next halfedge on the wire
             ihedg = get_next(brep, ihedg)
             !
             if ( ihedg(1) == ifirstedge ) then ! <------+
                ! the wire is complete                   !
                write (fedg,*) npf, ifirstpoint          !
                exit                                     !
             else ! -------------------------------------+
                write (fedg,*) npf, npf+1                !
             end if ! <----------------------------------+
             !           
          end do halfedges
       end do wires
       !
       ! close files
       close(fpts)
       close(fedg)
       !
       ! run mesher
       PRINT *,'MESH FACE #',IFACE
       !IF ( IFACE == 5 .OR. IFACE == 48 ) THEN
       !   PRINT *,'NPF =',NPF
       !   DO I = 1,NPF
       !      PRINT *,LOC2GLOB(I)
       !   END DO
       !END IF
       call system('/stck/bandrieu/Bureau/MeshGen/./meshgen.out &
            & tmp/c.cheb &
            & tmp/bpts.dat &
            & tmp/bedg.dat &
            & tmp/info.dat &
            & tmp/tri.dat &
            & tmp/uv.dat &
            & tmp/xyz.dat' )
       write (strnum,'(i3.3)') iface
       call system('cp tmp/c.cheb Jouke/meshgen/coeffs/c_'//strnum//'.cheb')
       call system('cp tmp/tri.dat Jouke/meshgen/brepmesh/uv/tri_'//strnum//'.dat')
       call system('cp tmp/uv.dat Jouke/meshgen/brepmesh/uv/uv_'//strnum//'.dat')
       !
       ! read face submesh                                                                     !
       call read_triangles( &                                                                  !
            'tmp/tri.dat', &                                                                   !
            trif, &                                                                            !
            ntrif )                                                                            !
       !
       call read_points( &                                                                     !
            'tmp/uv.dat', &                                                                    !
            uvi, &                                                                             !
            2, &                                                                               !
            np )                                                                               !
       if ( allocated(uvf) ) deallocate(uvf)
       allocate(uvf(2,2,np-npf))
       uvf(1:2,1,1:np-npf) = uvi(1:2,npf+1:np)
       !
       call read_points( &                                                                     !
            'tmp/xyz.dat', &                                                                   !
            xyzf, &                                                                            !
            3, &                                                                               !
            np )     
       !
       
       if ( size(loc2glob) < np ) call reallocate_list(loc2glob, np)
       loc2glob(npf+1:np) = mesh%nv + [(i, i=1,np-npf)]
       !                                                                                       !
       if ( size(idsf) < np ) call reallocate_list(idsf, np)                                   !
       idsf(npf+1:np) = iface                                                                  !
       !                                                                                       !
       if ( size(typf) < np ) call reallocate_list(typf, np)                                   !
       typf(npf+1:np) = 2                                                                      !
       !                                                                                       !
       ! ****************************
       call append_triangles( &
            meshvisu, &
            trif(1:3,1:ntrif) + meshvisu%nv, &
            ntrif )
       call append_vertices( &
            meshvisu, &
            xyzf(1:3,1:np), &
            spread(spread([0._fp,0._fp],dim=2,ncopies=2), dim=3, ncopies=np), &
            idsf(1:np), &
            typf(1:np), &
            np )
       ! ****************************
       
       do i = 1,3 ! <----------------------------------+                                       !
          trif(i,1:ntrif) = loc2glob(trif(i,1:ntrif))  !                                       !
       end do ! <--------------------------------------+                                       !
       ! add new triangles                                                                     !
       call append_triangles( &                                                                !
            mesh, &                                                                            !
            trif(1:3,1:ntrif), &                                                               !
            ntrif )                                                                            !
       !                                                                                       !
       ! add new mesh vertices                                                                 !
       call append_vertices( &                                                                 !
            mesh, &                                                                            !
            xyzf(1:3,npf+1:np), &                                                              !
            uvf(1:2,1:2,1:np-npf), &                                                           !
            idsf(npf+1:np), &                                                                  !
            typf(npf+1:np), &                                                                  !
            np-npf )
       !
       
       CALL INSERT_N_AFTER( &
            T2BF, &
            NT, &
            NTRIF, &
            SPREAD([IFACE], DIM=1, NCOPIES=NTRIF), &
            NT )
    end do faces

    ! create a path for each hyperedge
    do ihype = 1,nhe
       pathtmp%nv = 0
       pathtmp%hyperedge = ihype
       !
       if ( hyperedges(ihype)%verts(1) < 1 ) then ! <---+
          ihedg = hyperedges(ihype)%halfedges(:,1)      !
          ivert = get_orig(brep, ihedg)                 !
          hyperedges(ihype)%verts(1) = ivert            !
       end if ! <---------------------------------------+
       !
       do iedge = 1,hyperedges(ihype)%ne ! <---------------+
          ihedg = hyperedges(ihype)%halfedges(:,iedge)     !
          sens = 1 + mod(ihedg(2),2)                       !
          !                                                !
          call insert_after( &                             !
               pathtmp%verts, &                            !
               pathtmp%nv, &                               !
               bv2mv(get_orig(brep, ihedg)), &             !
               pathtmp%nv )                                !
          !                                                !
          np = be2mv(2,ihedg(1)) - be2mv(1,ihedg(1)) + 1   !
          if ( sens == 1 ) then ! <------+                 !
             ifirst = be2mv(2,ihedg(1))  !                 !
             ilast = be2mv(1,ihedg(1))   !                 !
          else ! ------------------------+                 !
             ifirst = be2mv(1,ihedg(1))  !                 !
             ilast = be2mv(2,ihedg(1))   !                 !
          end if ! <---------------------+                 !
          !                                                !
          call insert_n_after( &                           !
               pathtmp%verts, &                            !
               pathtmp%nv, &                               !
               np, &                                       !
               [(i, i=ifirst,ilast,(-1)**sens)], &         !
               pathtmp%nv )                                !
       end do ! <------------------------------------------+
       !
       if ( hyperedges(ihype)%verts(2) > 0 ) then ! <---+
          call insert_after( &                          !
               pathtmp%verts, &                         !
               pathtmp%nv, &                            !
               bv2mv(hyperedges(ihype)%verts(2)), &     !
               pathtmp%nv )                             !
       end if ! <---------------------------------------+
       !
       call append_path( &
            mesh, &
            pathtmp )
       !
       ! compute discrete curvilinear abscissa along the path
       allocate(mesh%paths(mesh%npaths)%s(mesh%paths(mesh%npaths)%nv))
       mesh%paths(mesh%npaths)%s(1) = 0._fp
       do i = 2,mesh%paths(mesh%npaths)%nv ! <-------------------------------+
          mesh%paths(mesh%npaths)%s(i) = mesh%paths(mesh%npaths)%s(i-1) + &  !
               norm2( &                                                      !
               mesh%xyz(:,mesh%paths(mesh%npaths)%verts(i)) - &              !
               mesh%xyz(:,mesh%paths(mesh%npaths)%verts(i-1)) )              !
       end do ! <------------------------------------------------------------+
       !
    end do


    call write_tecplot_mesh_iface( &
         meshvisu, &
         'Jouke/meshgen/brepmesh/brepmesh.dat', &
         'mesh_brep', &
         t2bf )
   
    
  end subroutine generate_brep_mesh








  subroutine read_triangles( &
       filename, &
       tri, &
       n )
    use mod_util
    implicit none
    character(*),         intent(in)    :: filename
    integer, allocatable, intent(inout) :: tri(:,:)
    integer,              intent(out)   :: n
    integer, allocatable                :: tmp(:,:)
    integer                             :: fid, io, t(3)

    call get_free_unit(fid)
    open(unit=fid, file=filename, action='read')
    n = 0
    do
       read (fid, *, iostat=io) t
       if ( io /= 0 ) exit
       if ( .not.allocated(tri) ) then
          allocate(tri(3,100))
       else
          if ( n + 1 > size(tri,2) ) then
             call move_alloc(from=tri, to=tmp)
             allocate(tri(3,n+100))
             tri(1:3,1:n) = tmp(1:3,1:n)
             deallocate(tmp)
          end if
       end if
       n = n + 1
       tri(1:3,n) = t(1:3)
    end do
    close(fid)

  end subroutine read_triangles




  subroutine read_points( &
       filename, &
       pts, &
       m, &
       n )
    use mod_util
    implicit none
    character(*),               intent(in)    :: filename
    real(kind=fp), allocatable, intent(inout) :: pts(:,:)
    integer,                    intent(in)    :: m
    integer,                    intent(out)   :: n
    real(kind=fp), allocatable                :: tmp(:,:)
    real(kind=fp)                             :: p(m)
    integer                                   :: fid, io

    call get_free_unit(fid)
    open(unit=fid, file=filename, action='read')
    n = 0
    do
       read (fid, *, iostat=io) p
       if ( io /= 0 ) exit
       if ( .not.allocated(pts) ) then
          allocate(pts(m,100))
       else
          if ( n + 1 > size(pts,2) ) then
             call move_alloc(from=pts, to=tmp)
             allocate(pts(m,n+100))
             pts(1:m,1:n) = tmp(1:m,1:n)
             deallocate(tmp)
          end if
       end if
       n = n + 1
       pts(1:m,n) = p(1:m)
    end do
    close(fid)

  end subroutine read_points






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
       write (fid,'(3I7)') mesh%tri(:,j)
    end do

    close(fid)
    
  end subroutine write_tecplot_mesh_iface


end program jouke
