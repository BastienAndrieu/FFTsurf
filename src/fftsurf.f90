program fftsurf

  use mod_util
  use mod_options
  use mod_math
  use mod_types_intersection
  use mod_types_brep
  use mod_brep
  use mod_hypergraph
  use mod_mesh
  use mod_optimmesh
  use mod_tolerances
  use mod_propagation
  use mod_polynomial
  use mod_intersection
  use mod_init

  implicit none

  ! --------------------------------------------------------------
  ! Parameters
  integer, parameter                    :: freq_checkpoint = 10
  integer, parameter                    :: PARAM_passmax = 30!
  real(kind=fp), parameter              :: PARAM_frac_conf1 = 1._fp
  real(kind=fp), parameter              :: PARAM_frac_conf2 = 0.7_fp
  integer, parameter                    :: PARAM_ipass1 = 0!5
  integer, parameter                    :: PARAM_ipass2 = 0!15
  ! --------------------------------------------------------------
  ! Options
  type(type_options)                    :: options
  integer                               :: narg, iarg, fid, irep, io
  character(99)                         :: fileoptions, arg
  ! --------------------------------------------------------------
  ! BREP
  type(type_surface), allocatable       :: surf(:)
  integer                               :: nsurf, isurf
  type(type_intersection_data), pointer :: interdata_old => null(), interdata_new => null()
  type(type_brep),              pointer :: brep_old => null(), brep_new => null()
  ! --------------------------------------------------------------
  ! Hypergraph
  type(type_hypergraph), pointer        :: hypergraph_old => null(), hypergraph_new => null()
  logical, allocatable                  :: feat_edge(:)
  logical, allocatable                  :: feat_vert(:)
  ! --------------------------------------------------------------
  ! Surface mesh
  type(type_surface_mesh)               :: mesh
  real(kind=fp), allocatable            :: xyzprev(:,:), xyztmp(:,:), dxyz(:,:)
  integer                               :: stat_corresp, stat_regen_path
  integer                               :: ivert, iface
  REAL(KIND=FP), ALLOCATABLE            :: htarget(:),tri_weights(:)
  ! --------------------------------------------------------------
  ! Physical time
  real                                  :: time
  integer                               :: instant
  ! --------------------------------------------------------------
  ! CPU time
  real                                  :: tic, toc
  ! --------------------------------------------------------------
  integer :: i, j
  character(3)                          :: strnum, strnum3
  integer :: nf_old

  !! Read argument (options file name)
  narg = command_argument_count()
  if ( narg < 1 ) then
     STOP 'fftsurf: user must provide an option file'
  else
     call get_command_argument(1, fileoptions)
  end if


  !! Parse additional options (override)
  iarg = 2
  do while(iarg <= narg)
     call get_command_argument(iarg, arg)
     print *,'arg =',arg
     select case(arg)
     case ('-rep')
        iarg = iarg + 1
        options%reprise = .true.
     case ('-msh')
        iarg = iarg + 1
        options%from_msh = .true.
     end select
  end do



  !! Initialization
  PRINT *,'INITIALIZATION...'
  call cpu_time(tic)
  allocate(interdata_old, brep_old, hypergraph_old)
  if ( options%from_msh ) then
     call init_from_mesh( &
          fileoptions, &
          options, &
          surf, &
          nsurf, &
          interdata_old, &
          brep_old, &
          hypergraph_old, &
          mesh )
  else
     call init_from_surfaces( &
          fileoptions, &
          options, &
          surf, &
          nsurf, &
          interdata_old, &
          brep_old, &
          hypergraph_old, &
          mesh )
  end if
  call cpu_time(toc)
  PRINT '(a8,1x,f8.3,1x,a1)','ELAPSED:', toc - tic, 's'



  if ( options%reprise ) then
     call get_free_unit(fid)
     open(unit=fid, file=trim(options%directory) // 'checkpoint/time.dat', action='read')
     read (fid,*) instant
     read (fid,*) time
     close(fid)

     open(unit=fid, file=trim(options%directory) // 'output/reprises.dat', action='readwrite')
     irep = 0
     do
        read (fid,*,iostat=io) iarg
        if ( io /= 0 ) exit
        irep = irep + 1
     end do
     backspace(fid)
     write (fid,*) instant
     close(fid)
     irep = irep + 1
  else
     time = 0.0
     instant = 0
     irep = 1
     call get_free_unit(fid)
     open(unit=fid, file=trim(options%directory) // 'output/reprises.dat', action='write')
     write (fid,*) instant
     close(fid)
  end if


  if ( options%mode > 1 ) then
     ! Export initial connectivity, face refs and vertex positions
     call write_connectivity( &
          trim(options%directory) // 'output/', &
          mesh, &
          irep )
     call write_face_ref( &
          trim(options%directory) // 'output/', &
          mesh, &
          irep )
     call write_xyz_positions( &
          trim(options%directory) // 'output/', &
          mesh, &
          instant )

     allocate(xyzprev(3,mesh%nv), xyztmp(3,mesh%nv), dxyz(3,mesh%nv))

     call write_vtk_mesh( &
          mesh, &
          '../debug/test.vtk' )
     call write_gmsh_mesh( &
          mesh, &
          '../debug/test.msh' )
     PAUSE
  end if


  !! Main loop
  PRINT *,'MAIN LOOP...'
  call cpu_time(tic)
  if ( options%mode > 0 ) allocate(interdata_new, brep_new)
  if ( options%mode > 1 ) allocate(hypergraph_new)
  if ( options%mode > 1 ) allocate(htarget(mesh%nv), tri_weights(mesh%nt))

  if ( options%mode > 1 ) nf_old = brep_old%nf
  main_loop : do
     ! checkpoint
     if ( mod(instant,freq_checkpoint) == 0 ) then
        call make_checkpoint( &
             options, &
             surf(1:nsurf), &
             nsurf, &
             interdata_old, &
             instant, &
             time )
     end if

     ! Propagate surfaces and advance one timestep
     do isurf = 1,nsurf
        if ( surf(isurf)%tag /= 1 ) cycle
        call propagation_step_RK4( &
             surf(isurf), &
             options%propagation_law, &
             options%timestep, &
             time )
     end do
     time = time + options%timestep
     instant = instant + 1

     PRINT *,''
     PRINT '(a7,1x,i0)','INSTANT',INSTANT
     PRINT *,'TIME:', time, '/', options%timespan



     !***************
     if ( .TRUE. ) then!options%mode == 0 ) then
        write (strnum, '(i3.3)') instant
        do isurf = 1,nsurf
           write (strnum3, '(i3.3)') isurf
           call write_polynomial( &
                surf(isurf)%x, &
                trim(options%directory) // 'output/debug/instant_' // strnum // 'c_' // strnum3 // '.cheb' )
        end do
     end if
     !***************


     if ( options%mode > 0 ) then
        ! reset brep data
        call free_brep(brep_new)
        !call free_intersection_data(interdata_new)
        nullify(interdata_new)
        allocate(interdata_new)

        call copy_tangent_intersections( &
             from = interdata_old, &
             to   = interdata_new )

        call update_intersection_curves(interdata_new)

        call make_brep( &
             surf, &
             nsurf, &
             interdata_new, &
             brep_new, &
             options%chord_err, &
             options%hmin, &
             options%hmax )
        print *,'nf_old =',nf_old

        ! debugging >> ..................
        call write_brep_files( &
             brep_new, &
             trim(options%directory) // 'output/debug/brep/verts.dat', &
             trim(options%directory) // 'output/debug/brep/edges.dat', &
             trim(options%directory) // 'output/debug/brep/faces.dat' )
        ! .................. <<


        if ( options%mode > 1 ) then
           call free_hypergraph(hypergraph_new)
           call make_hypergraph( &
                brep_new, &
                hypergraph_new, &
                feat_edge, &
                feat_vert )
           ! debugging >> ..................
           call get_free_unit( fid )
           open(unit=fid, file=trim(options%directory) // 'output/debug/brep/hyperfaces.dat', action='write')
           write (fid,*) hypergraph_new%nhf
           do i = 1,hypergraph_new%nhf
              write (fid,*) hypergraph_new%hyperfaces(i)%nf
              write (fid,*) hypergraph_new%hyperfaces(i)%faces(1:hypergraph_new%hyperfaces(i)%nf)
           end do
           close(fid)

           open(unit=fid, file=trim(options%directory) // 'output/debug/brep/hyperedges.dat', action='write')
           write (fid,*) hypergraph_new%nhe
           do i = 1,hypergraph_new%nhe
              write (fid,*) hypergraph_new%hyperedges(i)%ne
              write (fid,*) hypergraph_new%hyperedges(i)%verts
              write (fid,*) hypergraph_new%hyperedges(i)%hyperfaces
              do j = 1,hypergraph_new%hyperedges(i)%ne
                 write (fid,*) hypergraph_new%hyperedges(i)%halfedges(1:2,j)
              end do
           end do
           close(fid)
           ! .................. <<

           if ( brep_new%nf /= nf_old ) then
              stat_corresp = 1
           else
              call update_mesh_correspondance( &
                   brep_old, &
                   brep_new, &
                   hypergraph_old, &
                   hypergraph_new, &
                   mesh, &
                   stat_corresp )
           end if

           if ( stat_corresp > 0 ) then
              call make_checkpoint( &
                   options, &
                   surf(1:nsurf), &
                   nsurf, &
                   interdata_new, &
                   instant, &
                   time )
              PRINT *,'failed to update mesh correspondance, made a checkpoint'
              STOP
           end if

           nf_old = brep_new%nf
           !PAUSE

           !xyzprev(1:3,1:mesh%nv) = mesh%xyz(1:3,1:mesh%nv)
           call write_tecplot_mesh( &
                mesh, &
                '../debug/pre_deform0.dat', &
                'avant' )

           ! Regenerate features mesh
           ! uv, ids
           IF ( .false.) THEN
              call regenerate_feature_paths( &
                   brep_new, &
                   hypergraph_new, &
                   mesh, &
                   stat_regen_path )
           ELSE
              call new_regen_paths( &
                   brep_new, &
                   hypergraph_new, &
                   mesh, &
                   stat_regen_path )
              !PAUSE
           END IF
           if ( stat_regen_path > 0 ) then
              PRINT *,'failed to regenerate feature paths'
              PAUSE
           end if

           ! xyz
           dxyz(1:3,1:mesh%nv) = 0._fp
           do ivert = 1,mesh%nv ! <-------------------------------------------------+
              select case ( mesh%typ(ivert) ) ! <-------------------------------+   !
              case (0) ! -------------------------------------------------------+   !
                 !mesh%xyz(:,ivert) = brep_new%verts(mesh%ids(ivert))%point%xyz  !   !
                 dxyz(1:3,ivert) = brep_new%verts(mesh%ids(ivert))%point%xyz - &
                      mesh%xyz(1:3,ivert)
              case (1) ! -------------------------------------------------------+   !
                 iface = brep_new%edges(mesh%ids(ivert))%halfedges(2)%face      !   !
                 call eval( &                                                   !   !
                      dxyz(1:3,ivert), &                                        !   !
                      brep_new%faces(iface)%surface, &                          !   !
                      mesh%uv(:,1,ivert) )                                      !   !
                 dxyz(1:3,ivert) = dxyz(1:3,ivert) - mesh%xyz(1:3,ivert)        !   !
              end select ! <----------------------------------------------------+   !
           end do ! <---------------------------------------------------------------+

           ! (Pre-deform mesh to prevent inverted elements)
           call spring_displacement_smoothing( &
                mesh, &
                dxyz, &
                20 )
           xyztmp = mesh%xyz(1:3,1:mesh%nv) + dxyz
           xyzprev = mesh%xyz(1:3,1:mesh%nv) !*
           mesh%xyz(1:3,1:mesh%nv) = xyztmp!*mesh%xyz(1:3,1:mesh%nv) + dxyz !*
           ! ********
           !call write_xyz_positions( &
           !     '../debug/', &
           !     mesh, &
           !     0 )
           call write_tecplot_mesh( &
                mesh, &
                '../debug/pre_deform1.dat', &
                'smooth displacements' )
           !mesh%xyz(1:3,1:mesh%nv) = xyzprev
           ! ********

           !PAUSE
           do ivert = 1,mesh%nv ! <------------------------------------------+
              if ( mesh%typ(ivert) < 2 ) cycle                               !
              iface = mesh%ids(ivert)                                        !
              call eval( &                                                   !
                   mesh%xyz(:,ivert), &                                      !
                   brep_new%faces(iface)%surface, &                          !
                   mesh%uv(:,1,ivert) )                                      !
           end do ! <--------------------------------------------------------+

           call write_tecplot_mesh( &
                mesh, &
                '../debug/pre_deform2.dat', &
                'update f(uv)' )

           IF ( .true. ) THEN
              !dxyz(1:3,1:mesh%nv) = xyzprev(1:3,1:mesh%nv) - mesh%xyz(1:3,1:mesh%nv)
              dxyz(1:3,1:mesh%nv) = xyztmp(1:3,1:mesh%nv) - mesh%xyz(1:3,1:mesh%nv)
              !xyzprev = mesh%xyz(1:3,1:mesh%nv)!***
              !mesh%xyz(1:3,1:mesh%nv) = mesh%xyz(1:3,1:mesh%nv) + dxyz(1:3,1:mesh%nv)!***
              !call write_tecplot_mesh( &
              !  mesh, &
              !  '../debug/pre_deform4.dat', &
              !  'verif' )
              !mesh%xyz(1:3,1:mesh%nv) = xyzprev!***
              call pre_deformation( &
                   brep_new, &
                   mesh, &
                   dxyz )

              !do ivert = 1,mesh%nv ! <-------------------------------------------------+
              !select case ( mesh%typ(ivert) ) ! <-------------------------------+   !
              !case (0) ! -------------------------------------------------------+   !
              !   mesh%xyz(:,ivert) = brep_new%verts(mesh%ids(ivert))%point%xyz  !   !
              !case (1) ! -------------------------------------------------------+   !
              !   iface = brep_new%edges(mesh%ids(ivert))%halfedges(2)%face      !   !
              !   call eval( &                                                   !   !
              !        mesh%xyz(1:3,ivert), &                                        !   !
              !        brep_new%faces(iface)%surface, &                          !   !
              !        mesh%uv(:,1,ivert) )                                      !   !
              !end select ! <----------------------------------------------------+   !
              !end do ! <---------------------------------------------------------------+

              ! ********
              !call write_xyz_positions( &
              !     '../debug/', &
              !     mesh, &
              !     1 )
              call write_tecplot_mesh( &
                   mesh, &
                   '../debug/pre_deform3.dat', &
                   'pre-deformed' )
              ! ********
              !PAUSE
           END IF

           ! Mesh optimization
           call optim_jiao( &
                brep_new, &
                hypergraph_new%hyperedges(1:hypergraph_new%nhe), &
                hypergraph_new%nhe, &
                mesh, &
                PARAM_frac_conf1, &
                PARAM_frac_conf2, &
                PARAM_ipass1, &
                PARAM_ipass2, &
                PARAM_passmax, &
                options%hmin, &
                options%hmax, &
                .true. )

           call write_tecplot_mesh( &
                mesh, &
                '../debug/pre_deform4.dat', &
                'optimized' )
           !PAUSE
            IF ( .TRUE. ) THEN
               write (strnum, '(i3.3)') instant
               !call discrete_minimum_curvature_radius( &
               !   mesh, &
               !   htarget )
               !htarget = min(htarget, 10._fp*options%hmax)
               !htarget = max(htarget, 0.1_fp*options%hmin)
               !call Hc_correction( &
               !   mesh, &
               !   htarget, &
               !   2._fp, &
               !   10 )
               call target_edge_lengths(mesh, options%hmin, options%hmax, htarget, .true., brep_new)
               !call write_vtk_mesh_solv( &
               !   mesh, &
               !   htarget, &
               !   trim(options%directory) // 'output/debug/hve/instant_'//strnum//'.vtk' )
               call compute_triangle_weights(mesh, 4._fp, 1._fp, tri_weights)
               call write_vtk_mesh_sol( &
                  mesh, &
                  trim(options%directory) // 'output/debug/hve/instant_'//strnum//'.vtk', &
                  tri_weights, &
                  'weight', &
                  htarget, &
                  'hTarget' )
            END IF

           ! Export new positions
           call write_xyz_positions( &
                trim(options%directory) // 'output/', &
                mesh, &
                instant )

           call free_hypergraph(hypergraph_old)
           hypergraph_old => hypergraph_new
        end if

        call free_brep(brep_old)
        call free_intersection_data(interdata_old)


        interdata_old => interdata_new
        brep_old => brep_new

     end if

     if ( time >= options%timespan ) exit main_loop
     !PAUSE

  end do main_loop
  call cpu_time(toc)
  PRINT '(a8,1x,f8.3,1x,a1)','ELAPSED:', toc - tic, 's'

contains


  subroutine make_checkpoint( &
       options, &
       surf, &
       nsurf, &
       interdata, &
       instant, &
       time )
    use mod_util
    use mod_options
    use mod_diffgeom
    use mod_polynomial
    implicit none
    type(type_options),           intent(in) :: options
    integer,                      intent(in) :: nsurf
    type(type_surface), target,   intent(in) :: surf(nsurf)
    integer,                      intent(in) :: instant
    type(type_intersection_data), intent(in) :: interdata
    real,                         intent(in) :: time
    character(3)                             :: strnum3
    integer                                  :: isurf, jsurf, fid, icurv, ncurv, surfpair(2), i

    do isurf = 1,nsurf
       write (strnum3, '(i3.3)') isurf
       call write_polynomial( &
            surf(isurf)%x, &
            trim(options%directory) // 'checkpoint/coef/c_' // strnum3 // '.cheb' )
    end do

    call get_free_unit(fid)

    open( &
         unit=fid, &
         file=trim(options%directory) // 'checkpoint/surftag.dat', &
         action='write' )
    write (fid,*) nsurf
    do isurf = 1,nsurf
       write (fid,*) surf(isurf)%tag
    end do
    close(fid)

    open( &
         unit=fid, &
         file=trim(options%directory) // 'checkpoint/time.dat', &
         action='write' )
    write (fid,*) instant
    write (fid,*) time
    close(fid)

    !call system('cp '  //trim(options%directory) // 'init/tangent_curves.dat ' &
    !     & //trim(options%directory) // 'checkpoint/')
    ncurv = 0
    do icurv = 1,interdata%nc
       if ( interdata%curves(icurv)%smooth ) ncurv = ncurv + 1
    end do
    open( &
         unit=fid, &
         file=trim(options%directory) // 'checkpoint/tangent_curves.dat', &
         action='write' )
    write (fid,*) ncurv
    do icurv = 1,interdata%nc
       if ( interdata%curves(icurv)%smooth ) then
          surfpair(1:2) = 0
          do jsurf = 1,2
             do isurf = 1,nsurf
                if ( associated(interdata%curves(icurv)%surf(jsurf)%ptr, surf(isurf)) ) then
                   surfpair(jsurf) = isurf
                   exit
                end if
             end do
          end do
          write (fid,*) surfpair
          write (fid,*) interdata%curves(icurv)%polyline%np
          do i = 1,interdata%curves(icurv)%polyline%np
             write (fid,'(ES22.15,1X,ES22.15,1X,ES22.15)') interdata%curves(icurv)%polyline%xyz(1:3,i)
          end do
          do i = 1,interdata%curves(icurv)%polyline%np
             write (fid,'(ES22.15,1X,ES22.15,1X,ES22.15,1X,ES22.15)') interdata%curves(icurv)%polyline%uv(1:2,1:2,i)
          end do
       end if
    end do
    close(fid)

  end subroutine make_checkpoint









  subroutine regenerate_feature_paths( &
       brep, &
       hypergraph, &
       mesh, &
       stat )
    implicit none
    LOGICAL, PARAMETER :: DEBUG = .false.
    type(type_brep),         intent(in)    :: brep
    type(type_hypergraph),   intent(in)    :: hypergraph
    type(type_surface_mesh), intent(inout) :: mesh
    real(kind=fp), allocatable             :: uvpath(:,:,:)
    integer, allocatable                   :: idspath(:)
    integer                                :: stat
    real(kind=fp)                          :: xyz(3)
    integer                                :: ihype, ipath, ivert

    do ipath = 1,mesh%npaths
       ! compute curvilinear abscissa on mesh at previous time
       if ( allocated(mesh%paths(ipath)%s) ) deallocate(mesh%paths(ipath)%s)
       allocate(mesh%paths(ipath)%s(mesh%paths(ipath)%nv))
       !mesh%paths(ipath)%s(1) = 0._fp
       ihype = mesh%paths(ipath)%hyperedge
       ivert = hypergraph%hyperedges(ihype)%verts(1)
       mesh%paths(ipath)%s(1) = 0._fp!norm2(mesh%xyz(1:3,mesh%paths(ipath)%verts(1)) - brep%verts(ivert)%point%xyz)
       do ivert = 2,mesh%paths(ipath)%nv
          mesh%paths(ipath)%s(ivert) = mesh%paths(ipath)%s(ivert-1) + &
               norm2(mesh%xyz(:,mesh%paths(ipath)%verts(ivert)) - &
               mesh%xyz(:,mesh%paths(ipath)%verts(ivert-1)))
       end do
       ! normalize
       mesh%paths(ipath)%s = mesh%paths(ipath)%s / mesh%paths(ipath)%s(mesh%paths(ipath)%nv)

       IF ( DEBUG ) THEN
          CALL GET_FREE_UNIT(FID)
          OPEN(UNIT=FID, FILE='../debug/path_old.dat', ACTION='WRITE')
          DO IVERT = 1,mesh%paths(ipath)%nv
             WRITE (FID,*) MESH%XYZ(:,MESH%PATHS(IPATH)%VERTS(IVERT))
          END DO
          CLOSE(FID)
       END IF

       ! regenerate mesh at current time while keeping same normalized abscissa
       allocate(uvpath(2,2,mesh%paths(ipath)%nv), idspath(mesh%paths(ipath)%nv))
       call discretize_hyperedge( &
            brep, &
            hypergraph%hyperedges(ihype), &
            mesh%paths(ipath)%s(1:mesh%paths(ipath)%nv), &
            mesh%paths(ipath)%nv, &
            stat, &
            idspath, &
            uvpath )

       if ( stat > 0 ) then
          PRINT *,'regenerate_feature_paths: failed to regenerate feature path #',ipath
          RETURN
       end if

       mesh%uv(1:2,1:2,mesh%paths(ipath)%verts(2:mesh%paths(ipath)%nv-1)) = &
            uvpath(1:2,1:2,2:mesh%paths(ipath)%nv-1)
       mesh%ids(mesh%paths(ipath)%verts(2:mesh%paths(ipath)%nv-1)) = &
            idspath(2:mesh%paths(ipath)%nv-1)
       if ( hypergraph%hyperedges(ihype)%verts(1) == hypergraph%hyperedges(ihype)%verts(2) ) then
          mesh%uv(1:2,1:2,mesh%paths(ipath)%verts(1)) = uvpath(1:2,1:2,1)
          mesh%ids(mesh%paths(ipath)%verts(1)) = idspath(1)
       end if
       deallocate(uvpath, idspath)

       IF ( DEBUG ) THEN
          CALL GET_FREE_UNIT(FID)
          OPEN(UNIT=FID, FILE='../debug/path_new.dat', ACTION='WRITE')
          DO IVERT = 1,mesh%paths(ipath)%nv
             call eval( &
                  xyz, &
                  brep%edges(mesh%ids(mesh%paths(ipath)%verts(ivert)))%curve%surf(1)%ptr, &
                  mesh%uv(:,1,mesh%paths(ipath)%verts(ivert)) )
             WRITE (FID,*) XYZ
          END DO
          CLOSE(FID)
          PAUSE
       END IF
    end do

  end subroutine regenerate_feature_paths







  subroutine discretize_hyperedge( &
       brep, &
       hyperedge, &
       s, &
       nv, &
       stat, &
       ids, &
       uv )
    USE MOD_UTIL
    use mod_types_brep
    use mod_brep
    use mod_projection
    use mod_intersection
    use mod_diffgeom
    implicit none
    type(type_brep),      intent(in), target :: brep
    type(type_hyperedge), intent(in)         :: hyperedge
    integer,              intent(in)         :: nv
    real(kind=fp),        intent(in)         :: s(nv)
    integer,              intent(out)        :: stat
    integer,              intent(out)        :: ids(nv)
    real(kind=fp),        intent(out)        :: uv(2,2,nv)
    real(kind=fp)                            :: xyz(3), xyztmp(3)
    integer                                  :: stat_intersection
    real(kind=fp)                            :: duv_ds(2,2,2), dxyz_ds(3,2)
    real(kind=fp)                            :: stot, ds, duv(2,2), dxyz(3)
    integer                                  :: stat_proj
    integer                                  :: ihedg(2), ifirst, ilast, sens, np, jedge, kedge, ivert
    INTEGER :: FID

    IF ( .false. ) THEN      
       CALL GET_FREE_UNIT(FID)
       OPEN(UNIT=FID, FILE='../debug/hyperedge.dat', ACTION='WRITE')
       DO JEDGE = 1,HYPEREDGE%NE
          IHEDG = HYPEREDGE%HALFEDGES(:,JEDGE)
          PRINT *,'IHEDG =',IHEDG
          PRINT *,BREP%EDGES(IHEDG(1))%CURVE%POLYLINE%NP
          PRINT *,BREP%EDGES(IHEDG(1))%CURVE%ISPLIT(2,1:BREP%EDGES(IHEDG(1))%CURVE%NSPLIT)
          CALL get_polyline_endpoints( &
               brep, &
               ihedg, &
               ifirst, &
               ilast, &
               sens, &
               np )
          PRINT *,'IFIRST, ILAST =', ifirst, ilast
          DO IVERT = IFIRST,ILAST,(-1)**SENS
             WRITE (FID,*) BREP%EDGES(IHEDG(1))%CURVE%POLYLINE%XYZ(:,IVERT)
          END DO
       END DO
       CLOSE(FID)
    END IF

    ! compute the hyperedge's total curvilinear length
    stot = 0._fp
    do jedge = 1,hyperedge%ne
       ihedg = hyperedge%halfedges(:,jedge)
       call get_polyline_endpoints( &
            brep, &
            ihedg, &
            ifirst, &
            ilast, &
            sens, &
            np )
       if ( jedge == 1 ) xyz = brep%edges(ihedg(1))%curve%polyline%xyz(:,ifirst)
       do ivert = ifirst,ilast,(-1)**sens
          stot = stot + norm2(brep%edges(ihedg(1))%curve%polyline%xyz(:,ivert) - xyz)
          xyz = brep%edges(ihedg(1))%curve%polyline%xyz(:,ivert)
       end do
    end do

    stat = 0

    ihedg = hyperedge%halfedges(:,1)
    ids(1) = ihedg(1)
    kedge = 1
    call get_polyline_endpoints( &
         brep, &
         ihedg, &
         ifirst, &
         ilast, &
         sens, &
         np )
    uv(:,:,1) = brep%edges(ids(1))%curve%polyline%uv(:,:,ifirst)
    xyz = brep%edges(ids(1))%curve%polyline%xyz(:,ifirst)

    do ivert = 2,nv-1
       IF ( .false. ) THEN
          PRINT *,'IVERT =',IVERT
          PRINT *,'XYZprev =',xyz
       END IF
       call diffgeom_intersection( &
            brep%edges(ids(ivert-1))%curve%surf, &
            uv(:,:,ivert-1), &
            duv_ds, &
            dxyz_ds, &
            stat_intersection )
       if ( stat_intersection /= 0 ) then
          stat = 1
          PRINT *,'discretize_hyperedge: not a transversal intersection'
          PAUSE
          return
       end if

       ds = stot * (s(ivert) - s(ivert-1)) * real((-1)**sens, kind=fp)
       duv = ds * duv_ds(:,1,:)
       dxyz = ds * dxyz_ds(:,1)
       IF ( .FALSE. ) PRINT *,'DXYZ =',dxyz

       call projection_hyperedge( &
            brep, &
            hyperedge, &
            ids(ivert-1), &
            uv(:,:,ivert-1), &
            xyz, &
            duv, &
            dxyz, &
            ids(ivert), &
            uv(:,:,ivert), &
            xyztmp, &
            .false., &
            stat_proj )
       if ( stat_proj /= 0 ) then
          stat = 2
          PRINT *,'discretize_hyperedge: failed to project onto hyperedge'
          PRINT *,'DXYZ =',dxyz
          PAUSE
          return
       end if

       if ( ids(ivert) /= ids(ivert-1) ) then
          do jedge = kedge+1,hyperedge%ne
             if ( hyperedge%halfedges(1,jedge) == ids(ivert) ) then
                call get_polyline_endpoints( &
                     brep, &
                     hyperedge%halfedges(:,jedge), &
                     ifirst, &
                     ilast, &
                     sens, &
                     np )
                kedge = jedge
                exit
             end if
          end do
       end if

       call eval( &
            xyz, &
            brep%edges(ids(ivert))%curve%surf(1)%ptr, &
            uv(:,1,ivert) )
    end do

    ihedg = hyperedge%halfedges(:,1)
    ids(nv) = ihedg(1)
    call get_polyline_endpoints( &
         brep, &
         ihedg, &
         ifirst, &
         ilast, &
         sens, &
         np )
    uv(:,:,nv) = brep%edges(ids(nv))%curve%polyline%uv(:,:,ilast)

  end subroutine discretize_hyperedge
















  subroutine update_mesh_correspondance( &
       brep_old, &
       brep_new, &
       hypg_old, &
       hypg_new, &
       mesh, &
       stat )
    use mod_types_brep
    use mod_hypergraph
    use mod_mesh
    use mod_halfedge
    implicit none
    type(type_brep),         intent(in)    :: brep_old
    type(type_brep),         intent(in)    :: brep_new
    type(type_hypergraph),   intent(in)    :: hypg_old
    type(type_hypergraph),   intent(inout) :: hypg_new
    type(type_surface_mesh), intent(inout) :: mesh
    integer,                 intent(out)   :: stat
    integer                                :: f2f(brep_old%nf)
    integer                                :: hf2hf(hypg_old%nhf)
    integer                                :: v2v(brep_old%nv)
    logical, allocatable                   :: assigned(:)
    integer, allocatable                   :: faces_old(:), faces_new(:)
    integer                                :: nfaces_old, nfaces_new
    integer                                :: he2he(hypg_new%nhe)
    integer                                :: iface, jface
    integer                                :: ihypf
    integer                                :: ivert, jvert, kvert
    integer                                :: ihype, jhype
    integer                                :: ipath, nv
    integer                                :: sens
    integer                                :: iedge, jedge
    real(kind=fp)                          :: disti, distj, xyzv(3)

    ! BREP face ->  BREP face
    if ( brep_old%nf /= brep_new%nf ) then
       stat = 2
       return
    end if
    f2f(1:brep_old%nf) = [(iface, iface=1,brep_old%nf)] ! ***

    ! hyperface -> hyperface
    hf2hf(1:hypg_old%nhf) = [(ihypf, ihypf=1,hypg_old%nhf)] ! ***

    ! BREP vertex ->  BREP vertex
    v2v(1:brep_old%nv) = 0
    allocate(assigned(brep_new%nv))
    assigned(1:brep_new%nv) = .false.
    do ivert = 1,brep_old%nv ! <--------------------------------------------------+
       call get_v2f( &                                                            !
            brep_old, &                                                           !
            ivert, &                                                              !
            faces_old, &                                                          !
            nfaces_old )                                                          !
       faces_old(1:nfaces_old) = f2f(faces_old(1:nfaces_old))                     !
       !                                                                          !
       do jvert = 1,brep_new%nv ! <--------------------------------------------+  !
          if ( assigned(jvert) ) cycle                                         !  !
          call get_v2f( &                                                      !  !
               brep_new, &                                                     !  !
               jvert, &                                                        !  !
               faces_new, &                                                    !  !
               nfaces_new )                                                    !  !
          if ( nfaces_new /= nfaces_old ) cycle
          !                                                                    !  !
          do iface = 1,nfaces_new ! <---------------------------------------+  !  !
             if ( faces_new(iface) == faces_old(1) ) then ! <------------+  !  !  !
                assigned(jvert) = .true.                                 !  !  !  !
                do jface = 1,nfaces_new-1 ! <-------------------------+  !  !  !  !
                   if ( faces_new(1+mod(iface+jface-1,nfaces_new)) &  !  !  !  !  !
                        /= faces_old(jface+1) ) then ! <-----------+  !  !  !  !  !
                      assigned(jvert) = .false.                    !  !  !  !  !  !
                      exit                                         !  !  !  !  !  !
                   end if ! <--------------------------------------+  !  !  !  !  !
                end do ! <--------------------------------------------+  !  !  !  !
                if ( assigned(jvert) ) v2v(ivert) = jvert                !  !  !  !
                exit                                                     !  !  !  !
             end if ! <--------------------------------------------------+  !  !  !
          end do ! <--------------------------------------------------------+  !  !
          !                                                                    !  !
          if ( v2v(ivert) > 0 ) exit                                           !  !
       end do ! <--------------------------------------------------------------+  !
    end do ! <--------------------------------------------------------------------+
    if ( allocated(faces_old) ) deallocate(faces_old)
    if ( allocated(faces_new) ) deallocate(faces_new)

    PRINT *,'v2v:'
    DO IVERT = 1,BREP_OLD%NV
       PRINT *, IVERT, ' -->', V2V(IVERT)
    END DO

    do ivert = 1,mesh%nv ! <------------------------------------------------------+
       if ( mesh%typ(ivert) == 0 ) then ! <------------------------------------+  !
          if ( v2v(mesh%ids(ivert)) > 0 ) then ! <--------------------------+  !  !
             mesh%ids(ivert) = v2v(mesh%ids(ivert))                         !  !  !
          else ! -----------------------------------------------------------+  !  !
             ! ...                                                          !  !  !
             ! nouveau typ = 1 ou 2                                         !  !  !
          end if ! <--------------------------------------------------------+  !  !
       end if ! <--------------------------------------------------------------+  !
    end do ! <--------------------------------------------------------------------+


    ! hyperedge -> hyperedge
    deallocate(assigned)
    allocate(assigned(hypg_new%nhe))
    assigned(1:hypg_new%nhe) = .false.
    he2he(1:hypg_old%nhe) = 0
    do ihype = 1,hypg_old%nhe ! <---------------------------------------------------------------------+
       if ( hypg_old%hyperedges(ihype)%verts(1) /= hypg_old%hyperedges(ihype)%verts(2) ) then ! <--+  !
          ! open hyperedge (with 2 distinct endpoints)                                             !  !
          new_open_hyperedges : do jhype = 1,hypg_new%nhe ! <--------------------------------+     !  !
             if ( assigned(jhype) ) cycle                                                    !     !  !
             if ( hypg_new%hyperedges(jhype)%verts(1) == &                                   !     !  !
                  hypg_new%hyperedges(jhype)%verts(2) ) cycle                                !     !  !
             !                                                                               !     !  !
             do sens = 1,2 ! <------------------------------------------------------------+  !     !  !
                if ( v2v(hypg_old%hyperedges(ihype)%verts(1)) == &                        !  !     !  !
                     hypg_new%hyperedges(jhype)%verts(sens) .and. &                       !  !     !  !
                     v2v(hypg_old%hyperedges(ihype)%verts(2)) == &                        !  !     !  !
                     hypg_new%hyperedges(jhype)%verts(1+mod(sens,2)) .and. &              !  !     !  !
                     hf2hf(hypg_old%hyperedges(ihype)%hyperfaces(1)) == &                 !  !     !  !
                     hypg_new%hyperedges(jhype)%hyperfaces(sens) .and. &                  !  !     !  !
                     hf2hf(hypg_old%hyperedges(ihype)%hyperfaces(2)) == &                 !  !     !  !
                     hypg_new%hyperedges(jhype)%hyperfaces(1+mod(sens,2)) ) then ! <---+  !  !     !  !
                   assigned(jhype) = .true.                                            !  !  !     !  !
                   he2he(ihype) = -jhype*(-1)**sens                                    !  !  !     !  !
                   exit new_open_hyperedges                                            !  !  !     !  !
                end if ! <-------------------------------------------------------------+  !  !     !  !
             end do ! <-------------------------------------------------------------------+  !     !  !
          end do new_open_hyperedges ! <-----------------------------------------------------+     !  !
       else ! -------------------------------------------------------------------------------------+  !
          ! closed hyperedge                                                                       !  !
          new_closed_hyperedges : do jhype = 1,hypg_new%nhe ! <------------------------------+     !  !
             if ( assigned(jhype) ) cycle                                                    !     !  !
             if ( hypg_new%hyperedges(jhype)%verts(1) /= &                                   !     !  !
                  hypg_new%hyperedges(jhype)%verts(2) ) cycle                                !     !  !
             do sens = 1,2 ! <------------------------------------------------------------+  !     !  !
                if ( hf2hf(hypg_old%hyperedges(ihype)%hyperfaces(1)) == &                 !  !     !  !
                     hypg_new%hyperedges(jhype)%hyperfaces(sens) .and. &                  !  !     !  !
                     hf2hf(hypg_old%hyperedges(ihype)%hyperfaces(2)) == &                 !  !     !  !
                     hypg_new%hyperedges(jhype)%hyperfaces(1+mod(sens,2))) then ! <----+  !  !     !  !
                   assigned(jhype) = .true.                                            !  !  !     !  !
                   he2he(ihype) = -jhype*(-1)**sens                                    !  !  !     !  !
                   do iedge = 1,hypg_old%hyperedges(ihype)%ne ! <-------------------+  !  !  !     !  !
                      ivert = v2v(get_orig( &                                       !  !  !  !     !  !
                           brep_old, &                                              !  !  !  !     !  !
                           hypg_old%hyperedges(ihype)%halfedges(1:2,iedge) ))       !  !  !  !     !  !
                      do jedge = 1,hypg_new%hyperedges(jhype)%ne ! <-------------+  !  !  !  !     !  !
                         jvert = get_orig( &                                     !  !  !  !  !     !  !
                              brep_new, &                                           !  !  !  !  !     !  !
                              hypg_new%hyperedges(jhype)%halfedges(1:2,jedge) )     !  !  !  !  !     !  !
                         if ( ivert == jvert ) then ! <-----------------------+  !  !  !  !  !     !  !
                            hypg_new%hyperedges(jhype)%verts(1:2) = jvert     !  !  !  !  !  !     !  !
                            call circular_permutation( &                      !  !  !  !  !  !     !  !
                                 hypg_new%hyperedges(jhype)%halfedges&        !  !  !  !  !  !     !  !
                                 (1,1:hypg_new%hyperedges(jhype)%ne), &       !  !  !  !  !  !     !  !
                                 hypg_new%hyperedges(jhype)%ne, &             !  !  !  !  !  !     !  !
                                 jedge )                                      !  !  !  !  !  !     !  !
                            call circular_permutation( &                      !  !  !  !  !  !     !  !
                                 hypg_new%hyperedges(jhype)%halfedges&        !  !  !  !  !  !     !  !
                                 (2,1:hypg_new%hyperedges(jhype)%ne), &       !  !  !  !  !  !     !  !
                                 hypg_new%hyperedges(jhype)%ne, &             !  !  !  !  !  !     !  !
                                 jedge )                                      !  !  !  !  !  !     !  !
                            exit new_closed_hyperedges                        !  !  !  !  !  !     !  !
                         end if ! <-------------------------------------------+  !  !  !  !  !     !  !
                      end do ! <-------------------------------------------------+  !  !  !  !     !  !
                   end do ! <-------------------------------------------------------+  !  !  !     !  !
                   STOP 'update_mesh_correspondance: failed to find a common vertex'   !  !  !     !  !
                end if ! <-------------------------------------------------------------+  !  !     !  !
             end do ! <-------------------------------------------------------------------+  !     !  !
          end do new_closed_hyperedges ! <---------------------------------------------------+     !  !
       end if ! <----------------------------------------------------------------------------------+  !
    end do ! <----------------------------------------------------------------------------------------+

    PRINT *,'he2he:'
    DO IHYPE = 1,HYPG_OLD%NHE
       print *, IHYPE, ' -->', HE2HE(IHYPE)
    END DO

    if ( .not.all(assigned)) then
       stat = 1
       return
    end if


    stat = 0
    do ipath = 1,mesh%npaths ! <---------------------------------------+
       ihype = mesh%paths(ipath)%hyperedge                             !
       if ( he2he(ihype) /= 0 ) then ! <----------------------------+  !
          nv = mesh%paths(ipath)%nv                                 !  !
          mesh%paths(ipath)%hyperedge = abs(he2he(ihype))           !  !
          if ( he2he(ihype) < 0 ) then ! <-----------------------+  !  !
             mesh%paths(ipath)%verts(1:nv) = &                   !  !  !
                  mesh%paths(ipath)%verts(nv:1:-1)               !  !  !
          end if ! <---------------------------------------------+  !  !
          if ( hypg_old%hyperedges(ihype)%verts(1) == &             !  !
               hypg_old%hyperedges(ihype)%verts(2) ) then ! <----+  !  !
             ivert = 0                                           !  !  !
             disti = huge(1._fp)                                 !  !  !
             jvert = hypg_old%hyperedges(ihype)%verts(1)         !  !  !
             xyzv = brep_old%verts(jvert)%point%xyz              !  !  !
             do jvert = 1,mesh%paths(ipath)%nv-1 ! <----------+  !  !  !
                kvert = mesh%paths(ipath)%verts(jvert)        !  !  !  !
                distj = sum((mesh%xyz(1:3,kvert) - xyzv)**2)  !  !  !  !
                if ( distj < disti ) then ! <---+             !  !  !  !
                   ivert = jvert                !             !  !  !  !
                   disti = distj                !             !  !  !  !
                end if ! <----------------------+             !  !  !  !
             end do ! <---------------------------------------+  !  !  !
             kvert = mesh%paths(ipath)%verts(ivert)              !  !  !
             call circular_permutation( &                        !  !  !
                  mesh%paths(ipath)%verts(1:nv-1), &             !  !  !
                  nv-1, &                                        !  !  !
                  ivert )                                        !  !  !
             mesh%paths(ipath)%verts(nv) = kvert                 !  !  !
             !PRINT *,mesh%paths(ipath)%verts(1:nv)
             !PAUSE
          end if ! <---------------------------------------------+  !  !
       else ! ------------------------------------------------------+  !
          ! delete path and update ids & typ                        !  !
          ! ...   
       end if ! <---------------------------------------------------+  !
    end do ! <---------------------------------------------------------+

  end subroutine update_mesh_correspondance






  subroutine new_regen_paths( &
       brep, &
       hypg, &
       mesh, &
       stat )
    use mod_types_brep
    use mod_brep
    use mod_projection
    use mod_intersection
    use mod_diffgeom
    implicit none
    LOGICAL, PARAMETER :: DEBUG = .false.
    INTEGER            :: FID
    type(type_brep),         intent(in)    :: brep
    type(type_hypergraph),   intent(in)    :: hypg
    type(type_surface_mesh), intent(inout) :: mesh
    integer,                 intent(out)   :: stat
    real(kind=fp)                          :: dist, dist_min, uv(2,2), uv_min(2,2), xyz(3), xyz_min(3), xyztmp(3)
    integer                                :: iedge, iedge_min
    integer                                :: stat_intersection
    real(kind=fp)                          :: duv_ds(2,2,2), dxyz_ds(3,2)
    real(kind=fp)                          :: ds, duv(2,2), dxyz(3)
    integer                                :: stat_proj
    integer                                :: ipath, ihype, jvert, ivert, kvert

    IF ( DEBUG ) PRINT *,'>>>> NEW_REGEN_FEATURES'

    do ipath = 1,mesh%npaths
       IF ( DEBUG ) THEN
          PRINT *,''
          PRINT *,'PATH #', ipath

          CALL GET_FREE_UNIT(FID)
          OPEN(UNIT=FID, FILE='../debug/path_old.dat', ACTION='WRITE')
          DO IVERT = 1,mesh%paths(ipath)%nv
             WRITE (FID,*) MESH%XYZ(:,MESH%PATHS(IPATH)%VERTS(IVERT))
          END DO
          CLOSE(FID)

          CALL GET_FREE_UNIT(FID)
          OPEN(UNIT=FID, FILE='../debug/path_new.dat', ACTION='WRITE')
       END IF
       ihype = mesh%paths(ipath)%hyperedge

       do jvert = 1,mesh%paths(ipath)%nv
          ivert = mesh%paths(ipath)%verts(jvert)
          IF ( DEBUG ) THEN
             PRINT *,''
             PRINT *,'VERT #', IVERT, '( TYP =', mesh%typ(ivert), ')'
             PRINT *,'   XYZ =', mesh%xyz(1:3,ivert)
          END IF

          if ( mesh%typ(ivert) == 0 ) then
             ! ENDPOINT OF AN OPEN PATH (i.e. either jvert = 1 or mesh%paths(ipath)%nv)
             ! (?)
             cycle
          elseif ( mesh%typ(ivert) == 1 ) then

             ! find closest polyline point
             IF ( JVERT == 1 ) THEN
                dist_min = huge(1._fp)
                do iedge = 1,hypg%hyperedges(ihype)%ne
                   call distance_from_edge( &
                        brep, &
                        hypg%hyperedges(ihype)%halfedges(1:2,iedge), &
                        mesh%xyz(1:3,ivert), &
                        dist, &
                        uv, &
                        xyz )
                   !IF ( DEBUG ) PRINT *,'      DIST =', DIST
                   if ( dist < dist_min ) then
                      dist_min = dist
                      uv_min = uv
                      xyz_min = xyz
                      iedge_min = hypg%hyperedges(ihype)%halfedges(1,iedge)
                   end if
                end do
             ELSE
                !kvert = mesh%paths(ipath)%verts(jvert-1)
                iedge_min = mesh%ids(kvert)
                uv_min = mesh%uv(1:2,1:2,kvert)
                !xyz_min = mesh%xyz(1:3,kvert)
                call eval( &
                     xyz_min, &
                     brep%edges(iedge_min)%curve%surf(1)%ptr, &
                     uv_min(1:2,1) )
                dist_min = sum( (mesh%xyz(1:3,ivert) - xyz_min)**2 )
             END IF
             kvert = ivert

             IF ( DEBUG ) THEN
                PRINT *,'   DIST_MIN =', SQRT(DIST_MIN)
                PRINT *,'     UV_MIN =', UV_MIN
                !CYCLE
             END IF

             ! refine projection
             IF ( .FALSE.) THEN
                mesh%uv(1:2,1:2,ivert) = uv_min
                mesh%ids(ivert) = iedge_min
             ELSE
                call diffgeom_intersection( &
                     brep%edges(iedge_min)%curve%surf, &
                     uv_min, &
                     duv_ds, &
                     dxyz_ds, &
                     stat_intersection )
                if ( stat_intersection /= 0 ) then
                   stat = 1
                   PRINT *,'new_regen_paths: not a transversal intersection'
                   PRINT *,'STAT =', STAT
                   PAUSE
                   return
                else
                   IF ( DEBUG ) THEN
                      !PRINT *,'   DUV_DS  =', duv_ds(1:2,1,1:2)
                      PRINT *,'         DXYZ_DS =', dxyz_ds(1:3,1)
                      PRINT *,'   XYZ - XYZ_MIN =', mesh%xyz(1:3,ivert) - xyz_min
                   END IF
                end if

                ds = dot_product(dxyz_ds(1:3,1), mesh%xyz(1:3,ivert) - xyz_min)
                duv = ds * duv_ds(1:2,1,1:2)
                dxyz = ds * dxyz_ds(1:3,1)
                IF ( DEBUG ) THEN
                   PRINT *,'   |DS|   =', ABS(DS)
                   PRINT *,'   |DUV|  =', NORM2(DUV(1:2,1)), NORM2(DUV(1:2,2))
                   PRINT *,'   |DXYZ| =', NORM2(DXYZ)
                END IF

                !CYCLE
                call projection_hyperedge( &
                     brep, &
                     hypg%hyperedges(ihype), &
                     iedge_min, &
                     uv_min, &
                     xyz_min, &
                     duv, &
                     dxyz, &
                     mesh%ids(ivert), &
                     mesh%uv(1:2,1:2,ivert), &
                     xyztmp, &
                     .false., &
                     stat_proj )
                if ( stat_proj /= 0 ) then
                   stat = 2
                   PRINT *,'new_regen_paths: failed to project onto hyperedge'
                   PRINT *,'DXYZ =',dxyz
                   CLOSE(FID)
                   PAUSE
                   return
                end if
             END IF

             IF ( DEBUG ) THEN
                PRINT *,'   IDS:', IEDGE_MIN, ' ->', mesh%ids(ivert)
                PRINT *,'   UV =', mesh%uv(1:2,1:2,ivert)

                xyz_min = xyztmp
                !call eval( &
                !     xyz_min, &
                !     brep%edges(mesh%ids(ivert))%curve%surf(1)%ptr, &
                !     mesh%uv(1:2,1,ivert) )
                WRITE (FID,*) xyz_min
             END IF

          end if
       end do

       IF ( DEBUG ) THEN
          CLOSE(FID)
          PAUSE
       END IF
    end do

    IF ( DEBUG ) PRINT *,'NEW_REGEN_FEATURES >>>>'


  end subroutine new_regen_paths







  subroutine distance_from_edge( &
       brep, &
       ihedg, &
       xyz, &
       dist_min, &
       uv_min, &
       xyz_min )
    use mod_types_intersection
    implicit none
    type(type_brep), intent(in), target       :: brep
    integer,         intent(in)               :: ihedg(2)
    real(kind=fp),   intent(in)               :: xyz(3)
    real(kind=fp),   intent(out)              :: dist_min
    real(kind=fp),   intent(out)              :: uv_min(2,2)
    real(kind=fp),   intent(out)              :: xyz_min(3)
    type(type_intersection_polyline), pointer :: poly => null()
    real(kind=fp)                             :: dist
    integer                                   :: ihead, itail, sens, np, iplv

    call get_polyline_endpoints( &
         brep, &
         ihedg, &
         ihead, &
         itail, &
         sens, &
         np )

    dist_min = huge(1._fp)
    poly => brep%edges(ihedg(1))%curve%polyline
    do iplv = ihead,itail,sens
       dist = sum( (xyz - poly%xyz(1:3,iplv))**2 )
       if ( dist < dist_min ) then
          dist_min = dist
          uv_min = poly%uv(1:2,1:2,iplv)
          xyz_min = poly%xyz(1:3,iplv)
       end if
    end do

  end subroutine distance_from_edge

end program fftsurf
