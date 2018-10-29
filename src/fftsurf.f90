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

  implicit none

  ! --------------------------------------------------------------
  ! Parameters
  integer, parameter                    :: freq_checkpoint = 10
  integer, parameter                    :: PARAM_passmax = 20!40
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
  integer                               :: stat_corresp, stat_regen_path
  integer                               :: ivert, iface
  ! --------------------------------------------------------------
  ! Physical time
  real                                  :: time
  integer                               :: instant
  ! --------------------------------------------------------------
  ! CPU time
  real                                  :: tic, toc
  ! --------------------------------------------------------------
  integer :: i, j
  

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
     end select
  end do



  !! Initialization
  PRINT *,'INITIALIZATION...'
  call cpu_time(tic)
  !if ( options%mode > 0 ) allocate(interdata_old, brep_old)
  !if ( options%mode > 1 ) allocate(hypergraph_old)
  allocate(interdata_old, brep_old, hypergraph_old)
  call initialization( &
       fileoptions, &
       options, &
       surf, &
       nsurf, &
       interdata_old, &
       brep_old, &
       hypergraph_old, &
       mesh )
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
     PAUSE
  end if


  !! Main loop
  PRINT *,'MAIN LOOP...'
  call cpu_time(tic)
  if ( options%mode > 0 ) allocate(interdata_new, brep_new)
  if ( options%mode > 1 ) allocate(hypergraph_new)
  main_loop : do
     ! checkpoint
     if ( mod(instant,freq_checkpoint) == 0 ) then
        call make_checkpoint( &
             options, &
             surf(1:nsurf), &
             nsurf, &
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


     if ( options%mode > 0 ) then
        ! reset brep data
        call free_brep(brep_new)
        call free_intersection_data(interdata_new)

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

           call update_mesh_correspondance( &
                brep_old, &
                brep_new, &
                hypergraph_old, &
                hypergraph_new, &
                mesh, &
                stat_corresp )
           if ( stat_corresp > 0 ) then
              call make_checkpoint( &
                   options, &
                   surf(1:nsurf), &
                   nsurf, &
                   instant, &
                   time )
              PRINT *,'failed to update mesh correspondance, made a checkpoint'
              STOP
           end if
           !PAUSE

           ! Regenerate features mesh
           call regenerate_feature_paths( &
                brep_new, &
                hypergraph_new, &
                mesh, &
                stat_regen_path )
           if ( stat_regen_path > 0 ) then
              PRINT *,'failed to regenerate feature paths'
              PAUSE
           end if

           ! (Pre-deform mesh to prevent inverted elements)
           do ivert = 1,mesh%nv ! <-------------------------------------------------+
              select case ( mesh%typ(ivert) ) ! <-------------------------------+   !
              case (0) ! -------------------------------------------------------+   !
                 mesh%xyz(:,ivert) = brep_new%verts(mesh%ids(ivert))%point%xyz  !   !
              case (1) ! -------------------------------------------------------+   !
                 iface = brep_new%edges(mesh%ids(ivert))%halfedges(2)%face      !   !
                 call eval( &                                                   !   !
                      mesh%xyz(:,ivert), &                                      !   !
                      brep_new%faces(iface)%surface, &                          !   !
                      mesh%uv(:,1,ivert) )                                      !   !
              case (2) ! -------------------------------------------------------+   !
                 iface = mesh%ids(ivert)                                        !   !
                 call eval( &                                                   !   !
                      mesh%xyz(:,ivert), &                                      !   !
                      brep_new%faces(iface)%surface, &                          !   !
                      mesh%uv(:,1,ivert) )                                      !   !
              end select ! <----------------------------------------------------+   !
           end do ! <---------------------------------------------------------------+

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
                options%hmax )

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

  end do main_loop
  call cpu_time(toc)
  PRINT '(a8,1x,f8.3,1x,a1)','ELAPSED:', toc - tic, 's'

contains


  subroutine make_checkpoint( &
       options, &
       surf, &
       nsurf, &
       instant, &
       time )
    use mod_util
    use mod_options
    use mod_diffgeom
    use mod_polynomial
    implicit none
    type(type_options), intent(in) :: options
    integer,            intent(in) :: nsurf
    type(type_surface), intent(in) :: surf(nsurf)
    integer,            intent(in) :: instant
    real,               intent(in) :: time
    character(3)                   :: strnum3
    integer                        :: isurf, fid

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

    call system('cp '  //trim(options%directory) // 'init/tangent_curves.dat ' &
         & //trim(options%directory) // 'checkpoint/')

  end subroutine make_checkpoint


  


  subroutine initialization( &
       fileoptions, &
       options, &
       surf, &
       nsurf, &
       interdata, &
       brep, &
       hypergraph, &
       mesh )
    use mod_util
    use mod_options
    use mod_diffgeom
    use mod_types_intersection
    use mod_types_brep
    use mod_intersection
    use mod_brep
    use mod_mesh
    use mod_polynomial
    use mod_optimmesh
    use mod_tolerances
    ! mode = 0 : propagation des surfaces
    !        1 : regeneration BREP + hypergraphe
    !        2 : isomorphisme hypergrap
    implicit none
    character(*),                    intent(in)    :: fileoptions
    type(type_options),              intent(inout) :: options
    type(type_surface), allocatable, intent(inout) :: surf(:)
    integer,                         intent(out)   :: nsurf
    type(type_intersection_data),    intent(inout) :: interdata
    type(type_brep),                 intent(inout) :: brep
    type(type_hypergraph),           intent(inout) :: hypergraph
    type(type_surface_mesh),         intent(out)   :: mesh
    character(99)                                  :: dir
    character(3)                                   :: strnum3
    integer                                        :: fid
    integer                                        :: isurf, icurv, i, j

    ! Read options file
    call read_options( &
         fileoptions, &
         options )
    call print_options( &
         options)

    dir = trim(options%directory)
    if ( options%reprise ) then
       dir = trim(dir) // 'checkpoint/'
    else
       dir = trim(dir) // 'init/'
    end if
    print *,'dir =',trim(dir)

    call get_free_unit(fid)

    ! Import surfaces
    open( &
         unit = fid, &
         file = trim(dir) // 'surftag.dat', &
         action = 'read')
    read (fid,*) nsurf
    allocate(surf(nsurf))
    do isurf = 1,nsurf
       read (fid,*) surf(isurf)%tag
    end do
    close(fid)

    do isurf = 1,nsurf
       write (strnum3,'(i3.3)') isurf
       call read_polynomial( &
            surf(isurf)%x, &
            trim(dir) // 'coef/c_' // strnum3 // '.cheb', &
            nvar=2, &
            base=1 )

       call economize2( &
            surf(isurf)%x, &
            EPSmath )

       call compute_deriv1(surf(isurf))
       call compute_deriv2(surf(isurf))
       call compute_pseudonormal(surf(isurf))
       call economize2( &
            surf(isurf)%pn, &
            EPSmath )
    end do

    if ( options%mode > 0 ) then
       ! Read initial tangential intersection curves
       call read_intersection_curves( &
            trim(dir) // 'tangent_curves.dat', &
            surf, &
            interdata )
       do icurv = 1,interdata%nc
          interdata%curves(icurv)%smooth = .true.
       end do

       ! Generate initial BREP model
       call make_brep( &
            surf, &
            nsurf, &
            interdata, &
            brep, &
            options%chord_err, &
            options%hmin, &
            options%hmax )

       ! debugging >> ..................
       call write_brep_files( &
            brep, &
            trim(dir) // 'brep/verts.dat', &
            trim(dir) // 'brep/edges.dat', &
            trim(dir) // 'brep/faces.dat' )

       call write_intersection_data( &
            interdata, &
            trim(dir) // 'brep/intersection_points.dat', &
            trim(dir) // 'brep/intersection_curves.dat' )
       ! .................. <<

       ! Make hypergraph
       call make_hypergraph( &
            brep, &
            hypergraph, &
            feat_edge, &
            feat_vert )

       ! debugging >> ..................
       call get_free_unit( fid )
       open(unit=fid, file=trim(dir) // 'brep/hyperfaces.dat', action='write')
       write (fid,*) hypergraph%nhf
       do i = 1,hypergraph%nhf
          write (fid,*) hypergraph%hyperfaces(i)%nf
          write (fid,*) hypergraph%hyperfaces(i)%faces(1:hypergraph%hyperfaces(i)%nf)
       end do
       close(fid)

       open(unit=fid, file=trim(dir) // 'brep/hyperedges.dat', action='write')
       write (fid,*) hypergraph%nhe
       do i = 1,hypergraph%nhe
          write (fid,*) hypergraph%hyperedges(i)%ne
          write (fid,*) hypergraph%hyperedges(i)%verts
          write (fid,*) hypergraph%hyperedges(i)%hyperfaces
          do j = 1,hypergraph%hyperedges(i)%ne
             write (fid,*) hypergraph%hyperedges(i)%halfedges(1:2,j)
          end do
       end do
       close(fid)
       ! .................. <<

    end if

    if ( options%mode > 1 ) then
       !IF ( .false. ) THEN
       !   call init_mesh( &
       !        trim(options%directory) // 'init/mesh/xyzb.dat', &
       !        trim(options%directory) // 'init/mesh/trib.dat', &
       !        brep, &
       !        mesh )
       !   !***********
       !   call write_connectivity( &
       !        trim(options%directory) // 'output/', &
       !        mesh, &
       !        1 )
       !   call write_xyz_positions( &
       !        trim(options%directory) // 'output/', &
       !        mesh, &
       !        0 )
       !   !***********
       !ELSE
       ! Generate a first mesh that conforms to the BREP
       call generate_brep_mesh( &
            brep, &
            options%hmin, &
            options%hmax, &
            options%chord_err, &
            feat_edge(1:brep%ne), &
            feat_vert(1:brep%nv), &
            hypergraph%hyperedges(1:hypergraph%nhe), &
            hypergraph%nhe, &
            mesh )

       ! debugging >> ..................
       call write_connectivity( &
            '../debug/', &
            mesh, &
            0 )
       call write_xyz_positions( &
            '../debug/', &
            mesh, &
            0 )
       ! .................. <<

       ! Eliminate edges that are too short
       call contract_small_edges( &
            brep, &
            hypergraph%hyperedges(1:hypergraph%nhe), &
            hypergraph%nhe, &
            mesh, &
            0.25_fp*options%hmin )

       ! debugging >> ..................
       call write_connectivity( &
            '../debug/', &
            mesh, &
            1 )
       call write_xyz_positions( &
            '../debug/', &
            mesh, &
            1 )
       ! .................. <<

       ! Make mesh halfedges
       call make_halfedges( &
            mesh )
       !END IF

       ! debugging >> ..................
       call write_mesh_files( &
            mesh, &
            trim(dir) // 'mesh/tri.dat', &
            trim(dir) // 'mesh/xyz.dat', &
            trim(dir) // 'mesh/uv.dat', &
            trim(dir) // 'mesh/idstyp.dat', &
            trim(dir) // 'mesh/paths.dat' )

       call get_free_unit(fid)
       open(unit=fid, file=trim(dir) // 'mesh/mv2h.dat', action='write')
       do i = 1,mesh%nv
          write (fid,*) mesh%v2h(:,i)
       end do
       close(fid)
       open(unit=fid, file=trim(dir) // 'mesh/mtwin.dat', action='write')
       do i = 1,mesh%nt
          write (fid,*) mesh%twin(:,:,i)
       end do
       close(fid)
       ! .................. <<

       ! CHECK UVs
       IF ( .true. ) call check_uvs(brep, mesh)

       PAUSE

       ! Mesh smoothing
       call optim_jiao( &
            brep, &
            hypergraph%hyperedges(1:hypergraph%nhe), &
            hypergraph%nhe, &
            mesh, &
            1._fp, &!PARAM_frac_conf1, &
            0.7_fp, &!PARAM_frac_conf2, &
            5, &!PARAM_ipass1, &
            15, &!PARAM_ipass2, &
            30, &
            options%hmin, &
            options%hmax )


       ! CHECK UVs
       IF ( .true. ) call check_uvs(brep, mesh)

    end if

  end subroutine initialization






  subroutine check_uvs( &
       brep, &
       mesh )
    use mod_types_brep
    use mod_mesh
    use mod_diffgeom
    use mod_tolerances
    implicit none
    type(type_brep),         intent(in) :: brep
    type(type_surface_mesh), intent(in) :: mesh
    real(kind=fp)                       :: xyzverif(3,2)
    integer                             :: i, j, iface

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
          IF ( MAX(norm2(mesh%xyz(:,i) - xyzverif(:,1)), norm2(mesh%xyz(:,i) - xyzverif(:,2))) > EPSxyz ) THEN
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

  end subroutine check_uvs






  subroutine regenerate_feature_paths( &
       brep, &
       hypergraph, &
       mesh, &
       stat )
    implicit none
    type(type_brep),         intent(in)    :: brep
    type(type_hypergraph),   intent(in)    :: hypergraph
    type(type_surface_mesh), intent(inout) :: mesh
    real(kind=fp), allocatable             :: uvpath(:,:,:)
    integer, allocatable                   :: idspath(:)
    integer                                :: stat
    real(kind=fp)                          :: xyz(3)
    integer                                :: ipath, ivert

    do ipath = 1,mesh%npaths
       ! compute curvilinear abscissa on mesh at previous time
       if ( allocated(mesh%paths(ipath)%s) ) deallocate(mesh%paths(ipath)%s)
       allocate(mesh%paths(ipath)%s(mesh%paths(ipath)%nv))
       mesh%paths(ipath)%s(1) = 0._fp
       do ivert = 2,mesh%paths(ipath)%nv
          mesh%paths(ipath)%s(ivert) = mesh%paths(ipath)%s(ivert-1) + &
               norm2(mesh%xyz(:,mesh%paths(ipath)%verts(ivert)) - &
               mesh%xyz(:,mesh%paths(ipath)%verts(ivert-1)))
       end do
       ! normalize
       mesh%paths(ipath)%s = mesh%paths(ipath)%s / mesh%paths(ipath)%s(mesh%paths(ipath)%nv)

       IF ( .FALSE. ) THEN
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
            hypergraph%hyperedges(mesh%paths(ipath)%hyperedge), &
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
       deallocate(uvpath, idspath)

       IF ( .FALSE. ) THEN
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
    real(kind=fp)                            :: xyz(3)
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





  


  subroutine write_connectivity( &
       dir, &
       mesh, &
       irep )
    use mod_util
    use mod_mesh
    implicit none
    character(*),            intent(in) :: dir
    type(type_surface_mesh), intent(in) :: mesh
    integer,                 intent(in) :: irep
    character(2)                        :: strnum
    integer                             :: fid, i, j

    write (strnum,'(i2.2)') irep
    call get_free_unit(fid)
    open( &
         unit = fid, &
         file = dir // 'connect_' // strnum //'.dat', &
         action = 'write' )
    do j = 1,mesh%nt
       do i = 1,3
          write (fid, '(i0,1x)', advance='no') mesh%tri(i,j)
       end do
       write (fid,*)
    end do
    close(fid)

  end subroutine write_connectivity



  subroutine write_face_ref( &
       dir, &
       mesh, &
       irep )
    use mod_util
    use mod_mesh
    implicit none
    character(*),            intent(in) :: dir
    type(type_surface_mesh), intent(in) :: mesh
    integer,                 intent(in) :: irep
    character(2)                        :: strnum
    integer                             :: fid, i

    write (strnum,'(i2.2)') irep
    call get_free_unit(fid)
    open( &
         unit = fid, &
         file = dir // 'faceref_' // strnum //'.dat', &
         action = 'write' )
    do i = 1,mesh%nt
       write (fid, '(i0)') mesh%ihf(i)
    end do
    close(fid)

  end subroutine write_face_ref



  subroutine write_xyz_positions( &
       dir, &
       mesh, &
       instant )
    use mod_util
    use mod_mesh
    implicit none
    character(*),            intent(in) :: dir
    type(type_surface_mesh), intent(in) :: mesh
    integer,                 intent(in) :: instant
    character(3)                        :: str
    integer                             :: fid, i, j

    write (str,'(i3.3)') instant
    call get_free_unit(fid)
    open( &
         unit = fid, &
         file = dir // 'pos_' // str // '.dat', &
         action = 'write' )
    do j = 1,mesh%nv
       do i = 1,3
          write (fid, '(e22.15,1x)', advance='no') mesh%xyz(i,j)
       end do
       write (fid,*)
    end do
    close(fid)

  end subroutine write_xyz_positions








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
    integer                                :: iface, jface, ihypf, ivert, jvert, ihype, jhype, ipath

    ! BREP face ->  BREP face
    f2f(1:brep_old%nf) = [(iface, iface=1,brep_old%nf)] ! ***

    ! hyperface -> hyperface
    hf2hf(1:hypg_old%nhf) = [(ihypf, ihypf=1,hypg_old%nhf)] ! ***
    
    ! BREP vertex ->  BREP vertex
    v2v(1:brep_old%nv) = 0
    allocate(assigned(brep_new%nv))
    assigned(1:brep_new%nv) = .false.
    do ivert = 1,brep_old%nv
       call get_v2f( &
            brep_old, &
            ivert, &
            faces_old, &
            nfaces_old )
       faces_old(1:nfaces_old) = f2f(faces_old(1:nfaces_old))
       
       do jvert = 1,brep_new%nv
          if ( assigned(jvert) ) cycle
          call get_v2f( &
               brep_new, &
               jvert, &
               faces_new, &
               nfaces_new )

          do iface = 1,nfaces_new
             if ( faces_new(iface) == faces_old(1) ) then
                assigned(jvert) = .true.
                do jface = 1,nfaces_new-1
                   if ( faces_new(1+mod(iface+jface-1,nfaces_new)) /= faces_old(jface+1) ) then
                      assigned(jvert) = .false.
                      exit
                   end if
                end do
                if ( assigned(jvert) ) v2v(ivert) = jvert
                exit
             end if
          end do

          if ( v2v(ivert) > 0 ) exit
       end do
    end do
    if ( allocated(faces_old) ) deallocate(faces_old)
    if ( allocated(faces_new) ) deallocate(faces_new)

    print *,'v2v:'
    do ivert = 1,brep_old%nv
       print *, ivert, ' -->', v2v(ivert)
    end do

    ! hyperedge -> hyperedge
    deallocate(assigned)
    allocate(assigned(hypg_new%nhe))
    assigned(1:hypg_new%nhe) = .false.
    he2he(1:hypg_old%nhe) = 0
    do ihype = 1,hypg_old%nhe
       if ( hypg_old%hyperedges(ihype)%verts(1) /= hypg_old%hyperedges(ihype)%verts(2) ) then
          do jhype = 1,hypg_new%nhe
             if ( assigned(jhype) ) cycle
             if ( v2v(hypg_old%hyperedges(ihype)%verts(1)) == hypg_new%hyperedges(jhype)%verts(1) .and. &
                  v2v(hypg_old%hyperedges(ihype)%verts(2)) == hypg_new%hyperedges(jhype)%verts(2) .and. &
                  hf2hf(hypg_old%hyperedges(ihype)%hyperfaces(1)) == hypg_new%hyperedges(jhype)%hyperfaces(1) .and. &
                  hf2hf(hypg_old%hyperedges(ihype)%hyperfaces(2)) == hypg_new%hyperedges(jhype)%hyperfaces(2) ) then
                assigned(jhype) = .true.
                he2he(ihype) = jhype
                exit
             elseif &
                  ( v2v(hypg_old%hyperedges(ihype)%verts(1)) == hypg_new%hyperedges(jhype)%verts(2) .and. &
                  v2v(hypg_old%hyperedges(ihype)%verts(2)) == hypg_new%hyperedges(jhype)%verts(1) .and. &
                  hf2hf(hypg_old%hyperedges(ihype)%hyperfaces(1)) == hypg_new%hyperedges(jhype)%hyperfaces(2) .and. &
                  hf2hf(hypg_old%hyperedges(ihype)%hyperfaces(2)) == hypg_new%hyperedges(jhype)%hyperfaces(1) ) then
                assigned(jhype) = .true.
                he2he(ihype) = -jhype
                call reverse_hyperedge(hypg_new%hyperedges(jhype))
                exit
             end if
          end do
       else
          ! closed hyperedge (cycle) ...
       end if
    end do

    print *,'he2he:'
    do ihype = 1,hypg_old%nhe
       print *, ihype, ' -->', he2he(ihype)
    end do

    if ( .not.all(assigned)) then
       stat = 1
       return
    end if


    stat = 0
    do ipath = 1,mesh%npaths
       ihype = mesh%paths(ipath)%hyperedge
       if ( he2he(ihype) > 0 ) then
          mesh%paths(ipath)%hyperedge = he2he(ihype)
       else
          ! delete path ans update ids & typ
       end if
    end do
    
  end subroutine update_mesh_correspondance




  

end program fftsurf
