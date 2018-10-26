program fftsurf

  USE MOD_UTIL
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
  integer, parameter                    :: freq_checkpoint = 30
  integer, parameter                    :: PARAM_passmax = 20!40
  real(kind=fp), parameter              :: PARAM_frac_conf1 = 1._fp
  real(kind=fp), parameter              :: PARAM_frac_conf2 = 0.7_fp
  integer, parameter                    :: PARAM_ipass1 = 0!5
  integer, parameter                    :: PARAM_ipass2 = 0!15
  ! --------------------------------------------------------------
  ! Options
  type(type_options)                    :: options
  integer                               :: narg, iarg, irep, io
  character(99)                         :: fileoptions, arg
  ! --------------------------------------------------------------
  ! BREP
  type(type_surface), allocatable       :: surf(:)
  integer                               :: nsurf, isurf
  !type(type_intersection_data), pointer :: interdata => null()
  !type(type_brep),              pointer :: brep => null()
  type(type_intersection_data) :: interdata
  type(type_brep)              :: brep
  ! --------------------------------------------------------------
  ! Hypergraph
  !type(type_hypergraph), pointer        :: hypergraph => null()
  type(type_hypergraph)        :: hypergraph
  logical, allocatable                  :: feat_edge(:)
  logical, allocatable                  :: feat_vert(:)
  integer                               :: ihe, ihf, ied, ifa
  ! --------------------------------------------------------------
  ! Surface mesh
  type(type_surface_mesh)               :: mesh
  integer                               :: ivert, iface, ipath
  integer, allocatable                  :: idspath(:)
  real(kind=fp), allocatable            :: uvpath(:,:,:)
  real(kind=fp)                         :: xyzpath(3)
  integer                               :: stat_regenpath
  ! --------------------------------------------------------------
  ! Physical time
  real                                  :: time
  integer                               :: instant
  ! --------------------------------------------------------------
  ! CPU time
  real                                  :: tic, toc
  ! --------------------------------------------------------------
  character(3) :: strnum3, strnum
  INTEGER :: FID

  ! Lecture argument (nom du fichier configuration)
  narg = command_argument_count()
  if ( narg < 1 ) then
     ! fichier config par defaut
     fileoptions = 'default.opt'
  else
     call get_command_argument(1, fileoptions)
  end if


  ! Parse options supplementaires (override)
  iarg = 2
  do while(iarg <= narg)
     call get_command_argument(iarg, arg)
     print *,'arg =',arg
     select case(arg)
     case ('-rep')
        print *,'here'
        iarg = iarg + 1
        options%reprise = .true.
        print *,options%reprise
     end select
  end do


  

  !! Initialization
  PRINT *,'INITIALIZATION...'
  call cpu_time(tic)
  !if ( options%mode > 1 ) allocate(interdata, brep, hypergraph)
  call initialization( &
       fileoptions, &
       options, &
       surf, &
       nsurf, &
       interdata, &
       brep, &
       hypergraph, &
       mesh )
  call cpu_time(toc)
  PRINT '(a8,1x,f8.3,1x,a1)','ELAPSED:', toc - tic, 's'


  !! Main loop
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
     ! Export initial mesh
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




  
  PRINT *,'MAIN LOOP...'
  call cpu_time(tic)
  main_loop : do

     ! checkpoint
     if ( mod(instant,freq_checkpoint) == 0 ) then
        do isurf = 1,nsurf
           write (strnum3, '(i3.3)') isurf
           call write_polynomial( &
                surf(isurf)%x, &
                trim(options%directory) // 'checkpoint/coef/c_' // strnum3 // '.cheb' )           
        end do

        call get_free_unit(fid)

        !open(unit=fid, file=trim(options%directory) // 'checkpoint/nsurf.dat', action='write')
        open(unit=fid, file=trim(options%directory) // 'checkpoint/surftag.dat', action='write')
        write (fid,*) nsurf
        do isurf = 1,nsurf
           write (fid,*) surf(isurf)%tag
        end do
        close(fid)

        open(unit=fid, file=trim(options%directory) // 'checkpoint/time.dat', action='write')
        write (fid,*) instant
        write (fid,*) time
        close(fid)

        call system('cp '  //trim(options%directory) // 'init/tangent_curves.dat ' //trim(options%directory) // 'checkpoint/')
     end if
     


     
     ! Propagate surfaces
     do isurf = 1,nsurf
        !if ( all([4,5] /= isurf) ) cycle ! CAS 'COMPENSATEUR'
        if ( surf(isurf)%tag /= 1 ) cycle
        call propagation_step_RK4( &
             surf(isurf), &
             options%propagation_law, &
             options%timestep, &
             time )
     end do
     time = time + options%timestep
     instant = instant + 1
     !PRINT '(a5,1x,f5.3,1x,f8.3)','TIME:', time, '/', options%timespan
     PRINT *,''
     PRINT '(a7,1x,i0)','INSTANT',INSTANT
     PRINT *,'TIME:', time, '/', options%timespan

     !***************
     write (strnum, '(i3.3)') instant
     do isurf = 1,nsurf
        write (strnum3, '(i3.3)') isurf
        call write_polynomial( &
             surf(isurf)%x, &
             trim(options%directory) // 'output/debug/instant_' // strnum // 'c_' // strnum3 // '.cheb' )
     end do
     !***************

     if ( options%mode > 0 ) then
         ! update the tangent curves' polylines
        call update_intersection_curves(interdata)

     IF ( .TRUE. ) THEN ! CAS 'COMPENSATEUR'
        brep%nv = 0; brep%ne = 0; brep%nf = 0;
        deallocate(brep%verts, brep%edges, brep%faces)

        interdata%np = 0; interdata%nc = 0;
        deallocate(interdata%points, interdata%curves)

        call make_brep( &
             surf, &
             nsurf, &
             interdata, &
             brep, &
             options%chord_err, &
             options%hmin, &
             options%hmax )

        IF ( .TRUE. ) THEN
           call write_brep_files( &
                brep, &
                trim(options%directory) // 'output/debug/brep/verts.dat', &
                trim(options%directory) // 'output/debug/brep/edges.dat', &
                trim(options%directory) // 'output/debug/brep/faces.dat' )

           call write_intersection_data( &
                interdata, &
                trim(options%directory) // 'output/debug/brep/intersection_points.dat', &
                trim(options%directory) // 'output/debug/brep/intersection_curves.dat' )
        END IF

        ! set edge -> hyperedge
        do ihe = 1,hypergraph%nhe
           do ied = 1,hypergraph%hyperedges(ihe)%ne
              brep%edges(hypergraph%hyperedges(ihe)%halfedges(1,ied))%hyperedge = ihe
           end do
        end do

        ! set face -> hyperface
        do ihf = 1,hypergraph%nhf
           do ifa = 1,hypergraph%hyperfaces(ihf)%nf
              brep%faces(hypergraph%hyperfaces(ihf)%faces(ifa))%hyperface = ihf
           end do
        end do
        
     END IF
     end if
     
     if ( .FALSE. ) THEN!options%mode > 0 ) then
        !deallocate(interdata, brep)
        !allocate(interdata, brep)
        ! Regenerate BREP
        call make_brep( &
             surf, &
             nsurf, &
             interdata, &
             brep, &
             options%chord_err, &
             options%hmin, &
             options%hmax )

        ! Make hypergraph
        !deallocate(hypergraph)
        !allocate(hypergraph)
        call make_hypergraph( &
             brep, &
             hypergraph, &
             feat_edge, &
             feat_vert )
     end if

     if ( options%mode > 1 ) then
        ! (Hypergraph isomorphism)

        ! (Regenerate features mesh)
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
                stat_regenpath, &
                idspath, &
                uvpath )

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
                      xyzpath, &
                      brep%edges(mesh%ids(mesh%paths(ipath)%verts(ivert)))%curve%surf(1)%ptr, &
                      mesh%uv(:,1,mesh%paths(ipath)%verts(ivert)) )
                 WRITE (FID,*) XYZPATH
              END DO
              CLOSE(FID)
              PAUSE
           END IF
        end do

        ! (Pre-deform mesh to prevent inverted elements)
        do ivert = 1,mesh%nv ! <-------------------------------------------------+
           select case ( mesh%typ(ivert) ) ! <-------------------------------+   !
           case (0) ! -------------------------------------------------------+   !
              mesh%xyz(:,ivert) = brep%verts(mesh%ids(ivert))%point%xyz      !   !
           case (1) ! -------------------------------------------------------+   !
              iface = brep%edges(mesh%ids(ivert))%halfedges(2)%face          !   !
              call eval( &                                                   !   !
                   mesh%xyz(:,ivert), &                                      !   !
                   brep%faces(iface)%surface, &                              !   !
                   mesh%uv(:,1,ivert) )                                      !   !
           case (2) ! -------------------------------------------------------+   !
              iface = mesh%ids(ivert)                                        !   !
              call eval( &                                                   !   !
                   mesh%xyz(:,ivert), &                                      !   !
                   brep%faces(iface)%surface, &                              !   !
                   mesh%uv(:,1,ivert) )                                      !   !
           end select ! <----------------------------------------------------+   !
        end do ! <---------------------------------------------------------------+

        ! Mesh optimization
        IF ( .true. ) THEN
           call optim_jiao( &
                brep, &
                hypergraph%hyperedges(1:hypergraph%nhe), &
                hypergraph%nhe, &
                mesh, &
                PARAM_frac_conf1, &
                PARAM_frac_conf2, &
                PARAM_ipass1, &
                PARAM_ipass2, &
                PARAM_passmax, &
                options%hmin, &
                options%hmax )
        ELSE
           call laplacian_smoothing( &
                brep, &
                hypergraph%hyperedges(1:hypergraph%nhe), &
                hypergraph%nhe, &
                mesh, &
                PARAM_passmax )
           !call equilateral_triangle_smoothing( &
           !     brep, &
           !     hypergraph%hyperedges, &
           !     hypergraph%nhe, &
           !     mesh, &
           !     PARAM_passmax )
        END IF

        ! Export mesh (update positions)
        call write_xyz_positions( &
             trim(options%directory) // 'output/', &
             mesh, &
             instant )

     end if

     if ( time >= options%timespan ) exit main_loop

  end do main_loop

  call cpu_time(toc)
  PRINT '(a8,1x,f8.3,1x,a1)','ELAPSED:', toc - tic, 's'


  do isurf = 1,nsurf
     write (strnum3, '(i3.3)') isurf
     call write_polynomial( &
          surf(isurf)%x, &
          trim(options%directory) // 'output/final/coef/c_' // strnum3 // '.cheb' )
  end do
  





  
contains



  

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
    type(type_surface), allocatable, intent(out)   :: surf(:)
    integer,                         intent(out)   :: nsurf
    type(type_intersection_data),    intent(out)   :: interdata
    type(type_brep),                 intent(out)   :: brep
    type(type_hypergraph),           intent(out)   :: hypergraph
    type(type_surface_mesh),         intent(out)   :: mesh
    character(99)                                  :: dir
    character(3)                                   :: strnum3
    integer                                        :: fid
    integer                                        :: isurf, icurv, i, j
    real(kind=fp)                           :: xyzverif(3,2)
    
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
         file = trim(dir) // 'surftag.dat', &!file = trim(dir) // 'nsurf.dat', &
         action = 'read')
    read (fid,*) nsurf
    allocate(surf(nsurf))
    do isurf = 1,nsurf
       read (fid,*) surf(isurf)%tag
    end do
    close(fid)
    
    !allocate(surf(nsurf))
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
       
       call write_brep_files( &
            brep, &
            trim(dir) // 'brep/verts.dat', &
            trim(dir) // 'brep/edges.dat', &
            trim(dir) // 'brep/faces.dat' )

       call write_intersection_data( &
            interdata, &
            trim(dir) // 'brep/intersection_points.dat', &
            trim(dir) // 'brep/intersection_curves.dat' )
       
       ! Make hypergraph
       call make_hypergraph( &
            brep, &
            hypergraph, &
            feat_edge, &
            feat_vert )
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
       
    end if

    if ( options%mode > 1 ) then
       IF ( .false. ) THEN
          call init_mesh( &
               trim(options%directory) // 'init/mesh/xyzb.dat', &
               trim(options%directory) // 'init/mesh/trib.dat', &
               brep, &
               mesh )

          !***********
          call write_connectivity( &
               trim(options%directory) // 'output/', &
               mesh, &
               1 )
          call write_xyz_positions( &
               trim(options%directory) // 'output/', &
               mesh, &
               0 )
          !***********
          
       ELSE
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

          call write_connectivity( &
               '../debug/', &
               mesh, &
               0 )
          call write_xyz_positions( &
               '../debug/', &
               mesh, &
               0 )

          ! Eliminate edges that are too short
          call contract_small_edges( &
               brep, &
               hypergraph%hyperedges(1:hypergraph%nhe), &
               hypergraph%nhe, &
               mesh, &
               0.25_fp*options%hmin )

          call write_connectivity( &
               '../debug/', &
               mesh, &
               1 )
          call write_xyz_positions( &
               '../debug/', &
               mesh, &
               1 )

          !***********
          !call write_connectivity( &
          !     trim(options%directory) // 'output/', &
          !     mesh )
          !call write_xyz_positions( &
          !     trim(options%directory) // 'output/', &
          !     mesh, &
          !     0 )
          !***********
          
          ! Make mesh halfedges
          call make_halfedges( &
               mesh )
       END IF

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
                IF ( MAX(norm2(mesh%xyz(:,i) - xyzverif(:,1)), norm2(mesh%xyz(:,i) - xyzverif(:,2))) > 1.D-9 ) THEN
                   PRINT *,'I =', I,', ERR =', norm2(mesh%xyz(:,i) - xyzverif(:,1)), norm2(mesh%xyz(:,i) - xyzverif(:,2))
                   do j = 1,2
                      iface = brep%edges(mesh%ids(i))%halfedges(j)%face
                      call eval( &
                           xyzverif(:,j), &
                           brep%faces(iface)%surface, &
                           mesh%uv(:,j,i) )
                   end do
                   PRINT *,'         ', norm2(mesh%xyz(:,i) - xyzverif(:,1)), norm2(mesh%xyz(:,i) - xyzverif(:,2))
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
                IF ( MAX(norm2(mesh%xyz(:,i) - xyzverif(:,1)), norm2(mesh%xyz(:,i) - xyzverif(:,2))) > 1.D-9 ) THEN
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
  
    end if

  end subroutine initialization













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
    open(unit=finf, file='../tmp/info.dat', action='write')
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
       call write_polynomial(brep%faces(iface)%surface%x, '../tmp/c.cheb')
       !
       ! open boundary points & edges files
       call get_free_unit(fpts)
       open(unit=fpts, file='../tmp/bpts.dat', action='write')
       call get_free_unit(fedg)
       open(unit=fedg, file='../tmp/bedg.dat', action='write')
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
             ! polyligne uv dans le sens du contour de la face
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
                bv2mv(ivert) = mesh%nv + 1                                       !
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
                      PRINT *,'***',BV2MV(IVERT), sens, ihedg(1), 1, uvf(1:2,1:2,1)
                      mesh%ids(bv2mv(ivert)) = ihedg(1)
                      mesh%typ(bv2mv(ivert)) = 1
                      mesh%uv(1:2,1:2,bv2mv(ivert)) = uvf(1:2,[2,1],1)!uvf(1:2,1:2,1)
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
                write (fpts,*) uvf(:,1,i)                !
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
            & ../tmp/c.cheb &
            & ../tmp/bpts.dat &
            & ../tmp/bedg.dat &
            & ../tmp/info.dat &
            & ../tmp/tri.dat &
            & ../tmp/uv.dat &
            & ../tmp/xyz.dat' )
       write (strnum,'(i3.3)') iface
       !call system('cp tmp/c.cheb Jouke/meshgen/coeffs/c_'//strnum//'.cheb')
       !call system('cp tmp/tri.dat Jouke/meshgen/brepmesh/uv/tri_'//strnum//'.dat')
       !call system('cp tmp/uv.dat Jouke/meshgen/brepmesh/uv/uv_'//strnum//'.dat')
       !
       !PAUSE
       
       ! read face submesh                                                                     !
       call read_triangles( &                                                                  !
            '../tmp/tri.dat', &                                                                   !
            trif, &                                                                            !
            ntrif )                                                                            !
       !
       call read_points( &                                                                     !
            '../tmp/uv.dat', &                                                                    !
            uvi, &                                                                             !
            2, &                                                                               !
            np )                                                                               !
       if ( allocated(uvf) ) deallocate(uvf)
       allocate(uvf(2,2,np-npf))
       uvf(1:2,1,1:np-npf) = uvi(1:2,npf+1:np)
       !
       call read_points( &                                                                     !
            '../tmp/xyz.dat', &                                                                   !
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
            [(brep%faces(iface)%hyperface, i=1,ntrif)], &
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
            [(brep%faces(iface)%hyperface, i=1,ntrif)], &
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
         !'Jouke/meshgen/brepmesh/brepmesh.dat', &
         '../brepmesh.dat', &
         'mesh_brep', &
         t2bf )
   
    
  end subroutine generate_brep_mesh






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
























  subroutine init_mesh( &
       filexyz, &
       filetri, &
       brep, &
       mesh )
    implicit none
    character(*),            intent(in)  :: filexyz, filetri
    type(type_brep),         intent(in)  :: brep
    type(type_surface_mesh), intent(out) :: mesh
    logical, allocatable                 :: projected(:)
    integer                              :: stat_proj
    integer                              :: ivert

    call read_triangles( &
         filetri, &
         mesh%tri, &
         mesh%nt )
    
    call read_points( &
         filexyz, &
         mesh%xyz, &
         3, &
         mesh%nv )
    
    PRINT *,MESH%NV,' VERTICES,', MESH%NT,' TRIANGLES'
    call make_halfedges(mesh)
    
    allocate(mesh%uv(2,2,mesh%nv), mesh%ids(mesh%nv), mesh%typ(mesh%nv))
    mesh%typ(:) = 2

    ! project mesh vertices onto brep surface
    allocate(projected(mesh%nv))
    projected(:) = .false.
    PRINT *,'PROJECT FIRST VERTEX..'
    do ivert = 1,mesh%nv
       PRINT *,'  IVERT =',IVERT
       call project_first_vertex_init( &
            brep, &
            mesh%xyz(:,ivert), &
            stat_proj, &
            mesh%uv(:,1,ivert), &
            mesh%ids(ivert) )
       if ( stat_proj == 0 ) then
          projected(ivert) = .true.
          exit
       end if
    end do
    
    call project_mesh_vertices_init( &
         brep, &
         mesh, &
         ivert, &
         projected )

    PRINT *,'FAILED PROJECTIONS:',COUNT(.not.PROJECTED),'/',MESH%NV
    PAUSE

    do ivert = 1,mesh%nv
       if ( .not.projected(ivert) ) cycle
       call eval( &
            mesh%xyz(:,ivert), &
            brep%faces(mesh%ids(ivert))%surface, &
            mesh%uv(:,1,ivert) )
    end do
    
    deallocate(projected) 
    
  end subroutine init_mesh



  
  subroutine project_first_vertex_init( &
       brep, &
       xyz, &
       stat_proj, &
       uv, &
       iface )
    use mod_types_brep
    use mod_regiontree
    implicit none
    type(type_brep), target, intent(in)  :: brep
    real(kind=fp),           intent(in)  :: xyz(3)
    integer,                 intent(out) :: stat_proj
    real(kind=fp),           intent(out) :: uv(2)
    integer,                 intent(out) :: iface
    type(type_region)                    :: region
    integer                              :: i
    
    do iface = 1,brep%nf
       PRINT *,'    IFACE =',IFACE
       call init_region( &
            region, &
            2, &
            [( [-1._fp, 1._fp], i=1,2 )] )
       
       stat_proj = 1
       call recurs_proj( &
            brep%faces(iface)%surface, &
            region, &
            xyz, &
            stat_proj, &
            uv )
       call free_region_tree(region)
       
       if ( stat_proj == 0 ) return
    end do
    
  end subroutine project_first_vertex_init



  
  recursive subroutine recurs_proj( &
       surf, &
       region, &
       xyz, &
       stat_proj, &
       uv )
    use mod_diffgeom
    use mod_polynomial
    use mod_regiontree
    use mod_obb
    use mod_projection
    implicit none
    LOGICAL, PARAMETER :: DEBUG = .TRUE.
    real(kind=fp), parameter             :: TOLplansqr = (1.d-4)**2
    type(type_surface),    intent(in)    :: surf
    type(type_region),     intent(inout) :: region
    real(kind=fp),         intent(in)    :: xyz(3)
    integer,               intent(inout) :: stat_proj
    real(kind=fp),         intent(out)   :: uv(2)
    type(type_polynomial)                :: polyregion
    type(type_obb)                       :: box
    integer                              :: m, n
    real(kind=fp)                        :: nonplanarity, uv_subdiv(2), xyz_proj(3)
    real(kind=fp), dimension(2)          :: lowerb , upperb
    integer                              :: stat_subdiv
    integer                              :: ichild
    
    if ( stat_proj == 0 ) return

    IF ( DEBUG ) PRINT *,'      UVBOX =',REGION%UVBOX

    call chgvar2( &
       surf%x, &
       polyregion, &
       region%uvbox([1,3]), &
       region%uvbox([2,4]) )

    ! compute bounding box
    call chebOBB2( &
         polyregion%coef, &
         polyregion%degr, &
         box )

    ! check whether the point is inside the bounding box
    if ( .not.is_inside_OBB(xyz, box, EPSmath) ) return

    ! check whether the surface region is "planar enough"
    m = polyregion%degr(1)+1
    n = polyregion%degr(2)+1
    nonplanarity = sum( (&
         sum(abs(polyregion%coef(3:m,1,1:3)),1) + &
         sum(abs(polyregion%coef(1,3:n,1:3)),1) + &
         sum(sum(abs(polyregion%coef(2:m,2:n,1:3)),1),1) &
         )**2 )
    IF ( DEBUG ) PRINT *,'      NONPLAN =',SQRT(NONPLANARITY)
    if ( nonplanarity < TOLplansqr ) then
       lowerb = region%uvbox([1,3]) - EPSuv
       upperb = region%uvbox([2,4]) + EPSuv
       uv = 0.5_fp * (region%uvbox([1,3]) + region%uvbox([2,4]))
       call projection_surface( &
            surf, &
            xyz, &
            uv, &
            lowerb, &
            upperb, &
            stat_proj, &
            xyz_proj )
       return
    end if
    
    ! carry on recursion
    uv_subdiv = 0.5_fp * (region%uvbox([1,3]) + region%uvbox([2,4]))
    call subdiv_region( &
       region, &
       uv_subdiv, &
       stat_subdiv, &
       EPSregion )

    if ( stat_subdiv < 2 ) then
       do ichild = 1,size(region%child)
          call recurs_proj( &
               surf, &
               region%child(ichild), &
               xyz, &
               stat_proj, &
               uv )
          if ( stat_proj == 0 ) return
       end do
    end if
     
  end subroutine recurs_proj



  

  
  recursive subroutine project_mesh_vertices_init( &
       brep, &
       mesh, &
       ivert, &
       projected )
    use mod_types_brep
    use mod_mesh
    use mod_halfedge
    use mod_diffgeom
    use mod_linalg
    use mod_projection
    implicit none
    type(type_brep),         intent(in)    :: brep
    type(type_surface_mesh), intent(inout) :: mesh
    integer,                 intent(in)    :: ivert
    logical,                 intent(inout) :: projected(mesh%nv)
    real(kind=fp)                          :: dxyz(3), jac(3,2), duv(2)
    logical                                :: singular
    integer                                :: stat_proj
    real(kind=fp)                          :: xyzproj(3)
    integer                                :: ihedg(2), iface, jface, jvert, ivar

    ihedg = mesh%v2h(:,ivert)
    iface = get_face(ihedg)
    adjacent_verts : do
       jvert = get_dest(mesh, ihedg)
       if ( projected(jvert) ) then
          dxyz = mesh%xyz(:,ivert) - mesh%xyz(:,jvert)
          do ivar = 1,2
             call evald1( &
                  jac(:,ivar), &
                  brep%faces(mesh%ids(jvert))%surface, &
                  mesh%uv(:,1,jvert), &
                  ivar )
          end do
          call solve_NxN( &
               duv, &
               matmul(transpose(jac), jac), &
               matmul(transpose(jac), dxyz), &
               singular )

          call projection_hyperface( &
               brep, &
               mesh%ids(jvert), &
               mesh%uv(:,1,jvert), &
               mesh%xyz(:,jvert), &
               duv, &
               dxyz, &
               mesh%ids(ivert), &
               mesh%uv(:,1,ivert), &
               .false., &
               stat_proj )

          call projection_surface( &
               brep%faces(mesh%ids(ivert))%surface, &
               mesh%xyz(:,ivert), &
               mesh%uv(:,1,ivert), &
               spread(-[1._fp + EPSuv], dim=1, ncopies=2), &
               spread( [1._fp + EPSuv], dim=1, ncopies=2), &
               stat_proj, &
               xyzproj )

          if ( stat_proj == 0 ) then
             projected(ivert) = .true.
             exit adjacent_verts
          end if
       end if

       ihedg = get_prev(ihedg)
       ihedg = get_twin(mesh, ihedg)
       jface = get_face(ihedg)
       if ( jface < 1 .or. jface == iface ) exit adjacent_verts
    end do adjacent_verts

    if ( .not.projected(ivert) ) return

    ihedg = mesh%v2h(:,ivert)
    iface = get_face(ihedg)
    do
       jvert = get_dest(mesh, ihedg)
       if ( .not.projected(jvert) ) then
          call project_mesh_vertices_init( &
               brep, &
               mesh, &
               jvert, &
               projected )
       end if
       ihedg = get_prev(ihedg)
       ihedg = get_twin(mesh, ihedg)
       jface = get_face(ihedg)
       if ( jface < 1 .or. jface == iface ) exit
    end do

  end subroutine project_mesh_vertices_init

  








  




  







  !subroutine write_quadratic_mfem_mesh( &
  !     mesh, &
  !     filename )
  !  use mod_mesh
  !  use mod_util
  !  implicit none
  !  character(7), parameter             :: fmt = 'ES13.6'
  !  type(type_surface_mesh), intent(in) :: mesh
  !  character(*),            intent(in) :: filename
  !  integer                             :: fid
  !  integer                             :: nb
  !  logical                             :: explored(3,mesh%nt)
  !  integer                             :: it, ie, iv, jv, j
  !  
  !  call get_free_unit(fid)
  !
  !  open( &
  !       unit = fid, &
  !       file = filename, &
  !       action = 'write' )
  !
  !  write (fid,*) 'MFEM mesh v1.0'
  !  write (fid,*) 'dimension'
  !  write (fid,*) '2'
  !
  !  write (fid,*) 'elements'
  !  write (fid,*) mesh%nt
  !  do it = 1,mesh%nt
  !     write (fid,*) '1 2 ',mesh%tri(:,it)
  !  end do
  !
  !  nb = 0
  !  do it = 1,mesh%nt
  !     do ie = 1,3
  !        if ( mesh%twin(2,ie,it) < 1 ) nb = nb + 1
  !     end do
  !  end do
  !  write (fid,*) 'boundary'
  !  write (fid,*) nb
  !  if ( nb > 0 ) then
  !     do it = 1,mesh%nt
  !        do ie = 1,3
  !           if ( mesh%twin(2,ie,it) < 1 ) write (fid,*) mesh%tri([ie,1+mod(ie,3)],it)
  !        end do
  !     end do
  !  end if
  !  
  !  write (fid,*) 'vertices'
  !  write (fid,*) mesh%nv
  !
  !  write (fid,*) 'nodes'
  !  write (fid,*) 'FiniteElementSpace'
  !  write (fid,*) 'FiniteElementCollection: Quadratic'
  !  write (fid,*) 'VDim: 3'
  !  write (fid,*) 'Ordering: 1'
  !
  !  do iv = 1,mesh%nv
  !     do j = 1,3
  !        write ('('//fmt//',1x)', advance='no') mesh%xyz(j,iv)
  !     end do
  !  end do
  !  
  !  explored(:,:) = .false.
  !  do it = 1,mesh%nt
  !     do ie = 1,3
  !        if ( explored(ie,it) ) cycle
  !        iv = mesh%tri(ie,it)
  !        jv = mesh%tri(1+mod(ie,3),it)
  !     end do
  !  end do
  !
  !  ! midpoints...
  !  
  !  close(fid)
  !  
  !end subroutine write_quadratic_mfem_mesh





  !subroutine get_midpoint( &
  !     brep, &
  !     mesh, &
  !     verts, &
  !     ids, &
  !     uv, &
  !     xyz )
  !  implicit none
  !  type(type_brep),         intent(in)  :: brep
  !  type(type_surface_mesh), intent(in)  :: mesh
  !  integer,                 intent(in)  :: verts(2)
  !  integer,                 intent(out) :: ids
  !  real(kind=fp),           intent(out) :: uv(2,2)
  !  real(kind=fp),           intent(out) :: xyz(3)
  !  integer                              :: typ(2)
  !  typ = mesh%typ(verts)
  !  if ( all(typ == 2) ) then
  !     if ( mesh%ids(verts(1)) == mesh%ids(verts(2)) ) then
  !        ids = mesh%ids(verts(1))
  !        uv = 0.5_fp * (mesh%uv(:,1,verts(1)) + mesh%uv(:,1,verts(2)))
  !        call eval( &
  !             xyz, &
  !             brep%faces(ids)%surface, &
  !             uv )
  !     else
  !        ids = mesh%ids(verts(1))
  !        uv = mesh%uv(:,1,verts(1))
  !     end if
  !  end if
  !end subroutine get_midpoint




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
          CALL get_polyline_endpoints( &
               brep, &
               ihedg, &
               ifirst, &
               ilast, &
               sens, &
               np )
          DO IVERT = IFIRST,ILAST,(-1)**SENS
             WRITE (FID,*) BREP%EDGES(IHEDG(1))%CURVE%POLYLINE%XYZ(:,IVERT)
          END DO
       END DO
       CLOSE(FID)
    END IF

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
       IF ( .FALSE. ) THEN
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
  
  
end program fftsurf
