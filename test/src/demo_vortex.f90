program demo_vortex

  use mod_util
  use mod_options
  use mod_math
  use mod_types_intersection
  use mod_types_brep
  use mod_brep
  use mod_hypergraph
  use mod_mesh
  use mod_tolerances
  use mod_propagation
  use mod_polynomial
  use mod_intersection
  use mod_init

  implicit none

  character(5) :: pref = 'data/'!''
  integer, parameter              :: nsnap = 200!6
  
  !real                            :: ftsnap(nsnap) = [0., 1./8., 1./4., 3./8., 1./2., 1.]
  real                            :: ftsnap(nsnap)
  
  ! Options
  type(type_options)              :: options

  ! BRep
  type(type_surface), allocatable :: surf(:)
  integer                         :: nsurf, isurf, icurv, ihead, itail, sens, np
  type(type_intersection_data)    :: interdata
  type(type_brep)                 :: brep

  ! Surface mesh
  type(type_surface_mesh)         :: mesh

  ! Hypergraph
  type(type_hypergraph)           :: hypergraph
  logical, allocatable            :: feat_edge(:), feat_vert(:)

  ! Physical time
  real                            :: time
  integer                         :: instant

  integer                         :: itsnap(nsnap), isnap, jsnap, i, j, fid
  character(3)                    :: strnum, strnum2

  !ftsnap = [0., 1./8., 1./4., 3./8., 1./2., 1.]
  ftsnap = [(real(i)/real(nsnap-1), i=0,nsnap-1)]
  
  ! Initialization
  PRINT *,'INITIALISATION...'
  call init_from_surfaces( &
          '../cases/vortex/vortex_demo.opt', &
          options, &
          surf, &
          nsurf, &
          interdata, &
          brep, &
          hypergraph, &
          mesh )
  PRINT *,'OK'

  itsnap = nint(4.0*ftsnap/options%timestep)
  
  time = 0.0
  isnap = 0
  do instant = 0,ceiling(options%timespan/options%timestep)
     ! Snapshot
     do jsnap = isnap+1,size(ftsnap)
        if (instant == itsnap(jsnap)) then
           isnap = jsnap
           print *,'isnap =',isnap,', time =',time
           !exit
           do isurf = 1,nsurf
              surf(isurf)%x%degr(1) = size(surf(isurf)%x%coef,1) - 1
              surf(isurf)%x%degr(2) = size(surf(isurf)%x%coef,2) - 1
           end do
           ! refine curves
           call update_intersection_curves(interdata)
           do icurv = 1,interdata%nc
              if ( interdata%curves(icurv)%smooth ) then
                 call refine_intersection_polyline( &
                      interdata%curves(icurv), &
                      1.d0, &
                      options%chord_err )
              end if
           end do

           write (strnum,'(i3.3)') isnap
           ! write surfaces
           do isurf = 1,nsurf
              write (strnum2, '(i3.3)') isurf
              call write_polynomial( &
                   surf(isurf)%x, &
                   'demo_vortex/'//pref//'snap_' // strnum // 'c_' // strnum2 // '.cheb' )
           end do

           
           ! clear mesh
           call free_mesh(mesh)
           
           IF ( .FALSE. ) THEN
              ! Make hypergraph
              call make_hypergraph( &
                   brep, &
                   hypergraph, &
                   feat_edge, &
                   feat_vert )

              !do isurf = 1,nsurf
              !   call economize2( &
              !        surf(isurf)%x, &
              !        1.d-7 )
              !end do

              ! generate BRep conforming mesh
              call generate_brep_conforming_mesh( &
                   brep, &
                   options%hmin, &
                   options%hmax, &
                   options%chord_err, &
                   feat_edge(1:brep%ne), &
                   feat_vert(1:brep%nv), &
                   hypergraph%hyperedges(1:hypergraph%nhe), &
                   hypergraph%nhe, &
                   mesh )

              ! write mesh
              call write_obj_mesh( &
                   mesh, &
                   'demo_vortex/'//pref//'snap_'//strnum//'.obj' )
           END IF
        
           ! write BRep vertices and edges
           call get_free_unit(fid)
           open(unit=fid, file='demo_vortex/'//pref//'verts_'//strnum//'.dat', action='write')
           do i = 1,brep%nv
              write (fid,*) brep%verts(i)%point%xyz
           end do
           close(fid)

           open(unit=fid, file='demo_vortex/'//pref//'edges_'//strnum//'.dat', action='write')
           write (fid,*) brep%ne
           do i = 1,brep%ne
              call get_polyline_endpoints( &
                   brep, &
                   [i,1], &
                   ihead, &
                   itail, &
                   sens, &
                   np )
               write (fid,*) np
              do j = min(ihead,itail),max(ihead,itail)
                 write (fid,*) brep%edges(i)%curve%polyline%xyz(1:3,j)
              end do
           end do
           close(fid)
           
           exit
        end if
     end do

     if ( isnap >= size(ftsnap) ) return
     
     ! Propagate surfaces and advance one timestep
     do isurf = 1,nsurf
        call propagation_step_RK4( &
             surf(isurf), &
             2, &
             options%timestep, &
             time )
     end do
     time = time + options%timestep
     
  end do
  
end program demo_vortex
