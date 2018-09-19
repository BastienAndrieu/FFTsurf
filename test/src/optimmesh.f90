program optimmesh

  use mod_util
  use mod_mesh
  use mod_optimmesh

  implicit none

  integer, parameter         :: nv0 = 10000, nt0 = 2*nv0
  
  integer                    :: fid, err, narg
  character(99)              :: arg
  real(kind=fp)              :: ratio

  integer                    :: ip, iv
  integer                    :: tri(3)
  integer, allocatable       :: tritmp(:,:)
  real(kind=fp)              :: xyz(3)
  real(kind=fp), allocatable :: xyztmp(:,:)
    
  type(type_surface_mesh)    :: mesh
  real(kind=fp), allocatable :: rad(:)

  integer                    :: npass, ipass
  character(2)               :: strnum


  narg = command_argument_count()
  if ( narg < 1 ) then
     npass = 1
  else
     call get_command_argument(1, arg)
     read (arg,*) npass
  end if
  if ( narg < 2 ) then
     ratio = 1._fp
  else
     call get_command_argument(2, arg)
     read (arg,*) ratio
  end if

  ! read mesh from files
  call get_free_unit(fid)

  allocate(mesh%tri(3,nt0))
  PRINT *,'READ TRIANGLES...'
  !open(unit=fid, file='Jouke/meshgen/brepmesh/tri.dat', action='read')
  open(unit=fid, file='Jouke/meshgen/brepmesh/tri_ec.dat', action='read')
  do
     read (fid,*,iostat=err) tri
     if ( err /=0 ) exit
     if ( mesh%nt + 1 > size(mesh%tri,2) ) then
        call move_alloc(from=mesh%tri, to=tritmp)
        allocate(mesh%tri(3,mesh%nt + 1000))
        mesh%tri(1:3,1:mesh%nt) = tritmp(1:3,1:mesh%nt)
        deallocate(tritmp)
     end if
     mesh%nt = mesh%nt + 1
     mesh%tri(1:3,mesh%nt) = tri(1:3)
  end do
  close(fid)
  PRINT *,'...OK'

  allocate(mesh%xyz(3,nv0))
  PRINT *,'READ XYZ...'
  !open(unit=fid, file='Jouke/meshgen/brepmesh/xyz.dat', action='read')
  !open(unit=fid, file='Jouke/meshgen/brepmesh/xyz_ec.dat', action='read')
  open(unit=fid, file='Jouke/meshgen/brepmesh/xyz_smooth.dat', action='read')
  do
     read (fid,*,iostat=err) xyz
     if ( err /=0 ) exit
     if ( mesh%nv + 1 > size(mesh%xyz,2) ) then
        call move_alloc(from=mesh%xyz, to=xyztmp)
        allocate(mesh%xyz(3,mesh%nv + 1000))
        mesh%xyz(1:3,1:mesh%nv) = xyztmp(1:3,1:mesh%nv)
        deallocate(xyztmp)
     end if
     mesh%nv = mesh%nv + 1
     mesh%xyz(1:3,mesh%nv) = xyz(1:3)
  end do
  close(fid)
  PRINT *,'...OK'

  PRINT *,'READ PATHS...'
  open(unit=fid, file='Jouke/meshgen/brepmesh/paths_ec.dat', action='read')
  read (fid,*) mesh%npaths
  allocate(mesh%paths(mesh%npaths))
  do ip = 1,mesh%npaths
     read (fid,*) mesh%paths(ip)%hyperedge
     read (fid,*) mesh%paths(ip)%nv
     allocate(mesh%paths(ip)%verts(mesh%paths(ip)%nv))
     do iv = 1,mesh%paths(ip)%nv
        read (fid,*) mesh%paths(ip)%verts(iv)
     end do
  end do
  close(fid)
  PRINT *,'...OK'

  
  
  ! make halfedges
  PRINT *,'MAKE HALFEDGES...'
  call make_halfedges(mesh)
  PRINT *,'...OK'

  ! compute discrete minimum curvature radius
  allocate(rad(mesh%nv))
  PRINT *,'COMPUTE CURVATURE RADIUS...'
  call discrete_minimum_curvature_radius( &
       mesh, &
       rad )
  PRINT *,'...OK'
  !rad = min(1.d0, max(0._fp, rad))
  print *, minval(rad), minloc(rad,1)
  rad = rad/minval(rad)
  rad = min(ratio, rad)
  rad = (rad - minval(rad)) / (maxval(rad) - minval(rad))

  do ipass = 0,npass
     write (strnum,'(I2.2)') ipass 
     call write_tecplot_mesh_solv( &
          mesh, &
          'optimmesh/pass_'//strnum//'.dat', &
          'pass_'//strnum, &
          rad )
     if ( ipass == npass ) exit
     PRINT *,'PASS',IPASS+1,'/',NPASS
     call smooth_function_mesh( &
          mesh, &
          rad, &
          1 )
  end do


  
end program optimmesh
