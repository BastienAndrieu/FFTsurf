program LL_patch

  use mod_math
  use mod_diffgeom
  use mod_geometry
  use mod_chebyshev
  use mod_polynomial
  use mod_eos
  use mod_util

  implicit none
  character(99)                              :: argstr
  integer, parameter                         :: degr_default = 16
  integer                                    :: nplanes
  real(kind=fp), allocatable, dimension(:,:) :: cc_center, normal, corner
  real(kind=fp)                              :: vxyz(3)
  integer                                    :: degr
  type(type_surface)                         :: surf
  integer                                    :: i, fid

  ! READ (OPTIONAL) ARGUMENT
  if ( command_argument_count() < 1 ) then
     degr = degr_default
  else
     call get_command_argument(1, argstr)
     read (argstr, *) degr
  end if

  ! READ INPUT DATA
  call get_free_unit(fid)
  open(unit=fid, file='LL_patch/input.dat', action='read')
  read (fid, *) nplanes
  allocate(cc_center(3,nplanes), normal(3,nplanes), corner(3,nplanes))
  read (fid, *) vxyz(1:3)
  do i = 1,nplanes
     read (fid,*) cc_center(1:3,i)
  end do
  do i = 1,nplanes
     read (fid,*) normal(1:3,i)
  end do
  do i = 1,nplanes
     read (fid,*) corner(1:3,i)
  end do

  ! RUN 
  call LL_patch_from_arcs( &
       nplanes, &
       cc_center, &
       normal, &
       corner, &
       vxyz, &
       degr, &
       surf )

  ! EXPORT
  call write_polynomial(surf%x, 'LL_patch/surf.cheb')
  call system('cp ../debug/LL_patch_from_arcs/B.dat LL_patch/B.dat')
  
end program LL_patch
