program demo_intersection

  use mod_util
  use mod_math
  use mod_polynomial
  use mod_diffgeom
  use mod_types_intersection
  use mod_intersection

  implicit none

  call demo_simple()

  
contains

  subroutine demo_simple()
    implicit none
    logical, parameter                     :: ECONOMIZE = .true.
    type(type_surface)                     :: surf(2)
    type(type_intersection_data), target   :: interdata
    integer                                :: isurf
    character                              :: strnum
    integer, allocatable                   :: valence_out_in(:,:)
    integer                                :: icurv, isplit, ipoint
    type(type_intersection_curve), pointer :: curve => null()
    integer                                :: ncat, fxyz, fuv1, fuv2

    ! Read surfaces 
    do isurf = 1,2
       write (strnum,'(i1)') isurf

       call read_polynomial( &
          surf(isurf)%x, &
          'demo_intersection/simple/c_' // strnum // '.cheb', &
          nvar=2, &
          base=1 )

       call compute_deriv1(surf(isurf))
       call compute_deriv2(surf(isurf))
       call compute_pseudonormal(surf(isurf))
       if ( ECONOMIZE ) call economize2(surf(isurf)%pn, EPSmath)
    end do

    ! intersect
    call intersect_all_surfaces( &
       surf, &
       2, &
       interdata, &
       [.true., .true.], &
       5.d-3, &
       1.d-3, &
       1.d-2 )

    ! export intersection data
    call write_intersection_data( &
       interdata, &
       'demo_intersection/simple/intersection_points.dat', &
       'demo_intersection/simple/intersection_curves.dat' )

    ! concatenate curves
    allocate(valence_out_in(2,interdata%np))
    valence_out_in(1:2,1:interdata%np) = 0
    do icurv = 1,interdata%nc
       do isplit = 1,2
          ipoint = interdata%curves(icurv)%isplit(1,isplit)
          valence_out_in(isplit,ipoint) = valence_out_in(isplit,ipoint) + 1
       end do
    end do

    call print_mat(transpose(valence_out_in))

    ! get starting curve
    do icurv = 1,interdata%nc
       ipoint = interdata%curves(icurv)%isplit(1,1)
       if ( valence_out_in(1,ipoint) == 1 .and. &
            valence_out_in(2,ipoint) == 0 ) then
          print *,icurv
          curve => interdata%curves(icurv)
          exit
       end if
    end do

    ncat = 0
    call get_free_unit(fxyz)
    open(unit=fxyz, file='demo_intersection/simple/curve_xyz.dat', action='write')
    call get_free_unit(fuv1)
    open(unit=fuv1, file='demo_intersection/simple/curve_uv1.dat', action='write')
    call get_free_unit(fuv2)
    open(unit=fuv2, file='demo_intersection/simple/curve_uv2.dat', action='write')
    do 
       ncat = ncat + 1
       
       do ipoint = curve%isplit(2,1),curve%isplit(2,2)-1
          write (fxyz,*) curve%polyline%xyz(1:3,ipoint)
          write (fuv1,*) curve%polyline%uv(1:2,1,ipoint)
          write (fuv2,*) curve%polyline%uv(1:2,2,ipoint)
       end do
       ipoint = curve%isplit(1,2)
       
       ! find next curve
       do icurv = 1,interdata%nc
          if ( interdata%curves(icurv)%isplit(1,1) == ipoint ) then
             print *,icurv
             curve => interdata%curves(icurv)
             exit
          end if
       end do
       
       if ( ncat >= interdata%nc ) exit
    end do
    
    ! final point
    ipoint = curve%isplit(2,2)
    write (fxyz,*) curve%polyline%xyz(1:3,ipoint)
    write (fuv1,*) curve%polyline%uv(1:2,1,ipoint)
    write (fuv2,*) curve%polyline%uv(1:2,2,ipoint)
    
    close(fxyz); close(fuv1); close(fuv2)
    
  end subroutine demo_simple
  
end program demo_intersection
