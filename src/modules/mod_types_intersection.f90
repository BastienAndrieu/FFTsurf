module mod_types_intersection

  use mod_math
  use mod_polynomial
  use mod_diffgeom2
  use mod_obb

  implicit none

  type type_point_on_surface
     type(type_surface), pointer                :: surf => null()
     real(kind=fp)                              :: uv(2)
     type(type_point_on_surface), pointer       :: next => null()
  end type type_point_on_surface



  type type_intersection_point
     real(kind=fp)                              :: xyz(3)
     type(type_point_on_surface), pointer       :: pos => null()
     integer                                    :: npos = 0
  end type type_intersection_point



  type type_intersection_segment
     integer                                    :: endpoints(2)
     type(type_intersection_segment), pointer   :: parent => null()
     type(type_intersection_segment), pointer   :: child(:) => null()
  end type type_intersection_segment



  type type_intersection_polyline
     integer                                    :: np = 0    ! actual nb. of poyline points
     real(kind=fp), allocatable                 :: s(:)      ! approximate curvilinear abscissa
     real(kind=fp), allocatable                 :: uv(:,:,:) ! uv coords.  (u/v, #surface, #point)
     real(kind=fp), allocatable                 :: xyz(:,:)  ! xyz coords. (x/y/z, #point)
  end type type_intersection_polyline



  type type_intersection_curve
     logical                                    :: smooth = .false.
     type(ptr_surface)                          :: surf(2)
     real(kind=fp)                              :: uvbox(2,2,2) ! umin/max, vmin/max, #surf
     real(kind=fp)                              :: param_vector(3)
     real(kind=fp)                              :: w0
     integer, allocatable                       :: isplit(:,:) ! row 1: #intersection point, row 2: #polyline point
     integer                                    :: nsplit = 0
     type(type_intersection_segment), pointer   :: root => null()
     type(type_intersection_polyline), pointer  :: polyline => null()
  end type type_intersection_curve



  type type_intersection_data
     integer                                    :: np = 0
     type(type_intersection_point), allocatable :: points(:)
     integer                                    :: nc = 0
     type(type_intersection_curve), allocatable :: curves(:)
  end type type_intersection_data



  
  type ptr_intersection_curve
     type(type_intersection_curve), pointer     :: ptr => null()
  end type ptr_intersection_curve



contains


  subroutine free_intersection_data( &
      interdata )
    implicit none
    type(type_intersection_data), intent(inout) :: interdata
    type(type_point_on_surface), pointer        :: current, next
    integer                                     :: ip, ic
    
    if ( allocated(interdata%points) ) then
       do ip = 1,interdata%np
          interdata%points(ip)%npos = 0
          current => interdata%points(ip)%pos
          do while( associated(current) )
             next => current%next
             deallocate( current )
             nullify( current )
             current => next
          end do
       end do
       deallocate( interdata%points )
    end if
    interdata%np = 0

    if ( allocated(interdata%curves) ) then
       do ic = 1,interdata%nc
          nullify(interdata%curves(ic)%surf(1)%ptr, interdata%curves(ic)%surf(2)%ptr)

          if ( allocated(interdata%curves(ic)%isplit) ) deallocate(interdata%curves(ic)%isplit)

          if ( associated(interdata%curves(ic)%polyline) ) then
             if ( allocated(interdata%curves(ic)%polyline%s  ) ) deallocate(interdata%curves(ic)%polyline%s  )
             if ( allocated(interdata%curves(ic)%polyline%uv ) ) deallocate(interdata%curves(ic)%polyline%uv )
             if ( allocated(interdata%curves(ic)%polyline%xyz) ) deallocate(interdata%curves(ic)%polyline%xyz)
             interdata%curves(ic)%polyline%np = 0
             deallocate(interdata%curves(ic)%polyline)
          end if
          nullify(interdata%curves(ic)%polyline)

          if ( associated(interdata%curves(ic)%root) ) then
             call free_intersection_segment_tree(interdata%curves(ic)%root)
             deallocate(interdata%curves(ic)%root)
          end if
          nullify(interdata%curves(ic)%root)
       end do
       deallocate( interdata%curves )
    end if
    interdata%nc = 0

  end subroutine free_intersection_data





  recursive subroutine free_intersection_segment_tree( segment )
    implicit none
    type(type_intersection_segment), intent(inout) :: segment
    integer                                        :: ichild

    if ( associated( segment%child ) ) then
       do ichild = 1,size(segment%child)
          call free_intersection_segment_tree( segment%child(ichild) )
       end do
       deallocate( segment%child )
    end if

  end subroutine free_intersection_segment_tree




end module mod_types_intersection
