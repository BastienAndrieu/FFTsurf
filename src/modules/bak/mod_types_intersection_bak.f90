module mod_types_intersection
  
  use mod_math
  !use mod_list
  use mod_diffgeom
  use mod_obb

  implicit none

  ! PARAMETERS ==================================================================================================================
  integer, parameter                            :: param_init_interpt_n = 10                                                    !
  integer, parameter                            :: param_xtra_interpt_n = 10                                                    !
  integer, parameter                            :: param_init_polyline_np = 100                                                 !
  integer, parameter                            :: param_xtra_polyline_np = 20                                                  !
  
  real(kind=MATHpr), parameter                  :: EPSregion = real( 1.e-6, kind=MATHpr )                                       !
  ! =============================================================================================================================


  ! DERIVED TYPES ===============================================================================================================
  ! --------------------------------------------------------------------------------------------------------------------------- !
  ! Rectangular region of a parametric surface (contains no actual reference to a surface)                                      !
  type type_surface_region                                                                                                      !
     real(kind=MATHpr)                          :: uvbox(2,2)          ! [umin, vmin ; umax, vmax]                              !
     type(type_obb),            pointer         :: xyzbox => null()    ! bounding box in 3D-space                               !
     type(type_obb),            pointer         :: pnbox => null()     ! bounding box of pseudo-normal patch                    !
     type(type_surface_region), pointer         :: child(:) => null()  ! child regions in a tree                                !
     type(type_surface_region), pointer         :: parent => null()    ! parent region in a tree                                !
     !integer, allocatable                       :: ipoints(:)          ! list of inter_pts indices                              !
     !integer                                    :: npoints = 0         ! actual nb. of inter_pts                                !
     integer, allocatable                       :: ipts_ss(:)          ! indices of surface-surface intersection points         !
     integer, allocatable                       :: ipts_cs(:)          ! indices of curve-surface intersection points           !
     !type(type_list)                            :: ipts_cs             ! linked-list of indices of curve-surface intersect. pts !
     !integer                                    :: npts_cs = 0         ! actual nb. of curve-surface intersection points        !
  end type type_surface_region                                                                                                  !
  ! --------------------------------------------------------------------------------------------------------------------------- ! 


  ! --------------------------------------------------------------------------------------------------------------------------- !
  ! Workaround for making arrays of pointers to surface_regions                                                                 !
  type ptr_surface_region                                                                                                       !
     type(type_surface_region), pointer         :: ptr => null()                                                                !
  end type ptr_surface_region                                                                                                   !
  ! --------------------------------------------------------------------------------------------------------------------------- ! 


  ! --------------------------------------------------------------------------------------------------------------------------- !
  ! Region (segment) of a parametric curve (contains no actual reference to a curve)                                            !
  type type_curve_region                                                                                                        !
     real(kind=MATHpr)                          :: tbox(2)             ! [tmin, tmax]                                           !
     type(type_obb),          pointer           :: xyzbox => null()    ! bounding box in 3D-space                               !
     type(type_curve_region), pointer           :: child(:) => null()  ! child regions in a tree                                !
     integer, allocatable                       :: ipts_cs(:)          ! list of indices of curve-surface intersection points   !
     !type(type_list)                            :: ipts_cs             ! linked-list of indices of curve-surface intersect. pts !
     !integer                                    :: npts_cs = 0         ! actual nb. of curve-surface intersection points        !
  end type type_curve_region                                                                                                    !
  ! --------------------------------------------------------------------------------------------------------------------------- ! 



  ! --------------------------------------------------------------------------------------------------------------------------- !
  type type_point_on_surface                                                                                                    !
     type(type_parametric_surface), pointer     :: surface => null() ! surface the points lies on                               !
     real(kind=MATHpr)                          :: uv(2)             ! uv coordinates of the point                              !
     real(kind=MATHpr)                          :: uvtol             ! precision of the uv coords.                              !
     real(kind=MATHpr)                          :: xyz(3)            ! xyz coords. of the point                                 !
     type(type_point_on_surface), pointer       :: next => null()    ! next element in linked list (see intersection_point)     !
  end type type_point_on_surface                                                                                                !
  ! --------------------------------------------------------------------------------------------------------------------------- !
  

  ! --------------------------------------------------------------------------------------------------------------------------- !
  type type_intersection_point                                                                                                  !
     type(type_point_on_surface), pointer       :: head => null() ! head of linked-list of points_on_surface                    !
     real(kind=MATHpr)                          :: xyz(3)         ! xyz coords. of the point                                    !
  end type type_intersection_point                                                                                              !
  ! --------------------------------------------------------------------------------------------------------------------------- !


  ! --------------------------------------------------------------------------------------------------------------------------- !
  type type_intersection_polyline                                                                                               !
     integer                                    :: np = 0    ! actual nb. of poyline points                                     !
     real(kind=MATHpr), allocatable             :: s(:)      ! approximate curvilinear abscissa                                 !
     real(kind=MATHpr), allocatable             :: uv(:,:,:) ! uv coords.  (u/v, #surface, #point)                              !
     real(kind=MATHpr), allocatable             :: xyz(:,:)  ! xyz coords. (x/y/z, #point)                                      !
  end type type_intersection_polyline                                                                                           !
  ! --------------------------------------------------------------------------------------------------------------------------- !


  ! --------------------------------------------------------------------------------------------------------------------------- !
  type type_intersection_segment                                                                                                !
     integer                                    :: head = 0            ! first polyline point                                   !
     integer                                    :: tail = 0            ! last polyline point                                    !
     integer                                    :: endpoints(2)        ! indices of intersection_points at endpoints            !
     real(kind=MATHpr)                          :: endangles(2)        ! tangent direction at endpoints                         !
     type(type_intersection_segment), pointer   :: child(:) => null()  ! child segments in tree                                 !
     type(type_intersection_segment), pointer   :: parent(:) => null() ! parent segment in tree                                 !
  end type type_intersection_segment                                                                                            !
  ! --------------------------------------------------------------------------------------------------------------------------- !


  ! --------------------------------------------------------------------------------------------------------------------------- !
  ! When intersected with other curves, the curve is subdivided into segments that are stored into a tree data structure        !
  type type_intersection_curve                                                                                                  !
     type(ptr_parametric_surface)               :: surfaces(2)      ! pointers to intersected surfaces                          !
     type(ptr_surface_region)                   :: regions(2)       ! pointers to surface_regions                               !
     logical                                    :: smooth = .false. ! tangential intersection flag                              !
     type(type_intersection_polyline)           :: polyline         ! piecewise linear approximation                            !
     type(type_intersection_segment)            :: root             ! root segment in tree                                      !
     real(kind=MATHpr)                          :: param_vector(3)  ! parameterization vector                                   !
     type(type_intersection_curve), pointer     :: next => null()   ! next element in linked list                               !
  end type type_intersection_curve                                                                                              !
  ! --------------------------------------------------------------------------------------------------------------------------- !


  ! --------------------------------------------------------------------------------------------------------------------------- !
  ! Structure collecting all the intersection points and curves                                                                 !
  type type_intersection_data                                                                                                   !
     integer                                    :: np = 0    ! actual nb. of intersection_points                                !
     integer                                    :: nc = 0    ! actual nb. of intersection_curves                                !
     type(type_intersection_point), allocatable :: points(:) ! list of intersection_points                                      !
     type(type_intersection_curve), allocatable :: curves(:) ! list of intersection_curves                                      !
  end type type_intersection_data                                                                                               !
  ! --------------------------------------------------------------------------------------------------------------------------- !  
  ! =============================================================================================================================

end module mod_types_intersection
