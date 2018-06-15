module mod_regiontree

  use mod_math
  use mod_obb
  use mod_polynomial

  real(kind=fp), parameter          :: EPSregion = real( 1.e-6, kind=fp ) ! minimal range along any dimension

  type type_region
     integer                        :: dim = 0
     real(kind=fp), allocatable     :: uvbox(:) ! [ x_1_min, x_1_max, ..., x_dim_min, x_dim_max ]
     type(type_obb), pointer        :: xyzbox => null()
     type(type_polynomial), pointer :: poly => null()
     integer                        :: npts = 0
     integer, allocatable           :: ipts(:)
     type(type_region), pointer     :: parent => null()
     type(type_region), pointer     :: child(:) => null()
  end type type_region

  type ptr_region
     type(type_region), pointer     :: ptr => null()
  end type ptr_region
  


contains

  subroutine init_region( &
       region, &
       dim, &
       uvbox ) 
    implicit none
    type(type_region), intent(inout) :: region
    integer,           intent(in)    :: dim
    real(kind=fp),     intent(in)    :: uvbox(2*dim)

    if ( allocated(region%uvbox) ) then
       if ( size(region%uvbox) < 2*dim ) deallocate( region%uvbox )
    end if
    region%dim = dim
    if ( .not.allocated(region%uvbox) ) allocate( region%uvbox(2*dim) )
    region%uvbox(1:2*dim) = uvbox

    if ( associated(region%xyzbox) ) deallocate( region%xyzbox )
    if ( associated(region%child) )  deallocate( region%child )
    if ( allocated(region%ipts) )    deallocate( region%ipts )
    region%npts = 0

    if ( associated(region%poly) ) then
       call free_polynomial( region%poly )
       deallocate( region%poly )
    end if
    nullify( region%parent )

  end subroutine init_region





  subroutine subdiv_region( &
       region, &
       uvs, &
       stat )
    ! stat = -1  if region already has children
    !      = 0   if regular subdivision (2**dim children)
    !      = n>0 if n degenerate dimensions (2**(dim-n) children)
    !            except if n=dim (the region is not subdivided)
    implicit none
    type(type_region), intent(inout), target :: region
    real(kind=fp),     intent(in)            :: uvs(region%dim)
    integer,           intent(out)           :: stat
    real(kind=fp)                            :: uvbox(2*region%dim)
    logical                                  :: degenerate(region%dim)
    integer                                  :: idim, nchild, jdim, ichild

    if ( associated(region%child) ) then
       ! the region already has children
       stat = -1
       return
    end if


    do idim = 1,region%dim
       ! check for possibly "degenerate" dimensions
       degenerate(idim) = ( &
            abs( uvs(idim) - region%uvbox(2*idim-1) ) < EPSregion .or. &
            abs( uvs(idim) - region%uvbox(2*idim) ) < EPSregion )
    end do

    stat = count(degenerate)
    if ( stat >= region%dim ) return ! the subdivision point is located at one of the region's corners

    nchild = 2**( region%dim - stat ) ! number of children
    allocate( region%child(nchild) )

    ! get uvboxes of the children
    do ichild = 1,nchild ! loop over the children
       jdim = 0 ! reset the counter of non-degenerate dimensions
       do idim = 1,region%dim ! loop over all the dimensions
          if ( degenerate(idim) ) then ! degenerate dimension   ...
             uvbox(2*idim-[1,0]) = region%uvbox(2*idim-[1,0]) ! ...same bounds as the parent
          else
             ! non-degenerate dimension
             jdim = jdim + 1 ! increment the counter of non-degenerate dimensions
             if ( mod( (ichild-1)/(2**(jdim-1)), 2 ) == 0 ) then
                uvbox(2*idim-[1,0]) = [ region%uvbox(2*idim-1), uvs(idim) ]
             else
                uvbox(2*idim-[1,0]) = [ uvs(idim), region%uvbox(2*idim) ]
             end if
          end if
       end do

       call init_region( &
            region%child(ichild), &
            region%dim, &
            uvbox ) 
       region%child(ichild)%parent => region
    end do

  end subroutine subdiv_region




  recursive subroutine free_region_tree( &
       region )
    implicit none
    type(type_region), intent(inout) :: region
    integer                          :: ichild

    if ( associated(region%child) ) then
       do ichild = 1,size(region%child)
          call free_region_tree( region%child(ichild) )
       end do
       deallocate( region%child )
    end if

    if ( allocated(region%uvbox) )   deallocate( region%uvbox )
    if ( associated(region%xyzbox) ) deallocate( region%xyzbox )
    if ( allocated(region%ipts) )    deallocate( region%ipts )
    if ( associated(region%parent) ) then
       if ( associated(region%poly) ) then
          call free_polynomial( region%poly )
          deallocate( region%poly )
       end if
    end if

  end subroutine free_region_tree




  subroutine export_region_tree( &
       region, &
       filename )
    use mod_util
    implicit none
    type(type_region), intent(in) :: region
    character(*),      intent(in) :: filename
    integer                       :: file_unit

    call get_free_unit( file_unit )

    open( &
         unit = file_unit, &
         file = trim(filename), &
         action = 'write' )

    call write_region( &
         region, &
         file_unit )

    close( file_unit )

    PRINT *,'region tree written in ',trim(filename)

  end subroutine export_region_tree





  recursive subroutine write_region( &
       region, &
       file_unit )
    ! Writes the data contained in a (leaf) region.
    implicit none
    type(type_region), intent(in) :: region
    integer,           intent(in) :: file_unit
    integer                       :: ichild

    if ( associated(region%child) ) then
       ! the current region has child regions, carry on the recursion
       do ichild = 1,size(region%child)
          call write_region( &
               region%child(ichild), &
               file_unit )
       end do
    else
       ! the current region is a leaf, write its data
       write (file_unit,*) region%uvbox
    end if

  end subroutine write_region






  
  recursive subroutine add_point_bottom_up( &
       region, &
       ipt )
    use mod_util
    implicit none
    type(type_region), intent(inout) :: region
    integer,           intent(in)    :: ipt

    call append_list( &
         region%ipts, &
         region%npts, &
         ipt )
  end subroutine add_point_bottom_up


end module mod_regiontree
