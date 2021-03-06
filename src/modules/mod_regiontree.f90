module mod_regiontree

  use mod_math
  use mod_obb
  use mod_polynomial

  real(kind=fp), parameter          :: EPSregion = real( 1.e-7, kind=fp ) ! minimal range along any dimension

  type type_region
     integer                        :: dim = 0
     real(kind=fp), allocatable     :: uvbox(:) ! [ x_1_min, x_1_max, ..., x_dim_min, x_dim_max ]
     type(type_obb), pointer        :: xyzbox => null()
     !type(type_polynomial), pointer :: poly => null()
     type(ptr_polynomial), allocatable :: poly(:)
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
    integer                          :: i

    region%dim = dim
    if ( allocated(region%uvbox) ) deallocate( region%uvbox )
    if ( .not.allocated(region%uvbox) ) allocate( region%uvbox(2*dim) )
    region%uvbox(1:2*dim) = uvbox

    if ( associated(region%xyzbox) ) deallocate( region%xyzbox )
    if ( associated(region%child) )  deallocate( region%child )
    if ( allocated(region%ipts) )    deallocate( region%ipts )
    region%npts = 0

    if ( allocated(region%poly) ) then
       do i = 1,size(region%poly)
          if ( associated(region%poly(i)%ptr) ) then
             call free_polynomial( region%poly(i)%ptr )
             deallocate( region%poly(i)%ptr )
          end if
       end do
    end if

    nullify( region%parent )

  end subroutine init_region





  subroutine subdiv_region( &
       region, &
       uvs, &
       stat, &
       EPSopt )
    use mod_tolerances
    ! stat = -1  if region already has children
    !      = 0   if regular subdivision (2**dim children)
    !      = n>0 if n degenerate dimensions (2**(dim-n) children)
    !            except if n=dim (the region is not subdivided)
    implicit none
    type(type_region), intent(inout), target :: region
    real(kind=fp),     intent(in)            :: uvs(region%dim)
    integer,           intent(out)           :: stat
    real(kind=fp),     intent(in), optional  :: EPSopt
    real(kind=fp)                            :: EPS
    real(kind=fp)                            :: uvbox(2*region%dim)
    logical                                  :: degenerate(region%dim)
    integer                                  :: stat_alloc
    integer                                  :: idim, nchild, jdim, ichild

    if ( present(EPSopt) ) then
       EPS = EPSopt
    else
       EPS = EPSregion
    end if
    
    if ( associated(region%child) ) then
       ! the region already has children
       stat = -1
       return
    end if

    do idim = 1,region%dim
       ! check that uvs is not outside the region
       if ( .not.is_in_closed_interval( &
            uvs(idim), &
            region%uvbox(2*idim-1), &
            region%uvbox(2*idim), &
            tolerance=EPS ) ) then
          STOP 'subdiv_region : subdivision point outside the region'
       end if
       
       ! check for possibly "degenerate" dimensions
       degenerate(idim) = ( &
            abs(uvs(idim) - region%uvbox(2*idim-1)) < EPS .or. &
            abs(uvs(idim) - region%uvbox(2*idim)  ) < EPS )
            !abs(uvs(idim) - region%uvbox(2*idim-1)) < EPSuv .or. &
            !abs(uvs(idim) - region%uvbox(2*idim)  ) < EPSuv )
            !abs(uvs(idim) - region%uvbox(2*idim-1)) < EPSregion .or. &
            !abs(uvs(idim) - region%uvbox(2*idim)  ) < EPSregion )
    end do

    !PRINT *,'subdiv_region : degenerate =',degenerate

    stat = count(degenerate)
    if ( stat >= region%dim ) return ! the subdivision point is located at one of the region's corners

    nchild = 2**(region%dim - stat) ! number of children
    allocate(region%child(nchild), stat=stat_alloc)
    if ( stat_alloc /= 0 ) then
       STOP 'subdiv_region : could not allocate region%child'
    end if

    ! get uvboxes of the children
    do ichild = 1,nchild ! loop over the children
       jdim = 0 ! reset the counter of non-degenerate dimensions
       do idim = 1,region%dim ! loop over all the dimensions
          if ( degenerate(idim) ) then ! degenerate dimension   ...
             uvbox(2*idim-[1,0]) = region%uvbox(2*idim-[1,0]) ! ...same bounds as the parent
          else
             ! non-degenerate dimension
             jdim = jdim + 1 ! increment the counter of non-degenerate dimensions
             if ( mod((ichild-1)/(2**(jdim-1)), 2) == 0 ) then
                uvbox(2*idim-[1,0]) = [region%uvbox(2*idim-1), uvs(idim)]
             else
                uvbox(2*idim-[1,0]) = [uvs(idim), region%uvbox(2*idim)]
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



  recursive subroutine free_ipts( &
       region )
    type(type_region), intent(inout) :: region
    integer                          :: ichild

    if ( allocated(region%ipts) ) deallocate(region%ipts)
    region%npts = 0
    
    if ( associated(region%child) ) then
       do ichild = 1,size(region%child)
          call free_ipts( region%child(ichild) )
       end do
    end if

  end subroutine free_ipts
  




  recursive subroutine free_region_tree( &
       region )
    implicit none
    type(type_region), intent(inout) :: region
    integer                          :: ichild, i

    if ( allocated(region%uvbox) )   deallocate( region%uvbox )
    if ( associated(region%xyzbox) ) deallocate( region%xyzbox )
    if ( allocated(region%ipts) )    deallocate( region%ipts )

    !if ( associated(region%parent) ) then
    if ( allocated(region%poly) ) then
       do i = 1,size(region%poly)
          if ( associated(region%poly(i)%ptr) ) then
             call free_polynomial( region%poly(i)%ptr )
             deallocate( region%poly(i)%ptr )
          end if
       end do
       deallocate( region%poly )
    end if
    !end if

    if ( associated(region%child) ) then
       do ichild = 1,size(region%child)
          call free_region_tree( region%child(ichild) )
       end do
       deallocate( region%child )
    end if

  end subroutine free_region_tree



  
  subroutine copy_region( &
       from, &
       to )
    ! Creates a copy of a region (all attributes except 'child' are copied)
    implicit none
    type(type_region), intent(in)    :: from
    type(type_region), intent(inout) :: to
    integer                          :: ipoly

    call init_region( &
         to, &
         from%dim, &
         from%uvbox ) 

    to%xyzbox => from%xyzbox

    if ( allocated(from%poly) ) then
       allocate(to%poly(size(from%poly)))
       do ipoly = 1,size(from%poly)
          to%poly(ipoly)%ptr => from%poly(ipoly)%ptr
       end do
    end if

    to%npts = from%npts
    if ( allocated(from%ipts) ) then
       allocate(to%ipts(size(from%ipts)))
       to%ipts = from%ipts
    end if

    to%parent => from%parent

  end subroutine copy_region






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

    call append_n( &
         region%ipts, &
         region%npts, &
         [ipt], &
         1, &
         unique=.true. )
    
    if ( associated(region%parent) ) then
       call add_point_bottom_up( &
            region%parent, &
            ipt )
    end if
    
  end subroutine add_point_bottom_up




  recursive subroutine add_points_bottom_up( &
       region, &
       ipts, &
       npts )
    use mod_util
    implicit none
    type(type_region), intent(inout) :: region
    integer,           intent(in)    :: npts
    integer,           intent(in)    :: ipts(npts)

    call append_n( &
         region%ipts, &
         region%npts, &
         ipts(1:npts), &
         npts, &
         unique=.true. )

    IF ( .FALSE. ) THEN
       PRINT *,'ADDPTSB->U >>>-----------------------'
       PRINT *,'UBVOX =',REGION%UVBOX
       PRINT *,'IPTS ='
       IF ( REGION%NPTS < 1 ) THEN
          PRINT *,'N/A'
       ELSE
          PRINT *,REGION%IPTS(1:REGION%NPTS)
       END IF
       PRINT *,'----------------------------------<<<'
    END IF

    if ( associated(region%parent) ) then
       call add_points_bottom_up( &
            region%parent, &
            ipts(1:npts), &
            npts )
    end if

  end subroutine add_points_bottom_up






end module mod_regiontree
