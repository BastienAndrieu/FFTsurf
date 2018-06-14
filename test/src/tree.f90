program tree
  
  use mod_util
  use mod_obb

  implicit none

  integer,       parameter          :: fp = SELECTED_REAL_KIND(15,307)!kind(1.d0)
  real(kind=fp), parameter          :: EPSregion = real( 1.e-6, kind=fp )



  type type_polynomial
     integer                        :: nvar = 0
     integer                        :: base = 0
     integer                        :: degr(2) = -1
     integer                        :: dim = 0
     real(kind=fp), allocatable     :: coef(:,:,:)
  end type type_polynomial



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



! ########################################################################################



  integer                    :: fileunit, dim, npts
  real(kind=fp), allocatable :: points(:,:), uvboxroot(:)
  type(type_region)          :: root
  integer*8                  :: tic, toc, countrate
  integer                    :: ipt, idim

  call get_free_unit( fileunit )

  open( &
       unit=fileunit, &
       file='tree/points.dat', &
       action='read' )
  read (fileunit,*) npts, dim
  allocate( points(npts,dim), uvboxroot(2*dim) )
  do ipt = 1,npts
     read (fileunit,*) points(ipt,:)
  end do
  close( fileunit )


  do idim = 1,dim
     uvboxroot(2*idim-1) = minval( points(:,idim) )
     uvboxroot(2*idim)   = maxval( points(:,idim) )
  end do
  PRINT *,'UVBOXROOT =',UVBOXROOT

  call init_region( &
       root, &
       dim, &
       uvboxroot )
  deallocate( uvboxroot )
  
  root%npts = npts
  allocate( root%ipts(npts) )
  root%ipts = [ (ipt,ipt=1,npts) ]

  call system_clock( tic, countrate )
  call test_tree( &
       points, &
       npts, &
       dim, &
       root )
  call system_clock( toc )
  print *,'ELAPSED =',real(toc - tic) / real(countrate)
  

  call export_region_tree( &
       root, &
       'tree/tree.dat' )

  call free_region_tree( root )
  deallocate( points )


contains



  recursive subroutine test_tree( &
       points, &
       npts, &
       dim, &
       region )
    use mod_util
    implicit none
    integer,           intent(in)    :: npts, dim
    real(kind=fp),     intent(in)    :: points(npts,dim)
    type(type_region), intent(inout) :: region
    real(kind=fp), allocatable       :: x(:)
    real(kind=fp)                    :: uvm(dim)
    integer                          :: stat
    integer                          :: jpt, ipt, idim, ichild

    !RETURN
    if ( associated(region%parent) ) then
       outer : do jpt = 1,region%parent%npts
          ipt = region%parent%ipts(jpt)
          do idim = 1,dim
             if ( points(ipt,idim) < region%uvbox(2*idim-1) .or. &
                  points(ipt,idim) > region%uvbox(2*idim) ) cycle outer
          end do
          call append( &
               region%ipts, &
               ipt, &
               noduplicates=.true., &
               newlength=region%npts )
       end do outer
    end if

    if ( region%npts < 2 ) return

    allocate( x(region%npts) )
    IF (.true.) THEN
       ! pick median point along each dimension
       do idim = 1,dim
          x = points( region%ipts(1:region%npts), idim )
          call bubblesort_double( x )
          if ( mod(region%npts,2) == 0 ) then
             uvm(idim) = 0.5_fp * sum( x(region%npts/2+[0,1]) )
          else
             uvm(idim) =  x(region%npts/2+1)
          end if
       end do
    ELSE
       uvm = 0.5_fp * ( region%uvbox(1::2) + region%uvbox(2::2) )
    END IF
    deallocate( x )

    
    call subdiv_region( &
         region, &
         uvm, &
         stat )

    if ( stat >= 0 .and. stat < dim ) then
       do ichild = 1,size(region%child)
          call test_tree( &
               points, &
               npts, &
               dim, &
               region%child(ichild) )
       end do
    end if

  end subroutine test_tree























  subroutine free_polynomial( poly )
    implicit none
    type(type_polynomial), intent(inout) :: poly
    if ( allocated(poly%coef) ) deallocate( poly%coef )
  end subroutine free_polynomial









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











  subroutine print_region( region )
    implicit none
    type(type_region), intent(in) :: region
    
    print *,'-------------------------------------'
    print *,'dim =',region%dim
    if ( allocated( region%uvbox ) ) then
       print *,'uvbox =',region%uvbox
    else
       print *,'uvbox = N/A'
    end if
    if ( associated( region%child ) ) then
       print *,size(region%child),' children'
    else
       print *,0,' children'
    end if
    print *,'-------------------------------------'

  end subroutine print_region

  
end program tree
