subroutine inherit_points( &
     region, &
     coords, &
     npts )
  use mod_math
  use mod_regiontree
  implicit none
  type(type_region), intent(inout) :: region
  integer,           intent(in)    :: npts
  real(kind=fp),     intent(in)    :: coords(region%dim,npts)
  logical, allocatable             :: mask(:)
  integer                          :: idim, ipt, jpt

  if ( .not.associated(region%parent) ) return
  if ( region%parent%npts < 1 ) return

  allocate( mask(region%parent%npts) )
  mask(:) = .true.
  outer : do jpt = 1,region%parent%npts
     ipt = region%parent%ipts(jpt)
     do idim = 1,region%dim
        if ( .not.is_in_interval( &
             coords(idim,ipt), &
             region%uvbox(2*idim-1), &
             region%uvbox(2*idim), &
             tolerance=EPSregion) ) then
           mask(jpt) = .false.
           cycle outer
        end if
     end do
  end do outer

  call append_n( &
       region%ipts, &
       region%npts, &
       pack( region%parent%ipts(1:region%parent%npts), mask ), &
       count(mask), &
       unique=.true. )

  deallocate( mask )

end subroutine inherit_points
