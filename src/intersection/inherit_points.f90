subroutine inherit_points( &
     region, &
     coords, &
     npts )
  use mod_util
  use mod_math
  use mod_regiontree
  implicit none
  LOGICAL, PARAMETER :: DEBUG = ( GLOBALDEBUG .AND. .false. )
  type(type_region), intent(inout) :: region
  integer,           intent(in)    :: npts
  real(kind=fp),     intent(in)    :: coords(region%dim,npts)
  logical, allocatable             :: mask(:)
  integer                          :: idim, ipt, jpt

  if ( .not.associated(region%parent) ) return
  if ( region%parent%npts < 1 ) return

  allocate( mask(region%parent%npts) )
  mask(:) = .true.
  outer : do jpt = 1,region%parent%npts ! <-------------+
     ipt = region%parent%ipts(jpt)                      !
     IF ( IPT < 0 .OR. IPT > NPTS ) THEN
        PRINT *,'*** REGION%PARENT%IPTS =',REGION%PARENT%IPTS(1:REGION%PARENT%NPTS)
     END IF
     do idim = 1,region%dim ! <--------------------+    !
        if ( .not.is_in_closed_interval( &         !    !
             coords(idim,ipt), &                   !    !
             region%uvbox(2*idim-1), &             !    !
             region%uvbox(2*idim), &               !    !
             tolerance=EPSregion) ) then ! <--+    !    !
           mask(jpt) = .false.                !    !    !
           cycle outer                        !    !    !
        end if ! <----------------------------+    !    !
     end do ! <------------------------------------+    !
  end do outer ! <--------------------------------------+

  IF ( DEBUG ) THEN
     PRINT *,'   REGION%UVBOX =',REGION%UVBOX
     DO IPT = 1,region%parent%npts
        PRINT *,REGION%PARENT%IPTS(IPT), COORDS(:,REGION%PARENT%IPTS(IPT)), MASK(IPT)
     END DO
     PRINT *,'   INHERITS',pack(region%parent%ipts(1:region%parent%npts), mask)
     !DO IPT = 1,region%parent%npts
     !   IF ( MASK(IPT) ) PRINT *,COORDS(:,REGION%PARENT%IPTS(IPT))
     !END DO
  END IF
  call append_n( &
       region%ipts, &
       region%npts, &
       pack(region%parent%ipts(1:region%parent%npts), mask), &
       count(mask), &
       unique=.true. )

  deallocate( mask )

end subroutine inherit_points
