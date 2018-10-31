module mod_import

  use mod_math

contains

  subroutine read_msh( &
       filename, &
       xyz, &
       nv, &
       tri, &
       ref, &
       nt )
    use mod_util
    implicit none
    character(*),               intent(in)    :: filename
    real(kind=fp), allocatable, intent(inout) :: xyz(:,:)
    integer,                    intent(out)   :: nv
    integer, allocatable,       intent(inout) :: tri(:,:)
    integer, allocatable,       intent(inout) :: ref(:)
    integer,                    intent(out)   :: nt
    integer                                   :: fid, err
    character(100)                            :: line
    real(kind=fp)                             :: vert(4)
    integer                                   :: ne, ntags, elem(10)
    logical, allocatable                      :: keepv(:)
    integer, allocatable                      :: old2new(:)
    integer                                   :: ivert, ielem

    call get_free_unit(fid)
    open( &
         unit = fid, &
         file = filename, &
         action = 'read' )
    do ! <-------------------------------------------------------+
       read (fid, *, iostat=err) line                            !
       if ( err /= 0 ) exit                                      !
       line = trim(line)                                         !
       !                                                         !
       if ( line(1:6) == '$Nodes' ) then ! <-----------------+   !
          read (fid,*) nv                                    !   !
          if ( allocated(xyz) ) then ! <---------------+     !   !
             if ( size(xyz,2) < nv ) deallocate(xyz)   !     !   !
          end if ! <-----------------------------------+     !   !
          if ( .not.allocated(xyz) ) allocate(xyz(3,nv))     !   !
          !                                                  !   !
          do ivert = 1,nv ! <-------------+                  !   !
             read (fid,*) vert            !                  !   !
             xyz(1:3,ivert) = vert(2:4)   !                  !   !
          end do ! <----------------------+                  !   !
          !                                                  !   !
       elseif ( line(1:9) == '$Elements' ) then ! -----------+   !
          read (fid,*) ne                                    !   !
          if ( allocated(tri) ) then ! <-----------------+   !   !
             if ( size(tri,2) < ne ) deallocate(tri)     !   !   !
          end if ! <-------------------------------------+   !   !
          if ( .not.allocated(tri) ) allocate(tri(3,ne))     !   !
           if ( allocated(ref) ) then ! <----------------+   !   !
             if ( size(ref) < ne ) deallocate(ref)       !   !   !
          end if ! <-------------------------------------+   !   !
          if ( .not.allocated(ref) ) allocate(ref(ne))       !   !
          !                                                  !   !
          nt = 0                                             !   !
          do ielem = 1,ne ! <----------------------------+   !   !
             read (fid,*) elem(1:3)                      !   !   !
             if ( elem(2) == 2 ) then ! <------------+   !   !   !
                nt = nt + 1                          !   !   !   !
                ntags = elem(3)                      !   !   !   !
                backspace(fid)                       !   !   !   !
                read (fid,*) elem(1:ntags+6)         !   !   !   !
                tri(1:3,nt) = elem(ntags+4:ntags+6)  !   !   !   !
                ref(nt) = elem(ntags+3)              !   !   !   !
             end if ! <------------------------------+   !   !   !
          end do ! <-------------------------------------+   !   !
          !                                                  !   !
       end if ! <--------------------------------------------+   !
    end do ! <---------------------------------------------------+
    close(fid)


    ! fix
    allocate(keepv(nv))
    keepv(1:nv) = .false.
    do ielem = 1,nt
       keepv(tri(1:3,ielem)) = .true.
    end do

    allocate(old2new(nv))
    old2new(1:nv) = 0
    nv = 0
    do ivert = 1,size(old2new)
       if ( keepv(ivert) ) then
          nv = nv + 1
          old2new(ivert) = nv
       end if
    end do

    do ielem = 1,nt
       tri(1:3,ielem) = old2new(tri(1:3,ielem))
    end do

    xyz(1:3,1:nv) = reshape( &
         pack(xyz(1:3,1:size(keepv)), spread(keepv, dim=1, ncopies=3)), &
         [3,nv] )
    
    deallocate(keepv, old2new)

  end subroutine read_msh
  
end module mod_import
