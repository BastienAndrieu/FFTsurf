module mod_options

  use mod_math
  use mod_tolerances
  
  implicit none
  
  integer, parameter :: PARAM_stringlen = 99

  type type_options
     integer                    :: mode = 0
     real                       :: timespan = 0.0
     real                       :: timestep = 1.0
     character(PARAM_stringlen) :: directory = 'defaultcase/'
     integer                    :: propagation_law = 1
     real(kind=fp)              :: chord_err = real(5.d-3, kind=fp)!TOLchord
     real(kind=fp)              :: hmin = real(1.d-3, kind=fp)!PARAM_hmin
     real(kind=fp)              :: hmax = real(1.d-2, kind=fp)!PARAM_hmax
     logical                    :: reprise = .false.
     logical                    :: from_msh = .false.
  end type type_options

contains

  subroutine read_options( &
       fileoptions, &
       options )
    use mod_util
    implicit none
    character(*),       intent(in)    :: fileoptions
    type(type_options), intent(inout) :: options
    integer                           :: fid, stat
    character(PARAM_stringlen)        :: line, trimlin

    call get_free_unit(fid)

    open( &
         unit = fid, &
         file = fileoptions, &
         action = 'read' )

    do
       read (fid, "(a)", iostat=stat) line

       if (stat < 0) exit ! fin de fichier

       ! sauter commentaires et lignes vides
       if ((line(1:1) == "#") .or. (len_trim(line) == 0)) cycle

       trimlin = trim(line)

       if (trimlin(1:8) == "timestep") then
          ! timestep
          read (fid,*) options%timestep

       elseif (trimlin(1:8) == "timespan") then
          ! timespan
          read (fid,*) options%timespan

       elseif (trimlin(1:9) == "directory") then
          ! directory
          read (fid,"(a)") options%directory

       elseif (trimlin(1:11) == "propagation") then
          ! propagation law
          read (fid,*) options%propagation_law

       elseif (trimlin(1:4) == "mode") then
          ! mode
          read (fid,*) options%mode
          
       elseif (trimlin(1:5) == "chord") then
          ! chordal error
          read (fid,*) options%chord_err
          
       elseif (trimlin(1:4) == "hmin") then
          ! hmin
          read (fid,*) options%hmin
          
       elseif (trimlin(1:4) == "hmax") then
          ! hmax
          read (fid,*) options%hmax

       end if
    end do

    close(fid)

  end subroutine read_options






  subroutine print_options( &
       self)
    implicit none
    type(type_options), intent(in) :: self

    print *,''
    print *,'========= OPTIONS FFTsurf ========='
    print *,'timespan  = ', self%timespan
    print *,'timestep  = ', self%timestep
    print *,'directory = ', trim(self%directory)
    print *,'prop. law = ', self%propagation_law
    print *,'mode      = ', self%mode
    print *,'chord err = ', self%chord_err
    print *,'hmin      = ', self%hmin
    print *,'hmax      = ', self%hmax
    print *,'reprise?    ', self%reprise
    print *,'from msh?   ', self%from_msh
    print *,'==================================='
    print *,''

  end subroutine print_options




end module mod_options
