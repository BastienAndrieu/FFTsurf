module mod_constants
  
  implicit none
  
  integer,       parameter :: fp = SELECTED_REAL_KIND(15,307)!kind(1.d0)
  
  real(kind=fp), parameter :: CSTpi = real( 3.14159265358979323846, kind=fp )

end module mod_constants
