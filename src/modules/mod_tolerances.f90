module mod_tolerances

  use mod_math

  implicit none

  real(kind=fp), parameter :: EPSuv = real( 1.e-12, kind=fp )
  real(kind=fp), parameter :: EPSuvsqr = EPSuv**2

  real(kind=fp), parameter :: EPSxyz = real( 1.e-9, kind=fp )
  real(kind=fp), parameter :: EPSxyzsqr = EPSxyz**2

  real(kind=fp), parameter :: EPScollineal = real( 1.e-9, kind=fp )
  real(kind=fp), parameter :: EPScollinealsqr = EPScollineal**2

end module mod_tolerances
