module mod_tolerances

  use mod_math

  implicit none

  real(kind=fp), parameter :: EPSfp = epsilon( 1._fp )
  real(kind=fp), parameter :: EPSfpsqr = EPSfp**2

  real(kind=fp), parameter :: EPSuv = real( 1.d-12, kind=fp )
  real(kind=fp), parameter :: EPSuvsqr = EPSuv**2

  real(kind=fp), parameter :: EPSxyz = real( 1.d-9, kind=fp )
  real(kind=fp), parameter :: EPSxyzsqr = EPSxyz**2

  real(kind=fp), parameter :: EPScollineal = real( 4.d-9, kind=fp )
  real(kind=fp), parameter :: EPScollinealsqr = EPScollineal**2

end module mod_tolerances
