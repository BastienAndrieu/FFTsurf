module mod_tolerances

  use mod_math

  implicit none

  real(kind=MATHpr), parameter :: EPSuv = real( 1.e-12, kind=MATHpr )
  real(kind=MATHpr), parameter :: EPSuvsqr = EPSuv**2

  real(kind=MATHpr), parameter :: EPSxyz = real( 1.e-9, kind=MATHpr )
  real(kind=MATHpr), parameter :: EPSxyzsqr = EPSxyz**2

  real(kind=MATHpr), parameter :: EPScollineal = real( 1.e-9, kind=MATHpr )
  real(kind=MATHpr), parameter :: EPScollinealsqr = EPScollineal**2

end module mod_tolerances
