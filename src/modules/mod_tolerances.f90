module mod_tolerances

  use mod_math

  implicit none

  real(kind=fp), parameter :: EPSuv = real( 1.d-12, kind=fp )
  real(kind=fp), parameter :: EPSuvsqr = EPSuv**2

  real(kind=fp), parameter :: EPSxyz = real( 1.d-9, kind=fp )
  real(kind=fp), parameter :: EPSxyzsqr = EPSxyz**2

  real(kind=fp), parameter :: EPScollineal = real( 4.d-9, kind=fp )
  real(kind=fp), parameter :: EPScollinealsqr = EPScollineal**2


  ! Intersection polyline tracing
  real(kind=fp), parameter :: TOLchord = real(1.d-3, kind=fp)
  real(kind=fp), parameter :: FRACcurvature_radius = 2._fp*sqrt(TOLchord*(2._fp - TOLchord))
  real(kind=fp), parameter :: TOLh = real(1.d-2, kind=fp)
  real(kind=fp), parameter :: TOLhsqr = tolh**2

end module mod_tolerances
