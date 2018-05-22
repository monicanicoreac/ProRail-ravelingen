$ Declare
   SUBROUTINE CrackGrowthThreshold (DK0, a, f, A0, dKth)
!
! parameters
   USE PublicPar
   IMPLICIT NONE
!
   DOUBLE PRECISION, INTENT(IN)    :: DK0
   DOUBLE PRECISION, INTENT(IN)    :: a
   DOUBLE PRECISION, INTENT(IN)    :: f
   DOUBLE PRECISION, INTENT(IN)    :: A0
   DOUBLE PRECISION, INTENT(OUT)   :: dKth
!
! declarations
   DOUBLE PRECISION  :: Cth, Cthm, Cthp, Rcl, asc0
   DOUBLE PRECISION  :: sc, Rth
!
! material: steel in air (NasGro equations, parameters
! acc. to Fitnet and German dissertation)
!
! Constants and parameters
   Cthm  = 0.2D0
   Cthp  = 2.D0
   Rcl   = 0.75D0
   asc0  = 0.0381D0           ! Correction factor for short cracks
   sc = DSQRT(a/(a+asc0))     ! cracks
!
! Constant value for dKth in case of large R-values
   IF (R.GE.Rcl) THEN
	  Rth = Rcl
	  Cth = Cthp
   ELSE IF (R.LT.0.D0) THEN
      Rth = R
	  Cth = Cthm
   ELSE
      Rth = R
	  Cth = Cthp
   END IF
!
! Determine dKth
   dKth = DK0*sc*(((1-A0)*(1-Rth)/(1-f))**(1+Cth*Rth))
!
! End subroutine
   RETURN
   END

