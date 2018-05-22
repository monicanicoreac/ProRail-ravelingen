$ Declare
   SUBROUTINE SIMFactorMb (a, c, B, W, theta, Mb)
!
! Stress intensity magnification factor Mb
! Mb according to BS 7910:2005 M.3 Flat Plates
!
! Parameters
   USE PublicPar
   IMPLICIT NONE
!
   DOUBLE PRECISION, INTENT(IN)    :: a         ! FT = 2: Crack depth
   DOUBLE PRECISION, INTENT(IN)    :: c         ! FT = 2: Crack width
   DOUBLE PRECISION, INTENT(IN)    :: B         ! FT = 2: Plate thickness
   DOUBLE PRECISION, INTENT(IN)    :: W         ! FT = 2: Plate width
   DOUBLE PRECISION, INTENT(IN)    :: theta     ! FT = 2: Angle in radians (0 => c, pi/2 => a)
   DOUBLE PRECISION, INTENT(OUT)   :: Mb        ! Stress intensity magnification factor Mb
!
! Declarations
   DOUBLE PRECISION                :: Mm, Mmmax
   DOUBLE PRECISION                :: as, cs
   DOUBLE PRECISION                :: H
   DOUBLE PRECISION                :: H1
   DOUBLE PRECISION                :: H2
   DOUBLE PRECISION                :: G1
   DOUBLE PRECISION                :: G2
   DOUBLE PRECISION                :: q, g, radius
   DOUBLE PRECISION                :: adivD
!
! Constants
   DOUBLE PRECISION :: pi
   pi = 3.14159265358979D0
!
! M for cope
   IF (FlawType .EQ. 9) THEN
      Mb = 1.D0
!
! Mm for through-thickness deck plate crack
   ELSE IF (FlawType .EQ. 8) THEN
      Mb = 0.375D0
!
! Mb for edge flaws in plates BS 7910:2005 M.3.5
   ELSE IF (FlawType .EQ. 7) THEN
      CALL SIMFactorMm (a, c, B, W, theta, Mm, Mmmax)
      Mb = Mm
!
! Mb for surface flaws in bolts BS 7910:2005 M.6.3.2
   ELSE IF (FlawType .EQ. 6) THEN
      radius = B / 2.D0
	  adivD = a / (2.D0 * radius)
	  IF (adivD .GT. 0.65D0) THEN
		 STOP 'Incorrect input in Subroutine SIMFactorMb: flaw in bolt with a/2r > 0.65'
	  ELSE IF (adivD .LT. 0.004D0) THEN
		 STOP 'Incorrect input in Subroutine SIMFactorMb: flaw in bolt with a/2r < 0.004'
	  ELSE
		 Mb = 0.6301D0 + 0.03488D0 * adivD - 3.3365D0 * (adivD)**2 + 13.406D0 * (adivD)**3 - 6.0021D0 * (adivD)**4
	  END IF
!
! Mb for corner flaws in plates BS 7910:2005 M.3.6
   ELSE IF (FlawType .EQ. 5) THEN
	  IF ((a/c) .LT. 0.2D0) THEN
!         STOP 'Incorrect input in Subroutine SIMFactorMb: corner flaw with a/c < 0.2'
		 as = c * 0.2D0
		 cs = c
	  ELSE IF ((a/c) .GT. 2.D0) THEN
!         STOP 'Incorrect input in Subroutine SIMFactorMb: corner flaw with a/c > 2.0'
		 as = a
		 cs = a * 0.5D0
	  ELSE
	     as = a
		 cs = c
	  END IF
	  IF ((as/cs) .LE. 1.D0) THEN
		 G1 = -1.22D0 - 0.12D0 * (as / cs)
		 G2 = 0.64D0 - 1.05D0 * (as / cs)**0.75D0 + 0.47D0 * (as / cs)**1.5D0
		 H1 = 1.D0 - 0.34D0 * (as / B) - 0.11D0 * (as / cs) * (as / B)
		 q = 0.2D0 + (as / cs) + 0.6D0 * (as / B)
	  ELSE
	     G1 = -2.11D0 + 0.77D0 * (cs / as)
		 G2 = 0.64D0 - 0.72D0 * (cs / as)**0.75D0 + 0.14D0 * (cs / as)**1.5D0
		 H1 = 1.D0 - (0.04D0 + 0.41D0*(cs/as))*(as/B) + (0.55D0 - 1.93D0*(cs/as)**0.75D0 + 1.38D0*(cs/as)**1.5D0)*(as/B)**2.D0
		 q = 0.2D0 + (cs / as) + 0.6D0 * (as / B)
	  END IF
	  H2 = 1.D0 + G1 * (as / B) + G2 * (as / B)**2
	  H = H1 + (H2 - H1) * (DSIN(theta))**q
      CALL SIMFactorMm (a, c, B, W, theta, Mm, Mmmax)
	  Mb = H * Mm
!
! Mb for semi-circular surface flaw in round bar BS 7910:2005 M.3.x
   ELSE IF (FlawType .EQ. 4) THEN
	  radius = B / 2.D0
      g = 1.84D0 / pi * DSQRT (DTAN (pi * a / 4.D0 / radius) / (pi * a / 4.D0 / radius)) / DCOS (pi * a / 4.D0 / radius)
      IF (a/radius .LE. 1.2D0) THEN
		 Mb = g * (0.923D0 + 0.199D0 * (1.D0 - DSIN(pi * a / 4.D0 / radius))**4.D0)
	  ELSE
         STOP 'Incorrect input in Subroutine SIMFactorMb: round bar with a/2r > 0.6'
	  END IF
!
! Mb for through-width flaws in plates BS 7910:2005 M.3.3
   ELSE IF (FlawType .EQ. 3) THEN
      IF (a/B .LE. 0.6D0) THEN
		 Mb = 1.12D0 - 1.39D0 * (a/B) + 7.32D0 * (a/B)**2.D0 - 13.1D0 * (a/B)**3.D0 + 14.D0 * (a/B)**4.D0
	  ELSE
         STOP 'Incorrect input in Subroutine SIMFactorMb: long surface flaw with a/B > 0.6'
	  END IF
!
! Mb for through-thickness flaws in plates BS 7910:2005 M.3.1
   ELSE IF (FlawType .EQ. 1) THEN
      Mb = 1
!
! Mb for surface flaws in plates BS 7910:2005 M.3.2.3 
   ELSE IF (FlawType .EQ. 2) THEN
!
! Check input
      IF (a .LT. 0.0D0) STOP 'Incorrect input in Subroutine SIMFactorMb: a < 0'
      IF (c .LT. 0.0D0) STOP 'Incorrect input in Subroutine SIMFactorMb: c < 0'
      IF (a/c .LE. 0.2D0) THEN
!         IF (a/B .GE. 1.25D0 * (a/c + 0.6D0)) STOP 'a/B Outside validity conditions in Subroutine SIMFactorMb'
      ELSE IF (a/c .LE. 2.0D0) THEN
         IF (a/B .GE. 1.D0) STOP 'a/B Outside validity conditions in Subroutine SIMFactorMb'
      END IF
	  IF ((a/c) .GT. 2.D0) THEN
!         STOP 'Incorrect input in Subroutine SIMFactorMm: corner flaw with a/c > 2.0'
		 as = a
		 cs = a * 0.5D0
	  ELSE
	     as = a
		 cs = c
	  END IF
!
! Factors q, H, H1, H2, G1 and G2
      IF (as/cs .LE. 1.D0) THEN
	     q  = 0.2D0 + as / cs + 0.6D0 * as / B
         H1 = 1.D0 - 0.34D0 * as / B - 0.11D0 * as * as / cs / B
		 G1 = -1.22D0 - 0.12D0 * as / cs
		 G2 = 0.55D0 - 1.05D0 * ((as / cs)**0.75D0) + 0.47D0 * ((as / cs)**1.5D0)
	  ELSE IF (as/cs .LE. 2.D0) THEN
	     q  = 0.2D0 + cs / as + 0.6D0 * as / B
         H1 = 1.D0 - (0.04D0 + 0.41D0 * cs / as) * as / B + (0.55D0 - 1.93D0 * ((cs / as)**0.75D0) + 1.38D0 * ((cs / as)**1.5D0)) * as * as / B / B
	     G1 = -2.11D0 + 0.77D0 * cs / as
		 G2 = 0.55D0 - 0.72D0 * ((cs / as)**0.75D0) + 0.14D0 * ((cs / as)**1.5D0)
	  END IF
	  H2 = 1.D0 + G1 * as / B + G2 * as * as / B / B
      H  = H1 + (H2 - H1) * (DSIN(theta)**q)
      CALL SIMFactorMm (as, cs, B, W, theta, Mm, Mmmax)
      Mb = H * Mm
   ELSE
      STOP 'Incorrect input in Subroutine SIMFactorMb: Unknown flawtype'
   END IF
   IF (Mb .LT. 0.D0) Mb = 0.D0

!
! End subroutine
   RETURN
   END