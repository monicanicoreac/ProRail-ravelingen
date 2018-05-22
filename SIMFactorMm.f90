$ Declare
   SUBROUTINE SIMFactorMm (a, c, B, W, theta, Mm, Mmmax)
!
! Stress intensity magnification factor Mm
! Mm according to BS 7910:2005 M.3 Flat Plates   
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
   DOUBLE PRECISION, INTENT(OUT)   :: Mm		 ! Stress intensity magnification factor Mm
   DOUBLE PRECISION, INTENT(OUT)   :: Mmmax		 ! Stress intensity magnification factor Mm for maximum stress
!
! Declarations
   DOUBLE PRECISION                :: M1, M2, M3
   DOUBLE PRECISION                :: as, cs
   DOUBLE PRECISION                :: g, radius, adivD, g1, g2
   DOUBLE PRECISION                :: ftheta
   DOUBLE PRECISION                :: fi
!
! Constants
   DOUBLE PRECISION :: pi
   pi = 3.14159265358979D0
!
! Mm for cope
   IF (FlawType .EQ. 9) THEN
      !M1 = 0.076071 * DEXP (-0.16676*a) + 0.088289 * DEXP (-0.010558*a)  ! Positive stress
	  !M2 = -0.31589 * DEXP (-0.14504*a) + -0.24447 * DEXP (0.00052707*a) ! Negative stress
      Mmmax = 1.12D0*(0.5372D0 * DEXP (-0.01056D0*a) + 0.4628D0 * DEXP (-0.1668D0*a)) ! Range
	  Mm = 1.12D0*(0.4344D0 * DEXP (-0.0006855D0*a) + 0.5608D0 * DEXP (-0.1384D0*a))
!
! Mm for through-thickness deck plate crack
   ELSE IF (FlawType .EQ. 8) THEN
      Mm = 1.D0
	  Mmmax = 1.D0
!
! Mm for edge flaws in plates BS 7910:2005 M.3.5
   ELSE IF (FlawType .EQ. 7) THEN
      Mm = 1.12D0 - 0.23D0 * (a/W) + 10.6D0 * (a/W)**2 - 21.7D0 * (a/W)**3 + 30.4D0 * (a/W)**4
	  Mmmax = Mm
!
! Mm for surface flaws in bolts BS 7910:2005 M.6.3.2
   ELSE IF (FlawType .EQ. 6) THEN
      radius = B / 2.D0
	  adivD = a / (2.D0 * radius)
	  IF (adivD .GT. 0.65D0) THEN
		 STOP 'Incorrect input in Subroutine SIMFactorMm: flaw in bolt with a/2r > 0.65'
	  ELSE IF (adivD .LT. 0.004D0) THEN
		 STOP 'Incorrect input in Subroutine SIMFactorMm: flaw in bolt with a/2r < 0.004'
	  ELSE
		Mm = 2.043D0 * DEXP (-31.332D0 * adivD) + 0.6507D0 + 0.5367D0 * adivD + 3.0469D0 * (adivD)**2 - 19.504D0 * (adivD)**3 + 45.647 * (adivD)**4
		Mmmax = Mm
	  END IF
!
! Mm for corner flaws in plates BS 7910:2005 M.3.6
   ELSE IF (FlawType .EQ. 5) THEN
	  IF ((a/c) .LT. 0.2D0) THEN
!         STOP 'Incorrect input in Subroutine SIMFactorMm: corner flaw with a/c < 0.2'
		 as = c * 0.2D0
		 cs = c
	  ELSE IF ((a/c) .GT. 2.D0) THEN
!         STOP 'Incorrect input in Subroutine SIMFactorMm: corner flaw with a/c > 2.0'
		 as = a
		 cs = a * 0.5D0
	  ELSE
	     as = a
		 cs = c
	  END IF
	  IF ((a/B) .GT. 1.D0) THEN
         STOP 'Incorrect input in Subroutine SIMFactorMm: corner flaw with a/B > 1.0'
	  ELSE IF ((as/cs) .LE. 1.D0) THEN
		 M1 = 1.08D0 - 0.03D0 * (as / cs)
		 M2 = (1.06D0 / (0.3D0 + (as / cs))) - 0.44D0
		 M3 = -0.5D0 + 0.25D0 * (as / cs) + 14.8D0 * (1.D0 - (as / cs))**15.D0
		 g1 = 1.D0 + (0.08D0 + 0.4D0 * (as/B)*(as/B)) * (1.D0 - DSIN(theta))**3.D0
		 g2 = 1.D0 + (0.08D0 + 0.15D0 * (as/B)*(as/B)) * (1.D0 - DCOS(theta))**3.D0
	     ftheta = ((as / cs * DCOS(theta))**2 + DSIN(theta)**2)**0.25D0
	     fi = DSQRT(1 + 1.464D0 * ((as / cs)**1.65D0))
	  ELSE
	     M1 = (1.08D0 - 0.03D0 * (cs / as)) * (cs / as)**0.5D0
		 M2 = 0.375 * (cs / as) * (cs / as)
		 M3 = -0.25D0 * (cs / as) * (cs / as)
		 g1 = 1.D0 + (0.08D0 + 0.4D0 * (cs/B)*(cs/B)) * (1.D0 - DSIN(theta))**3.D0
		 g2 = 1.D0 + (0.08D0 + 0.15D0 * (cs/B)*(cs/B)) * (1.D0 - DCOS(theta))**3.D0
	     ftheta = ((cs / as * DSIN(theta))**2 + DCOS(theta)**2)**0.25D0
	     fi = DSQRT(1 + 1.464D0 * ((cs / as)**1.65D0))
	  END IF
      Mm = (M1 + M2 * (as/B)**2 + M3 * (as/B)**4) * g1 * g2 * ftheta / fi
	  Mmmax = Mm
!
! Mm for semi-circular surface flaw in round bar BS 7910:2005 M.3.x
   ELSE IF (FlawType .EQ. 4) THEN
	  radius = B / 2.D0
      g = 1.84D0 / pi * DSQRT (DTAN (pi * a / 4.D0 / radius) / (pi * a / 4.D0 / radius)) / DCOS (pi * a / 4.D0 / radius)
      IF (a/radius .LE. 1.2D0) THEN
		 Mm = g * (0.752D0 + 2.02D0 * (a/2.D0/radius) + 0.37D0 * (1.D0 - DSIN(pi * a / 4.D0 / radius))**3.D0)
		 Mmmax = Mm
	  ELSE
         STOP 'Incorrect input in Subroutine SIMFactorMm: round bar with a/2r > 0.6'
	  END IF
	  !
! Mm for through-width flaws in plates BS 7910:2005 M.3.3
   ELSE IF (FlawType .EQ. 3) THEN
      IF (a/B .LE. 0.6D0) THEN
		 Mm = 1.12D0 - 0.23D0 * (a/B) + 10.6D0 * (a/B)**2.D0 - 21.7D0 * (a/B)**3.D0 + 30.4D0 * (a/B)**4.D0
		 Mmmax = Mm
	  ELSE
         STOP 'Incorrect input in Subroutine SIMFactorMm: long surface flaw with a/B > 0.6'
	  END IF
	  
!
! Mm for through-thickness flaws in plates BS 7910:2005 M.3.1
   ELSE IF (FlawType .EQ. 1) THEN
      Mm = 1.D0
	  Mmmax = Mm
!
! Mm for surface flaws in plates BS 7910:2005 M.3.2.2
   ELSE IF (FlawType .EQ. 2) THEN
!
! Check input
!      write (*,*) a, c, B
      IF (a .LT. 0.0D0) STOP 'Incorrect input in Subroutine SIMFactorMm: a < 0'
      IF (c .LT. 0.0D0) STOP 'Incorrect input in Subroutine SIMFactorMm: c < 0'
!      IF (a/c .LE. 0.2D0) THEN
!         IF (a/B .GE. 1.25D0 * (a/c + 0.6D0)) STOP 'a/B Outside validity conditions in Subroutine SIMFactorMm'
!      END IF
      IF (a/B .GE. 1.D0) STOP 'a/B Outside validity conditions in Subroutine SIMFactorMm'
	  IF ((a/c) .GT. 2.D0) THEN
!         STOP 'Incorrect input in Subroutine SIMFactorMm: corner flaw with a/c > 2.0'
		 as = a
		 cs = a * 0.5D0
	  ELSE
	     as = a
		 cs = c
	  END IF
!
! Factors M1, M2, M3, g, ftheta and fi
! Remark: The functions of M2 and M3 are discontinuing at a/c = 1
      IF (as/cs .LE. 1.D0) THEN
         M1     = 1.13D0 - 0.09D0 * as / cs
	     M2     = 0.89D0 / (0.2D0 + as / cs) - 0.54D0
	     M3     = 0.5D0 - 1.D0 / (0.65D0 + as / cs) + 14.D0 * ((1.D0 - (as / cs))**24.D0)
         g      = 1.D0 + (0.1D0 + 0.35D0 * as * as / B / B) * ((1.D0 - DSIN(theta))**2)
	     ftheta = ((as / cs * DCOS(theta))**2 + DSIN(theta)**2)**0.25D0
	     fi     = DSQRT(1 + 1.464D0 * ((as / cs)**1.65D0))
      ELSE IF (as / cs .LE. 2.D0) THEN
         M1     = DSQRT(cs / as) * (1.D0 + 0.04D0 * cs / as)
	     M2     = 0.2D0 * ((cs / as)**4)
	     M3     = -0.11D0 * ((cs / as)**4)
         g      = 1.D0 + (0.1D0 + 0.35D0 * cs * as / B / B) * ((1.D0 - DSIN(theta))**2)
   	     ftheta = ((cs / as * DSIN(theta))**2 + DCOS(theta)**2)**0.25D0
	     fi     = DSQRT(1 + 1.464D0 * ((cs / as)**1.65D0))
      END IF
!
! Mm
      Mm = (M1 + (M2 + M3 * as * as / B / B) * as * as / B / B) * g * ftheta / fi
	  Mmmax = Mm
	  !Mmmax = 1.D0
   ELSE
      STOP 'Incorrect input in Subroutine SIMFactorMm: Unknown flawtype'
   END IF
   IF (Mm.LT.0.D0) Mm = 0.D0
!
! End subroutine
   RETURN
   END