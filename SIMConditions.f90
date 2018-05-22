$ Declare
   SUBROUTINE SIMConditions (a, c, B, W, err)
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
   INTEGER, INTENT(OUT)            :: err       ! Error message, 1 = error
!
! Declarations
   DOUBLE PRECISION                :: M1, M2, M3
   DOUBLE PRECISION                :: as, cs
   DOUBLE PRECISION                :: g, radius, adivD, g1, g2
   DOUBLE PRECISION                :: ftheta
   DOUBLE PRECISION                :: fi
!
! Constants
   err = 0
!
!  cope
   IF (FlawType .EQ. 9) THEN
!
!  through-thickness deck plate crack
   ELSE IF (FlawType .EQ. 8) THEN
!
!  edge flaws in plates BS 7910:2005 M.3.5
   ELSE IF (FlawType .EQ. 7) THEN
!
!  surface flaws in bolts BS 7910:2005 M.6.3.2
   ELSE IF (FlawType .EQ. 6) THEN
      radius = B / 2.D0
	  adivD = a / (2.D0 * radius)
	  IF (adivD .GT. 0.65D0) THEN
		 err = 1
	  ELSE IF (adivD .LT. 0.004D0) THEN
		 err = 1
	  END IF
!
!  corner flaws in plates BS 7910:2005 M.3.6
   ELSE IF (FlawType .EQ. 5) THEN
	  IF ((a/c) .LT. 0.2D0) THEN
		 err = 1
	  ELSE IF ((a/c) .GT. 2.D0) THEN
		 err = 1
	  ELSE IF ((a/B) .GT. 1.D0) THEN
         err = 1
	  END IF
!
!  semi-circular surface flaw in round bar BS 7910:2005 M.3.x
   ELSE IF (FlawType .EQ. 4) THEN
      IF (a/radius .GT. 1.2D0) THEN
		 err = 1
	  END IF
!
!  through-width flaws in plates BS 7910:2005 M.3.3
   ELSE IF (FlawType .EQ. 3) THEN
      IF (a/B .GT. 0.6D0) THEN
	     err = 1
	  END IF
!
!  through-thickness flaws in plates BS 7910:2005 M.3.1
   ELSE IF (FlawType .EQ. 1) THEN
!
!  surface flaws in plates BS 7910:2005 M.3.2.2
   ELSE IF (FlawType .EQ. 2) THEN
      IF (a .LT. 0.0D0) THEN
	     err = 1
      ELSEIF (c .LT. 0.0D0) THEN
	     err = 1
      ELSEIF (a/B .GE. 1.D0) THEN
         err = 1
	  ELSEIF ((a/c) .GT. 2.D0) THEN
         err = 1
	  ELSEIF (c/W .GT. 0.4D0) THEN
	     err = 1
	  END IF
!
! Unknown flaw type
   ELSE
      STOP 'Incorrect input in Subroutine SIMConditions: Unknown flawtype'
   END IF
!
! End subroutine
   RETURN
   END