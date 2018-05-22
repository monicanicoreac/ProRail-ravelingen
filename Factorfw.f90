$ Declare
   SUBROUTINE Factorfw (a, c, B, W, fw)
!
! Correction term fw in SIF for eliptical flaws
! Fw according to BS 7910:2005 M.3 Flat plates
!
! Parameters
   USE PublicPar
   IMPLICIT NONE
!
   DOUBLE PRECISION, INTENT(IN)    :: a         ! Crack depth
   DOUBLE PRECISION, INTENT(IN)    :: c         ! FT = 2: Crack width
   DOUBLE PRECISION, INTENT(IN)    :: B         ! FT = 2: Plate thickness
   DOUBLE PRECISION, INTENT(IN)    :: W         ! Plate width
   DOUBLE PRECISION, INTENT(OUT)   :: fw        ! Correction term fw
!
! Declarations
   DOUBLE PRECISION                :: lambda
!
! Constants
   DOUBLE PRECISION                :: pi
   pi = 3.14159265358979D0
!
! Fw for copes
   IF (FlawType .EQ. 9) THEN
      fw = DSQRT(1.D0 / DCOS(pi * a / W ))
!
! Fw for edge flaws in plates BS 7910:2005 M.3.5
   ELSE IF (FlawType .EQ. 7) THEN
      fw = 1.D0
!
! Fw for surface flaws in bolts BS 7910:2005 M.6.3.2
   ELSE IF (FlawType .EQ. 6) THEN
      fw = 1.D0
!
! Fw for corner flaws in plates BS 7910:2005 M.3.6
   ELSE IF (FlawType .EQ. 5) THEN
	  lambda = (c / W) * DSQRT(a / B)
	  fw = 1.D0 - 0.2D0 * lambda + 9.4D0 * lambda**2 - 19.4D0 * lambda**3 + 27.1D0 * lambda**4
!
! Fw for semi-circular surface flaw in round bar BS 7910:2005 M.3.X
   ELSE IF (FlawType .EQ. 4) THEN
	  fw = 1.D0
!
! Fw for through-width flaws in plates BS 7910:2005 M.3.3
   ELSE IF (FlawType .EQ. 3) THEN
	  fw = 1.D0
!
! Fw for through-thickness flaws in plates BS 7910:2005 M.3.1
   ELSE IF (FlawType .EQ. 1) THEN
      fw = DSQRT(1.D0 / DCOS(pi * a / W ))
!
! Fw for surface flaws in plates BS 7910:2005 M.3.2.2
   ELSE IF (FlawType .EQ. 2) THEN
      IF (c/W .LE. 0.41D0) THEN
         fw = DSQRT(1.D0 / DCOS(pi * c / W * DSQRT(a / B)))
      ELSE
!         STOP 'Error calculating f_w in Subroutine GeometryFactorY: 2 c/W > 0.4'
         STOP 'Error calculating f_w in Subroutine GeometryFactorY: 2 c/W > 0.41'
      END IF
   ELSE
      STOP 'Error in Subroutine Factorfw: Unknown flawtype'
   END IF
!
! End subroutine
   RETURN
   END