$ Declare
   SUBROUTINE EllipticShapeFactorQ (a, c, Q)
!
! According to JCSS Probabilistic model code Part 3: Resistance models
!
! Parameters
   DOUBLE PRECISION, INTENT(IN)    :: a  ! Minimum value of the cyclic component of the total stress intensity
   DOUBLE PRECISION, INTENT(IN)    :: c  ! Maximum value of the cyclic component of the total stress intensity
   DOUBLE PRECISION, INTENT(OUT)   :: Q  ! Crack closure factor
!
! Declarations
   DOUBLE PRECISION                :: Rs
!
! Elliptic shape factor Q
   IF (a/c .LE. 1.D0) THEN
      Q = 1.D0 + 1.464D0*((a/c)**1.65D0)
   ELSE
      STOP 'Error calculating shape factor in Subroutine EllipticShapeFactorQ: a/c > 1'
   END IF
!
! Einde subroutine
   RETURN
   END