$ Declare
   SUBROUTINE GeometryFactorY (StressType, a, c, b, W, theta, Y)
!
! According to BS 7910:2005 M.3 Flat plates
! but: g is different en the complete elliptic integral of the second kind is missing.
!
! Parameters
   INTEGER,          INTENT(IN)    :: StressType  ! 1 = Surface flaws membrane loading
                                                  ! 2 = Surface flaws bending loading
   DOUBLE PRECISION, INTENT(IN)    :: a           ! Crack depth
   DOUBLE PRECISION, INTENT(IN)    :: c           ! Crack width
   DOUBLE PRECISION, INTENT(IN)    :: B           ! Plate thickness
   DOUBLE PRECISION, INTENT(IN)    :: W           ! Plate width
   DOUBLE PRECISION, INTENT(IN)    :: theta       ! Angle in radians (0 => c, pi/2 => a)
   DOUBLE PRECISION, INTENT(OUT)   :: Y           ! Geometry factor
!
! Declarations
   DOUBLE PRECISION                :: fw
   DOUBLE PRECISION                :: M
!
! Constants
   DOUBLE PRECISION                :: pi
   pi = 3.14159265358979D0
!
! Calculate stress intensity magnification factor M (is Mm or Mb)
   CALL SIMFactorM (StressType, a, c, B, theta, M)
!
! Factor Fw
   IF (StressType .EQ. 1 .OR. StressType .EQ. 2) THEN
      IF (c/W .LE. 0.4D0) THEN
         fw = DSQRT(1.D0 / DCOS(pi * c / W * DSQRT(a / B)))
      ELSE
         STOP 'Error calculating f_w in Subroutine GeometryFactorY: c/b > 0.4'
      END IF
   ELSE
      STOP 'Error in Subroutine GeometryFactorY: Unknown flaw and stress type'
   END IF
!
! Geometry factor Y 
   Y = M * fw
!
! End subroutine
   RETURN
   END