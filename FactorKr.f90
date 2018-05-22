$ Declare
   SUBROUTINE FactorKr ( a, c, B, W, theta, Ym, Yb, Spm, Spb, Ssm, Ssb, K1C, Lr, Kr )
!
! Parameters
   IMPLICIT NONE
   DOUBLE PRECISION, INTENT(IN)    :: a            ! Crack depth
   DOUBLE PRECISION, INTENT(IN)    :: c            ! Crack width
   DOUBLE PRECISION, INTENT(IN)    :: B            ! Plate thickness
   DOUBLE PRECISION, INTENT(IN)    :: W            ! Plate width
   DOUBLE PRECISION, INTENT(IN)    :: theta        ! FT = 2: Angle in radians (0 => c, pi/2 => a)
   DOUBLE PRECISION, INTENT(IN)    :: Ym           ! Geometric Correction Factor for membrane stress
   DOUBLE PRECISION, INTENT(IN)    :: Yb           ! Geometric Correction Factor for bending stress
   DOUBLE PRECISION, INTENT(IN)    :: Spm          ! primary membrane stress
   DOUBLE PRECISION, INTENT(IN)    :: Spb          ! primary bending stress
   DOUBLE PRECISION, INTENT(IN)	   :: Ssm          ! secondary membrane stress
   DOUBLE PRECISION, INTENT(IN)	   :: Ssb          ! secondary bending stress
   DOUBLE PRECISION, INTENT(IN)	   :: K1C          ! fracture toughness
   DOUBLE PRECISION, INTENT(IN)    :: Lr		   ! plasticity factor
   DOUBLE PRECISION, INTENT(OUT)   :: Kr		   ! fracture factor
!
! Declarations
   DOUBLE PRECISION                :: KIp, KIs, KI, chi, rho1, rho
   DOUBLE PRECISION                :: Mm, Mb, Mmmax
!
! Constants
   DOUBLE PRECISION :: pi
   pi = 3.14159265358979D0
!
! Calculate Kr
	  CALL SIMFactorMm ( a, c, B, W, theta, Mm, Mmmax)
	  CALL SIMFactorMb ( a, c, B, W, theta, Mb)
   	  KIp = (Ym * Spm + Yb * Spb) * DSQRT(pi * a)
      KIs = (Mmmax * Ssm + Mb * Ssb) * DSQRT(pi * a)
	  KI = KIp + KIs
	  chi = (KIs * Lr / KIp)
	  rho1 = 0.1D0 * (chi)**0.714D0 - 0.007D0 * chi * chi + 3.D-5 * (chi)**5
	  IF (chi .GT. 4.D0) THEN
		 STOP 'Error in Subroutine FactorKr: KIs*Lr/KIp >4'
	  ELSE IF (Lr .LE. 0.8D0) THEN
		 rho = rho1
	  ELSE IF (Lr .LT. 1.05D0) THEN
		 rho = 4.D0 * rho1 * (1.05D0 - Lr)
	  ELSE
		 rho = 0.D0
	  END IF
	  Kr = KI/K1C + rho
!
! End subroutine
 RETURN
 END