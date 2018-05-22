$ Declare
   SUBROUTINE SIGrefC ( a, c, B, W, Spm, Spb, Sref)
!
! Parameters
   USE PublicPar
   IMPLICIT NONE
!
   DOUBLE PRECISION, INTENT(IN)    :: a     ! Crack depth
   DOUBLE PRECISION, INTENT(IN)    :: c     ! Crack width
   DOUBLE PRECISION, INTENT(IN)    :: B     ! Plate thickness
   DOUBLE PRECISION, INTENT(IN)    :: W     ! Plate width
   DOUBLE PRECISION, INTENT(IN)    :: Spm   ! primary membrane stress
   DOUBLE PRECISION, INTENT(IN)    :: Spb   ! primary bending stress
   DOUBLE PRECISION, INTENT(OUT)   :: Sref  ! reference stress
!
! Declarations
   DOUBLE PRECISION                :: alphadacc
!
! Calculate reference stress
   IF (FlawType .EQ. 2) THEN
      IF (W .GE. 2*(c+B)) THEN
		 alphadacc = (a/B)/(1+(B/c))
	  ELSE
		 alphadacc = (2*a/B)*(c/W)
	  END IF 
	  Sref = (Spb + DSQRT(Spb*Spb+9.D0*Spm*Spm*(1.D0-alphadacc)*(1.D0-alphadacc))) / (3.D0*(1.D0-alphadacc)*(1.D0-alphadacc))
   ELSE
	  Sref = (Spb + DSQRT(Spb*Spb+9.D0*Spm*Spm)) / (3.D0 * (1.D0 - (2.D0*a/W)))
   END IF
!
! End subroutine
   RETURN
   END