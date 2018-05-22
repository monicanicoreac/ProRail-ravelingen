$ Declare
   SUBROUTINE CrackClosureFactorU (Kmincyc, Kmaxcyc, Kres, U)
!
! Crack closure factor U
! (Shetty and Baker, 1990)
!
! Parameters
   DOUBLE PRECISION, INTENT(IN)    :: Kmincyc  ! Minimum value of the cyclic component of the total stress intensity
   DOUBLE PRECISION, INTENT(IN)    :: Kmaxcyc  ! Maximum value of the cyclic component of the total stress intensity
   DOUBLE PRECISION, INTENT(IN)    :: Kres     ! Mean level of the total stress intensity, due to welding residual stresses
   DOUBLE PRECISION, INTENT(OUT)   :: U        ! Crack closure factor
!
! Declarations
   DOUBLE PRECISION                :: Rs
!
! Check input
   IF (Kmincyc .GT. 0.0D0) STOP 'Incorrect input in Subroutine CrackClosureFactorU: Kmincyc > 0'
   IF (Kmaxcyc .LT. 0.0D0) STOP 'Incorrect input in Subroutine CrackClosureFactorU: Kmaxcyc < 0'
!
! Stress ratio Rs 
   Rs = (Kmincyc + Kres) / (Kmaxcyc + Kres)
!
! Crack closure factor U
   IF (Rs .GE. 0.D0) THEN
      U = 1.D0
   ELSE
      U = 1.D0 / (1.D0 - Rs)
   END IF
!
! Einde subroutine
   RETURN
   END