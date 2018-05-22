$ Declare
   SUBROUTINE CrackGrowthLaw (A1, m1, A2, m2, DK0, K1C, fy, fu, Smax, dK, a, dadn)
!
! parameters
   USE PublicPar
   IMPLICIT NONE
!
   DOUBLE PRECISION, INTENT(IN)    :: A1
   DOUBLE PRECISION, INTENT(IN)    :: m1
   DOUBLE PRECISION, INTENT(IN)    :: A2
   DOUBLE PRECISION, INTENT(IN)    :: m2
   DOUBLE PRECISION, INTENT(IN)    :: DK0
   DOUBLE PRECISION, INTENT(IN)    :: K1C
   DOUBLE PRECISION, INTENT(IN)    :: fy
   DOUBLE PRECISION, INTENT(IN)    :: fu
   DOUBLE PRECISION, INTENT(IN)    :: Smax
   DOUBLE PRECISION, INTENT(IN)    :: dK
   DOUBLE PRECISION, INTENT(IN)    :: a
   DOUBLE PRECISION, INTENT(OUT)   :: dadn
!
! declarations
   DOUBLE PRECISION  :: ep, eq
   DOUBLE PRECISION  :: f, A0, dKth, dKeff, Kmax, curve
!
! Fatigue crack growth curves (N,mm)
! material: steel in air (Modified from BS 7910)
!
   ep = 2
   eq = 2
!
! Crack growth laws
   IF (law.EQ.2) THEN
!
! Following 2-stage law BS7910:2013
      IF (DK.LE.DK0) THEN
	     dadn = 0.D0
	  ELSE
	     dadn = MIN ((A2 * (DK**m2)),(A1 * (DK**m1)))
      END IF
   ELSE IF (law.EQ.1) THEN
!
! Following Forman Mettu
      Kmax = dK / (1.D0 - R)
      CALL Parameter_f (Smax, fy, fu, f, A0)
      CALL CrackGrowthThreshold (DK0, a, f, A0, dKth)
      dKeff = (1.D0-f)/(1.D0-R)*dK
      IF (dK.LE.dKth) THEN
         dadn = 0.D0
      ELSE IF (Kmax.GE.K1C) THEN
         dadn = 1.D0   ! Hele grote waarde gegeven. Feitelijk treedt breuk op
      ELSE
         curve = ((1.D0-dKth/dK)**ep)/((1.D0-Kmax/K1C)**eq)
         dadn = A1*curve*(dKeff)**m1
      END IF
   END IF
!
! End subroutine
   RETURN
   END

