$ Declare
   SUBROUTINE Parameter_f (Smax, fy, fu, f, A0)
!
! parameters
   USE PublicPar
   IMPLICIT NONE
!
   DOUBLE PRECISION, INTENT(IN)    :: Smax, fy, fu
   DOUBLE PRECISION, INTENT(OUT)   :: f, A0
!
! declarations
   DOUBLE PRECISION  :: pi
   DOUBLE PRECISION  :: Sflow
   DOUBLE PRECISION  :: A1, A2, A3, f1, f2, fK
   DOUBLE PRECISION  :: ratio
!
! material: steel in air (NasGro equations, parameters
! acc. to Fitnet and German dissertation)
!
! Constants and parameters
   pi    = 3.14159265358979D0
!
   Sflow = (fy + fu) / (2*fy)
   ratio = MIN (Smax / Sflow, 0.5D0)
   A0 = (0.825D0-0.34D0*alpha+0.05D0*alpha**2D0)*(DCOS(pi/2.D0*ratio))**(1.D0/alpha)
   A1 = (0.415D0-0.071D0*alpha)*ratio
   A3 = 2.D0*A0+A1-1.D0
   A2 = 1.D0-A0-A1-A3
   f1 = A0 + A1*R + A2*R**2 + A3*R**3
   f2 = A0 + A1*R
!
   IF (R.GE.0.) THEN
     fK = MAX (R,f1)
   ELSE IF (R.GE.-2.) THEN
     fK = f2
   ELSE IF (R.LT.-2.) THEN
     fK = A0 - 2*A1
   END IF
!
   IF (fK.GE.0.) THEN
     f = fK
   ELSE
     f = 0.
   END IF	   
!
! End subroutine
   RETURN
   END