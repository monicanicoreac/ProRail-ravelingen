$ Declare
   SUBROUTINE SICFactorY (LoadingType, a, c, B, W, L, tw, phi, rho, theta, Y, err)
!
! Stress intensity correction factor Y according to BS 7910:2005 M.3 Flat plates
!
! Parameters
   USE PublicPar
   IMPLICIT NONE
!
   INTEGER, INTENT(IN)             :: LoadingType  ! 1 = Membrane loading
                                                   ! 2 = Bending loading
   DOUBLE PRECISION, INTENT(IN)    :: a            ! Crack depth
   DOUBLE PRECISION, INTENT(IN)    :: c            ! Crack width
   DOUBLE PRECISION, INTENT(IN)    :: B            ! Plate thickness
   DOUBLE PRECISION, INTENT(IN)    :: W            ! Plate width
   DOUBLE PRECISION, INTENT(IN)    :: L            ! Attachment length [mm]
   DOUBLE PRECISION, INTENT(IN)    :: tw           ! WG = 2: Weld throat thickness [mm]
   DOUBLE PRECISION, INTENT(IN)	   :: phi		   ! Weld toe angle in degrees
   DOUBLE PRECISION, INTENT(IN)	   :: rho		   ! Weld toe radius in mm
   DOUBLE PRECISION, INTENT(IN)    :: theta        ! Crack angle in radians (0 => c, pi/2 => a)
   DOUBLE PRECISION, INTENT(OUT)   :: Y            ! Geometry factor
   INTEGER, INTENT(OUT)            :: err          ! out of validity for one of the input parameters
!
! Declarations
   DOUBLE PRECISION                :: fw
   DOUBLE PRECISION                :: M
   DOUBLE PRECISION                :: Mscf
   DOUBLE PRECISION                :: Mm
   DOUBLE PRECISION                :: Mmmax
   DOUBLE PRECISION                :: Mb
   DOUBLE PRECISION                :: Mk
   DOUBLE PRECISION                :: z
!
! Voorwaarden
   CALL SIMConditions (a, c, B, W, err)
   IF (err.EQ.0) THEN
!
! Calculate stress intensity magnification factor M 
      CALL SIMFactorM (M)
!
! Factor Fw
      CALL Factorfw (a, c, B, W, fw)
!
! Calculate weld geometry factor Mk
      IF (WeldGeometry .EQ. 0 ) THEN
         Mk = 1.D0
      ELSE
         IF (theta .EQ. 0 ) THEN
	        z = 0.15D0
         ELSE
	        z = a
         END IF
         CALL WeldGeometryFactorMk (LoadingType, z, B, L, tw, phi, rho, Mk)
      END IF
!
! Calculate stress intensity magnification factor Mm or Mb
! Stress intensity correction factor Y = M fw Mm or Y = M fw Mb 
      IF (LoadingType .EQ. 1) THEN
         CALL SIMFactorMm (a, c, B, W, theta, Mm, Mmmax)
	     Y = M * fw * Mm * Mk
      ELSE IF (LoadingType .EQ. 2) THEN
         CALL SIMFactorMb (a, c, B, W, theta, Mb)
	     Y = M * fw * Mb * Mk
      ELSE
         STOP 'Error in Subroutine SICFactorY: Unknown loadingtype'
      END IF
   END IF
!
! End subroutine
   RETURN
   END