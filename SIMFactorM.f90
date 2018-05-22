$ Declare
   SUBROUTINE SIMFactorM (M)
!
! Stress intensity magnification factor M
! M according to BS 7910:2005 M.3 Flat plates
!
! Parameters
   USE PublicPar
   IMPLICIT NONE
   DOUBLE PRECISION, INTENT(OUT)   :: M         ! Stress intensity magnification factor M
!
! M for cope
   IF (FlawType .EQ. 9) THEN
      M = 1.D0
!
! M for through-thickness deck plate crack
   ELSE IF (FlawType .EQ. 8) THEN
      M = 1.D0
!
! M for edge flaws in plates BS 7910:2005 M.3.5
   ELSE IF (FlawType .EQ. 7) THEN
      M = 1.D0
!
! M for surface flaws in bolts BS 7910:2005 M.6.3.2
   ELSE IF (FlawType .EQ. 6) THEN
      M = 1.D0
!
! M for corner flaws in plates
   ELSE IF (FlawType .EQ. 5) THEN
      M = 1.D0
!
! M for surface flaws in round bars BS 7910:2005 M.6
   ELSE IF (FlawType .EQ. 4) THEN
      M = 1.D0
!
! M for long surface flaws in plates BS 7910:2005 M.3.3
   ELSE IF (FlawType .EQ. 3) THEN
      M = 1.D0
!
! M for through-thickness flaws in plates BS 7910:2005 M.3.1
   ELSE IF (FlawType .EQ. 1) THEN
      M = 1.D0
!
! M for surface flaws in plates
   ELSE IF (FlawType .EQ. 2) THEN
      M = 1.D0
   ELSE
      STOP 'Incorrect input in Subroutine SIMFactorM: Unknown flawtype'
   END IF
!
! End subroutine
   RETURN
   END