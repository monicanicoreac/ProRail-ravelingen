$ Declare
   SUBROUTINE StressFact (DetailType, c, Sfact)
!
! Stress intensity correction factor Y according to BS 7910:2005 M.3 Flat plates
!
! Parameters
	IMPLICIT NONE
   INTEGER,          INTENT(IN)    :: DetailType   ! 1 = Detail 11b
                                                   ! 2 = Detail 12
												   ! 3 = Detail 4a
   DOUBLE PRECISION, INTENT(IN)    :: c            ! half crack width
   DOUBLE PRECISION, INTENT(OUT)   :: Sfact        ! Factor with which the stress needs to be multiplied
!
! Declarations
   DOUBLE PRECISION                :: SCF
!
! Set stress concentration factor for flaw type 4 (detail 12)
   SCF = 1.5D0
!
! Calculate stress factor
   IF (DetailType .EQ. 2) THEN
	  Sfact = 1.D0 / SCF + (1.D0 - 1.D0 / SCF) * DEXP(-1.D0 * c)
   ELSE
      Sfact = 1.D0
   END IF
!
! End subroutine
   RETURN
   END