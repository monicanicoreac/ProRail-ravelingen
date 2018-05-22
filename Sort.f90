$ Declare
   SUBROUTINE Sort (loc, Cunc)
!
! Subroutine sorteren van een array in oplopende volgorde
!
! Parameters
	IMPLICIT NONE
   INTEGER,          INTENT(IN)    :: loc
   DOUBLE PRECISION, INTENT(INOUT) :: Cunc
!
! Declaraties
   INTEGER			:: I, J, L, IR
   DOUBLE PRECISION	:: RRA
!
! Initialisatie
   DIMENSION Cunc(loc)
IF (loc.nE.1) THEN
   L = loc / 2 + 1
   IR = loc
!
!
10	CONTINUE
	IF (L.GT.1) then
		L = L - 1
		RRA = Cunc(L)
	ELSE
		RRA = Cunc(IR)
		Cunc(IR) = Cunc(1)
		IR=IR - 1
		IF (IR.EQ.1) THEN
			Cunc(1) = RRA
			RETURN
		END IF
	END IF
	I = L
	J = L + L
20	IF (J.LE.IR) THEN
		IF (J.LT.IR) THEN
			IF (Cunc(J).LT.Cunc(J+1)) J = J + 1
		END IF
		IF (RRA.LT.Cunc(J)) THEN
			Cunc(I) = Cunc(J)
			I = J
			J = J + J
		ELSE
			J = IR + 1
		END IF
	    GO TO 20
	END IF
	Cunc(I) = RRA
    GO TO 10
END IF
!
! Einde routine
   RETURN
   END